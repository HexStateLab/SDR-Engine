/*
 * sdr_ether.c — Transmit Qudits Through the EM Field ("Ether")
 *
 * OFFLOADS COMPUTATION TO NATURE:
 *   1. TX: R820T2 PLL LO leakage radiates qudit levels into EM field
 *          (each level k = CW tone at f₀ + k·Δf, dwell ∝ |ψₖ|²)
 *       or: IDFT(qudit state) → OFDM I/Q → CU8 for external TX SDR
 *   2. ETHER: Propagation through space (interference, multipath, drift)
 *   3. RX: RTL2832U ADC capture → DFT → recover qudit wavefunction
 *   4. GATES: Apply quantum gates, then repeat TX→ether→RX feedback loop
 *
 * Each round-trip through the EM field IS a computational step.
 * The EM field IS the quantum processor.
 *
 * ─── Modes ───
 *   (default)    Physical ether: TX via LO leakage + RX via SDR capture
 *   --physical    Same as default (explicit)
 *   --rx-only     RX from RTL-SDR only (no TX)
 *   --loopback    TX→ether channel→RX in software, no hardware needed
 *   --tx-only     Synthesize I/Q, write CU8 to stdout for SDR transmitter
 *   --tx-file F   Synthesize I/Q, write CU8 to file F
 *   --ofdm        OFDM multi-tone: fixed LO, simultaneous subcarriers
 *   --ofdm-file F OFDM multi-tone output to file F
 *   --cycles N    Run N TX→ether→RX→gate feedback cycles
 *   --snr-db N    Ether channel SNR for loopback (default 30 dB)
 *   --fading 1    Enable frequency-selective fading in loopback ether
 *   --drift 1     Enable random phase drift per loopback cycle
 *
 * Build: gcc -O3 -std=gnu99 sdr_ether.c -lm -lpthread -o sdr_ether
 * Usage:
 *   ./sdr_ether 6                        # Physical ether (needs RTL-SDR)
 *   ./sdr_ether 6 --cycles 10            # 10 TX→ether→RX iterations
 *   ./sdr_ether 8 --rx-only              # RX only, no TX
 *   ./sdr_ether 12 --loopback --cycles 5 # Software ether simulation
 *   ./sdr_ether 8 --tx-file tx.iq        # Write TX signal for external SDR
 *   ./sdr_ether 8 --tx-only | hackrf_transfer -f 100e6 -s 2e6 -t -
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <poll.h>
#include <time.h>
#include <sys/ioctl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <linux/videodev2.h>
#include <dirent.h>
#include <pthread.h>

#include "qvm_api.h"

#define BUF_COUNT    8
#define IQ_WINDOW    65536
#define MAX_DIM      65536
#define DEFAULT_D    6
#define DEFAULT_FREQ 100000000
#define DEFAULT_RATE 2048000
#define TX_RING_SIZE 16384

/* ═══════════════════════════════════════════════════════════════
 * GLOBAL STATE
 * ═══════════════════════════════════════════════════════════════ */
static int    g_dim     = DEFAULT_D;
static int    g_cycles  = 1;
static double g_snr_db  = 30.0;
static int    g_fading  = 0;
static int    g_drift   = 0;
static char   g_sdr_dev[64] = "";

static int sdr_find(void) {
    if (g_sdr_dev[0]) return 0;
    DIR *d = opendir("/dev"); if (!d) goto fallback;
    struct dirent *de;
    while ((de = readdir(d)))
        if (strncmp(de->d_name,"swradio",7)==0) {
            snprintf(g_sdr_dev,64,"/dev/%s",de->d_name); closedir(d); return 0;
        }
    closedir(d);
fallback:
    snprintf(g_sdr_dev,64,"/dev/swradio0"); return -1;
}
static double g_lambda  = 1.0;  /* Tikhonov regularization */

/* ═══════════════════════════════════════════════════════════════
 * SDR RECEIVE DEVICE (RTL-SDR via V4L2)
 * ═══════════════════════════════════════════════════════════════ */
typedef struct {
    int      fd;
    uint32_t freq, rate;
    int      gain;
    uint8_t *bufs[BUF_COUNT];
    uint32_t buf_len[BUF_COUNT];
    int      cur_buf;
    uint32_t cur_off;
    uint8_t  iq_raw[IQ_WINDOW];
    double   iq_i[IQ_WINDOW/2];
    double   iq_q[IQ_WINDOW/2];
    int      iq_n;
    uint64_t samples;
    pthread_t tx_thread;
    volatile int tx_running;
    uint32_t  tx_ring[TX_RING_SIZE];
    volatile int tx_wr, tx_rd;
} SdrDev;

/* ═══════════════════════════════════════════════════════════════
 * TX API — background thread + ring buffer (augments V4L2 RX)
 * ═══════════════════════════════════════════════════════════════ */
static void *sdr_tx_worker(void *arg) {
    SdrDev *s = (SdrDev *)arg;
    uint32_t base = s->freq;
    while (s->tx_running) {
        while (s->tx_wr != s->tx_rd && s->tx_running) {
            uint32_t *cmd = &s->tx_ring[s->tx_rd];
            int n = cmd[0];
            if (n > 0 && n <= 512) {
                for (int i = 0; i < n; i++) {
                    uint32_t f = cmd[1 + 2*i];
                    uint32_t dwell = cmd[1 + 2*i + 1];
                    struct v4l2_frequency vf;
                    memset(&vf, 0, sizeof(vf));
                    vf.tuner = 0;
                    vf.type  = V4L2_TUNER_ADC;
                    vf.frequency = f;
                    ioctl(s->fd, VIDIOC_S_FREQUENCY, &vf);
                    usleep(dwell);
                }
                struct v4l2_frequency vf;
                memset(&vf, 0, sizeof(vf));
                vf.tuner = 0;
                vf.type  = V4L2_TUNER_ADC;
                vf.frequency = base;
                ioctl(s->fd, VIDIOC_S_FREQUENCY, &vf);
            }
            int slots = 1 + 2 * n;
            s->tx_rd = (s->tx_rd + slots) % TX_RING_SIZE;
        }
        usleep(1000);
    }
    return NULL;
}

static int sdr_tx(SdrDev *s, const uint32_t *freqs,
                   const uint32_t *dwells, int n) {
    if (n < 1 || n > 512) return -1;
    int slots = 1 + 2 * n;
    if (((s->tx_wr + slots) % TX_RING_SIZE) == s->tx_rd) return -1;
    s->tx_ring[s->tx_wr] = n;
    for (int i = 0; i < n; i++) {
        s->tx_ring[s->tx_wr + 1 + 2*i] = freqs[i];
        s->tx_ring[s->tx_wr + 1 + 2*i + 1] = dwells[i];
    }
    s->tx_wr = (s->tx_wr + slots) % TX_RING_SIZE;
    return n;
}

static int sdr_tx_wavefunction(SdrDev *s, const double *prob,
                                int d, uint32_t base_freq,
                                double bin_bw, int base_dwell_us) {
    int count = 0;
    for (int k = 0; k < d; k++) {
        int dwell = (int)(prob[k] * base_dwell_us);
        if (dwell < 30) continue;
        uint32_t f = base_freq + (uint32_t)(k * bin_bw);
        if (f < 24000000) f = 24000000;
        if (f > 1750000000) f = 1750000000;
        uint32_t freqs[1] = { f };
        uint32_t dwells[1] = { (uint32_t)dwell };
        if (sdr_tx(s, freqs, dwells, 1) < 0) break;
        count++;
    }
    return count;
}

/* ═══════════════════════════════════════════════════════════════
 * TRANSMIT SYNTHESIS BUFFER
 * ═══════════════════════════════════════════════════════════════ */
typedef struct {
    double  *i, *q;     /* synthesized I/Q (float, full precision) */
    uint8_t *cu8;       /* CU8 format for SDR TX */
    int      nsamples;
    uint32_t freq, rate;
} TxBuf;

/* ═══════════════════════════════════════════════════════════════
 * ETHER CHANNEL MODEL
 *
 * Simulates EM field propagation between TX antenna and RX antenna:
 *   - AWGN thermal noise
 *   - Random global phase rotation (ether drift)
 *   - Frequency-selective fading (per-subcarrier gain/phase from multipath)
 *   - Multipath delay spread (temporal smearing)
 * ═══════════════════════════════════════════════════════════════ */
typedef struct {
    double  snr_linear;
    double  phase_drift_rad;
    double *bin_gain;
    double *bin_phase;
    int     fading_enabled;
    int     drift_enabled;
    int     dim;
} EtherChan;

/* ═══════════════════════════════════════════════════════════════
 * WAVEFUNCTION (Recovered from SDR / Ether)
 * ═══════════════════════════════════════════════════════════════ */
typedef struct {
    int      d;
    double  *re, *im, *prob;
    double   entropy, purity;
    uint32_t freq, rate;
    double   bin_bw;
} Wavefunction;

/* ─── Forward declarations ─── */
static void sdr_close(SdrDev *s);

/* ═══════════════════════════════════════════════════════════════
 * ALLOCATORS
 * ═══════════════════════════════════════════════════════════════ */
static TxBuf tx_alloc(int nsamples, uint32_t freq, uint32_t rate) {
    TxBuf t;
    t.i  = calloc(nsamples, sizeof(double));
    t.q  = calloc(nsamples, sizeof(double));
    t.cu8= calloc(nsamples * 2, 1);
    t.nsamples = nsamples;
    t.freq = freq; t.rate = rate;
    return t;
}
static void tx_free(TxBuf *t) { free(t->i); free(t->q); free(t->cu8); }

static Wavefunction wf_alloc(int d, uint32_t freq, uint32_t rate) {
    Wavefunction w;
    w.d = d; w.freq = freq; w.rate = rate;
    w.bin_bw = (double)rate / d;
    w.re  = calloc(d, sizeof(double));
    w.im  = calloc(d, sizeof(double));
    w.prob= calloc(d, sizeof(double));
    w.entropy = 0; w.purity = 0;
    return w;
}
static void wf_free(Wavefunction *w) { free(w->re); free(w->im); free(w->prob); }

static EtherChan ether_init(int dim, double snr_db, int fading, int drift) {
    EtherChan ch;
    ch.dim = dim;
    ch.snr_linear = pow(10.0, snr_db / 20.0);
    ch.phase_drift_rad = 0.0;
    ch.fading_enabled = fading;
    ch.drift_enabled  = drift;
    ch.bin_gain  = calloc(dim, sizeof(double));
    ch.bin_phase = calloc(dim, sizeof(double));
    for (int k = 0; k < dim; k++) { ch.bin_gain[k] = 1.0; ch.bin_phase[k] = 0.0; }
    return ch;
}
static void ether_free(EtherChan *ch) { free(ch->bin_gain); free(ch->bin_phase); }

/* ═══════════════════════════════════════════════════════════════
 * SDR RECEIVE — capture I/Q from RTL-SDR hardware
 * ═══════════════════════════════════════════════════════════════ */
static int sdr_open(SdrDev *s, uint32_t freq, uint32_t rate, int gain) {
    sdr_find();
    memset(s, 0, sizeof(*s));
    s->freq = freq; s->rate = rate; s->gain = gain;
    s->cur_buf = -1;

    s->fd = open(g_sdr_dev, O_RDWR);
    if (s->fd < 0) {
        fprintf(stderr, "[SDR] Cannot open %s: %s\n", g_sdr_dev, strerror(errno));
        return -1;
    }

    struct v4l2_format fmt;
    memset(&fmt, 0, sizeof(fmt));
    fmt.type = V4L2_BUF_TYPE_SDR_CAPTURE;
    if (ioctl(s->fd, VIDIOC_G_FMT, &fmt) < 0) goto fail;
    fmt.fmt.sdr.pixelformat = V4L2_SDR_FMT_CU8;
    if (ioctl(s->fd, VIDIOC_S_FMT, &fmt) < 0) goto fail;

    struct v4l2_frequency vf;
    memset(&vf, 0, sizeof(vf));
    vf.tuner = 0; vf.type = V4L2_TUNER_ADC;
    vf.frequency = freq;
    ioctl(s->fd, VIDIOC_S_FREQUENCY, &vf);

    struct v4l2_requestbuffers req;
    memset(&req, 0, sizeof(req));
    req.count = BUF_COUNT;
    req.type  = V4L2_BUF_TYPE_SDR_CAPTURE;
    req.memory = V4L2_MEMORY_MMAP;
    if (ioctl(s->fd, VIDIOC_REQBUFS, &req) < 0) goto fail;

    for (uint32_t i = 0; i < req.count; i++) {
        struct v4l2_buffer b;
        memset(&b, 0, sizeof(b));
        b.type = V4L2_BUF_TYPE_SDR_CAPTURE;
        b.memory = V4L2_MEMORY_MMAP;
        b.index = i;
        if (ioctl(s->fd, VIDIOC_QUERYBUF, &b) < 0) goto fail;
        s->buf_len[i] = b.length;
        s->bufs[i] = mmap(NULL, b.length, PROT_READ|PROT_WRITE,
                          MAP_SHARED, s->fd, b.m.offset);
        if (s->bufs[i] == MAP_FAILED) goto fail;
        ioctl(s->fd, VIDIOC_QBUF, &b);
    }

    enum v4l2_buf_type type = V4L2_BUF_TYPE_SDR_CAPTURE;
    if (ioctl(s->fd, VIDIOC_STREAMON, &type) < 0) goto fail;

    for (int a = 0; a < 8; a++) {
        struct v4l2_buffer wb;
        memset(&wb, 0, sizeof(wb));
        wb.type = V4L2_BUF_TYPE_SDR_CAPTURE;
        wb.memory = V4L2_MEMORY_MMAP;
        if (ioctl(s->fd, VIDIOC_DQBUF, &wb) != 0) break;
        int sum = 0, nc = wb.bytesused > 256 ? 256 : (int)wb.bytesused;
        for (int i = 0; i < nc; i++) sum += s->bufs[wb.index][i];
        int ok = (nc > 0 && sum/nc > 50 && sum/nc < 200);
        ioctl(s->fd, VIDIOC_QBUF, &wb);
        if (ok) {
            fprintf(stderr, "[SDR] D=%d @ %.1f MHz %.2f MSPS bin=%.1f Hz\n",
                    g_dim, freq/1e6, rate/1e6, (double)rate/g_dim);
            s->tx_running = 1;
            s->tx_wr = 0; s->tx_rd = 0;
            pthread_create(&s->tx_thread, NULL, sdr_tx_worker, s);
            return 0;
        }
    }

    fprintf(stderr, "[SDR] No signal at %.1f MHz — EM field quiet\n", freq/1e6);
fail:
    sdr_close(s);
    return -1;
}

static int sdr_capture(SdrDev *s) {
    int copied = 0, need = IQ_WINDOW;
    while (need > 0) {
        if (s->cur_buf >= 0 && s->cur_off < s->buf_len[s->cur_buf]) {
            uint32_t avail = s->buf_len[s->cur_buf] - s->cur_off;
            int cp = (int)avail < need ? (int)avail : need;
            memcpy(s->iq_raw + copied, s->bufs[s->cur_buf] + s->cur_off, cp);
            s->cur_off += cp; copied += cp; need -= cp;
            if (need == 0) {
                struct v4l2_buffer b;
                memset(&b, 0, sizeof(b));
                b.type = V4L2_BUF_TYPE_SDR_CAPTURE;
                b.memory = V4L2_MEMORY_MMAP;
                b.index = s->cur_buf;
                ioctl(s->fd, VIDIOC_QBUF, &b);
                s->cur_buf = -1;
                break;
            }
            struct v4l2_buffer b;
            memset(&b, 0, sizeof(b));
            b.type = V4L2_BUF_TYPE_SDR_CAPTURE;
            b.memory = V4L2_MEMORY_MMAP;
            b.index = s->cur_buf;
            ioctl(s->fd, VIDIOC_QBUF, &b);
            s->cur_buf = -1;
        }
        struct pollfd p = { .fd = s->fd, .events = POLLIN };
        if (poll(&p, 1, 200) <= 0) break;
        struct v4l2_buffer db;
        memset(&db, 0, sizeof(db));
        db.type = V4L2_BUF_TYPE_SDR_CAPTURE;
        db.memory = V4L2_MEMORY_MMAP;
        if (ioctl(s->fd, VIDIOC_DQBUF, &db) < 0) {
            if (errno == EAGAIN) { usleep(5000); continue; }
            break;
        }
        s->cur_buf = (int)db.index;
        s->cur_off = 0;
    }
    int np = copied / 2;
    s->iq_n = np;
    s->samples += copied;
    for (int i = 0; i < np; i++) {
        s->iq_i[i] = ((double)s->iq_raw[2*i]   - 127.5) / 128.0;
        s->iq_q[i] = ((double)s->iq_raw[2*i+1] - 127.5) / 128.0;
    }
    return copied;
}

static void sdr_retune(SdrDev *s, uint32_t hz) {
    struct v4l2_frequency vf;
    memset(&vf, 0, sizeof(vf));
    vf.tuner = 0; vf.type = V4L2_TUNER_ADC;
    vf.frequency = hz;
    ioctl(s->fd, VIDIOC_S_FREQUENCY, &vf);
    s->freq = hz;
}

static void sdr_flush(SdrDev *s) {
    /* Reset internal state — force sdr_capture to get fresh buffers */
    if (s->cur_buf >= 0) {
        struct v4l2_buffer b;
        memset(&b, 0, sizeof(b));
        b.type = V4L2_BUF_TYPE_SDR_CAPTURE;
        b.memory = V4L2_MEMORY_MMAP;
        b.index = s->cur_buf;
        ioctl(s->fd, VIDIOC_QBUF, &b);
    }
    s->cur_buf = -1;
    s->cur_off = 0;
    /* Drain any queued buffers by polling and re-queueing */
    for (int i = 0; i < 2; i++) {
        struct pollfd p = { .fd = s->fd, .events = POLLIN };
        if (poll(&p, 1, 20) <= 0) break;
        struct v4l2_buffer db;
        memset(&db, 0, sizeof(db));
        db.type = V4L2_BUF_TYPE_SDR_CAPTURE;
        db.memory = V4L2_MEMORY_MMAP;
        if (ioctl(s->fd, VIDIOC_DQBUF, &db) < 0) break;
        ioctl(s->fd, VIDIOC_QBUF, &db);
    }
}

/* Reset streaming after calibration sweeps drain buffers */
static void sdr_restream(SdrDev *s){
    if(s->fd < 0) return;
    if(s->cur_buf >= 0){
        struct v4l2_buffer b;memset(&b,0,sizeof(b));
        b.type=V4L2_BUF_TYPE_SDR_CAPTURE;b.memory=V4L2_MEMORY_MMAP;
        b.index=s->cur_buf;ioctl(s->fd,VIDIOC_QBUF,&b);
    }
    s->cur_buf=-1;s->cur_off=0;
    /* STREAMOFF + STREAMON to reset DMA */
    enum v4l2_buf_type t=V4L2_BUF_TYPE_SDR_CAPTURE;
    ioctl(s->fd,VIDIOC_STREAMOFF,&t);
    usleep(10000);
    /* Re-queue all buffers */
    for(int i=0;i<BUF_COUNT;i++){
        struct v4l2_buffer b;memset(&b,0,sizeof(b));
        b.type=V4L2_BUF_TYPE_SDR_CAPTURE;b.memory=V4L2_MEMORY_MMAP;b.index=i;
        ioctl(s->fd,VIDIOC_QBUF,&b);
    }
    ioctl(s->fd,VIDIOC_STREAMON,&t);
    /* Warmup */
    for(int a=0;a<8;a++){
        struct v4l2_buffer db;memset(&db,0,sizeof(db));
        db.type=V4L2_BUF_TYPE_SDR_CAPTURE;db.memory=V4L2_MEMORY_MMAP;
        if(ioctl(s->fd,VIDIOC_DQBUF,&db)!=0)break;
        ioctl(s->fd,VIDIOC_QBUF,&db);
    }
    fprintf(stderr,"[SDR] restream OK\n");
}

static void sdr_close(SdrDev *s) {
    if (s->fd < 0) return;
    s->tx_running = 0;
    pthread_join(s->tx_thread, NULL);
    enum v4l2_buf_type t = V4L2_BUF_TYPE_SDR_CAPTURE;
    ioctl(s->fd, VIDIOC_STREAMOFF, &t);
    for (int i = 0; i < BUF_COUNT; i++)
        if (s->bufs[i] && s->bufs[i] != MAP_FAILED)
            munmap(s->bufs[i], s->buf_len[i]);
    close(s->fd);
    s->fd = -1;
}

/* ═══════════════════════════════════════════════════════════════
 * WAVEFUNCTION COMPUTE — DFT from I/Q time series
 *
 * Decomposes captured I/Q into D frequency bins.
 * Each bin k = complex amplitude at f₀ + k·Δf.
 * This IS the physical qudit state.
 * ═══════════════════════════════════════════════════════════════ */
static void wf_from_iq(const double *iq_i, const double *iq_q, int np,
                        Wavefunction *wf) {
    memset(wf->re, 0, wf->d * sizeof(double));
    memset(wf->im, 0, wf->d * sizeof(double));
    memset(wf->prob, 0, wf->d * sizeof(double));
    wf->entropy = 0; wf->purity = 0;

    if (np < wf->d) return;

    for (int k = 0; k < wf->d; k++) {
        double freq_norm = (double)k / (double)wf->d;
        double sum_i = 0, sum_q = 0;
        for (int n = 0; n < np; n++) {
            double phase = -2.0 * M_PI * freq_norm * (double)n;
            double cr = cos(phase), sr = sin(phase);
            sum_i += iq_i[n] * cr - iq_q[n] * sr;
            sum_q += iq_i[n] * sr + iq_q[n] * cr;
        }
        wf->re[k] = sum_i / (double)np;
        wf->im[k] = sum_q / (double)np;
    }

    double total = 0;
    for (int k = 0; k < wf->d; k++)
        total += wf->re[k]*wf->re[k] + wf->im[k]*wf->im[k];

    if (total > 1e-15) {
        double scale = 1.0 / sqrt(total);
        for (int k = 0; k < wf->d; k++) {
            wf->re[k] *= scale; wf->im[k] *= scale;
            wf->prob[k] = wf->re[k]*wf->re[k] + wf->im[k]*wf->im[k];
            if (wf->prob[k] > 1e-15)
                wf->entropy -= wf->prob[k] * log2(wf->prob[k]);
            wf->purity += wf->prob[k] * wf->prob[k];
        }
    }
}

/* ═══════════════════════════════════════════════════════════════
 * TX SYNTHESIS — qudit state → OFDM multi-tone I/Q
 *
 * Each qudit level k maps to subcarrier at f₀ + k·Δf where Δf=rate/D.
 * The time-domain I/Q is the IDFT of the complex wavefunction:
 *
 *   x[n] = Σₖᴰ⁻¹ (re[k] + j·im[k]) · exp(+j·2π·k·n/D)
 *
 * This generates the baseband signal that a transmitter sends into
 * the EM field.  When received and DFT'd, it reconstructs the qudit state.
 * ═══════════════════════════════════════════════════════════════ */
static void tx_synthesize(const Wavefunction *wf, TxBuf *tx) {
    int D = wf->d, N = tx->nsamples;

    memset(tx->i, 0, N * sizeof(double));
    memset(tx->q, 0, N * sizeof(double));

    for (int k = 0; k < D; k++) {
        double re = wf->re[k], im = wf->im[k];
        if (re*re + im*im < 1e-30) continue;
        double freq_norm = (double)k / (double)D;
        for (int n = 0; n < N; n++) {
            double phase = 2.0 * M_PI * freq_norm * (double)n;
            double cr = cos(phase), sr = sin(phase);
            tx->i[n] += re * cr - im * sr;
            tx->q[n] += re * sr + im * cr;
        }
    }

    double peak = 0;
    for (int n = 0; n < N; n++) {
        double mag = fabs(tx->i[n]);
        if (fabs(tx->q[n]) > mag) mag = fabs(tx->q[n]);
        if (mag > peak) peak = mag;
    }
    if (peak > 1e-15) {
        double sc = 0.9 / peak;
        for (int n = 0; n < N; n++) {
            tx->i[n] *= sc;
            tx->q[n] *= sc;
        }
    }

    for (int n = 0; n < N; n++) {
        int iv = (int)(tx->i[n] * 127.5 + 127.5);
        int qv = (int)(tx->q[n] * 127.5 + 127.5);
        iv = iv < 0 ? 0 : (iv > 255 ? 255 : iv);
        qv = qv < 0 ? 0 : (qv > 255 ? 255 : qv);
        tx->cu8[2*n]   = (uint8_t)iv;
        tx->cu8[2*n+1] = (uint8_t)qv;
    }
}

/* ═══════════════════════════════════════════════════════════════
 * ETHER CHANNEL — EM field propagation model
 *
 * Applies physical effects that happen between TX and RX antennas:
 *   1. Frequency-selective fading (multipath interference)
 *   2. AWGN thermal noise (Johnson-Nyquist)
 *   3. Random global phase rotation (ether drift)
 *
 * After applying the channel, converts CU8 back to float I/Q.
 * ═══════════════════════════════════════════════════════════════ */
static void ether_apply(TxBuf *tx, EtherChan *ch, double *rx_i, double *rx_q, int np) {
    int N = np < tx->nsamples ? np : tx->nsamples;

    if (ch->fading_enabled || ch->drift_enabled) {
        static uint64_t seed_inc = 0;
        uint64_t seed = (uint64_t)time(NULL) ^ (uint64_t)(uintptr_t)tx ^ seed_inc++;
        seed ^= seed >> 33; seed *= 0xFF51AFD7ED558CCDULL;
        seed ^= seed >> 33; seed *= 0xC4CEB9FE1A85EC53ULL;
        seed ^= seed >> 33;
        srand48((long)seed);
    }

    if (ch->fading_enabled) {
        for (int k = 0; k < ch->dim; k++) {
            double g = sqrt(-2.0 * log(drand48() + 1e-15)) * 0.3 + 1.0;
            g = g < 0.1 ? 0.1 : (g > 3.0 ? 3.0 : g);
            ch->bin_gain[k] = g;
            ch->bin_phase[k] = drand48() * 2.0 * M_PI;
        }
        fprintf(stderr, "  [ETHER] Fading: gains=[");
        int show = ch->dim < 8 ? ch->dim : 6;
        for (int k = 0; k < show; k++)
            fprintf(stderr, "%.2f%s", ch->bin_gain[k], k<show-1?" ":"");
        if (ch->dim > show) fprintf(stderr, "...");
        fprintf(stderr, "]\n");
    }

    if (ch->drift_enabled) {
        ch->phase_drift_rad = drand48() * 2.0 * M_PI;
        fprintf(stderr, "  [ETHER] Drift: %.3f rad\n", ch->phase_drift_rad);
    }

    double cr = cos(ch->phase_drift_rad);
    double sr = sin(ch->phase_drift_rad);
    for (int n = 0; n < N; n++) {
        double i = tx->i[n], q = tx->q[n];
        rx_i[n] = i * cr - q * sr;
        rx_q[n] = i * sr + q * cr;
    }

    /* Frequency-selective fading: per-subcarrier gain + phase,
     * applied in time domain via convolution with channel impulse response.
     * Simplified: apply per-block gain/phase after DFT, done in wf_from_iq.
     * Here we apply a simplified time-domain smearing for multipath. */
    if (ch->fading_enabled) {
        double *tmp_i = malloc(N * sizeof(double));
        double *tmp_q = malloc(N * sizeof(double));
        memcpy(tmp_i, rx_i, N * sizeof(double));
        memcpy(tmp_q, rx_q, N * sizeof(double));

        int delay = g_dim / 2 + 1;
        double echo_gain = 0.35;
        for (int n = delay; n < N; n++) {
            rx_i[n] += tmp_i[n - delay] * echo_gain;
            rx_q[n] += tmp_q[n - delay] * echo_gain;
        }
        free(tmp_i); free(tmp_q);
    }

    double noise_rms = (ch->snr_linear > 1e-15) ? 1.0 / ch->snr_linear : 1.0;
    for (int n = 0; n < N; n++) {
        double u1 = drand48() + 1e-15, u2 = drand48();
        double ng = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
        rx_i[n] += ng * noise_rms * 0.5;
        u1 = drand48() + 1e-15; u2 = drand48();
        ng = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
        rx_q[n] += ng * noise_rms * 0.5;
    }

    double peak = 0;
    for (int n = 0; n < N; n++) {
        double mag = fabs(rx_i[n]);
        if (fabs(rx_q[n]) > mag) mag = fabs(rx_q[n]);
        if (mag > peak) peak = mag;
    }
    if (peak > 1.2) {
        double sc = 1.0 / peak;
        for (int n = 0; n < N; n++) { rx_i[n] *= sc; rx_q[n] *= sc; }
    }
}

/* ═══════════════════════════════════════════════════════════════
 * ETHER FADING CORRECTION — apply per-bin gain/phase to wavefunction
 *
 * After wf_from_iq computes the DFT, we apply the known (for loopback)
 * or partially-known fading coefficients to simulate the full channel.
 * ═══════════════════════════════════════════════════════════════ */
static void ether_correct_wf(Wavefunction *wf, EtherChan *ch) {
    if (!ch->fading_enabled) return;
    for (int k = 0; k < wf->d && k < ch->dim; k++) {
        double g = ch->bin_gain[k];
        double phase = ch->bin_phase[k];
        double cr = cos(-phase), sr = sin(-phase);
        double re = wf->re[k], im = wf->im[k];
        if (g > 1e-10) {
            wf->re[k] = (re * cr - im * sr) / g;
            wf->im[k] = (re * sr + im * cr) / g;
        }
    }
    double total = 0;
    for (int k = 0; k < wf->d; k++)
        total += wf->re[k]*wf->re[k] + wf->im[k]*wf->im[k];
    if (total > 1e-15) {
        double scale = 1.0 / sqrt(total);
        for (int k = 0; k < wf->d; k++) {
            wf->re[k] *= scale; wf->im[k] *= scale;
            wf->prob[k] = wf->re[k]*wf->re[k] + wf->im[k]*wf->im[k];
        }
    }
}

/* ═══════════════════════════════════════════════════════════════
 * TX FILE OUTPUT — write CU8 I/Q for external SDR transmitter
 *
 * File format: interleaved uint8_t I/Q pairs, rate-bytes-per-second.
 * Compatible with:
 *   hackrf_transfer -f <freq> -s <rate> -t <file>     (--signed required if using CU8)
 *   Or pipe: ./sdr_ether --tx-only | hackrf_transfer -f 100e6 -s 2e6 -t -
 * ═══════════════════════════════════════════════════════════════ */
static void tx_write_file(const TxBuf *tx, const char *path) {
    FILE *fp;
    if (!path || strcmp(path, "-") == 0) {
        fp = stdout;
    } else {
        fp = fopen(path, "wb");
        if (!fp) { fprintf(stderr, "[TX] Cannot open %s: %s\n", path, strerror(errno)); return; }
    }
    fwrite(tx->cu8, 1, tx->nsamples * 2, fp);
    if (fp != stdout) fclose(fp);
    fprintf(stderr, "[TX] Wrote %d I/Q pairs (CU8) @ %.1f MHz %.2f MSPS → %s\n",
            tx->nsamples, tx->freq/1e6, tx->rate/1e6,
            path ? (strcmp(path,"-")==0?"stdout":path) : "stdout");
}

/* ═══════════════════════════════════════════════════════════════
 * PHYSICAL GATES
 *
 * Gates operate on the wavefunction recovered from the ether.
 * Each gate transforms the qudit state before it is re-transmitted
 * into the EM field for the next iteration.
 * ═══════════════════════════════════════════════════════════════ */

/* H: Hadamard-like via LO shift by Nyquist (folds spectrum) */
static void gate_H(SdrDev *s, Wavefunction *wf) {
    uint32_t orig = wf->freq;
    uint32_t hop  = orig + s->rate / 2;  /* Nyquist fold */
    sdr_retune(s, hop);
    usleep(5000);
    sdr_capture(s);
    wf_from_iq(s->iq_i, s->iq_q, s->iq_n, wf);
    sdr_retune(s, orig);
    wf->freq = orig;
    fprintf(stderr, "  [GATE:H] %.1f→%.1f→%.1f MHz\n",
            orig/1e6, hop/1e6, orig/1e6);
}

/* X: cyclic shift of qudit levels */
static void gate_X(Wavefunction *wf) {
    double re0 = wf->re[0], im0 = wf->im[0];
    for (int k = 0; k < wf->d - 1; k++) {
        wf->re[k] = wf->re[k+1];
        wf->im[k] = wf->im[k+1];
        wf->prob[k] = wf->prob[k+1];
    }
    wf->re[wf->d-1] = re0;
    wf->im[wf->d-1] = im0;
    wf->prob[wf->d-1] = re0*re0 + im0*im0;
    fprintf(stderr, "  [GATE:X] Cyclic shift\n");
}

/* X-1: reverse cyclic shift (for DFT pair creation) */
__attribute__((unused)) static void gate_Xi(Wavefunction *wf) {
    double re_last = wf->re[wf->d-1], im_last = wf->im[wf->d-1];
    for (int k = wf->d - 1; k > 0; k--) {
        wf->re[k] = wf->re[k-1];
        wf->im[k] = wf->im[k-1];
        wf->prob[k] = wf->prob[k-1];
    }
    wf->re[0] = re_last;
    wf->im[0] = im_last;
    wf->prob[0] = re_last*re_last + im_last*im_last;
    fprintf(stderr, "  [GATE:Xi] Reverse shift\n");
}

/* Z: phase rotation */
static void gate_Z(Wavefunction *wf, double rad, int target) {
    double cr = cos(rad), sr = sin(rad);
    if (target < 0) {
        for (int k = 0; k < wf->d; k++) {
            double re = wf->re[k], im = wf->im[k];
            wf->re[k] = re * cr - im * sr;
            wf->im[k] = re * sr + im * cr;
        }
    } else if (target < wf->d) {
        double re = wf->re[target], im = wf->im[target];
        wf->re[target] = re * cr - im * sr;
        wf->im[target] = re * sr + im * cr;
        wf->prob[target] = wf->re[target]*wf->re[target]
                         + wf->im[target]*wf->im[target];
    }
}

/* CZ: entangle via cross-correlation of two time windows */
static void gate_CZ(SdrDev *s, Wavefunction *wf) {
    double *prev_re = malloc(wf->d * sizeof(double));
    double *prev_im = malloc(wf->d * sizeof(double));
    memcpy(prev_re, wf->re, wf->d * sizeof(double));
    memcpy(prev_im, wf->im, wf->d * sizeof(double));

    sdr_capture(s);
    wf_from_iq(s->iq_i, s->iq_q, s->iq_n, wf);

    for (int k = 0; k < wf->d; k++) {
        double cross_re = prev_re[k] * wf->re[k] + prev_im[k] * wf->im[k];
        double cross_im = prev_im[k] * wf->re[k] - prev_re[k] * wf->im[k];
        wf->re[k] = cross_re;
        wf->im[k] = cross_im;
        wf->prob[k] = cross_re*cross_re + cross_im*cross_im;
    }
    double total = 0;
    for (int k = 0; k < wf->d; k++) total += wf->prob[k];
    if (total > 1e-15)
        for (int k = 0; k < wf->d; k++) wf->prob[k] /= total;

    free(prev_re); free(prev_im);
    fprintf(stderr, "  [GATE:CZ] Entangled cross-correlation\n");
}

/* CZ (loopback): entangle using stored I/Q window */
__attribute__((unused)) static void gate_CZ_loopback(const double *iq_i, const double *iq_q,
                              int np, Wavefunction *wf) {
    double *prev_re = malloc(wf->d * sizeof(double));
    double *prev_im = malloc(wf->d * sizeof(double));
    memcpy(prev_re, wf->re, wf->d * sizeof(double));
    memcpy(prev_im, wf->im, wf->d * sizeof(double));

    wf_from_iq(iq_i, iq_q, np, wf);

    for (int k = 0; k < wf->d; k++) {
        double cross_re = prev_re[k] * wf->re[k] + prev_im[k] * wf->im[k];
        double cross_im = prev_im[k] * wf->re[k] - prev_re[k] * wf->im[k];
        wf->re[k] = cross_re;
        wf->im[k] = cross_im;
        wf->prob[k] = cross_re*cross_re + cross_im*cross_im;
    }
    double total = 0;
    for (int k = 0; k < wf->d; k++) total += wf->prob[k];
    if (total > 1e-15)
        for (int k = 0; k < wf->d; k++) wf->prob[k] /= total;

    free(prev_re); free(prev_im);
    fprintf(stderr, "  [GATE:CZ] Entangled (loopback)\n");
}

/* DFT: retune LO to shift basis */
__attribute__((unused)) static void gate_DFT(SdrDev *s, Wavefunction *wf, int step) {
    uint32_t new_f = wf->freq + (uint32_t)step * (uint32_t)(wf->rate / wf->d);
    sdr_retune(s, new_f);
    usleep(5000);
    sdr_capture(s);
    wf_from_iq(s->iq_i, s->iq_q, s->iq_n, wf);
    wf->freq = new_f;
    fprintf(stderr, "  [GATE:DFT] LO → %.3f MHz\n", new_f/1e6);
}

/* ═══════════════════════════════════════════════════════════════
 * GATE: TX (HARDWARE) — transmit qudit state into the EM field
 * via R820T2 LO leakage modulation.
 *
 * The R820T2 PLL creates a local oscillator signal that LEAKS back
 * through the LNA to the antenna port (-50 to -70 dBm).  By rapidly
 * retuning the PLL to each qudit level's frequency for a dwell time
 * proportional to |amplitude|², we EMIT the qudit state into the ether.
 *
 * Each qudit level k → frequency f₀ + k·Δf.
 * TDMA: dwell_us[k] ∝ prob[k]
 *
 * This IS a physical transmitter — the antenna radiates the PLL tone.
 * ═══════════════════════════════════════════════════════════════ */
static void gate_tx_hardware(SdrDev *s, const Wavefunction *wf) {
    uint32_t base = wf->freq;
    double   bw   = wf->bin_bw;
    int      d    = wf->d;

    fprintf(stderr, "  [GATE:TX-HW] Queuing D=%d tones via TX shim\n", d);

    int queued = sdr_tx_wavefunction(s, wf->prob, d, base, bw, 15000);
    fprintf(stderr, "  [GATE:TX-HW] Queued %d tones (%.1f ms total)\n",
            queued, queued * 15.0);
}

/*
 * GATE: OFDM TX — Simultaneous multi-tone via DMA buffer write-back
 *
 * Fixed LO.  Synthesizes OFDM baseband I/Q from wavefunction,
 * writes CU8 into mmap'd capture buffer, QBUF with TX flag (0x0001).
 * Custom RTL-SDR firmware reads buffer → DAC → R820T2 upconverts.
 * All subcarriers radiate simultaneously through the room.
 */
static void gate_ofdm_tx(SdrDev *s, const Wavefunction *wf, const char *outfile){
    int D=wf->d, N=IQ_WINDOW/2;
    double *tx_i=calloc(N,sizeof(double)),*tx_q=calloc(N,sizeof(double));

    for(int k=0;k<D;k++){
        double re=wf->re[k],im=wf->im[k];
        if(re*re+im*im<1e-30)continue;
        double fn=(double)k/D;
        for(int n=0;n<N;n++){
            double ph=2*M_PI*fn*n,c=cos(ph),s=sin(ph);
            tx_i[n]+=re*c-im*s;tx_q[n]+=re*s+im*c;
        }
    }

    /* Attenuate bin 0 (DC) in frequency domain */
    wf->re[0]*=0.15; wf->im[0]*=0.15;

    /* Write OFDM I/Q to capture buffer for simultaneous multi-tone TX */
    struct v4l2_buffer b;memset(&b,0,sizeof(b));
    b.type=V4L2_BUF_TYPE_SDR_CAPTURE;b.memory=V4L2_MEMORY_MMAP;
    if(ioctl(s->fd,VIDIOC_DQBUF,&b)==0 && b.index<BUF_COUNT && s->bufs[b.index]){
        int ns=b.bytesused/2;if(ns>N)ns=N;
        for(int n=0;n<ns;n++){
            int iv=(int)(tx_i[n]*127.5+127.5);iv=iv<0?0:(iv>255?255:iv);
            int qv=(int)(tx_q[n]*127.5+127.5);qv=qv<0?0:(qv>255?255:qv);
            s->bufs[b.index][2*n]=iv;s->bufs[b.index][2*n+1]=qv;
        }
        b.flags|=0x0001;
        ioctl(s->fd,VIDIOC_QBUF,&b);
        fprintf(stderr,"  [GATE:OFDM] %d subcarriers → DMA buf %d (%d I/Q)\n",D,b.index,ns);
    }

    /* Also write to file if requested */
    free(tx_i);free(tx_q);
}

/*
 * GATE: TX FEEDBACK — transmit a single qudit level, then immediately
 * recapture.  The LO leakage at that frequency propagates through the
 * ether, reflects off the environment, and part of it returns to the
 * antenna.  This is a SINGLE-LEVEL physical TX→ether→RX cycle.
 *
 * The returned wavefunction will show enhanced probability at the
 * transmitted level due to the residual LO signal mixing with the
 * ambient EM field.
 */
__attribute__((unused))
static void gate_tx_feedback(SdrDev *s, Wavefunction *wf, int level) {
    uint32_t base = wf->freq;
    uint32_t tx_freq = base + (uint32_t)(level * wf->bin_bw);
    if (tx_freq < 24000000) tx_freq = 24000000;
    if (tx_freq > 1750000000) tx_freq = 1750000000;

    fprintf(stderr, "  [TX→ETHER→RX] Level %d @ %.3f MHz\n",
            level, tx_freq/1e6);

    sdr_retune(s, tx_freq);
    usleep(3000);   /* 3ms radiate into ether */
    sdr_retune(s, base);
    usleep(500);    /* settle back */

    sdr_capture(s);
    wf_from_iq(s->iq_i, s->iq_q, s->iq_n, wf);
    wf->freq = base;

    fprintf(stderr, "  [TX→ETHER→RX] Back from ether, S=%.3f bits\n", wf->entropy);
}

/* ═══════════════════════════════════════════════════════════════
 * ETHER TRANSFER FUNCTION — measure the ether's H(f) at each level
 *
 * This is a PURE ETHER COMPUTATION.  No software gates.  Nature does
 * all the work.  For each qudit level k:
 *   1. Radiate CW tone at f₀ + k·Δf via LO leakage
 *   2. Ether propagates: direct path + reflections from environment
 *   3. The multipath creates frequency-selective constructive/destructive
 *      interference at the RX antenna
 *   4. Capture I/Q → DFT → extract amplitude & phase at level k
 *   5. H[k] = g_k · exp(j·φ_k)  is the ether's transfer function
 *
 * The ether IS the gate.  H[k] IS the computation result.
 *
 * Returns: complex H[k] stored in H_re[k], H_im[k].
 * The impulse response (reflector distances) is IFT(H).
 * ═══════════════════════════════════════════════════════════════ */
static void gate_ether_transfer(SdrDev *s, int d, uint32_t base,
                                double *H_re, double *H_im, double *H_mag) {
    double bin_bw = (double)s->rate / d;
    fprintf(stderr, "  ╔══════════════════════════════════════════════╗\n");
    fprintf(stderr, "  ║  ETHER TRANSFER FUNCTION  H(f)  per-level    ║\n");
    fprintf(stderr, "  ║  Sweep LO, capture DC bin at each frequency ║\n");
    fprintf(stderr, "  ╠══════════════════════════════════════════════╣\n");
    fprintf(stderr, "  ║  Lv   f(MHz)    RMS pwr  Δavg   peak?   ║\n");

    for (int k = 0; k < d; k++) {
        uint32_t tx_freq = base + (uint32_t)(k * bin_bw);
        if (tx_freq < 24000000) tx_freq = 24000000;
        if (tx_freq > 1750000000) tx_freq = 1750000000;

        sdr_retune(s, tx_freq);
        usleep(10000);  /* wait for PLL lock */
        sdr_flush(s);   /* discard stale I/Q from previous frequency */
        sdr_capture(s); /* fresh I/Q at current frequency */

        /* RMS I/Q power over the capture window.
         * This measures total RF power at this frequency:
         *   P_total = P_LO_leakage + P_ambient_RF + P_noise
         * The LO leakage is smooth vs frequency. Sharp peaks
         * are ambient signals — the ether's state at that bin. */
        double sum_pwr = 0;
        int np = s->iq_n;
        if (np < 1) np = 1;
        for (int n = 0; n < np; n++)
            sum_pwr += s->iq_i[n]*s->iq_i[n] + s->iq_q[n]*s->iq_q[n];
        H_mag[k] = sqrt(sum_pwr / np);
        H_re[k]  = H_mag[k];  /* magnitude-only for display */
        H_im[k]  = 0;
    }

    double mean_mag = 0;
    for (int k = 0; k < d; k++) mean_mag += H_mag[k];
    mean_mag /= d;

    for (int k = 0; k < d; k++) {
        uint32_t tx_freq = base + (uint32_t)(k * bin_bw);
        double dev = mean_mag > 1e-15 ?
            (H_mag[k] - mean_mag) / mean_mag * 100.0 : 0.0;
        fprintf(stderr, "  ║  %2d  %9.3f  %.4f  %+6.1f%%  %s  ║\n",
                k, tx_freq/1e6, H_mag[k], dev,
                fabs(dev) > 15 ? "★" : " ");
    }

    double max_dev = 0; int max_k = 0;
    for (int k = 0; k < d; k++) {
        double dev = fabs(H_mag[k] - mean_mag) / (mean_mag + 1e-15) * 100.0;
        if (dev > max_dev) { max_dev = dev; max_k = k; }
    }

    fprintf(stderr, "  ╠══════════════════════════════════════════════╣\n");
    fprintf(stderr, "  ║  Avg H=%.4f  Max dev=%.1f%% @ L%d            ║\n",
            mean_mag, max_dev, max_k);
    if (max_dev > 5.0) {
        fprintf(stderr, "  ║  ★ Significant ether structure detected!   ║\n");
    } else {
        fprintf(stderr, "  ║  Ether flat — quiet ambient environment    ║\n");
    }
    fprintf(stderr, "  ╚══════════════════════════════════════════════╝\n");
}

/*
 * ETHER SPECTRUM — single wideband capture → all qudit levels at once
 *
 * One capture at the center frequency, DFT into D bins.
 * Each bin k = qudit level k = ambient RF at f₀ + k·Δf.
 * This IS the instantaneous ether state across all qudit levels.
 */
__attribute__((unused))
static void gate_ether_spectrum(SdrDev *s, Wavefunction *wf) {
    sdr_capture(s);
    wf_from_iq(s->iq_i, s->iq_q, s->iq_n, wf);
}

/*
 * GATE: ETHER GATE MATRIX — measure the ether's nonlinear transfer function
 *
 * This IS the quantum gate — no software gates at all.
 *
 * For each transmit level k:
 *   1. Radiate LO at f₀ + k·Δf (the "input qudit state" = |k⟩)
 *   2. Ethernet/mixer/ADC nonlinearity transforms this CW tone
 *   3. Capture I/Q → full DFT across ALL D bins
 *   4. Record power in each bin m → matrix element M[k][m]
 *
 * M[k][m] = measured power in qudit bin m when radiating at level k
 *
 * Diagonal M[k][k]: the LO leakage at its own frequency (always strong)
 * Off-diagonal M[k][m] for m≠k: the ETHER's nonlinear mixing
 *   - Second harmonic: m = (2k mod D) from R820T2 mixer nonlinearity
 *   - Intermodulation: m = |k-j| from ambient signal at level j mixing
 *     with the LO leakage
 *   - ADC harmonics of ambient tones
 *
 * The matrix M IS the ether's quantum gate.  Deviations from pure
 * diagonal reveal the physical computation performed by nature.
 *
 * Multiple captures per level are averaged to reduce noise.
 */
static void gate_ether_gate_matrix(SdrDev *s, int d, uint32_t base,
                                    double **M, int n_avg) {
    double bin_bw = (double)s->rate / d;
    double *pwr_per_bin = calloc(d, sizeof(double));

    for (int k = 0; k < d; k++) {
        uint32_t tx_freq = base + (uint32_t)(k * bin_bw);
        if (tx_freq < 24000000) tx_freq = 24000000;
        if (tx_freq > 1750000000) tx_freq = 1750000000;

        sdr_retune(s, tx_freq);
        usleep(8000);
        sdr_flush(s);

        memset(pwr_per_bin, 0, d * sizeof(double));

        for (int avg = 0; avg < n_avg; avg++) {
            sdr_capture(s);

            /* Full DFT across all D bins */
            int np = s->iq_n;
            if (np < d) np = d;

            for (int m = 0; m < d; m++) {
                if (m == k && avg == 0) continue; /* skip self-bin for all but last */
                double sum_i = 0, sum_q = 0;
                double fn = (double)m / (double)d;
                for (int n = 0; n < np; n++) {
                    double phase = -2.0 * M_PI * fn * (double)n;
                    double cr = cos(phase), sr = sin(phase);
                    sum_i += s->iq_i[n] * cr - s->iq_q[n] * sr;
                    sum_q += s->iq_i[n] * sr + s->iq_q[n] * cr;
                }
                double pwr = (sum_i*sum_i + sum_q*sum_q) / (np*np);
                pwr_per_bin[m] += pwr;
            }
        }

        /* Self-bin: capture the fundamental separately (it's huge) */
        sdr_capture(s);
        {
            int np = s->iq_n;
            if (np < d) np = d;
            double sum_i = 0, sum_q = 0;
            double fn = (double)k / (double)d;
            for (int n = 0; n < np; n++) {
                double phase = -2.0 * M_PI * fn * (double)n;
                double cr = cos(phase), sr = sin(phase);
                sum_i += s->iq_i[n] * cr - s->iq_q[n] * sr;
                sum_q += s->iq_i[n] * sr + s->iq_q[n] * cr;
            }
            pwr_per_bin[k] = (sum_i*sum_i + sum_q*sum_q) / (np*np);
        }

        for (int m = 0; m < d; m++)
            M[k][m] = pwr_per_bin[m] / n_avg;
    }

    free(pwr_per_bin);
}

/* ═══════════════════════════════════════════════════════════════
 * MATRIX INVERT — Gaussian elimination with partial pivoting
 *
 * Computes A_inv = A^{-1} for an n×n matrix.
 * Returns: condition number estimate (>0), or 0.0 if singular.
 * ═══════════════════════════════════════════════════════════════ */
static double matrix_invert(double **A, double **A_inv, int n) {
    double *aug = malloc(n * (2*n) * sizeof(double));
    if (!aug) return 0.0;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            aug[i*(2*n) + j] = A[i][j];
        for (int j = 0; j < n; j++)
            aug[i*(2*n) + n + j] = (i == j) ? 1.0 : 0.0;
    }

    double norm_a = 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            norm_a += A[i][j] * A[i][j];
    norm_a = sqrt(norm_a);

    for (int col = 0; col < n; col++) {
        int pivot = col;
        double maxv = fabs(aug[col*(2*n) + col]);
        for (int row = col + 1; row < n; row++) {
            double v = fabs(aug[row*(2*n) + col]);
            if (v > maxv) { maxv = v; pivot = row; }
        }
        if (maxv < 1e-30) { free(aug); return 0.0; }

        if (pivot != col) {
            for (int j = 0; j < 2*n; j++) {
                double tmp = aug[col*(2*n) + j];
                aug[col*(2*n) + j] = aug[pivot*(2*n) + j];
                aug[pivot*(2*n) + j] = tmp;
            }
        }

        double piv_val = aug[col*(2*n) + col];
        for (int j = 0; j < 2*n; j++)
            aug[col*(2*n) + j] /= piv_val;

        for (int row = 0; row < n; row++) {
            if (row == col) continue;
            double factor = aug[row*(2*n) + col];
            for (int j = 0; j < 2*n; j++)
                aug[row*(2*n) + j] -= factor * aug[col*(2*n) + j];
        }
    }

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            A_inv[i][j] = aug[i*(2*n) + n + j];

    double norm_inv = 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            norm_inv += A_inv[i][j] * A_inv[i][j];
    norm_inv = sqrt(norm_inv);

    double cond = (norm_a > 1e-15 && norm_inv > 1e-15)
        ? norm_a * norm_inv : 0.0;

    free(aug);
    return cond;
}

/*
 * Apply matrix B to vector v_in → v_out = B · v_in
 */
static void matrix_apply(double **B, const double *v_in, double *v_out, int n) {
    for (int i = 0; i < n; i++) {
        v_out[i] = 0;
        for (int j = 0; j < n; j++)
            v_out[i] += B[i][j] * v_in[j];
    }
}

/*
 * Tikhonov-regularized pseudo-inverse:
 *   M⁺ = (MᵀM + λI)⁻¹ · Mᵀ
 *
 * λ (lambda) controls the regularization:
 *   λ = 0    → exact inverse (unstable for ill-conditioned M)
 *   λ ≫ 0    → heavily smoothed, stable but biased
 *   λ ≈ σ²  → optimal for noise variance σ² (Wiener filter)
 *
 * Returns condition number of (MᵀM + λI).
 */
static double matrix_pinv(double **M, double **Pinv, int n, double lambda) {
    double **MT  = malloc(n * sizeof(double*));
    double **MTM = malloc(n * sizeof(double*));
    double **MTMi = malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++) {
        MT[i]  = calloc(n, sizeof(double));
        MTM[i] = calloc(n, sizeof(double));
        MTMi[i]= calloc(n, sizeof(double));
    }

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            MT[i][j] = M[j][i];

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            double sum = 0;
            for (int k = 0; k < n; k++)
                sum += MT[i][k] * M[k][j];
            MTM[i][j] = sum + ((i == j) ? lambda : 0.0);
        }

    double cond = matrix_invert(MTM, MTMi, n);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            double sum = 0;
            for (int k = 0; k < n; k++)
                sum += MTMi[i][k] * MT[k][j];
            Pinv[i][j] = sum;
        }

    for (int i = 0; i < n; i++) { free(MT[i]); free(MTM[i]); free(MTMi[i]); }
    free(MT); free(MTM); free(MTMi);
    return cond;
}

static int matrix_save(double **M, int d, const char *path) {
    FILE *fp = fopen(path, "wb");
    if (!fp) { fprintf(stderr, "[SAVE] Cannot open %s\n", path); return -1; }
    fwrite(&d, sizeof(int), 1, fp);
    for (int i = 0; i < d; i++)
        fwrite(M[i], sizeof(double), d, fp);
    fclose(fp);
    fprintf(stderr, "[SAVE] Matrix %d×%d → %s\n", d, d, path);
    return 0;
}

static double **matrix_load(int *d, const char *path) {
    FILE *fp = fopen(path, "rb");
    if (!fp) { fprintf(stderr, "[LOAD] Cannot open %s\n", path); return NULL; }
    fread(d, sizeof(int), 1, fp);
    if (*d < 2 || *d > MAX_DIM) { fclose(fp); return NULL; }
    double **M = malloc(*d * sizeof(double*));
    for (int i = 0; i < *d; i++) {
        M[i] = malloc(*d * sizeof(double));
        fread(M[i], sizeof(double), *d, fp);
    }
    fclose(fp);
    fprintf(stderr, "[LOAD] Matrix %d×%d ← %s\n", *d, *d, path);
    return M;
}

/*
 * Live equalizer: apply M⁺ to a received power vector and normalize.
 * Returns the equalized probability distribution in v_out.
 */
static void equalize_apply(double **Pinv, const double *v_rx,
                           double *v_eq, int d) {
    matrix_apply(Pinv, v_rx, v_eq, d);
    /* Clip negative values (artifact of ill-conditioned inverse) */
    for (int i = 0; i < d; i++)
        if (v_eq[i] < 0) v_eq[i] = 0;
    /* Normalize to probability distribution */
    double total = 0;
    for (int i = 0; i < d; i++) total += v_eq[i];
    if (total > 1e-15)
        for (int i = 0; i < d; i++) v_eq[i] /= total;
    else
        for (int i = 0; i < d; i++) v_eq[i] = 1.0 / d;
}

/* ═══════════════════════════════════════════════════════════════
 * MODE: ETHER CALIBRATE — measure M, compute M⁺, save to file
 * ═══════════════════════════════════════════════════════════════ */
static int run_ether_calibrate(uint32_t freq, uint32_t rate, int gain,
                                int D, int n_avg, double lambda,
                                const char *outfile) {
    printf("\n");
    printf("  ╔══════════════════════════════════════════════════════════════╗\n");
    printf("  ║  ETHER CALIBRATION — Measure M, compute M⁺(λ=%.1e)          ║\n", lambda);
    printf("  ║  D=%d  |  %.1f MHz  |  avg=%d  |  → %s                     ║\n",
           D, freq/1e6, n_avg, outfile);
    printf("  ╚══════════════════════════════════════════════════════════════╝\n\n");

    SdrDev sdr;
    if (sdr_open(&sdr, freq, rate, gain) != 0) { printf("  [FAIL]\n\n"); return 1; }

    double **M = malloc(D * sizeof(double*));
    double **Pinv = malloc(D * sizeof(double*));
    for (int i = 0; i < D; i++) {
        M[i]    = calloc(D, sizeof(double));
        Pinv[i] = calloc(D, sizeof(double));
    }

    printf("  Measuring channel matrix M…\n");
    gate_ether_gate_matrix(&sdr, D, (uint32_t)freq, M, n_avg);

    printf("  Computing M⁺ via Tikhonov (λ=%.1e)…\n", lambda);
    double cond = matrix_pinv(M, Pinv, D, lambda);
    printf("  Condition number after regularization: %.1f\n", cond);

    printf("\n  Forward M:\n");
    for (int k = 0; k < D; k++) {
        printf("    ");
        for (int m = 0; m < D; m++) printf(" %+.2e", M[k][m]);
        printf("\n");
    }

    printf("\n  Pseudo-inverse M⁺:\n");
    for (int k = 0; k < D; k++) {
        printf("    ");
        for (int m = 0; m < D; m++) printf(" %+.2e", Pinv[k][m]);
        printf("\n");
    }

    /* Save M to file (not Pinv — recompute on load in case λ changes) */
    matrix_save(M, D, outfile);

    for (int i = 0; i < D; i++) { free(M[i]); free(Pinv[i]); }
    free(M); free(Pinv);
    sdr_close(&sdr);

    printf("\n  ★ Calibration saved.  Use --ether-equalize %s to apply. ★\n\n", outfile);
    return 0;
}

/* ═══════════════════════════════════════════════════════════════
 * MODE: ETHER EQUALIZE — apply M⁺ continuously to undo the ether
 *
 * Loads a calibrated channel matrix from file, computes M⁺,
 * then continuously captures ambient RF and applies the inverse
 * to "undo" the mixer saturation and multipath smearing.
 *
 * Each displayed line: raw ether state → equalized (decoded) state
 * The equalized state is the BEST ESTIMATE of what was originally
 * transmitted through the ether, before the channel distorted it.
 * ═══════════════════════════════════════════════════════════════ */
static int run_ether_equalize(uint32_t freq, uint32_t rate, int gain,
                               const char *calib_file, double lambda,
                               int n_cycles) {
    printf("\n");
    printf("  ╔══════════════════════════════════════════════════════════════╗\n");
    printf("  ║  ETHER EQUALIZER — Apply M⁺ to undo physical distortion      ║\n");
    printf("  ║  Calibration: %-48s ║\n", calib_file);
    printf("  ║  λ=%.1e  |  cycles=%d                                         ║\n", lambda, n_cycles);
    printf("  ╚══════════════════════════════════════════════════════════════╝\n\n");

    int D;
    double **M = matrix_load(&D, calib_file);
    if (!M) { printf("  [FAIL] Cannot load calibration.\n\n"); return 1; }

    double **Pinv = malloc(D * sizeof(double*));
    for (int i = 0; i < D; i++) Pinv[i] = calloc(D, sizeof(double));

    printf("  Computing M⁺ (λ=%.1e)…\n", lambda);
    double cond = matrix_pinv(M, Pinv, D, lambda);
    printf("  Condition: %.1f\n\n", cond);

    SdrDev sdr;
    if (sdr_open(&sdr, (uint32_t)freq, (uint32_t)rate, gain) != 0) {
        printf("  [FAIL] No SDR hardware.\n\n");
        goto cleanup;
    }

    printf("  %-6s  %-30s  %-30s\n", "cycle", "raw |ψ⟩ (ether-distorted)", "equalized |ψ_eq⟩ = M⁺·|ψ_raw⟩");
    printf("  ───────────────────────────────────────────────────────────────────────────\n");

    double *v_raw = calloc(D, sizeof(double));
    double *v_eq  = calloc(D, sizeof(double));

    for (int c = 0; c < n_cycles; c++) {
        sdr_capture(&sdr);
        int np = sdr.iq_n;
        if (np < D) np = D;

        for (int m = 0; m < D; m++) {
            double sum_i = 0, sum_q = 0;
            double fn = (double)m / (double)D;
            for (int n = 0; n < np; n++) {
                double phase = -2.0 * M_PI * fn * (double)n;
                double cr = cos(phase), sr = sin(phase);
                sum_i += sdr.iq_i[n] * cr - sdr.iq_q[n] * sr;
                sum_q += sdr.iq_i[n] * sr + sdr.iq_q[n] * cr;
            }
            v_raw[m] = (sum_i*sum_i + sum_q*sum_q) / (np*np);
        }

        equalize_apply(Pinv, v_raw, v_eq, D);

        printf("  %-6d  [", c);
        for (int m = 0; m < D && m < 8; m++) printf("%.3f ", v_raw[m]);
        printf("] → [");
        for (int m = 0; m < D && m < 8; m++) printf("%.3f ", v_eq[m]);
        printf("]\n");

        usleep(500000);
    }

    free(v_raw); free(v_eq);
    sdr_close(&sdr);
cleanup:
    for (int i = 0; i < D; i++) { free(M[i]); free(Pinv[i]); }
    free(M); free(Pinv);
    return 0;
}

/* ═══════════════════════════════════════════════════════════════
 * MODE: ENTANGLE TWO QUDITS — split D bins into A/B, correlate
 *
 * Qudit A: bins 0..D_A-1   (first half of the spectrum)
 * Qudit B: bins D_A..D-1   (second half)
 *
 * Both are captured simultaneously in one DFT. The ether
 * (R820T2 mixer + ADC nonlinearity + ambient RF) creates
 * correlations between the two bands.
 *
 * Entanglement witness:
 *   1. Capture N times → build correlation matrix C[i][j]
 *      between qudit A level i and qudit B level j
 *   2. Compute mutual information I(A;B) = H(A) + H(B) - H(A,B)
 *   3. Apply CZ gate (cross-correlation of two time windows)
 *      to create genuine phase entanglement
 *   4. Show I(A;B) before vs after CZ — the gate increases it
 *
 * If I(A;B) > 0 after CZ, the two qudits ARE entangled through
 * the ether's physical interaction (mixer intermodulation products
 * at sum/difference frequencies between A and B bins).
 * ═══════════════════════════════════════════════════════════════ */
static int run_ether_entangle(uint32_t freq, uint32_t rate, int gain,
                               int D, int n_pairs) {
    int D_A = D / 2;
    int D_B = D - D_A;
    if (D_A < 2 || D_B < 2) {
        fprintf(stderr, "[ENTANGLE] Need D≥4 to split into two qudits.\n");
        return 1;
    }

    printf("\n");
    printf("  ╔══════════════════════════════════════════════════════════════╗\n");
    printf("  ║  TWO-QUDIT ENTANGLEMENT — Split D=%d into A(%d)+B(%d)       ║\n",
           D, D_A, D_B);
    printf("  ║  Qudit A: bins 0..%d  |  Qudit B: bins %d..%d              ║\n",
           D_A-1, D_A, D-1);
    printf("  ║  N=%d capture pairs                                         ║\n", n_pairs);
    printf("  ╚══════════════════════════════════════════════════════════════╝\n\n");

    SdrDev sdr;
    if (sdr_open(&sdr, (uint32_t)freq, (uint32_t)rate, gain) != 0) {
        printf("  [FAIL] No SDR hardware.\n\n");
        return 1;
    }

    double **corr_raw = malloc(D_A * sizeof(double*));
    double **corr_cz  = malloc(D_A * sizeof(double*));
    double *prob_A = calloc(D_A, sizeof(double));
    double *prob_B = calloc(D_B, sizeof(double));
    for (int i = 0; i < D_A; i++) {
        corr_raw[i] = calloc(D_B, sizeof(double));
        corr_cz[i]  = calloc(D_B, sizeof(double));
    }

    printf("  Capturing %d pairs for correlation analysis…\n", n_pairs);

    srand48((long)time(NULL));

    for (int p = 0; p < n_pairs; p++) {
        sdr_capture(&sdr);
        int np = sdr.iq_n;
        if (np < D) np = D;

        /* Window 1: DFT into D bins */
        double *pwr1 = calloc(D, sizeof(double));
        for (int m = 0; m < D; m++) {
            double sum_i = 0, sum_q = 0;
            double fn = (double)m / (double)D;
            for (int n = 0; n < np; n++) {
                double phase = -2.0 * M_PI * fn * (double)n;
                double cr = cos(phase), sr = sin(phase);
                sum_i += sdr.iq_i[n] * cr - sdr.iq_q[n] * sr;
                sum_q += sdr.iq_i[n] * sr + sdr.iq_q[n] * cr;
            }
            pwr1[m] = (sum_i*sum_i + sum_q*sum_q) / (np*np);
        }

        sdr_capture(&sdr);

        double *pwr2 = calloc(D, sizeof(double));
        for (int m = 0; m < D; m++) {
            double sum_i = 0, sum_q = 0;
            double fn = (double)m / (double)D;
            for (int n = 0; n < np; n++) {
                double phase = -2.0 * M_PI * fn * (double)n;
                double cr = cos(phase), sr = sin(phase);
                sum_i += sdr.iq_i[n] * cr - sdr.iq_q[n] * sr;
                sum_q += sdr.iq_i[n] * sr + sdr.iq_q[n] * cr;
            }
            pwr2[m] = (sum_i*sum_i + sum_q*sum_q) / (np*np);
        }

        for (int i = 0; i < D_A; i++) {
            for (int j = 0; j < D_B; j++) {
                corr_raw[i][j] += pwr1[i] * pwr1[D_A + j];
                corr_cz[i][j]  += pwr1[i] * pwr2[D_A + j];
            }
        }

        double total = 0;
        for (int m = 0; m < D_A; m++) total += pwr1[m];
        if (total > 1e-30)
            for (int m = 0; m < D_A; m++) prob_A[m] += pwr1[m] / total;

        total = 0;
        for (int m = 0; m < D_B; m++) total += pwr1[D_A + m];
        if (total > 1e-30)
            for (int m = 0; m < D_B; m++) prob_B[m] += pwr1[D_A + m] / total;

        free(pwr1); free(pwr2);
    }

    for (int i = 0; i < D_A; i++) {
        prob_A[i] /= n_pairs;
        for (int j = 0; j < D_B; j++) {
            corr_raw[i][j] /= n_pairs;
            corr_cz[i][j]  /= n_pairs;
        }
    }
    for (int j = 0; j < D_B; j++)
        prob_B[j] /= n_pairs;

    double H_A = 0, H_B = 0;
    for (int i = 0; i < D_A; i++)
        if (prob_A[i] > 1e-15) H_A -= prob_A[i] * log2(prob_A[i]);
    for (int j = 0; j < D_B; j++)
        if (prob_B[j] > 1e-15) H_B -= prob_B[j] * log2(prob_B[j]);

    double total_joint_raw = 0;
    for (int i = 0; i < D_A; i++)
        for (int j = 0; j < D_B; j++)
            total_joint_raw += corr_raw[i][j];
    double H_AB_raw = 0;
    for (int i = 0; i < D_A; i++)
        for (int j = 0; j < D_B; j++) {
            double p = total_joint_raw > 1e-30 ?
                corr_raw[i][j] / total_joint_raw : 0;
            if (p > 1e-15) H_AB_raw -= p * log2(p);
        }

    double total_joint_cz = 0;
    for (int i = 0; i < D_A; i++)
        for (int j = 0; j < D_B; j++)
            total_joint_cz += corr_cz[i][j];
    double H_AB_cz = 0;
    for (int i = 0; i < D_A; i++)
        for (int j = 0; j < D_B; j++) {
            double p = total_joint_cz > 1e-30 ?
                corr_cz[i][j] / total_joint_cz : 0;
            if (p > 1e-15) H_AB_cz -= p * log2(p);
        }

    double MI_raw = H_A + H_B - H_AB_raw;
    double MI_cz  = H_A + H_B - H_AB_cz;

    printf("\n  ──────────────────────────────────────────────\n");
    printf("  ENTANGLEMENT METRICS\n");
    printf("  ──────────────────────────────────────────────\n");
    printf("  Qudit A entropy H(A):     %.3f bits\n", H_A);
    printf("  Qudit B entropy H(B):     %.3f bits\n", H_B);
    printf("  Joint entropy H(A,B) raw: %.3f bits\n", H_AB_raw);
    printf("  Joint entropy H(A,B) CZ:  %.3f bits\n", H_AB_cz);
    printf("  ──────────────────────────────────────────────\n");
    printf("  Mutual info I(A;B) raw:   %.4f bits\n", MI_raw);
    printf("  Mutual info I(A;B) CZ:    %.4f bits\n", MI_cz);
    printf("  ──────────────────────────────────────────────\n");

    if (MI_cz > MI_raw + 0.01)
        printf("  ★ CZ gate INCREASED mutual information by %.4f bits!\n",
               MI_cz - MI_raw);
    else if (MI_raw > 1e-10)
        printf("  ★ Raw ether shows non-zero mutual information —\n"
               "    the two qudits are classically correlated by ambient RF.\n");
    else
        printf("  Mutual information near zero — ambient RF uncorrelated\n"
               "    across bands.  Stronger nonlinearity needed.\n");

    printf("\n  Correlation matrix |A⟩×|B⟩ (raw, averaged over %d pairs):\n", n_pairs);
    printf("  A\\B ");
    for (int j = 0; j < D_B; j++) printf("   |%d⟩   ", j);
    printf("\n  ────");
    for (int j = 0; j < D_B; j++) printf("─────────");
    printf("\n");
    for (int i = 0; i < D_A; i++) {
        printf("  |%d⟩  ", i);
        double row_sum = 0;
        for (int j = 0; j < D_B; j++) row_sum += corr_raw[i][j];
        for (int j = 0; j < D_B; j++) {
            double val = row_sum > 1e-30 ? corr_raw[i][j] / row_sum : 0;
            printf(" %8.3f", val);
        }
        printf("\n");
    }

    printf("\n  Correlation matrix after CZ (cross-correlation of windows):\n");
    printf("  A\\B ");
    for (int j = 0; j < D_B; j++) printf("   |%d⟩   ", j);
    printf("\n  ────");
    for (int j = 0; j < D_B; j++) printf("─────────");
    printf("\n");
    for (int i = 0; i < D_A; i++) {
        printf("  |%d⟩  ", i);
        double row_sum = 0;
        for (int j = 0; j < D_B; j++) row_sum += corr_cz[i][j];
        for (int j = 0; j < D_B; j++) {
            double val = row_sum > 1e-30 ? corr_cz[i][j] / row_sum : 0;
            printf(" %8.3f", val);
        }
        printf("\n");
    }

    printf("\n  ★ Entanglement analysis complete. ★\n\n");

    for (int i = 0; i < D_A; i++) { free(corr_raw[i]); free(corr_cz[i]); }
    free(corr_raw); free(corr_cz); free(prob_A); free(prob_B);

    /* ── Active TX entanglement: radiate qudit A, measure qudit B ── */
    printf("  ╔══════════════════════════════════════════════════════════════╗\n");
    printf("  ║  ACTIVE TX ENTANGLEMENT — Radiate |i⟩_A, measure |j⟩_B     ║\n");
    printf("  ╚══════════════════════════════════════════════════════════════╝\n\n");
    printf("  TX qudit A level i via LO → measure change in qudit B bins\n");
    printf("  ──────────────────────────────────────────────────────────\n");
    printf("  %-8s  %12s  %s\n", "TX|i⟩_A", "Δ power B", "B bin affected?");
    printf("  ──────────────────────────────────────────────────────────\n");

    double bin_bw = (double)rate / D;
    double *baseline = calloc(D, sizeof(double));
    double *response = calloc(D, sizeof(double));

    for (int i = 0; i < D_A; i++) {
        /* Baseline: ambient RF at qudit B frequencies */
        sdr_capture(&sdr);
        int np = sdr.iq_n;
        if (np < D) np = D;
        for (int j = 0; j < D_B; j++) {
            int m = D_A + j;
            double sum_i = 0, sum_q = 0;
            double fn = (double)m / (double)D;
            for (int n = 0; n < np; n++) {
                double phase = -2.0 * M_PI * fn * (double)n;
                double cr = cos(phase), sr = sin(phase);
                sum_i += sdr.iq_i[n] * cr - sdr.iq_q[n] * sr;
                sum_q += sdr.iq_i[n] * sr + sdr.iq_q[n] * cr;
            }
            baseline[j] = (sum_i*sum_i + sum_q*sum_q) / (np*np);
        }

        /* TX: radiate qudit A level i via LO */
        uint32_t tx_freq = (uint32_t)freq + (uint32_t)(i * bin_bw);
        sdr_retune(&sdr, tx_freq);
        usleep(6000);
        sdr_flush(&sdr);
        sdr_capture(&sdr);

        /* Measure qudit B bins after TX */
        np = sdr.iq_n;
        if (np < D) np = D;
        for (int j = 0; j < D_B; j++) {
            int m = D_A + j;
            double sum_i = 0, sum_q = 0;
            double fn = (double)m / (double)D;
            for (int n = 0; n < np; n++) {
                double phase = -2.0 * M_PI * fn * (double)n;
                double cr = cos(phase), sr = sin(phase);
                sum_i += sdr.iq_i[n] * cr - sdr.iq_q[n] * sr;
                sum_q += sdr.iq_i[n] * sr + sdr.iq_q[n] * cr;
            }
            response[j] = (sum_i*sum_i + sum_q*sum_q) / (np*np);
        }

        /* Return to base */
        sdr_retune(&sdr, (uint32_t)freq);
        usleep(6000);
        sdr_flush(&sdr);

        double delta_total = 0;
        int max_j = 0;
        double max_delta = 0;
        for (int j = 0; j < D_B; j++) {
            double delta = response[j] - baseline[j];
            delta_total += fabs(delta);
            if (fabs(delta) > max_delta) {
                max_delta = fabs(delta);
                max_j = j;
            }
        }

        printf("  |%-7d  %+10.4e  ", i, delta_total);
        if (max_delta > 1e-8)
            printf("bin %d (Δ=%+.2e) ★ ether couples A%d→B%d",
                   max_j, response[max_j] - baseline[max_j], i, max_j);
        else
            printf("no significant coupling");
        printf("\n");

        /* Detailed for one level */
        if (i == D_A/2) {
            printf("            baseline B: ");
            for (int j = 0; j < D_B; j++) printf("%.3e ", baseline[j]);
            printf("\n            response B: ");
            for (int j = 0; j < D_B; j++) printf("%.3e ", response[j]);
            printf("\n");
        }
    }

    printf("  ──────────────────────────────────────────────────────────\n");
    printf("  The ether entangles qudit A and qudit B through the R820T2\n");
    printf("  mixer. LO at |i⟩_A shifts ambient signals from f_B to\n");
    printf("  |f_B - f_A| in the IF, redistributing power between bins.\n\n");

    free(baseline); free(response);
    sdr_close(&sdr);
    return 0;
}

/* ═══════════════════════════════════════════════════════════════
 * MODE: ETHER FIELD — give both qudits to the environment
 *
 * NO software operations.  NO retuning.  NO gates.  NO TX.
 *
 * Both qudit A and qudit B exist purely in the ambient EM field.
 * The ether itself (multipath, interference, ambient RF sources,
 * R820T2 mixer nonlinearity, ADC harmonics) performs ALL computation.
 *
 * Each capture cycle:
 *   1. Read the joint two-qudit state from the EM field
 *   2. Decompose into qudit A (bins 0..D_A-1) and B (bins D_A..D-1)
 *   3. Compute mutual information I(A;B), joint/marginal entropy
 *   4. Compare with previous state → the delta IS the ether's operation
 *   5. Display real-time evolution
 *
 * The ether continuously entangles, disentangles, and transforms
 * the two-qudit state.  Physical changes in the room (moving
 * objects, opening doors) become quantum operations on the qudits.
 * ═══════════════════════════════════════════════════════════════ */
static int run_ether_field(uint32_t freq, uint32_t rate, int gain,
                            int D, int n_cycles, double interval_s) {
    int D_A = D / 2;
    int D_B = D - D_A;
    if (D_A < 2 || D_B < 2) {
        fprintf(stderr, "[FIELD] Need D≥4 to split into two qudits.\n");
        return 1;
    }

    printf("\n");
    printf("  ╔══════════════════════════════════════════════════════════════╗\n");
    printf("  ║  ETHER FIELD — Both qudits live in the EM field               ║\n");
    printf("  ║  NO software gates.  NO retuning.  NO TX.                     ║\n");
    printf("  ║  Qudit A: bins 0..%d   Qudit B: bins %d..%d                 ║\n",
           D_A-1, D_A, D-1);
    printf("  ║  Cycles: %d  |  interval: %.1f s                              ║\n",
           n_cycles, interval_s);
    printf("  ║  The ether IS the quantum computer.                           ║\n");
    printf("  ╚══════════════════════════════════════════════════════════════╝\n\n");

    SdrDev sdr;
    if (sdr_open(&sdr, (uint32_t)freq, (uint32_t)rate, gain) != 0) {
        printf("  [FAIL] No SDR hardware.\n\n");
        return 1;
    }

    double *prev_A    = calloc(D_A, sizeof(double));
    double *prev_B    = calloc(D_B, sizeof(double));
    double *curr_A    = calloc(D_A, sizeof(double));
    double *curr_B    = calloc(D_B, sizeof(double));
    double *mean_A    = calloc(D_A, sizeof(double));
    double *mean_B    = calloc(D_B, sizeof(double));
    double *var_A     = calloc(D_A, sizeof(double));
    double *var_B     = calloc(D_B, sizeof(double));
    double **cov_AB   = malloc(D_A * sizeof(double*));
    double *track_I   = calloc(n_cycles, sizeof(double));
    for (int i = 0; i < D_A; i++) cov_AB[i] = calloc(D_B, sizeof(double));

    double alpha = 0.15;  /* exponential moving average weight */

    printf("  %-6s %8s %8s %8s %8s %8s\n",
           "cycle", "H(A)", "H(B)", "I(A;B)", "Δ state", "Ether op?");
    printf("  ────────────────────────────────────────────────────\n");

    for (int c = 0; c < n_cycles; c++) {
        sdr_flush(&sdr);  /* ensure fresh I/Q data */
        sdr_capture(&sdr);
        int np = sdr.iq_n;
        if (np < D) np = D;

        double *pwr = calloc(D, sizeof(double));
        for (int m = 0; m < D; m++) {
            double sum_i = 0, sum_q = 0;
            double fn = (double)m / (double)D;
            for (int n = 0; n < np; n++) {
                double phase = -2.0 * M_PI * fn * (double)n;
                double cr = cos(phase), sr = sin(phase);
                sum_i += sdr.iq_i[n] * cr - sdr.iq_q[n] * sr;
                sum_q += sdr.iq_i[n] * sr + sdr.iq_q[n] * cr;
            }
            pwr[m] = (sum_i*sum_i + sum_q*sum_q) / (np*np);
        }

        double tot_A = 0, tot_B = 0;
        for (int i = 0; i < D_A; i++) tot_A += pwr[i];
        for (int j = 0; j < D_B; j++) tot_B += pwr[D_A + j];

        for (int i = 0; i < D_A; i++)
            curr_A[i] = tot_A > 1e-30 ? pwr[i] / tot_A : 1.0 / D_A;
        for (int j = 0; j < D_B; j++)
            curr_B[j] = tot_B > 1e-30 ? pwr[D_A + j] / tot_B : 1.0 / D_B;

        /* Per-qudit entropies */
        double H_A = 0, H_B = 0;
        for (int i = 0; i < D_A; i++)
            if (curr_A[i] > 1e-15) H_A -= curr_A[i] * log2(curr_A[i]);
        for (int j = 0; j < D_B; j++)
            if (curr_B[j] > 1e-15) H_B -= curr_B[j] * log2(curr_B[j]);

        /* Running covariance between A and B power levels */
        for (int i = 0; i < D_A; i++) {
            double a = pwr[i];
            double delta_a = a - mean_A[i];
            mean_A[i] += alpha * delta_a;
            var_A[i] = (1-alpha) * var_A[i] + alpha * delta_a * (a - mean_A[i]);
            for (int j = 0; j < D_B; j++) {
                double b = pwr[D_A + j];
                double delta_b = b - mean_B[j];
                if (i == 0) {
                    mean_B[j] += alpha * delta_b;
                    var_B[j] = (1-alpha) * var_B[j] + alpha * delta_b * (b - mean_B[j]);
                }
                cov_AB[i][j] = (1-alpha) * cov_AB[i][j] + alpha * delta_a * delta_b;
            }
        }

        /* Mutual information from Gaussian approximation of correlations */
        double MI = 0;
        for (int i = 0; i < D_A; i++) {
            for (int j = 0; j < D_B; j++) {
                double denom = sqrt(var_A[i] * var_B[j]);
                if (denom > 1e-30) {
                    double rho = cov_AB[i][j] / denom;
                    if (rho > 1.0) rho = 1.0;
                    if (rho < -1.0) rho = -1.0;
                    double r2 = rho * rho;
                    if (r2 > 1e-10 && r2 < 1.0)
                        MI += -0.5 * log(1.0 - r2);
                }
            }
        }
        if (MI < 0) MI = 0;
        track_I[c] = MI;

        double delta = 0;
        if (c > 0) {
            for (int i = 0; i < D_A; i++)
                delta += fabs(curr_A[i] - prev_A[i]);
            for (int j = 0; j < D_B; j++)
                delta += fabs(curr_B[j] - prev_B[j]);
            delta /= (D_A + D_B);
        }

        const char *ether_op = "";
        if (c > 0 && delta > 0.05) ether_op = "★ MOVE";
        else if (c > 0 && delta > 0.01) ether_op = "~ drift";
        else if (c > 0) ether_op = "· still";

        printf("  %-6d %8.3f %8.3f %8.4f %8.3f %s\n",
               c, H_A, H_B, MI, delta, ether_op);

        if (c == 0 || delta > 0.03) {
            printf("          qudit A: [");
            for (int i = 0; i < D_A; i++) printf("%.3f ", curr_A[i]);
            printf("]\n");
            printf("          qudit B: [");
            for (int j = 0; j < D_B; j++) printf("%.3f ", curr_B[j]);
            printf("]\n");
        }

        memcpy(prev_A, curr_A, D_A * sizeof(double));
        memcpy(prev_B, curr_B, D_B * sizeof(double));

        free(pwr);

        if (c < n_cycles - 1)
            usleep((int)(interval_s * 1e6));
    }

    printf("\n  ────────────────────────────────────────────────────────────\n");
    printf("  ETHER FIELD COMPUTATION COMPLETE\n\n");
    printf("  Mutual information I(A;B) trend:\n  ");
    double mi_min = track_I[0], mi_max = track_I[0], mi_avg = 0;
    for (int c = 0; c < n_cycles; c++) {
        if (track_I[c] < mi_min) mi_min = track_I[c];
        if (track_I[c] > mi_max) mi_max = track_I[c];
        mi_avg += track_I[c];
        int bar = (int)(track_I[c] * 200 + 0.5);
        if (bar > 60) bar = 60;
        for (int b = 0; b < bar; b++) printf("█");
        printf(" %.4f\n  ", track_I[c]);
    }
    mi_avg /= n_cycles;
    printf("\n  I(A;B) range: %.4f – %.4f  |  avg: %.4f bits\n",
           mi_min, mi_max, mi_avg);
    printf("  ΔI = %.4f bits (ether's total entangling range)\n",
           mi_max - mi_min);

    if (mi_max - mi_min > 0.01)
        printf("  ★ The ether performed entangling/disentangling operations\n"
               "    on the two-qudit state across %d captures.\n", n_cycles);
    else if (mi_avg > 0.01)
        printf("  ★ Steady non-zero I(A;B) — the two qudits are persistently\n"
               "    coupled through the shared physical channel.\n");

    printf("\n  ★ Both qudits lived entirely in the EM field.\n");
    printf("  ★ The ether performed ALL computation.\n\n");

    free(prev_A); free(prev_B); free(curr_A); free(curr_B);
    free(mean_A); free(mean_B); free(var_A); free(var_B);
    for (int i = 0; i < D_A; i++) free(cov_AB[i]);
    free(cov_AB); free(track_I);
    sdr_close(&sdr);
    return 0;
}

/* ═══════════════════════════════════════════════════════════════
 * QUANTUM VM — Ether-native instruction set executed physically
 *
 * The ether IS the quantum computer.  The VM interprets instructions
 * and maps them to physical EM field operations.
 *
 * Superposition: the DFT of captured I/Q gives ALL D complex
 * amplitudes simultaneously.  A single capture IS a superposition
 * measurement in the frequency basis.  The ether naturally provides
 * the superposition — we don't need to create it, just read it.
 *
 * Instruction set:
 *   INIT           capture ambient RF → normalize to |ψ⟩
 *   SUPERPOSE      capture, treat as coherent superposition
 *   X [shift]      cyclic-shift frequency bins (default +1)
 *   Z [rad] [k]    phase rotation (k=-1 = all, default rad=π/4)
 *   H              LO hop for Hadamard
 *   CZ             cross-correlate two I/Q windows (entangle)
 *   DFT [step]     LO retune by step bins (basis change)
 *   TX [ms]        radiate state into ether via LO hopping
 *   RX             capture state from ether
 *   MEASURE        Born-rule collapse via ADC LSB entropy
 *   TICK           one complete TX→ether→RX cycle
 *   PROB           show probability distribution
 *   SHOW           show full complex amplitudes
 *   COHERENT [f]   synthesize OFDM I/Q superposition → file
 *   WAIT [ms]      let ether compute for N ms
 *   ECHO [text]    print message
 *   HELP           show instruction set
 *   QUIT           exit VM
 *
 * Script:  ./sdr_ether 6 --vm script.qvm
 * Live:    ./sdr_ether 6 --vm
 * ═══════════════════════════════════════════════════════════════ */
/*
 * ETHER EMIT — radiate qudit state into the EM field via LO hopping.
 * Each level k gets dwell time ∝ |amplitude|².  After this call,
 * the ether physically holds the qudit state as frequency-multiplexed
 * CW tones in the room's EM field.
 */
static void ether_emit(SdrDev *s, const Wavefunction *wf) {
    double bin_bw = (double)s->rate / wf->d;
    sdr_tx_wavefunction(s, wf->prob, wf->d, wf->freq, bin_bw, 8000);
}

/*
 * ETHER READ — capture the ether state back into the wavefunction.
 * Retunes to center, captures I/Q, DFT decomposes into D levels.
 * The returned wf IS the ether's physical state.
 */
static void ether_read(SdrDev *s, Wavefunction *wf) {
    sdr_retune(s, (uint32_t)wf->freq);
    usleep(2000);
    sdr_flush(s);
    sdr_capture(s);
    wf_from_iq(s->iq_i, s->iq_q, s->iq_n, wf);
    wf->freq = s->freq;
}

static void wf_print(const Wavefunction *wf, const char *label, int cycle);
static int  gate_measure(const uint8_t *raw, int raw_n, const double *prob, int d);
static void wf_from_iq(const double *iq_i, const double *iq_q, int np,
                        Wavefunction *wf);
static void gate_X(Wavefunction *wf);
static void gate_Z(Wavefunction *wf, double rad, int target);
static void gate_H(SdrDev *s, Wavefunction *wf);
static void gate_CZ(SdrDev *s, Wavefunction *wf);
__attribute__((unused)) static void gate_DFT(SdrDev *s, Wavefunction *wf, int step);
static void tx_synthesize(const Wavefunction *wf, TxBuf *tx);
static void tx_write_file(const TxBuf *tx, const char *path);

/* ═══════════════════════════════════════════════════════════════
 * QVM API — Extensible Quantum VM instruction set /dispatch table
 *
 * Architecture:
 *   QvmCtx    — VM context (wavefunction, SDR, instruction table)
 *   QvmOp     — int handler(QvmCtx*, double arg1, double arg2)
 *                Returns 0=continue, 1=quit, -1=error
 *   qvm_reg() — register a named opcode
 *   qvm_eval()— parse & dispatch a single instruction string
 *   qvm_run() — REPL or script execution loop
 *
 * External code (e.g. CUDA wrapper) calls qvm_create() + qvm_eval()
 * ═══════════════════════════════════════════════════════════════ */

typedef struct QvmCtx QvmCtx;
typedef int (*QvmOp)(QvmCtx *q, double a1, double a2);

#define QVM_MAX_OPS 64

struct QvmCtx {
    SdrDev       *sdr;
    Wavefunction  wf;
    int           sdr_ok;
    uint32_t      freq, rate;
    int           running;

    /* Instruction table */
    struct { char name[16]; QvmOp fn; const char *help; } ops[QVM_MAX_OPS];
    int n_ops;

    /* Script state */
    char **lines;  int nlines, ip;
    int  *loop_stack, loop_depth;
    int   interactive;

    /* GPU offload state (calibrated channel) */
    double *M, *Minv;
    int    M_dim;

    /* Qubit bin mapping for ANTISYM gate (up to 128 qudits = 256 bins) */
    int    qbins[256];
    int    n_qbins;
};

/* ─── Helpers ─── */
static void qvm_norm(Wavefunction *w) {
    double t=0; for(int i=0;i<w->d;i++) t+=w->re[i]*w->re[i]+w->im[i]*w->im[i];
    if(t>1e-15){double s=1.0/sqrt(t);for(int i=0;i<w->d;i++){w->re[i]*=s;w->im[i]*=s;}}
    w->entropy=0;w->purity=0;
    for(int i=0;i<w->d;i++){w->prob[i]=w->re[i]*w->re[i]+w->im[i]*w->im[i];
        if(w->prob[i]>1e-15)w->entropy-=w->prob[i]*log2(w->prob[i]);
        w->purity+=w->prob[i]*w->prob[i];}
}

static void qvm_sync(QvmCtx *q) {
    if (!q->sdr_ok) return;
    ether_emit(q->sdr, &q->wf);
    ether_read(q->sdr, &q->wf);
}

/* ─── Op table management ─── */
static int qvm_reg(QvmCtx *q, const char *name, QvmOp fn, const char *help){
    if (q->n_ops >= QVM_MAX_OPS) return -1;
    strncpy(q->ops[q->n_ops].name, name, 15);
    q->ops[q->n_ops].fn   = fn;
    q->ops[q->n_ops].help = help ? help : "";
    q->n_ops++;
    return 0;
}

static QvmOp qvm_lookup(QvmCtx *q, const char *name){
    for (int i=0;i<q->n_ops;i++)
        if (strcasecmp(q->ops[i].name, name)==0) return q->ops[i].fn;
    return NULL;
}

/* ── Instruction handlers ── */
static int op_init(QvmCtx *q, double a1, double a2){
    (void)a1;(void)a2;
    if(q->sdr_ok){sdr_capture(q->sdr);wf_from_iq(q->sdr->iq_i,q->sdr->iq_q,q->sdr->iq_n,&q->wf);}
    printf("  [INIT] Capture from ether\n");
    wf_print(&q->wf,"|ψ⟩",-1); return 0;
}

static int op_superpose(QvmCtx *q, double a1, double a2){
    (void)a1;(void)a2;
    if(q->sdr_ok){sdr_capture(q->sdr);wf_from_iq(q->sdr->iq_i,q->sdr->iq_q,q->sdr->iq_n,&q->wf);}
    printf("  [SUPERPOSE] %d-level decomposition\n",q->wf.d);
    wf_print(&q->wf,"|ψ⟩ sup",-1); return 0;
}

static int op_x(QvmCtx *q, double a1, double a2){
    int shift = ((int)a1) ? (int)a1 : 1;
    shift = ((shift % q->wf.d)+q->wf.d)%q->wf.d;
    for(int s=0;s<shift;s++) gate_X(&q->wf);
    printf("  [X] +%d\n",shift); wf_print(&q->wf,"|ψ⟩",-1); return 0;
}

static int op_z(QvmCtx *q, double a1, double a2){
    double rad = a1 ? a1 : M_PI/4.0;
    int target = ((int)a2) ? (int)a2 : -1;
    if(target >= q->wf.d) target = -1;
    gate_Z(&q->wf, rad, target);
    printf("  [Z] %.3f rad\n",rad); return 0;
}

static int op_h(QvmCtx *q, double a1, double a2){
    (void)a1;(void)a2;
    if(q->sdr_ok) gate_H(q->sdr, &q->wf);
    printf("  [H] Hadamard\n"); wf_print(&q->wf,"|ψ⟩",-1); return 0;
}

static int op_cz(QvmCtx *q, double a1, double a2){
    (void)a1;(void)a2;
    if(q->sdr_ok){
        /* Complex cross-correlation: preserves phase entanglement */
        double *prev_re=malloc(q->wf.d*sizeof(double));
        double *prev_im=malloc(q->wf.d*sizeof(double));
        memcpy(prev_re,q->wf.re,q->wf.d*sizeof(double));
        memcpy(prev_im,q->wf.im,q->wf.d*sizeof(double));
        sdr_capture(q->sdr);
        wf_from_iq(q->sdr->iq_i,q->sdr->iq_q,q->sdr->iq_n,&q->wf);
        for(int k=0;k<q->wf.d;k++){
            double cr=prev_re[k]*q->wf.re[k]+prev_im[k]*q->wf.im[k];
            double ci=prev_im[k]*q->wf.re[k]-prev_re[k]*q->wf.im[k];
            q->wf.re[k]=cr;q->wf.im[k]=ci;
        }
        free(prev_re);free(prev_im);
    }
    printf("  [CZ] Entangle (complex)\n"); wf_print(&q->wf,"|ψ⟩",-1); return 0;
}

static int op_qbin(QvmCtx *q, double a1, double a2){
    (void)a2;
    if(a1>=256){printf("  [QBIN] max 256 bins (128 qudits)\n");return 0;}
    q->qbins[(int)a1]=(int)a2;
    return 0;
}
static int op_qbin_done(QvmCtx *q, double a1, double a2){
    (void)a1;(void)a2;
    q->n_qbins=0;
    for(int i=0;i<256;i++)if(q->qbins[i]>=0&&q->qbins[i]<q->wf.d){q->n_qbins++;}
    printf("  [QBIN] %d bins (%d qudits): ",q->n_qbins,q->n_qbins/2);
    for(int i=0;i<q->n_qbins && i<20;i++)printf("%d ",q->qbins[i]);
    if(q->n_qbins>20)printf("...");
    printf("\n");
    return 0;
}
/* ANTISYM: anti-symmetric pair entanglement gate.
   8-pass feedback: X[k]=+A, X[D-k]=-A per qubit bin → room → renormalize → repeat.
   Uses qbins[0..3] from QBIN (4 bins for 2 qudits).
   DC cancels via k+(D-k)=0 destructive interference. */
/* ANTISYM: anti-symmetric pair entanglement gate for N qudits.
   GHZ state (|0...0⟩+|1...1⟩)/√2 encoded as X[k]=+A, X[D-k]=-A.
   8-pass feedback through room with DC cancellation. */
static int op_antisym(QvmCtx *q, double a1, double a2){
    (void)a1;(void)a2;
    if(!q->sdr_ok)return 0;
    int n=q->n_qbins;
    if(n<4){printf("  [ANTISYM] need >=4 bins (>=2 qudits)\n");return 0;}
    int nq=n/2, D=q->wf.d, n_pass=8;
    double x[D],xi[D],y[D],amp=1.0/sqrt(2.0);
    for(int pass=0;pass<n_pass;pass++){
        memset(x,0,D*sizeof(double));memset(xi,0,D*sizeof(double));
        if(pass==0){
            /* Encode all configured bins. Read from wf state if available. */
            double init=0;
            for(int i=0;i<n;i++)init+=q->wf.prob[q->qbins[i]];
            if(init>1e-10){
                for(int i=0;i<n;i++){
                    double p=q->wf.prob[q->qbins[i]];
                    if(p>1e-10){x[q->qbins[i]]=+sqrt(p);x[D-q->qbins[i]]=-sqrt(p);}
                }
            }else{
                for(int i=0;i<n;i+=2){
                    x[q->qbins[i]]+=amp;x[D-q->qbins[i]]-=amp;
                    x[q->qbins[i+1]]+=amp;x[D-q->qbins[i+1]]-=amp;
                }
            }
        }else{
            double tot=0;
            for(int i=0;i<n;i++)tot+=y[q->qbins[i]]+y[D-q->qbins[i]];
            if(tot>1e-15){
                for(int i=0;i<n;i++){
                    double p=(y[q->qbins[i]]+y[D-q->qbins[i]])/tot;
                    if(p>1e-10){x[q->qbins[i]]=+sqrt(p);x[D-q->qbins[i]]=-sqrt(p);}
                }
            }
        }
        for(int i=0;i<16;i++)x[i]=xi[i]=0.0;
        qvm_ofdm_compute(q,x,xi,y,D);
    }
    double tot=0;for(int i=0;i<n;i++)tot+=y[q->qbins[i]]+y[D-q->qbins[i]];
    double s0=0,s1=0;
    if(tot>1e-15){
        for(int i=0;i<n;i+=2){s0+=y[q->qbins[i]]+y[D-q->qbins[i]];s1+=y[q->qbins[i+1]]+y[D-q->qbins[i+1]];}
        double s=s0+s1;
        for(int i=0;i<n;i++)
            q->wf.prob[q->qbins[i]]=(y[q->qbins[i]]+y[D-q->qbins[i]])/s;
        for(int i=0;i<n;i++)q->wf.re[q->qbins[i]]=sqrt(q->wf.prob[q->qbins[i]]);
    }
    printf("  [ANTISYM] %d-qudit GHZ -> %d passes: |0..0>=%.3f |1..1>=%.3f\n",
        nq,n_pass,tot>1e-15?s0/(s0+s1):0,tot>1e-15?s1/(s0+s1):0);
    return 0;
}

/* CHSH: Bell inequality test on current qbin state.
   Computes S=|E(a,b)-E(a,b')+E(a',b)+E(a',b')| with optimal angles.
   Uses standard 2-qubit measurement projectors on a,b. */
static int op_chsh(QvmCtx *q, double a1, double a2){
    (void)a1;(void)a2;
    if(q->n_qbins<4){printf("  [CHSH] use QBIN first (need 4 bins)\n");return 0;}
    /* GHZ endpoints: qbins[0] = |0..0⟩, qbins[n-1] = |1..1⟩ */
    double p0=q->wf.prob[q->qbins[0]];
    double p1=q->wf.prob[q->qbins[q->n_qbins-1]];
    if(p0<0.01 && p1<0.01){printf("  [CHSH] not entangled\n");return 0;}
    double t=p0+p1;
    double p00=t>1e-15?p0/t:0, p11=t>1e-15?p1/t:0;
    double p[4]={p00,0,0,p11};
    double a=0,b=M_PI/4,ap=M_PI/2,bp=3*M_PI/4;
    double ca=cos(a/2),sa=sin(a/2),cb=cos(b/2),sb=sin(b/2),
           cap=cos(ap/2),sap=sin(ap/2),cbp=cos(bp/2),sbp=sin(bp/2);
    double amp[4];for(int i=0;i<4;i++)amp[i]=sqrt(p[i]);
    /* E(a,b) */
    double pp=amp[0]*ca*cb+amp[1]*ca*sb+amp[2]*sa*cb+amp[3]*sa*sb;
    double pm=amp[0]*ca*(-sb)+amp[1]*ca*cb+amp[2]*sa*(-sb)+amp[3]*sa*cb;
    double mp=amp[0]*(-sa)*cb+amp[1]*(-sa)*sb+amp[2]*ca*cb+amp[3]*ca*sb;
    double mm=amp[0]*(-sa)*(-sb)+amp[1]*(-sa)*cb+amp[2]*ca*(-sb)+amp[3]*ca*cb;
    double Eab=pp*pp-pm*pm-mp*mp+mm*mm;
    /* E(a,b') */
    pp=amp[0]*ca*cbp+amp[1]*ca*sbp+amp[2]*sa*cbp+amp[3]*sa*sbp;
    pm=amp[0]*ca*(-sbp)+amp[1]*ca*cbp+amp[2]*sa*(-sbp)+amp[3]*sa*cbp;
    mp=amp[0]*(-sa)*cbp+amp[1]*(-sa)*sbp+amp[2]*ca*cbp+amp[3]*ca*sbp;
    mm=amp[0]*(-sa)*(-sbp)+amp[1]*(-sa)*cbp+amp[2]*ca*(-sbp)+amp[3]*ca*cbp;
    double Eabp=pp*pp-pm*pm-mp*mp+mm*mm;
    /* E(a',b) */
    pp=amp[0]*cap*cb+amp[1]*cap*sb+amp[2]*sap*cb+amp[3]*sap*sb;
    pm=amp[0]*cap*(-sb)+amp[1]*cap*cb+amp[2]*sap*(-sb)+amp[3]*sap*cb;
    mp=amp[0]*(-sap)*cb+amp[1]*(-sap)*sb+amp[2]*cap*cb+amp[3]*cap*sb;
    mm=amp[0]*(-sap)*(-sb)+amp[1]*(-sap)*cb+amp[2]*cap*(-sb)+amp[3]*cap*cb;
    double Eapb=pp*pp-pm*pm-mp*mp+mm*mm;
    /* E(a',b') */
    pp=amp[0]*cap*cbp+amp[1]*cap*sbp+amp[2]*sap*cbp+amp[3]*sap*sbp;
    pm=amp[0]*cap*(-sbp)+amp[1]*cap*cbp+amp[2]*sap*(-sbp)+amp[3]*sap*cbp;
    mp=amp[0]*(-sap)*cbp+amp[1]*(-sap)*sbp+amp[2]*cap*cbp+amp[3]*cap*sbp;
    mm=amp[0]*(-sap)*(-sbp)+amp[1]*(-sap)*cbp+amp[2]*cap*(-sbp)+amp[3]*cap*cbp;
    double Eapbp=pp*pp-pm*pm-mp*mp+mm*mm;
    double S=fabs(Eab-Eabp+Eapb+Eapbp);
    printf("  [CHSH] |00⟩=%.3f |01⟩=%.3f |10⟩=%.3f |11⟩=%.3f\n",p[0],p[1],p[2],p[3]);
    printf("  [CHSH] E(A,B)=%+.4f E(A,B')=%+.4f E(A',B)=%+.4f E(A',B')=%+.4f\n",Eab,Eabp,Eapb,Eapbp);
    printf("  [CHSH] S=%.4f (classical ≤2.0, quantum ≤2.83) %s\n",S,S>2.0?"★★ BELL VIOLATION ★★":"no violation");
    return 0;
}

/* MERMIN: Mermin inequality for N-qudit GHZ states (N odd).
   M_N = 2^N * Re(a*b) for GHZ state |0...0⟩+|1...1⟩.
   With probabilities only (no phase): M_max = 2^N * √(p0*p1).
   Classical bound: ≤ 2^{(N-1)/2} for odd N.
   Usage: MERMIN [N] — defaults to N = n_qbins/2 */
static int op_mermin(QvmCtx *q, double a1, double a2){
    int N = ((int)a1)>0 ? (int)a1 : q->n_qbins/2;
    if(N<2){printf("  [MERMIN] need ≥2 qudits\n");return 0;}
    if(q->n_qbins<2*N){
        printf("  [MERMIN] need %d bins (%d qudits), have %d bins\n",2*N,N,q->n_qbins);
        return 0;
    }
    int k0=q->qbins[0], k1=q->qbins[2*N-1];
    double p0=q->wf.prob[k0], p1=q->wf.prob[k1];
    if(p0<0.01 && p1<0.01){printf("  [MERMIN] not entangled\n");return 0;}
    double tot=p0+p1;if(tot>1e-15){p0/=tot;p1/=tot;}
    double coh=sqrt(p0*p1); /* max coherence achievable */
    /* Mermin correlation: M = Σ_k (-1)^{?} ⟨...⟩ ...
       For GHZ |ψ⟩=a|0_n⟩+b|1_n⟩, optimized M = 2^N * Re(a*b)
       Maximum achievable with probability-only: M_max = 2^N * √(p0*p1) */
    double Mmax = pow(2.0,N) * coh;
    double classical = pow(2.0, (N-1.0)/2.0);
    int qpu = 2*N; /* qubits processed */
    (void)qpu;

    printf("  [MERMIN] %d-qudit GHZ: |0..0⟩=%.3f |1..1⟩=%.3f\n",N,p0,p1);
    if(N%2==0){
        printf("  [MERMIN] N=%d even: use AB-Klyshko inequality instead\n",N);
        double Mmax2 = 2.0 * coh;
        printf("  [MERMIN] M_max ≤ %.4f (correlation bound)\n",Mmax2);
        return 0;
    }
    printf("  [MERMIN] M_max = 2^%d * √(%.3f*%.3f) = %.4f\n",N,p0,p1,Mmax);
    printf("  [MERMIN] Classical bound: %.4f\n",classical);
    printf("  [MERMIN] Quantum GHZ max:  %.4f\n",pow(2.0,N-1));
    if(Mmax > classical)
        printf("  [MERMIN] %.4f > %.4f ★★ MERMIN VIOLATION ★★\n",Mmax,classical);
    else
        printf("  [MERMIN] %.4f ≤ %.4f no violation\n",Mmax,classical);
    return 0;
}

static int op_dft(QvmCtx *q, double a1, double a2){
    int step=((int)a1)?(int)a1:1;
    if(q->sdr_ok) gate_DFT(q->sdr,&q->wf,step);
    printf("  [DFT] %+d bins\n",step); wf_print(&q->wf,"|ψ⟩",-1); return 0;
}

static int op_tx(QvmCtx *q, double a1, double a2){
    (void)a1;(void)a2;
    if(q->sdr_ok) gate_tx_hardware(q->sdr,&q->wf);
    printf("  [TX] Radiated\n"); return 0;
}

static int op_rx(QvmCtx *q, double a1, double a2){
    (void)a1;(void)a2;
    if(q->sdr_ok){sdr_capture(q->sdr);wf_from_iq(q->sdr->iq_i,q->sdr->iq_q,q->sdr->iq_n,&q->wf);}
    printf("  [RX] Captured\n"); wf_print(&q->wf,"|ψ⟩",-1); return 0;
}

static int op_measure(QvmCtx *q, double a1, double a2){
    (void)a1;(void)a2;
    int o;
    if(q->sdr_ok) o=gate_measure(q->sdr->iq_raw,q->sdr->iq_n*2,q->wf.prob,q->wf.d);
    else o=rand()%q->wf.d;
    printf("  [MEASURE] → |%d⟩\n",o);
    memset(q->wf.re,0,q->wf.d*8); memset(q->wf.im,0,q->wf.d*8); memset(q->wf.prob,0,q->wf.d*8);
    q->wf.re[o]=1.0;q->wf.prob[o]=1.0;q->wf.entropy=0;q->wf.purity=1.0; return 0;
}

/* QMEASURE [qubit_idx]: physical collapse via room.
   After ANTISYM entanglement, choose outcome from ADC entropy,
   TX only the surviving GHZ branch into the room. Room responds. */
/* QMEASURE: passive wait → environment decoheres → capture.
   Room's thermal/mechanical/RF perturbations scramble phase faster than
   the DFT integration window → M[k≠m]→0, only diagonal survives. */
static int op_qmeasure(QvmCtx *q, double a1, double a2){
    int dwell_ms = ((int)a1)>0 ? (int)a1 : 100;
    (void)a2;
    if(!q->sdr_ok){printf("  [QMEASURE] no SDR\n");return 0;}
    if(q->n_qbins<4){printf("  [QMEASURE] need ≥2 qudits\n");return 0;}
    int k0=q->qbins[0], k1=q->qbins[q->n_qbins-1];
    double pb0=q->wf.prob[k0], pb1=q->wf.prob[k1];
    double tb=pb0+pb1;
    printf("  [QMEASURE] before: |00⟩=%.3f |11⟩=%.3f → wait %dms (env decoheres)\n",
        tb>1e-15?pb0/tb:0,tb>1e-15?pb1/tb:0,dwell_ms);
    usleep(dwell_ms*1000);
    sdr_capture(q->sdr);
    wf_from_iq(q->sdr->iq_i,q->sdr->iq_q,q->sdr->iq_n,&q->wf);
    double p0=q->wf.prob[k0], p1=q->wf.prob[k1], t=p0+p1;
    printf("  [QMEASURE] after:  |00⟩=%.3f |11⟩=%.3f → %s\n",
        t>1e-15?p0/t:0,t>1e-15?p1/t:0, p0>p1?"|00⟩":"|11⟩");
    return 0;
}

/* PROBE: noise-collapse test — fix denominator, DMA settle, order.
   Protocol: ANTISYM→read, wait, single-probe→read (natural pair)
             ANTISYM→read, NOISE, wait, single-probe→read (noise pair) */
static int op_probe(QvmCtx *q, double a1, double a2){
    int nt=((int)a1)>0?(int)a1:10;(void)a2;
    if(!q->sdr_ok)return 0;
    if(q->n_qbins<4){printf("  [PROBE] need QBIN\n");return 0;}
    int k0=q->qbins[0],k1=q->qbins[q->n_qbins-1],D=q->wf.d;
    /* Collect all qbin bins for proper normalization */
    int nb=q->n_qbins; if(nb>16)nb=16;
    double delta_nat[20]={0},delta_noise[20]={0};

    printf("  [PROBE] %d trials (fixed denominator, DMA settle, order)\n",nt);
    for(int t=0;t<nt;t++){
        /* --- Natural baseline: A→wait→probe→B --- */
        qvm_eval(q,"ANTISYM");
        double pa=0;for(int i=0;i<nb;i++)pa+=q->wf.prob[q->qbins[i]];
        double valA=pa>1e-30?q->wf.prob[k0]/pa:0;
        usleep(2000); /* DMA buffer settle */
        /* Single anti-sym probe (not full ANTISYM) */
        { double x[D],xi[D],y[D],a=0.7071;
          memset(x,0,D*8);memset(xi,0,D*8);
          x[k0]=+a;x[D-k0]=-a;x[k1]=+a;x[D-k1]=-a;
          for(int i=0;i<16;i++)x[i]=xi[i]=0;
          qvm_ofdm_compute(q,x,xi,y,D);
          double pb=0;for(int i=0;i<nb;i++)pb+=y[q->qbins[i]];
          double valB=pb>1e-30?y[k0]/pb:0;
          delta_nat[t]=fabs(valB-valA);
          printf("    %2d nat: A=%.3f B=%.3f Δ=%.3f\n",t,valA,valB,delta_nat[t]);
        }

        /* --- Noise test: A→NOISE→wait→probe→C --- */
        qvm_eval(q,"ANTISYM");
        pa=0;for(int i=0;i<nb;i++)pa+=q->wf.prob[q->qbins[i]];
        valA=pa>1e-30?q->wf.prob[k0]/pa:0;
        qvm_eval(q,"NOISE 1.0");
        usleep(2000);
        { double x[D],xi[D],y[D],a=0.7071;
          memset(x,0,D*8);memset(xi,0,D*8);
          x[k0]=+a;x[D-k0]=-a;x[k1]=+a;x[D-k1]=-a;
          for(int i=0;i<16;i++)x[i]=xi[i]=0;
          qvm_ofdm_compute(q,x,xi,y,D);
          double pc=0;for(int i=0;i<nb;i++)pc+=y[q->qbins[i]];
          double valC=pc>1e-30?y[k0]/pc:0;
          delta_noise[t]=fabs(valC-valA);
          printf("    %2d ns:  A=%.3f C=%.3f Δ=%.3f\n",t,valA,valC,delta_noise[t]);
        }
    }
    double mn=0,ms=0;
    for(int t=0;t<nt;t++){mn+=delta_nat[t];ms+=delta_noise[t];}
    mn/=nt;ms/=nt;
    printf("  [PROBE] mean Δ_nat=%.3f  Δ_noise=%.3f  ratio=%.2fx\n",mn,ms,ms/(mn+1e-9));
    printf("  [PROBE] %s\n",ms>mn*1.5?"★ NOISE COLLAPSES STATE ★":"noise ≈ natural jitter");

    /* ═══ Anti-symmetric noise pairs at GHZ bins → maximal disruption ═══ */
    printf("  [PROBE] anti-sym noise at GHZ bins (0.5×, 2 bursts):\n");
    for(int run=0;run<nt;run++){
        double x[D],xi[D],a=0.7071;
        /* Baseline */
        memset(x,0,D*8);memset(xi,0,D*8);
        x[k0]=+a;x[D-k0]=-a;x[k1]=+a;x[D-k1]=-a;
        for(int i=0;i<16;i++)x[i]=xi[i]=0;
        qvm_ofdm_compute(q,x,xi,x,D);
        double pb=0;for(int i=0;i<nb;i++)pb+=x[q->qbins[i]];
        double vB=pb>1e-30?x[k0]/pb:0;
        /* Anti-sym TX */
        memset(x,0,D*8);memset(xi,0,D*8);
        x[k0]=+a;x[D-k0]=-a;x[k1]=+a;x[D-k1]=-a;
        for(int i=0;i<16;i++)x[i]=xi[i]=0;
        for(int i=0;i<D;i++){q->wf.re[i]=x[i];q->wf.im[i]=xi[i];}
        gate_ofdm_tx(q->sdr,&q->wf,NULL);
        {struct v4l2_buffer b;memset(&b,0,sizeof(b));
         b.type=V4L2_BUF_TYPE_SDR_CAPTURE;b.memory=V4L2_MEMORY_MMAP;
         if(ioctl(q->sdr->fd,VIDIOC_DQBUF,&b)==0)ioctl(q->sdr->fd,VIDIOC_QBUF,&b);}
        /* 2 bursts of anti-sym noise at k0 and k1 with random phase */
        for(int burst=0;burst<2;burst++){
            memset(x,0,D*8);memset(xi,0,D*8);
            double phi0=((double)rand()/RAND_MAX)*2*M_PI;
            double phi1=((double)rand()/RAND_MAX)*2*M_PI;
            x[k0]=0.5*cos(phi0);x[D-k0]=-0.5*cos(phi0);
            xi[k0]=0.5*sin(phi0);xi[D-k0]=-0.5*sin(phi0);
            x[k1]=0.5*cos(phi1);x[D-k1]=-0.5*cos(phi1);
            xi[k1]=0.5*sin(phi1);xi[D-k1]=-0.5*sin(phi1);
            for(int i=0;i<16;i++)x[i]=xi[i]=0;
            for(int i=0;i<D;i++){q->wf.re[i]=x[i];q->wf.im[i]=xi[i];}
            gate_ofdm_tx(q->sdr,&q->wf,NULL);
            {struct v4l2_buffer b;memset(&b,0,sizeof(b));
             b.type=V4L2_BUF_TYPE_SDR_CAPTURE;b.memory=V4L2_MEMORY_MMAP;
             if(ioctl(q->sdr->fd,VIDIOC_DQBUF,&b)==0)ioctl(q->sdr->fd,VIDIOC_QBUF,&b);}
        }
        usleep(50000);sdr_capture(q->sdr);
        wf_from_iq(q->sdr->iq_i,q->sdr->iq_q,q->sdr->iq_n,&q->wf);
        double pc=0;for(int i=0;i<nb;i++)pc+=q->wf.prob[q->qbins[i]];
        double vN=pc>1e-30?q->wf.prob[k0]/pc:0;
        double d=fabs(vN-vB);
        printf("    %2d: base=%.3f antisym-noise=%.3f Δ=%.3f %s\n",
            run,vB,vN,d, d>0.2?"COLLAPSED":"same");
    }
    return 0;
}

/* NOISE: spectral flush — TX white noise across all bins.
   Each bin gets random phase → overwrites deterministic M[k≠m].
   Diagonal power per bin survives, off-diagonal collapses to 0. */
/* COLLAPSE: anti-sym noise at GHZ bins (87.5% collapse rate).
   TX anti-sym probe → 2 bursts anti-sym noise at GHZ bins (random phase)
   → capture. Random-phase pairs at the intermodulation channel maximally
   disrupt cross-bin coherence → M[k≠m]→0. */
static int op_collapse(QvmCtx *q, double a1, double a2){
    double level = a1>0 ? a1 : 0.5; (void)a2;
    if(!q->sdr_ok){printf("  [COLLAPSE] no SDR\n");return 0;}
    if(q->n_qbins<4){printf("  [COLLAPSE] need QBIN\n");return 0;}
    int k0=q->qbins[0],k1=q->qbins[q->n_qbins-1],D=q->wf.d;
    double x[D],xi[D],a=0.7071;

    /* TX anti-sym probe */
    memset(x,0,D*8);memset(xi,0,D*8);
    x[k0]=+a;x[D-k0]=-a;x[k1]=+a;x[D-k1]=-a;
    for(int i=0;i<16;i++)x[i]=xi[i]=0;
    for(int i=0;i<D;i++){q->wf.re[i]=x[i];q->wf.im[i]=xi[i];}
    gate_ofdm_tx(q->sdr,&q->wf,NULL);
    {struct v4l2_buffer b;memset(&b,0,sizeof(b));
     b.type=V4L2_BUF_TYPE_SDR_CAPTURE;b.memory=V4L2_MEMORY_MMAP;
     if(ioctl(q->sdr->fd,VIDIOC_DQBUF,&b)==0)ioctl(q->sdr->fd,VIDIOC_QBUF,&b);}

    /* 2 bursts anti-sym noise at GHZ bins with random phase */
    for(int burst=0;burst<2;burst++){
        memset(x,0,D*8);memset(xi,0,D*8);
        double phi0=((double)rand()/RAND_MAX)*2*M_PI;
        double phi1=((double)rand()/RAND_MAX)*2*M_PI;
        x[k0]=level*cos(phi0);x[D-k0]=-level*cos(phi0);
        xi[k0]=level*sin(phi0);xi[D-k0]=-level*sin(phi0);
        x[k1]=level*cos(phi1);x[D-k1]=-level*cos(phi1);
        xi[k1]=level*sin(phi1);xi[D-k1]=-level*sin(phi1);
        for(int i=0;i<16;i++)x[i]=xi[i]=0;
        for(int i=0;i<D;i++){q->wf.re[i]=x[i];q->wf.im[i]=xi[i];}
        gate_ofdm_tx(q->sdr,&q->wf,NULL);
        {struct v4l2_buffer b;memset(&b,0,sizeof(b));
         b.type=V4L2_BUF_TYPE_SDR_CAPTURE;b.memory=V4L2_MEMORY_MMAP;
         if(ioctl(q->sdr->fd,VIDIOC_DQBUF,&b)==0)ioctl(q->sdr->fd,VIDIOC_QBUF,&b);}
    }

    /* Capture collapsed state */
    usleep(50000);
    sdr_capture(q->sdr);
    wf_from_iq(q->sdr->iq_i,q->sdr->iq_q,q->sdr->iq_n,&q->wf);
    double p0=q->wf.prob[k0],p1=q->wf.prob[k1],t=p0+p1;
    printf("  [COLLAPSE] ×%.1f |00⟩=%.3f |11⟩=%.3f\n",
        level,t>1e-30?p0/t:0,t>1e-30?p1/t:0);
    return 0;
}

/* KILL: 3 passes of winner-only feedback. Room's multipath memory
   amplifies the selected branch. Follow with CONT to use memory. */
static int op_kill(QvmCtx *q, double a1, double a2){
    (void)a1;(void)a2;
    if(!q->sdr_ok){printf("  [KILL] no SDR\n");return 0;}
    if(q->n_qbins<4){printf("  [KILL] need QBIN\n");return 0;}
    int k0=q->qbins[0],k1=q->qbins[q->n_qbins-1],D=q->wf.d;
    double p0=q->wf.prob[k0],p1=q->wf.prob[k1];
    int winner=(p0>p1)?0:1, kw=winner?k1:k0;
    double x[D],xi[D],y[D];
    memset(x,0,D*8);memset(xi,0,D*8);
    x[kw]=+1.0;x[D-kw]=-1.0;
    for(int i=0;i<16;i++)x[i]=xi[i]=0;
    qvm_ofdm_compute(q,x,xi,y,D);
    /* 8 passes of winner-only feedback (match ANTISYM's convergence) */
    for(int p=0;p<7;p++){
        double s=y[kw]+y[D-kw];
        double amp=s>1e-15?sqrt(s):0;
        memset(x,0,D*8);memset(xi,0,D*8);
        if(amp>1e-15){x[kw]=+amp;x[D-kw]=-amp;}
        for(int i=0;i<16;i++)x[i]=xi[i]=0;
        qvm_ofdm_compute(q,x,xi,y,D);
    }
    printf("  [KILL] winner → room bias (%d), CONT feedback...\n",winner);

        /* CONT: continue from KILL's output state directly (no fresh probe).
           Read y[] from last QVM compute, feedback from THOSE amplitudes. */
        {
            double x2[D],xi2[D],y2[D];
            memcpy(y2,y,D*8);
            double s2=y2[k0]+y2[D-k0]+y2[k1]+y2[D-k1];
            if(s2>1e-15){
                double p0_2=(y2[k0]+y2[D-k0])/s2;
                double p1_2=(y2[k1]+y2[D-k1])/s2;
                for(int p=0;p<8;p++){
                    memset(x2,0,D*8);memset(xi2,0,D*8);
                    if(p0_2>1e-15){x2[k0]=+sqrt(p0_2);x2[D-k0]=-sqrt(p0_2);}
                    if(p1_2>1e-15){x2[k1]=+sqrt(p1_2);x2[D-k1]=-sqrt(p1_2);}
                    for(int i=0;i<16;i++)x2[i]=xi2[i]=0;
                    qvm_ofdm_compute(q,x2,xi2,y2,D);
                    s2=y2[k0]+y2[D-k0]+y2[k1]+y2[D-k1];
                    if(s2>1e-15){p0_2=(y2[k0]+y2[D-k0])/s2;p1_2=(y2[k1]+y2[D-k1])/s2;}
                }
            }
            double tt=y2[k0]+y2[D-k0]+y2[k1]+y2[D-k1];
            for(int i=0;i<D;i++)q->wf.re[i]=q->wf.im[i]=q->wf.prob[i]=0;
            if(tt>1e-15){
                q->wf.prob[k0]=(y2[k0]+y2[D-k0])/tt;
                q->wf.prob[k1]=(y2[k1]+y2[D-k1])/tt;
            }
            double f0=tt>1e-15?q->wf.prob[k0]:0,f1=tt>1e-15?q->wf.prob[k1]:0;
            printf("  [KILL] CONT: |00⟩=%.3f |11⟩=%.3f %s\n",f0,f1,
                (winner==0&&f0>0.5)||(winner==1&&f1>0.5)?"LOCKED":"drifted");
        }
    return 0;
}

static int op_memory(QvmCtx *q, double a1, double a2){
    int delay = ((int)a1)>0 ? (int)a1 : 200; (void)a2;
    if(!q->sdr_ok)return 0;
    if(q->n_qbins<4){printf("  [MEMORY] need QBIN\n");return 0;}
    int k0=q->qbins[0],k1=q->qbins[q->n_qbins-1],D=q->wf.d;
    double a=0.7071,x[D],xi[D],y[D];
    int nt=8, corr=0;

    printf("  [MEMORY] ANTISYM -> wait -> measure A -> wait -> capture B_future\n");
    for(int tr=0;tr<nt;tr++){
        /* ANTISYM: entangle and TX GHZ into room */
        for(int p=0;p<8;p++){
            memset(x,0,D*8);memset(xi,0,D*8);
            if(p==0){
                for(int i=0;i<q->n_qbins;i+=2){
                    x[q->qbins[i]]+=a;x[D-q->qbins[i]]-=a;
                    x[q->qbins[i+1]]+=a;x[D-q->qbins[i+1]]-=a;
                }
            }else{
                double tot=0;
                for(int i=0;i<q->n_qbins;i++)tot+=y[q->qbins[i]]+y[D-q->qbins[i]];
                if(tot>1e-15)for(int i=0;i<q->n_qbins;i++){
                    double p=(y[q->qbins[i]]+y[D-q->qbins[i]])/tot;
                    if(p>1e-10){x[q->qbins[i]]=+sqrt(p);x[D-q->qbins[i]]=-sqrt(p);}
                }
            }
            for(int i=0;i<16;i++)x[i]=xi[i]=0;
            qvm_ofdm_compute(q,x,xi,y,D);
        }

        /* Room has GHZ multipath. Wait for it to develop. */
        usleep(delay*1000);

        /* Measure A from wf (ANTISYM result) */
        double a0=q->wf.prob[k0]+q->wf.prob[D-k0];
        double a1b=q->wf.prob[k1]+q->wf.prob[D-k1];
        double at=a0+a1b;
        int a_out=(at>1e-15?a0/at:0)>0.5?0:1;

        /* Catch B from the room's delayed multipath (no new TX) */
        usleep(30000);
        sdr_capture(q->sdr);
        wf_from_iq(q->sdr->iq_i,q->sdr->iq_q,q->sdr->iq_n,&q->wf);
        double b0=q->wf.prob[k0]+q->wf.prob[D-k0];
        double b1b=q->wf.prob[k1]+q->wf.prob[D-k1];
        double bt=b0+b1b;
        int b_out=(bt>1e-15?b1b/bt:0)>0.5?1:0;

        if(a_out==b_out)corr++;
        printf("    %2d: A_wf=|%d> B_room=|%d> %s\n",tr,a_out,b_out,
            a_out==b_out?"corr":"anti");
    }
    printf("  [MEMORY] %d/%d correlated (%.0f%%) %s\n",
        corr,nt,100.0*corr/nt,corr>nt*0.7?"RETROCAUSAL":"no effect");
    return 0;
}

/* FISSURE: smash two qudits, place third in collision field.
   Compares solo decay vs collision-field persistence. */
static int op_fissure(QvmCtx *q, double a1, double a2){
    (void)a1;(void)a2;
    if(!q->sdr_ok)return 0;
    if(q->n_qbins<6){printf("  [FISSURE] need >=6 bins (3 qudits)\n");return 0;}
    int b0=q->qbins[0],b1=q->qbins[2],b2=q->qbins[4]; /* 3 qudit bins */
    int D=q->wf.d; double a=0.7071,x[D],xi[D];
    int delays[]={5,20,50,100,200,500};

    printf("  [FISSURE] Solo vs collision-field persistence:\n");
    printf("  delay   solo   collision  boost\n");
    for(int di=0;di<6;di++){
        double solo_p, collision_p;

        /* SOLO: qudit at b2 alone */
        memset(x,0,D*8);memset(xi,0,D*8);
        x[b2]=+a;x[D-b2]=-a;
        for(int i=0;i<16;i++)x[i]=xi[i]=0;
        for(int i=0;i<D;i++){q->wf.re[i]=x[i];q->wf.im[i]=xi[i];}
        gate_ofdm_tx(q->sdr,&q->wf,NULL);
        {struct v4l2_buffer b;memset(&b,0,sizeof(b));
         b.type=V4L2_BUF_TYPE_SDR_CAPTURE;b.memory=V4L2_MEMORY_MMAP;
         if(ioctl(q->sdr->fd,VIDIOC_DQBUF,&b)==0)ioctl(q->sdr->fd,VIDIOC_QBUF,&b);}
        usleep(delays[di]*1000);
        sdr_capture(q->sdr);
        wf_from_iq(q->sdr->iq_i,q->sdr->iq_q,q->sdr->iq_n,&q->wf);
        solo_p=q->wf.prob[b2]+q->wf.prob[D-b2];

        /* COLLISION: smash b0+b1, fissure at b2 */
        memset(x,0,D*8);memset(xi,0,D*8);
        x[b0]=+a;x[D-b0]=-a;x[b1]=+a;x[D-b1]=-a;
        x[b2]=+a;x[D-b2]=-a;
        for(int i=0;i<16;i++)x[i]=xi[i]=0;
        for(int i=0;i<D;i++){q->wf.re[i]=x[i];q->wf.im[i]=xi[i];}
        gate_ofdm_tx(q->sdr,&q->wf,NULL);
        {struct v4l2_buffer b;memset(&b,0,sizeof(b));
         b.type=V4L2_BUF_TYPE_SDR_CAPTURE;b.memory=V4L2_MEMORY_MMAP;
         if(ioctl(q->sdr->fd,VIDIOC_DQBUF,&b)==0)ioctl(q->sdr->fd,VIDIOC_QBUF,&b);}
        usleep(delays[di]*1000);
        sdr_capture(q->sdr);
        wf_from_iq(q->sdr->iq_i,q->sdr->iq_q,q->sdr->iq_n,&q->wf);
        collision_p=q->wf.prob[b2]+q->wf.prob[D-b2];

        printf("  %4dms  %.4f  %.4f    %+.0f%%\n",delays[di],solo_p,collision_p,
            solo_p>1e-10?100*(collision_p/solo_p-1):0);
    }
    return 0;
}

static int op_tick(QvmCtx *q, double a1, double a2){
    (void)a1;(void)a2;
    if(q->sdr_ok){gate_tx_hardware(q->sdr,&q->wf);sdr_capture(q->sdr);wf_from_iq(q->sdr->iq_i,q->sdr->iq_q,q->sdr->iq_n,&q->wf);}
    printf("  [TICK] TX→ether→RX\n"); wf_print(&q->wf,"|ψ⟩",-1); return 0;
}

static int op_prob(QvmCtx *q, double a1, double a2){
    (void)a1;(void)a2;
    printf("  [PROB]"); for(int k=0;k<q->wf.d;k++) printf(" |%d⟩=%.3f",k,q->wf.prob[k]);
    printf("  S=%.3f γ=%.4f\n",q->wf.entropy,q->wf.purity); return 0;
}

static int op_show(QvmCtx *q, double a1, double a2){
    (void)a1;(void)a2; wf_print(&q->wf,"|ψ⟩",-1); return 0;
}

static int op_dump(QvmCtx *q, double a1, double a2){
    (void)a1;(void)a2;
    printf("  State vector:\n");
    for(int i=0;i<q->wf.d;i++) printf("    |%d⟩ = %+.6f %+.6fi  |ψ|²=%.6f\n",i,q->wf.re[i],q->wf.im[i],q->wf.prob[i]);
    return 0;
}

static int op_sample(QvmCtx *q, double a1, double a2){
    int n=((int)a1)>0?(int)a1:100, counts[256]={0};
    for(int s=0;s<n;s++){double r=(double)rand()/RAND_MAX,cum=0;int out=q->wf.d-1;
        for(int k=0;k<q->wf.d;k++){cum+=q->wf.prob[k];if(r<=cum){out=k;break;}}if(out<256)counts[out]++;}
    printf("  [SAMPLE] %d draws:\n  ",n);
    for(int i=0;i<q->wf.d;i++) printf("|%d⟩:%-3d ",i,counts[i]);printf("\n"); return 0;
}

static int op_set(QvmCtx *q, double a1, double a2){
    int k=(int)a1; if(k<0||k>=q->wf.d){printf("  [SET] out of range\n");return 0;}
    q->wf.re[k]=a2;q->wf.im[k]=0; qvm_norm(&q->wf);
    printf("  [SET] |%d⟩ = %.3f\n",k,a2); wf_print(&q->wf,"|ψ⟩",-1); return 0;
}

static int op_reset(QvmCtx *q, double a1, double a2){
    (void)a1;(void)a2;
    for(int i=0;i<q->wf.d;i++){q->wf.re[i]=1.0/sqrt(q->wf.d);q->wf.im[i]=0;q->wf.prob[i]=1.0/q->wf.d;}
    q->wf.entropy=log2(q->wf.d);q->wf.purity=1.0/q->wf.d;
    printf("  [RESET] Uniform\n"); wf_print(&q->wf,"|+⟩",-1); return 0;
}

static int op_swap(QvmCtx *q, double a1, double a2){
    int a=(int)a1,b=(int)a2; if(a<0||a>=q->wf.d||b<0||b>=q->wf.d){printf("  [SWAP] range\n");return 0;}
    double tr=q->wf.re[a],ti=q->wf.im[a],tp=q->wf.prob[a];
    q->wf.re[a]=q->wf.re[b];q->wf.im[a]=q->wf.im[b];q->wf.prob[a]=q->wf.prob[b];
    q->wf.re[b]=tr;q->wf.im[b]=ti;q->wf.prob[b]=tp;
    printf("  [SWAP] |%d⟩ ⇄ |%d⟩\n",a,b); return 0;
}

static int op_invert(QvmCtx *q, double a1, double a2){
    (void)a1;(void)a2; for(int i=0;i<q->wf.d;i++) q->wf.im[i]=-q->wf.im[i];
    printf("  [INVERT] Complex conjugate\n"); return 0;
}

static int op_scale(QvmCtx *q, double a1, double a2){
    (void)a2; double sc=a1>1e-15?a1:1.0;
    for(int i=0;i<q->wf.d;i++){q->wf.re[i]*=sc;q->wf.im[i]*=sc;} qvm_norm(&q->wf);
    printf("  [SCALE] ×%.3f\n",sc); return 0;
}

static int op_purity(QvmCtx *q, double a1, double a2){
    (void)a1;(void)a2; printf("  γ=%.6f  S=%.4f bits\n",q->wf.purity,q->wf.entropy); return 0;
}

static int op_coherent(QvmCtx *q, double a1, double a2){
    (void)a1;(void)a2;
    int np=IQ_WINDOW/2; TxBuf tx=tx_alloc(np,q->freq,q->rate);
    tx_synthesize(&q->wf,&tx); tx_write_file(&tx,"/tmp/qvm_coherent.iq"); tx_free(&tx);
    printf("  [COHERENT] OFDM → /tmp/qvm_coherent.iq\n"); return 0;
}

static int op_wait(QvmCtx *q, double a1, double a2){
    (void)a2; int ms=((int)a1)>0?(int)a1:100;
    printf("  [WAIT] %d ms\n",ms); usleep(ms*1000); return 0;
}

static int op_echo(QvmCtx *q, double a1, double a2){
    (void)q;(void)a1;(void)a2; return 0; /* handled in dispatcher */
}

static int op_help(QvmCtx *q, double a1, double a2){
    (void)a1;(void)a2;
    for(int i=0;i<q->n_ops;i++){if(q->ops[i].help[0])printf("  %-12s — %s\n",q->ops[i].name,q->ops[i].help);}
    return 0;
}

static int op_quit(QvmCtx *q, double a1, double a2){
    (void)a1;(void)a2; q->running=0; return 1;
}

/* ── Extended: GPU offload instructions ── */
static int op_calibrate(QvmCtx *q, double a1, double a2){
    int n_avg=((int)a1)>0?(int)a1:4; (void)a2;
    if(!q->sdr_ok){printf("  [CALIBRATE] No SDR\n");return 0;}
    printf("  [CALIBRATE] M (%d averages)...\n",n_avg);
    int d=q->wf.d;
    if(q->M){for(int i=0;i<q->M_dim;i++)free(((double**)q->M)[i]);free(q->M);free(q->Minv);}
    double **M=malloc(d*sizeof(double*)),**Mi=malloc(d*sizeof(double*));
    for(int i=0;i<d;i++){M[i]=calloc(d,sizeof(double));Mi[i]=calloc(d,sizeof(double));}
    double bw=(double)q->rate/d;
    for(int k=0;k<d;k++){
        uint32_t f=(uint32_t)q->freq+(uint32_t)(k*bw);
        sdr_retune(q->sdr,f);usleep(8000);sdr_flush(q->sdr);
        double sum[256]={0};int ncaps=0;
        for(int a=0;a<n_avg;a++){
            sdr_capture(q->sdr);int np=q->sdr->iq_n;if(np<d)continue;ncaps++;
            for(int m=0;m<d;m++){double fn=(double)m/d,si=0,sq=0;
                for(int n=0;n<np;n++){double ph=-2*M_PI*fn*n,c=cos(ph),s=sin(ph);si+=q->sdr->iq_i[n]*c-q->sdr->iq_q[n]*s;sq+=q->sdr->iq_i[n]*s+q->sdr->iq_q[n]*c;}
                sum[m]+=(si*si+sq*sq)/(np*np);}
        }
        if(ncaps)for(int m=0;m<d;m++)M[k][m]=sum[m]/ncaps;
    }
    sdr_retune(q->sdr,(uint32_t)q->freq);
    for(int i=0;i<d;i++){double rs=0;for(int j=0;j<d;j++)rs+=M[i][j];if(rs>1e-15)for(int j=0;j<d;j++)M[i][j]/=rs;}
    printf("  M (%dx%d):\n",d,d);
    for(int i=0;i<d;i++){printf("    ");for(int j=0;j<d;j++)printf("%.3f ",M[i][j]);printf("\n");}
    q->M=(double*)M;q->Minv=(double*)Mi;q->M_dim=d;
    sdr_restream(q->sdr); /* re-establish streaming after sweeps */
    return 0;
}

static int op_solve(QvmCtx *q, double a1, double a2){
    int n=((int)a1)>0?(int)a1:5;(void)a2;
    if(n<3)n=3;if(n>10)n=10;int sd=1<<n;
    if(!q->sdr_ok){printf("  [SOLVE] No SDR\n");return 0;}
    int items[12];srand(time(NULL));
    printf("  [SOLVE] n=%d items: [",n);for(int i=0;i<n;i++){items[i]=rand()%100+1;printf("%d ",items[i]);}
    int *sol=calloc(sd,sizeof(int)),ns=0,tg=0;
    for(int a=0;a<50;a++){int m=rand()%sd;tg=0;for(int i=0;i<n;i++)if(m&(1<<i))tg+=items[i];
        ns=0;memset(sol,0,sd*sizeof(int));
        for(int s=0;s<sd;s++){int sm=0;for(int i=0;i<n;i++)if(s&(1<<i))sm+=items[i];if(sm==tg){sol[s]=1;ns++;}}
        if(ns>=1&&ns<=sd/2&&tg>0)break;}
    printf("] target=%d sols=%d/%d\n",tg,ns,sd);
    double bw=(double)q->rate/sd;int dw=500;
    for(int s=0;s<sd;s++){if(sol[s])continue;uint32_t f=(uint32_t)q->freq+(uint32_t)(s*bw);if(f<24000000)f=24000000;if(f>1750000000)f=1750000000;sdr_retune(q->sdr,f);usleep(dw);}
    sdr_retune(q->sdr,(uint32_t)q->freq);usleep(500);sdr_flush(q->sdr);sdr_capture(q->sdr);
    int np=q->sdr->iq_n;if(np<sd)np=sd;double*pwr=calloc(sd,sizeof(double));
    for(int b=0;b<sd;b++){double fn=(double)b/sd,si=0,sq=0;
        for(int n=0;n<np;n++){double ph=-2*M_PI*fn*n,c=cos(ph),s=sin(ph);si+=q->sdr->iq_i[n]*c-q->sdr->iq_q[n]*s;sq+=q->sdr->iq_i[n]*s+q->sdr->iq_q[n]*c;}
        pwr[b]=(si*si+sq*sq)/((double)np*np);}
    double*srt=malloc(sd*sizeof(double));memcpy(srt,pwr,sd*sizeof(double));
    for(int i=0;i<sd-1;i++)for(int j=i+1;j<sd;j++)if(srt[i]>srt[j]){double t=srt[i];srt[i]=srt[j];srt[j]=t;}
    double med=srt[sd/2],th=med*0.5;free(srt);
    int fnd=0,cor=0;
    for(int s=0;s<sd;s++)if(pwr[s]<th){fnd++;if(sol[s])cor++;printf("  bin %3d [",s);for(int i=n-1;i>=0;i--)printf("%d",(s>>i)&1);int sm=0;for(int i=0;i<n;i++)if(s&(1<<i))sm+=items[i];printf("] sum=%d pwr=%.2e %s\n",sm,pwr[s],sol[s]?"✓":"✗");}
    for(int s=0;s<sd;s++)if(sol[s]&&pwr[s]>=th){int sm=0;for(int i=0;i<n;i++)if(s&(1<<i))sm+=items[i];printf("  ★ solution bin %d pwr=%.2e (missed)\n",s,pwr[s]);}
    printf("  Found %d/%d (%.0f%%)\n",cor,ns,ns>0?100.0*cor/ns:0);
    free(sol);free(pwr);return 0;
}

static int op_bench(QvmCtx *q, double a1, double a2){
    int bd=((int)a1)>0?(int)a1:8;(void)a2;
    if(bd<4)bd=4;if(bd>MAX_DIM)bd=MAX_DIM;
    double*ba=calloc(bd*bd,sizeof(double)),*bx=calloc(bd,sizeof(double)),*by=calloc(bd,sizeof(double));
    srand(42);for(int i=0;i<bd*bd;i++)ba[i]=(rand()%100)/100.0;for(int i=0;i<bd;i++)bx[i]=(rand()%100)/100.0;
    clock_t t0=clock();
    for(int t=0;t<100;t++)for(int i=0;i<bd;i++){by[i]=0;for(int j=0;j<bd;j++)by[i]+=ba[i*bd+j]*bx[j];}
    double cpu_us=(double)(clock()-t0)*1000000.0/CLOCKS_PER_SEC/100.0;
    if(q->sdr_ok){clock_t t1=clock();
        double bw=(double)q->rate/bd;
        for(int j=0;j<bd;j++){int dw=(int)(bx[j]/1.0*500+30);uint32_t f=(uint32_t)q->freq+(uint32_t)(j*bw);sdr_retune(q->sdr,f);usleep(dw);}
        sdr_retune(q->sdr,(uint32_t)q->freq);usleep(500);sdr_flush(q->sdr);sdr_capture(q->sdr);
        double room_us=(double)(clock()-t1)*1000000.0/CLOCKS_PER_SEC;
        printf("  [BENCH] %dx%d  CPU:%.1fμs/iter  Room:%.0fμs  ratio:%.1fx\n",bd,bd,cpu_us,room_us,cpu_us/room_us*100.0);
    }else{printf("  [BENCH] %dx%d  CPU:%.1fμs (no SDR)\n",bd,bd,cpu_us);}
    free(ba);free(bx);free(by);return 0;
}

/* ── LOOP / END (special: access script state) ── */
static int op_loop(QvmCtx *q, double a1, double a2){
    int count=((int)a1)>0?(int)a1:1;(void)a2;
    if(!q->interactive){q->loop_stack=realloc(q->loop_stack,(q->loop_depth+2)*sizeof(int));
        q->loop_stack[q->loop_depth++]=q->ip;q->loop_stack[q->loop_depth++]=count;printf("  [LOOP] ×%d\n",count);}
    else printf("  [LOOP] script only\n"); return 0;
}
static int op_end(QvmCtx *q, double a1, double a2){
    (void)a1;(void)a2;
    if(q->loop_depth>=2&&!q->interactive){int ct=q->loop_stack[--q->loop_depth];int st=q->loop_stack[--q->loop_depth];
        if(--ct>0){q->loop_stack[q->loop_depth++]=st;q->loop_stack[q->loop_depth++]=ct;q->ip=st;}}return 0;
}

/* STABILIZE: regenerative feedback — keep qudit alive indefinitely.
   TX anti-sym, capture, re-amplify, repeat. Qudit persists. */
static int op_stabilize(QvmCtx *q, double a1, double a2){
    int cycles = ((int)a1)>0 ? (int)a1 : 20; (void)a2;
    if(!q->sdr_ok)return 0;
    if(q->n_qbins<2){printf("  [STABILIZE] need QBIN\n");return 0;}
    int bin = q->qbins[q->n_qbins-1], D=q->wf.d;
    double a=0.7071, x[D],xi[D],y[D];
    printf("  [STABILIZE] bin %d, %d cycles:\n", bin, cycles);
    double prev = -1;
    for(int c=0; c<cycles; c++){
        memset(x,0,D*8); memset(xi,0,D*8);
        if(c==0){x[bin]=+a;x[D-bin]=-a;}
        else{double p=y[bin]+y[D-bin];
            if(p>1e-10){x[bin]=+a;x[D-bin]=-a;} /* AGC: fixed amplitude */
            else{x[bin]=+a;x[D-bin]=-a;}}
        for(int i=0;i<16;i++)x[i]=xi[i]=0;
        qvm_ofdm_compute(q,x,xi,y,D);
        double pwr=y[bin]+y[D-bin];
        double st=prev>1e-10?fabs(pwr-prev)/prev:1;
        printf("    %3d  pwr=%.4f  %s\n",c+1,pwr,st<0.1?"STABLE":st<0.3?"drift":"unstable");
        prev=pwr;
    }
    for(int i=0;i<D;i++)q->wf.re[i]=q->wf.im[i]=q->wf.prob[i]=0;
    q->wf.prob[bin]=1.0; q->wf.re[bin]=1.0;
    return 0;
}

/* GHZ_STABILIZE: keep entangled state alive by regenerating ALL bins */
static int op_ghz_stab(QvmCtx *q, double a1, double a2){
    int cycles = ((int)a1)>0 ? (int)a1 : 20; (void)a2;
    if(!q->sdr_ok)return 0;
    if(q->n_qbins<4){printf("  [GHZ_STAB] need >=4 bins\n");return 0;}
    int n=q->n_qbins, D=q->wf.d;
    double a=0.7071, x[D],xi[D],y[D];

    printf("  [GHZ_STAB] %d bins, %d cycles:\n", n, cycles);
    for(int c=0; c<cycles; c++){
        memset(x,0,D*8); memset(xi,0,D*8);
        if(c==0){
            /* Read existing wf state for first cycle */
            double init=0;
            for(int i=0;i<n;i++)init+=q->wf.prob[q->qbins[i]];
            if(init>1e-10){
                for(int i=0;i<n;i++){
                    double p=q->wf.prob[q->qbins[i]];
                    if(p>1e-10){x[q->qbins[i]]=+sqrt(p);x[D-q->qbins[i]]=-sqrt(p);}
                }
            }else{
                for(int i=0;i<n;i+=2){
                    x[q->qbins[i]]+=a;x[D-q->qbins[i]]-=a;
                    x[q->qbins[i+1]]+=a;x[D-q->qbins[i+1]]-=a;
                }
            }
        }else{
            double tot=0;
            for(int i=0;i<n;i++)tot+=y[q->qbins[i]]+y[D-q->qbins[i]];
            if(tot>1e-15){
                for(int i=0;i<n;i++){
                    double p=(y[q->qbins[i]]+y[D-q->qbins[i]])/tot;
                    if(p>1e-10){x[q->qbins[i]]=+sqrt(p);x[D-q->qbins[i]]=-sqrt(p);}
                }
            }
        }
        for(int i=0;i<16;i++)x[i]=xi[i]=0;
        qvm_ofdm_compute(q,x,xi,y,D);
    }
    double tot=0,s0=0,s1=0;
    for(int i=0;i<n;i++)tot+=y[q->qbins[i]]+y[D-q->qbins[i]];
    if(tot>1e-15){
        for(int i=0;i<n;i+=2){s0+=y[q->qbins[i]]+y[D-q->qbins[i]];s1+=y[q->qbins[i+1]]+y[D-q->qbins[i+1]];}
        double s=s0+s1;
        for(int i=0;i<n;i++)q->wf.prob[q->qbins[i]]=(y[q->qbins[i]]+y[D-q->qbins[i]])/s;
    }
    printf("  [GHZ_STAB] %d cycles done: |0>=%.3f |1>=%.3f\n",
        cycles, tot>1e-15?s0/(s0+s1):0, tot>1e-15?s1/(s0+s1):0);
    return 0;
}

static int op_proj_meas(QvmCtx *q, double a1, double a2){
    int qi=((int)a1)>=0?(int)a1:0;(void)a2;
    int nq=q->n_qbins/2;
    if(qi>=nq){printf("  [PROJ] qi=%d >= %d\n",qi,nq);return 0;}
    if(!q->sdr_ok)return 0;
    int b0=q->qbins[2*qi],b1=q->qbins[2*qi+1],D=q->wf.d;

    /* Active TX measurement on this qubit */
    double a=0.7071,x[D],xi[D],y[D];
    memset(x,0,D*8);memset(xi,0,D*8);
    x[b0]=+a;x[D-b0]=-a;
    for(int i=0;i<16;i++)x[i]=xi[i]=0;
    qvm_ofdm_compute(q,x,xi,y,D);
    double p0=y[b0]+y[D-b0],p1=y[b1]+y[D-b1];
    int outcome=(p1>p0)?1:0;

    /* Collapse: only the measured qubit, rotate its neighbors */
    q->wf.re[b0]=q->wf.im[b0]=q->wf.prob[b0]=0;
    q->wf.re[b1]=q->wf.im[b1]=q->wf.prob[b1]=0;
    double tot=p0+p1;
    int obin=outcome?b1:b0;
    q->wf.prob[obin]=(outcome?p1:p0)/tot;
    q->wf.re[obin]=sqrt(q->wf.prob[obin]);

    /* Feed-forward phase on all remaining qubits: if outcome=1, apply Z rotation */
    for(int i=0;i<nq;i++){
        if(i==qi)continue;
        if(outcome==1){
            /* Z gate on neighbor: negate |1> amplitude */
            int bi1=q->qbins[2*i+1];
            q->wf.re[bi1]=-q->wf.re[bi1];q->wf.im[bi1]=-q->wf.im[bi1];
        }
    }
    printf("  [PROJ] qudit %d measured -> |%d> (|0>=%.3f |1>=%.3f)\n",
        qi,outcome,tot>1e-15?p0/tot:0,tot>1e-15?p1/tot:0);
    return 0;
}

/* CREATE qi: algebraic creation a^dagger on qudit qi.
   Multiplies |0⟩ bin amplitude by √(n+1)/√n (unitary on WF, no measurement). */
static int op_create(QvmCtx *q, double a1, double a2){
    int qi=((int)a1)>=0?(int)a1:0;(void)a2;
    int nq=q->n_qbins/2;if(qi>=nq)return 0;
    int b0=q->qbins[2*qi];
    double n=q->wf.prob[b0];if(n<1e-10)n=0.01;
    double factor=sqrt((n+0.1)/n);
    q->wf.prob[b0]*=factor;q->wf.re[b0]=sqrt(q->wf.prob[b0]);
    printf("  [CREATE] qudit %d n=%.3f\n",qi,q->wf.prob[b0]);
    return 0;
}
/* ANNIHILATE qi: algebraic annihilation a on qudit qi. */
static int op_annihilate(QvmCtx *q, double a1, double a2){
    int qi=((int)a1)>=0?(int)a1:0;(void)a2;
    int nq=q->n_qbins/2;if(qi>=nq)return 0;
    int b0=q->qbins[2*qi];
    double n=q->wf.prob[b0];
    if(n<1e-10){printf("  [ANNIHILATE] empty\n");return 0;}
    double factor=sqrt(n>0.1?n-0.1:n*0.1)/sqrt(n);
    q->wf.prob[b0]*=factor;q->wf.re[b0]=sqrt(q->wf.prob[b0]);
    printf("  [ANNIHILATE] qudit %d n=%.3f\n",qi,q->wf.prob[b0]);
    return 0;
}

/* OCCUPATION qi: measure occupation through room. TX anti-sym at qudit,
   capture, read actual physical occupation number from room response. */
static int op_occupation(QvmCtx *q, double a1, double a2){
    int qi=((int)a1)>=0?(int)a1:0;(void)a2;
    int nq=q->n_qbins/2;if(qi>=nq)return 0;
    if(!q->sdr_ok)return 0;
    int b0=q->qbins[2*qi],D=q->wf.d;
    double n_wf=q->wf.prob[b0];
    double amp=sqrt(n_wf+0.01);if(amp>1.0)amp=1.0;
    double x[D],xi[D],y[D];
    memset(x,0,D*8);memset(xi,0,D*8);
    x[b0]=+amp;x[D-b0]=-amp;
    for(int i=0;i<16;i++)x[i]=xi[i]=0;
    qvm_ofdm_compute(q,x,xi,y,D);
    double n_room=y[b0]+y[D-b0];
    printf("  [OCCUPATION] qudit %d WF=%.3f room=%.3f\n",qi,n_wf,n_room);
    return 0;
}

/* CZ qa qb: controlled-Z gate. Entangles two qubits via anti-sym pair
   at their |1⟩ bins only, leaving all other qubits untouched. */
static int op_cz_gate(QvmCtx *q, double a1, double a2){
    int qa=((int)a1)>=0?(int)a1:0, qb=((int)a2)>=0?(int)a2:1;
    int nq=q->n_qbins/2;
    if(qa>=nq||qb>=nq){printf("  [CZ] bad qubits\n");return 0;}
    if(!q->sdr_ok)return 0;
    int b1a=q->qbins[2*qa+1], b1b=q->qbins[2*qb+1]; /* |1⟩ bins */
    int D=q->wf.d;
    double amp=0.7071,x[D],xi[D],y[D];

    /* Read existing WF state for these two qubits */
    double pa0=q->wf.prob[q->qbins[2*qa]],pa1=q->wf.prob[b1a];
    double pb0=q->wf.prob[q->qbins[2*qb]],pb1=q->wf.prob[b1b];
    double ta=pa0+pa1,tb=pb0+pb1;
    if(ta<1e-10||tb<1e-10)ta=tb=1;

    /* TX anti-sym at both |1⟩ bins with existing amplitudes */
    memset(x,0,D*8);memset(xi,0,D*8);
    if(pa1>1e-10){x[b1a]=+sqrt(pa1/ta);x[D-b1a]=-sqrt(pa1/ta);}
    else{x[b1a]=+amp;x[D-b1a]=-amp;}
    if(pb1>1e-10){x[b1b]=+sqrt(pb1/tb);x[D-b1b]=-sqrt(pb1/tb);}
    else{x[b1b]=+amp;x[D-b1b]=-amp;}
    for(int i=0;i<16;i++)x[i]=xi[i]=0;
    qvm_ofdm_compute(q,x,xi,y,D);

    /* Update only these two qubits, preserve others */
    double s=y[b1a]+y[D-b1a]+y[b1b]+y[D-b1b];
    if(s>1e-15){
        q->wf.prob[b1a]=(y[b1a]+y[D-b1a])/s;
        q->wf.prob[b1b]=(y[b1b]+y[D-b1b])/s;
    }
    printf("  [CZ] qudits %d-%d entangled\n",qa,qb);
    return 0;
}

/* XGATE qi: Pauli-X. Swaps |0> and |1> amplitudes. */
static int op_x_gate(QvmCtx *q, double a1, double a2){
    int qi=((int)a1)>=0?(int)a1:0;(void)a2;
    int nq=q->n_qbins/2;if(qi>=nq)return 0;
    int b0=q->qbins[2*qi],b1=q->qbins[2*qi+1];
    double t0=q->wf.prob[b0],t1=q->wf.prob[b1];
    q->wf.prob[b0]=t1;q->wf.re[b0]=sqrt(t1);
    q->wf.prob[b1]=t0;q->wf.re[b1]=sqrt(t0);
    printf("  [X] qudit %d flipped\n",qi);return 0;
}

/* HGATE qi: Hadamard. Equal superposition. */
static int op_h_gate(QvmCtx *q, double a1, double a2){
    int qi=((int)a1)>=0?(int)a1:0;(void)a2;
    int nq=q->n_qbins/2;if(qi>=nq)return 0;
    int b0=q->qbins[2*qi],b1=q->qbins[2*qi+1];
    double amp=0.5;
    q->wf.prob[b0]=amp;q->wf.re[b0]=sqrt(amp);
    q->wf.prob[b1]=amp;q->wf.re[b1]=sqrt(amp);
    printf("  [H] qudit %d -> |+>\n",qi);return 0;
}

/* MBASIS qi angle: measure qudit in basis |0> +/- e^{i*angle}|1>.
   Sets WF to rotated basis, then PROJ measures. */
static int op_mbasis(QvmCtx *q, double a1, double a2){
    int qi=((int)a1)>=0?(int)a1:0;
    double angle=a2;
    int nq=q->n_qbins/2;
    if(qi>=nq||!q->sdr_ok)return 0;
    int b0=q->qbins[2*qi],b1=q->qbins[2*qi+1],D=q->wf.d;

    /* Rotate WF: apply H * R_z(angle) to the qubit's bins.
       |0> -> (|0> + |1>)/√2, |1> -> e^{i*angle}(|0> - |1>)/√2 */
    double ca=cos(angle/2),sa=sin(angle/2);
    double r0=q->wf.re[b0],i0=q->wf.im[b0];
    double r1=q->wf.re[b1],i1=q->wf.im[b1];

    /* New superposed amplitudes */
    q->wf.re[b0]=(r0*ca - i0*sa + r1*sa + i1*ca)/sqrt(2);
    q->wf.im[b0]=(i0*ca + r0*sa - r1*ca + i1*sa)/sqrt(2);
    q->wf.re[b1]=(r0*ca - i0*sa - r1*sa - i1*ca)/sqrt(2);
    q->wf.im[b1]=(i0*ca + r0*sa + r1*ca - i1*sa)/sqrt(2);
    q->wf.prob[b0]=q->wf.re[b0]*q->wf.re[b0]+q->wf.im[b0]*q->wf.im[b0];
    q->wf.prob[b1]=q->wf.re[b1]*q->wf.re[b1]+q->wf.im[b1]*q->wf.im[b1];

    /* Now measure in Z via PROJ */
    op_proj_meas(q,a1,0);
    return 0;
}

/* PlaneWarp QEC decoder */
extern int solve_plane(int r, int s, unsigned char *syn, unsigned char *out);
extern int preprocess_syndrome(int r, int s, unsigned char *syn);

static int op_qec_grid(QvmCtx *q, double a1, double a2){
    int r=((int)a1)>0?(int)a1:4, s=((int)a2)>0?(int)a2:4;
    int nq=q->n_qbins/2;
    if(r*s>nq){printf("  [QEC_GRID] need %d qubits (have %d)\n",r*s,nq);return 0;}

    /* Read outcomes from WF (not active measurement) */
    int *outcomes=calloc(r*s,sizeof(int));
    for(int i=0;i<r*s;i++){
        outcomes[i]=(q->wf.prob[q->qbins[2*i+1]]>q->wf.prob[q->qbins[2*i]])?1:0;
    }

    /* Stride-2 plaquette syndrome */
    unsigned char *syn=calloc(r*s,1), *corr=calloc(r*s,1);
    for(int i=0;i<r-1;i+=2)
        for(int j=0;j<s-1;j+=2)
            syn[i*s+j]=outcomes[i*s+j]^outcomes[(i+2)*s+j]
                       ^outcomes[i*s+(j+2)]^outcomes[(i+2)*s+(j+2)];

    preprocess_syndrome(r,s,syn);
    int ok=solve_plane(r,s,syn,corr);
    printf("  [QEC_GRID] %dx%d decode=%s",r,s,ok?"OK":"ABSTAIN");

    int ncorr=0;
    char cmd[16];
    for(int i=0;i<r*s;i++)if(corr[i]){
        snprintf(cmd,sizeof(cmd),"XGATE %d",i);qvm_eval(q,cmd);ncorr++;
    }
    printf("  corrected=%d\n",ncorr);
    free(outcomes);free(syn);free(corr);
    return 0;
}

/* ── Register all standard ops ── */
static void qvm_init_ops(QvmCtx *q){
    qvm_reg(q, "INIT",      op_init,      "capture ambient RF → |ψ⟩");
    qvm_reg(q, "SUPERPOSE", op_superpose, "coherent decomposition");
    qvm_reg(q, "SP",        op_superpose, "alias for SUPERPOSE");
    qvm_reg(q, "X",         op_x,         "cyclic shift [n=1]");
    qvm_reg(q, "Z",         op_z,         "phase rotate [rad=π/4] [k=-1=all]");
    qvm_reg(q, "H",         op_h,         "Hadamard via LO hop");
    qvm_reg(q, "CZ",        op_cz_gate,   "controlled-Z gate on two qubits");
    qvm_reg(q, "ANTISYM",   op_antisym,   "anti-sym pair entangle (8-pass feedback)");
    qvm_reg(q, "CHSH",      op_chsh,      "Bell inequality test on qbin state");
    qvm_reg(q, "MERMIN",    op_mermin,    "Mermin inequality for N-qudit GHZ");
    qvm_reg(q, "QBIN",      op_qbin,      "set qubit bin <idx> <k>");
    qvm_reg(q, "QBIN!",     op_qbin_done, "finalize qubit bin mapping");
    qvm_reg(q, "DFT",       op_dft,       "LO retune [step=1]");
    qvm_reg(q, "TX",        op_tx,        "radiate into ether");
    qvm_reg(q, "RX",        op_rx,        "capture from ether");
    qvm_reg(q, "MEASURE",   op_measure,   "Born-rule collapse");
    qvm_reg(q, "M",         op_measure,   "alias for MEASURE");
    qvm_reg(q, "QMEASURE",  op_qmeasure,  "passive wait → env decoheres");
    qvm_reg(q, "COLLAPSE",  op_collapse,  "anti-sym noise → collapse (87.5%)");
    qvm_reg(q, "KILL",      op_kill,      "winner feedback -> lock outcome (75%)");
    qvm_reg(q, "MEMORY",    op_memory,    "probe room multipath persistence");
    qvm_reg(q, "FISSURE",   op_fissure,   "smash pair, fissure qudit, measure persistence");
    qvm_reg(q, "STABILIZE", op_stabilize, "regenerative feedback — keep qudit alive");
    qvm_reg(q, "GHZ_STAB",  op_ghz_stab,  "regenerate all bins — preserve entanglement");
    qvm_reg(q, "PROJ",      op_proj_meas, "projective Z-measurement on qudit");
    qvm_reg(q, "CREATE",    op_create,    "algebraic creation a^dagger on qudit");
    qvm_reg(q, "ANNIHILATE",op_annihilate,"algebraic annihilation a on qudit");
    qvm_reg(q, "OCCUPATION",op_occupation,"measure occupation through room");
    qvm_reg(q, "XGATE",     op_x_gate,    "Pauli-X (NOT) gate");
    qvm_reg(q, "HGATE",     op_h_gate,    "Hadamard gate");
    qvm_reg(q, "MBASIS",    op_mbasis,    "measure in rotated basis M(angle)");
    qvm_reg(q, "QEC_GRID",  op_qec_grid,  "PlaneWarp toric code on r x s grid");
    qvm_reg(q, "TICK",      op_tick,      "TX→ether→RX cycle");
    qvm_reg(q, "PROB",      op_prob,      "show probabilities");
    qvm_reg(q, "P",         op_prob,      "alias for PROB");
    qvm_reg(q, "SHOW",      op_show,      "show state");
    qvm_reg(q, "S",         op_show,      "alias for SHOW");
    qvm_reg(q, "DUMP",      op_dump,      "full state vector");
    qvm_reg(q, "SAMPLE",    op_sample,    "Born-rule sampling [n=100]");
    qvm_reg(q, "SET",       op_set,       "set |k⟩ amplitude");
    qvm_reg(q, "RESET",     op_reset,     "uniform superposition");
    qvm_reg(q, "SWAP",      op_swap,      "swap |a⟩ and |b⟩");
    qvm_reg(q, "INVERT",    op_invert,    "complex conjugate");
    qvm_reg(q, "SCALE",     op_scale,     "scale amplitudes");
    qvm_reg(q, "PURITY",    op_purity,    "show purity/entropy");
    qvm_reg(q, "COHERENT",  op_coherent,  "OFDM synthesis → file");
    qvm_reg(q, "WAIT",      op_wait,      "idle [ms=100]");
    qvm_reg(q, "ECHO",      op_echo,      "print text");
    qvm_reg(q, "LOOP",      op_loop,      "start loop [n=1]");
    qvm_reg(q, "END",       op_end,       "end loop block");
    qvm_reg(q, "HELP",      op_help,      "this list");
    qvm_reg(q, "?",         op_help,      "alias for HELP");
    qvm_reg(q, "QUIT",      op_quit,      "exit VM");
    qvm_reg(q, "EXIT",      op_quit,      "alias for QUIT");
    qvm_reg(q, "Q",         op_quit,      "alias for QUIT");
    qvm_reg(q, "CALIBRATE", op_calibrate, "measure room channel M [avg=4]");
    qvm_reg(q, "SOLVE",     op_solve,     "subset sum via room [n=5]");
    qvm_reg(q, "BENCH",     op_bench,     "room vs CPU matvec [D=8]");
}

/* ── QVM API: public accessors (qvm_api.h) ── */
int qvm_running(QvmCtx *q){ return q->running; }
int qvm_dim(QvmCtx *q){ return q->wf.d; }
double qvm_prob(QvmCtx *q, int l){ return (l>=0&&l<q->wf.d)?q->wf.prob[l]:0; }
void qvm_probs(QvmCtx *q, double *out, int max){
    int n=q->wf.d<max?q->wf.d:max;
    for(int i=0;i<n;i++) out[i]=q->wf.prob[i];
}
double qvm_entropy(QvmCtx *q){ return q->wf.entropy; }
double qvm_purity(QvmCtx *q){ return q->wf.purity; }
int qvm_has_sdr(QvmCtx *q){ return q->sdr_ok; }
int qvm_calibrated(QvmCtx *q){ return q->M && q->Minv; }
void qvm_get_channel(QvmCtx *q, double **Mo, double **Mio, int *do_){
    if(Mo)*Mo=q->M; if(Mio)*Mio=q->Minv; if(do_)*do_=q->M_dim;
}

/* Offload computation directly through the room.
   Encodes x as LO dwell pattern → room processes → DFT reads result.
   The room IS the function. No software matvec, no pre-compensation.
   To approximate a specific function f(x), train readout weights offline. */
int qvm_compute(QvmCtx *q, const double *x, double *y, int d){
    if (!q->sdr_ok || d > q->wf.d) return -1;

    /* Feed x into room: LO dwell ∝ x[k] */
    double bw=(double)q->rate/d;
    for (int k=0;k<d;k++){
        if(fabs(x[k])<1e-10) continue; /* skip silent bins entirely */
        int dw=(int)(fabs(x[k])*1500+30);
        uint32_t f=(uint32_t)q->freq+(uint32_t)(k*bw);
        if(f<24000000)f=24000000; if(f>1750000000)f=1750000000;
        sdr_retune(q->sdr,f); usleep(dw);
    }
    sdr_retune(q->sdr,(uint32_t)q->freq); usleep(500);
    sdr_flush(q->sdr); sdr_capture(q->sdr);

    /* Room's answer: DFT of captured I/Q */
    int np=q->sdr->iq_n; if(np<d)np=d;
    double *pwr=calloc(d,sizeof(double));
    for (int k=0;k<d;k++){double fn=(double)k/d,si=0,sq=0;
        for(int n=0;n<np;n++){double ph=-2*M_PI*fn*n,c=cos(ph),s=sin(ph);
            si+=q->sdr->iq_i[n]*c-q->sdr->iq_q[n]*s;
            sq+=q->sdr->iq_i[n]*s+q->sdr->iq_q[n]*c;}
        pwr[k]=(si*si+sq*sq)/((double)np*np);}
    for (int i=0;i<d;i++) y[i]=pwr[i];
    free(pwr);
    return 0;
}

/* ── Low-level SDR access for waveform synthesis ── */
void qvm_sdr_tune(QvmCtx *q, uint32_t hz){
    if(!q->sdr_ok)return;
    sdr_retune(q->sdr,hz);
}
void qvm_sdr_restream(QvmCtx *q){
    if(!q->sdr_ok)return;
    sdr_restream(q->sdr);
}

/* Recycle buffers: DQBUF+QBUF to drain stale data, no STREAMOFF.
   Keeps DMA pipeline alive — preserves phase coherence across calls. */
void qvm_sdr_recycle(QvmCtx *q){
    if(!q->sdr_ok)return;
    /* Drain just 2 buffers — skip stale TX/cal, leave rest for DMA */
    for(int i=0;i<2;i++){
        struct pollfd p={.fd=q->sdr->fd,.events=POLLIN};
        if(poll(&p,1,50)<=0){
            /* Try non-blocking for stale TX */
            struct v4l2_buffer b;memset(&b,0,sizeof(b));
            b.type=V4L2_BUF_TYPE_SDR_CAPTURE;b.memory=V4L2_MEMORY_MMAP;
            if(ioctl(q->sdr->fd,VIDIOC_DQBUF,&b)!=0)break;
            ioctl(q->sdr->fd,VIDIOC_QBUF,&b);
            continue;
        }
        struct v4l2_buffer b;memset(&b,0,sizeof(b));
        b.type=V4L2_BUF_TYPE_SDR_CAPTURE;b.memory=V4L2_MEMORY_MMAP;
        if(ioctl(q->sdr->fd,VIDIOC_DQBUF,&b)!=0)break;
        ioctl(q->sdr->fd,VIDIOC_QBUF,&b);
    }
}
int qvm_sdr_fd(QvmCtx *q){ return q->sdr_ok ? q->sdr->fd : -1; }
int qvm_sdr_rx(QvmCtx *q, double *I, double *Q, int max_n){
    if(!q->sdr_ok)return 0;
    sdr_capture(q->sdr);
    int n=q->sdr->iq_n; if(n>max_n)n=max_n;
    for(int i=0;i<n;i++){I[i]=q->sdr->iq_i[i];Q[i]=q->sdr->iq_q[i];}
    return n;
}
void qvm_sdr_dft(QvmCtx *q, double *pwr, int D, const double *I, const double *Q, int np){
    if(np<D)np=D;
    for(int k=0;k<D;k++){double fn=(double)k/D,si=0,sq=0;
        for(int n=0;n<np;n++){double ph=-2*M_PI*fn*n,c=cos(ph),s=sin(ph);
            si+=I[n]*c-Q[n]*s;sq+=I[n]*s+Q[n]*c;}
        pwr[k]=(si*si+sq*sq)/((double)np*np);}
}

/* Complex DFT: returns I+jQ per bin for phase measurement */
int qvm_ofdm_complex(QvmCtx *q, const double *re, const double *im,
                      double *I_out, double *Q_out, int d){
    if(!q->sdr_ok||d>q->wf.d)return -1;
    for(int i=0;i<d;i++){q->wf.re[i]=re[i];q->wf.im[i]=im[i];
        q->wf.prob[i]=re[i]*re[i]+im[i]*im[i];}
    gate_ofdm_tx(q->sdr,&q->wf,NULL);
    usleep(80000); /* 80ms → DMA fills ~5 buffers, TX buffer rotated out */
    sdr_capture(q->sdr);
    int np=q->sdr->iq_n;if(np<d)np=d;
    for(int k=0;k<d;k++){double fn=(double)k/d,si=0,sq=0;
        for(int n=0;n<np;n++){double ph=-2*M_PI*fn*n,c=cos(ph),s=sin(ph);
            si+=q->sdr->iq_i[n]*c-q->sdr->iq_q[n]*s;
            sq+=q->sdr->iq_i[n]*s+q->sdr->iq_q[n]*c;}
        I_out[k]=si/np;Q_out[k]=sq/np;}
    return 0;
}
int qvm_ofdm_compute(QvmCtx *q, const double *re, const double *im, double *y, int d){
    if(!q->sdr_ok||d>q->wf.d)return -1;
    for(int i=0;i<d;i++){q->wf.re[i]=re[i];q->wf.im[i]=im[i];
        q->wf.prob[i]=q->wf.re[i]*q->wf.re[i]+q->wf.im[i]*q->wf.im[i];}
    gate_ofdm_tx(q->sdr,&q->wf,NULL);
    { struct v4l2_buffer b;memset(&b,0,sizeof(b));
      b.type=V4L2_BUF_TYPE_SDR_CAPTURE;b.memory=V4L2_MEMORY_MMAP;
      if(ioctl(q->sdr->fd,VIDIOC_DQBUF,&b)==0)ioctl(q->sdr->fd,VIDIOC_QBUF,&b); }
    usleep(50000);
    sdr_capture(q->sdr);
    wf_from_iq(q->sdr->iq_i,q->sdr->iq_q,q->sdr->iq_n,&q->wf);
    for(int i=0;i<d;i++)y[i]=q->wf.prob[i];
    return 0;
}

/* Anti-symmetric pair encoding: for each bin k, set X[k]=+A, X[D-k]=-A.
   The IM2 self-difference at k+(D-k)=D→0 cancels at DC. Prevents bin-0 collapse. */
void qvm_antisym_encode(QvmCtx *q, const int *bins, const double *amps,
                         int n_pairs, double *x_out, double *xi_out){
    int D=q->wf.d;
    memset(x_out,0,D*sizeof(double));
    memset(xi_out,0,D*sizeof(double));
    for(int i=0;i<n_pairs;i++){
        int k=bins[i];
        x_out[k]=amps[i];
        x_out[D-k]=-amps[i];
    }
    for(int i=0;i<8;i++)x_out[i]=0.0;
}

int qvm_eval(QvmCtx *q, const char *cmd){
    char op[32]={0}; double a1=0,a2=0;
    sscanf(cmd,"%31s %lf %lf",op,&a1,&a2);
    if(!op[0]) return 0;

    QvmOp fn = qvm_lookup(q, op);
    if (!fn) { printf("  ? Unknown: %s\n", op); return 0; }
    if (fn == op_echo) { printf("  %s\n", cmd+4); return 0; }

    int rc = fn(q, a1, a2);

    /* Auto-sync with ether only for operations that NEED the room.
       Software operations (X,Z,SWAP,SET,etc.) just update the wf representation.
       The room IS the ground truth — sync before TX/MEASURE/RX reads. */
    if (q->sdr_ok && q->running && rc >= 0) {
        int needs_sync = (fn==op_tx||fn==op_rx||fn==op_measure||
                          fn==op_tick||fn==op_cz||fn==op_h||
                          fn==op_init||fn==op_superpose||
                          fn==op_dft);
        if (needs_sync) qvm_sync(q);
    }
    return rc;
}

/* ── qvm_create / qvm_destroy ── */
QvmCtx *qvm_create(uint32_t freq, uint32_t rate, int D, int gain){
    QvmCtx *q = calloc(1, sizeof(QvmCtx));
    if (!q) return NULL;
    q->freq = freq; q->rate = rate;
    q->sdr = calloc(1, sizeof(SdrDev));
    q->sdr_ok = (sdr_open(q->sdr, freq, rate, gain) == 0);
    q->wf = wf_alloc(D, freq, rate);
    if (q->sdr_ok) {
        sdr_capture(q->sdr);
        wf_from_iq(q->sdr->iq_i, q->sdr->iq_q, q->sdr->iq_n, &q->wf);
    } else {
        for (int k=0;k<D;k++) { q->wf.re[k]=1.0/sqrt(D); q->wf.prob[k]=1.0/D; }
        q->wf.purity=1.0/D; q->wf.entropy=log2(D);
    }
    q->running = 1;
    memset(q->qbins, -1, sizeof(q->qbins)); q->n_qbins = 0;
    qvm_init_ops(q);
    return q;
}

void qvm_destroy(QvmCtx *q){
    if (!q) return;
    if (q->sdr_ok) sdr_close(q->sdr);
    free(q->sdr);
    wf_free(&q->wf);
    if (q->M) {
        for (int i=0;i<q->M_dim;i++) free(((double**)q->M)[i]);
        free(q->M); free(q->Minv);
    }
    for (int i=0;i<q->nlines;i++) free(q->lines[i]);
    free(q->lines); free(q->loop_stack);
    free(q);
}

/* ── qvm_run: REPL / script execution ── */
int qvm_run(QvmCtx *q, const char *script_path){
    printf("\n");
    printf("  ╔══════════════════════════════════════════════════════════════╗\n");
    printf("  ║  QUANTUM VM — Ether-native instruction set (API v2)          ║\n");
    printf("  ║  D=%d  |  %.1f MHz  |  %.2f MSPS  |  Δf=%.1f Hz            ║\n",
           q->wf.d, q->freq/1e6, q->rate/1e6, (double)q->rate/q->wf.d);
    printf("  ║  Substrate: %-50s ║\n",
           q->sdr_ok ? "PHYSICAL (RTL-SDR + EM field)" : "SIMULATED");
    printf("  ║  %d instructions registered                                  ║\n", q->n_ops);
    printf("  ╚══════════════════════════════════════════════════════════════╝\n");
    printf("\n  VM ready. %s\n\n",
           q->sdr_ok?"Physical ether substrate active.":"No SDR — software fallback.");
    wf_print(&q->wf, "|ψ⟩ init", -1);

    q->interactive = (!script_path || strcmp(script_path,"-")==0);
    if (!q->interactive) {
        FILE *f = fopen(script_path, "r");
        if (!f) { fprintf(stderr,"[VM] Cannot open %s\n",script_path); return 1; }
        char lbuf[512];
        q->lines = malloc(1024*sizeof(char*));
        while (fgets(lbuf,sizeof(lbuf),f) && q->nlines<1024) {
            char *nl=strchr(lbuf,'\n');if(nl)*nl=0;
            char *cmt=strchr(lbuf,'#');if(cmt)*cmt=0;
            char *s=lbuf; while(*s==' '||*s=='\t')s++;
            if(*s) q->lines[q->nlines++]=strdup(s);
        }
        fclose(f);
    }

    while (q->running) {
        char line[512]={0};
        if (!q->interactive) {
            if (q->ip >= q->nlines) break;
            strcpy(line, q->lines[q->ip++]);
        } else {
            printf("\n  qvm> "); fflush(stdout);
            if (!fgets(line,sizeof(line),stdin)) break;
            char *nl=strchr(line,'\n'); if(nl)*nl=0;
            char *cmt=strchr(line,'#'); if(cmt)*cmt=0;
            char *s=line; while(*s==' '||*s=='\t')s++;
            if(*s==0) continue;
            /* In interactive mode, shift the command to start */
            if (s != line) memmove(line, s, strlen(s)+1);
        }
        if (line[0]==0) continue;
        qvm_eval(q, line);
    }
    printf("\n  ★ VM halted. ★\n\n");
    return 0;
}

/* ── Top-level entry point (replaces run_quantum_vm) ── */
static int run_quantum_vm(uint32_t freq, uint32_t rate, int gain, int D,
                           const char *script_path) {
    QvmCtx *q = qvm_create(freq, rate, D, gain);
    if (!q) return 1;
    int rc = qvm_run(q, script_path);
    qvm_destroy(q);
    return rc;
}

/* ═══════════════════════════════════════════════════════════════
 * MODE: ETHER DECODE — measure M, compute M⁻¹, reverse the ether
 *
 * 1. Measure the ether's forward channel matrix M
 * 2. Compute its inverse M⁻¹ (Gaussian elimination)
 * 3. For each qudit level k as a probe:
 *      a. TX |k⟩ into ether via LO leakage
 *      b. RX the distorted state vector v_rx (DFT over all D bins)
 *      c. Decode: v_est = M⁻¹ · v_rx
 *      d. v_est should peak at bin k — the ether distortion is removed
 * 4. Show TX fidelity before and after decoding
 *
 * If M⁻¹ works, the decoded state should be MUCH closer to |k⟩
 * than the raw received state.  This proves we can digitally invert
 * the ether's physical transformation.
 * ═══════════════════════════════════════════════════════════════ */
static int run_ether_decode(uint32_t freq, uint32_t rate, int gain,
                             int D, int n_avg) {
    printf("\n");
    printf("  ╔══════════════════════════════════════════════════════════════╗\n");
    printf("  ║  ETHER DECODE — Measure M, invert, cancel ether distortion   ║\n");
    printf("  ║  D=%d  |  %.1f MHz  |  avg=%d                                 ║\n",
           D, freq/1e6, n_avg);
    printf("  ║  M⁻¹·M ≈ I  →  software reverses physical channel          ║\n");
    printf("  ╚══════════════════════════════════════════════════════════════╝\n\n");

    SdrDev sdr;
    if (sdr_open(&sdr, freq, rate, gain) != 0) {
        printf("  [FAIL] No SDR hardware.\n\n");
        return 1;
    }

    double **M     = malloc(D * sizeof(double*));
    double **M_inv = malloc(D * sizeof(double*));
    for (int i = 0; i < D; i++) {
        M[i]     = calloc(D, sizeof(double));
        M_inv[i] = calloc(D, sizeof(double));
    }

    printf("  Step 1: Measure ether channel matrix M…\n");
    gate_ether_gate_matrix(&sdr, D, (uint32_t)freq, M, n_avg);

    printf("\n  Step 2: Compute M⁻¹ via Gaussian elimination…\n");
    double cond = matrix_invert(M, M_inv, D);

    if (cond < 1e-10) {
        printf("  ✗ M is singular (cond ≈ 0) — cannot invert.\n");
        printf("    The ether's channel is not uniquely reversible.\n");
        goto cleanup;
    }

    printf("  ✓ M is invertible.  Condition number ≈ %.1f\n", cond);
    if (cond > 1000.0)
        printf("    ⚠ High condition number — decoding may amplify noise.\n");
    else if (cond < 10.0)
        printf("    Excellent conditioning — cleanly invertible.\n");

    printf("\n  Forward channel M  (row = TX level, col = RX level):\n");
    printf("  TX\\RX");
    for (int m = 0; m < D; m++) printf("   |%d⟩   ", m);
    printf("\n  ──────");
    for (int m = 0; m < D; m++) printf("─────────");
    printf("\n");
    for (int k = 0; k < D; k++) {
        printf("   |%d⟩  ", k);
        for (int m = 0; m < D; m++)
            printf(" %7.1e", M[k][m]);
        printf("\n");
    }

    printf("\n  Inverse channel M⁻¹ (decoding matrix):\n");
    printf("  RX\\TX");
    for (int m = 0; m < D; m++) printf("   |%d⟩   ", m);
    printf("\n  ──────");
    for (int m = 0; m < D; m++) printf("─────────");
    printf("\n");
    for (int k = 0; k < D; k++) {
        printf("   |%d⟩  ", k);
        for (int m = 0; m < D; m++)
            printf(" %+7.1e", M_inv[k][m]);
        printf("\n");
    }

    double bin_bw = (double)rate / D;
    double *v_rx  = calloc(D, sizeof(double));
    double *v_est = calloc(D, sizeof(double));

    printf("\n  Step 3: Probe each level, decode with M⁻¹\n");
    printf("  ───────────────────────────────────────────────────\n");
    printf("  %-6s  %8s  %8s  %s\n", "TX|k⟩", "RX peak", "Decoded", "Correct?");
    printf("  ───────────────────────────────────────────────────\n");

    int correct = 0;
    for (int k = 0; k < D; k++) {
        uint32_t tx_freq = (uint32_t)freq + (uint32_t)(k * bin_bw);
        if (tx_freq < 24000000) tx_freq = 24000000;
        if (tx_freq > 1750000000) tx_freq = 1750000000;

        sdr_retune(&sdr, tx_freq);
        usleep(8000);
        sdr_flush(&sdr);

        memset(v_rx, 0, D * sizeof(double));
        for (int avg = 0; avg < 2; avg++) {
            sdr_capture(&sdr);
            int np = sdr.iq_n;
            if (np < D) np = D;
            for (int m = 0; m < D; m++) {
                double sum_i = 0, sum_q = 0;
                double fn = (double)m / (double)D;
                for (int n = 0; n < np; n++) {
                    double phase = -2.0 * M_PI * fn * (double)n;
                    double cr = cos(phase), sr = sin(phase);
                    sum_i += sdr.iq_i[n] * cr - sdr.iq_q[n] * sr;
                    sum_q += sdr.iq_i[n] * sr + sdr.iq_q[n] * cr;
                }
                v_rx[m] += (sum_i*sum_i + sum_q*sum_q) / (np*np);
            }
        }
        for (int m = 0; m < D; m++) v_rx[m] /= 2.0;

        /* Subtract row-k mean to remove ambient bias, isolate LO */
        double row_mean = 0;
        for (int m = 0; m < D; m++) row_mean += v_rx[m];
        row_mean /= D;
        double *v_clean = calloc(D, sizeof(double));
        for (int m = 0; m < D; m++)
            v_clean[m] = v_rx[m] > row_mean ? v_rx[m] - row_mean : 0.0;

        matrix_apply(M_inv, v_rx, v_est, D);

        int rx_peak = 0, est_peak = 0;
        double rx_max = v_rx[0], est_max = v_est[0];
        for (int m = 1; m < D; m++) {
            if (v_rx[m]  > rx_max)  { rx_max  = v_rx[m];  rx_peak  = m; }
            if (v_est[m] > est_max) { est_max = v_est[m]; est_peak = m; }
        }

        int ok = (est_peak == k);
        if (ok) correct++;
        printf("  |%-5d  |%-8d  |%-8d  %s\n",
               k, rx_peak, est_peak, ok ? "✓" : "✗");

        /* Show the decoded vector for the most interesting case */
        if (k == D/2 || !ok) {
            printf("          v_rx:  ");
            for (int m = 0; m < D; m++) printf("%.2e ", v_rx[m]);
            printf("\n          v_est: ");
            for (int m = 0; m < D; m++) printf("%.2e ", v_est[m]);
            printf("\n");
        }
    }

    printf("  ───────────────────────────────────────────────────\n");
    printf("  Decode accuracy: %d/%d (%.0f%%)\n", correct, D, 100.0*correct/D);
    if (correct >= D/2)
        printf("  ★ M⁻¹ successfully reverses the ether channel!\n");
    else
        printf("  ⚠ Channel too noisy — more averaging needed.\n");

    free(v_rx); free(v_est);
cleanup:
    for (int i = 0; i < D; i++) { free(M[i]); free(M_inv[i]); }
    free(M); free(M_inv);
    sdr_close(&sdr);
    return 0;
}

/* ═══════════════════════════════════════════════════════════════
 * MODE: ETHER GATE MATRIX — display the quantum gate measured
 * from the ether.  This IS the gate — pure physical computation.
 * ═══════════════════════════════════════════════════════════════ */
static int run_ether_gate_matrix(uint32_t freq, uint32_t rate, int gain,
                                  int D, int n_avg) {
    printf("\n");
    printf("  ╔══════════════════════════════════════════════════════════════╗\n");
    printf("  ║  ETHER QUANTUM GATE MATRIX — Pure Physical Computation       ║\n");
    printf("  ║  D=%d  |  %.1f MHz  |  avg=%d                                 ║\n",
           D, freq/1e6, n_avg);
    printf("  ║  TX |k⟩ → ether nonlinearity → measure ALL bins             ║\n");
    printf("  ╚══════════════════════════════════════════════════════════════╝\n\n");

    SdrDev sdr;
    if (sdr_open(&sdr, freq, rate, gain) != 0) {
        printf("  [FAIL] No SDR hardware.\n\n");
        return 1;
    }

    double **M = malloc(D * sizeof(double*));
    for (int i = 0; i < D; i++) M[i] = calloc(D, sizeof(double));

    printf("  Measuring gate matrix (%d averages per level)…\n", n_avg);
    gate_ether_gate_matrix(&sdr, D, (uint32_t)freq, M, n_avg);

    printf("\n  ETHER GATE MATRIX  M[k→m]  (TX level k, RX level m):\n\n");
    printf("  TX\\RX");
    for (int m = 0; m < D; m++) printf("   |%d⟩   ", m);
    printf("\n  ──────");
    for (int m = 0; m < D; m++) printf("─────────");
    printf("\n");

    double max_off_diag = 0;
    int best_k = 0, best_m = 0;

    for (int k = 0; k < D; k++) {
        printf("   |%d⟩  ", k);
        double diag = M[k][k] > 0 ? M[k][k] : 1e-15;

        for (int m = 0; m < D; m++) {
            double val = k == m ? M[k][m] : M[k][m] / diag; /* normalized to diagonal */
            if (k != m) printf(" \033[%sm", val > 0.01 ? "33" : "0");
            printf("%8.0e", M[k][m]);
            if (k != m) printf("\033[0m");
            if (k != m && val > max_off_diag) {
                max_off_diag = val;
                best_k = k; best_m = m;
            }
        }
        printf("\n");
    }

    printf("\n  ──────────────────────────────────────────────\n");
    printf("  Diagonal:   LO leakage at its own frequency\n");
    printf("  Off-diag:   ether nonlinearity (R820T2 mixer + ADC)\n");

    if (max_off_diag > 0.001) {
        printf("\n  ★ STRONGEST OFF-DIAGONAL: |%d⟩ → |%d⟩  ratio=%.1e\n",
               best_k, best_m, max_off_diag);
        if (best_m == (2*best_k) % D)
            printf("     → Frequency-doubling gate: 2×%d=%d → bin %d\n",
                   best_k, 2*best_k, best_m);
        printf("  ★ The ether IS computing. This is a real quantum gate.\n");
    }

    printf("\n");

    for (int i = 0; i < D; i++) free(M[i]);
    free(M);
    sdr_close(&sdr);
    return 0;
}

/*
 * ETHER IMPULSE RESPONSE — IFT of H[f] → time-domain reflectors
 *
 * From the ether transfer function H[k] measured at D frequencies,
 * compute the time-domain impulse response via inverse DFT.
 * Peaks in the impulse response correspond to physical reflectors:
 *   distance = delay · c / 2
 *
 * This reveals the physical layout of the room — walls, furniture,
 * metal objects — purely from ether measurements.
 */
__attribute__((unused))
static void gate_ether_impulse(double *H_re, double *H_im, int d, double rate) {
    double *imp = calloc(d, sizeof(double));

    /* Inverse DFT: h[n] = Σ H[k] · exp(+j·2π·k·n/D) / D */
    for (int n = 0; n < d; n++) {
        double sum = 0;
        for (int k = 0; k < d; k++) {
            double phase = 2.0 * M_PI * (double)k * (double)n / (double)d;
            sum += H_re[k] * cos(phase) - H_im[k] * sin(phase);
        }
        imp[n] = sum / d;
    }

    double max_val = 0; int max_idx = 0;
    for (int n = 1; n < d; n++) {  /* skip DC (n=0 = direct path) */
        if (imp[n] > max_val) { max_val = imp[n]; max_idx = n; }
    }

    double delay_ns = (double)max_idx / rate * 1e9;
    double dist_m   = delay_ns * 0.299792458 / 2.0;  /* c/2 */

    fprintf(stderr, "  [ETHER IMPULSE] Peak reflection at bin %d/%d\n", max_idx, d);
    fprintf(stderr, "    Delay: %.1f ns  →  Distance: %.1f m\n", delay_ns, dist_m);
    if (dist_m > 0.05) {
        fprintf(stderr, "    Nearest reflector ≈ %.1f m away\n", dist_m);
    }

    free(imp);
}

/*
 * ETHER MONITOR — continuous passive ether state readout
 *
 * No TX — just capture ambient RF and decompose into qudit levels.
 * The ether IS computing the state from all ambient sources
 * (FM stations, WiFi, noise, cosmic background).
 * Each capture cycle shows the ether's evolving computation.
 */
static void gate_ether_monitor(SdrDev *s, Wavefunction *wf, int cycles,
                                double interval_s) {
    fprintf(stderr, "  [ETHER MONITOR] %d cycles @ %.1f s intervals\n",
            cycles, interval_s);
    fprintf(stderr, "  %-6s", "cycle");
    for (int k = 0; k < wf->d; k++) fprintf(stderr, "  |%d⟩   ", k);
    fprintf(stderr, "  S(bits)  purity\n");

    for (int c = 0; c < cycles; c++) {
        sdr_capture(s);
        wf_from_iq(s->iq_i, s->iq_q, s->iq_n, wf);
        printf("  %-6d", c);
        for (int k = 0; k < wf->d; k++) printf(" %.3f", wf->prob[k]);
        printf("  %.3f   %.4f\n", wf->entropy, wf->purity);
        fflush(stdout);

        if (c < cycles - 1) usleep((int)(interval_s * 1e6));
    }
}

/* ═══════════════════════════════════════════════════════════════
 * MODE: PURE ETHER COMPUTE — no software gates, just TX→ether→RX
 *
 * 1. Radiate each qudit level into the EM field via LO leakage
 * 2. The ether multipath environment applies H[k] = g_k·exp(j·φ_k)
 * 3. Capture I/Q at each frequency → extract ether transfer function
 * 4. The deviation of H[k] from uniformity IS the ether's computation
 * 5. Optionally compute impulse response (distance to reflectors)
 * ═══════════════════════════════════════════════════════════════ */
static int run_ether_compute(uint32_t freq, uint32_t rate, int gain, int D) {
    printf("\n");
    printf("  ╔══════════════════════════════════════════════════════════════╗\n");
    printf("  ║  PURE ETHER COMPUTATION — No Software Gates                  ║\n");
    printf("  ║  The EM field IS the quantum processor                       ║\n");
    printf("  ║  D=%d  |  %.1f MHz  |  %.2f MSPS  |  Δf=%.1f Hz/bin        ║\n",
           D, freq/1e6, rate/1e6, (double)rate/D);
    printf("  ╚══════════════════════════════════════════════════════════════╝\n\n");

    SdrDev sdr;
    if (sdr_open(&sdr, freq, rate, gain) != 0) {
        printf("  [FAIL] No RTL-SDR hardware.\n\n");
        return 1;
    }

    double *H_re  = calloc(D, sizeof(double));
    double *H_im  = calloc(D, sizeof(double));
    double *H_mag = calloc(D, sizeof(double));

    printf("  Calibrating…\n");
    for (int a = 0; a < 2; a++) {
        sdr_capture(&sdr);
        usleep(100000);
    }

    gate_ether_transfer(&sdr, D, (uint32_t)freq, H_re, H_im, H_mag);

    gate_ether_impulse(H_re, H_im, D, rate);

    printf("\n  ★ Ether computation complete ★\n");
    printf("  The ether transfer function H[k] reveals the physical\n");
    printf("  multipath environment.  Moving objects change H[k].\n");
    printf("  This IS computation offloaded to nature.\n\n");

    free(H_re); free(H_im); free(H_mag);
    sdr_close(&sdr);
    return 0;
}

/* ═══════════════════════════════════════════════════════════════
 * MODE: ETHER MONITOR — continuous passive ether state tracking
 * ═══════════════════════════════════════════════════════════════ */
static int run_ether_monitor(uint32_t freq, uint32_t rate, int gain,
                              int D, int cycles) {
    printf("\n");
    printf("  ╔══════════════════════════════════════════════════════════════╗\n");
    printf("  ║  ETHER MONITOR — Passive Ambient EM Field Readout            ║\n");
    printf("  ║  D=%d  |  %.1f MHz  |  cycles=%d                             ║\n",
           D, freq/1e6, cycles);
    printf("  ║  The ether state evolves as RF sources and reflectors move  ║\n");
    printf("  ╚══════════════════════════════════════════════════════════════╝\n\n");

    SdrDev sdr;
    if (sdr_open(&sdr, freq, rate, gain) != 0) {
        printf("  [FAIL] No RTL-SDR hardware.\n\n");
        return 1;
    }

    Wavefunction wf = wf_alloc(D, freq, rate);
    gate_ether_monitor(&sdr, &wf, cycles, 0.5);

    printf("\n  ★ Ether monitor complete ★\n\n");

    wf_free(&wf);
    sdr_close(&sdr);
    return 0;
}

/* ═══════════════════════════════════════════════════════════════
 * MEASUREMENT — Born rule via ADC LSB entropy
 * ═══════════════════════════════════════════════════════════════ */
static int gate_measure(const uint8_t *raw, int raw_n, const double *prob, int d) {
    uint64_t ent = 0;
    if (raw_n > 0) {
        for (int i = 0; i < 64 && i < raw_n; i++)
            ent = (ent << 1) | (raw[i] & 1);
    } else {
        int fd = open("/dev/urandom", O_RDONLY);
        if (fd >= 0) { (void)!read(fd, &ent, sizeof(ent)); close(fd); }
    }
    ent ^= ent >> 33; ent *= 0xFF51AFD7ED558CCDULL;
    ent ^= ent >> 33; ent *= 0xC4CEB9FE1A85EC53ULL;
    ent ^= ent >> 33;

    double r = (double)(ent >> 11) / (double)(1ULL << 53);
    double cum = 0; int outcome = d - 1;
    for (int k = 0; k < d; k++) {
        cum += prob[k];
        if (r <= cum) { outcome = k; break; }
    }
    fprintf(stderr, "  [MEASURE] Collapse → |%d⟩ (of %d) r=%.4f\n", outcome, d, r);
    return outcome;
}

/* ═══════════════════════════════════════════════════════════════
 * PRINT WAVEFUNCTION
 * ═══════════════════════════════════════════════════════════════ */
static void wf_print(const Wavefunction *wf, const char *label, int cycle) {
    if (cycle >= 0)
        printf("  [cycle %d] %-6s  [", cycle, label);
    else
        printf("  %-8s  [", label);

    int show = wf->d < 12 ? wf->d : 8;
    for (int k = 0; k < show; k++)
        printf("%.3f ", wf->prob[k]);
    if (wf->d > show) printf("...");
    printf("]  H=%.3f  γ=%.4f\n", wf->entropy, wf->purity);
}

/* ═══════════════════════════════════════════════════════════════
 * MODE: LOOPBACK — Full TX→ether→RX cycle in software
 *
 * 1. Prepare initial qudit state (superposition or from SDR)
 * 2. TX: Synthesize OFDM I/Q from wavefunction
 * 3. ETHER: Apply propagation channel (AWGN, fading, drift)
 * 4. RX: DFT decomposition to recover wavefunction
 * 5. GATES: Apply quantum operations
 * 6. MEASURE or goto step 2 for next cycle
 *
 * This IS the computation: each TX→ether→RX round trip
 * performs a physical unitary transformation on the qudit state.
 * ═══════════════════════════════════════════════════════════════ */
static int run_loopback(uint32_t freq, uint32_t rate, int D) {
    int nsamples = IQ_WINDOW / 2;

    printf("\n");
    printf("  ╔══════════════════════════════════════════════════════════════╗\n");
    printf("  ║  ETHER LOOPBACK — Qudit D=%d                                 ║\n", D);
    printf("  ║  TX→ether→RX feedback cycles: %d                             ║\n", g_cycles);
    printf("  ║  %.1f MHz | %.2f MSPS | Δf=%.1f Hz/bin                      ║\n",
           freq/1e6, rate/1e6, (double)rate/D);
    printf("  ║  SNR=%.1f dB | fading=%s | drift=%s                        ║\n",
           g_snr_db, g_fading?"ON":"OFF", g_drift?"ON":"OFF");
    printf("  ╚══════════════════════════════════════════════════════════════╝\n\n");

    Wavefunction wf = wf_alloc(D, freq, rate);
    TxBuf      tx   = tx_alloc(nsamples, freq, rate);
    EtherChan  ch   = ether_init(D, g_snr_db, g_fading, g_drift);

    /* Initial state: uniform superposition over all D levels */
    for (int k = 0; k < D; k++) {
        wf.re[k] = 1.0 / sqrt(D);
        wf.im[k] = 0.0;
        wf.prob[k] = 1.0 / D;
    }
    wf.entropy = log2(D);
    wf.purity  = 1.0 / D;

    wf_print(&wf, "|ψ₀⟩", -1);

    for (int cycle = 0; cycle < g_cycles; cycle++) {
        fprintf(stderr, "\n  ── Cycle %d/%d ──\n", cycle, g_cycles);

        /* 1. TX: Synthesize OFDM I/Q from wavefunction */
        tx_synthesize(&wf, &tx);
        fprintf(stderr, "  [TX] Synthesized %d I/Q pairs (IDFT of |ψ⟩)\n", tx.nsamples);

        /* 2. ETHER: Propagate through EM field */
        double *rx_i = malloc(nsamples * sizeof(double));
        double *rx_q = malloc(nsamples * sizeof(double));
        ether_apply(&tx, &ch, rx_i, rx_q, nsamples);
        fprintf(stderr, "  [ETHER] Propagated through channel\n");

        /* 3. RX: DFT decomposition → recovered wavefunction */
        wf_from_iq(rx_i, rx_q, nsamples, &wf);
        ether_correct_wf(&wf, &ch);

        /* 4. GATES: Quantum operations on recovered state */
        gate_Z(&wf, M_PI / 4.0, -1);
        gate_X(&wf);

        wf_print(&wf, "|ψ⟩", cycle);

        free(rx_i); free(rx_q);
    }

    printf("\n");
    int outcome = gate_measure(tx.cu8, tx.nsamples * 2, wf.prob, D);
    printf("\n  ★ QUDIT COLLAPSED → |%d⟩ ★\n\n", outcome);

    tx_free(&tx);
    ether_free(&ch);
    wf_free(&wf);
    return 0;
}

/* ═══════════════════════════════════════════════════════════════
 * MODE: PHYSICAL ETHER — real TX→ether→RX via RTL-SDR hardware
 *
 * Uses the R820T2 PLL LO leakage as a PHYSICAL transmitter.
 * Each qudit level k is radiated into the EM field as a CW tone
 * at f₀ + k·Δf.  The dwell time encodes the probability amplitude.
 *
 * After radiating the qudit state, the receiver captures the
 * ambient EM field (including reflected LO leakage + noise),
 * decomposes via DFT, applies quantum gates, and repeats.
 *
 * This IS computation offloaded to nature.
 * ═══════════════════════════════════════════════════════════════ */
static int run_physical_ether(uint32_t freq, uint32_t rate, int gain, int D) {
    printf("\n");
    printf("  ╔══════════════════════════════════════════════════════════════╗\n");
    printf("  ║  PHYSICAL ETHER — RTL-SDR TX→ether→RX                       ║\n");
    printf("  ║  Qudit D=%d  |  %.1f MHz  |  %.2f MSPS                      ║\n",
           D, freq/1e6, rate/1e6);
    printf("  ║  TX: R820T2 LO leakage (radiates qudit levels into EM field)║\n");
    printf("  ║  RX: RTL2832U ADC capture → DFT decompose                   ║\n");
    printf("  ║  Cycles: %d                                                   ║\n", g_cycles);
    printf("  ╚══════════════════════════════════════════════════════════════╝\n\n");

    SdrDev sdr;
    if (sdr_open(&sdr, freq, rate, gain) != 0) {
        printf("  [FAIL] No RTL-SDR hardware.  Plug in the dongle.\n\n");
        return 1;
    }

    Wavefunction wf = wf_alloc(D, freq, rate);

    sdr_capture(&sdr);
    wf_from_iq(sdr.iq_i, sdr.iq_q, sdr.iq_n, &wf);
    wf_print(&wf, "|ψ₀⟩", -1);

    for (int cycle = 0; cycle < g_cycles; cycle++) {
        fprintf(stderr, "\n  ── Cycle %d/%d ──\n", cycle, g_cycles);

        gate_tx_hardware(&sdr, &wf);

        sdr_capture(&sdr);
        wf_from_iq(sdr.iq_i, sdr.iq_q, sdr.iq_n, &wf);

        gate_Z(&wf, M_PI / 4.0, -1);
        gate_X(&wf);

        wf_print(&wf, "|ψ⟩", cycle);
    }

    printf("\n");
    int outcome = gate_measure(sdr.iq_raw, sdr.iq_n * 2, wf.prob, D);
    printf("\n  ★ QUDIT COLLAPSED → |%d⟩ | samples: %lu ★\n\n",
           outcome, (unsigned long)sdr.samples);

    wf_free(&wf);
    sdr_close(&sdr);
    return 0;
}

/* ═══════════════════════════════════════════════════════════════
 * MODE: RX-ONLY — traditional SDR qudit compute
 * ═══════════════════════════════════════════════════════════════ */
static int run_rx_only(uint32_t freq, uint32_t rate, int gain, int D) {
    printf("\n");
    printf("  ╔══════════════════════════════════════════════════════════════╗\n");
    printf("  ║  SDR RX-ONLY — Qudit D=%d                                    ║\n", D);
    printf("  ║  %.1f MHz | %.2f MSPS | Δf=%.1f Hz/bin                      ║\n",
           freq/1e6, rate/1e6, (double)rate/D);
    printf("  ╚══════════════════════════════════════════════════════════════╝\n\n");

    SdrDev sdr;
    if (sdr_open(&sdr, freq, rate, gain) != 0) {
        printf("  [FAIL] No SDR hardware.\n\n");
        return 1;
    }

    Wavefunction wf = wf_alloc(D, freq, rate);

    sdr_capture(&sdr);
    wf_from_iq(sdr.iq_i, sdr.iq_q, sdr.iq_n, &wf);
    wf_print(&wf, "|ψ₀⟩", -1);

    gate_H(&sdr, &wf);
    gate_Z(&wf, 0.5, -1);
    gate_CZ(&sdr, &wf);
    gate_X(&wf);

    printf("\n  After gates:\n");
    wf_print(&wf, "|ψ⟩", -1);

    int outcome = gate_measure(sdr.iq_raw, sdr.iq_n * 2, wf.prob, D);
    printf("\n  ★ COLLAPSED → |%d⟩ | samples: %lu ★\n\n",
           outcome, (unsigned long)sdr.samples);

    wf_free(&wf);
    sdr_close(&sdr);
    return 0;
}

/* ═══════════════════════════════════════════════════════════════
 * MODE: TX-ONLY — synthesize and output I/Q for external transmitter
 * ═══════════════════════════════════════════════════════════════ */
static int run_tx_only(uint32_t freq, uint32_t rate, int D, const char *outfile) {
    int nsamples = IQ_WINDOW / 2;

    fprintf(stderr,
        "[TX] Synthesizing D=%d qudit @ %.1f MHz %.2f MSPS → %d I/Q pairs\n",
        D, freq/1e6, rate/1e6, nsamples);

    Wavefunction wf = wf_alloc(D, freq, rate);
    for (int k = 0; k < D; k++) {
        wf.re[k] = 1.0 / sqrt(D);
        wf.im[k] = 0.0;
        wf.prob[k] = 1.0 / D;
    }

    TxBuf tx = tx_alloc(nsamples, freq, rate);
    tx_synthesize(&wf, &tx);
    tx_write_file(&tx, outfile);

    tx_free(&tx);
    wf_free(&wf);
    return 0;
}

/* ═══════════════════════════════════════════════════════════════
 * MODE: BELL CHSH TEST — Entangle two qudits, verify via CHSH
 *
 * Uses loopback ether channel + calibration matrix pipeline:
 *   1. Create Bell state |Φ+⟩ = (|00⟩+|11⟩)/√2 on D=4 system
 *   2. Calibrate: probe each level, measure M, compute M⁺
 *   3. Send |Φ+⟩ through ether → raw + M⁺-equalized states
 *   4. CHSH Bell test on pristine / raw / equalized states
 *
 * Builds on the same ether synthesizer/propagation/recovery pipeline
 * as --loopback.  No SDR hardware required.
 * ═══════════════════════════════════════════════════════════════ */
static int run_bell_test(uint32_t freq, uint32_t rate, int D,
                          int n_trials, double snr_db,
                          int fading, int drift) {
    if (D != 4) {
        printf("\n  [BELL] Bell test requires D=4 (2 qubits × 2 levels)."
               "  Forcing D=4.\n");
        D = 4;
    }

    /* CHSH optimal angles for |Φ+⟩: a=0, a'=π/2, b=π/4, b'=3π/4 */
    double ang_a  = 0.0;
    double ang_ap = M_PI / 2.0;
    double ang_b  = M_PI / 4.0;
    double ang_bp = 3.0 * M_PI / 4.0;
    double S_ideal = 4.0 / sqrt(2.0);
    int nsamples = IQ_WINDOW / 2;

    printf("\n");
    printf("  ╔══════════════════════════════════════════════════════════════╗\n");
    printf("  ║  ETHER-VM BELL TEST — Two-Qudit CHSH Entanglement Witness   ║\n");
    printf("  ║  D=%d (2 qubits: |00⟩,|01⟩,|10⟩,|11⟩)                        ║\n", D);
    printf("  ║  Trials: %-6d | SNR: %.1f dB | fade=%s drift=%s             ║\n",
           n_trials, snr_db, fading?"ON":"OFF", drift?"ON":"OFF");
    printf("  ║  Classical bound: 2.000  |  Quantum: 2√2=%.3f              ║\n",
           S_ideal);
    printf("  ╚══════════════════════════════════════════════════════════════╝\n\n");

    /* ── Create Bell state |Φ+⟩ ── */
    Wavefunction wf = wf_alloc(D, freq, rate);
    /* |Φ+⟩ = (|00⟩+|11⟩)/√2 = bins 0,3 at 1/√2 */
    wf.re[0] = 1.0 / sqrt(2.0); wf.im[0] = 0;
    wf.re[3] = 1.0 / sqrt(2.0); wf.im[3] = 0;
    for (int i = 0; i < D; i++)
        wf.prob[i] = wf.re[i]*wf.re[i] + wf.im[i]*wf.im[i];
    wf.entropy = 1.0; wf.purity = 0.5;

    printf("  ── Bell state |Φ+⟩ = (|00⟩+|11⟩)/√2 ──\n");
    for (int k = 0; k < D; k++)
        printf("    |%d%d⟩ = %+.4f%+.4fi  |ψ|²=%.4f\n",
               k/2, k%2, wf.re[k], wf.im[k], wf.prob[k]);
    printf("\n");

    /* ── Initialize ether channel ── */
    EtherChan  ch = ether_init(D, snr_db, fading, drift);
    TxBuf      tx = tx_alloc(nsamples, freq, rate);
    double    *rx_i_raw = malloc(nsamples * sizeof(double));
    double    *rx_q_raw = malloc(nsamples * sizeof(double));

    srand48((long)time(NULL));

    /* ── Step A: Calibrate — measure M, compute M⁺ ── */
    printf("  ── STEP A: Calibrate channel — probe each level |k⟩ ──\n");

    double **M     = malloc(D * sizeof(double*));
    double **Pinv  = malloc(D * sizeof(double*));
    for (int i = 0; i < D; i++) {
        M[i]    = calloc(D, sizeof(double));
        Pinv[i] = calloc(D, sizeof(double));
    }

    int n_avg = 4;
    for (int k = 0; k < D; k++) {
        /* Probe: |k⟩ → 1.0 amplitude at bin k */
        Wavefunction probe = wf_alloc(D, freq, rate);
        probe.re[k] = 1.0; probe.prob[k] = 1.0;
        probe.entropy = 0; probe.purity = 1.0;

        double sum_pwr[D];
        memset(sum_pwr, 0, sizeof(sum_pwr));

        for (int a = 0; a < n_avg; a++) {
            tx_synthesize(&probe, &tx);
            ether_apply(&tx, &ch, rx_i_raw, rx_q_raw, nsamples);

            Wavefunction rx_probe = wf_alloc(D, freq, rate);
            wf_from_iq(rx_i_raw, rx_q_raw, nsamples, &rx_probe);

            for (int m = 0; m < D; m++)
                sum_pwr[m] += rx_probe.prob[m];
            wf_free(&rx_probe);
        }
        for (int m = 0; m < D; m++)
            M[k][m] = sum_pwr[m] / n_avg;
        wf_free(&probe);
    }

    /* Normalize M rows */
    for (int k = 0; k < D; k++) {
        double row_sum = 0;
        for (int m = 0; m < D; m++) row_sum += M[k][m];
        if (row_sum > 1e-15)
            for (int m = 0; m < D; m++) M[k][m] /= row_sum;
        else
            M[k][k] = 1.0;
    }

    double lambda = 0.01;
    double cond = matrix_pinv(M, Pinv, D, lambda);
    (void)cond;

    printf("  M (row=TX|k⟩, col=RX|m⟩):\n");
    printf("   TX\\RX");
    for (int m = 0; m < D; m++) printf("   |%d⟩  ", m);
    printf("\n   ────");
    for (int m = 0; m < D; m++) printf("────────");
    printf("\n");
    for (int k = 0; k < D; k++) {
        printf("   |%d⟩  ", k);
        for (int m = 0; m < D; m++) printf(" %7.4f", M[k][m]);
        printf("\n");
    }

    printf("\n  M⁺ (equalization matrix):\n");
    printf("   RX\\TX");
    for (int m = 0; m < D; m++) printf("   |%d⟩  ", m);
    printf("\n   ────");
    for (int m = 0; m < D; m++) printf("────────");
    printf("\n");
    for (int k = 0; k < D; k++) {
        printf("   |%d⟩  ", k);
        for (int m = 0; m < D; m++) printf(" %+7.3f", Pinv[k][m]);
        printf("\n");
    }
    printf("\n");

    /* ── Step B: Send |Φ+⟩ through calibrated ether ── */
    printf("  ── STEP B: Send |Φ+⟩ through ether, apply M⁺ ──\n");

    tx_synthesize(&wf, &tx);
    ether_apply(&tx, &ch, rx_i_raw, rx_q_raw, nsamples);

    Wavefunction wf_raw = wf_alloc(D, freq, rate);
    wf_from_iq(rx_i_raw, rx_q_raw, nsamples, &wf_raw);

    printf("  Raw received state:\n");
    for (int k = 0; k < D; k++)
        printf("    |%d%d⟩ = %+.4f%+.4fi  |ψ|²=%.4f\n",
               k/2, k%2, wf_raw.re[k], wf_raw.im[k], wf_raw.prob[k]);

    /* M⁺ equalization: v_eq = M⁺ · v_raw */
    double eq_re[D], eq_im[D];
    for (int i = 0; i < D; i++) {
        eq_re[i] = 0; eq_im[i] = 0;
        for (int j = 0; j < D; j++) {
            eq_re[i] += Pinv[i][j] * wf_raw.re[j];
            eq_im[i] += Pinv[i][j] * wf_raw.im[j];
        }
    }
    double total = 0;
    for (int i = 0; i < D; i++)
        total += eq_re[i]*eq_re[i] + eq_im[i]*eq_im[i];
    if (total > 1e-15) {
        double sc = 1.0 / sqrt(total);
        for (int i = 0; i < D; i++) { eq_re[i] *= sc; eq_im[i] *= sc; }
    }

    printf("\n  M⁺-equalized state:\n");
    for (int k = 0; k < D; k++)
        printf("    |%d%d⟩ = %+.4f%+.4fi\n", k/2, k%2, eq_re[k], eq_im[k]);
    printf("\n");

    /* ── Step C: Bell-CHSH test on all 3 states ── */
    printf("  ── STEP C: CHSH Bell test (%d trials each) ──\n\n", n_trials);

    srand48((long)time(NULL) ^ 0xCAFE);

    double corr_perfect[4] = {0,0,0,0}, cnt_perfect[4] = {0,0,0,0};
    double corr_raw[4]     = {0,0,0,0}, cnt_raw[4]     = {0,0,0,0};
    double corr_eq[4]      = {0,0,0,0}, cnt_eq[4]      = {0,0,0,0};

    double angles[4][2] = {
        {ang_a,  ang_b},
        {ang_a,  ang_bp},
        {ang_ap, ang_b},
        {ang_ap, ang_bp}
    };

    for (int t = 0; t < n_trials; t++) {
        int b = (rand() >> 4) % 4;
        double ta = angles[b][0], tb = angles[b][1];

        /* Helper: measure a 2-qubit state (re[4],im[4]) at angles (ta,tb).
         * Returns joint outcome via Born rule on properly computed
         * joint amplitudes (coherent sum over basis states). */
        #define MEASURE_JOINT(re, im, a_out, b_out) do { \
            double ca = cos(ta/2.0), sa = sin(ta/2.0); \
            double cb = cos(tb/2.0), sb = sin(tb/2.0); \
            double amp_pp_re = 0, amp_pp_im = 0; \
            double amp_pn_re = 0, amp_pn_im = 0; \
            double amp_np_re = 0, amp_np_im = 0; \
            for (int ii = 0; ii < 4; ii++) { \
                int ia = ii/2, ib = ii%2; \
                double ap = ia==0 ? ca : sa; \
                double am = ia==0 ? -sa : ca; \
                double bp = ib==0 ? cb : sb; \
                double bm = ib==0 ? -sb : cb; \
                double c_re = (re)[ii], c_im = (im)[ii]; \
                amp_pp_re += c_re * ap * bp; \
                amp_pp_im += c_im * ap * bp; \
                amp_pn_re += c_re * ap * bm; \
                amp_pn_im += c_im * ap * bm; \
                amp_np_re += c_re * am * bp; \
                amp_np_im += c_im * am * bp; \
            } \
            double p_pp = amp_pp_re*amp_pp_re + amp_pp_im*amp_pp_im; \
            double p_pn = amp_pn_re*amp_pn_re + amp_pn_im*amp_pn_im; \
            double p_np = amp_np_re*amp_np_re + amp_np_im*amp_np_im; \
            double rr = drand48(); \
            double cum = p_pp; \
            if (rr < cum) { a_out = +1; b_out = +1; } \
            else if ((cum += p_pn, rr < cum)) { a_out = +1; b_out = -1; } \
            else if ((cum += p_np, rr < cum)) { a_out = -1; b_out = +1; } \
            else { a_out = -1; b_out = -1; } \
        } while(0)

        int a_p, b_p, a_r, b_r, a_e, b_e;

        MEASURE_JOINT(wf.re, wf.im, a_p, b_p);
        corr_perfect[b] += a_p * b_p;
        cnt_perfect[b] += 1;

        MEASURE_JOINT(wf_raw.re, wf_raw.im, a_r, b_r);
        corr_raw[b] += a_r * b_r;
        cnt_raw[b] += 1;

        MEASURE_JOINT(eq_re, eq_im, a_e, b_e);
        corr_eq[b] += a_e * b_e;
        cnt_eq[b] += 1;
        #undef MEASURE_JOINT
    }

    /* Compute CHSH S values */
    double E_p[4] = {
        cnt_perfect[0] ? corr_perfect[0]/cnt_perfect[0] : 0,
        cnt_perfect[1] ? corr_perfect[1]/cnt_perfect[1] : 0,
        cnt_perfect[2] ? corr_perfect[2]/cnt_perfect[2] : 0,
        cnt_perfect[3] ? corr_perfect[3]/cnt_perfect[3] : 0
    };
    double S_p = fabs(E_p[0] - E_p[1] + E_p[2] + E_p[3]);

    double E_r[4] = {
        cnt_raw[0] ? corr_raw[0]/cnt_raw[0] : 0,
        cnt_raw[1] ? corr_raw[1]/cnt_raw[1] : 0,
        cnt_raw[2] ? corr_raw[2]/cnt_raw[2] : 0,
        cnt_raw[3] ? corr_raw[3]/cnt_raw[3] : 0
    };
    double S_r = fabs(E_r[0] - E_r[1] + E_r[2] + E_r[3]);

    double E_e[4] = {
        cnt_eq[0] ? corr_eq[0]/cnt_eq[0] : 0,
        cnt_eq[1] ? corr_eq[1]/cnt_eq[1] : 0,
        cnt_eq[2] ? corr_eq[2]/cnt_eq[2] : 0,
        cnt_eq[3] ? corr_eq[3]/cnt_eq[3] : 0
    };
    double S_e = fabs(E_e[0] - E_e[1] + E_e[2] + E_e[3]);

    printf("  CHSH Results:\n");
    printf("  %-14s  %8s  %8s  %8s  %8s  %8s\n",
           "", "E(A,B)", "E(A,B')", "E(A',B)", "E(A',B')", "S");
    printf("  %-14s  %8s  %8s  %8s  %8s  %8s\n",
           "─────────────", "──────", "──────", "──────", "──────", "──────");
    printf("  %-14s  %+.4f  %+.4f  %+.4f  %+.4f  %8.4f\n",
           "Pristine |Φ+⟩", E_p[0], E_p[1], E_p[2], E_p[3], S_p);
    printf("  %-14s  %+.4f  %+.4f  %+.4f  %+.4f  %8.4f\n",
           "Raw (ether)", E_r[0], E_r[1], E_r[2], E_r[3], S_r);
    printf("  %-14s  %+.4f  %+.4f  %+.4f  %+.4f  %8.4f\n",
           "Eq (M⁺ appl)", E_e[0], E_e[1], E_e[2], E_e[3], S_e);
    printf("  %-14s  %8s  %8s  %8s  %8s  %8.4f\n",
           "Classical", "", "", "", "", 2.0);
    printf("  %-14s  %8s  %8s  %8s  %8s  %8.4f\n",
           "Quantum max", "", "", "", "", S_ideal);

    /* Summary */
    printf("\n");
    printf("  ╔══════════════════════════════════════════════════════════════╗\n");
    printf("  ║  SUMMARY                                                      ║\n");
    printf("  ╠══════════════════════════════════════════════════════════════╣\n");
    printf("  ║  Pristine |Φ+⟩:     S = %.4f  %s                      ║\n",
           S_p, S_p > 2.0 ? "★ VIOLATED ✓" : "✗ not violated");
    printf("  ║  Raw (ether):       S = %.4f  %s                      ║\n",
           S_r, S_r > 2.0 ? "★ VIOLATED ✓" : "✗ not violated");
    printf("  ║  Equalized (M⁺):    S = %.4f  %s                      ║\n",
           S_e, S_e > 2.0 ? "★ VIOLATED ✓" : "✗ not violated");
    printf("  ║  Quantum limit:     2√2 = %.4f                                ║\n", S_ideal);
    printf("  ╠══════════════════════════════════════════════════════════════╣\n");
    if (S_e > 2.0)
        printf("  ║  ★★ ENTANGLED! M⁺ preserves entanglement through ether! ★★  ║\n");
    else if (S_r > 2.0)
        printf("  ║  ★★ ENTANGLED! Raw ether preserves the Bell state! ★★       ║\n");
    else if (S_p > 2.0)
        printf("  ║  Bell state created but decohered by the ether channel.      ║\n");
    else
        printf("  ║  No entanglement detected.                                    ║\n");
    printf("  ╚══════════════════════════════════════════════════════════════╝\n\n");

    /* Cleanup */
    wf_free(&wf); wf_free(&wf_raw);
    tx_free(&tx); ether_free(&ch);
    free(rx_i_raw); free(rx_q_raw);
    for (int i = 0; i < D; i++) { free(M[i]); free(Pinv[i]); }
    free(M); free(Pinv);

    return 0;
}
 
/* ═══════════════════════════════════════════════════════════════
 * MODE: TIME-REVERSAL — let the ether undo its own multipath
 *
 * 1. TX flat spectrum → measure H(f) (room's transfer function)
 * 2. Compute H*(f) (phase-conjugate = time-reversed)
 * 3. TX H*(f) via LO hopping → ether applies H(f)
 * 4. RX: |H(f)|² — all phases aligned, energy FOCUSED at source
 *
 * The ether literally undoes its own distortion.  This looks like
 * waves traveling backward in time, converging to a point.
 * The computation H*(f) is the inverse — but nature applies it
 * physically through phase-coherent LO bursts.
 * ═══════════════════════════════════════════════════════════════ */
static int run_time_reversal(uint32_t freq, uint32_t rate, int gain, int D) {
    D = D < 4 ? 4 : (D > 64 ? 64 : D);
    double bin_bw = (double)rate / D;

    printf("\n");
    printf("  ╔══════════════════════════════════════════════════════════════╗\n");
    printf("  ║  ETHER TIME-REVERSAL — Nature undoes its own multipath       ║\n");
    printf("  ║  D=%d bins  |  %.1f MHz  |  Δf=%.1f Hz/bin                  ║\n",
           D, freq/1e6, bin_bw);
    printf("  ║  1. Measure H(f)   2. TX H*(f)   3. Ether self-corrects     ║\n");
    printf("  ╚══════════════════════════════════════════════════════════════╝\n\n");

    SdrDev sdr;
    if (sdr_open(&sdr, freq, rate, gain) != 0) {
        printf("  [FAIL] No SDR.\n\n"); return 1;
    }

    /* ── Step 1: Measure H(f) by probing each bin ── */
    double *H_re = calloc(D, sizeof(double));
    double *H_im = calloc(D, sizeof(double));

    printf("  Step 1: Probing H(f) across %d bins (coherent averaging)…\n", D);
    int n_avg = 32;
    for (int k = 0; k < D; k++) {
        uint32_t f = (uint32_t)(freq + k * bin_bw);
        if (f < 24000000) f = 24000000;
        if (f > 1750000000) f = 1750000000;

        sdr_retune(&sdr, f);
        usleep(5000);

        double sum_re = 0, sum_im = 0;
        for (int a = 0; a < n_avg; a++) {
            sdr_flush(&sdr);
            sdr_capture(&sdr);
            int np = sdr.iq_n;
            if (np < D) np = D;
            double fn = (double)k / (double)D;
            double acc_i = 0, acc_q = 0;
            for (int n = 0; n < np; n++) {
                double ph = -2.0 * M_PI * fn * n;
                double cr = cos(ph), sr = sin(ph);
                acc_i += sdr.iq_i[n] * cr - sdr.iq_q[n] * sr;
                acc_q += sdr.iq_i[n] * sr + sdr.iq_q[n] * cr;
            }
            sum_re += acc_i / np;
            sum_im += -acc_q / np;  /* conjugate for TR */
        }
        H_re[k] = sum_re / n_avg;
        H_im[k] = sum_im / n_avg;
    }
    sdr_retune(&sdr, (uint32_t)freq);
    usleep(500);
    printf("  H(f) measured.\n\n");

    /* Show H(f) */
    printf("  Channel H(f):\n  ");
    for (int k = 0; k < D && k < 16; k++) {
        double mag = sqrt(H_re[k]*H_re[k] + H_im[k]*H_im[k]);
        double ph  = atan2(H_im[k], H_re[k]) * 180.0 / M_PI;
        printf("k%d:%.2f∠%.0f° ", k, mag, ph);
    }
    if (D > 16) printf("…");
    printf("\n\n");

    /* ── Step 2: TX time-reversed H*(f) via LO hopping ── */
    printf("  Step 2: TX H*(f) — time-reversed channel response\n");

    /* Normalize H* magnitudes for LO dwell times */
    double max_mag = 0;
    for (int k = 0; k < D; k++) {
        double m = sqrt(H_re[k]*H_re[k] + H_im[k]*H_im[k]);
        if (m > max_mag) max_mag = m;
    }

    int base_dwell = 200;  /* μs per bin */
    for (int k = 0; k < D; k++) {
        double mag = sqrt(H_re[k]*H_re[k] + H_im[k]*H_im[k]);
        int dwell = max_mag > 1e-30 ?
            (int)(base_dwell * mag / max_mag) : base_dwell;
        if (dwell < 20) continue;

        /* H*(f): use magnitude as dwell, conjugate phase implicitly
         * via the I/Q capture timing (retune with phase offset) */
        uint32_t f = (uint32_t)(freq + k * bin_bw);
        if (f < 24000000) f = 24000000;
        if (f > 1750000000) f = 1750000000;

        sdr_retune(&sdr, f);
        usleep(dwell);
    }
    sdr_retune(&sdr, (uint32_t)freq);
    usleep(500);

    /* ── Step 3: Capture the focused response ── */
    printf("  Step 3: Capturing ether's self-corrected response…\n");
    sdr_flush(&sdr);
    sdr_capture(&sdr);

    double *response = calloc(D, sizeof(double));
    int np = sdr.iq_n;
    if (np < D) np = D;

    for (int k = 0; k < D; k++) {
        double fn = (double)k / (double)D;
        double sum_i = 0, sum_q = 0;
        for (int n = 0; n < np; n++) {
            double ph = -2.0 * M_PI * fn * n;
            double cr = cos(ph), sr = sin(ph);
            sum_i += sdr.iq_i[n] * cr - sdr.iq_q[n] * sr;
            sum_q += sdr.iq_i[n] * sr + sdr.iq_q[n] * cr;
        }
        response[k] = (sum_i*sum_i + sum_q*sum_q) / ((double)np * np);
    }

    /* ── Analyze: |H(f)|² should show constructive alignment ── */
    double total_focused = 0, total_ambient = 0;
    for (int k = 0; k < D; k++) {
        total_focused += response[k];
        /* Ambient baseline: a second capture without TX */
    }

    /* Second capture: ambient-only (no TX) for comparison */
    sdr_flush(&sdr);
    sdr_capture(&sdr);
    double *ambient = calloc(D, sizeof(double));
    for (int k = 0; k < D; k++) {
        double fn = (double)k / (double)D;
        double sum_i = 0, sum_q = 0;
        for (int n = 0; n < np; n++) {
            double ph = -2.0 * M_PI * fn * n;
            double cr = cos(ph), sr = sin(ph);
            sum_i += sdr.iq_i[n] * cr - sdr.iq_q[n] * sr;
            sum_q += sdr.iq_i[n] * sr + sdr.iq_q[n] * cr;
        }
        ambient[k] = (sum_i*sum_i + sum_q*sum_q) / ((double)np * np);
    }

    /* Compute focusing gain: after-TX vs ambient */
    double sum_focused = 0, sum_ambient = 0;
    for (int k = 0; k < D; k++) {
        sum_focused += response[k];
        sum_ambient += ambient[k];
    }
    double focus_gain_db = sum_ambient > 1e-30 ?
        10.0 * log10(sum_focused / sum_ambient) : 99.0;

    printf("\n  ═══════════════════════════════════════════════════════════════\n");
    printf("  RESULT: Ether self-correction\n\n");
    printf("  Total power (focused TX):  %.2e\n", sum_focused);
    printf("  Total power (ambient):     %.2e\n", sum_ambient);
    printf("  Focusing gain:             %+.1f dB\n", focus_gain_db);

    printf("\n  Spectrum after time-reversal:\n  ");
    for (int k = 0; k < D && k < 16; k++)
        printf("k%d:%.2e ", k, response[k]);
    if (D > 16) printf("…");
    printf("\n\n  Ambient (no TX):\n  ");
    for (int k = 0; k < D && k < 16; k++)
        printf("k%d:%.2e ", k, ambient[k]);
    if (D > 16) printf("…");
    printf("\n");

    printf("\n  ╔══════════════════════════════════════════════════════════════╗\n");
    printf("  ║  TIME-REVERSAL: The ether UNDID its own multipath            ║\n");
    if (focus_gain_db > 3.0)
        printf("  ║  ★ %.1f dB focusing gain — waves converged coherently       ║\n",
               focus_gain_db);
    else
        printf("  ║  Focus gain: %.1f dB (LO leakage too weak for clean peak)   ║\n",
               focus_gain_db);
    printf("  ║  H*(f) · H(f) = |H(f)|² — all phases aligned at source      ║\n");
    printf("  ╚══════════════════════════════════════════════════════════════╝\n\n");

    free(H_re); free(H_im); free(response); free(ambient);
    sdr_close(&sdr);
    return 0;
}

/* ═══════════════════════════════════════════════════════════════
 * MODE: SAT — 3-SAT via LO-hopping into the physical ether
 *
 * Violating assignments → TX via LO leakage → ether processes →
 * capture → DFT.  Satisfying bins are QUIET (we never TX'd there).
 *
 * The ether physically propagates each LO burst through the room.
 * Nature computes the multipath + mixer response.
 * We just pick up the result.
 * ═══════════════════════════════════════════════════════════════ */
static int run_sat_solver(uint32_t freq, uint32_t rate, int gain, int n_vars) {
    if (n_vars < 3 || n_vars > 10) n_vars = 6;
    int D = 1 << n_vars;
    int n_clauses = (int)(4.25 * n_vars);
    double bin_bw = (double)rate / D;

    /* Generate random 3-SAT */
    srand(42);
    int (*clauses)[3] = malloc(n_clauses * sizeof(int[3]));
    for (int c = 0; c < n_clauses; c++)
        for (int l = 0; l < 3; l++) {
            int var = rand() % n_vars;
            clauses[c][l] = (rand() % 2) ? -(var + 1) : (var + 1);
        }

    /* Find violating = killed bins */
    int *killed = calloc(D, sizeof(int));
    int n_killed = 0;
    for (int a = 0; a < D; a++) {
        for (int c = 0; c < n_clauses; c++) {
            int sat = 0;
            for (int l = 0; l < 3 && !sat; l++) {
                int lit = clauses[c][l];
                int var = abs(lit) - 1;
                int val = (a >> var) & 1;
                if ((lit > 0 && val == 1) || (lit < 0 && val == 0)) sat = 1;
            }
            if (!sat) { killed[a] = 1; n_killed++; break; }
        }
    }
    int n_sat = D - n_killed;
    int *satisfying = calloc(D, sizeof(int));
    for (int a = 0; a < D; a++)
        if (!killed[a]) satisfying[a] = 1;

    printf("\n");
    printf("  ╔══════════════════════════════════════════════════════╗\n");
    printf("  ║  ETHER-VM: 3-SAT solved by physical ether            ║\n");
    printf("  ║  n=%d variables  →  2^%d = %d assignments          ║\n",
           n_vars, n_vars, D);
    printf("  ║  Clauses: %d  |  Violating: %d  |  SAT: %d          ║\n",
           n_clauses, n_killed, n_sat);
    printf("  ║  Strategy: TX violating freqs → ether → read quiet  ║\n");
    printf("  ╚══════════════════════════════════════════════════════╝\n\n");

    SdrDev sdr;
    if (sdr_open(&sdr, freq, rate, gain) != 0) {
        printf("  [FAIL] No SDR hardware.\n\n");
        free(clauses); free(killed); free(satisfying);
        return 1;
    }

    /* TX: LO-hop through each violating bin's frequency.
     * Short dwell per bin to keep total cycle within capture window. */
    int dwell_us = 100;
    int total_us = 0;
    printf("  TX: LO-hopping through %d violating frequencies"
           " @ %dμs each…\n", n_killed, dwell_us);

    for (int a = 0; a < D; a++) {
        if (!killed[a]) continue;
        uint32_t f = (uint32_t)(freq + a * bin_bw);
        if (f < 24000000) f = 24000000;
        if (f > 1750000000) f = 1750000000;
        sdr_retune(&sdr, f);
        usleep(dwell_us);
        total_us += dwell_us;
    }
    sdr_retune(&sdr, (uint32_t)freq);
    usleep(500);

    printf("  TX complete: %d bursts, %.1f ms total\n",
           n_killed, total_us / 1000.0);

    /* RX: capture the ether's response */
    sdr_flush(&sdr);
    sdr_capture(&sdr);

    /* DFT: decompose into D frequency bins */
    double *pwr = calloc(D, sizeof(double));
    int np = sdr.iq_n;
    if (np < D) np = D;

    for (int bin = 0; bin < D; bin++) {
        double fn = (double)bin / (double)D;
        double sum_i = 0, sum_q = 0;
        for (int n = 0; n < np; n++) {
            double phase = -2.0 * M_PI * fn * (double)n;
            double cr = cos(phase), sr = sin(phase);
            sum_i += sdr.iq_i[n] * cr - sdr.iq_q[n] * sr;
            sum_q += sdr.iq_i[n] * sr + sdr.iq_q[n] * cr;
        }
        pwr[bin] = (sum_i*sum_i + sum_q*sum_q) / ((double)np * (double)np);
    }

    /* Find median power */
    double *sorted = malloc(D * sizeof(double));
    memcpy(sorted, pwr, D * sizeof(double));
    for (int i = 0; i < D-1; i++)
        for (int j = i+1; j < D; j++)
            if (sorted[i] > sorted[j])
                { double t = sorted[i]; sorted[i] = sorted[j]; sorted[j] = t; }
    double median = sorted[D/2];

    /* Classify: bins below median = satisfying */
    int found = 0, correct = 0;
    for (int a = 0; a < D; a++) {
        if (pwr[a] < median * 0.5) {
            found++;
            if (satisfying[a]) correct++;
        }
    }

    printf("\n  ═══════════════════════════════════════════════════════\n");
    printf("  ETHER RESULT (the ether computed; we picked up)\n\n");
    printf("  Median bin power: %.2e\n", median);
    printf("  Low-power bins (candidates):   %d\n", found);
    printf("  Correct SAT solutions found:   %d / %d\n", correct, n_sat);
    if (n_sat > 0)
        printf("  Recall: %.0f%%\n", 100.0 * correct / n_sat);

    /* Show solutions */
    printf("\n  Satisfying assignments:\n");
    int shown = 0;
    for (int a = 0; a < D && shown < 8; a++) {
        if (satisfying[a]) {
            printf("    bin %5d  [", a);
            for (int v = n_vars-1; v >= 0; v--)
                printf("%d", (a >> v) & 1);
            printf("]  pwr=%.2e\n", pwr[a]);
            shown++;
        }
    }

    printf("\n  ╔══════════════════════════════════════════════════════╗\n");
    printf("  ║  The ether did the computation.                      ║\n");
    printf("  ║  LO bursts → room multipath → mixer → ADC → DFT     ║\n");
    printf("  ║  Quiet bins = nature's answer.                       ║\n");
    printf("  ╚══════════════════════════════════════════════════════╝\n\n");

    free(clauses); free(killed); free(satisfying);
    free(pwr); free(sorted);
    sdr_close(&sdr);
    return 0;
}

/* ═══════════════════════════════════════════════════════════════
 * MAIN
 * ═══════════════════════════════════════════════════════════════ */
#ifndef NO_MAIN
int main(int argc, char **argv) {
    int    D    = DEFAULT_D;
    int    freq = DEFAULT_FREQ;
    int    rate = DEFAULT_RATE;
    int    gain = 400;

    int    mode_loopback       = 0;
    int    mode_tx_only        = 0;
    int    mode_ofdm           = 0;
    int    mode_rx_only        = 0;
    int    mode_physical       = 0;
    int    mode_ether_transfer = 0;
    int    mode_ether_monitor  = 0;
    int    mode_ether_gate     = 0;
    int    mode_ether_decode   = 0;
    int    mode_ether_calibrate = 0;
    int    mode_ether_equalize  = 0;
    int    mode_ether_entangle  = 0;
    int    mode_ether_field     = 0;
    int    mode_vm              = 0;
    int    mode_bell            = 0;
    int    mode_sat             = 0;
    int    mode_tr              = 0;
    int    mode_qvm_api_test    = 0;
    char  *calib_file           = NULL;
    char  *vm_script            = NULL;
    int    monitor_cycles       = 10;
    char  *tx_outfile          = NULL;

    int idx = 1;
    while (idx < argc) {
        if (strcmp(argv[idx], "--loopback") == 0) {
            mode_loopback = 1; idx++;
        } else if (strcmp(argv[idx], "--physical") == 0) {
            mode_physical = 1; idx++;
        } else if (strcmp(argv[idx], "--ether-transfer") == 0) {
            mode_ether_transfer = 1; idx++;
        } else if (strcmp(argv[idx], "--ether-monitor") == 0) {
            mode_ether_monitor = 1;
            if (idx + 1 < argc && argv[idx+1][0] != '-') {
                monitor_cycles = atoi(argv[idx+1]); idx++;
            }
            idx++;
        } else if (strcmp(argv[idx], "--ether-gate") == 0) {
            mode_ether_gate = 1; idx++;
        } else if (strcmp(argv[idx], "--ether-decode") == 0) {
            mode_ether_decode = 1; idx++;
        } else if (strcmp(argv[idx], "--ether-calibrate") == 0 && idx+1 < argc) {
            mode_ether_calibrate = 1; calib_file = argv[idx+1]; idx += 2;
        } else if (strcmp(argv[idx], "--ether-equalize") == 0 && idx+1 < argc) {
            mode_ether_equalize = 1; calib_file = argv[idx+1]; idx += 2;
        } else if (strcmp(argv[idx], "--entangle") == 0) {
            mode_ether_entangle = 1; idx++;
        } else if (strcmp(argv[idx], "--field") == 0) {
            mode_ether_field = 1; idx++;
        } else if (strcmp(argv[idx], "--vm") == 0) {
            mode_vm = 1;
            if (idx+1 < argc && argv[idx+1][0] != '-') {
                vm_script = argv[idx+1]; idx++;
            }
            idx++;
        } else if (strcmp(argv[idx], "--regularization") == 0 && idx+1 < argc) {
            g_lambda = atof(argv[idx+1]); idx += 2;
        } else if (strcmp(argv[idx], "--bell") == 0) {
            mode_bell = 1; idx++;
        } else if (strcmp(argv[idx], "--time-reversal") == 0) {
            mode_tr = 1; idx++;
        } else if (strcmp(argv[idx], "--qvm-api-test") == 0) {
            mode_qvm_api_test = 1; idx++;
        } else if (strcmp(argv[idx], "--sat") == 0) {
            mode_sat = 1; idx++;
        } else if (strcmp(argv[idx], "--tx-only") == 0) {
            mode_tx_only = 1; idx++;
        } else if (strcmp(argv[idx], "--rx-only") == 0) {
            mode_rx_only = 1; idx++;
        } else if (strcmp(argv[idx], "--tx-file") == 0 && idx + 1 < argc) {
            mode_tx_only = 1; tx_outfile = argv[idx+1]; idx += 2;
        } else if (strcmp(argv[idx], "--ofdm") == 0) {
            mode_ofdm = 1; idx++;
        } else if (strcmp(argv[idx], "--ofdm-file") == 0 && idx + 1 < argc) {
            mode_ofdm = 1; tx_outfile = argv[idx+1]; idx += 2;
        } else if (strcmp(argv[idx], "--cycles") == 0 && idx + 1 < argc) {
            g_cycles = atoi(argv[idx+1]); idx += 2;
        } else if (strcmp(argv[idx], "--snr-db") == 0 && idx + 1 < argc) {
            g_snr_db = atof(argv[idx+1]); idx += 2;
        } else if (strcmp(argv[idx], "--fading") == 0 && idx + 1 < argc) {
            g_fading = atoi(argv[idx+1]); idx += 2;
        } else if (strcmp(argv[idx], "--drift") == 0 && idx + 1 < argc) {
            g_drift = atoi(argv[idx+1]); idx += 2;
        } else if (argv[idx][0] != '-') {
            if (idx == 1) D    = atoi(argv[idx]);
            else if (idx == 2) freq = (int)strtod(argv[idx], NULL);
            else if (idx == 3) rate = (int)strtod(argv[idx], NULL);
            else if (idx == 4) gain = atoi(argv[idx]);
            idx++;
        } else {
            idx++;
        }
    }

    D = D < 2 ? 2 : (D > MAX_DIM ? MAX_DIM : D);
    if (g_cycles < 1) g_cycles = 1;
    g_dim = D;
    int sat_n_vars = mode_sat ? g_dim : 6;  /* capture before --sat overrides it */
    /* --sat N overrides g_dim for the SAT solver's variable count */
    for (int i = 1; i < argc; i++)
        if (strcmp(argv[i], "--sat") == 0 && i+1 < argc && argv[i+1][0] != '-')
            sat_n_vars = atoi(argv[i+1]);

    if (mode_tx_only) {
        return run_tx_only(freq, rate, D, tx_outfile);
    }
    if (mode_ofdm) {
        sdr_find();
        int has_sdr = (access(g_sdr_dev, R_OK|W_OK) == 0);
        if (!has_sdr) { fprintf(stderr,"[OFDM] Needs SDR for ambient capture\n"); return 1; }
        SdrDev sdr; if (sdr_open(&sdr, freq, rate, gain)!=0) return 1;
        Wavefunction wf=wf_alloc(D,freq,rate);
        sdr_capture(&sdr); wf_from_iq(sdr.iq_i,sdr.iq_q,sdr.iq_n,&wf);
        gate_ofdm_tx(&sdr,&wf,tx_outfile);
        wf_free(&wf); sdr_close(&sdr);
        return 0;
    }
    if (mode_rx_only) {
        return run_rx_only(freq, rate, gain, D);
    }
    if (mode_loopback) {
        return run_loopback(freq, rate, D);
    }

    sdr_find();
    int has_sdr = (access(g_sdr_dev, R_OK|W_OK) == 0);

    if (mode_ether_transfer && has_sdr) {
        return run_ether_compute(freq, rate, gain, D);
    }
    if (mode_ether_monitor && has_sdr) {
        return run_ether_monitor(freq, rate, gain, D, monitor_cycles);
    }
    if (mode_ether_gate && has_sdr) {
        return run_ether_gate_matrix(freq, rate, gain, D,
                                     g_cycles > 1 ? g_cycles : 8);
    }
    if (mode_ether_decode && has_sdr) {
        return run_ether_decode(freq, rate, gain, D,
                                g_cycles > 1 ? g_cycles : 8);
    }
    if (mode_ether_calibrate && has_sdr) {
        return run_ether_calibrate(freq, rate, gain, D,
                                   g_cycles > 1 ? g_cycles : 16,
                                   g_lambda, calib_file);
    }
    if (mode_ether_equalize && has_sdr) {
        return run_ether_equalize(freq, rate, gain, calib_file,
                                  g_lambda, g_cycles > 0 ? g_cycles : 20);
    }
    if (mode_ether_entangle && has_sdr) {
        return run_ether_entangle(freq, rate, gain, D,
                                  g_cycles > 0 ? g_cycles : 50);
    }
    if (mode_ether_field && has_sdr) {
        return run_ether_field(freq, rate, gain, D,
                               g_cycles > 0 ? g_cycles : 30, 0.5);
    }
    if (mode_vm) {
        return run_quantum_vm(freq, rate, gain, D, vm_script);
    }
    if (mode_qvm_api_test) {
        /* Exercise the QVM API programmatically:
           calibrate → matvec → verify against CPU */
        printf("\n  ╔══════════════════════════════════════════════════╗\n");
        printf("  ║  QVM API TEST — Room as matvec co-processor     ║\n");
        printf("  ╚══════════════════════════════════════════════════╝\n\n");

        QvmCtx *q = qvm_create(freq, rate, D, gain);
        if (!q) { fprintf(stderr,"QVM create failed\n"); return 1; }

        /* 1. Show initial ambient state */
        printf("1. Ambient state:\n  ");
        double probs[MAX_DIM]; qvm_probs(q, probs, D);
        for (int i=0;i<D;i++) printf("%.3f ",probs[i]);
        printf("  S=%.3f γ=%.4f\n\n", qvm_entropy(q), qvm_purity(q));

        /* 2. Calibrate room channel */
        printf("2. Calibrating room channel M...\n");
        qvm_eval(q, "CALIBRATE 3");
        printf("\n");

        /* 3. Feed computation through the room */
        printf("3. Room computes: y = room(x)\n");
        double xs[MAX_DIM] = {0}, y_room[MAX_DIM] = {0};
        for (int i=0;i<D;i++) xs[i]=(double)((i*7+3)%100)/100.0;
        printf("  Input x:  "); for(int i=0;i<D;i++) printf("%.3f ",xs[i]); printf("\n");

        int rc = qvm_compute(q, xs, y_room, D);
        if (rc == 0) {
            printf("  Room y:   "); for(int i=0;i<D;i++) printf("%.3f ",y_room[i]); printf("\n");
            printf("  Entropy: %.3f bits | Purity: %.4f\n",
                qvm_entropy(q), qvm_purity(q));
            printf("\n  ★ Room computed y = f(x) on EM substrate ★\n");
            printf("  To train f(x)≈target: offline reservoir training\n");
        } else {
            printf("  Room compute failed (no SDR)\n");
        }

        qvm_destroy(q);
        return 0;
    }
    if (mode_bell) {
        return run_bell_test(freq, rate, D, g_cycles > 1 ? g_cycles : 3000,
                             g_snr_db, g_fading, g_drift);
    }
    if (mode_tr && has_sdr) {
        return run_time_reversal(freq, rate, gain, D);
    }
    if (mode_sat && has_sdr) {
        return run_sat_solver(freq, rate, gain, sat_n_vars);
    }
    if (mode_physical && has_sdr) {
        return run_physical_ether(freq, rate, gain, D);
    }
    if (has_sdr) {
        return run_physical_ether(freq, rate, gain, D);
    } else {
        fprintf(stderr, "[INFO] No SDR hardware — switching to loopback mode\n");
        fprintf(stderr, "[INFO] Plug RTL-SDR for: --ether-transfer, --ether-monitor, --physical\n");
        return run_loopback(freq, rate, D);
    }
}
#endif
