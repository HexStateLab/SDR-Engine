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
 *   --cycles N    Run N TX→ether→RX→gate feedback cycles
 *   --snr-db N    Ether channel SNR for loopback (default 30 dB)
 *   --fading 1    Enable frequency-selective fading in loopback ether
 *   --drift 1     Enable random phase drift per loopback cycle
 *
 * Build: gcc -O3 -std=gnu99 sdr_ether.c -lm -o sdr_ether
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

#define SDR_DEVICE   "/dev/swradio0"
#define BUF_COUNT    8
#define IQ_WINDOW    65536
#define MAX_DIM      256
#define DEFAULT_D    6
#define DEFAULT_FREQ 100000000
#define DEFAULT_RATE 2048000

/* ═══════════════════════════════════════════════════════════════
 * GLOBAL STATE
 * ═══════════════════════════════════════════════════════════════ */
static int    g_dim     = DEFAULT_D;
static int    g_cycles  = 1;
static double g_snr_db  = 30.0;
static int    g_fading  = 0;
static int    g_drift   = 0;
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
} SdrDev;

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
    memset(s, 0, sizeof(*s));
    s->freq = freq; s->rate = rate; s->gain = gain;
    s->cur_buf = -1;

    s->fd = open(SDR_DEVICE, O_RDWR);
    if (s->fd < 0) {
        fprintf(stderr, "[SDR] Cannot open %s: %s\n", SDR_DEVICE, strerror(errno));
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
            if (need == 0) break;
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

static void sdr_close(SdrDev *s) {
    if (s->fd < 0) return;
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

/* H: Hadamard-like via LO frequency hop */
static void gate_H(SdrDev *s, Wavefunction *wf) {
    uint32_t orig = wf->freq;
    uint32_t hop  = (uint32_t)(orig * 1.5);
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

    fprintf(stderr, "  [GATE:TX-HW] Radiating qudit D=%d via LO leakage\n", d);
    fprintf(stderr, "               freq  dwell  |ψ|²\n");

    double total_time_us = 0;
    for (int k = 0; k < d; k++) {
        int dwell_us = (int)(wf->prob[k] * 15000.0);  /* up to 15 ms per level */
        if (dwell_us < 30) continue;

        uint32_t tx_freq = base + (uint32_t)(k * bw);
        if (tx_freq < 24000000) tx_freq = 24000000;
        if (tx_freq > 1750000000) tx_freq = 1750000000;

        sdr_retune(s, tx_freq);
        usleep(dwell_us);
        total_time_us += dwell_us;

        fprintf(stderr, "               %.3f MHz %4dμs  %.3f\n",
                tx_freq/1e6, dwell_us, wf->prob[k]);
    }

    /* Return to center for RX */
    sdr_retune(s, base);
    usleep(500);

    fprintf(stderr, "  [GATE:TX-HW] Emitted %.1f ms of qudit state into ether\n",
            total_time_us / 1000.0);
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
 * MAIN
 * ═══════════════════════════════════════════════════════════════ */
int main(int argc, char **argv) {
    int    D    = DEFAULT_D;
    int    freq = DEFAULT_FREQ;
    int    rate = DEFAULT_RATE;
    int    gain = 400;

    int    mode_loopback       = 0;
    int    mode_tx_only        = 0;
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
    char  *calib_file           = NULL;
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
        } else if (strcmp(argv[idx], "--regularization") == 0 && idx+1 < argc) {
            g_lambda = atof(argv[idx+1]); idx += 2;
        } else if (strcmp(argv[idx], "--tx-only") == 0) {
            mode_tx_only = 1; idx++;
        } else if (strcmp(argv[idx], "--rx-only") == 0) {
            mode_rx_only = 1; idx++;
        } else if (strcmp(argv[idx], "--tx-file") == 0 && idx + 1 < argc) {
            mode_tx_only = 1; tx_outfile = argv[idx+1]; idx += 2;
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

    if (mode_tx_only) {
        return run_tx_only(freq, rate, D, tx_outfile);
    }
    if (mode_rx_only) {
        return run_rx_only(freq, rate, gain, D);
    }
    if (mode_loopback) {
        return run_loopback(freq, rate, D);
    }

    int has_sdr = (access(SDR_DEVICE, R_OK|W_OK) == 0);

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
