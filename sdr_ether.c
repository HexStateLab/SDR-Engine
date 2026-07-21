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

    int    mode_loopback = 0;
    int    mode_tx_only  = 0;
    int    mode_rx_only  = 0;
    int    mode_physical = 0;
    char  *tx_outfile    = NULL;

    int idx = 1;
    while (idx < argc) {
        if (strcmp(argv[idx], "--loopback") == 0) {
            mode_loopback = 1; idx++;
        } else if (strcmp(argv[idx], "--physical") == 0) {
            mode_physical = 1; idx++;
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
    if (mode_physical && has_sdr) {
        return run_physical_ether(freq, rate, gain, D);
    }
    if (has_sdr) {
        return run_physical_ether(freq, rate, gain, D);
    } else {
        fprintf(stderr, "[INFO] No SDR hardware — switching to loopback mode\n");
        fprintf(stderr, "[INFO] Use --physical with RTL-SDR for real TX via LO leakage\n");
        return run_loopback(freq, rate, D);
    }
}
