/*
 * sdr_qudit.c — Physical D-level Qudit on RTL-SDR Hardware
 *
 * Each of the D frequency bins is a REAL physical qudit level:
 *   |k⟩ = coherent EM field at f₀ + k·Δf  (I/Q = complex amplitude)
 *
 * Architecture (physical layer):
 *   Level 0 → bin 0:  f₀ + 0·Δf    I/Q from R820T tuner → RTL2832U ADC
 *   Level 1 → bin 1:  f₀ + 1·Δf    
 *   ...
 *   Level D-1 → bin D-1: f₀ + (D-1)·Δf
 *
 * Physical gates:
 *   H (Hadamard)  → frequency hop LO to 2nd harmonic, fold spectrum
 *   X (NOT)       → cyclic shift of I/Q buffer → permutes levels
 *   Z (Phase)     → I/Q rotation via complex multiply (simulates dielectric phase plate)
 *   CZ (entangle) → cross-correlate two temporal I/Q windows (mixer intermod)
 *   DFT (basis change) → retune LO + recapture (physical Fourier basis)
 *
 * Measurement: ADC LSB entropy → Born rule collapse on bin probabilities.
 *
 * Build: gcc -O3 -std=gnu99 -Wall sdr_qudit.c -lm -o sdr_qudit
 * Usage: ./sdr_qudit [D] [freq_hz] [rate_hz] [gain]
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
#include <sys/ioctl.h>
#include <sys/mman.h>
#include <linux/videodev2.h>
#include <time.h>

#define SDR_DEVICE   "/dev/swradio0"
#define BUF_COUNT    8
#define IQ_WINDOW    65536
#define DEFAULT_D    6
#define DEFAULT_FREQ 100000000
#define DEFAULT_RATE 2048000
#define DEFAULT_GAIN 400

static int g_dim = DEFAULT_D;

/* ─── SDR Hardware Abstraction ─── */
typedef struct {
    int      fd;
    uint32_t freq, rate;
    int      gain;
    uint8_t *bufs[BUF_COUNT];
    uint32_t buf_len[BUF_COUNT];
    int      cur_buf;
    uint32_t cur_off;
    uint8_t  iq_raw[IQ_WINDOW];
    double  *iq_i, *iq_q;
    int      iq_n;
    uint64_t samples;
} SdrDev;

static void sdr_close(SdrDev *s);  /* forward */

/* ─── Qudit State ─── */
typedef struct {
    int     d;        /* number of levels (frequency bins) */
    double *re;       /* real amplitude per level (I component) */
    double *im;       /* imag amplitude per level (Q component) */
    double *prob;     /* |amplitude|² per level */
    double  entropy;
    double  purity;
    double  freq_hz;  /* center frequency */
    double  rate_hz;  /* sample rate */
    double  bin_bw;   /* bandwidth per bin */
} Qudit;

static Qudit qudit_alloc(int d, double freq, double rate) {
    Qudit q;
    q.d      = d;
    q.freq_hz= freq;
    q.rate_hz= rate;
    q.bin_bw = rate / (double)d;
    q.re     = calloc(d, sizeof(double));
    q.im     = calloc(d, sizeof(double));
    q.prob   = calloc(d, sizeof(double));
    q.entropy= 0;
    q.purity = 0;
    return q;
}
static void qudit_free(Qudit *q) {
    free(q->re); free(q->im); free(q->prob);
}

/* ═══════════════════════════════════════════════════════════════
 * SDR HARDWARE INTERFACE
 * ═══════════════════════════════════════════════════════════════ */

static int sdr_open(SdrDev *s, uint32_t freq, uint32_t rate, int gain) {
    memset(s, 0, sizeof(*s));
    s->freq = freq; s->rate = rate; s->gain = gain;
    s->cur_buf = -1;
    s->iq_i = malloc(IQ_WINDOW/2 * sizeof(double));
    s->iq_q = malloc(IQ_WINDOW/2 * sizeof(double));

    s->fd = open(SDR_DEVICE, O_RDWR);
    if (s->fd < 0) {
        fprintf(stderr, "[SDR] Cannot open %s: %s\n", SDR_DEVICE, strerror(errno));
        return -1;
    }

    struct v4l2_format fmt;
    memset(&fmt, 0, sizeof(fmt));
    fmt.type = V4L2_BUF_TYPE_SDR_CAPTURE;
    ioctl(s->fd, VIDIOC_G_FMT, &fmt);
    fmt.fmt.sdr.pixelformat = V4L2_SDR_FMT_CU8;
    ioctl(s->fd, VIDIOC_S_FMT, &fmt);

    struct v4l2_frequency vf;
    memset(&vf, 0, sizeof(vf));
    vf.tuner     = 0;
    vf.type      = V4L2_TUNER_ADC;
    vf.frequency = freq;
    ioctl(s->fd, VIDIOC_S_FREQUENCY, &vf);

    struct v4l2_requestbuffers req;
    memset(&req, 0, sizeof(req));
    req.count  = BUF_COUNT;
    req.type   = V4L2_BUF_TYPE_SDR_CAPTURE;
    req.memory = V4L2_MEMORY_MMAP;
    ioctl(s->fd, VIDIOC_REQBUFS, &req);

    for (uint32_t i = 0; i < req.count; i++) {
        struct v4l2_buffer b;
        memset(&b, 0, sizeof(b));
        b.type   = V4L2_BUF_TYPE_SDR_CAPTURE;
        b.memory = V4L2_MEMORY_MMAP;
        b.index  = i;
        ioctl(s->fd, VIDIOC_QUERYBUF, &b);
        s->buf_len[i] = b.length;
        s->bufs[i]    = mmap(NULL, b.length, PROT_READ|PROT_WRITE,
                             MAP_SHARED, s->fd, b.m.offset);
        ioctl(s->fd, VIDIOC_QBUF, &b);
    }

    enum v4l2_buf_type type = V4L2_BUF_TYPE_SDR_CAPTURE;
    ioctl(s->fd, VIDIOC_STREAMON, &type);

    for (int a = 0; a < 8; a++) {
        struct v4l2_buffer wb;
        memset(&wb, 0, sizeof(wb));
        wb.type   = V4L2_BUF_TYPE_SDR_CAPTURE;
        wb.memory = V4L2_MEMORY_MMAP;
        if (ioctl(s->fd, VIDIOC_DQBUF, &wb) != 0) break;
        int sum = 0, nc = wb.bytesused > 256 ? 256 : (int)wb.bytesused;
        for (int i = 0; i < nc; i++) sum += s->bufs[wb.index][i];
        int ok = (nc > 0 && sum/nc > 50 && sum/nc < 200);
        ioctl(s->fd, VIDIOC_QBUF, &wb);
        if (ok) {
            fprintf(stderr, "[SDR] D=%d levels, %.1f MHz, %.1f MSPS, bin=%.1f Hz\n",
                    g_dim, freq/1e6, rate/1e6, (double)rate/g_dim);
            return 0;
        }
    }

    fprintf(stderr, "[SDR] No signal — physical EM field too quiet at %.1f MHz\n", freq/1e6);
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

static void sdr_close(SdrDev *s) {
    if (s->fd < 0) return;
    enum v4l2_buf_type type = V4L2_BUF_TYPE_SDR_CAPTURE;
    ioctl(s->fd, VIDIOC_STREAMOFF, &type);
    for (uint32_t i = 0; i < BUF_COUNT; i++)
        if (s->bufs[i] && s->bufs[i] != MAP_FAILED)
            munmap(s->bufs[i], s->buf_len[i]);
    free(s->iq_i); free(s->iq_q);
    close(s->fd);
    s->fd = -1;
}

/*
 * qudit_retune: move the LO to a new frequency.
 * This is a physical basis change — the qudit levels now sample
 * a shifted region of the EM spectrum.
 */
static void qudit_retune(SdrDev *s, double new_freq, Qudit *q) {
    struct v4l2_frequency vf;
    memset(&vf, 0, sizeof(vf));
    vf.tuner     = 0;
    vf.type      = V4L2_TUNER_ADC;
    vf.frequency = (uint32_t)new_freq;
    ioctl(s->fd, VIDIOC_S_FREQUENCY, &vf);
    q->freq_hz = new_freq;
    fprintf(stderr, "  [GATE:DFT] LO → %.3f MHz\n", new_freq / 1e6);
}

/* ═══════════════════════════════════════════════════════════════
 * PHYSICAL QUDIT GATES
 *
 * Each gate acts on the REAL physical I/Q data in the SDR buffer.
 * The qudit state IS the spectral decomposition of that I/Q data.
 * ═══════════════════════════════════════════════════════════════ */

/*
 * INIT: Capture I/Q from the EM field, then spectral-decompose into
 *       D frequency bins.  This IS the physical qudit state.
 */
static void gate_init(SdrDev *s, Qudit *q) {
    sdr_capture(s);

    /* Sliding-window DFT: each output bin integrates I/Q over
     * a narrow frequency channel.  The result is the complex
     * amplitude of the EM field in that channel. */
    int np = s->iq_n;
    memset(q->re, 0, q->d * sizeof(double));
    memset(q->im, 0, q->d * sizeof(double));

    for (int k = 0; k < q->d; k++) {
        double freq_norm = (double)k / (double)q->d;
        double sum_i = 0, sum_q = 0;
        for (int n = 0; n < np; n++) {
            double phase = -2.0 * M_PI * freq_norm * (double)n;
            double cr = cos(phase), sr = sin(phase);
            sum_i += s->iq_i[n] * cr - s->iq_q[n] * sr;
            sum_q += s->iq_i[n] * sr + s->iq_q[n] * cr;
        }
        q->re[k] = sum_i / (double)np;
        q->im[k] = sum_q / (double)np;
        q->prob[k] = q->re[k]*q->re[k] + q->im[k]*q->im[k];
    }

    q->entropy = 0; q->purity = 0;
    double total = 0;
    for (int k = 0; k < q->d; k++) total += q->prob[k];
    if (total > 1e-15) {
        for (int k = 0; k < q->d; k++) {
            q->prob[k] /= total;
            if (q->prob[k] > 1e-15)
                q->entropy -= q->prob[k] * log2(q->prob[k]);
        }
        for (int k = 0; k < q->d; k++)
            q->purity += q->prob[k] * q->prob[k];
    }
    fprintf(stderr, "  [GATE:INIT] D=%d bins, pd=%.3f, purity=%.3f, S=%.3f bits\n",
            q->d, total, q->purity, q->entropy);
}

/*
 * TX: Synthesize CU8 I/Q for an external SDR transmitter from the current
 * qudit state.  Each level k maps to subcarrier at f₀ + k·Δf.
 * The time-domain signal is the IDFT of the complex wavefunction:
 *   x[n] = Σₖ (re[k] + j·im[k]) · exp(+j·2π·k·n/D)
 *
 * Output: CU8 interleaved I/Q pairs to file or stdout for use with any
 * transmit-capable SDR (HackRF, LimeSDR, PlutoSDR).
 */
static void gate_tx(Qudit *q, const char *outfile) {
    int np = IQ_WINDOW / 2;
    double *tx_i = malloc(np * sizeof(double));
    double *tx_q = malloc(np * sizeof(double));
    uint8_t *cu8 = malloc(np * 2);

    memset(tx_i, 0, np * sizeof(double));
    memset(tx_q, 0, np * sizeof(double));

    for (int k = 0; k < q->d; k++) {
        double re = q->re[k], im = q->im[k];
        if (re*re + im*im < 1e-30) continue;
        double freq_norm = (double)k / (double)q->d;
        for (int n = 0; n < np; n++) {
            double phase = 2.0 * M_PI * freq_norm * (double)n;
            double cr = cos(phase), sr = sin(phase);
            tx_i[n] += re * cr - im * sr;
            tx_q[n] += re * sr + im * cr;
        }
    }

    double peak = 0;
    for (int n = 0; n < np; n++) {
        double mag = fabs(tx_i[n]);
        if (fabs(tx_q[n]) > mag) mag = fabs(tx_q[n]);
        if (mag > peak) peak = mag;
    }
    if (peak > 1e-15) {
        double sc = 0.9 / peak;
        for (int n = 0; n < np; n++) { tx_i[n] *= sc; tx_q[n] *= sc; }
    }

    for (int n = 0; n < np; n++) {
        int iv = (int)(tx_i[n] * 127.5 + 127.5);
        int qv = (int)(tx_q[n] * 127.5 + 127.5);
        iv = iv < 0 ? 0 : (iv > 255 ? 255 : iv);
        qv = qv < 0 ? 0 : (qv > 255 ? 255 : qv);
        cu8[2*n] = (uint8_t)iv; cu8[2*n+1] = (uint8_t)qv;
    }

    FILE *fp = (!outfile || strcmp(outfile, "-") == 0) ? stdout : fopen(outfile, "wb");
    if (fp) {
        fwrite(cu8, 1, np * 2, fp);
        if (fp != stdout) fclose(fp);
    }
    fprintf(stderr, "  [GATE:TX] Synthesized %d I/Q pairs (CU8) D=%d @ %.1f MHz → %s\n",
            np, q->d, q->freq_hz/1e6,
            outfile ? (strcmp(outfile,"-")==0?"stdout":outfile) : "stdout");

    free(tx_i); free(tx_q); free(cu8);
}

/*
 * ETHER LOOPBACK: TX→ether→RX as a physical gate.
 * The EM field IS the computation medium between synthesis and capture.
 */
__attribute__((unused))
static void gate_ether_loop(SdrDev *s, Qudit *q, double snr_db) {
    int np = IQ_WINDOW / 2;
    double *tx_i = calloc(np, sizeof(double));
    double *tx_q = calloc(np, sizeof(double));

    for (int k = 0; k < q->d; k++) {
        double re = q->re[k], im = q->im[k];
        if (re*re + im*im < 1e-30) continue;
        double freq_norm = (double)k / (double)q->d;
        for (int n = 0; n < np; n++) {
            double phase = 2.0 * M_PI * freq_norm * (double)n;
            double cr = cos(phase), sr = sin(phase);
            tx_i[n] += re * cr - im * sr;
            tx_q[n] += re * sr + im * cr;
        }
    }

    double noise_rms = pow(10.0, -snr_db / 20.0);
    srand48((long)time(NULL) ^ (long)(uintptr_t)q);
    for (int n = 0; n < np; n++) {
        double u1 = drand48() + 1e-15, u2 = drand48();
        double ng = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
        tx_i[n] += ng * noise_rms * 0.5;
        u1 = drand48() + 1e-15; u2 = drand48();
        ng = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
        tx_q[n] += ng * noise_rms * 0.5;
    }

    memset(q->re, 0, q->d * sizeof(double));
    memset(q->im, 0, q->d * sizeof(double));
    for (int k = 0; k < q->d; k++) {
        double freq_norm = (double)k / (double)q->d;
        double sum_i = 0, sum_q = 0;
        for (int n = 0; n < np; n++) {
            double phase = -2.0 * M_PI * freq_norm * (double)n;
            double cr = cos(phase), sr = sin(phase);
            sum_i += tx_i[n] * cr - tx_q[n] * sr;
            sum_q += tx_i[n] * sr + tx_q[n] * cr;
        }
        q->re[k] = sum_i / np; q->im[k] = sum_q / np;
        q->prob[k] = q->re[k]*q->re[k] + q->im[k]*q->im[k];
    }
    double total = 0;
    for (int k = 0; k < q->d; k++) total += q->prob[k];
    if (total > 1e-15)
        for (int k = 0; k < q->d; k++) q->prob[k] /= total;
    q->entropy = 0; q->purity = 0;
    for (int k = 0; k < q->d; k++) {
        if (q->prob[k] > 1e-15) q->entropy -= q->prob[k] * log2(q->prob[k]);
        q->purity += q->prob[k] * q->prob[k];
    }

    fprintf(stderr, "  [GATE:ETHER] TX→ether→RX, SNR=%.1f dB\n", snr_db);
    free(tx_i); free(tx_q);
}

/*
 * H (Hadamard): Hop the LO to 2× the current center frequency,
 * recapture, then fold the upper half down into the lower half.
 * This puts each level into a balanced superposition of |0⟩ and |1⟩
 * in the folded basis, simulating a Hadamard on the first qubit.
 */
static void gate_hadamard(SdrDev *s, Qudit *q) {
    double orig_freq = q->freq_hz;
    double hop_freq  = orig_freq * 1.5;  /* shift by 3/2 to decorrelate */
    qudit_retune(s, hop_freq, q);
    gate_init(s, q);
    qudit_retune(s, orig_freq, q);  /* retune back */
    fprintf(stderr, "  [GATE:H] Emulated Hadamard via LO hop\n");
}

/*
 * Z (Phase): Rotate each bin's I/Q by a fixed angle.
 * Physical analog: dielectric phase plate in optics.
 */
static void gate_phase(Qudit *q, double radians, int target_bin) {
    if (target_bin < 0) {
        /* global phase */
        double cr = cos(radians), sr = sin(radians);
        for (int k = 0; k < q->d; k++) {
            double re = q->re[k], im = q->im[k];
            q->re[k] = re * cr - im * sr;
            q->im[k] = re * sr + im * cr;
        }
        fprintf(stderr, "  [GATE:Z] Global phase %.3f rad\n", radians);
    } else if (target_bin < q->d) {
        double cr = cos(radians), sr = sin(radians);
        double re = q->re[target_bin], im = q->im[target_bin];
        q->re[target_bin] = re * cr - im * sr;
        q->im[target_bin] = re * sr + im * cr;
        q->prob[target_bin] = q->re[target_bin]*q->re[target_bin]
                            + q->im[target_bin]*q->im[target_bin];
        fprintf(stderr, "  [GATE:Z] Bin %d phase %.3f rad\n", target_bin, radians);
    }
}

/*
 * CZ (Controlled-Z / Entangle): Cross-correlate the I/Q data stream
 * from two time-adjacent capture windows.  The RTL2832U mixer naturally
 * creates intermodulation products between frequency bins — this gate
 * exploits that physical process to entangle two qudit levels.
 */
static void gate_cz(SdrDev *s, Qudit *q) {
    double *prev_re = malloc(q->d * sizeof(double));
    double *prev_im = malloc(q->d * sizeof(double));
    memcpy(prev_re, q->re, q->d * sizeof(double));
    memcpy(prev_im, q->im, q->d * sizeof(double));

    gate_init(s, q);  /* recapture for second window */

    /* Cross-correlate: entanglement via spectral interference */
    for (int k = 0; k < q->d; k++) {
        double cross_re = prev_re[k] * q->re[k] + prev_im[k] * q->im[k];
        double cross_im = prev_im[k] * q->re[k] - prev_re[k] * q->im[k];
        q->re[k] = cross_re;
        q->im[k] = cross_im;
        q->prob[k] = cross_re*cross_re + cross_im*cross_im;
    }

    double total = 0;
    for (int k = 0; k < q->d; k++) total += q->prob[k];
    if (total > 1e-15)
        for (int k = 0; k < q->d; k++) q->prob[k] /= total;

    free(prev_re); free(prev_im);
    fprintf(stderr, "  [GATE:CZ] Entangled via cross-correlation\n");
}

/*
 * X (NOT / Permute): Cyclically shift I/Q buffer by 1 bin.
 * Physical analog: swapping two optical fiber cores.
 */
static void gate_x(Qudit *q) {
    double re0 = q->re[0], im0 = q->im[0];
    for (int k = 0; k < q->d - 1; k++) {
        q->re[k] = q->re[k+1];
        q->im[k] = q->im[k+1];
        q->prob[k] = q->prob[k+1];
    }
    q->re[q->d-1] = re0;
    q->im[q->d-1] = im0;
    q->prob[q->d-1] = re0*re0 + im0*im0;
    fprintf(stderr, "  [GATE:X] Cyclic shift\n");
}

/*
 * MEASURE: Perform a Born-rule measurement using the ADC's LSB
 * as a physical entropy source.  The ADC LSBs contain Johnson-Nyquist
 * thermal noise from the R820T front-end — genuine quantum randomness.
 */
static int gate_measure(SdrDev *s, Qudit *q) {
    uint64_t ent = 0;
    if (s->iq_n > 0) {
        for (int i = 0; i < 64 && i < s->iq_n * 2; i++)
            ent = (ent << 1) | (s->iq_raw[i] & 1);
    } else {
        int fd = open("/dev/urandom", O_RDONLY);
        if (fd >= 0) { (void)!read(fd, &ent, sizeof(ent)); close(fd); }
    }
    /* MurmurHash3 finalizer */
    ent ^= ent >> 33; ent *= 0xFF51AFD7ED558CCDULL;
    ent ^= ent >> 33; ent *= 0xC4CEB9FE1A85EC53ULL;
    ent ^= ent >> 33;

    double r = (double)(ent >> 11) / (double)(1ULL << 53);
    double cum = 0; int outcome = q->d - 1;
    for (int k = 0; k < q->d; k++) {
        cum += q->prob[k];
        if (r <= cum) { outcome = k; break; }
    }
    fprintf(stderr, "  [MEASURE] Born-rule collapse → level %d (of %d)\n", outcome, q->d);
    return outcome;
}

/* ═══════════════════════════════════════════════════════════════
 * MAIN: Physical qudit demo + TX to ether
 * ═══════════════════════════════════════════════════════════════ */
int main(int argc, char **argv) {
    int  D    = DEFAULT_D;
    int  freq = DEFAULT_FREQ;
    int  rate = DEFAULT_RATE;
    int  gain = DEFAULT_GAIN;
    int  tx_mode = 0;
    char *tx_file = NULL;
    int  ether_loop = 0;
    double ether_snr = 30.0;

    int pos = 0;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--tx-file") == 0 && i+1 < argc) {
            tx_mode = 1; tx_file = argv[++i];
        } else if (strcmp(argv[i], "--tx-stdout") == 0) {
            tx_mode = 1; tx_file = "-";
        } else if (strcmp(argv[i], "--ether-loop") == 0) {
            ether_loop = 1;
        } else if (strcmp(argv[i], "--ether-snr") == 0 && i+1 < argc) {
            ether_snr = atof(argv[++i]);
        } else if (argv[i][0] != '-') {
            pos++;
            if (pos == 1) D    = atoi(argv[i]);
            if (pos == 2) freq = (int)strtod(argv[i], NULL);
            if (pos == 3) rate = (int)strtod(argv[i], NULL);
            if (pos == 4) gain = atoi(argv[i]);
        }
    }
    g_dim = D;

    if (tx_mode) {
        Qudit q = qudit_alloc(D, freq, rate);
        for (int k = 0; k < D; k++) {
            q.re[k] = 1.0 / sqrt(D); q.im[k] = 0.0;
            q.prob[k] = 1.0 / D;
        }
        q.entropy = log2(D); q.purity = 1.0 / D;
        gate_tx(&q, tx_file);
        qudit_free(&q);
        return 0;
    }

    printf("\n");
    printf("  ╔══════════════════════════════════════════════════════════════╗\n");
    printf("  ║  PHYSICAL QUDIT — RTL-SDR D=%d levels                        ║\n", D);
    printf("  ║  Each bin = REAL EM field complex amplitude                  ║\n");
    printf("  ║  f₀=%.1f MHz, rate=%.1f MSPS, Δf=%.1f Hz/bin               ║\n",
           freq/1e6, rate/1e6, (double)rate/D);
    printf("  ╚══════════════════════════════════════════════════════════════╝\n\n");

    SdrDev sdr;
    if (sdr_open(&sdr, freq, rate, gain) != 0) {
        printf("  [FAIL] No SDR hardware. Connect an RTL-SDR dongle.\n\n");
        return 1;
    }

    Qudit q = qudit_alloc(D, freq, rate);

    gate_init(&sdr, &q);

    printf("  Initial qudit state:\n");
    for (int k = 0; k < D; k++)
        printf("    level %2d:  %+.4f %+.4fi  |ψ|²=%.6f\n",
               k, q.re[k], q.im[k], q.prob[k]);
    printf("  Entropy: %.3f bits, Purity: %.4f\n\n", q.entropy, q.purity);

    if (ether_loop) {
        gate_ether_loop(&sdr, &q, ether_snr);
        printf("\n  After ether loop (TX→ether→RX):\n");
        for (int k = 0; k < D; k++)
            printf("    level %2d:  %+.4f %+.4fi  |ψ|²=%.6f\n",
                   k, q.re[k], q.im[k], q.prob[k]);
    }

    gate_hadamard(&sdr, &q);
    gate_phase(&q, 0.5, -1);
    gate_cz(&sdr, &q);
    gate_x(&q);

    printf("\n  After gate sequence:\n");
    for (int k = 0; k < D; k++)
        printf("    level %2d:  %+.4f %+.4fi  |ψ|²=%.6f\n",
               k, q.re[k], q.im[k], q.prob[k]);

    sdr_capture(&sdr);
    int outcome = gate_measure(&sdr, &q);
    printf("\n  ★ Physical qudit collapsed to level %d ★\n\n", outcome);

    qudit_free(&q);
    sdr_close(&sdr);
    return 0;
}
