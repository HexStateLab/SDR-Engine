/*
 * sdr_quantum_engine.c — RTL-SDR Quantum Engine (runtime-D)
 *
 * Usage: ./sdr_qe [D] [freq_hz] [rate_hz] [gain]
 *        ./sdr_qe 6                          # D=6, defaults
 *        ./sdr_qe 12 100e6 2048000 400       # D=12, custom
 *        ./sdr_qe 256 50e6                   # D=256, 50 MHz
 *
 * Build: gcc -O3 -std=gnu99 sdr_quantum_engine.c -lm -o sdr_qe
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
#define IQ_WINDOW    8192    /* 4096 I/Q pairs — 16/bin at D=256 */
#define MAX_DIM      256

/* ═══════════════════════════════════════════
 * GLOBAL STATE — runtime dimension
 * ═══════════════════════════════════════════ */
static int g_dim = 6;

/* ═══════════════════════════════════════════
 * SDR DEVICE
 * ═══════════════════════════════════════════ */
typedef struct {
    int fd; uint32_t freq, rate; int gain;
    uint8_t **bufs; uint32_t buf_len[BUF_COUNT];
    int cur_buf; uint32_t cur_off; uint64_t samples;
    uint8_t iq_raw[IQ_WINDOW];
    double iq_i[IQ_WINDOW/2], iq_q[IQ_WINDOW/2];
    int iq_n;
} SdrDev;

/* ═══════════════════════════════════════════
 * WAVEFUNCTION
 * ═══════════════════════════════════════════ */
typedef struct {
    double *re, *im, *prob;
    double entropy, purity;
} Wavefunction;

static Wavefunction wf_alloc(int dim) {
    Wavefunction w;
    w.re = calloc(dim, sizeof(double));
    w.im = calloc(dim, sizeof(double));
    w.prob = calloc(dim, sizeof(double));
    w.entropy = 0; w.purity = 0;
    return w;
}
static void wf_free(Wavefunction *w) { free(w->re); free(w->im); free(w->prob); }

/* ═══════════════════════════════════════════
 * SDR LIFECYCLE
 * ═══════════════════════════════════════════ */
static int sdr_open(SdrDev *s, uint32_t f, uint32_t r, int g) {
    memset(s, 0, sizeof(*s)); s->fd = -1; s->freq = f; s->rate = r; s->gain = g;
    s->cur_buf = -1;
    s->fd = open(SDR_DEVICE, O_RDWR);  /* blocking — wait for USB data */
    if (s->fd < 0) { fprintf(stderr, "[SDR] open: %s\n", strerror(errno)); return -1; }

    struct v4l2_format fmt; memset(&fmt, 0, sizeof(fmt));
    fmt.type = V4L2_BUF_TYPE_SDR_CAPTURE;
    if (ioctl(s->fd, VIDIOC_G_FMT, &fmt) < 0) goto fail;
    fmt.fmt.sdr.pixelformat = V4L2_SDR_FMT_CU8;
    if (ioctl(s->fd, VIDIOC_S_FMT, &fmt) < 0) goto fail;

    { struct v4l2_frequency vf; memset(&vf,0,sizeof(vf));
      vf.tuner=0; vf.type=V4L2_TUNER_ADC; vf.frequency=f; ioctl(s->fd,VIDIOC_S_FREQUENCY,&vf); }

    struct v4l2_requestbuffers req; memset(&req,0,sizeof(req));
    req.count=BUF_COUNT; req.type=V4L2_BUF_TYPE_SDR_CAPTURE; req.memory=V4L2_MEMORY_MMAP;
    if (ioctl(s->fd,VIDIOC_REQBUFS,&req)<0) goto fail;

    s->bufs = calloc(req.count, sizeof(uint8_t*));
    for (uint32_t i=0; i<req.count; i++) {
        struct v4l2_buffer b; memset(&b,0,sizeof(b));
        b.type=V4L2_BUF_TYPE_SDR_CAPTURE; b.memory=V4L2_MEMORY_MMAP; b.index=i;
        if (ioctl(s->fd,VIDIOC_QUERYBUF,&b)<0) goto fail;
        s->buf_len[i]=b.length;
        s->bufs[i]=mmap(NULL,b.length,PROT_READ|PROT_WRITE,MAP_SHARED,s->fd,b.m.offset);
    }
    for (uint32_t i=0; i<req.count; i++) {
        struct v4l2_buffer b; memset(&b,0,sizeof(b));
        b.type=V4L2_BUF_TYPE_SDR_CAPTURE; b.memory=V4L2_MEMORY_MMAP; b.index=i;
        ioctl(s->fd,VIDIOC_QBUF,&b);
    }

    { enum v4l2_buf_type t=V4L2_BUF_TYPE_SDR_CAPTURE;
      if (ioctl(s->fd,VIDIOC_STREAMON,&t)<0) goto fail; }

    /* Warmup: drain until USB data flows */
    for (int a=0; a<16; a++) {
        struct v4l2_buffer wb; memset(&wb,0,sizeof(wb));
        wb.type=V4L2_BUF_TYPE_SDR_CAPTURE; wb.memory=V4L2_MEMORY_MMAP;
        if (ioctl(s->fd,VIDIOC_DQBUF,&wb)!=0) break;
        int sum=0, nc=wb.bytesused>256?256:(int)wb.bytesused;
        for (int i=0;i<nc;i++) sum+=s->bufs[wb.index][i];
        int ok=(sum/nc>50&&sum/nc<200);
        ioctl(s->fd,VIDIOC_QBUF,&wb);
        if (ok) { fprintf(stderr,"[SDR] warm after %d DQBUFs\n",a+1); break; }
    }
    fprintf(stderr,"[SDR] D=%d @ %.1f MHz %.3f MSPS g=%d\n",g_dim,f/1e6,r/1e6,g);
    return 0;
fail:
    fprintf(stderr,"[SDR] init fail: %s\n",strerror(errno));
    close(s->fd); s->fd=-1; return -1;
}

static int sdr_read(SdrDev *s) {
    int copied=0, need=IQ_WINDOW;
    while (need>0) {
        if (s->cur_buf>=0 && s->cur_off < s->buf_len[s->cur_buf]) {
            uint32_t av=s->buf_len[s->cur_buf]-s->cur_off;
            int cp=(int)av<need?(int)av:need;
            memcpy(s->iq_raw+copied, s->bufs[s->cur_buf]+s->cur_off, cp);
            s->cur_off+=cp; copied+=cp; need-=cp;
            if (need==0) break;
            struct v4l2_buffer b; memset(&b,0,sizeof(b));
            b.type=V4L2_BUF_TYPE_SDR_CAPTURE; b.memory=V4L2_MEMORY_MMAP;
            b.index=s->cur_buf; ioctl(s->fd,VIDIOC_QBUF,&b);
            s->cur_buf=-1;
        }
        struct pollfd p={.fd=s->fd,.events=POLLIN};
        if (poll(&p,1,200)<=0) break;
        struct v4l2_buffer db; memset(&db,0,sizeof(db));
        db.type=V4L2_BUF_TYPE_SDR_CAPTURE; db.memory=V4L2_MEMORY_MMAP;
        if (ioctl(s->fd,VIDIOC_DQBUF,&db)<0){if(errno==EAGAIN){usleep(5000);continue;}break;}
        s->cur_buf=(int)db.index; s->cur_off=0;
    }
    int np=copied/2; s->iq_n=np; s->samples+=copied;
    for (int i=0;i<np;i++) {
        s->iq_i[i]=((double)s->iq_raw[2*i]-127.5)/128.0;
        s->iq_q[i]=((double)s->iq_raw[2*i+1]-127.5)/128.0;
    }
    return copied;
}

static void sdr_retune(SdrDev *s, uint32_t hz) {
    struct v4l2_frequency vf; memset(&vf,0,sizeof(vf));
    vf.tuner=0; vf.type=V4L2_TUNER_ADC; vf.frequency=hz;
    ioctl(s->fd,VIDIOC_S_FREQUENCY,&vf); s->freq=hz;
}

static void sdr_close(SdrDev *s) {
    if (s->fd<0) return;
    enum v4l2_buf_type t=V4L2_BUF_TYPE_SDR_CAPTURE; ioctl(s->fd,VIDIOC_STREAMOFF,&t);
    for (int i=0;i<BUF_COUNT;i++) if (s->bufs&&s->bufs[i]&&s->bufs[i]!=MAP_FAILED)
        munmap(s->bufs[i],s->buf_len[i]);
    free(s->bufs); close(s->fd); s->fd=-1;
}

/* ═══════════════════════════════════════════
 * WAVEFUNCTION
 * ═══════════════════════════════════════════ */
static void wf_compute(SdrDev *s, Wavefunction *wf) {
    memset(wf->re,0,g_dim*8); memset(wf->im,0,g_dim*8);
    memset(wf->prob,0,g_dim*8); wf->entropy=0; wf->purity=0;
    sdr_read(s); int np=s->iq_n;
    if (np<g_dim) { wf->prob[0]=1.0; wf->re[0]=1.0; wf->purity=1.0; return; }
    int blk=np/g_dim;
    for (int k=0;k<g_dim;k++) {
        int st=k*blk, en=(k+1)*blk;
        for (int i=st;i<en&&i<np;i++) { wf->re[k]+=s->iq_i[i]; wf->im[k]+=s->iq_q[i]; }
    }
    double total=0;
    for (int k=0;k<g_dim;k++) total+=wf->re[k]*wf->re[k]+wf->im[k]*wf->im[k];
    if (total<1e-30) { wf->prob[0]=1.0; wf->re[0]=1.0; wf->purity=1.0; return; }
    double sc=1.0/sqrt(total);
    for (int k=0;k<g_dim;k++) {
        wf->re[k]*=sc; wf->im[k]*=sc;
        wf->prob[k]=wf->re[k]*wf->re[k]+wf->im[k]*wf->im[k];
        if (wf->prob[k]>1e-14) wf->entropy-=wf->prob[k]*log2(wf->prob[k]);
        wf->purity+=wf->prob[k]*wf->prob[k];
    }
}

static void wf_print(const Wavefunction *wf, const char *label) {
    printf("  %-8s [",label);
    int show=g_dim<12?g_dim:8;
    for (int k=0;k<show;k++) printf("%.3f ",wf->prob[k]);
    if (g_dim>show) printf("...");
    printf("]  H=%.3f γ=%.4f\n",wf->entropy,wf->purity);
}

/* ═══════════════════════════════════════════
 * PHYSICAL GATES — use runtime g_dim
 * ═══════════════════════════════════════════ */
static void gate_H(SdrDev *s) {
    uint32_t o=s->freq; sdr_retune(s,o*2); usleep(3000); sdr_retune(s,o); usleep(3000);
    fprintf(stderr,"  [H] %.1f→%.1f→%.1f MHz\n",o/1e6,o*2/1e6,o/1e6);
}
static void gate_DFT(SdrDev *s, int step) {
    uint32_t n=s->freq+(uint32_t)step*(s->rate/g_dim);
    sdr_retune(s,n); usleep(3000);
    fprintf(stderr,"  [DFT] %.3f MHz\n",n/1e6);
}
static void gate_Phase(SdrDev *s, int ch, double rad) {
    sdr_read(s); if (s->iq_n<g_dim) return;
    int blk=s->iq_n/g_dim; double cr=cos(rad),sr=sin(rad);
    int st=ch*blk, en=(ch+1)*blk;
    for (int i=st;i<en&&i<s->iq_n;i++) {
        double I=s->iq_i[i], Q=s->iq_q[i];
        s->iq_i[i]=I*cr-Q*sr; s->iq_q[i]=I*sr+Q*cr;
        s->iq_raw[2*i]=(uint8_t)(s->iq_i[i]*128.0+127.5);
        s->iq_raw[2*i+1]=(uint8_t)(s->iq_q[i]*128.0+127.5);
    }
    fprintf(stderr,"  [Phase] ch=%d %.2f rad\n",ch,rad);
}
static void gate_CZ(SdrDev *s) {
    sdr_read(s); int n1=s->iq_n; if (n1<g_dim) return;
    double *sv=malloc(n1*2*8);
    memcpy(sv,s->iq_i,n1*8); memcpy(sv+n1,s->iq_q,n1*8);
    sdr_read(s); int n=n1<s->iq_n?n1:s->iq_n;
    for (int i=0;i<n;i++) {
        double ar=sv[i],ai=sv[i+n1],br=s->iq_i[i],bi=s->iq_q[i];
        s->iq_i[i]=ar*br-ai*bi; s->iq_q[i]=ar*bi+ai*br;
        s->iq_raw[2*i]=(uint8_t)(s->iq_i[i]*128.0+127.5);
        s->iq_raw[2*i+1]=(uint8_t)(s->iq_q[i]*128.0+127.5);
    }
    free(sv); fprintf(stderr,"  [CZ] %d pairs\n",n);
}
static void gate_X(SdrDev *s, int shift) {
    sdr_read(s); if (s->iq_n<g_dim) return;
    int blk=s->iq_n/g_dim;
    for (int b=0;b<blk;b++) {
        double *ti=malloc(g_dim*8), *tq=malloc(g_dim*8);
        int base=b*g_dim;
        for (int k=0;k<g_dim;k++){ti[k]=s->iq_i[base+k];tq[k]=s->iq_q[base+k];}
        for (int k=0;k<g_dim;k++) {
            int src=(k-shift+g_dim)%g_dim;
            s->iq_i[base+k]=ti[src]; s->iq_q[base+k]=tq[src];
            s->iq_raw[2*(base+k)]=(uint8_t)(ti[src]*128.0+127.5);
            s->iq_raw[2*(base+k)+1]=(uint8_t)(tq[src]*128.0+127.5);
        }
        free(ti); free(tq);
    }
    fprintf(stderr,"  [X] shift %+d\n",shift);
}

/*
 * GATE:TX — Synthesize CU8 I/Q for external SDR transmitter
 * Uses IDFT of wavefunction to generate multi-tone OFDM baseband.
 * Output to file or stdout for use with HackRF, LimeSDR, etc.
 */
static void gate_TX(Wavefunction *wf, uint32_t freq, uint32_t rate,
                    const char *outfile) {
    int np = IQ_WINDOW / 2;
    double *tx_i = calloc(np, sizeof(double));
    double *tx_q = calloc(np, sizeof(double));
    uint8_t *cu8 = malloc(np * 2);

    for (int k = 0; k < g_dim; k++) {
        double re = wf->re[k], im = wf->im[k];
        if (re*re + im*im < 1e-30) continue;
        double fn = (double)k / (double)g_dim;
        for (int n = 0; n < np; n++) {
            double ph = 2.0 * M_PI * fn * (double)n;
            double cr = cos(ph), sr = sin(ph);
            tx_i[n] += re * cr - im * sr;
            tx_q[n] += re * sr + im * cr;
        }
    }
    double peak = 0;
    for (int n = 0; n < np; n++) {
        double m = fabs(tx_i[n]); if (fabs(tx_q[n]) > m) m = fabs(tx_q[n]);
        if (m > peak) peak = m;
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
    if (fp) { fwrite(cu8, 1, np * 2, fp); if (fp != stdout) fclose(fp); }
    fprintf(stderr, "  [TX] %d I/Q pairs CU8 @ %.1f MHz → %s\n",
            np, freq/1e6,
            outfile ? (strcmp(outfile,"-")==0?"stdout":outfile) : "stdout");
    free(tx_i); free(tx_q); free(cu8);
}

static int gate_Measure(SdrDev *s, Wavefunction *wf) {
    uint64_t ent=0;
    for (int i=0;i<64&&i<s->iq_n*2;i++) ent=(ent<<1)|(s->iq_raw[i]&1);
    ent^=ent>>33; ent*=0xFF51AFD7ED558CCDULL; ent^=ent>>33;
    ent*=0xC4CEB9FE1A85EC53ULL; ent^=ent>>33;
    double r=(double)(ent>>11)/(double)(1ULL<<53);
    double cum=0; int o=g_dim-1;
    for (int k=0;k<g_dim;k++){cum+=wf->prob[k];if(r<=cum){o=k;break;}}
    fprintf(stderr,"  [MEASURE] → |%d⟩ r=%.4f\n",o,r);
    return o;
}

/* ═══════════════════════════════════════════
 * MAIN
 * ═══════════════════════════════════════════ */
int main(int argc, char **argv) {
    g_dim = 6;
    uint32_t freq = 100000000;
    uint32_t rate = 2048000;
    int    gain  = 400;
    int    tx_mode = 0;
    char  *tx_file = NULL;
    int    pos = 0;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--tx-file") == 0 && i+1 < argc) {
            tx_mode = 1; tx_file = argv[++i];
        } else if (strcmp(argv[i], "--tx-stdout") == 0) {
            tx_mode = 1; tx_file = "-";
        } else if (argv[i][0] != '-') {
            pos++;
            if (pos == 1) g_dim = atoi(argv[i]);
            if (pos == 2) freq = (uint32_t)atof(argv[i]);
            if (pos == 3) rate = (uint32_t)atof(argv[i]);
            if (pos == 4) gain = atoi(argv[i]);
        }
    }
    g_dim = g_dim < 2 ? 2 : (g_dim > MAX_DIM ? MAX_DIM : g_dim);

    if (tx_mode) {
        Wavefunction wf = wf_alloc(g_dim);
        for (int k=0;k<g_dim;k++) wf.re[k]=1.0/sqrt(g_dim);
        for (int k=0;k<g_dim;k++) wf.prob[k]=1.0/g_dim;
        wf.purity=1.0/g_dim; wf.entropy=log2(g_dim);
        gate_TX(&wf, freq, rate, tx_file);
        wf_free(&wf);
        return 0;
    }

    SdrDev sdr;
    int sdr_ok = (sdr_open(&sdr, freq, rate, gain) == 0);

    printf("\n");
    printf("  ╔══════════════════════════════════════════════════════════════╗\n");
    printf("  ║  RTL-SDR QUANTUM ENGINE — D=%-3d                              ║\n", g_dim);
    printf("  ║  %.1f MHz | %.3f MSPS | gain=%.1f dB                          ║\n",
           freq/1e6, rate/1e6, gain/10.0);
    printf("  ║  Source: %-50s ║\n", sdr_ok ? "PHYSICAL (RTL-SDR)" : "SIMULATED");
    printf("  ╚══════════════════════════════════════════════════════════════╝\n\n");

    Wavefunction wf = wf_alloc(g_dim);

    if (sdr_ok) { wf_compute(&sdr, &wf); }
    else { for (int k=0;k<g_dim;k++) wf.re[k]=1.0/sqrt(g_dim);
           for (int k=0;k<g_dim;k++) wf.prob[k]=1.0/g_dim;
           wf.purity=1.0/g_dim; wf.entropy=log2(g_dim); }
    wf_print(&wf, "|ψ₀⟩");

    if (sdr_ok) { gate_H(&sdr); wf_compute(&sdr, &wf); } wf_print(&wf, "H");
    if (sdr_ok) { gate_DFT(&sdr,1); wf_compute(&sdr, &wf); } wf_print(&wf, "DFT");
    if (sdr_ok) gate_DFT(&sdr,-1);
    if (sdr_ok) { gate_Phase(&sdr,2%g_dim,M_PI/4); wf_compute(&sdr, &wf); } wf_print(&wf, "Phase");
    if (sdr_ok) { gate_CZ(&sdr); wf_compute(&sdr, &wf); } wf_print(&wf, "CZ");
    if (sdr_ok) { gate_X(&sdr,1); wf_compute(&sdr, &wf); } wf_print(&wf, "X");

    int outcome;
    if (sdr_ok) outcome = gate_Measure(&sdr, &wf);
    else outcome = rand() % g_dim;

    printf("\n  ★ COLLAPSED → |%d⟩  |  samples: %lu ★\n\n",
           outcome, (unsigned long)sdr.samples);

    wf_free(&wf);
    if (sdr_ok) sdr_close(&sdr);
    return 0;
}
