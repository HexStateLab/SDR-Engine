/*
 * qvm_api.h — Public API for the Ether Quantum VM
 *
 * Include this header to drive the room as a computational substrate
 * from external programs (CUDA wrappers, Python bindings, test suites).
 *
 * Usage:
 *   QvmCtx *q = qvm_create(100e6, 2.048e6, 256, 496);
 *   qvm_eval(q, "INIT");
 *   qvm_eval(q, "CALIBRATE 4");
 *   qvm_eval(q, "MEASURE");
 *   qvm_destroy(q);
 */

#ifndef QVM_API_H
#define QVM_API_H

#include <stdint.h>
#include <stdio.h>

/* ─── Opaque context (caller treats as handle) ─── */
typedef struct QvmCtx QvmCtx;

/* ─── Lifecycle ─── */
QvmCtx *qvm_create(uint32_t freq_hz, uint32_t rate_hz, int dim, int gain);
void    qvm_destroy(QvmCtx *q);

/* ─── Execution ─── */
int  qvm_eval(QvmCtx *q, const char *instruction);
int  qvm_run(QvmCtx *q, const char *script_path);
int  qvm_running(QvmCtx *q);

/* ─── State access ─── */
int    qvm_dim(QvmCtx *q);
double qvm_prob(QvmCtx *q, int level);
void   qvm_probs(QvmCtx *q, double *out, int max);
double qvm_entropy(QvmCtx *q);
double qvm_purity(QvmCtx *q);

/* ─── Hardware info ─── */
int qvm_has_sdr(QvmCtx *q);

/* ─── Extended: GPU offload ─── */
int    qvm_calibrated(QvmCtx *q);
void   qvm_get_channel(QvmCtx *q, double **M_out, double **Minv_out, int *dim_out);
int    qvm_compute(QvmCtx *q, const double *x, double *y, int d);

/* ─── Low-level SDR access (for advanced waveform synthesis) ─── */
void   qvm_sdr_tune(QvmCtx *q, uint32_t freq_hz);
void   qvm_sdr_flush(QvmCtx *q);
int    qvm_sdr_rx(QvmCtx *q, double *I_out, double *Q_out, int max_samples);
void   qvm_sdr_dft(QvmCtx *q, double *pwr, int D, const double *I, const double *Q, int n_samples);

#endif /* QVM_API_H */
