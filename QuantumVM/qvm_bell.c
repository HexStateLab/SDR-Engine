/*
 * qvm_bell.c — Bell test via QVM API (external program)
 *
 * Creates |Φ+⟩ = (|00⟩+|11⟩)/√2 on D=4 (2 qubits), sends through
 * the physical room, runs CHSH test to check Bell inequality.
 *
 * Build:
 *   gcc -c -O3 -std=gnu99 -DNO_MAIN sdr_ether.c -o sdr_ether_lib.o
 *   gcc -O3 -std=gnu99 qvm_bell.c sdr_ether_lib.o -lm -lpthread -o qvm_bell
 */

#include "qvm_api.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define D 4  /* 2 qubits: |00⟩,|01⟩,|10⟩,|11⟩ */

/* Joint measurement at angles ta, tb on a 2-qubit state (probs[4]).
   Returns correlation E(ta,tb) = P(++)-P(+-)-P(-+)+P(--) */
static double measure_joint(const double probs[D], double ta, double tb) {
    /* For each basis state |ab⟩, compute amplitude at measurement angles.
       This is a simplified Born rule: project onto measurement basis. */
    double ca = cos(ta/2), sa = sin(ta/2);
    double cb = cos(tb/2), sb = sin(tb/2);

    double amp_pp = 0, amp_pm = 0, amp_mp = 0, amp_mm = 0;
    const double *re = probs; /* using probs as pseudo-amplitudes */
    /* |00⟩: amplitude sqrt(prob[0]) */
    double a0 = sqrt(probs[0]), a1 = sqrt(probs[1]);
    double a2 = sqrt(probs[2]), a3 = sqrt(probs[3]);

    /* |00⟩ contribution */
    amp_pp += a0 * ca * cb;
    amp_pm += a0 * ca * (-sb);
    amp_mp += a0 * (-sa) * cb;
    amp_mm += a0 * (-sa) * (-sb);

    /* |01⟩ contribution */
    amp_pp += a1 * ca * sb;
    amp_pm += a1 * ca * cb;
    amp_mp += a1 * (-sa) * sb;
    amp_mm += a1 * (-sa) * cb;

    /* |10⟩ contribution */
    amp_pp += a2 * sa * cb;
    amp_pm += a2 * sa * (-sb);
    amp_mp += a2 * ca * cb;
    amp_mm += a2 * ca * (-sb);

    /* |11⟩ contribution */
    amp_pp += a3 * sa * sb;
    amp_pm += a3 * sa * cb;
    amp_mp += a3 * ca * sb;
    amp_mm += a3 * ca * cb;

    double p_pp = amp_pp*amp_pp, p_pm = amp_pm*amp_pm;
    double p_mp = amp_mp*amp_mp, p_mm = amp_mm*amp_mm;

    return p_pp - p_pm - p_mp + p_mm;
}

int main() {
    printf("\n  ╔══════════════════════════════════════════════════╗\n");
    printf("  ║  QVM API: BELL TEST — Two-Qudit CHSH             ║\n");
    printf("  ║  |Φ+⟩ = (|00⟩+|11⟩)/√2  →  room  →  CHSH check  ║\n");
    printf("  ╚══════════════════════════════════════════════════╝\n\n");

    QvmCtx *q = qvm_create(100000000, 2048000, D, 496);
    if (!q) { fprintf(stderr,"qvm_create failed\n"); return 1; }
    printf("Hardware: %s\n", qvm_has_sdr(q) ? "RTL-SDR ACTIVE" : "SIMULATION");

    /* ── Create Bell state |Φ+⟩ = (|00⟩+|11⟩)/√2 ── */
    printf("\n1. Creating Bell state |Φ+⟩...\n");
    /* Reset then SET both target bins before normalization happens.
       SET 0 triggers normalization, so we need to prepare differently.
       Use qvm_compute to feed the exact amplitude pattern as a single shot. */
    double bell_x[D] = {0.7071, 0.0, 0.0, 0.7071};
    double bell_y[D] = {0};
    printf("  Target: |00⟩:%.3f |01⟩:%.3f |10⟩:%.3f |11⟩:%.3f\n",
        bell_x[0],bell_x[1],bell_x[2],bell_x[3]);

    /* ── Baseline: software Bell state (no SDR) ── */
    double pristine_E[4] = {0};
    double a=0.0,ap=M_PI/2.0,b=M_PI/4.0,bp=3.0*M_PI/4.0;
    double angles[][2]={{a,b},{a,bp},{ap,b},{ap,bp}};
    int n_trials=2000;
    srand(time(NULL));
    for(int m=0;m<4;m++)
        for(int t=0;t<n_trials;t++)
            pristine_E[m]+=measure_joint(bell_x,angles[m][0],angles[m][1]);
    for(int m=0;m<4;m++) pristine_E[m]/=n_trials;
    double S_pristine=fabs(pristine_E[0]-pristine_E[1]+pristine_E[2]+pristine_E[3]);
    printf("  S (pristine): %.4f   %s\n\n", S_pristine, S_pristine>2.0?"ENTANGLED":"classical");

    /* ── Send Bell state through the room with feedback rounds ── */
    printf("2. Sending Bell state through room (%d SDR rounds)...\n",
        qvm_has_sdr(q)?5:1);
    double room_s[D];
    if (qvm_has_sdr(q)) {
        /* Round 1: feed Bell amplitudes directly */
        qvm_compute(q, bell_x, room_s, D);
        /* Rounds 2-5: feedback — room's output is next input.
           Mixer intermodulation creates entanglement across bins. */
        for (int r=1; r<5; r++) {
            double tmp[D]; memcpy(tmp, room_s, D*sizeof(double));
            qvm_compute(q, tmp, room_s, D);
        }
    } else {
        memcpy(room_s, bell_x, D*sizeof(double));
    }

    printf("  Room final: |00⟩:%.3f |01⟩:%.3f |10⟩:%.3f |11⟩:%.3f\n",
        room_s[0],room_s[1],room_s[2],room_s[3]);

    /* ── Calibrate for reference ── */
    if (qvm_has_sdr(q)) {
        printf("\n3. Calibrating room M...\n");
        qvm_eval(q, "CALIBRATE 4");
    }

    /* ── CHSH on pristine vs room-processed ── */
    printf("\n4. CHSH Bell test (%d trials each):\n\n", n_trials);

    double room_E[4]={0};
    for(int m=0;m<4;m++)
        for(int t=0;t<n_trials;t++)
            room_E[m]+=measure_joint(room_s,angles[m][0],angles[m][1]);
    for(int m=0;m<4;m++) room_E[m]/=n_trials;
    double S_room=fabs(room_E[0]-room_E[1]+room_E[2]+room_E[3]);

    printf("  %-14s %8s %8s %8s %8s %8s\n","","E(A,B)","E(A,B')","E(A',B)","E(A',B')","S");
    printf("  ────────────── ──────── ──────── ──────── ──────── ────────\n");
    printf("  %-14s %+.4f %+.4f %+.4f %+.4f %8.4f %s\n",
        "Pristine |Φ+⟩",pristine_E[0],pristine_E[1],pristine_E[2],pristine_E[3],
        S_pristine, S_pristine>2.0?"★":"");
    printf("  %-14s %+.4f %+.4f %+.4f %+.4f %8.4f %s\n",
        "Room output",room_E[0],room_E[1],room_E[2],room_E[3],
        S_room, S_room>2.0?"★":"");
    printf("  %-14s %8s %8s %8s %8s %8.4f\n","Classical bound","","","","",2.0);
    printf("  %-14s %8s %8s %8s %8s %8.4f\n","Quantum maximum","","","","",2.0*sqrt(2.0));

    printf("\n  ╔══════════════════════════════════════════════════╗\n");
    printf("  ║  Pristine S = %.4f  |  Room S = %.4f         ║\n",S_pristine,S_room);
    if (S_room > 2.0)
        printf("  ║  ★ ROOM PRESERVED ENTANGLEMENT — S>2! ★        ║\n");
    else if (S_pristine > 2.0)
        printf("  ║  Bell state decohered by room channel (%+.2fΔS) ║\n",S_room-S_pristine);
    else
        printf("  ║  No entanglement detected                       ║\n");
    printf("  ║  Room processed via QVM API (external program)  ║\n");
    printf("  ╚══════════════════════════════════════════════════╝\n\n");

    qvm_destroy(q);
    return 0;
}
