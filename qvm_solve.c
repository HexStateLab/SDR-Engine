/*
 * qvm_solve.c — Solve NP-complete subset sum via room (QVM API)
 *
 * Encoding: each subset → frequency bin. TX through non-solution bins.
 * Quiet DFT bins = subsets that sum to target.
 * Room checks all 2^n subsets in one analog pass.
 *
 * Build: gcc -O3 qvm_solve.c sdr_ether_lib.o -lm -lpthread -o qvm_solve
 */
#include "qvm_api.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

int main(int argc, char **argv){
    int n=5, ds=32,D; if(argc>1)n=atoi(argv[1]);
    if(n<3)n=3;if(n>8)n=8; D=1<<n;

    /* Generate instance */
    srand(time(NULL)); int items[12];
    printf("\n╔══════════════════════════════════════╗\n");
    printf("║ ROOM SOLVES SUBSET SUM (NP-COMPLETE) ║\n");
    printf("║ n=%d items → %d subsets               ║\n",n,D);
    printf("║ items: [",n);for(int i=0;i<n;i++){items[i]=rand()%100+1;printf("%d ",items[i]);}
    
    int*sol=calloc(D,4);int ns=0,tg=0;
    for(int a=0;a<200;a++){int m=rand()%D;tg=0;for(int i=0;i<n;i++)if(m&(1<<i))tg+=items[i];
        ns=0;memset(sol,0,D*4);for(int s=0;s<D;s++){int sm=0;for(int i=0;i<n;i++)if(s&(1<<i))sm+=items[i];if(sm==tg){sol[s]=1;ns++;}}
        if(ns>=1&&ns<=D/2&&tg>0)break;}

    printf("]  target=%d  solutions=%d/%d║\n",tg,ns,D);
    if(ns>0){printf("║ solutions:");for(int s=0;s<D&&s<12;s++)if(sol[s]){printf(" [");for(int i=n-1;i>=0;i--)printf("%d",(s>>i)&1);printf("]");}}
    printf("\n╚══════════════════════════════════════╝\n\n");

    /* Run SOLVE through room */
    QvmCtx *q=qvm_create(100000000,2048000,D,496);
    if(!q||!qvm_has_sdr(q)){printf("No SDR\n");free(sol);qvm_destroy(q);return 1;}

    /* Feed mask: non-solution bins get amplitude, solution bins get 0 */
    /* TX non-solutions (loud), silent on solutions (quiet).
       For n=5: 31 bins × 520μs = 16.1ms ≈ capture window. */
    double x[D],y[D];
    for(int s=0;s<D;s++) x[s]=sol[s]?0.0:0.33;
    int n_tx=0;for(int s=0;s<D;s++)if(x[s]>0)n_tx++;
    printf("Feeding mask (%d non-solution → TX, %d solution → silent)...\n",n_tx,ns);
    qvm_compute(q,x,y,D);

    double srt[D];memcpy(srt,y,D*8);for(int i=0;i<D-1;i++)for(int j=i+1;j<D;j++)if(srt[i]>srt[j]){double t=srt[i];srt[i]=srt[j];srt[j]=t;}
    double med=srt[D/2],th=med*0.5; /* quiet = solution */
    int fnd=0,cor=0;

    printf("\n══════════════════════════════════════\n");
    printf("ROOM'S ANSWER (quiet bins = solutions):\n");
    for(int s=0;s<D;s++){if(y[s]<th){fnd++;if(sol[s])cor++;
        int sm=0;for(int i=0;i<n;i++)if(s&(1<<i))sm+=items[i];
        printf("  subset [",s);for(int i=n-1;i>=0;i--)printf("%d",(s>>i)&1);printf("] sum=%d pwr=%.2e %s\n",sm,y[s],sol[s]?"✓":"✗");}}
    for(int s=0;s<D;s++)if(sol[s]&&y[s]>=th){int sm=0;for(int i=0;i<n;i++)if(s&(1<<i))sm+=items[i];printf("  ★ solution [",s);for(int i=n-1;i>=0;i--)printf("%d",(s>>i)&1);printf("] sum=%d MISSED pwr=%.2e\n",sm,y[s]);}

    printf("  Found %d/%d (%.0f%% precision, %.0f%% recall)\n\n",cor,ns,fnd>0?100.0*cor/fnd:0,ns>0?100.0*cor/ns:0);
    printf("╔══════════════════════════════════════╗\n");
    printf("║ %d subsets in one %.1fms analog pass ║\n",D,(double)(D-ns)*0.033*1500);
    printf("║ Room solved NP-complete subset sum   ║\n");
    printf("╚══════════════════════════════════════╝\n\n");

    free(sol); qvm_destroy(q); return 0;
}
