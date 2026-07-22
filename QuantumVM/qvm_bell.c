/*
 * qvm_bell.c — Calibrated + M⁺-precompensated Bell test
 * Build: gcc -O3 qvm_bell.c sdr_ether_lib.o -lm -lpthread -o qvm_bf
 */
#include "qvm_api.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define D 4

static void matvec(double **A, double *x, double *y, int n){
    for(int i=0;i<n;i++){y[i]=0;for(int j=0;j<n;j++)y[i]+=A[i][j]*x[j];}
}

static double chsh(const double *p, double a, double b){
    double ca=cos(a/2),sa=sin(a/2),cb=cos(b/2),sb=sin(b/2);
    double amp[4];for(int i=0;i<4;i++)amp[i]=sqrt(p[i]);
    double pp=amp[0]*ca*cb+amp[1]*ca*sb+amp[2]*sa*cb+amp[3]*sa*sb;
    double pm=amp[0]*ca*(-sb)+amp[1]*ca*cb+amp[2]*sa*(-sb)+amp[3]*sa*cb;
    double mp=amp[0]*(-sa)*cb+amp[1]*(-sa)*sb+amp[2]*ca*cb+amp[3]*ca*sb;
    double mm=amp[0]*(-sa)*(-sb)+amp[1]*(-sa)*cb+amp[2]*ca*(-sb)+amp[3]*ca*cb;
    return pp*pp-pm*pm-mp*mp+mm*mm;
}

int main(){
    printf("\n╔══════════════════════════════════════════════╗\n");
    printf("║  CALIBRATED BELL — M⁺ pre-comp each round    ║\n");
    printf("╚══════════════════════════════════════════════╝\n\n");

    QvmCtx *q=qvm_create(100000000,2048000,D,496);
    if(!q||!qvm_has_sdr(q)){printf("No SDR\n");qvm_destroy(q);return 1;}

    /* Calibrate */
    printf("Calibrating M...\n");
    qvm_eval(q,"CALIBRATE 12");
    double *Mraw,*Miraw;int Md;
    qvm_get_channel(q,&Mraw,&Miraw,&Md);
    double **M=(double**)Mraw;

    /* Compute M⁺ via Tikhonov */
    double MTM[4][4]={{0}},MTMi[4][4]={{0}},Mpinv[4][4]={{0}};
    for(int i=0;i<4;i++)for(int j=0;j<4;j++)for(int k=0;k<4;k++)MTM[i][j]+=M[k][i]*M[k][j];
    for(int i=0;i<4;i++)MTM[i][i]+=0.001; /* sharper inverse */
    double aug[4][8];memset(aug,0,sizeof(aug));
    for(int i=0;i<4;i++){for(int j=0;j<4;j++)aug[i][j]=MTM[i][j];aug[i][4+i]=1;}
    for(int c=0;c<4;c++){
        int piv=c;double mv=fabs(aug[c][c]);
        for(int r=c+1;r<4;r++){double v=fabs(aug[r][c]);if(v>mv){mv=v;piv=r;}}
        if(piv!=c)for(int j=0;j<8;j++){double t=aug[c][j];aug[c][j]=aug[piv][j];aug[piv][j]=t;}
        double pv=aug[c][c];for(int j=0;j<8;j++)aug[c][j]/=pv;
        for(int r=0;r<4;r++){if(r==c)continue;double f=aug[r][c];for(int j=0;j<8;j++)aug[r][j]-=f*aug[c][j];}
    }
    for(int i=0;i<4;i++)for(int j=0;j<4;j++)MTMi[i][j]=aug[i][4+j];
    double **Mp=(double**)malloc(4*sizeof(double*));
    for(int i=0;i<4;i++){Mp[i]=calloc(4,8);
        for(int j=0;j<4;j++)for(int k=0;k<4;k++)Mp[i][j]+=MTMi[i][k]*M[j][k];}

    printf("M⁺ computed. Starting Bell test...\n\n");

    /* Bell state */
    double bell[4]={0.5,0,0,0.5};
    printf("Target: |00⟩:%.3f |11⟩:%.3f\n\n",bell[0],bell[3]);

    /* 5 rounds: pre-comp xp = M⁺·target, room computes M·xp */
    double state[4];memcpy(state,bell,32);
    int rounds=8;
    for(int r=0;r<rounds;r++){
        double xp[4]={0};matvec(Mp,state,xp,4);
        for(int i=0;i<4;i++)if(xp[i]<0)xp[i]=0;
        double t=0;for(int i=0;i<4;i++)t+=xp[i];
        if(t>1e-15)for(int i=0;i<4;i++)xp[i]/=t;

        double y[4];
        qvm_compute(q,xp,y,D);

        /* Normalize room output */
        double s=0;for(int i=0;i<4;i++)s+=y[i];
        if(s>1e-15)for(int i=0;i<4;i++)y[i]/=s;
        memcpy(state,y,32);

        printf("  R%d: |00⟩:%.4f |01⟩:%.4f |10⟩:%.4f |11⟩:%.4f\n",
            r+1,state[0],state[1],state[2],state[3]);
    }

    /* H gate to redistribute before CHSH */
    qvm_eval(q,"H");
    qvm_probs(q,state,D);
    printf("  H:   |00⟩:%.4f |01⟩:%.4f |10⟩:%.4f |11⟩:%.4f\n",
        state[0],state[1],state[2],state[3]);

    /* CHSH */
    double a=0,ap=M_PI/2,b=M_PI/4,bp=3*M_PI/4;
    double E[4]={chsh(state,a,b),chsh(state,a,bp),chsh(state,ap,b),chsh(state,ap,bp)};
    double S=fabs(E[0]-E[1]+E[2]+E[3]);
    double Eb[4]={chsh(bell,a,b),chsh(bell,a,bp),chsh(bell,ap,b),chsh(bell,ap,bp)};
    double Sb=fabs(Eb[0]-Eb[1]+Eb[2]+Eb[3]);

    printf("\nCHSH: Room S=%.4f  Bell S=%.4f  Classical=2.0\n",S,Sb);
    printf("%s\n\n",S>2.0?"★★ VIOLATED ★★":"decohered");

    for(int i=0;i<4;i++)free(Mp[i]);free(Mp);
    qvm_destroy(q);
    return 0;
}
