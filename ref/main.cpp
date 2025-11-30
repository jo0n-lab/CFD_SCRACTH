#include <iostream>
#define COMMON_DEF
#include "common.h" 


int main()
{
    NI = 22; DX0 = 0.05; NJ = 22; DY0 = 0.05; U0 = 0.; V0 = 0.; P0 = 0.; T0 = 0.;
    IR = 0; JUMP = 50, GRAVY = 0; BETA = 0.e-3;
    K0 = 1.; ROCP0 = 1.e+6; QDOT = 0.; RHOP0 = 1.e3; VIS0 = 1.e-3;
    ITMAX = 1000; ERRMAX = 1.e-5; RELAX = 1.;  RELAXUV = 0.5; RELAXP = 1.; //RELAXT = 1.;

    //IBCW = IWAL; QBCW = 0.;
    //IBCE = IWQQ; QBCE = 0.;
    //IBCS = IWQQ; TBCS = 0.;
    //IBCN = IWAL; TBCN = 1.;

    //IBCW = IOUT; TBCW = 0.; QBCW = 0.; UBCW = 0.; VBCW = 0.; PBCW = 0.;
    //IBCE = IOUT; TBCE = 0.; QBCE = 0.; UBCE = 0.; VBCE = 0.; PBCE = 0.;
    //IBCS = IINF; TBCS = 0.; QBCS = 0.; UBCS = 0.; VBCS = 0.; PBCS = 0.;
    //IBCN = IINF; TBCN = 0.; QBCN = 0.; UBCN = 1.; VBCN = 0.; PBCN = 0.;
     //4-1-1


    //IBCW = IINF; TBCE = 0.; QBCE = 0.; UBCE = 1.; VBCE = 0.; PBCE = 0.;
    //IBCS = ISYM; TBCW = 0.; QBCW = 0.; UBCW = 0.; VBCW = 0.; PBCW = 0.;
    //IBCE = IINF; TBCS = 0.; QBCS = 0.; UBCS = 0.; VBCS = 0.; PBCS = 0.;
    //IBCN = ISYM; TBCN = 0.; QBCN = 0.; UBCN = 0.; VBCN = 0.; PBCN = 0.;
    //4-1-2


    //IBCW = ISYM; TBCW = 0.; QBCW = 0.; UBCW = 0.; VBCW = 0.; PBCW = 0.;
    //IBCE = ISYM; TBCE = 0.; QBCE = 0.; UBCE = 0.; VBCE = 0.; PBCE = 0.;
    //IBCS = IINF; TBCS = 0.; QBCS = 0.; UBCS = 0.; VBCS = 0.; PBCS = 0.;
    //IBCN = IINF; TBCN = 0.; QBCN = 0.; UBCN = 0.; VBCN = 1.; PBCN = 0.;
    //4-2


    //IBCW = IWAL; TBCW = 0.; QBCW = 0.; UBCW = 0.; VBCW = 0.; PBCW = 0.;
    //IBCE = IWAL; TBCE = 0.; QBCE = 0.; UBCE = 0.; VBCE = 0.1; PBCE = 0.;
    //IBCS = IINF; TBCS = 0.; QBCS = 0.; UBCS = 0.; VBCS = 0.; PBCS = 0.;
    //IBCN = IOUT; TBCN = 0.; QBCN = 0.; UBCN = 0.; VBCN = 0.; PBCN = 0.;
    //4-3

    IBCW = IWAL; TBCW = 0.; QBCW = 0.; UBCW = 0.; VBCW = 0.; PBCW = 0.;
    IBCE = IWAL; TBCE = 0.; QBCE = 0.; UBCE = 0.; VBCE = 0.1; PBCE = 0.;
    IBCS = IINF; TBCS = 0.; QBCS = 0.; UBCS = 0.; VBCS = 0.; PBCS = 0.;
    IBCN = IPRE; TBCN = 0.; QBCN = 0.; UBCN = 0.; VBCN = 0.; PBCN = 0.;
    //4-4

    INIT(); //PROP();

    for (i = 1; i < NI; i++) { V[i][1] = 0.1 * XP[i]; } // V에 대한 BC인듯 v = 0.1x case:4-3 & 4-4


    printf("%4d %11.3e %11.3e %11.3e %11.3e %11.3e %11.3e %11.3e \n",
       0, U[10][10], V[10][10], P[10][10], 0., 0., 0., 0.);


    for (int iter = 0; iter <= ITMAX; iter++) {

        VELMAX = 1.e-20;
        for (j = 1; j <= NJ; j++) {
            for (i = 1; i <= NI; i++) {
                BU[i][j] = U[i][j]; BV[i][j] = V[i][j]; BT[i][j] = T[i][j];
                VELMAX = fmax(VELMAX, U[i][j]); VELMAX = fmax(VELMAX, V[i][j]);

                //UU[i][j] = 0.; VV[i][j] = 0.1 * XP[i]; V[i][j] = 0.1 * XP[i];
            }
        }


        PXY_GET(); U_SOLVE(); V_SOLVE(); P_SOLVE();

        double uvsmax = fmax(USMAX, VSMAX); double uvcmax = fmax(UCMAX, VCMAX);
        SSMAX = fmax(uvsmax, PSMAX); FCMAX = fmax(uvcmax, PCMAX);
        //SSMAX = fmax(TSMAX, SSMAX); FCMAX = fmax(TCMAX, FCMAX);

        if (iter == iter / JUMP * JUMP) {
            printf("%5d %11.3e %11.3e %11.3e %11.3e %11.3e %11.3e %11.3e \n",
                iter, U[10][10], V[10][10], P[10][10], uvsmax, PSMAX, uvcmax, PCMAX);
        }

        if (iter > 0 && SSMAX < ERRMAX && FCMAX < ERRMAX) {
            printf("%5d %11.3e %11.3e %11.3e %11.3e %11.3e %11.3e %11.3e \n",
                iter, U[10][10], V[10][10], P[10][10], uvsmax, PSMAX, uvcmax, PCMAX);
            break;
        }
    }

    FILE* fp = fopen("out_x", "w");

    double Nuav = 0.;
    for (i = 1; i <= NI; i++) {
        //double Pe = RHOP0 * U0 / VIS0; // / VIS0 * fmax(fabs(U0), fabs(V0));
        //double Uex = (exp(Pe * XP[i]) - 1.) / (exp(Pe) - 1.);
        double Vex = 0.1 * XP[i];
        fprintf(fp, "%4d %15.7e %15.7e %15.7e %15.7e %15.7e\n", i, XP[i], U[i][7], V[i][7], P[i][7], Vex);
    }
    fclose(fp);

    FILE* fp1 = fopen("out_y", "w");
    for (j = 1; j <= NJ; j++) {
       double Pe = RHOP0 * fmax(U0, V0) / VIS0; // / VIS0 * fmax(fabs(U0), fabs(V0));
       double Uex = (exp(Pe * YP[j]) - 1.) / (exp(Pe) - 1.);

        fprintf(fp1, "%4d %15.7e %15.7e %15.7e %15.7e\n", j, YP[j], U[7][j], P[7][j], Uex);
    }
    fclose(fp1);

    FILE* fp2 = fopen("out_p", "w");
    for (j = 1; j <= NJ; j++) {
       for (i = 1; i <= NI; i++) {
            fprintf(fp2,"%15.7e ", P[i][j]);
       }
       fprintf(fp2,"\n");
    }
    fclose(fp2);

    

    return 0;
}

