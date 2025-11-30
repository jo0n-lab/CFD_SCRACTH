#include <iostream>
#include "common.h"

double TKK_INTP(double a, double b) { return 0.5 * (a + b); }
double VIS_INTP(double a, double b) { return 0.5 * (a + b); }

void INIT() {
    NIM = NI - 1; NJM = NJ - 1;

    for (i = 2; i <= NI; i++) XU[i] = DX0 * (i - 2.);
    for (i = 2; i <= NIM; i++) XP[i] = 0.5 * (XU[i + 1] + XU[i]); XP[1] = XU[2]; XP[NI] = XU[NI];
    for (i = 2; i <= NIM; i++) DXP[i] = XU[i + 1] - XU[i];
    for (i = 2; i <= NI; i++) DXU[i] = XP[i] - XP[i - 1];

    for (i = 2; i <= NIM; i++) DXUE[i] = 0.5 * DXP[i] / DXU[i]; DXUE[NI] = 0.;
    for (i = 2; i <= NI; i++) DXUW[i] = 1. - DXUE[i];

    for (j = 2; j <= NJ; j++) YV[j] = DY0 * (j - 2.);
    for (j = 2; j <= NJM; j++) YP[j] = 0.5 * (YV[j + 1] + YV[j]); YP[1] = YV[2]; YP[NJ] = YV[NJ];
    for (j = 2; j <= NJM; j++) DYP[j] = YV[j + 1] - YV[j];
    for (j = 2; j <= NJ; j++) DYV[j] = YP[j] - YP[j - 1];

    for (j = 2; j <= NJM; j++) DYVN[j] = 0.5 * DYP[j] / DYV[j]; DYVN[NJ] = 0.;
    for (j = 2; j <= NJ; j++) DYVS[j] = 1. - DYVN[j];

    for (i = 2; i <= NI; i++) RU[i] = 1.;
    for (i = 1; i <= NI; i++) RP[i] = 1.;

    if (IR == 1) {
        for (i = 2; i <= NI; i++) RU[i] = fabs(XU[i]);
        for (i = 1; i <= NI; i++) RP[i] = fabs(XP[i]);
    }

    for (i = 2; i <= NIM; i++) {
        for (j = 2; j <= NJM; j++) {
            VOL[i][j] = RP[i] * DXP[i] * DYP[j];
        }
    }

    for (i = 1; i <= NI; i++) {
        for (j = 1; j <= NJ; j++) {
            TKKP[i][j] = K0; ROCP[i][j] = ROCP0; VISP[i][j] = VIS0; RHOP[i][j] = RHOP0;
        }
    }

    for (i = 1; i <= NI; i++) { for (j = 1; j <= NJ; j++) { T[i][j] = T0; P[i][j] = P0; } }
    for (i = 1; i <= NI; i++) { for (j = 1; j <= NJ; j++) { U[i][j] = U0; V[i][j] = V0; } }
    //for (i = 1; i <= NI; i++) { for (j = 1; j <= NJ; j++) { U[i][j] = 0.; V[i][j] = 0.; } }

    for (j = 1; j <= NJ; j++) { T[1][j] = TBCW; T[NI][j] = TBCE; }
    for (i = 1; i <= NI; i++) { T[i][1] = TBCS; T[i][NJ] = TBCN; }
    for (j = 1; j <= NJ; j++) { U[1][j] = UBCW; U[NI][j] = UBCE; }
    for (i = 1; i <= NI; i++) { U[i][1] = UBCS; U[i][NJ] = UBCN; }
    for (j = 1; j <= NJ; j++) { V[1][j] = VBCW; V[NI][j] = VBCE; }
    for (i = 1; i <= NI; i++) { V[i][1] = VBCS; V[i][NJ] = VBCN; }
}

void PROP() {
    for (i = 1; i <= NI; i++) {
        for (j = 1; j <= NJ; j++) {
            TKKP[i][j] = K0; ROCP[i][j] = ROCP0; VISP[i][j] = VIS0; RHOP[i][j] = RHOP0;
        }
    }
}



void PXY_GET() {
    for (j = 2; j <= NJM; j++) {
        for (i = 2; i <= NI; i++) {
            PXU[i][j] = (P[i][j] - P[i - 1][j]) / DXU[i];
        }
        PXU[1][j] = PXU[2][j]; PXU[NI + 1][j] = PXU[NI][j];
    }

    for (i = 2; i <= NIM; i++) {
        for (j = 2; j <= NJ; j++) {
            PYV[i][j] = (P[i][j] - P[i][j - 1]) / DYV[j];
        }
        PYV[i][1] = PYV[i][2]; PYV[i][NJ + 1] = PYV[i][NJ];
    }
}

void U_SOLVE()
{
    PXY_GET();

    for (j = 2; j <= NJM; j++) {
        for (i = 2; i <= NIM; i++) {


            double visw = VIS_INTP(VISP[i][j], VISP[i - 1][j]);
            double vise = VIS_INTP(VISP[i][j], VISP[i + 1][j]);
            double viss = VIS_INTP(VISP[i][j], VISP[i][j - 1]);
            double visn = VIS_INTP(VISP[i][j], VISP[i][j + 1]);


            AW[i][j] = visw * RU[i] * DYP[j] / DXU[i];
            AE[i][j] = vise * RU[i + 1] * DYP[j] / DXU[i + 1];
            AS[i][j] = viss * RP[i] * DXP[i] / DYV[j];
            AN[i][j] = visn * RP[i] * DXP[i] / DYV[j + 1];


            double flw = RHOP[i][j] * DYP[j] * RU[i] * UU[i][j];
            double fle = RHOP[i][j] * DYP[j] * RU[i + 1] * UU[i + 1][j];
            double fls = RHOP[i][j] * DXP[i] * RP[i] * VV[i][j];
            double fln = RHOP[i][j] * DXP[i] * RP[i] * VV[i][j + 1];


            AW[i][j] += fmax(0., flw);
            AE[i][j] += fmax(0., -fle);
            AS[i][j] += fmax(0., fls);
            AN[i][j] += fmax(0., -fln);


            USOR[i][j] = 0.;


            double ac = 0.;
            AP[i][j] = AW[i][j] + AE[i][j] + AS[i][j] + AN[i][j] + ac;
            //USOR[i][j] += 0.;
            //if (IR == 1) AP[i][j] += VOL[i][j] * VISP[i][j] / RP[i] / RP[i];
            USOR[i][j] -= VOL[i][j] * 0.5 * (PXU[i][j] + PXU[i + 1][j]);

            AP[i][j] = AP[i][j] / RELAXUV;
            USOR[i][j] += (1. - RELAXUV) * AP[i][j] * U[i][j];
            APC[i][j] = AP[i][j];

            APU_P[i][j] = VOL[i][j] / AP[i][j];

        }
    }

    //------------------------------------ X-BC
    for (j = 2; j <= NJM; j++) {
        if (IBCW == IPRE || IBCW == IOUT) AP[2][j] -= AW[2][j];
        if (IBCE == IPRE || IBCE == IOUT) AP[NIM][j] -= AE[NIM][j];

        APU_P[1][j] = APU_P[2][j]; APU_P[NI][j] = APU_P[NIM][j];
        
    }

    //------------------------------------- Y-BC
    for (i = 2; i <= NIM; i++) {
        if (IBCS == ISYM) AP[i][2] -= AS[i][2];
        if (IBCN == ISYM) AP[i][NJM] -= AN[i][NJM];
        
    }


    for (j = 2; j <= NJM; j++) {
        for (i = 2; i <= NIM; i++) {
            SS[i][j] = USOR[i][j] - APC[i][j] * U[i][j]
                + AW[i][j] * U[i - 1][j] + AE[i][j] * U[i + 1][j]
                + AS[i][j] * U[i][j - 1] + AN[i][j] * U[i][j + 1];
        }
    }


    LNTDMA(2, NIM, 2, NJM);


    for (j = 2; j <= NJM; j++) {
        for (i = 2; i <= NIM; i++) {
            U[i][j] = U[i][j] + CC[i][j];
        }
    }

    //---------------------------- X_BC
    for (j = 2; j <= NJM; j++) {
        if (IBCW == IPRE || IBCW == IOUT) U[1][j] = U[2][j];
        if (IBCE == IPRE || IBCE == IOUT) U[NI][j] = U[NIM][j];
    }

    //---------------------------- Y-BC
    for (i = 2; i <= NIM; i++) {
        if (IBCS == ISYM) U[i][1] = U[i][2];
        if (IBCN == ISYM) U[i][NJ] = U[i][NJM];
    }

    //--------------------------- SS
    USMAX = 0.; UFMAX = 0.; UCMAX = 0.;
    for (j = 2; j <= NJM; j++) {
        for (i = 2; i <= NIM; i++) {
            SS[i][j] = USOR[i][j] - APC[i][j] * U[i][j]
                + AW[i][j] * U[i - 1][j] + AE[i][j] * U[i + 1][j]
                + AS[i][j] * U[i][j - 1] + AN[i][j] * U[i][j + 1];
            USMAX = fmax(USMAX, fabs(SS[i][j]) / APC[i][j]);
            UFMAX = fmax(UFMAX, fabs(U[i][j]));
            UCMAX = fmax(UCMAX, fabs(CC[i][j]));
        }
    }

    USMAX = USMAX / VELMAX; UCMAX = UCMAX / VELMAX;
}

void V_SOLVE()
{
    PXY_GET();

    for (j = 2; j <= NJM; j++) {
        for (i = 2; i <= NIM; i++) {


            double visw = VIS_INTP(VISP[i][j], VISP[i - 1][j]);
            double vise = VIS_INTP(VISP[i][j], VISP[i + 1][j]);
            double viss = VIS_INTP(VISP[i][j], VISP[i][j - 1]);
            double visn = VIS_INTP(VISP[i][j], VISP[i][j + 1]);


            AW[i][j] = visw * RU[i] * DYP[j] / DXU[i];
            AE[i][j] = vise * RU[i + 1] * DYP[j] / DXU[i + 1];
            AS[i][j] = viss * RP[i] * DXP[i] / DYV[j];
            AN[i][j] = visn * RP[i] * DXP[i] / DYV[j + 1];


            double flw = RHOP[i][j] * DYP[j] * RU[i] * UU[i][j];
            double fle = RHOP[i][j] * DYP[j] * RU[i + 1] * UU[i + 1][j];
            double fls = RHOP[i][j] * DXP[i] * RP[i] * VV[i][j];
            double fln = RHOP[i][j] * DXP[i] * RP[i] * VV[i][j + 1];


            AW[i][j] += fmax(0., flw);
            AE[i][j] += fmax(0., -fle);
            AS[i][j] += fmax(0., fls);
            AN[i][j] += fmax(0., -fln);


            VSOR[i][j] = 0.;

            double ac = 0.;
            AP[i][j] = AW[i][j] + AE[i][j] + AS[i][j] + AN[i][j] + ac;
            VSOR[i][j] -= VOL[i][j] * 0.5 * (PYV[i][j] + PYV[i][j + 1]);
            
            //if (IR == 1) AP[i][j] += VOL[i][j] * VISP[i][j] / RP[i] / RP[i];
            //USOR[i][j] -= VOL[i][j] * 0.5 * (PXU[i][j] + PXU[i + 1][j]);

            AP[i][j] = AP[i][j] / RELAXUV;
            VSOR[i][j] += (1. - RELAXUV) * AP[i][j] * V[i][j];
            APC[i][j] = AP[i][j];

            APV_P[i][j] = VOL[i][j] / AP[i][j];

        }
    }

    //------------------------------------ X-BC
    for (j = 2; j <= NJM; j++) {
        if (IBCW == ISYM) AP[2][j] -= AW[2][j];
        if (IBCE == ISYM) AP[NIM][j] -= AE[NIM][j];

    }

    //------------------------------------- Y-BC
    for (i = 2; i <= NIM; i++) {
        if (IBCS == IPRE || IBCS == IOUT) AP[i][2] -= AS[i][2];
        if (IBCN == IPRE || IBCN == IOUT) AP[i][NJM] -= AN[i][NJM];

        APV_P[i][1] = APV_P[i][2]; APV_P[i][NJ] = APV_P[i][NJM];
    }


    for (j = 2; j <= NJM; j++) {
        for (i = 2; i <= NIM; i++) {
            SS[i][j] = VSOR[i][j] - APC[i][j] * V[i][j]
                + AW[i][j] * V[i - 1][j] + AE[i][j] * V[i + 1][j]
                + AS[i][j] * V[i][j - 1] + AN[i][j] * V[i][j + 1];
        }
    }


    LNTDMA(2, NIM, 2, NJM);


    for (j = 2; j <= NJM; j++) {
        for (i = 2; i <= NIM; i++) {
            V[i][j] = V[i][j] + CC[i][j];
        }
    }

    //---------------------------- X_BC
    for (j = 2; j <= NJM; j++) {
        if (IBCW == ISYM) V[1][j] = V[2][j];
        if (IBCE == ISYM) V[NI][j] = V[NIM][j];
    }

    //---------------------------- Y-BC
    for (i = 2; i <= NIM; i++) {
        if (IBCS == IPRE || IBCS == IOUT) V[i][1] = V[i][2];
        if (IBCN == IPRE || IBCN == IOUT) V[i][NJ] = V[i][NJM];
    }

    //--------------------------- SS
    VSMAX = 0.; VFMAX = 0.; VCMAX = 0.;
    for (j = 2; j <= NJM; j++) {
        for (i = 2; i <= NIM; i++) {
            SS[i][j] = VSOR[i][j] - APC[i][j] * V[i][j]
                + AW[i][j] * V[i - 1][j] + AE[i][j] * V[i + 1][j]
                + AS[i][j] * V[i][j - 1] + AN[i][j] * V[i][j + 1];
            VSMAX = fmax(VSMAX, fabs(SS[i][j]) / APC[i][j]);
            VFMAX = fmax(VFMAX, fabs(V[i][j]));
            VCMAX = fmax(UCMAX, fabs(CC[i][j]));
        }
    }

    VSMAX = VSMAX / VELMAX; VCMAX = VCMAX / VELMAX;
}



void T_SOLVE()
{
    for (j = 1; j <= NJ; j++) { for (i = 1; i <= NI; i++) { BT[i][j] = T[i][j]; } }
    //
    for (j = 2; j <= NJM; j++) {
        for (i = 2; i <= NIM; i++) {
            double kw = TKK_INTP(TKKP[i][j], TKKP[i - 1][j]);
            double ke = TKK_INTP(TKKP[i][j], TKKP[i + 1][j]);
            double ks = TKK_INTP(TKKP[i][j], TKKP[i][j - 1]);
            double kn = TKK_INTP(TKKP[i][j], TKKP[i][j + 1]);
            AW[i][j] = kw * RU[i] * DYP[j] / DXU[i];
            AE[i][j] = ke * RU[i + 1] * DYP[j] / DXU[i + 1];
            AS[i][j] = ks * RP[i] * DXP[i] / DYV[j];
            AN[i][j] = kn * RP[i] * DXP[i] / DYV[j + 1];

            double flw = ROCP[i][j] * RU[i] * DYP[j] * U[i][j];
            double fle = ROCP[i][j] * RU[i + 1] * DYP[j] * U[i + 1][j];
            double fls = ROCP[i][j] * RP[i] * DXP[i] * V[i][j];
            double fln = ROCP[i][j] * RP[i] * DXP[i] * V[i][j + 1];
            AW[i][j] += fmax(0., flw);
            AE[i][j] += fmax(0., -fle);
            AS[i][j] += fmax(0., fls);
            AN[i][j] += fmax(0., -fln);

            double ac = 0.;
            AP[i][j] = AW[i][j] + AE[i][j] + AS[i][j] + AN[i][j] + ac;
            TSOR[i][j] = QDOT * RP[i] * DXP[i] * DYP[j] + ac * BT[i][j];
            APC[i][j] = AP[i][j];
        }
    }

    //---------------------------- X-BC
    for (j = 2; j <= NJM; j++) {
        if (IBCW == IWQQ) AP[2][j] -= AW[2][j];
        if (IBCE == IWQQ) AP[NIM][j] -= AE[NIM][j];

    }

    //---------------------------- Y-BC
    for (i = 2; i <= NIM; i++) {
        if (IBCS == IWQQ) AP[i][2] -= AS[i][2];
        if (IBCN == IWQQ) AP[i][NJM] -= AN[i][NJM];
    }

    for (j = 2; j <= NJM; j++) {
        for (i = 2; i <= NIM; i++) {
            SS[i][j] = TSOR[i][j] - APC[i][j] * T[i][j] + AW[i][j] * T[i - 1][j] + AE[i][j] * T[i + 1][j]
                + AS[i][j] * T[i][j - 1] + AN[i][j] * T[i][j + 1];
        }
    }
    //--------------------------------------------- iteration
    for (int iter = 0; iter <= ITMAX; iter++) {

        if (iter == iter / 1 * 1) {
            printf(" iter= %4d T= %11.3e ResMax= %11.3e Fcmax= %11.3e  \n", iter, T[10][10], SSMAX, FCMAX);
        }

        if (iter > 0 && SSMAX < ERRMAX && FCMAX < ERRMAX) {
            printf(" iter= %4d T= %11.3e ResMax= %11.3e Fcmax= %11.3e  \n", iter, T[10][10], SSMAX, FCMAX);
            break;
        }

        LNTDMA(2, NIM, 2, NJM);

        for (j = 2; j <= NJM; j++) {
            for (i = 2; i <= NIM; i++) {
                T[i][j] = T[i][j] + RELAX * CC[i][j];
            }
        }

        //---------------------------- X-BC
        for (j = 2; j <= NJM; j++) {
            if (IBCW == IWQQ) T[1][j] = T[2][j]; //+ QBCW * RU[2] * DYP[j] / AW[2][j];
            if (IBCE == IWQQ) T[NI][j] = T[NIM][j]; //+ QBCE * RU[NI] * DYP[j] / AE[NIM][j];
        }

        //---------------------------- Y-BC
        for (i = 2; i <= NIM; i++) {
            if (IBCS == IWQQ) T[i][1] = T[i][2]; //+ QBCS * RP[i] * DXP[i] / AS[i][2];
            if (IBCN == IWQQ) T[i][NJ] = T[i][NJM]; //+ QBCN * RP[i] * DXP[i] / AS[i][NJM];
        }
        //--------------------------------------- SS
        SSMAX = 0.; FFMAX = 0.; FCMAX = 0.;
        for (j = 2; j <= NJM; j++) {
            for (i = 2; i <= NIM; i++) {
                SS[i][j] = TSOR[i][j] - APC[i][j] * T[i][j] + AW[i][j] * T[i - 1][j] + AE[i][j] * T[i + 1][j]
                    + AS[i][j] * T[i][j - 1] + AN[i][j] * T[i][j + 1];
                SSMAX = fmax(SSMAX, fabs(SS[i][j]) / APC[i][j]);
                FFMAX = fmax(FFMAX, fabs(T[i][j]));
                FCMAX = fmax(FCMAX, fabs(CC[i][j]));
            }
        }
        FFMAX = fmax(1.e-30, FFMAX); SSMAX = SSMAX / FFMAX; FCMAX = FCMAX / FFMAX;
    }

}

void P_SOLVE()
{
    for (j = 2; j <= NJM; j++) {
        for (i = 1; i <= NI; i++) {
            U[i][j] += 0.5 * APU_P[i][j] * (PXU[i][j] + PXU[i + 1][j]); //dpdx ������ �κ�
        }
    }

    for (i = 2; i <= NIM; i++) {
        for (j = 1; j <= NJ; j++) {
            V[i][j] += 0.5 * APV_P[i][j] * (PYV[i][j] + PYV[i][j+1]); //dpdy ������ �κ�
        }
    }

    for (j = 2; j <= NJM; j++) {
        for (i = 1; i <= NI; i++) {
            APU_U[i][j] = APU_P[i][j] * DXUW[i] + APU_P[i - 1][j] * DXUE[i]; // interpolation
            UU[i][j] = U[i][j] * DXUW[i] + U[i - 1][j] * DXUE[i]; // velocity interpolation
        }
    }

    for (i = 2; i <= NIM; i++) {
        for (j = 1; j <= NJ; j++) {
            APV_V[i][j] = APV_P[i][j] * DYVS[j] + APV_V[i][j-1] * DYVN[j]; // interpolation
            VV[i][j] = V[i][j] * DYVS[j] + V[i][j-1] * DYVN[j]; //velocity interpolation
        }
    }

    for (j = 2; j <= NJM; j++) {
        for (i = 2; i <= NIM; i++) {


            AW[i][j] = APU_U[i][j] * RU[i] * DYP[j] / DXU[i];
            AE[i][j] = APU_U[i+1][j] * RU[i + 1] * DYP[j] / DXU[i + 1];
            AS[i][j] = APV_V[i][j] * RP[i] * DXP[i] / DYV[j];
            AN[i][j] = APV_V[i][j+1] * RP[i] * DXP[i] / DYV[j + 1];
            AP[i][j] = AW[i][j] + AE[i][j] + AS[i][j] + AN[i][j];


            double flw = RU[i    ] * DYP[j] * UU[i][j];
            double fle = RU[i + 1] * DYP[j] * UU[i + 1][j];
            double fls = RP[i    ] * DXP[i] * VV[i][j];
            double fln = RP[i    ] * DXP[i] * VV[i][j + 1];

            PSOR[i][j] = flw - fle + fls - fln;
            PC[i][j] = 0.;
            APC[i][j] = AP[i][j];
                            
        }
    }


    //------------------------------------ X-BC
    for (j = 2; j <= NJM; j++) {
        if (IBCW != IPRE) AP[  2][j] -= AW[  2][j];
        if (IBCE != IPRE) AP[NIM][j] -= AE[NIM][j];

    }

    //------------------------------------- Y-BC
    for (i = 2; i <= NIM; i++) {
        if (IBCS != IPRE) AP[i][  2] -= AS[i][2];
        if (IBCN != IPRE) AP[i][NJM] -= AN[i][NJM];

        APV_P[i][1] = APV_P[i][2]; APV_P[i][NJ] = APV_P[i][NJM];

    }


    for (j = 2; j <= NJM; j++) {
        for (i = 2; i <= NIM; i++) {
            SS[i][j] = PSOR[i][j] - APC[i][j] * P[i][j]
                + AW[i][j] * P[i - 1][j] + AE[i][j] * P[i + 1][j]
                + AS[i][j] * P[i][j - 1] + AN[i][j] * P[i][j + 1];
        }
    }

    // -------------------------DEBUG-------------------------
    // std::cout<<"APU_U: "<<std::endl;
    // for(int i=2; i<=NIM; i++) {
    //     for(int j=2; j<=NJM; j++) {
    //         std::cout<<APU_U[i][j]<<" ";
    //     }
    //     std::cout<<std::endl;
    // }
    // std::cout<<"APU_V: "<<std::endl;
    // for(int i=2; i<=NIM; i++) {
    //     for(int j=2; j<=NJM; j++) {
    //         std::cout<<APV_V[i][j]<<" ";
    //     }
    //     std::cout<<std::endl;
    // }
    // std::cout<<"AP: "<<std::endl;
    // for(int i=2; i<=NIM; i++) {
    //     for(int j=2; j<=NJM; j++) {
    //         std::cout<<AP[i][j]<<" ";
    //     }
    //     std::cout<<std::endl;
    // }
    // std::cout<<"SS: "<<std::endl;
    // for(int i=2; i<=NIM; i++) {
    //     for(int j=2; j<=NJM; j++) {
    //         std::cout<<SS[i][j]<<" ";
    //     }
    //     std::cout<<std::endl;
    // }
    // std::cout<<"V: "<<std::endl;
    // for(int i=1; i<=NI; i++) {
    //     for(int j=1; j<=NJ; j++) {
    //         std::cout<<V[i][j]<<" ";
    //     }
    //     std::cout<<std::endl;
    // }
    std::cout.precision(3);
    std::cout<<"P: "<<std::endl;
    for(int i=1; i<=NI; i++) {
        for(int j=1; j<=NJM; j++) {
            std::cout<<P[i][j]<<" ";
        }
        std::cout<<std::endl;
    }
    std::cout.precision(3);
    std::cout<<"APV_V: "<<std::endl;
    for(int i=1; i<=NI; i++) {
        for(int j=1; j<=NJ; j++) {
            std::cout<<APV_V[i][j]<<" ";
        }
        std::cout<<std::endl;
    }
    std::cout.precision(3);
    std::cout<<"APU_U: "<<std::endl;
    for(int i=1; i<=NI; i++) {
        for(int j=1; j<=NJ; j++) {
            std::cout<<APU_U[i][j]<<" ";
        }
        std::cout<<std::endl;
    }
    // -------------------------DEBUG-------------------------


    LNTDMA(2, NIM, 2, NJM);


    for (j = 2; j <= NJM; j++) {
        for (i = 2; i <= NIM; i++) {
            PC[i][j] = PC[i][j] + CC[i][j];
        }
    }

    //---------------------------- X_BC
    for (j = 2; j <= NJM; j++) {
        if (IBCW != IPRE) PC[ 1][j] = PC[  2][j];
        if (IBCE != IPRE) PC[NI][j] = PC[NIM][j];
    }

    //---------------------------- Y-BC
    for (i = 2; i <= NIM; i++) {
        if (IBCS != IPRE) PC[i][ 1] = PC[i][  2];
        if (IBCN != IPRE) PC[i][NJ] = PC[i][NJM];
    }

    for (j = 1; j <= NJ; j++) {
        for (i = 1; i <= NI; i++) {
            P[i][j] += RELAXP * PC[i][j];
        }
    }

    //--------------------------- SS
    PSMAX = 0.; PFMAX = 0.; PCMAX = 0.;
    for (j = 2; j <= NJM; j++) {
        for (i = 2; i <= NIM; i++) {
            SS[i][j] = PSOR[i][j] - APC[i][j] * P[i][j]
                + AW[i][j] * P[i - 1][j] + AE[i][j] * P[i + 1][j]
                + AS[i][j] * P[i][j - 1] + AN[i][j] * P[i][j + 1];
            PSMAX = fmax(PSMAX, fabs(SS[i][j]) / APC[i][j]);
            PFMAX = fmax(PFMAX, fabs(P[i][j]));
            PCMAX = fmax(PCMAX, fabs(CC[i][j]));
        }
    }

    PFMAX = fmax(1., PFMAX); PSMAX = PSMAX / PFMAX; PCMAX = PCMAX / PFMAX;

    PXY_GET();

    for (j = 2; j <= NJM; j++) {
        for (i = 2; i <= NI; i++) {
            UU[i][j] -=  APU_U[i][j] * PXU[i][j]; // UU update
        }
    }

    for (i = 2; i <= NIM; i++) {
        for (j = 2; j <= NJ; j++) {
            VV[i][j] -= APV_V[i][j] * PYV[i][j]; // VV update
        }
    }

    for (j = 2; j <= NJM; j++) {
        for (i = 1; i <= NI; i++) {
            U[i][j] -= 0.5 * APU_P[i][j] * (PXU[i][j] + PXU[i + 1][j]); // U update
        }
    }

    for (i = 2; i <= NIM; i++) {
        for (j = 1; j <= NJ; j++) {
            V[i][j] -= 0.5 * APV_P[i][j] * (PYV[i][j] + PYV[i][j + 1]); // V update
        }
    }



}


//double TKK_INTP(double a, double b) { return 0.5 * (a + b); }

/*
void INIT()
{
   NIM=NI-1;NJM=NJ-1;

   for(i=2; i<=NI;  i++) XU[i]=DX0*(i-2.);
   for(i=2; i<=NIM; i++) XP[i]=0.5*(XU[i+1]+XU[i]); XP[1]=XU[2]; XP[NI]=XU[NI];
   for(i=2; i<=NIM; i++) DXP[i]=XU[i+1]-XU[i];
   for(i=2; i<=NI;  i++) DXU[i]=XP[i]-XP[i-1];

   for(j=2; j<=NJ;  j++) YV[j]=DY0*(j-2.);
   for(j=2; j<=NJM; j++) YP[j]=0.5*(YV[j+1]+YV[j]); YP[j]=YV[2]; YP[NJ]=YV[NJ];
   for(j=2; j<=NJM; j++) DYP[j]=YV[j+1]-YV[j];
   for(j=2; j<=NJ;  j++) DYV[j]=YP[j]-YP[j-1];

   for(i=2; i<=NI; i++) RU[i]=1.;
   for(i=1; i<=NI; i++) RP[i]=1.;

   if(IR==1){
      for(i=2; i<=NI; i++) RU[i]=fabs(XU[i]);
      for(i=1; i<=NI; i++) RP[i]=fabs(XP[i]);
   }

   for(i=1; i<=NI; i++){
      for(j=1; j<=NJ; j++){
         TKKP[i][j]=K0; ROCP[i][j]=ROCP0;
      }
   }

} */

void LNTDMA(int ist, int iend, int jst, int jend)
{
    double beta, sss, gama[IJDIM];

    for (i = ist - 1; i <= iend + 1; i++) { for (j = jst - 1; j <= jend + 1; j++) { CC[i][j] = 0.; } }
    //
    for (int nt = 1; nt <= 2; nt++) {
        //-------------------------------------- x-TDMA
        for (j = jst; j <= jend; j++) {
            i = ist; sss = SS[i][j] + AS[i][j] * CC[i][j - 1] + AN[i][j] * CC[i][j + 1];
            beta = AP[i][j]; CC[i][j] = sss / beta;
            for (i = ist + 1; i <= iend; i++) {
                sss = SS[i][j] + AS[i][j] * CC[i][j - 1] + AN[i][j] * CC[i][j + 1];
                gama[i] = -AE[i - 1][j] / beta; beta = AP[i][j] + AW[i][j] * gama[i];
                CC[i][j] = (sss + AW[i][j] * CC[i - 1][j]) / beta;
            }
            for (i = iend - 1; i >= ist; i--) { CC[i][j] = CC[i][j] - gama[i + 1] * CC[i + 1][j]; }
        }
        //------------------------------------- y-TDMA
        for (i = ist; i <= iend; i++) {
            j = jst; sss = SS[i][j] + AW[i][j] * CC[i - 1][j] + AE[i][j] * CC[i + 1][j];
            beta = AP[i][j]; CC[i][j] = sss / beta;
            for (j = jst + 1; j <= jend; j++) {
                sss = SS[i][j] + AW[i][j] * CC[i - 1][j] + AE[i][j] * CC[i + 1][j];
                gama[j] = -AN[i][j - 1] / beta; beta = AP[i][j] + AS[i][j] * gama[j];
                CC[i][j] = (sss + AS[i][j] * CC[i][j - 1]) / beta;
            }
            for (j = jend - 1; j >= jst; j--) { CC[i][j] = CC[i][j] - gama[j + 1] * CC[i][j + 1]; }
        }
    }
}

void Plot_Contour(FILE* fp, double f0, int ibeg, int iend, int jbeg, int jend, double x[IDIM], double y[JDIM], double f[IDIM][JDIM])
{
    int i0, j0, i1, j1, n;
    double x1[4], y1[4], df00, df10, df11, df01;
    for (j0 = jbeg; j0 < jend; j0++) {
        for (i0 = ibeg; i0 < iend; i0++) {
            j1 = j0 + 1; i1 = i0 + 1; n = -1;
            df00 = f0 - f[i0][j0]; df10 = f0 - f[i1][j0];
            df11 = f0 - f[i1][j1]; df01 = f0 - f[i0][j1];

            if ((df00 * df10 <= 0.) && (df00 != df10)) {
                n++; x1[n] = x[i0] + df00 / (df00 - df10) * (x[i1] - x[i0]); y1[n] = y[j0];
            }

            if ((df10 * df11 <= 0.) && (df10 != df11)) {
                n++; x1[n] = x[i1]; y1[n] = y[j0] + df10 / (df10 - df11) * (y[j1] - y[j0]);
            }

            if ((df11 * df01 <= 0.) && (df11 != df01)) {
                n++; x1[n] = x[i0] + df01 / (df01 - df11) * (x[i1] - x[i0]); y1[n] = y[j1];
            }

            if ((df01 * df00 <= 0.) && (df01 != df00)) {
                n++; x1[n] = x[i0]; y1[n] = y[j0] + df00 / (df00 - df01) * (y[j1] - y[j0]);
            }

            if (n >= 1) {
                for (int nn = 0; nn < n; nn++) {
                    fprintf(fp, "  %9.5f %9.5f \n", x1[nn], y1[nn]);
                    fprintf(fp, "  %9.5f %9.5f \n \n", x1[nn + 1], y1[nn + 1]);
                }
                if (n != 1) {
                    fprintf(fp, "  %9.5f %9.5f \n", x1[n], y1[n]);
                    fprintf(fp, "  %9.5f %9.5f \n \n", x1[0], y1[0]);
                }
            }
        }
    }
}

