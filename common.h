#pragma once
#include <math.h>

#define IDIM 100
#define JDIM 100
#define IJDIM 100

#define IINF 1
#define ISYM 2
#define IPRE 3
#define IWAL 4
#define IWQQ 5

#ifdef COMMON_DEF
#define EXT
#else
#define EXT extern
#endif 

EXT int IR;
EXT int NI, NIM, NJ, NJM;
EXT int i, j, k;

EXT double DX0, XP[IDIM], XU[IDIM], RP[IDIM], RU[IDIM], DXP[IDIM], DXU[IDIM];
EXT double DY0, YP[JDIM], YV[JDIM], DYP[JDIM], DYV[JDIM];
EXT double T[IDIM][JDIM], TSOR[IDIM][JDIM], BT[IDIM][JDIM];
EXT double K0, ROCP0, QDOT, TKKP[IDIM][JDIM], ROCP[IDIM][JDIM];


EXT double U0, U[IDIM][JDIM], USOR[IDIM][JDIM];
EXT double V0, V[IDIM][JDIM], VSOR[IDIM][JDIM];

EXT int IBCW, IBCE, IBCS, IBCN;
EXT double TBCW, TBCE, TBCS, TBCN;
EXT double QBCW, QBCE, QBCS, CBCN;

EXT int LTIME, ITIME, JUMP, ITMAX;
EXT double FFMAX, FCMAX, SSMAX, ERRMAX, RELAX;

EXT double CC[IDIM][JDIM], SS[IDIM][JDIM];
EXT double AP[IDIM][JDIM], AW[IDIM][JDIM], AE[IDIM][JDIM], AS[IDIM][JDIM], AN[IDIM][JDIM], APC[IDIM][JDIM];

//EXT double gama[IJDIM];


void LNTDMA(int ist, int iend, int jst, int jend);
void INIT();
void T_SOLVE();
double TKK_INTP(double a, double b);

void Plot_Contour(FILE* fp, double f0, int ibeg, int iend, int jbeg, int jend, double x[IDIM], double y[JDIM], double f[IDIM][JDIM]);

