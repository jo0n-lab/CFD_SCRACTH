#include <math.h>

#define IDIM 200
#define JDIM 200
#define IJDIM 200

#define IINF 1
#define ISYM 2
#define IPRE 3
#define IWAL 4
#define IWQQ 5
#define IOUT 6

#ifdef COMMON_DEF
#define EXT
#else
#define EXT extern
#endif 

EXT int IR;
EXT int NI, NIM, NJ, NJM;
EXT int i, j, k;

EXT double DX0, XP[IDIM], XU[IDIM], DXP[IDIM], DXU[IDIM], DXUW[IDIM], DXUE[IDIM];
EXT double DY0, YP[JDIM], YV[JDIM], DYP[JDIM], DYV[JDIM], DYVS[JDIM], DYVN[JDIM];
EXT double RP[IDIM], RU[IDIM], VOL[IDIM][JDIM];

EXT double T0, T[IDIM][JDIM], TSOR[IDIM][JDIM], BT[IDIM][JDIM];
EXT double U0, U[IDIM][JDIM], USOR[IDIM][JDIM], BU[IDIM][JDIM], UU[IDIM][JDIM], PXU[IDIM][JDIM];
EXT double V0, V[IDIM][JDIM], VSOR[IDIM][JDIM], BV[IDIM][JDIM], VV[IDIM][JDIM], PYV[IDIM][JDIM];
EXT double P0, P[IDIM][JDIM], PSOR[IDIM][JDIM], PC[IDIM][JDIM];

EXT double CONV_SOR[IDIM][JDIM], SOR[IDIM][JDIM];

EXT double PCXU[IDIM][JDIM], UC[IDIM][JDIM], HUCP[IDIM][JDIM], HUCU[IDIM][JDIM];
EXT double PCYV[IDIM][JDIM], VC[IDIM][JDIM], HVCP[IDIM][JDIM], HVCV[IDIM][JDIM];

EXT double APU_P[IDIM][JDIM], APV_P[IDIM][JDIM], APU_U[IDIM][JDIM], APV_V[IDIM][JDIM];

EXT double K0, ROCP0, QDOT, TKKP[IDIM][JDIM], ROCP[IDIM][JDIM];
EXT double RHOP0, RHOP[IDIM][JDIM];
EXT double VIS0, VISP[IDIM][JDIM];

EXT double K0_SOL, ROCP0_SOL, VIS0_SOL;

EXT double GRAVY, BETA;

EXT double SF[IDIM][JDIM];


EXT int LTIME, ITIME, JUMP, ITMAX;
EXT double FFMAX, FCMAX, SSMAX, ERRMAX, RELAX;
EXT double UFMAX, UCMAX, USMAX, RELAXUV;
EXT double VFMAX, VCMAX, VSMAX;
EXT double PFMAX, PCMAX, PSMAX, RELAXP;
EXT double TFMAX, TCMAX, TSMAX, RELAXT;
EXT double VELMAX;

EXT int IBCW, IBCE, IBCS, IBCN;
EXT double TBCW, TBCE, TBCS, TBCN;
EXT double QBCW, QBCE, QBCS, QBCN;
EXT double UBCW, UBCE, UBCS, UBCN;
EXT double VBCW, VBCE, VBCS, VBCN;
EXT double PBCW, PBCE, PBCS, PBCN;

EXT double CC[IDIM][JDIM], SS[IDIM][JDIM];
EXT double AP[IDIM][JDIM], AW[IDIM][JDIM], AE[IDIM][JDIM], AS[IDIM][JDIM], AN[IDIM][JDIM], APC[IDIM][JDIM];

//EXT double gama[IJDIM];


void LNTDMA(int ist, int iend, int jst, int jend);
void INIT();
void PXY_GET();
void T_SOLVE(); void T_SOLVE1();
void U_SOLVE(); void V_SOLVE(); void P_SOLVE();
void SF_Get();
void PROP();

double TKK_INTP(double a, double b);
double VIS_INTP(double a, double b);

void Plot_Contour(FILE* fp, double f0, int ibeg, int iend, int jbeg, int jend, double x[IDIM], double y[JDIM], double f[IDIM][JDIM]);

