#include <iostream>
#define COMMON_DEF 
#include "common.h"


int main()
{
	NI = 42; DX0 = 0.1 / (NI - 2.); NJ = 42; DY0 = 0.1 / (NJ - 2.); U0 = 0.; V0 = 0.; P0 = 0.; T0 = 0.;
	IR = 0.; JUMP = 100.; GRAVY = 9.81; BETA = 1.e-3;
	DTIME = 1.; I_NC = 1.; I_STEADY = 0.; LTIME = 200.;		// unsteady
	K0 = 1.e-2; ROCP0 = 1.e2; QDOT = 0.; RHOP0 = 1.e0; VIS0 = 0.71e-4;
	K0_S = 10.; ROCP0_S = 1.e6; RHOP0_S = 1.e3; VIS0_S = 10.;

	ITMAX = 10000.; ERRMAX = 1.e-5; RELAX = 1.; RELAXT = 1.; RELAXUV = 0.8; RELAXP = 1.;

	double Ra = 1.e4, DT = Ra / 1.382e3;

	IBCW = IWAL; TBCW = 0.5 * DT; QBCW = 0.; UBCW = 0.; VBCW = 0.; PBCW = 0.;
	IBCE = IWAL; TBCE = -0.5 * DT; QBCE = 0.; UBCE = 0.; VBCE = 0.; PBCE = 0.;
	IBCS = IWQQ; TBCS = 0.; QBCS = 0.; UBCS = 0.; VBCS = 0.; PBCS = 0.;
	IBCN = IWQQ; TBCN = 0.; QBCN = 0.; UBCN = 0.; VBCN = 0.; PBCN = 0.;
	//9-1

	INIT(); //F_Get(0.03, 0.03, 0.01, 0.07, 0.07, 0.01); H_Get(DX0); PROP();

	// 초기 속도 초기화 및 설정
	//for (i = 2; i < NI; i++) { for (j = 1; j < NJ; j++) { U[i][j] = U0 * XU[i]; } };
	//for (j = 1; j < NJ; j++) { U[1][j] = 0; V[1][j] = 0; };

	printf(" iter       U           V           P           T         uvsMAX      PSMAX       TSMAX       uvCMAX       PCMAX       TCMAX \n");
	printf(" %4d %11.3e %11.3e %11.3e %11.3e %11.3e %11.3e %11.3e  %11.3e %11.3e %11.3e \n\n",
		0., U[10][10], V[10][10], P[10][10], T[10][10], 0., 0., 0., 0., 0., 0.);

	int last = LTIME;
	if (I_STEADY == 1) last = ITMAX;

	FILE* fpgnu = _fsopen("plt", "w", _SH_DENYNO);
	fprintf(fpgnu, " set xr [0:0.1] \n");
	fprintf(fpgnu, " set yr [0:0.1] \n");
	fprintf(fpgnu, " set grid \n");

	for (int itime = 0; itime <= last; itime++) {
		VELMAX = 1.e-20;

		for (j = 1; j <= NJ; j++) {
			for (i = 1; i <= NI; i++) {
				BU[i][j] = U[i][j]; BV[i][j] = V[i][j]; BT[i][j] = T[i][j];
				VELMAX = fmax(VELMAX, U[i][j]); VELMAX = fmax(VELMAX, V[i][j]);
				//UU[i][j] = 0.; VV[i][j] = 0.1 * XP[i];
				//V[i][j] = 0.1 * XP[i];
			}
		}

		PXY_GET(); U_SOLVE(); V_SOLVE(); P_SOLVE();

		if (I_NC == 1 || I_STEADY == 0) T_SOLVE();

		double uvSMAX = fmax(USMAX, VSMAX); double uvCMAX = fmax(UCMAX, VCMAX);
		SSMAX = fmax(uvSMAX, PSMAX); FCMAX = fmax(uvCMAX, PCMAX);
		SSMAX = fmax(SSMAX, TSMAX); FCMAX = fmax(FCMAX, TCMAX);

		int jump = 1;
		if (I_STEADY == 1) jump = JUMP;
		if (itime == itime / jump * jump) {
			printf(" iter       U           V           P           T         uvsMAX      PSMAX       TSMAX       uvCMAX       PCMAX       TCMAX \n");
			printf(" %4d %11.3e %11.3e %11.3e %11.3e %11.3e %11.3e %11.3e  %11.3e %11.3e %11.3e \n\n",
				itime, U[10][10], V[10][10], P[10][10], T[10][10], uvSMAX, PSMAX, TSMAX, uvCMAX, PCMAX, TCMAX);
		}

		if (I_STEADY == 1 && itime > 0 && SSMAX < ERRMAX && FCMAX < ERRMAX) {
			printf(" iter       U           V           P           T         uvsMAX      PSMAX       TSMAX       uvCMAX       PCMAX       TCMAX \n");
			printf(" %4d %11.3e %11.3e %11.3e %11.3e %11.3e %11.3e %11.3e  %11.3e %11.3e %11.3e \n\n",
				itime, U[10][10], V[10][10], P[10][10], T[10][10], uvSMAX, PSMAX, TSMAX, uvCMAX, PCMAX, TCMAX);
			break;
		}

		if (itime == itime / 10 * 10) {
			char Text[5];
			char TextS[10] = "Scont"; strcat(TextS, _itoa(itime / 10, Text, 10));
			char TextT[10] = "Tcont"; strcat(TextT, _itoa(itime / 10, Text, 10));
			fprintf(fpgnu, " plot \"%s\" w l, \"%s\" w l lt 3\n", TextS, TextT);
			fprintf(fpgnu, " pause -1 \n");

			FILE* fpS = _fsopen(TextS, "w", _SH_DENYNO);
			FILE* fpT = _fsopen(TextT, "w", _SH_DENYNO);

			SF_Get();

			//double theta[IDIM][JDIM];

			double SFmax = -1.e10, SFmin = 1.e10;
			for (j = 2; j <= NJ; j++) {
				for (i = 2; i <= NI; i++) {
					SFmax = fmax(SFmax, SF[i][j]);
					SFmin = fmin(SFmin, SF[i][j]);
				}
			}

			for (i = 1; i <= 19; i++) {
				double SFout = SFmin + (SFmax - SFmin) / 20. * i;
				Plot_Contour(fpS, SFout, 2, NI, 2, NJ, XU, YV, SF);
			}

			for (i = 1; i <= 19; i++) {
				double Tout = -0.5 * DT + DT / 20. * i;
				Plot_Contour(fpT, Tout, 1, NI, 1, NJ, XP, YP, T);
			}

			fclose(fpS); fclose(fpT);
		}
	}

	if (I_NC == 0 && I_STEADY == 1) T_SOLVE();

	
	FILE* fp = _fsopen("out", "w", _SH_DENYNO);
	FILE* fp2 = _fsopen("out1", "w", _SH_DENYNO);

	double qsum = 0., Nu = 0.;

	for (j = 2; j <= NJM; j++) {
		double q = K0 * (T[1][j] - T[2][j]) / DXU[2];
		qsum += q * DYP[j];

		fprintf(fp2, " %4d %15.7e %15.7e  %15.7e %15.7e  %15.7e\n",
			j, YP[j], U[NI / 2][j], U[NIM][j], T[NI / 2][j], T[NIM][j]);
	}

	
	
	fclose(fp2);

	Nu = qsum / DT / K0;

	SF_Get();

	//double theta[IDIM][JDIM];

	double SFmax = -1.e10, SFmin = 1.e10;
	for (j = 2; j <= NJ; j++) {
		for (i = 2; i <= NI; i++) {
			SFmax = fmax(SFmax, SF[i][j]);
			SFmin = fmin(SFmin, SF[i][j]);
		}
	}

	printf("\n  Nuav = %15.7e  SFmax = %15.7e  SFmin = %15.7e \n", Nu, SFmax, SFmin);

	/*
	for (i = 1; i <= NI; i++) {
		double Rex = RHOP0 * U0 * XP[i] / VIS0;
		double tau = VIS0 * (U[i][2] - U[i][1]) / DYV[2];
		double Cfx = tau * 2. / RHOP0 / U0 / U0;
		double Cfxex = 0.664 / sqrt(Rex);

		double h = QBCS / (T[i][1] - T0);
		double Nux = h * XP[i] / K0;
		double Nuxex = 0.453 * sqrt(Rex);

		fprintf(fp, " %4d	%15.7e	%15.7e	%15.7e	%15.7e	%15.7e	%15.7e\n",
			i, XP[i], V[i][NJM], Cfx, Cfxex, Nux, Nuxex);
	}

	printf("file name: %s", NAME_I);
	printf("            i          XP[i]          V[i][NJM]          Cfx         Cfxex         Nux         Nuxex\n");
	*/
	fclose(fp);


	// 등온선 contour plot
	FILE* fp1 = _fsopen("Tcont", "w", _SH_DENYNO);
	for (i = 1; i <= 19; i++) {
		double Tout = -0.5 * DT + DT / 20. * i;
		Plot_Contour(fp1, Tout, 1, NI, 1, NJ, XP, YP, T);
	}

	fclose(fp1);

	FILE* fp3 = _fsopen("SFcont", "w", _SH_DENYNO);
	for (i = 1; i <= 19; i++) {
		double SFout = SFmin + (SFmax - SFmin) / 20. * i;
		Plot_Contour(fp3, SFout, 2, NI, 2, NJ, XU, YV, SF);
	}

	fclose(fp3);

	FILE* fp4 = _fsopen("Fcont", "w", _SH_DENYNO);
	Plot_Contour(fp4, 0, 1, NI, 1, NJ, XP, YP, F);
	fclose(fp4);
	return 0;
}

