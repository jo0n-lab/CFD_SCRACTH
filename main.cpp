#include <iostream>
#define COMMON_DEF
#include "common.h" 


int main()
{
	NI = 22; DX0 = 0.05; NJ = 22; DY0 = 0.05; U0 = 0.; V0 = -1.e-6;
	IR = 0; JUMP = 1;
	K0 = 1.; ROCP0 = 1.e+6; QDOT = 0.;
	ITMAX = 1000; ERRMAX = 1.e-5; RELAX = 1.;

	//IBCW = IWAL; QBCW = 0.;
	//IBCE = IWQQ; QBCE = 0.;
	//IBCS = IWQQ; TBCS = 0.;
	//IBCN = IWAL; TBCN = 1.;

	IBCW = IWAL; TBCW = 0.;
	IBCE = IWQQ; QBCE = 0.;
	IBCS = IWQQ; QBCS = 0.; 
	IBCN = IWAL; TBCN = 1.;


	for (i = 1; i < NI; i++) { for (j = 1; j < NJ; j++) { T[i][j] = 0.; } };
	

	//for (i = 1; i < NI; i++) { T[i][NJ] = 1; } //case2
	//for (j = 1; j < NJ; j++) { T[NI][j] = 1; } //case1
	if (IBCW == IWAL) for (j = 1; j < NJ; j++) { T[1][j] = TBCW; }
	if (IBCE == IWAL) for (j = 1; j < NJ; j++) { T[NI][j] = TBCE; }
	if (IBCS == IWAL) for (i = 1; i < NI; i++) { T[i][1] = TBCS; }
	if (IBCN == IWAL) for (i = 1; i < NI; i++) { T[i][NJ] = TBCN; }

	INIT(); 

	for (i = 2; i < NI; i++) { for (j = 1; j < NJ; j++) { U[i][j] = 1.e-5 * XU[i]; } };
	for (i = 1; i < NI; i++) { for (j = 2; j < NJ; j++) { V[i][j] = -1.e-5 * YV[j]; } };

	T_SOLVE();

	FILE *fp = fopen("out", "w");

	

	for (j = 1; j <= NJ; j++) {
		double Pe = ROCP0 * V0;
		double Tex = (exp(Pe * YP[j]) - 1.) / (exp(Pe) - 1.);
		fprintf(fp, " %4d %15.7e %15.7e  %15.7e %15.7e \n", j, YP[j], T[7][j], T[17][j], Tex);
	}
		

	fclose(fp);

	FILE *fp1 = fopen("out", "w");
	//for(i=0; i <= 10; i++) {}
	Plot_Contour(fp1, 0.9, 1, NI, 1, NJ, XP, YP, T);
	Plot_Contour(fp1, 0.8, 1, NI, 1, NJ, XP, YP, T);
	Plot_Contour(fp1, 0.7, 1, NI, 1, NJ, XP, YP, T);
	Plot_Contour(fp1, 0.6, 1, NI, 1, NJ, XP, YP, T);
	Plot_Contour(fp1, 0.5, 1, NI, 1, NJ, XP, YP, T);
	Plot_Contour(fp1, 0.4, 1, NI, 1, NJ, XP, YP, T);
	Plot_Contour(fp1, 0.3, 1, NI, 1, NJ, XP, YP, T);
	Plot_Contour(fp1, 0.2, 1, NI, 1, NJ, XP, YP, T);
	Plot_Contour(fp1, 0.1, 1, NI, 1, NJ, XP, YP, T);
	fclose(fp1);


	FILE *fp2 = fopen("T_out", "w");
	for (i = 1; i <= NI; i++) {
		for (j = 1; j <= NJ; j++) {
			fprintf(fp, "%15.7e ", T[i][j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp2);
	

	return 0;
}

