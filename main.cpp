#include <iostream>
#define COMMON_DEF
#include "common.h"               

int main()
{
	NI = 22, DX0 = 0.1/(NI-2), NJ = 22, DY0 = 0.1/(NJ-2), T0 = 0., P0 = 0., U0 = 0., V0 = 0., P0 = 0.; //NI=82
	INC=1;
	IR = 0, JUMP = 100, GRAVY = 9.81, BETA = 1.e-2;
	K0 = 1.e-2, ROCP0 = 1.e2, QDOT = 0., RHOP0 = 1.e0, VIS0 = 0.71e-4;
	double Ra = 1.e4, DTemp = Ra / 1.382e4;
	ITMAX = 2000, ERRMAX = 1.e-5, RELAX = 1., RELAXUV = 0.9, RELAXP = 1.0; RELAXT = 0.9;
	

	IBCW = IWAL, TBCW =  0.5 * DTemp, 	QBCW = 0., UBCW = 0., VBCW = 0.,  PBCW = 0.;
	IBCE = IWAL, TBCE = -0.5 * DTemp, 	QBCE = 0., UBCE = 0., VBCE = 0.,  PBCE = 0.;
	IBCS = IWQQ, TBCS =  0., 				QBCS = 0., UBCS = 0., VBCS = 0.,  PBCS = 0.;
	IBCN = IWQQ, TBCN =  0., 				QBCN = 0., UBCN = 0., VBCN = 0.,  PBCN = 0.;

	INIT();

	printf("DTemp= %10.3e T0= %10.3e Ra= %10.3e TBCW= %10.3e TBCE= %10.3e\n", DTemp, T0, Ra, TBCW, TBCE);
	printf("%4d %10.2e %10.2e %10.2e %10.2e\n", 0, U[10][10], V[10][10], P[10][10], T[10][10]);

	for (int iter = 0; iter <= ITMAX; iter++)
	{
		VELMAX = 1.e-20;
		for (j = 1; j <= NJ; j++)
		{
			for (i = 1; i <= NI; i++)
			{
				// continuity?
				BU[i][j] = U[i][j];
				BV[i][j] = V[i][j];
				BT[i][j] = T[i][j];
				VELMAX = fmax(VELMAX, U[i][j]);
				VELMAX = fmax(VELMAX, V[i][j]);
			}
		}

		if(INC==1) T_SOLVE();
		U_SOLVE(); V_SOLVE(); P_SOLVE();

		// ex 4.1 main.cpp
		double uvsmax=fmax(USMAX,VSMAX);
		double uvcmax=fmax(UCMAX,VCMAX);
		SSMAX = fmax(uvsmax, PSMAX);
		FCMAX = fmax(uvcmax, PCMAX);

		SSMAX = fmax(SSMAX, TSMAX);
		FCMAX = fmax(FCMAX, TCMAX);

		if (iter == iter / JUMP * JUMP)
		{
			printf("%4d %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e\n", iter, U[10][10], V[10][10], P[10][10], T[10][10] ,uvcmax, PCMAX, TCMAX, uvsmax, PSMAX, TSMAX);

		}
		if (iter > 0 && SSMAX < ERRMAX && FCMAX < ERRMAX)
		{
			printf("%4d %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e\n", iter, U[10][10], V[10][10], P[10][10], T[10][10] ,uvcmax, PCMAX, TCMAX, uvsmax, PSMAX, TSMAX);
			break;
		}
	}

	if(INC==0) T_SOLVE();

	FILE *fpT = fopen("out_T", "w");
	for (int i=1;i<=NI;i++){
		for (int j=1;j<=NJ;j++){
			fprintf(fpT,"%15.7e ",T[i][j]);
		}
		fprintf(fpT,"\n");
	}
	fclose(fpT);


	// for(j=1;j<=NJ;j++){
	// 	double vrsum=0;
	// 	double vrTsum=0;

	// 	for(i=2;i<=NIM;i++){
	// 		vrsum+=V[i][j]*RP[i]*DXP[i];
	// 		vrTsum+=V[i][j]*RP[i]*DXP[i]*T[i][j];
	// 	}
	// 	Tm[j]=vrTsum/vrsum;
	// }


	// FILE *fp1 = fopen("out_x", "w");
	// for (i = 2; i <= NIM; i++)
	// {
	// 	// double Vex = 0.2*(1-XP[i]*XP[i]/XP[NI]/XP[NI]);
	// 	// Ts=TBCE
	// 	// double theta = (T[NI][NJM] - T[i][NJM])/(T[NI][NJM] - Tm[NJM]);
	// 	// double eta = XP[i] / XP[NI];
	// 	// double thetaex = 4.36 * (3/8. - pow(eta,2)/2 + pow(eta,4)/8);
	// 	// double Rex = XP[i] / VIS0;
	// 	double Rex = RHOP0 * U0 * XP[i] / VIS0;
	// 	double Cfex = 0.664 * pow(Rex, -0.5);
	// 	// double Nuex = 0.332 * pow(Rex, 0.5);
	// 	double Nuex = 0.454 * sqrt(Rex);
	// 	double tau = VIS0 * (U[i][2] - U[i][1]) / DYV[2];
	// 	double Cf = tau * 2;
	// 	double q = K0 * (T[i][1] - T[i][2]) / DYV[2];
	// 	double h = QBCS / (T[i][1] - T0);
	// 	// double Nu = q * XP[i] / K0;
	// 	double Nu = h * XP[i] / K0;
	// 	fprintf(fp1, " %4d %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e\n", i, XP[i], V[i][NJ], Cf, Cfex, Nu, Nuex);
	// }
	// fclose(fp1);

	FILE *fp2 = fopen("out_y", "w");

	double qsum = 0;
	for (j = 2; j <= NJM; j++)
	{
		// double h = QBCE / (T[NI][j] - Tm[j]);
		// double Nu = h * 2 * XU[NI] / K0;
		
		// fprintf(fp2, " %4d %15.7e %15.7e %15.7e %15.7e %15.7e\n", j, YP[j], U[NI/2][j], U[NIM][j], T[NI/2][j], T[NIM][j]);

		qsum+=(T[1][j]-T[2][j])/DXU[2]*DYP[j];
	}
	double Nu=qsum/DTemp;
	

	SF_GET();
	double SFmax=-1.e20;
	double SFmin=1.e20;
	for(j=2; j<=NJ;j++){
		for(i=2; i<=NI;i++){
			SFmax=fmax(SFmax,SF[i][j]);
			SFmin = fmin(SFmin, SF[i][j]);
		}
	}

	printf("NU= %15.7e SFmax= %15.7e SFmin= %15.7e\n",Nu, SFmax,SFmin);
	fclose(fp2);

	FILE *fp3 = fopen("cont", "w");
	for (int n=1;n<10;n++){
		double Tout=-0.5*DTemp+DTemp/10.*n;
		Plot_Contour(fp3, Tout, 1, NI, 1, NJ, XP, YP, T);
	}
	fclose(fp3);

	FILE *fp4 = fopen("stream", "w");
	for (int n=1;n<10;n++){
		double fout=SFmin + (SFmax-SFmin)/10.*n;
		Plot_Contour(fp4, fout, 2, NI, 2, NJ, XU, YV, SF);
	}
	fclose(fp4);

	FILE *fpU = fopen("out_U", "w");
	for (int i=1;i<=NI;i++){
		for (int j=1;j<=NJ;j++){
			fprintf(fpU,"%15.7e ",U[i][j]);
		}
		fprintf(fpU,"\n");
	}
	fclose(fpU);

	FILE *fpV = fopen("out_V", "w");
	for (int i=1;i<=NI;i++){
		for (int j=1;j<=NJ;j++){
			fprintf(fpV,"%15.7e ",V[i][j]);
		}
		fprintf(fpV,"\n");
	}
	fclose(fpV);

	return 0;
}
