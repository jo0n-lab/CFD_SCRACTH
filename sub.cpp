#include <iostream>
#include "common.h" 

void T_SOLVE()
{
	for(j=1; j<=NJ;  j++){ for(i=1; i<=NI; i++){ BT[i][j]=T[i][j]; } }
//
	for(j=2; j<=NJM;  j++){
		for(i=2; i<=NIM; i++){
			double kw = TKK_INTP(TKKP[i][j], TKKP[i - 1][j]);
			double ke = TKK_INTP(TKKP[i][j], TKKP[i + 1][j]);
			double ks = TKK_INTP(TKKP[i][j], TKKP[i][j - 1]);
			double kn = TKK_INTP(TKKP[i][j], TKKP[i][j + 1]);
			AW[i][j]=kw*RU[i  ]*DYP[j]/DXU[i  ];
			AE[i][j]=ke*RU[i+1]*DYP[j]/DXU[i+1];
			AS[i][j]=ks*RP[i  ]*DXP[i]/DYV[j  ];
			AN[i][j]=kn*RP[i  ]*DXP[i]/DYV[j+1];

			double flw = ROCP[i][j] * RU[i] * DYP[j] * U[i][j];
			double fle = ROCP[i][j] * RU[i + 1] * DYP[j] * U[i + 1][j];
			double fls = ROCP[i][j] * RP[i] * DXP[i] * V[i][j];
			double fln = ROCP[i][j] * RP[i] * DXP[i] * V[i][j+1];
			AW[i][j] += fmax(0., flw);
			AE[i][j] += fmax(0., -fle);
			AS[i][j] += fmax(0., fls);
			AN[i][j] += fmax(0., -fln);

			double ac=0.;
			AP[i][j]=AW[i][j]+AE[i][j]+AS[i][j]+AN[i][j]+ac;
			TSOR[i][j]=QDOT*RP[i]*DXP[i]*DYP[j]+ac*BT[i][j];
			APC[i][j]=AP[i][j];
		}
	}

//---------------------------- X-BC
	for(j=2; j<=NJM; j++){
		if (IBCW == IWQQ) AP[2][j] -= AW[2][j];
		if (IBCE == IWQQ) AP[NIM][j] -= AE[NIM][j];

	}

//---------------------------- Y-BC
	for (i = 2; i <= NIM; i++) {
		if (IBCS == IWQQ) AP[i][2] -= AS[i][2];
		if (IBCN == IWQQ) AP[i][NJM] -= AN[i][NJM];
	}
	
	for(j=2; j<=NJM;  j++){
		for(i=2; i<=NIM; i++){
			SS[i][j]=TSOR[i][j]-APC[i][j]*T[i][j]+AW[i][j]*T[i-1][j]+AE[i][j]*T[i+1][j]
												 +AS[i][j]*T[i][j-1]+AN[i][j]*T[i][j+1];
		}
	}
//--------------------------------------------- iteration
	for(int iter=0; iter<=ITMAX; iter++){

		if(iter==iter/1*1){
			printf(" iter= %4d T= %11.3e ResMax= %11.3e Fcmax= %11.3e  \n", iter, T[10][10], SSMAX, FCMAX);
		}

		if(iter>0 && SSMAX<ERRMAX && FCMAX<ERRMAX){
			printf(" iter= %4d T= %11.3e ResMax= %11.3e Fcmax= %11.3e  \n", iter, T[10][10], SSMAX, FCMAX);
			break;
		}

		LNTDMA(2,NIM,2,NJM);

		for(j=2; j<=NJM;  j++){
			for(i=2; i<=NIM; i++){
				T[i][j]=T[i][j]+RELAX*CC[i][j];
			}
		}

//---------------------------- X-BC
		for(j=1; j<=NJ; j++){
			if (IBCW == IWQQ) T[1][j] = T[2][j]; //+ QBCW * RU[2] * DYP[j] / AW[2][j];
			if (IBCE == IWQQ) T[NI][j] = T[NIM][j]; //+ QBCE * RU[NI] * DYP[j] / AE[NIM][j];
		}

//---------------------------- Y-BC
		for(i=2; i<=NIM;  i++){
			if (IBCS == IWQQ) T[i][1] = T[i][2]; //+ QBCS * RP[i] * DXP[i] / AS[i][2];
			if (IBCN == IWQQ) T[i][NJ] = T[i][NJM]; //+ QBCN * RP[i] * DXP[i] / AS[i][NJM];
		}
//--------------------------------------- SS
		SSMAX=0.; FFMAX=0.; FCMAX=0.;
		for (j = 2; j <= NJM; j++) {
			for (i = 2; i <= NIM; i++) {
				SS[i][j] = TSOR[i][j] - APC[i][j] * T[i][j] + AW[i][j] * T[i - 1][j] + AE[i][j] * T[i + 1][j]
					+ AS[i][j] * T[i][j - 1] + AN[i][j] * T[i][j + 1];
				SSMAX = fmax(SSMAX, fabs(SS[i][j]) / APC[i][j]);
				FFMAX = fmax(FFMAX, fabs(T[i][j]));
				FCMAX = fmax(FCMAX, fabs(CC[i][j]));
			}
		}
		FFMAX=fmax(1.e-30,FFMAX); SSMAX=SSMAX/FFMAX; FCMAX=FCMAX/FFMAX;
	}

}

double TKK_INTP(double a, double b) { return 0.5 * (a + b); }

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
	
}

void LNTDMA(int ist, int iend, int jst, int jend)
{
	double beta, sss, gama[IJDIM];

	for(i=ist-1; i<=iend+1;  i++){ for(j=jst-1; j<=jend+1; j++){CC[i][j]=0.;} }
//
	for(int nt=1; nt<=2; nt++){
//-------------------------------------- x-TDMA
		for(j=jst; j<=jend; j++){
			i=ist; sss=SS[i][j]+AS[i][j]*CC[i][j-1]+AN[i][j]*CC[i][j+1];
			beta=AP[i][j]; CC[i][j]=sss/beta;
			for(i=ist+1; i<=iend; i++){
				sss=SS[i][j]+AS[i][j]*CC[i][j-1]+AN[i][j]*CC[i][j+1];
				gama[i]=-AE[i-1][j]/beta; beta=AP[i][j]+AW[i][j]*gama[i];
				CC[i][j]=(sss+AW[i][j]*CC[i-1][j])/beta;
			}
			for(i=iend-1; i>=ist; i--){ CC[i][j]=CC[i][j]-gama[i+1]*CC[i+1][j]; }
		}
//------------------------------------- y-TDMA
		for(i=ist; i<=iend; i++){
			j=jst; sss=SS[i][j]+AW[i][j]*CC[i-1][j]+AE[i][j]*CC[i+1][j];
			beta=AP[i][j]; CC[i][j]=sss/beta;
			for(j=jst+1; j<=jend; j++){
				sss=SS[i][j]+AW[i][j]*CC[i-1][j]+AE[i][j]*CC[i+1][j];
				gama[j]=-AN[i][j-1]/beta; beta=AP[i][j]+AS[i][j]*gama[j];
				CC[i][j]=(sss+AS[i][j]*CC[i][j-1])/beta;
			}
			for(j=jend-1; j>=jst; j--){ CC[i][j]=CC[i][j]-gama[j+1]*CC[i][j+1]; }
		}
	}	
}

void Plot_Contour(FILE *fp, double f0, int ibeg, int iend, int jbeg, int jend, double x[IDIM], double y[JDIM], double f[IDIM][JDIM])
{
	int i0,j0,i1,j1,n;
	double x1[4],y1[4],df00,df10,df11,df01;
	for(j0=jbeg; j0<jend; j0++){
		for(i0=ibeg; i0<iend; i0++){
			j1=j0+1; i1=i0+1; n=-1;
			df00=f0-f[i0][j0]; df10=f0-f[i1][j0];
			df11=f0-f[i1][j1]; df01=f0-f[i0][j1];

			if((df00*df10 <= 0.)&&(df00 != df10)){
				n++; x1[n]=x[i0]+df00/(df00-df10)*(x[i1]-x[i0]); y1[n]=y[j0];
			}

			if((df10*df11 <= 0.)&&(df10 != df11)){
				n++; x1[n]=x[i1]; y1[n]=y[j0]+df10/(df10-df11)*(y[j1]-y[j0]); 
			}

			if((df11*df01 <= 0.)&&(df11 != df01)){
				n++; x1[n]=x[i0]+df01/(df01-df11)*(x[i1]-x[i0]); y1[n]=y[j1]; 
			}

			if((df01*df00 <= 0.)&&(df01 != df00)){
				n++; x1[n]=x[i0]; y1[n]=y[j0]+df00/(df00-df01)*(y[j1]-y[j0]); 
			}

			if(n >= 1){
				for(int nn=0; nn<n; nn++){
					fprintf(fp,"  %9.5f %9.5f \n",    x1[nn  ], y1[nn  ]);
					fprintf(fp,"  %9.5f %9.5f \n \n", x1[nn+1], y1[nn+1]);
				}
				if(n != 1){
					fprintf(fp,"  %9.5f %9.5f \n",    x1[n], y1[n]);
					fprintf(fp,"  %9.5f %9.5f \n \n", x1[0], y1[0]);
				}
			}
		}
	}
}

