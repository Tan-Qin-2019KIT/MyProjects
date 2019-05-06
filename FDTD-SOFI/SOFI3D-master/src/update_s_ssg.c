/*------------------------------------------------------------------------
 * Copyright (C) 2011 For the list of authors, see file AUTHORS.
 *
 * This file is part of SOFI3D.
 * 
 * SOFI3D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 * 
 * SOFI3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with SOFI3D. See file COPYING and/or 
  * <http://www.gnu.org/licenses/gpl-2.0.html>.
--------------------------------------------------------------------------*/
/*------------------------------------------------------------------------
 *   updating stress at gridpoints [nx1...nx2][ny1...ny2][nz1...nz2]
 *   by a staggered grid finite difference scheme of various order accuracy in space
 *   and second order accuracy in time
 *   viscoelastic version
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

double update_s(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2,  int nt,
		float *** vx, float *** vy, float *** vz,
		float *** sxx, float *** syy, float *** szz, float *** sxy, float *** syz, float *** sxz, float *** rxx, float *** ryy,
		float *** rzz, float *** rxy, float *** ryz, float *** rxz, float ***  pi, float ***  u, float ***  uipjp, float ***  ujpkp, float ***  uipkp,
                float  ***  taus, float  ***  tausipjp, float  ***  tausjpkp, float  ***  tausipkp, float  ***  taup, float *  eta,
                float *** vxyyx,float *** vyzzy,float *** vxzzx,float *** vxxyyzz,float *** vyyzz,float *** vxxzz,float *** vxxyy,
                float *** vxyyx_2,float *** vyzzy_2,float *** vxzzx_2,float *** vxxyyzz_2,float *** vyyzz_2,float *** vxxzz_2,float *** vxxyy_2,
                float *** vxyyx_3,float *** vyzzy_3,float *** vxzzx_3,float *** vxxyyzz_3,float *** vyyzz_3,float *** vxxzz_3,float *** vxxyy_3,
                float *** vxyyx_4,float *** vyzzy_4,float *** vxzzx_4,float *** vxxyyzz_4,float *** vyyzz_4,float *** vxxzz_4,float *** vxxyy_4,float *** rxx_2, float *** ryy_2,
                float *** rzz_2, float *** rxy_2, float *** ryz_2, float *** rxz_2,float *** rxx_3, float *** ryy_3,
                float *** rzz_3, float *** rxy_3, float *** ryz_3, float *** rxz_3,float *** rxx_4, float *** ryy_4,
                float *** rzz_4, float *** rxy_4, float *** ryz_4, float *** rxz_4){

	extern float DT, DX, DY, DZ;
	extern int L, MYID, FDORDER,FDORDER_TIME, LOG, FDCOEFF;
	extern FILE *FP;
	extern int OUTNTIMESTEPINFO;

	int i, j, k, l=1;
	double time=0.0, time1=0.0, time2=0.0;
	float vxx,vxy,vxz,vyx,vyy,vyz,vzx,vzy,vzz;
	float b,c,e,g,d,f,fipjp,fjpkp,fipkp,dipjp,djpkp,dipkp;
    float vdiag;
    float n1,n2,n3,n4,n5;
    
    float vxyyx_T2,vyzzy_T2,vxzzx_T2,vxxyyzz_T2,vyyzz_T2,vxxzz_T2,vxxyy_T2; /* for FDORDER_TIME=2 */
    float c1, c2, c3, c4; /* Coefficients for Adam Bashforth */
    
    float **vxyyx_j,**vyzzy_j,**vxzzx_j,**vxxyyzz_j,**vyyzz_j,**vxxzz_j,**vxxyy_j;
    float **vxyyx_j_2,**vyzzy_j_2,**vxzzx_j_2,**vxxyyzz_j_2,**vyyzz_j_2,**vxxzz_j_2,**vxxyy_j_2;
    float **vxyyx_j_3,**vyzzy_j_3,**vxzzx_j_3,**vxxyyzz_j_3,**vyyzz_j_3,**vxxzz_j_3,**vxxyy_j_3;
    float **vxyyx_j_4,**vyzzy_j_4,**vxzzx_j_4,**vxxyyzz_j_4,**vyyzz_j_4,**vxxzz_j_4,**vxxyy_j_4;
    float *vxyyx_j_i,*vyzzy_j_i,*vxzzx_j_i,*vxxyyzz_j_i,*vyyzz_j_i,*vxxzz_j_i,*vxxyy_j_i;
    float *vxyyx_j_i_2,*vyzzy_j_i_2,*vxzzx_j_i_2,*vxxyyzz_j_i_2,*vyyzz_j_i_2,*vxxzz_j_i_2,*vxxyy_j_i_2;
    float *vxyyx_j_i_3,*vyzzy_j_i_3,*vxzzx_j_i_3,*vxxyyzz_j_i_3,*vyyzz_j_i_3,*vxxzz_j_i_3,*vxxyy_j_i_3;
    float *vxyyx_j_i_4,*vyzzy_j_i_4,*vxzzx_j_i_4,*vxxyyzz_j_i_4,*vyyzz_j_i_4,*vxxzz_j_i_4,*vxxyy_j_i_4;
    
    float **rxx_j,**ryy_j,**rzz_j,**rxy_j,**ryz_j,**rxz_j;
    float **rxx_j_2,**ryy_j_2,**rzz_j_2,**rxy_j_2,**ryz_j_2,**rxz_j_2;
    float **rxx_j_3,**ryy_j_3,**rzz_j_3,**rxy_j_3,**ryz_j_3,**rxz_j_3;
    float **rxx_j_4,**ryy_j_4,**rzz_j_4,**rxy_j_4,**ryz_j_4,**rxz_j_4;
    
    float *rxx_j_i,*ryy_j_i,*rzz_j_i,*rxy_j_i,*ryz_j_i,*rxz_j_i;
    float *rxx_j_i_2,*ryy_j_i_2,*rzz_j_i_2,*rxy_j_i_2,*ryz_j_i_2,*rxz_j_i_2;
    float *rxx_j_i_3,*ryy_j_i_3,*rzz_j_i_3,*rxy_j_i_3,*ryz_j_i_3,*rxz_j_i_3;
    float *rxx_j_i_4,*ryy_j_i_4,*rzz_j_i_4,*rxy_j_i_4,*ryz_j_i_4,*rxz_j_i_4;
    
	float sumrxy,sumryz,sumrxz,sumrxx,sumryy,sumrzz;
	float b1, b2, b3, b4, b5, b6, dthalbe;

	dthalbe = DT/2.0;

	if (LOG)
		if ((MYID==0) && ((nt+(OUTNTIMESTEPINFO-1))%OUTNTIMESTEPINFO)==0) time1=MPI_Wtime();

    
    
    switch (FDORDER_TIME) {
            
        case 2:
            switch (FDORDER){
                    
                case 2 :
                    
                    
                    
                    for (j=ny1;j<=ny2;j++){
                        for (i=nx1;i<=nx2;i++){
                            for (k=nz1;k<=nz2;k++){
                                
                                /* spatial derivatives of the components of the velocities
                                 are computed */
                                
                                vxx = (vx[j][i][k]-vx[j][i-1][k])/DX;
                                vxy = (vx[j+1][i][k]-vx[j][i][k])/DY;
                                vxz = (vx[j][i][k+1]-vx[j][i][k])/DZ;
                                vyx = (vy[j][i+1][k]-vy[j][i][k])/DX;
                                vyy = (vy[j][i][k]-vy[j-1][i][k])/DY;
                                vyz = (vy[j][i][k+1]-vy[j][i][k])/DZ;
                                vzx = (vz[j][i+1][k]-vz[j][i][k])/DX;
                                vzy = (vz[j+1][i][k]-vz[j][i][k])/DY;
                                vzz = (vz[j][i][k]-vz[j][i][k-1])/DZ;
                                
                                
                                /* computing sums of the old memory variables in this version only one mechanism is possi
                                 ble */
                                
                                sumrxy=rxy[j][i][k];
                                sumryz=ryz[j][i][k];
                                sumrxz=rxz[j][i][k];
                                sumrxx=rxx[j][i][k];
                                sumryy=ryy[j][i][k];
                                sumrzz=rzz[j][i][k];
                                
                                
                                /* updating components of the stress tensor, partially */
                                fipjp=uipjp[j][i][k]*DT*(1.0+L*tausipjp[j][i][k]);
                                fjpkp=ujpkp[j][i][k]*DT*(1.0+L*tausjpkp[j][i][k]);
                                fipkp=uipkp[j][i][k]*DT*(1.0+L*tausipkp[j][i][k]);
                                g=pi[j][i][k]*(1.0+L*taup[j][i][k]);
                                f=2.0*u[j][i][k]*(1.0+L*taus[j][i][k]);
                                
                                vxyyx_T2=vxy+vyx;
                                vyzzy_T2=vyz+vzy;
                                vxzzx_T2=vxz+vzx;
                                vxxyyzz_T2=vxx+vyy+vzz;
                                vyyzz_T2=vyy+vzz;
                                vxxzz_T2=vxx+vzz;
                                vxxyy_T2=vxx+vyy;
                                
                                sxy[j][i][k]+=(fipjp*vxyyx_T2)+(dthalbe*sumrxy);
                                syz[j][i][k]+=(fjpkp*vyzzy_T2)+(dthalbe*sumryz);
                                sxz[j][i][k]+=(fipkp*vxzzx_T2)+(dthalbe*sumrxz);
                                sxx[j][i][k]+=DT*((g*vxxyyzz_T2)-(f*vyyzz_T2))+(dthalbe*sumrxx);
                                syy[j][i][k]+=DT*((g*vxxyyzz_T2)-(f*vxxzz_T2))+(dthalbe*sumryy);
                                szz[j][i][k]+=DT*((g*vxxyyzz_T2)-(f*vxxyy_T2))+(dthalbe*sumrzz);
                                
                                /* in the case of L=1 update
                                 the memory-variables without suming*/
                                l=1;
                                b=1.0/(1.0+(eta[l]*0.5));
                                c=1.0-(eta[l]*0.5);
                                dipjp=uipjp[j][i][k]*eta[l]*tausipjp[j][i][k];
                                djpkp=ujpkp[j][i][k]*eta[l]*tausjpkp[j][i][k];
                                dipkp=uipkp[j][i][k]*eta[l]*tausipkp[j][i][k];
                                d=2.0*u[j][i][k]*eta[l]*taus[j][i][k];
                                e=pi[j][i][k]*eta[l]*taup[j][i][k];
                                rxy[j][i][k]=b*(rxy[j][i][k]*c-(dipjp*vxyyx_T2));
                                ryz[j][i][k]=b*(ryz[j][i][k]*c-(djpkp*vyzzy_T2));
                                rxz[j][i][k]=b*(rxz[j][i][k]*c-(dipkp*vxzzx_T2));
                                rxx[j][i][k]=b*(rxx[j][i][k]*c-(e*vxxyyzz_T2)+(d*vyyzz_T2));
                                ryy[j][i][k]=b*(ryy[j][i][k]*c-(e*vxxyyzz_T2)+(d*vxxzz_T2));
                                rzz[j][i][k]=b*(rzz[j][i][k]*c-(e*vxxyyzz_T2)+(d*vxxyy_T2));
                                sumrxy=rxy[j][i][k];
                                sumryz=ryz[j][i][k];
                                sumrxz=rxz[j][i][k];
                                sumrxx=rxx[j][i][k];
                                sumryy=ryy[j][i][k];
                                sumrzz=rzz[j][i][k];
                                
                                
                                
                                /* and now the components of the stress tensor are
                                 completely updated */
                                sxy[j][i][k]+=(dthalbe*sumrxy);
                                syz[j][i][k]+=(dthalbe*sumryz);
                                sxz[j][i][k]+=(dthalbe*sumrxz);
                                sxx[j][i][k]+=(dthalbe*sumrxx);
                                syy[j][i][k]+=(dthalbe*sumryy);
                                szz[j][i][k]+=(dthalbe*sumrzz);
                                
                                
                                
                                
                                
                                
                            }
                        }
                    }
                    break;
                    
                case 4 :
                    
                    b1=9.0/8.0; b2=-1.0/24.0; /* Taylor coefficients */
                    if(FDCOEFF==2){
                        b1=1.1382; b2=-0.046414;} /* Holberg coefficients E=0.1 %*/
                    
                    for (j=ny1;j<=ny2;j++){
                        for (i=nx1;i<=nx2;i++){
                            for (k=nz1;k<=nz2;k++){
                                
                                /* spatial derivatives of the components of the velocities
                                 are computed */
                                
                                vxx = (b1*(vx[j][i][k]-vx[j][i-1][k])+b2*(vx[j][i+1][k]-vx[j][i-2][k]))/DX;
                                vxy = (b1*(vx[j+1][i][k]-vx[j][i][k])+b2*(vx[j+2][i][k]-vx[j-1][i][k]))/DY;
                                vxz = (b1*(vx[j][i][k+1]-vx[j][i][k])+b2*(vx[j][i][k+2]-vx[j][i][k-1]))/DZ;
                                vyx = (b1*(vy[j][i+1][k]-vy[j][i][k])+b2*(vy[j][i+2][k]-vy[j][i-1][k]))/DX;
                                vyy = (b1*(vy[j][i][k]-vy[j-1][i][k])+b2*(vy[j+1][i][k]-vy[j-2][i][k]))/DY;
                                vyz = (b1*(vy[j][i][k+1]-vy[j][i][k])+b2*(vy[j][i][k+2]-vy[j][i][k-1]))/DZ;
                                vzx = (b1*(vz[j][i+1][k]-vz[j][i][k])+b2*(vz[j][i+2][k]-vz[j][i-1][k]))/DX;
                                vzy = (b1*(vz[j+1][i][k]-vz[j][i][k])+b2*(vz[j+2][i][k]-vz[j-1][i][k]))/DY;
                                vzz = (b1*(vz[j][i][k]-vz[j][i][k-1])+b2*(vz[j][i][k+1]-vz[j][i][k-2]))/DZ;
                                
                                
                                /* computing sums of the old memory variables in this version only one mechanism is possi
                                 ble */
                                
                                sumrxy=rxy[j][i][k];
                                sumryz=ryz[j][i][k];
                                sumrxz=rxz[j][i][k];
                                sumrxx=rxx[j][i][k];
                                sumryy=ryy[j][i][k];
                                sumrzz=rzz[j][i][k];
                                
                                
                                /* updating components of the stress tensor, partially */
                                fipjp=uipjp[j][i][k]*DT*(1.0+L*tausipjp[j][i][k]);
                                fjpkp=ujpkp[j][i][k]*DT*(1.0+L*tausjpkp[j][i][k]);
                                fipkp=uipkp[j][i][k]*DT*(1.0+L*tausipkp[j][i][k]);
                                g=pi[j][i][k]*(1.0+L*taup[j][i][k]);
                                f=2.0*u[j][i][k]*(1.0+L*taus[j][i][k]);
                                
                                
                                vxyyx_T2=vxy+vyx;
                                vyzzy_T2=vyz+vzy;
                                vxzzx_T2=vxz+vzx;
                                vxxyyzz_T2=vxx+vyy+vzz;
                                vyyzz_T2=vyy+vzz;
                                vxxzz_T2=vxx+vzz;
                                vxxyy_T2=vxx+vyy;
                                
                                sxy[j][i][k]+=(fipjp*vxyyx_T2)+(dthalbe*sumrxy);
                                syz[j][i][k]+=(fjpkp*vyzzy_T2)+(dthalbe*sumryz);
                                sxz[j][i][k]+=(fipkp*vxzzx_T2)+(dthalbe*sumrxz);
                                sxx[j][i][k]+=DT*((g*vxxyyzz_T2)-(f*vyyzz_T2))+(dthalbe*sumrxx);
                                syy[j][i][k]+=DT*((g*vxxyyzz_T2)-(f*vxxzz_T2))+(dthalbe*sumryy);
                                szz[j][i][k]+=DT*((g*vxxyyzz_T2)-(f*vxxyy_T2))+(dthalbe*sumrzz);
                                
                                /* in the case of L=1 update
                                 the memory-variables without suming*/
                                l=1;
                                b=1.0/(1.0+(eta[l]*0.5));
                                c=1.0-(eta[l]*0.5);
                                dipjp=uipjp[j][i][k]*eta[l]*tausipjp[j][i][k];
                                djpkp=ujpkp[j][i][k]*eta[l]*tausjpkp[j][i][k];
                                dipkp=uipkp[j][i][k]*eta[l]*tausipkp[j][i][k];
                                d=2.0*u[j][i][k]*eta[l]*taus[j][i][k];
                                e=pi[j][i][k]*eta[l]*taup[j][i][k];
                                rxy[j][i][k]=b*(rxy[j][i][k]*c-(dipjp*vxyyx_T2));
                                ryz[j][i][k]=b*(ryz[j][i][k]*c-(djpkp*vyzzy_T2));
                                rxz[j][i][k]=b*(rxz[j][i][k]*c-(dipkp*vxzzx_T2));
                                rxx[j][i][k]=b*(rxx[j][i][k]*c-(e*vxxyyzz_T2)+(d*vyyzz_T2));
                                ryy[j][i][k]=b*(ryy[j][i][k]*c-(e*vxxyyzz_T2)+(d*vxxzz_T2));
                                rzz[j][i][k]=b*(rzz[j][i][k]*c-(e*vxxyyzz_T2)+(d*vxxyy_T2));
                                sumrxy=rxy[j][i][k];
                                sumryz=ryz[j][i][k];
                                sumrxz=rxz[j][i][k];
                                sumrxx=rxx[j][i][k];
                                sumryy=ryy[j][i][k];
                                sumrzz=rzz[j][i][k];
                                
                                
                                
                                /* and now the components of the stress tensor are
                                 completely updated */
                                sxy[j][i][k]+=(dthalbe*sumrxy);
                                syz[j][i][k]+=(dthalbe*sumryz);
                                sxz[j][i][k]+=(dthalbe*sumrxz);
                                sxx[j][i][k]+=(dthalbe*sumrxx);
                                syy[j][i][k]+=(dthalbe*sumryy);
                                szz[j][i][k]+=(dthalbe*sumrzz);
                                
                                
                                
                                
                            }
                        }
                    }
                    break;
                    
                case 6 :
                    
                    b1=75.0/64.0; b2=-25.0/384.0; b3=3.0/640.0; /* Taylor coefficients */
                    if(FDCOEFF==2){
                        b1=1.1965; b2=-0.078804; b3=0.0081781;}   /* Holberg coefficients E=0.1 %*/
                    
                    for (j=ny1;j<=ny2;j++){
                        for (i=nx1;i<=nx2;i++){
                            for (k=nz1;k<=nz2;k++){
                                
                                /* spatial derivatives of the components of the velocities
                                 are computed */
                                
                                vxx = (b1*(vx[j][i][k]-vx[j][i-1][k])+
                                       b2*(vx[j][i+1][k]-vx[j][i-2][k])+
                                       b3*(vx[j][i+2][k]-vx[j][i-3][k]))/DX;
                                
                                vxy = (b1*(vx[j+1][i][k]-vx[j][i][k])+
                                       b2*(vx[j+2][i][k]-vx[j-1][i][k])+
                                       b3*(vx[j+3][i][k]-vx[j-2][i][k]))/DY;
                                
                                vxz = (b1*(vx[j][i][k+1]-vx[j][i][k])+
                                       b2*(vx[j][i][k+2]-vx[j][i][k-1])+
                                       b3*(vx[j][i][k+3]-vx[j][i][k-2]))/DZ;
                                
                                vyx = (b1*(vy[j][i+1][k]-vy[j][i][k])+
                                       b2*(vy[j][i+2][k]-vy[j][i-1][k])+
                                       b3*(vy[j][i+3][k]-vy[j][i-2][k]))/DX;
                                
                                vyy = (b1*(vy[j][i][k]-vy[j-1][i][k])+
                                       b2*(vy[j+1][i][k]-vy[j-2][i][k])+
                                       b3*(vy[j+2][i][k]-vy[j-3][i][k]))/DY;
                                
                                vyz = (b1*(vy[j][i][k+1]-vy[j][i][k])+
                                       b2*(vy[j][i][k+2]-vy[j][i][k-1])+
                                       b3*(vy[j][i][k+3]-vy[j][i][k-2]))/DZ;
                                
                                vzx = (b1*(vz[j][i+1][k]-vz[j][i][k])+
                                       b2*(vz[j][i+2][k]-vz[j][i-1][k])+
                                       b3*(vz[j][i+3][k]-vz[j][i-2][k]))/DX;
                                
                                vzy = (b1*(vz[j+1][i][k]-vz[j][i][k])+
                                       b2*(vz[j+2][i][k]-vz[j-1][i][k])+
                                       b3*(vz[j+3][i][k]-vz[j-2][i][k]))/DY;
                                
                                vzz = (b1*(vz[j][i][k]-vz[j][i][k-1])+
                                       b2*(vz[j][i][k+1]-vz[j][i][k-2])+
                                       b3*(vz[j][i][k+2]-vz[j][i][k-3]))/DZ;
                                
                                
                                /* computing sums of the old memory variables in this version only one mechanism is possi
                                 ble */
                                
                                sumrxy=rxy[j][i][k];
                                sumryz=ryz[j][i][k];
                                sumrxz=rxz[j][i][k];
                                sumrxx=rxx[j][i][k];
                                sumryy=ryy[j][i][k];
                                sumrzz=rzz[j][i][k];
                                
                                
                                /* updating components of the stress tensor, partially */
                                fipjp=uipjp[j][i][k]*DT*(1.0+L*tausipjp[j][i][k]);
                                fjpkp=ujpkp[j][i][k]*DT*(1.0+L*tausjpkp[j][i][k]);
                                fipkp=uipkp[j][i][k]*DT*(1.0+L*tausipkp[j][i][k]);
                                g=pi[j][i][k]*(1.0+L*taup[j][i][k]);
                                f=2.0*u[j][i][k]*(1.0+L*taus[j][i][k]);
                                
                                
                                vxyyx_T2=vxy+vyx;
                                vyzzy_T2=vyz+vzy;
                                vxzzx_T2=vxz+vzx;
                                vxxyyzz_T2=vxx+vyy+vzz;
                                vyyzz_T2=vyy+vzz;
                                vxxzz_T2=vxx+vzz;
                                vxxyy_T2=vxx+vyy;
                                
                                sxy[j][i][k]+=(fipjp*vxyyx_T2)+(dthalbe*sumrxy);
                                syz[j][i][k]+=(fjpkp*vyzzy_T2)+(dthalbe*sumryz);
                                sxz[j][i][k]+=(fipkp*vxzzx_T2)+(dthalbe*sumrxz);
                                sxx[j][i][k]+=DT*((g*vxxyyzz_T2)-(f*vyyzz_T2))+(dthalbe*sumrxx);
                                syy[j][i][k]+=DT*((g*vxxyyzz_T2)-(f*vxxzz_T2))+(dthalbe*sumryy);
                                szz[j][i][k]+=DT*((g*vxxyyzz_T2)-(f*vxxyy_T2))+(dthalbe*sumrzz);
                                
                                /* in the case of L=1 update
                                 the memory-variables without suming*/
                                l=1;
                                b=1.0/(1.0+(eta[l]*0.5));
                                c=1.0-(eta[l]*0.5);
                                dipjp=uipjp[j][i][k]*eta[l]*tausipjp[j][i][k];
                                djpkp=ujpkp[j][i][k]*eta[l]*tausjpkp[j][i][k];
                                dipkp=uipkp[j][i][k]*eta[l]*tausipkp[j][i][k];
                                d=2.0*u[j][i][k]*eta[l]*taus[j][i][k];
                                e=pi[j][i][k]*eta[l]*taup[j][i][k];
                                rxy[j][i][k]=b*(rxy[j][i][k]*c-(dipjp*vxyyx_T2));
                                ryz[j][i][k]=b*(ryz[j][i][k]*c-(djpkp*vyzzy_T2));
                                rxz[j][i][k]=b*(rxz[j][i][k]*c-(dipkp*vxzzx_T2));
                                rxx[j][i][k]=b*(rxx[j][i][k]*c-(e*vxxyyzz_T2)+(d*vyyzz_T2));
                                ryy[j][i][k]=b*(ryy[j][i][k]*c-(e*vxxyyzz_T2)+(d*vxxzz_T2));
                                rzz[j][i][k]=b*(rzz[j][i][k]*c-(e*vxxyyzz_T2)+(d*vxxyy_T2));
                                sumrxy=rxy[j][i][k];
                                sumryz=ryz[j][i][k];
                                sumrxz=rxz[j][i][k];
                                sumrxx=rxx[j][i][k];
                                sumryy=ryy[j][i][k];
                                sumrzz=rzz[j][i][k];
                                
                                
                                
                                /* and now the components of the stress tensor are
                                 completely updated */
                                sxy[j][i][k]+=(dthalbe*sumrxy);
                                syz[j][i][k]+=(dthalbe*sumryz);
                                sxz[j][i][k]+=(dthalbe*sumrxz);
                                sxx[j][i][k]+=(dthalbe*sumrxx);
                                syy[j][i][k]+=(dthalbe*sumryy);
                                szz[j][i][k]+=(dthalbe*sumrzz);
                                
                                
                                
                                
                                
                            }
                        }
                    }
                    break;
                    
                case 8 :
                    
                    b1=1225.0/1024.0; b2=-245.0/3072.0; b3=49.0/5120.0; b4=-5.0/7168.0; /* Taylor coefficients */
                    if(FDCOEFF==2){
                        b1=1.2257; b2=-0.099537; b3=0.018063; b4=-0.0026274;} /* Holberg coefficients E=0.1 %*/
                    
                    for (j=ny1;j<=ny2;j++){
                        for (i=nx1;i<=nx2;i++){
                            for (k=nz1;k<=nz2;k++){
                                
                                /* spatial derivatives of the components of the velocities
                                 are computed */
                                
                                vxx = (b1*(vx[j][i][k]-vx[j][i-1][k])+
                                       b2*(vx[j][i+1][k]-vx[j][i-2][k])+
                                       b3*(vx[j][i+2][k]-vx[j][i-3][k])+
                                       b4*(vx[j][i+3][k]-vx[j][i-4][k]))/DX;
                                
                                vxy = (b1*(vx[j+1][i][k]-vx[j][i][k])+
                                       b2*(vx[j+2][i][k]-vx[j-1][i][k])+
                                       b3*(vx[j+3][i][k]-vx[j-2][i][k])+
                                       b4*(vx[j+4][i][k]-vx[j-3][i][k]))/DY;
                                
                                vxz = (b1*(vx[j][i][k+1]-vx[j][i][k])+
                                       b2*(vx[j][i][k+2]-vx[j][i][k-1])+
                                       b3*(vx[j][i][k+3]-vx[j][i][k-2])+
                                       b4*(vx[j][i][k+4]-vx[j][i][k-3]))/DZ;
                                
                                vyx = (b1*(vy[j][i+1][k]-vy[j][i][k])+
                                       b2*(vy[j][i+2][k]-vy[j][i-1][k])+
                                       b3*(vy[j][i+3][k]-vy[j][i-2][k])+
                                       b4*(vy[j][i+4][k]-vy[j][i-3][k]))/DX;
                                
                                vyy = (b1*(vy[j][i][k]-vy[j-1][i][k])+
                                       b2*(vy[j+1][i][k]-vy[j-2][i][k])+
                                       b3*(vy[j+2][i][k]-vy[j-3][i][k])+
                                       b4*(vy[j+3][i][k]-vy[j-4][i][k]))/DY;
                                
                                vyz = (b1*(vy[j][i][k+1]-vy[j][i][k])+
                                       b2*(vy[j][i][k+2]-vy[j][i][k-1])+
                                       b3*(vy[j][i][k+3]-vy[j][i][k-2])+
                                       b4*(vy[j][i][k+4]-vy[j][i][k-3]))/DZ;
                                
                                vzx = (b1*(vz[j][i+1][k]-vz[j][i][k])+
                                       b2*(vz[j][i+2][k]-vz[j][i-1][k])+
                                       b3*(vz[j][i+3][k]-vz[j][i-2][k])+
                                       b4*(vz[j][i+4][k]-vz[j][i-3][k]))/DX;
                                
                                vzy = (b1*(vz[j+1][i][k]-vz[j][i][k])+
                                       b2*(vz[j+2][i][k]-vz[j-1][i][k])+
                                       b3*(vz[j+3][i][k]-vz[j-2][i][k])+
                                       b4*(vz[j+4][i][k]-vz[j-3][i][k]))/DY;
                                
                                vzz = (b1*(vz[j][i][k]-vz[j][i][k-1])+
                                       b2*(vz[j][i][k+1]-vz[j][i][k-2])+
                                       b3*(vz[j][i][k+2]-vz[j][i][k-3])+
                                       b4*(vz[j][i][k+3]-vz[j][i][k-4]))/DZ;
                                
                                
                                /* computing sums of the old memory variables in this version only one mechanism is possi
                                 ble */
                                
                                sumrxy=rxy[j][i][k];
                                sumryz=ryz[j][i][k];
                                sumrxz=rxz[j][i][k];
                                sumrxx=rxx[j][i][k];
                                sumryy=ryy[j][i][k];
                                sumrzz=rzz[j][i][k];
                                
                                
                                /* updating components of the stress tensor, partially */
                                fipjp=uipjp[j][i][k]*DT*(1.0+L*tausipjp[j][i][k]);
                                fjpkp=ujpkp[j][i][k]*DT*(1.0+L*tausjpkp[j][i][k]);
                                fipkp=uipkp[j][i][k]*DT*(1.0+L*tausipkp[j][i][k]);
                                g=pi[j][i][k]*(1.0+L*taup[j][i][k]);
                                f=2.0*u[j][i][k]*(1.0+L*taus[j][i][k]);
                                
                                
                                vxyyx_T2=vxy+vyx;
                                vyzzy_T2=vyz+vzy;
                                vxzzx_T2=vxz+vzx;
                                vxxyyzz_T2=vxx+vyy+vzz;
                                vyyzz_T2=vyy+vzz;
                                vxxzz_T2=vxx+vzz;
                                vxxyy_T2=vxx+vyy;
                                
                                sxy[j][i][k]+=(fipjp*vxyyx_T2)+(dthalbe*sumrxy);
                                syz[j][i][k]+=(fjpkp*vyzzy_T2)+(dthalbe*sumryz);
                                sxz[j][i][k]+=(fipkp*vxzzx_T2)+(dthalbe*sumrxz);
                                sxx[j][i][k]+=DT*((g*vxxyyzz_T2)-(f*vyyzz_T2))+(dthalbe*sumrxx);
                                syy[j][i][k]+=DT*((g*vxxyyzz_T2)-(f*vxxzz_T2))+(dthalbe*sumryy);
                                szz[j][i][k]+=DT*((g*vxxyyzz_T2)-(f*vxxyy_T2))+(dthalbe*sumrzz);
                                
                                /* in the case of L=1 update
                                 the memory-variables without suming*/
                                l=1;
                                b=1.0/(1.0+(eta[l]*0.5));
                                c=1.0-(eta[l]*0.5);
                                dipjp=uipjp[j][i][k]*eta[l]*tausipjp[j][i][k];
                                djpkp=ujpkp[j][i][k]*eta[l]*tausjpkp[j][i][k];
                                dipkp=uipkp[j][i][k]*eta[l]*tausipkp[j][i][k];
                                d=2.0*u[j][i][k]*eta[l]*taus[j][i][k];
                                e=pi[j][i][k]*eta[l]*taup[j][i][k];
                                rxy[j][i][k]=b*(rxy[j][i][k]*c-(dipjp*vxyyx_T2));
                                ryz[j][i][k]=b*(ryz[j][i][k]*c-(djpkp*vyzzy_T2));
                                rxz[j][i][k]=b*(rxz[j][i][k]*c-(dipkp*vxzzx_T2));
                                rxx[j][i][k]=b*(rxx[j][i][k]*c-(e*vxxyyzz_T2)+(d*vyyzz_T2));
                                ryy[j][i][k]=b*(ryy[j][i][k]*c-(e*vxxyyzz_T2)+(d*vxxzz_T2));
                                rzz[j][i][k]=b*(rzz[j][i][k]*c-(e*vxxyyzz_T2)+(d*vxxyy_T2));
                                sumrxy=rxy[j][i][k];
                                sumryz=ryz[j][i][k];
                                sumrxz=rxz[j][i][k];
                                sumrxx=rxx[j][i][k];
                                sumryy=ryy[j][i][k];
                                sumrzz=rzz[j][i][k];
                                
                                
                                
                                /* and now the components of the stress tensor are
                                 completely updated */
                                sxy[j][i][k]+=(dthalbe*sumrxy);
                                syz[j][i][k]+=(dthalbe*sumryz);
                                sxz[j][i][k]+=(dthalbe*sumrxz);
                                sxx[j][i][k]+=(dthalbe*sumrxx);
                                syy[j][i][k]+=(dthalbe*sumryy);
                                szz[j][i][k]+=(dthalbe*sumrzz);
                                
                                
                                
                                
                                
                                
                            }
                        }
                    }
                    break;
                    
                case 10 :
                    
                    b1=19845.0/16384.0; b2=-735.0/8192.0; b3=567.0/40960.0; b4=-405.0/229376.0; b5=35.0/294912.0; /* Taylor coefficients */
                    if(FDCOEFF==2){
                        b1=1.2415; b2=-0.11231; b3=0.026191; b4=-0.0064682; b5=0.001191;} /* Holberg coefficients E=0.1 %*/
                    
                    for (j=ny1;j<=ny2;j++){
                        for (i=nx1;i<=nx2;i++){
                            for (k=nz1;k<=nz2;k++){
                                
                                /* spatial derivatives of the components of the velocities
                                 are computed */
                                
                                vxx = (b1*(vx[j][i][k]-vx[j][i-1][k])+
                                       b2*(vx[j][i+1][k]-vx[j][i-2][k])+
                                       b3*(vx[j][i+2][k]-vx[j][i-3][k])+
                                       b4*(vx[j][i+3][k]-vx[j][i-4][k])+
                                       b5*(vx[j][i+4][k]-vx[j][i-5][k]))/DX;
                                
                                vxy = (b1*(vx[j+1][i][k]-vx[j][i][k])+
                                       b2*(vx[j+2][i][k]-vx[j-1][i][k])+
                                       b3*(vx[j+3][i][k]-vx[j-2][i][k])+
                                       b4*(vx[j+4][i][k]-vx[j-3][i][k])+
                                       b5*(vx[j+5][i][k]-vx[j-4][i][k]))/DY;
                                
                                vxz = (b1*(vx[j][i][k+1]-vx[j][i][k])+
                                       b2*(vx[j][i][k+2]-vx[j][i][k-1])+
                                       b3*(vx[j][i][k+3]-vx[j][i][k-2])+
                                       b4*(vx[j][i][k+4]-vx[j][i][k-3])+
                                       b5*(vx[j][i][k+5]-vx[j][i][k-4]))/DZ;
                                
                                vyx = (b1*(vy[j][i+1][k]-vy[j][i][k])+
                                       b2*(vy[j][i+2][k]-vy[j][i-1][k])+
                                       b3*(vy[j][i+3][k]-vy[j][i-2][k])+
                                       b4*(vy[j][i+4][k]-vy[j][i-3][k])+
                                       b5*(vy[j][i+5][k]-vy[j][i-4][k]))/DX;
                                
                                vyy = (b1*(vy[j][i][k]-vy[j-1][i][k])+
                                       b2*(vy[j+1][i][k]-vy[j-2][i][k])+
                                       b3*(vy[j+2][i][k]-vy[j-3][i][k])+
                                       b4*(vy[j+3][i][k]-vy[j-4][i][k])+
                                       b5*(vy[j+4][i][k]-vy[j-5][i][k]))/DY;
                                
                                vyz = (b1*(vy[j][i][k+1]-vy[j][i][k])+
                                       b2*(vy[j][i][k+2]-vy[j][i][k-1])+
                                       b3*(vy[j][i][k+3]-vy[j][i][k-2])+
                                       b4*(vy[j][i][k+4]-vy[j][i][k-3])+
                                       b5*(vy[j][i][k+5]-vy[j][i][k-4]))/DZ;
                                
                                vzx = (b1*(vz[j][i+1][k]-vz[j][i][k])+
                                       b2*(vz[j][i+2][k]-vz[j][i-1][k])+
                                       b3*(vz[j][i+3][k]-vz[j][i-2][k])+
                                       b4*(vz[j][i+4][k]-vz[j][i-3][k])+
                                       b5*(vz[j][i+5][k]-vz[j][i-4][k]))/DX;
                                
                                vzy = (b1*(vz[j+1][i][k]-vz[j][i][k])+
                                       b2*(vz[j+2][i][k]-vz[j-1][i][k])+
                                       b3*(vz[j+3][i][k]-vz[j-2][i][k])+
                                       b4*(vz[j+4][i][k]-vz[j-3][i][k])+
                                       b5*(vz[j+5][i][k]-vz[j-4][i][k]))/DY;
                                
                                vzz = (b1*(vz[j][i][k]-vz[j][i][k-1])+
                                       b2*(vz[j][i][k+1]-vz[j][i][k-2])+
                                       b3*(vz[j][i][k+2]-vz[j][i][k-3])+
                                       b4*(vz[j][i][k+3]-vz[j][i][k-4])+
                                       b5*(vz[j][i][k+4]-vz[j][i][k-5]))/DZ;
                                
                                
                                /* computing sums of the old memory variables in this version only one mechanism is possi
                                 ble */
                                
                                sumrxy=rxy[j][i][k];
                                sumryz=ryz[j][i][k];
                                sumrxz=rxz[j][i][k];
                                sumrxx=rxx[j][i][k];
                                sumryy=ryy[j][i][k];
                                sumrzz=rzz[j][i][k];
                                
                                
                                /* updating components of the stress tensor, partially */
                                fipjp=uipjp[j][i][k]*DT*(1.0+L*tausipjp[j][i][k]);
                                fjpkp=ujpkp[j][i][k]*DT*(1.0+L*tausjpkp[j][i][k]);
                                fipkp=uipkp[j][i][k]*DT*(1.0+L*tausipkp[j][i][k]);
                                g=pi[j][i][k]*(1.0+L*taup[j][i][k]);
                                f=2.0*u[j][i][k]*(1.0+L*taus[j][i][k]);
                                
                                
                                vxyyx_T2=vxy+vyx;
                                vyzzy_T2=vyz+vzy;
                                vxzzx_T2=vxz+vzx;
                                vxxyyzz_T2=vxx+vyy+vzz;
                                vyyzz_T2=vyy+vzz;
                                vxxzz_T2=vxx+vzz;
                                vxxyy_T2=vxx+vyy;
                                
                                sxy[j][i][k]+=(fipjp*vxyyx_T2)+(dthalbe*sumrxy);
                                syz[j][i][k]+=(fjpkp*vyzzy_T2)+(dthalbe*sumryz);
                                sxz[j][i][k]+=(fipkp*vxzzx_T2)+(dthalbe*sumrxz);
                                sxx[j][i][k]+=DT*((g*vxxyyzz_T2)-(f*vyyzz_T2))+(dthalbe*sumrxx);
                                syy[j][i][k]+=DT*((g*vxxyyzz_T2)-(f*vxxzz_T2))+(dthalbe*sumryy);
                                szz[j][i][k]+=DT*((g*vxxyyzz_T2)-(f*vxxyy_T2))+(dthalbe*sumrzz);
                                
                                /* in the case of L=1 update
                                 the memory-variables without suming*/
                                l=1;
                                b=1.0/(1.0+(eta[l]*0.5));
                                c=1.0-(eta[l]*0.5);
                                dipjp=uipjp[j][i][k]*eta[l]*tausipjp[j][i][k];
                                djpkp=ujpkp[j][i][k]*eta[l]*tausjpkp[j][i][k];
                                dipkp=uipkp[j][i][k]*eta[l]*tausipkp[j][i][k];
                                d=2.0*u[j][i][k]*eta[l]*taus[j][i][k];
                                e=pi[j][i][k]*eta[l]*taup[j][i][k];
                                rxy[j][i][k]=b*(rxy[j][i][k]*c-(dipjp*vxyyx_T2));
                                ryz[j][i][k]=b*(ryz[j][i][k]*c-(djpkp*vyzzy_T2));
                                rxz[j][i][k]=b*(rxz[j][i][k]*c-(dipkp*vxzzx_T2));
                                rxx[j][i][k]=b*(rxx[j][i][k]*c-(e*vxxyyzz_T2)+(d*vyyzz_T2));
                                ryy[j][i][k]=b*(ryy[j][i][k]*c-(e*vxxyyzz_T2)+(d*vxxzz_T2));
                                rzz[j][i][k]=b*(rzz[j][i][k]*c-(e*vxxyyzz_T2)+(d*vxxyy_T2));
                                sumrxy=rxy[j][i][k];
                                sumryz=ryz[j][i][k];
                                sumrxz=rxz[j][i][k];
                                sumrxx=rxx[j][i][k];
                                sumryy=ryy[j][i][k];
                                sumrzz=rzz[j][i][k];
                                
                                
                                
                                /* and now the components of the stress tensor are
                                 completely updated */
                                sxy[j][i][k]+=(dthalbe*sumrxy);
                                syz[j][i][k]+=(dthalbe*sumryz);
                                sxz[j][i][k]+=(dthalbe*sumrxz);
                                sxx[j][i][k]+=(dthalbe*sumrxx);
                                syy[j][i][k]+=(dthalbe*sumryy);
                                szz[j][i][k]+=(dthalbe*sumrzz);
                                
                                
                                
                                
                                
                                
                            }
                        }
                    }
                    break;
                    
                case 12 :
                    
                    /* Taylor coefficients */
                    b1=160083.0/131072.0; b2=-12705.0/131072.0; b3=22869.0/1310720.0;
                    b4=-5445.0/1835008.0; b5=847.0/2359296.0; b6=-63.0/2883584;
                    
                    /* Holberg coefficients E=0.1 %*/
                    if(FDCOEFF==2){
                        b1=1.2508; b2=-0.12034; b3=0.032131; b4=-0.010142; b5=0.0029857; b6=-0.00066667;}
                    
                    for (j=ny1;j<=ny2;j++){
                        for (i=nx1;i<=nx2;i++){
                            for (k=nz1;k<=nz2;k++){
                                
                                /* spatial derivatives of the components of the velocities
                                 are computed */
                                
                                vxx = (b1*(vx[j][i][k]-vx[j][i-1][k])+
                                       b2*(vx[j][i+1][k]-vx[j][i-2][k])+
                                       b3*(vx[j][i+2][k]-vx[j][i-3][k])+
                                       b4*(vx[j][i+3][k]-vx[j][i-4][k])+
                                       b5*(vx[j][i+4][k]-vx[j][i-5][k])+
                                       b6*(vx[j][i+5][k]-vx[j][i-6][k]))/DX;
                                
                                vxy = (b1*(vx[j+1][i][k]-vx[j][i][k])+
                                       b2*(vx[j+2][i][k]-vx[j-1][i][k])+
                                       b3*(vx[j+3][i][k]-vx[j-2][i][k])+
                                       b4*(vx[j+4][i][k]-vx[j-3][i][k])+
                                       b5*(vx[j+5][i][k]-vx[j-4][i][k])+
                                       b6*(vx[j+6][i][k]-vx[j-5][i][k]))/DY;
                                
                                vxz = (b1*(vx[j][i][k+1]-vx[j][i][k])+
                                       b2*(vx[j][i][k+2]-vx[j][i][k-1])+
                                       b3*(vx[j][i][k+3]-vx[j][i][k-2])+
                                       b4*(vx[j][i][k+4]-vx[j][i][k-3])+
                                       b5*(vx[j][i][k+5]-vx[j][i][k-4])+
                                       b6*(vx[j][i][k+6]-vx[j][i][k-5]))/DZ;
                                
                                vyx = (b1*(vy[j][i+1][k]-vy[j][i][k])+
                                       b2*(vy[j][i+2][k]-vy[j][i-1][k])+
                                       b3*(vy[j][i+3][k]-vy[j][i-2][k])+
                                       b4*(vy[j][i+4][k]-vy[j][i-3][k])+
                                       b5*(vy[j][i+5][k]-vy[j][i-4][k])+
                                       b6*(vy[j][i+6][k]-vy[j][i-5][k]))/DX;
                                
                                vyy = (b1*(vy[j][i][k]-vy[j-1][i][k])+
                                       b2*(vy[j+1][i][k]-vy[j-2][i][k])+
                                       b3*(vy[j+2][i][k]-vy[j-3][i][k])+
                                       b4*(vy[j+3][i][k]-vy[j-4][i][k])+
                                       b5*(vy[j+4][i][k]-vy[j-5][i][k])+
                                       b6*(vy[j+5][i][k]-vy[j-6][i][k]))/DY;
                                
                                vyz = (b1*(vy[j][i][k+1]-vy[j][i][k])+
                                       b2*(vy[j][i][k+2]-vy[j][i][k-1])+
                                       b3*(vy[j][i][k+3]-vy[j][i][k-2])+
                                       b4*(vy[j][i][k+4]-vy[j][i][k-3])+
                                       b5*(vy[j][i][k+5]-vy[j][i][k-4])+
                                       b6*(vy[j][i][k+6]-vy[j][i][k-5]))/DZ;
                                
                                vzx = (b1*(vz[j][i+1][k]-vz[j][i][k])+
                                       b2*(vz[j][i+2][k]-vz[j][i-1][k])+
                                       b3*(vz[j][i+3][k]-vz[j][i-2][k])+
                                       b4*(vz[j][i+4][k]-vz[j][i-3][k])+
                                       b5*(vz[j][i+5][k]-vz[j][i-4][k])+
                                       b6*(vz[j][i+6][k]-vz[j][i-5][k]))/DX;
                                
                                vzy = (b1*(vz[j+1][i][k]-vz[j][i][k])+
                                       b2*(vz[j+2][i][k]-vz[j-1][i][k])+
                                       b3*(vz[j+3][i][k]-vz[j-2][i][k])+
                                       b4*(vz[j+4][i][k]-vz[j-3][i][k])+
                                       b5*(vz[j+5][i][k]-vz[j-4][i][k])+
                                       b6*(vz[j+6][i][k]-vz[j-5][i][k]))/DY;
                                
                                vzz = (b1*(vz[j][i][k]-vz[j][i][k-1])+
                                       b2*(vz[j][i][k+1]-vz[j][i][k-2])+
                                       b3*(vz[j][i][k+2]-vz[j][i][k-3])+
                                       b4*(vz[j][i][k+3]-vz[j][i][k-4])+
                                       b5*(vz[j][i][k+4]-vz[j][i][k-5])+
                                       b6*(vz[j][i][k+5]-vz[j][i][k-6]))/DZ;
                                
                                
                                /* computing sums of the old memory variables in this version only one mechanism is possi
                                 ble */
                                
                                sumrxy=rxy[j][i][k];
                                sumryz=ryz[j][i][k];
                                sumrxz=rxz[j][i][k];
                                sumrxx=rxx[j][i][k];
                                sumryy=ryy[j][i][k];
                                sumrzz=rzz[j][i][k];
                                
                                
                                /* updating components of the stress tensor, partially */
                                fipjp=uipjp[j][i][k]*DT*(1.0+L*tausipjp[j][i][k]);
                                fjpkp=ujpkp[j][i][k]*DT*(1.0+L*tausjpkp[j][i][k]);
                                fipkp=uipkp[j][i][k]*DT*(1.0+L*tausipkp[j][i][k]);
                                g=pi[j][i][k]*(1.0+L*taup[j][i][k]);
                                f=2.0*u[j][i][k]*(1.0+L*taus[j][i][k]);
                                
                                
                                vxyyx_T2=vxy+vyx;
                                vyzzy_T2=vyz+vzy;
                                vxzzx_T2=vxz+vzx;
                                vxxyyzz_T2=vxx+vyy+vzz;
                                vyyzz_T2=vyy+vzz;
                                vxxzz_T2=vxx+vzz;
                                vxxyy_T2=vxx+vyy;
                                
                                sxy[j][i][k]+=(fipjp*vxyyx_T2)+(dthalbe*sumrxy);
                                syz[j][i][k]+=(fjpkp*vyzzy_T2)+(dthalbe*sumryz);
                                sxz[j][i][k]+=(fipkp*vxzzx_T2)+(dthalbe*sumrxz);
                                sxx[j][i][k]+=DT*((g*vxxyyzz_T2)-(f*vyyzz_T2))+(dthalbe*sumrxx);
                                syy[j][i][k]+=DT*((g*vxxyyzz_T2)-(f*vxxzz_T2))+(dthalbe*sumryy);
                                szz[j][i][k]+=DT*((g*vxxyyzz_T2)-(f*vxxyy_T2))+(dthalbe*sumrzz);
                                
                                /* in the case of L=1 update
                                 the memory-variables without suming*/
                                l=1;
                                b=1.0/(1.0+(eta[l]*0.5));
                                c=1.0-(eta[l]*0.5);
                                dipjp=uipjp[j][i][k]*eta[l]*tausipjp[j][i][k];
                                djpkp=ujpkp[j][i][k]*eta[l]*tausjpkp[j][i][k];
                                dipkp=uipkp[j][i][k]*eta[l]*tausipkp[j][i][k];
                                d=2.0*u[j][i][k]*eta[l]*taus[j][i][k];
                                e=pi[j][i][k]*eta[l]*taup[j][i][k];
                                rxy[j][i][k]=b*(rxy[j][i][k]*c-(dipjp*vxyyx_T2));
                                ryz[j][i][k]=b*(ryz[j][i][k]*c-(djpkp*vyzzy_T2));
                                rxz[j][i][k]=b*(rxz[j][i][k]*c-(dipkp*vxzzx_T2));
                                rxx[j][i][k]=b*(rxx[j][i][k]*c-(e*vxxyyzz_T2)+(d*vyyzz_T2));
                                ryy[j][i][k]=b*(ryy[j][i][k]*c-(e*vxxyyzz_T2)+(d*vxxzz_T2));
                                rzz[j][i][k]=b*(rzz[j][i][k]*c-(e*vxxyyzz_T2)+(d*vxxyy_T2));
                                sumrxy=rxy[j][i][k];
                                sumryz=ryz[j][i][k];
                                sumrxz=rxz[j][i][k];
                                sumrxx=rxx[j][i][k];
                                sumryy=ryy[j][i][k];
                                sumrzz=rzz[j][i][k];
                                
                                
                                
                                /* and now the components of the stress tensor are
                                 completely updated */
                                sxy[j][i][k]+=(dthalbe*sumrxy);
                                syz[j][i][k]+=(dthalbe*sumryz);
                                sxz[j][i][k]+=(dthalbe*sumrxz);
                                sxx[j][i][k]+=(dthalbe*sumrxx);
                                syy[j][i][k]+=(dthalbe*sumryy);
                                szz[j][i][k]+=(dthalbe*sumrzz);
                                
                                
                                
                                
                                
                            }
                        }
                    }
                    break;
                    
            }
            break; /* Break FDORDER_TIME=2 */
            
        case 3:
            c1=25.0/24.0; c2=-1.0/12.0; c3=1.0/24.0;
            switch (FDORDER){
                    
                case 2 :
                    
                    for (j=ny1;j<=ny2;j++){
                        
                        rxx_j=*(rxx+j);ryy_j=*(ryy+j);rzz_j=*(rzz+j);rxy_j=*(rxy+j);ryz_j=*(ryz+j);rxz_j=*(rxz+j);
                        rxx_j_2=*(rxx_2+j);ryy_j_2=*(ryy_2+j);rzz_j_2=*(rzz_2+j);rxy_j_2=*(rxy_2+j);ryz_j_2=*(ryz_2+j); rxz_j_2=*(rxz_2+j);
                        rxx_j_3=*(rxx_3+j);ryy_j_3=*(ryy_3+j);rzz_j_3=*(rzz_3+j);rxy_j_3=*(rxy_3+j);ryz_j_3=*(ryz_3+j);rxz_j_3=*(rxz_3+j);
                        vxyyx_j=*(vxyyx+j);vyzzy_j=*(vyzzy+j);vxzzx_j=*(vxzzx+j);vxxyyzz_j=*(vxxyyzz+j);vyyzz_j=*(vyyzz+j);vxxzz_j=*(vxxzz+j);vxxyy_j=*(vxxyy+j);
                        vxyyx_j_2=*(vxyyx_2+j);vyzzy_j_2=*(vyzzy_2+j); vxzzx_j_2=*(vxzzx_2+j);vxxyyzz_j_2=*(vxxyyzz_2+j);vyyzz_j_2=*(vyyzz_2+j);vxxzz_j_2=*(vxxzz_2+j);vxxyy_j_2=*(vxxyy_2+j);
                        vxyyx_j_3=*(vxyyx_3+j);vyzzy_j_3=*(vyzzy_3+j);vxzzx_j_3=*(vxzzx_3+j);vxxyyzz_j_3=*(vxxyyzz_3+j);vyyzz_j_3=*(vyyzz_3+j);vxxzz_j_3=*(vxxzz_3+j);vxxyy_j_3=*(vxxyy_3+j);
                        for (i=nx1;i<=nx2;i++){
                            
                            rxx_j_i=*(rxx_j+i);ryy_j_i=*(ryy_j+i);rzz_j_i=*(rzz_j+i);rxy_j_i=*(rxy_j+i);ryz_j_i=*(ryz_j+i);rxz_j_i=*(rxz_j+i);
                            rxx_j_i_2=*(rxx_j_2+i);ryy_j_i_2=*(ryy_j_2+i);rzz_j_i_2=*(rzz_j_2+i);rxy_j_i_2=*(rxy_j_2+i);ryz_j_i_2=*(ryz_j_2+i);rxz_j_i_2=*(rxz_j_2+i);
                            rxx_j_i_3=*(rxx_j_3+i);ryy_j_i_3=*(ryy_j_3+i);rzz_j_i_3=*(rzz_j_3+i);rxy_j_i_3=*(rxy_j_3+i);ryz_j_i_3=*(ryz_j_3+i);rxz_j_i_3=*(rxz_j_3+i);
                            vxyyx_j_i=*(vxyyx_j+i);vyzzy_j_i=*(vyzzy_j+i);vxzzx_j_i=*(vxzzx_j+i);vxxyyzz_j_i=*(vxxyyzz_j+i);vyyzz_j_i=*(vyyzz_j+i);vxxzz_j_i=*(vxxzz_j+i);vxxyy_j_i=*(vxxyy_j+i);
                            vxyyx_j_i_2=*(vxyyx_j_2+i);vyzzy_j_i_2=*(vyzzy_j_2+i);vxzzx_j_i_2=*(vxzzx_j_2+i);vxxyyzz_j_i_2=*(vxxyyzz_j_2+i);vyyzz_j_i_2=*(vyyzz_j_2+i);vxxzz_j_i_2=*(vxxzz_j_2+i);vxxyy_j_i_2=*(vxxyy_j_2+i);
                            vxyyx_j_i_3=*(vxyyx_j_3+i);vyzzy_j_i_3=*(vyzzy_j_3+i);vxzzx_j_i_3=*(vxzzx_j_3+i); vxxyyzz_j_i_3=*(vxxyyzz_j_3+i);vyyzz_j_i_3=*(vyyzz_j_3+i);vxxzz_j_i_3=*(vxxzz_j_3+i); vxxyy_j_i_3=*(vxxyy_j_3+i);
                            
                            for (k=nz1;k<=nz2;k++){
                                /* spatial derivatives of the components of the velocities
                                 are computed */
                                
                                vxx = (vx[j][i][k]-vx[j][i-1][k])/DX;
                                vxy = (vx[j+1][i][k]-vx[j][i][k])/DY;
                                vxz = (vx[j][i][k+1]-vx[j][i][k])/DZ;
                                vyx = (vy[j][i+1][k]-vy[j][i][k])/DX;
                                vyy = (vy[j][i][k]-vy[j-1][i][k])/DY;
                                vyz = (vy[j][i][k+1]-vy[j][i][k])/DZ;
                                vzx = (vz[j][i+1][k]-vz[j][i][k])/DX;
                                vzy = (vz[j+1][i][k]-vz[j][i][k])/DY;
                                vzz = (vz[j][i][k]-vz[j][i][k-1])/DZ;
                                
                                /* computing sums of the old memory variables in this version only one mechanism is possi
                                 ble */
                                
                                sumrxy=c1*(*(rxy_j_i+k))+c2*(*(rxy_j_i_2+k))+c3*(*(rxy_j_i_3+k));
                                sumryz=c1*(*(ryz_j_i+k))+c2*(*(ryz_j_i_2+k))+c3*(*(ryz_j_i_3+k));
                                sumrxz=c1*(*(rxz_j_i+k))+c2*(*(rxz_j_i_2+k))+c3*(*(rxz_j_i_3+k));
                                sumrxx=c1*(*(rxx_j_i+k))+c2*(*(rxx_j_i_2+k))+c3*(*(rxx_j_i_3+k));
                                sumryy=c1*(*(ryy_j_i+k))+c2*(*(ryy_j_i_2+k))+c3*(*(ryy_j_i_3+k));
                                sumrzz=c1*(*(rzz_j_i+k))+c2*(*(rzz_j_i_2+k))+c3*(*(rzz_j_i_3+k));
                                
                                /* updating components of the stress tensor, partially */
                                fipjp=uipjp[j][i][k]*DT*(1.0+L*tausipjp[j][i][k]);
                                fjpkp=ujpkp[j][i][k]*DT*(1.0+L*tausjpkp[j][i][k]);
                                fipkp=uipkp[j][i][k]*DT*(1.0+L*tausipkp[j][i][k]);
                                g=pi[j][i][k]*(1.0+L*taup[j][i][k]);
                                f=2.0*u[j][i][k]*(1.0+L*taus[j][i][k]);
                                
                                /* Save spatial derivations */
                                *(vxyyx_j_i+k)=vxy+vyx;
                                *(vyzzy_j_i+k)=vyz+vzy;
                                *(vxzzx_j_i+k)=vxz+vzx;
                                *(vxxyyzz_j_i+k)=vxx+vyy+vzz;
                                *(vyyzz_j_i+k)=vyy+vzz;
                                *(vxxzz_j_i+k)=vxx+vzz;
                                *(vxxyy_j_i+k)=vxx+vyy;
                                
                                 /* Do the Adams Bashforth stuff and save it, we need it for the memory variables again */
                                vxyyx_T2=(c1*(*(vxyyx_j_i+k))+c2*(*(vxyyx_j_i_2+k))+c3*(*(vxyyx_j_i_3+k)));
                                vyzzy_T2=(c1*(*(vyzzy_j_i+k))+c2*(*(vyzzy_j_i_2+k))+c3*(*(vyzzy_j_i_3+k)));
                                vxzzx_T2=(c1*(*(vxzzx_j_i+k))+c2*(*(vxzzx_j_i_2+k))+c3*(*(vxzzx_j_i_3+k)));
                                vyyzz_T2=(c1*(*(vyyzz_j_i+k))+c2*(*(vyyzz_j_i_2+k))+c3*(*(vyyzz_j_i_3+k)));
                                vxxzz_T2=(c1*(*(vxxzz_j_i+k))+c2*(*(vxxzz_j_i_2+k))+c3*(*(vxxzz_j_i_3+k)));
                                vxxyy_T2=(c1*(*(vxxyy_j_i+k))+c2*(*(vxxyy_j_i_2+k))+c3*(*(vxxyy_j_i_3+k)));
                                vdiag=((c1*(*(vxxyyzz_j_i+k))+c2*(*(vxxyyzz_j_i_2+k))+c3*(*(vxxyyzz_j_i_3+k))));

                                /* Update the stress tensor */
                                sxy[j][i][k]+=(fipjp*vxyyx_T2)+(dthalbe*sumrxy);
                                syz[j][i][k]+=(fjpkp*vyzzy_T2)+(dthalbe*sumryz);
                                sxz[j][i][k]+=(fipkp*vxzzx_T2)+(dthalbe*sumrxz);
                                sxx[j][i][k]+=DT*((g*vdiag)-(f*vyyzz_T2))+(dthalbe*sumrxx);
                                syy[j][i][k]+=DT*((g*vdiag)-(f*vxxzz_T2))+(dthalbe*sumryy);
                                szz[j][i][k]+=DT*((g*vdiag)-(f*vxxyy_T2))+(dthalbe*sumrzz);
  
                                /* Update the memory variables */
                                
                                /* in the case of L=1 update the memory-variables without suming*/
                                l=1;
                                n1=1.0/(1.0+(c1*eta[l]*0.5)); /* 1/(1+DT*c1/(2*tau_sigl)) */
                                n2=1.0-(c1*eta[l]*0.5);
                                n3=c2*eta[l]/2.0;
                                n4=c3*eta[l]/2.0;
                                dipjp=uipjp[j][i][k]*eta[l]*tausipjp[j][i][k];
                                djpkp=ujpkp[j][i][k]*eta[l]*tausjpkp[j][i][k];
                                dipkp=uipkp[j][i][k]*eta[l]*tausipkp[j][i][k];
                                d=2.0*u[j][i][k]*eta[l]*taus[j][i][k];
                                e=pi[j][i][k]*eta[l]*taup[j][i][k];
                                
                                
                                *(rxy_j_i_3+k)=n1*((*(rxy_j_i+k))*n2-n3*((*(rxy_j_i+k))+(*(rxy_j_i_2+k)))-n4*((*(rxy_j_i_2+k))+(*(rxy_j_i_3+k)))-(dipjp*vxyyx_T2));
                                *(ryz_j_i_3+k)=n1*((*(ryz_j_i+k))*n2-n3*((*(ryz_j_i+k))+(*(ryz_j_i_2+k)))-n4*((*(ryz_j_i_2+k))+(*(ryz_j_i_3+k)))-(djpkp*vyzzy_T2));
                                *(rxz_j_i_3+k)=n1*((*(rxz_j_i+k))*n2-n3*((*(rxz_j_i+k))+(*(rxz_j_i_2+k)))-n4*((*(rxz_j_i_2+k))+(*(rxz_j_i_3+k)))-(dipkp*vxzzx_T2));
                                *(rxx_j_i_3+k)=n1*((*(rxx_j_i+k))*n2-n3*((*(rxx_j_i+k))+(*(rxx_j_i_2+k)))-n4*((*(rxx_j_i_2+k))+(*(rxx_j_i_3+k)))-(e*vdiag)+(d*vyyzz_T2));
                                *(ryy_j_i_3+k)=n1*((*(ryy_j_i+k))*n2-n3*((*(ryy_j_i+k))+(*(ryy_j_i_2+k)))-n4*((*(ryy_j_i_2+k))+(*(ryy_j_i_3+k)))-(e*vdiag)+(d*vxxzz_T2));
                                *(rzz_j_i_3+k)=n1*((*(rzz_j_i+k))*n2-n3*((*(rzz_j_i+k))+(*(rzz_j_i_2+k)))-n4*((*(rzz_j_i_2+k))+(*(rzz_j_i_3+k)))-(e*vdiag)+(d*vxxyy_T2));
                                
                                
                                sumrxy=c1*(*(rxy_j_i_3+k))+c2*(*(rxy_j_i+k))+c3*(*(rxy_j_i_2+k));
                                sumryz=c1*(*(ryz_j_i_3+k))+c2*(*(ryz_j_i+k))+c3*(*(ryz_j_i_2+k));
                                sumrxz=c1*(*(rxz_j_i_3+k))+c2*(*(rxz_j_i+k))+c3*(*(rxz_j_i_2+k));
                                sumrxx=c1*(*(rxx_j_i_3+k))+c2*(*(rxx_j_i+k))+c3*(*(rxx_j_i_2+k));
                                sumryy=c1*(*(ryy_j_i_3+k))+c2*(*(ryy_j_i+k))+c3*(*(ryy_j_i_2+k));
                                sumrzz=c1*(*(rzz_j_i_3+k))+c2*(*(rzz_j_i+k))+c3*(*(rzz_j_i_2+k));

                                
                                /* and now the components of the stress tensor are
                                 completely updated */
                                sxy[j][i][k]+=(dthalbe*sumrxy);
                                syz[j][i][k]+=(dthalbe*sumryz);
                                sxz[j][i][k]+=(dthalbe*sumrxz);
                                sxx[j][i][k]+=(dthalbe*sumrxx);
                                syy[j][i][k]+=(dthalbe*sumryy);
                                szz[j][i][k]+=(dthalbe*sumrzz);
                                
                                
                            }
                        }
                    }
                    break;
                    
                case 4 :
                    
                    b1=9.0/8.0; b2=-1.0/24.0; /* Taylor coefficients */
                    if(FDCOEFF==2){
                        b1=1.1382; b2=-0.046414;} /* Holberg coefficients E=0.1 %*/
                    
                    for (j=ny1;j<=ny2;j++){
                        
                        rxx_j=*(rxx+j);ryy_j=*(ryy+j);rzz_j=*(rzz+j);rxy_j=*(rxy+j);ryz_j=*(ryz+j);rxz_j=*(rxz+j);
                        rxx_j_2=*(rxx_2+j);ryy_j_2=*(ryy_2+j);rzz_j_2=*(rzz_2+j);rxy_j_2=*(rxy_2+j);ryz_j_2=*(ryz_2+j); rxz_j_2=*(rxz_2+j);
                        rxx_j_3=*(rxx_3+j);ryy_j_3=*(ryy_3+j);rzz_j_3=*(rzz_3+j);rxy_j_3=*(rxy_3+j);ryz_j_3=*(ryz_3+j);rxz_j_3=*(rxz_3+j);
                        vxyyx_j=*(vxyyx+j);vyzzy_j=*(vyzzy+j);vxzzx_j=*(vxzzx+j);vxxyyzz_j=*(vxxyyzz+j);vyyzz_j=*(vyyzz+j);vxxzz_j=*(vxxzz+j);vxxyy_j=*(vxxyy+j);
                        vxyyx_j_2=*(vxyyx_2+j);vyzzy_j_2=*(vyzzy_2+j); vxzzx_j_2=*(vxzzx_2+j);vxxyyzz_j_2=*(vxxyyzz_2+j);vyyzz_j_2=*(vyyzz_2+j);vxxzz_j_2=*(vxxzz_2+j);vxxyy_j_2=*(vxxyy_2+j);
                        vxyyx_j_3=*(vxyyx_3+j);vyzzy_j_3=*(vyzzy_3+j);vxzzx_j_3=*(vxzzx_3+j);vxxyyzz_j_3=*(vxxyyzz_3+j);vyyzz_j_3=*(vyyzz_3+j);vxxzz_j_3=*(vxxzz_3+j);vxxyy_j_3=*(vxxyy_3+j);
                        for (i=nx1;i<=nx2;i++){
                            
                            rxx_j_i=*(rxx_j+i);ryy_j_i=*(ryy_j+i);rzz_j_i=*(rzz_j+i);rxy_j_i=*(rxy_j+i);ryz_j_i=*(ryz_j+i);rxz_j_i=*(rxz_j+i);
                            rxx_j_i_2=*(rxx_j_2+i);ryy_j_i_2=*(ryy_j_2+i);rzz_j_i_2=*(rzz_j_2+i);rxy_j_i_2=*(rxy_j_2+i);ryz_j_i_2=*(ryz_j_2+i);rxz_j_i_2=*(rxz_j_2+i);
                            rxx_j_i_3=*(rxx_j_3+i);ryy_j_i_3=*(ryy_j_3+i);rzz_j_i_3=*(rzz_j_3+i);rxy_j_i_3=*(rxy_j_3+i);ryz_j_i_3=*(ryz_j_3+i);rxz_j_i_3=*(rxz_j_3+i);
                            vxyyx_j_i=*(vxyyx_j+i);vyzzy_j_i=*(vyzzy_j+i);vxzzx_j_i=*(vxzzx_j+i);vxxyyzz_j_i=*(vxxyyzz_j+i);vyyzz_j_i=*(vyyzz_j+i);vxxzz_j_i=*(vxxzz_j+i);vxxyy_j_i=*(vxxyy_j+i);
                            vxyyx_j_i_2=*(vxyyx_j_2+i);vyzzy_j_i_2=*(vyzzy_j_2+i);vxzzx_j_i_2=*(vxzzx_j_2+i);vxxyyzz_j_i_2=*(vxxyyzz_j_2+i);vyyzz_j_i_2=*(vyyzz_j_2+i);vxxzz_j_i_2=*(vxxzz_j_2+i);vxxyy_j_i_2=*(vxxyy_j_2+i);
                            vxyyx_j_i_3=*(vxyyx_j_3+i);vyzzy_j_i_3=*(vyzzy_j_3+i);vxzzx_j_i_3=*(vxzzx_j_3+i); vxxyyzz_j_i_3=*(vxxyyzz_j_3+i);vyyzz_j_i_3=*(vyyzz_j_3+i);vxxzz_j_i_3=*(vxxzz_j_3+i); vxxyy_j_i_3=*(vxxyy_j_3+i);
                            
                            for (k=nz1;k<=nz2;k++){

                                /* spatial derivatives of the components of the velocities
                                 are computed */
                                
                                vxx = (b1*(vx[j][i][k]-vx[j][i-1][k])+b2*(vx[j][i+1][k]-vx[j][i-2][k]))/DX;
                                vxy = (b1*(vx[j+1][i][k]-vx[j][i][k])+b2*(vx[j+2][i][k]-vx[j-1][i][k]))/DY;
                                vxz = (b1*(vx[j][i][k+1]-vx[j][i][k])+b2*(vx[j][i][k+2]-vx[j][i][k-1]))/DZ;
                                vyx = (b1*(vy[j][i+1][k]-vy[j][i][k])+b2*(vy[j][i+2][k]-vy[j][i-1][k]))/DX;
                                vyy = (b1*(vy[j][i][k]-vy[j-1][i][k])+b2*(vy[j+1][i][k]-vy[j-2][i][k]))/DY;
                                vyz = (b1*(vy[j][i][k+1]-vy[j][i][k])+b2*(vy[j][i][k+2]-vy[j][i][k-1]))/DZ;
                                vzx = (b1*(vz[j][i+1][k]-vz[j][i][k])+b2*(vz[j][i+2][k]-vz[j][i-1][k]))/DX;
                                vzy = (b1*(vz[j+1][i][k]-vz[j][i][k])+b2*(vz[j+2][i][k]-vz[j-1][i][k]))/DY;
                                vzz = (b1*(vz[j][i][k]-vz[j][i][k-1])+b2*(vz[j][i][k+1]-vz[j][i][k-2]))/DZ;
                                
                                
                                /* computing sums of the old memory variables in this version only one mechanism is possi
                                 ble */
                                
                                sumrxy=c1*(*(rxy_j_i+k))+c2*(*(rxy_j_i_2+k))+c3*(*(rxy_j_i_3+k));
                                sumryz=c1*(*(ryz_j_i+k))+c2*(*(ryz_j_i_2+k))+c3*(*(ryz_j_i_3+k));
                                sumrxz=c1*(*(rxz_j_i+k))+c2*(*(rxz_j_i_2+k))+c3*(*(rxz_j_i_3+k));
                                sumrxx=c1*(*(rxx_j_i+k))+c2*(*(rxx_j_i_2+k))+c3*(*(rxx_j_i_3+k));
                                sumryy=c1*(*(ryy_j_i+k))+c2*(*(ryy_j_i_2+k))+c3*(*(ryy_j_i_3+k));
                                sumrzz=c1*(*(rzz_j_i+k))+c2*(*(rzz_j_i_2+k))+c3*(*(rzz_j_i_3+k));
                                
                                /* updating components of the stress tensor, partially */
                                fipjp=uipjp[j][i][k]*DT*(1.0+L*tausipjp[j][i][k]);
                                fjpkp=ujpkp[j][i][k]*DT*(1.0+L*tausjpkp[j][i][k]);
                                fipkp=uipkp[j][i][k]*DT*(1.0+L*tausipkp[j][i][k]);
                                g=pi[j][i][k]*(1.0+L*taup[j][i][k]);
                                f=2.0*u[j][i][k]*(1.0+L*taus[j][i][k]);
                                
                                /* Save spatial derivations */
                                *(vxyyx_j_i+k)=vxy+vyx;
                                *(vyzzy_j_i+k)=vyz+vzy;
                                *(vxzzx_j_i+k)=vxz+vzx;
                                *(vxxyyzz_j_i+k)=vxx+vyy+vzz;
                                *(vyyzz_j_i+k)=vyy+vzz;
                                *(vxxzz_j_i+k)=vxx+vzz;
                                *(vxxyy_j_i+k)=vxx+vyy;
                                
                                /* Do the Adams Bashforth stuff and save it, we need it for the memory variables again */
                                vxyyx_T2=(c1*(*(vxyyx_j_i+k))+c2*(*(vxyyx_j_i_2+k))+c3*(*(vxyyx_j_i_3+k)));
                                vyzzy_T2=(c1*(*(vyzzy_j_i+k))+c2*(*(vyzzy_j_i_2+k))+c3*(*(vyzzy_j_i_3+k)));
                                vxzzx_T2=(c1*(*(vxzzx_j_i+k))+c2*(*(vxzzx_j_i_2+k))+c3*(*(vxzzx_j_i_3+k)));
                                vyyzz_T2=(c1*(*(vyyzz_j_i+k))+c2*(*(vyyzz_j_i_2+k))+c3*(*(vyyzz_j_i_3+k)));
                                vxxzz_T2=(c1*(*(vxxzz_j_i+k))+c2*(*(vxxzz_j_i_2+k))+c3*(*(vxxzz_j_i_3+k)));
                                vxxyy_T2=(c1*(*(vxxyy_j_i+k))+c2*(*(vxxyy_j_i_2+k))+c3*(*(vxxyy_j_i_3+k)));
                                vdiag=((c1*(*(vxxyyzz_j_i+k))+c2*(*(vxxyyzz_j_i_2+k))+c3*(*(vxxyyzz_j_i_3+k))));
                                
                                /* Update the stress tensor */
                                sxy[j][i][k]+=(fipjp*vxyyx_T2)+(dthalbe*sumrxy);
                                syz[j][i][k]+=(fjpkp*vyzzy_T2)+(dthalbe*sumryz);
                                sxz[j][i][k]+=(fipkp*vxzzx_T2)+(dthalbe*sumrxz);
                                sxx[j][i][k]+=DT*((g*vdiag)-(f*vyyzz_T2))+(dthalbe*sumrxx);
                                syy[j][i][k]+=DT*((g*vdiag)-(f*vxxzz_T2))+(dthalbe*sumryy);
                                szz[j][i][k]+=DT*((g*vdiag)-(f*vxxyy_T2))+(dthalbe*sumrzz);
                                
                                /* Update the memory variables */
                                
                                /* in the case of L=1 update the memory-variables without suming*/
                                l=1;
                                n1=1.0/(1.0+(c1*eta[l]*0.5)); /* 1/(1+DT*c1/(2*tau_sigl)) */
                                n2=1.0-(c1*eta[l]*0.5);
                                n3=c2*eta[l]/2.0;
                                n4=c3*eta[l]/2.0;
                                dipjp=uipjp[j][i][k]*eta[l]*tausipjp[j][i][k];
                                djpkp=ujpkp[j][i][k]*eta[l]*tausjpkp[j][i][k];
                                dipkp=uipkp[j][i][k]*eta[l]*tausipkp[j][i][k];
                                d=2.0*u[j][i][k]*eta[l]*taus[j][i][k];
                                e=pi[j][i][k]*eta[l]*taup[j][i][k];
                                
                                
                                *(rxy_j_i_3+k)=n1*((*(rxy_j_i+k))*n2-n3*((*(rxy_j_i+k))+(*(rxy_j_i_2+k)))-n4*((*(rxy_j_i_2+k))+(*(rxy_j_i_3+k)))-(dipjp*vxyyx_T2));
                                *(ryz_j_i_3+k)=n1*((*(ryz_j_i+k))*n2-n3*((*(ryz_j_i+k))+(*(ryz_j_i_2+k)))-n4*((*(ryz_j_i_2+k))+(*(ryz_j_i_3+k)))-(djpkp*vyzzy_T2));
                                *(rxz_j_i_3+k)=n1*((*(rxz_j_i+k))*n2-n3*((*(rxz_j_i+k))+(*(rxz_j_i_2+k)))-n4*((*(rxz_j_i_2+k))+(*(rxz_j_i_3+k)))-(dipkp*vxzzx_T2));
                                *(rxx_j_i_3+k)=n1*((*(rxx_j_i+k))*n2-n3*((*(rxx_j_i+k))+(*(rxx_j_i_2+k)))-n4*((*(rxx_j_i_2+k))+(*(rxx_j_i_3+k)))-(e*vdiag)+(d*vyyzz_T2));
                                *(ryy_j_i_3+k)=n1*((*(ryy_j_i+k))*n2-n3*((*(ryy_j_i+k))+(*(ryy_j_i_2+k)))-n4*((*(ryy_j_i_2+k))+(*(ryy_j_i_3+k)))-(e*vdiag)+(d*vxxzz_T2));
                                *(rzz_j_i_3+k)=n1*((*(rzz_j_i+k))*n2-n3*((*(rzz_j_i+k))+(*(rzz_j_i_2+k)))-n4*((*(rzz_j_i_2+k))+(*(rzz_j_i_3+k)))-(e*vdiag)+(d*vxxyy_T2));
                                
                                
                                sumrxy=c1*(*(rxy_j_i_3+k))+c2*(*(rxy_j_i+k))+c3*(*(rxy_j_i_2+k));
                                sumryz=c1*(*(ryz_j_i_3+k))+c2*(*(ryz_j_i+k))+c3*(*(ryz_j_i_2+k));
                                sumrxz=c1*(*(rxz_j_i_3+k))+c2*(*(rxz_j_i+k))+c3*(*(rxz_j_i_2+k));
                                sumrxx=c1*(*(rxx_j_i_3+k))+c2*(*(rxx_j_i+k))+c3*(*(rxx_j_i_2+k));
                                sumryy=c1*(*(ryy_j_i_3+k))+c2*(*(ryy_j_i+k))+c3*(*(ryy_j_i_2+k));
                                sumrzz=c1*(*(rzz_j_i_3+k))+c2*(*(rzz_j_i+k))+c3*(*(rzz_j_i_2+k));
                                
                                
                                /* and now the components of the stress tensor are
                                 completely updated */
                                sxy[j][i][k]+=(dthalbe*sumrxy);
                                syz[j][i][k]+=(dthalbe*sumryz);
                                sxz[j][i][k]+=(dthalbe*sumrxz);
                                sxx[j][i][k]+=(dthalbe*sumrxx);
                                syy[j][i][k]+=(dthalbe*sumryy);
                                szz[j][i][k]+=(dthalbe*sumrzz);
                                
                                
                                
                            }
                        }
                    }
                    break;
                    
                case 6 :
                    
                    b1=75.0/64.0; b2=-25.0/384.0; b3=3.0/640.0; /* Taylor coefficients */
                    if(FDCOEFF==2){
                        b1=1.1965; b2=-0.078804; b3=0.0081781;}   /* Holberg coefficients E=0.1 %*/
                    
                    
                    for (j=ny1;j<=ny2;j++){
                        
                        rxx_j=*(rxx+j);ryy_j=*(ryy+j);rzz_j=*(rzz+j);rxy_j=*(rxy+j);ryz_j=*(ryz+j);rxz_j=*(rxz+j);
                        rxx_j_2=*(rxx_2+j);ryy_j_2=*(ryy_2+j);rzz_j_2=*(rzz_2+j);rxy_j_2=*(rxy_2+j);ryz_j_2=*(ryz_2+j); rxz_j_2=*(rxz_2+j);
                        rxx_j_3=*(rxx_3+j);ryy_j_3=*(ryy_3+j);rzz_j_3=*(rzz_3+j);rxy_j_3=*(rxy_3+j);ryz_j_3=*(ryz_3+j);rxz_j_3=*(rxz_3+j);
                        vxyyx_j=*(vxyyx+j);vyzzy_j=*(vyzzy+j);vxzzx_j=*(vxzzx+j);vxxyyzz_j=*(vxxyyzz+j);vyyzz_j=*(vyyzz+j);vxxzz_j=*(vxxzz+j);vxxyy_j=*(vxxyy+j);
                        vxyyx_j_2=*(vxyyx_2+j);vyzzy_j_2=*(vyzzy_2+j); vxzzx_j_2=*(vxzzx_2+j);vxxyyzz_j_2=*(vxxyyzz_2+j);vyyzz_j_2=*(vyyzz_2+j);vxxzz_j_2=*(vxxzz_2+j);vxxyy_j_2=*(vxxyy_2+j);
                        vxyyx_j_3=*(vxyyx_3+j);vyzzy_j_3=*(vyzzy_3+j);vxzzx_j_3=*(vxzzx_3+j);vxxyyzz_j_3=*(vxxyyzz_3+j);vyyzz_j_3=*(vyyzz_3+j);vxxzz_j_3=*(vxxzz_3+j);vxxyy_j_3=*(vxxyy_3+j);
                        for (i=nx1;i<=nx2;i++){
                            
                            rxx_j_i=*(rxx_j+i);ryy_j_i=*(ryy_j+i);rzz_j_i=*(rzz_j+i);rxy_j_i=*(rxy_j+i);ryz_j_i=*(ryz_j+i);rxz_j_i=*(rxz_j+i);
                            rxx_j_i_2=*(rxx_j_2+i);ryy_j_i_2=*(ryy_j_2+i);rzz_j_i_2=*(rzz_j_2+i);rxy_j_i_2=*(rxy_j_2+i);ryz_j_i_2=*(ryz_j_2+i);rxz_j_i_2=*(rxz_j_2+i);
                            rxx_j_i_3=*(rxx_j_3+i);ryy_j_i_3=*(ryy_j_3+i);rzz_j_i_3=*(rzz_j_3+i);rxy_j_i_3=*(rxy_j_3+i);ryz_j_i_3=*(ryz_j_3+i);rxz_j_i_3=*(rxz_j_3+i);
                            vxyyx_j_i=*(vxyyx_j+i);vyzzy_j_i=*(vyzzy_j+i);vxzzx_j_i=*(vxzzx_j+i);vxxyyzz_j_i=*(vxxyyzz_j+i);vyyzz_j_i=*(vyyzz_j+i);vxxzz_j_i=*(vxxzz_j+i);vxxyy_j_i=*(vxxyy_j+i);
                            vxyyx_j_i_2=*(vxyyx_j_2+i);vyzzy_j_i_2=*(vyzzy_j_2+i);vxzzx_j_i_2=*(vxzzx_j_2+i);vxxyyzz_j_i_2=*(vxxyyzz_j_2+i);vyyzz_j_i_2=*(vyyzz_j_2+i);vxxzz_j_i_2=*(vxxzz_j_2+i);vxxyy_j_i_2=*(vxxyy_j_2+i);
                            vxyyx_j_i_3=*(vxyyx_j_3+i);vyzzy_j_i_3=*(vyzzy_j_3+i);vxzzx_j_i_3=*(vxzzx_j_3+i); vxxyyzz_j_i_3=*(vxxyyzz_j_3+i);vyyzz_j_i_3=*(vyyzz_j_3+i);vxxzz_j_i_3=*(vxxzz_j_3+i); vxxyy_j_i_3=*(vxxyy_j_3+i);
                            
                            for (k=nz1;k<=nz2;k++){
                                
                                /* spatial derivatives of the components of the velocities
                                 are computed */
                                
                                vxx = (b1*(vx[j][i][k]-vx[j][i-1][k])+
                                       b2*(vx[j][i+1][k]-vx[j][i-2][k])+
                                       b3*(vx[j][i+2][k]-vx[j][i-3][k]))/DX;
                                
                                vxy = (b1*(vx[j+1][i][k]-vx[j][i][k])+
                                       b2*(vx[j+2][i][k]-vx[j-1][i][k])+
                                       b3*(vx[j+3][i][k]-vx[j-2][i][k]))/DY;
                                
                                vxz = (b1*(vx[j][i][k+1]-vx[j][i][k])+
                                       b2*(vx[j][i][k+2]-vx[j][i][k-1])+
                                       b3*(vx[j][i][k+3]-vx[j][i][k-2]))/DZ;
                                
                                vyx = (b1*(vy[j][i+1][k]-vy[j][i][k])+
                                       b2*(vy[j][i+2][k]-vy[j][i-1][k])+
                                       b3*(vy[j][i+3][k]-vy[j][i-2][k]))/DX;
                                
                                vyy = (b1*(vy[j][i][k]-vy[j-1][i][k])+
                                       b2*(vy[j+1][i][k]-vy[j-2][i][k])+
                                       b3*(vy[j+2][i][k]-vy[j-3][i][k]))/DY;
                                
                                vyz = (b1*(vy[j][i][k+1]-vy[j][i][k])+
                                       b2*(vy[j][i][k+2]-vy[j][i][k-1])+
                                       b3*(vy[j][i][k+3]-vy[j][i][k-2]))/DZ;
                                
                                vzx = (b1*(vz[j][i+1][k]-vz[j][i][k])+
                                       b2*(vz[j][i+2][k]-vz[j][i-1][k])+
                                       b3*(vz[j][i+3][k]-vz[j][i-2][k]))/DX;
                                
                                vzy = (b1*(vz[j+1][i][k]-vz[j][i][k])+
                                       b2*(vz[j+2][i][k]-vz[j-1][i][k])+
                                       b3*(vz[j+3][i][k]-vz[j-2][i][k]))/DY;
                                
                                vzz = (b1*(vz[j][i][k]-vz[j][i][k-1])+
                                       b2*(vz[j][i][k+1]-vz[j][i][k-2])+
                                       b3*(vz[j][i][k+2]-vz[j][i][k-3]))/DZ;
                                
                                
                                /* computing sums of the old memory variables in this version only one mechanism is possi
                                 ble */
                                
                                sumrxy=c1*(*(rxy_j_i+k))+c2*(*(rxy_j_i_2+k))+c3*(*(rxy_j_i_3+k));
                                sumryz=c1*(*(ryz_j_i+k))+c2*(*(ryz_j_i_2+k))+c3*(*(ryz_j_i_3+k));
                                sumrxz=c1*(*(rxz_j_i+k))+c2*(*(rxz_j_i_2+k))+c3*(*(rxz_j_i_3+k));
                                sumrxx=c1*(*(rxx_j_i+k))+c2*(*(rxx_j_i_2+k))+c3*(*(rxx_j_i_3+k));
                                sumryy=c1*(*(ryy_j_i+k))+c2*(*(ryy_j_i_2+k))+c3*(*(ryy_j_i_3+k));
                                sumrzz=c1*(*(rzz_j_i+k))+c2*(*(rzz_j_i_2+k))+c3*(*(rzz_j_i_3+k));
                                
                                /* updating components of the stress tensor, partially */
                                fipjp=uipjp[j][i][k]*DT*(1.0+L*tausipjp[j][i][k]);
                                fjpkp=ujpkp[j][i][k]*DT*(1.0+L*tausjpkp[j][i][k]);
                                fipkp=uipkp[j][i][k]*DT*(1.0+L*tausipkp[j][i][k]);
                                g=pi[j][i][k]*(1.0+L*taup[j][i][k]);
                                f=2.0*u[j][i][k]*(1.0+L*taus[j][i][k]);
                                
                                /* Save spatial derivations */
                                *(vxyyx_j_i+k)=vxy+vyx;
                                *(vyzzy_j_i+k)=vyz+vzy;
                                *(vxzzx_j_i+k)=vxz+vzx;
                                *(vxxyyzz_j_i+k)=vxx+vyy+vzz;
                                *(vyyzz_j_i+k)=vyy+vzz;
                                *(vxxzz_j_i+k)=vxx+vzz;
                                *(vxxyy_j_i+k)=vxx+vyy;
                                
                                /* Do the Adams Bashforth stuff and save it, we need it for the memory variables again */
                                vxyyx_T2=(c1*(*(vxyyx_j_i+k))+c2*(*(vxyyx_j_i_2+k))+c3*(*(vxyyx_j_i_3+k)));
                                vyzzy_T2=(c1*(*(vyzzy_j_i+k))+c2*(*(vyzzy_j_i_2+k))+c3*(*(vyzzy_j_i_3+k)));
                                vxzzx_T2=(c1*(*(vxzzx_j_i+k))+c2*(*(vxzzx_j_i_2+k))+c3*(*(vxzzx_j_i_3+k)));
                                vyyzz_T2=(c1*(*(vyyzz_j_i+k))+c2*(*(vyyzz_j_i_2+k))+c3*(*(vyyzz_j_i_3+k)));
                                vxxzz_T2=(c1*(*(vxxzz_j_i+k))+c2*(*(vxxzz_j_i_2+k))+c3*(*(vxxzz_j_i_3+k)));
                                vxxyy_T2=(c1*(*(vxxyy_j_i+k))+c2*(*(vxxyy_j_i_2+k))+c3*(*(vxxyy_j_i_3+k)));
                                vdiag=((c1*(*(vxxyyzz_j_i+k))+c2*(*(vxxyyzz_j_i_2+k))+c3*(*(vxxyyzz_j_i_3+k))));
                                
                                /* Update the stress tensor */
                                sxy[j][i][k]+=(fipjp*vxyyx_T2)+(dthalbe*sumrxy);
                                syz[j][i][k]+=(fjpkp*vyzzy_T2)+(dthalbe*sumryz);
                                sxz[j][i][k]+=(fipkp*vxzzx_T2)+(dthalbe*sumrxz);
                                sxx[j][i][k]+=DT*((g*vdiag)-(f*vyyzz_T2))+(dthalbe*sumrxx);
                                syy[j][i][k]+=DT*((g*vdiag)-(f*vxxzz_T2))+(dthalbe*sumryy);
                                szz[j][i][k]+=DT*((g*vdiag)-(f*vxxyy_T2))+(dthalbe*sumrzz);
                                
                                /* Update the memory variables */
                                
                                /* in the case of L=1 update the memory-variables without suming*/
                                l=1;
                                n1=1.0/(1.0+(c1*eta[l]*0.5)); /* 1/(1+DT*c1/(2*tau_sigl)) */
                                n2=1.0-(c1*eta[l]*0.5);
                                n3=c2*eta[l]/2.0;
                                n4=c3*eta[l]/2.0;
                                dipjp=uipjp[j][i][k]*eta[l]*tausipjp[j][i][k];
                                djpkp=ujpkp[j][i][k]*eta[l]*tausjpkp[j][i][k];
                                dipkp=uipkp[j][i][k]*eta[l]*tausipkp[j][i][k];
                                d=2.0*u[j][i][k]*eta[l]*taus[j][i][k];
                                e=pi[j][i][k]*eta[l]*taup[j][i][k];
                                
                                
                                *(rxy_j_i_3+k)=n1*((*(rxy_j_i+k))*n2-n3*((*(rxy_j_i+k))+(*(rxy_j_i_2+k)))-n4*((*(rxy_j_i_2+k))+(*(rxy_j_i_3+k)))-(dipjp*vxyyx_T2));
                                *(ryz_j_i_3+k)=n1*((*(ryz_j_i+k))*n2-n3*((*(ryz_j_i+k))+(*(ryz_j_i_2+k)))-n4*((*(ryz_j_i_2+k))+(*(ryz_j_i_3+k)))-(djpkp*vyzzy_T2));
                                *(rxz_j_i_3+k)=n1*((*(rxz_j_i+k))*n2-n3*((*(rxz_j_i+k))+(*(rxz_j_i_2+k)))-n4*((*(rxz_j_i_2+k))+(*(rxz_j_i_3+k)))-(dipkp*vxzzx_T2));
                                *(rxx_j_i_3+k)=n1*((*(rxx_j_i+k))*n2-n3*((*(rxx_j_i+k))+(*(rxx_j_i_2+k)))-n4*((*(rxx_j_i_2+k))+(*(rxx_j_i_3+k)))-(e*vdiag)+(d*vyyzz_T2));
                                *(ryy_j_i_3+k)=n1*((*(ryy_j_i+k))*n2-n3*((*(ryy_j_i+k))+(*(ryy_j_i_2+k)))-n4*((*(ryy_j_i_2+k))+(*(ryy_j_i_3+k)))-(e*vdiag)+(d*vxxzz_T2));
                                *(rzz_j_i_3+k)=n1*((*(rzz_j_i+k))*n2-n3*((*(rzz_j_i+k))+(*(rzz_j_i_2+k)))-n4*((*(rzz_j_i_2+k))+(*(rzz_j_i_3+k)))-(e*vdiag)+(d*vxxyy_T2));
                                
                                
                                sumrxy=c1*(*(rxy_j_i_3+k))+c2*(*(rxy_j_i+k))+c3*(*(rxy_j_i_2+k));
                                sumryz=c1*(*(ryz_j_i_3+k))+c2*(*(ryz_j_i+k))+c3*(*(ryz_j_i_2+k));
                                sumrxz=c1*(*(rxz_j_i_3+k))+c2*(*(rxz_j_i+k))+c3*(*(rxz_j_i_2+k));
                                sumrxx=c1*(*(rxx_j_i_3+k))+c2*(*(rxx_j_i+k))+c3*(*(rxx_j_i_2+k));
                                sumryy=c1*(*(ryy_j_i_3+k))+c2*(*(ryy_j_i+k))+c3*(*(ryy_j_i_2+k));
                                sumrzz=c1*(*(rzz_j_i_3+k))+c2*(*(rzz_j_i+k))+c3*(*(rzz_j_i_2+k));
                                
                                
                                /* and now the components of the stress tensor are
                                 completely updated */
                                sxy[j][i][k]+=(dthalbe*sumrxy);
                                syz[j][i][k]+=(dthalbe*sumryz);
                                sxz[j][i][k]+=(dthalbe*sumrxz);
                                sxx[j][i][k]+=(dthalbe*sumrxx);
                                syy[j][i][k]+=(dthalbe*sumryy);
                                szz[j][i][k]+=(dthalbe*sumrzz);
                                
                            }
                        }
                    }
                    break;
                    
                case 8 :
                    
                    b1=1225.0/1024.0; b2=-245.0/3072.0; b3=49.0/5120.0; b4=-5.0/7168.0; /* Taylor coefficients */
                    if(FDCOEFF==2){
                        b1=1.2257; b2=-0.099537; b3=0.018063; b4=-0.0026274;} /* Holberg coefficients E=0.1 %*/
                    
                    
                    for (j=ny1;j<=ny2;j++){
                        
                        rxx_j=*(rxx+j);ryy_j=*(ryy+j);rzz_j=*(rzz+j);rxy_j=*(rxy+j);ryz_j=*(ryz+j);rxz_j=*(rxz+j);
                        rxx_j_2=*(rxx_2+j);ryy_j_2=*(ryy_2+j);rzz_j_2=*(rzz_2+j);rxy_j_2=*(rxy_2+j);ryz_j_2=*(ryz_2+j); rxz_j_2=*(rxz_2+j);
                        rxx_j_3=*(rxx_3+j);ryy_j_3=*(ryy_3+j);rzz_j_3=*(rzz_3+j);rxy_j_3=*(rxy_3+j);ryz_j_3=*(ryz_3+j);rxz_j_3=*(rxz_3+j);
                        vxyyx_j=*(vxyyx+j);vyzzy_j=*(vyzzy+j);vxzzx_j=*(vxzzx+j);vxxyyzz_j=*(vxxyyzz+j);vyyzz_j=*(vyyzz+j);vxxzz_j=*(vxxzz+j);vxxyy_j=*(vxxyy+j);
                        vxyyx_j_2=*(vxyyx_2+j);vyzzy_j_2=*(vyzzy_2+j); vxzzx_j_2=*(vxzzx_2+j);vxxyyzz_j_2=*(vxxyyzz_2+j);vyyzz_j_2=*(vyyzz_2+j);vxxzz_j_2=*(vxxzz_2+j);vxxyy_j_2=*(vxxyy_2+j);
                        vxyyx_j_3=*(vxyyx_3+j);vyzzy_j_3=*(vyzzy_3+j);vxzzx_j_3=*(vxzzx_3+j);vxxyyzz_j_3=*(vxxyyzz_3+j);vyyzz_j_3=*(vyyzz_3+j);vxxzz_j_3=*(vxxzz_3+j);vxxyy_j_3=*(vxxyy_3+j);
                        for (i=nx1;i<=nx2;i++){
                            
                            rxx_j_i=*(rxx_j+i);ryy_j_i=*(ryy_j+i);rzz_j_i=*(rzz_j+i);rxy_j_i=*(rxy_j+i);ryz_j_i=*(ryz_j+i);rxz_j_i=*(rxz_j+i);
                            rxx_j_i_2=*(rxx_j_2+i);ryy_j_i_2=*(ryy_j_2+i);rzz_j_i_2=*(rzz_j_2+i);rxy_j_i_2=*(rxy_j_2+i);ryz_j_i_2=*(ryz_j_2+i);rxz_j_i_2=*(rxz_j_2+i);
                            rxx_j_i_3=*(rxx_j_3+i);ryy_j_i_3=*(ryy_j_3+i);rzz_j_i_3=*(rzz_j_3+i);rxy_j_i_3=*(rxy_j_3+i);ryz_j_i_3=*(ryz_j_3+i);rxz_j_i_3=*(rxz_j_3+i);
                            vxyyx_j_i=*(vxyyx_j+i);vyzzy_j_i=*(vyzzy_j+i);vxzzx_j_i=*(vxzzx_j+i);vxxyyzz_j_i=*(vxxyyzz_j+i);vyyzz_j_i=*(vyyzz_j+i);vxxzz_j_i=*(vxxzz_j+i);vxxyy_j_i=*(vxxyy_j+i);
                            vxyyx_j_i_2=*(vxyyx_j_2+i);vyzzy_j_i_2=*(vyzzy_j_2+i);vxzzx_j_i_2=*(vxzzx_j_2+i);vxxyyzz_j_i_2=*(vxxyyzz_j_2+i);vyyzz_j_i_2=*(vyyzz_j_2+i);vxxzz_j_i_2=*(vxxzz_j_2+i);vxxyy_j_i_2=*(vxxyy_j_2+i);
                            vxyyx_j_i_3=*(vxyyx_j_3+i);vyzzy_j_i_3=*(vyzzy_j_3+i);vxzzx_j_i_3=*(vxzzx_j_3+i); vxxyyzz_j_i_3=*(vxxyyzz_j_3+i);vyyzz_j_i_3=*(vyyzz_j_3+i);vxxzz_j_i_3=*(vxxzz_j_3+i); vxxyy_j_i_3=*(vxxyy_j_3+i);
                            
                            for (k=nz1;k<=nz2;k++){
                                /* spatial derivatives of the components of the velocities
                                 are computed */
                                
                                vxx = (b1*(vx[j][i][k]-vx[j][i-1][k])+
                                       b2*(vx[j][i+1][k]-vx[j][i-2][k])+
                                       b3*(vx[j][i+2][k]-vx[j][i-3][k])+
                                       b4*(vx[j][i+3][k]-vx[j][i-4][k]))/DX;
                                
                                vxy = (b1*(vx[j+1][i][k]-vx[j][i][k])+
                                       b2*(vx[j+2][i][k]-vx[j-1][i][k])+
                                       b3*(vx[j+3][i][k]-vx[j-2][i][k])+
                                       b4*(vx[j+4][i][k]-vx[j-3][i][k]))/DY;
                                
                                vxz = (b1*(vx[j][i][k+1]-vx[j][i][k])+
                                       b2*(vx[j][i][k+2]-vx[j][i][k-1])+
                                       b3*(vx[j][i][k+3]-vx[j][i][k-2])+
                                       b4*(vx[j][i][k+4]-vx[j][i][k-3]))/DZ;
                                
                                vyx = (b1*(vy[j][i+1][k]-vy[j][i][k])+
                                       b2*(vy[j][i+2][k]-vy[j][i-1][k])+
                                       b3*(vy[j][i+3][k]-vy[j][i-2][k])+
                                       b4*(vy[j][i+4][k]-vy[j][i-3][k]))/DX;
                                
                                vyy = (b1*(vy[j][i][k]-vy[j-1][i][k])+
                                       b2*(vy[j+1][i][k]-vy[j-2][i][k])+
                                       b3*(vy[j+2][i][k]-vy[j-3][i][k])+
                                       b4*(vy[j+3][i][k]-vy[j-4][i][k]))/DY;
                                
                                vyz = (b1*(vy[j][i][k+1]-vy[j][i][k])+
                                       b2*(vy[j][i][k+2]-vy[j][i][k-1])+
                                       b3*(vy[j][i][k+3]-vy[j][i][k-2])+
                                       b4*(vy[j][i][k+4]-vy[j][i][k-3]))/DZ;
                                
                                vzx = (b1*(vz[j][i+1][k]-vz[j][i][k])+
                                       b2*(vz[j][i+2][k]-vz[j][i-1][k])+
                                       b3*(vz[j][i+3][k]-vz[j][i-2][k])+
                                       b4*(vz[j][i+4][k]-vz[j][i-3][k]))/DX;
                                
                                vzy = (b1*(vz[j+1][i][k]-vz[j][i][k])+
                                       b2*(vz[j+2][i][k]-vz[j-1][i][k])+
                                       b3*(vz[j+3][i][k]-vz[j-2][i][k])+
                                       b4*(vz[j+4][i][k]-vz[j-3][i][k]))/DY;
                                
                                vzz = (b1*(vz[j][i][k]-vz[j][i][k-1])+
                                       b2*(vz[j][i][k+1]-vz[j][i][k-2])+
                                       b3*(vz[j][i][k+2]-vz[j][i][k-3])+
                                       b4*(vz[j][i][k+3]-vz[j][i][k-4]))/DZ;
                                
                                
                                /* computing sums of the old memory variables in this version only one mechanism is possi
                                 ble */
                                
                                sumrxy=c1*(*(rxy_j_i+k))+c2*(*(rxy_j_i_2+k))+c3*(*(rxy_j_i_3+k));
                                sumryz=c1*(*(ryz_j_i+k))+c2*(*(ryz_j_i_2+k))+c3*(*(ryz_j_i_3+k));
                                sumrxz=c1*(*(rxz_j_i+k))+c2*(*(rxz_j_i_2+k))+c3*(*(rxz_j_i_3+k));
                                sumrxx=c1*(*(rxx_j_i+k))+c2*(*(rxx_j_i_2+k))+c3*(*(rxx_j_i_3+k));
                                sumryy=c1*(*(ryy_j_i+k))+c2*(*(ryy_j_i_2+k))+c3*(*(ryy_j_i_3+k));
                                sumrzz=c1*(*(rzz_j_i+k))+c2*(*(rzz_j_i_2+k))+c3*(*(rzz_j_i_3+k));
                                
                                /* updating components of the stress tensor, partially */
                                fipjp=uipjp[j][i][k]*DT*(1.0+L*tausipjp[j][i][k]);
                                fjpkp=ujpkp[j][i][k]*DT*(1.0+L*tausjpkp[j][i][k]);
                                fipkp=uipkp[j][i][k]*DT*(1.0+L*tausipkp[j][i][k]);
                                g=pi[j][i][k]*(1.0+L*taup[j][i][k]);
                                f=2.0*u[j][i][k]*(1.0+L*taus[j][i][k]);
                                
                                /* Save spatial derivations */
                                *(vxyyx_j_i+k)=vxy+vyx;
                                *(vyzzy_j_i+k)=vyz+vzy;
                                *(vxzzx_j_i+k)=vxz+vzx;
                                *(vxxyyzz_j_i+k)=vxx+vyy+vzz;
                                *(vyyzz_j_i+k)=vyy+vzz;
                                *(vxxzz_j_i+k)=vxx+vzz;
                                *(vxxyy_j_i+k)=vxx+vyy;
                                
                                /* Do the Adams Bashforth stuff and save it, we need it for the memory variables again */
                                vxyyx_T2=(c1*(*(vxyyx_j_i+k))+c2*(*(vxyyx_j_i_2+k))+c3*(*(vxyyx_j_i_3+k)));
                                vyzzy_T2=(c1*(*(vyzzy_j_i+k))+c2*(*(vyzzy_j_i_2+k))+c3*(*(vyzzy_j_i_3+k)));
                                vxzzx_T2=(c1*(*(vxzzx_j_i+k))+c2*(*(vxzzx_j_i_2+k))+c3*(*(vxzzx_j_i_3+k)));
                                vyyzz_T2=(c1*(*(vyyzz_j_i+k))+c2*(*(vyyzz_j_i_2+k))+c3*(*(vyyzz_j_i_3+k)));
                                vxxzz_T2=(c1*(*(vxxzz_j_i+k))+c2*(*(vxxzz_j_i_2+k))+c3*(*(vxxzz_j_i_3+k)));
                                vxxyy_T2=(c1*(*(vxxyy_j_i+k))+c2*(*(vxxyy_j_i_2+k))+c3*(*(vxxyy_j_i_3+k)));
                                vdiag=((c1*(*(vxxyyzz_j_i+k))+c2*(*(vxxyyzz_j_i_2+k))+c3*(*(vxxyyzz_j_i_3+k))));
                                
                                /* Update the stress tensor */
                                sxy[j][i][k]+=(fipjp*vxyyx_T2)+(dthalbe*sumrxy);
                                syz[j][i][k]+=(fjpkp*vyzzy_T2)+(dthalbe*sumryz);
                                sxz[j][i][k]+=(fipkp*vxzzx_T2)+(dthalbe*sumrxz);
                                sxx[j][i][k]+=DT*((g*vdiag)-(f*vyyzz_T2))+(dthalbe*sumrxx);
                                syy[j][i][k]+=DT*((g*vdiag)-(f*vxxzz_T2))+(dthalbe*sumryy);
                                szz[j][i][k]+=DT*((g*vdiag)-(f*vxxyy_T2))+(dthalbe*sumrzz);
                                
                                /* Update the memory variables */
                                
                                /* in the case of L=1 update the memory-variables without suming*/
                                l=1;
                                n1=1.0/(1.0+(c1*eta[l]*0.5)); /* 1/(1+DT*c1/(2*tau_sigl)) */
                                n2=1.0-(c1*eta[l]*0.5);
                                n3=c2*eta[l]/2.0;
                                n4=c3*eta[l]/2.0;
                                dipjp=uipjp[j][i][k]*eta[l]*tausipjp[j][i][k];
                                djpkp=ujpkp[j][i][k]*eta[l]*tausjpkp[j][i][k];
                                dipkp=uipkp[j][i][k]*eta[l]*tausipkp[j][i][k];
                                d=2.0*u[j][i][k]*eta[l]*taus[j][i][k];
                                e=pi[j][i][k]*eta[l]*taup[j][i][k];
                                
                                
                                *(rxy_j_i_3+k)=n1*((*(rxy_j_i+k))*n2-n3*((*(rxy_j_i+k))+(*(rxy_j_i_2+k)))-n4*((*(rxy_j_i_2+k))+(*(rxy_j_i_3+k)))-(dipjp*vxyyx_T2));
                                *(ryz_j_i_3+k)=n1*((*(ryz_j_i+k))*n2-n3*((*(ryz_j_i+k))+(*(ryz_j_i_2+k)))-n4*((*(ryz_j_i_2+k))+(*(ryz_j_i_3+k)))-(djpkp*vyzzy_T2));
                                *(rxz_j_i_3+k)=n1*((*(rxz_j_i+k))*n2-n3*((*(rxz_j_i+k))+(*(rxz_j_i_2+k)))-n4*((*(rxz_j_i_2+k))+(*(rxz_j_i_3+k)))-(dipkp*vxzzx_T2));
                                *(rxx_j_i_3+k)=n1*((*(rxx_j_i+k))*n2-n3*((*(rxx_j_i+k))+(*(rxx_j_i_2+k)))-n4*((*(rxx_j_i_2+k))+(*(rxx_j_i_3+k)))-(e*vdiag)+(d*vyyzz_T2));
                                *(ryy_j_i_3+k)=n1*((*(ryy_j_i+k))*n2-n3*((*(ryy_j_i+k))+(*(ryy_j_i_2+k)))-n4*((*(ryy_j_i_2+k))+(*(ryy_j_i_3+k)))-(e*vdiag)+(d*vxxzz_T2));
                                *(rzz_j_i_3+k)=n1*((*(rzz_j_i+k))*n2-n3*((*(rzz_j_i+k))+(*(rzz_j_i_2+k)))-n4*((*(rzz_j_i_2+k))+(*(rzz_j_i_3+k)))-(e*vdiag)+(d*vxxyy_T2));
                                
                                
                                sumrxy=c1*(*(rxy_j_i_3+k))+c2*(*(rxy_j_i+k))+c3*(*(rxy_j_i_2+k));
                                sumryz=c1*(*(ryz_j_i_3+k))+c2*(*(ryz_j_i+k))+c3*(*(ryz_j_i_2+k));
                                sumrxz=c1*(*(rxz_j_i_3+k))+c2*(*(rxz_j_i+k))+c3*(*(rxz_j_i_2+k));
                                sumrxx=c1*(*(rxx_j_i_3+k))+c2*(*(rxx_j_i+k))+c3*(*(rxx_j_i_2+k));
                                sumryy=c1*(*(ryy_j_i_3+k))+c2*(*(ryy_j_i+k))+c3*(*(ryy_j_i_2+k));
                                sumrzz=c1*(*(rzz_j_i_3+k))+c2*(*(rzz_j_i+k))+c3*(*(rzz_j_i_2+k));
                                
                                
                                /* and now the components of the stress tensor are
                                 completely updated */
                                sxy[j][i][k]+=(dthalbe*sumrxy);
                                syz[j][i][k]+=(dthalbe*sumryz);
                                sxz[j][i][k]+=(dthalbe*sumrxz);
                                sxx[j][i][k]+=(dthalbe*sumrxx);
                                syy[j][i][k]+=(dthalbe*sumryy);
                                szz[j][i][k]+=(dthalbe*sumrzz);
                                
                                
                                
                            }
                        }
                    }
                    break;
                    
                case 10 :
                    
                    b1=19845.0/16384.0; b2=-735.0/8192.0; b3=567.0/40960.0; b4=-405.0/229376.0; b5=35.0/294912.0; /* Taylor coefficients */
                    if(FDCOEFF==2){
                        b1=1.2415; b2=-0.11231; b3=0.026191; b4=-0.0064682; b5=0.001191;} /* Holberg coefficients E=0.1 %*/
                    
                    
                    for (j=ny1;j<=ny2;j++){
                        
                        rxx_j=*(rxx+j);ryy_j=*(ryy+j);rzz_j=*(rzz+j);rxy_j=*(rxy+j);ryz_j=*(ryz+j);rxz_j=*(rxz+j);
                        rxx_j_2=*(rxx_2+j);ryy_j_2=*(ryy_2+j);rzz_j_2=*(rzz_2+j);rxy_j_2=*(rxy_2+j);ryz_j_2=*(ryz_2+j); rxz_j_2=*(rxz_2+j);
                        rxx_j_3=*(rxx_3+j);ryy_j_3=*(ryy_3+j);rzz_j_3=*(rzz_3+j);rxy_j_3=*(rxy_3+j);ryz_j_3=*(ryz_3+j);rxz_j_3=*(rxz_3+j);
                        vxyyx_j=*(vxyyx+j);vyzzy_j=*(vyzzy+j);vxzzx_j=*(vxzzx+j);vxxyyzz_j=*(vxxyyzz+j);vyyzz_j=*(vyyzz+j);vxxzz_j=*(vxxzz+j);vxxyy_j=*(vxxyy+j);
                        vxyyx_j_2=*(vxyyx_2+j);vyzzy_j_2=*(vyzzy_2+j); vxzzx_j_2=*(vxzzx_2+j);vxxyyzz_j_2=*(vxxyyzz_2+j);vyyzz_j_2=*(vyyzz_2+j);vxxzz_j_2=*(vxxzz_2+j);vxxyy_j_2=*(vxxyy_2+j);
                        vxyyx_j_3=*(vxyyx_3+j);vyzzy_j_3=*(vyzzy_3+j);vxzzx_j_3=*(vxzzx_3+j);vxxyyzz_j_3=*(vxxyyzz_3+j);vyyzz_j_3=*(vyyzz_3+j);vxxzz_j_3=*(vxxzz_3+j);vxxyy_j_3=*(vxxyy_3+j);
                        for (i=nx1;i<=nx2;i++){
                            
                            rxx_j_i=*(rxx_j+i);ryy_j_i=*(ryy_j+i);rzz_j_i=*(rzz_j+i);rxy_j_i=*(rxy_j+i);ryz_j_i=*(ryz_j+i);rxz_j_i=*(rxz_j+i);
                            rxx_j_i_2=*(rxx_j_2+i);ryy_j_i_2=*(ryy_j_2+i);rzz_j_i_2=*(rzz_j_2+i);rxy_j_i_2=*(rxy_j_2+i);ryz_j_i_2=*(ryz_j_2+i);rxz_j_i_2=*(rxz_j_2+i);
                            rxx_j_i_3=*(rxx_j_3+i);ryy_j_i_3=*(ryy_j_3+i);rzz_j_i_3=*(rzz_j_3+i);rxy_j_i_3=*(rxy_j_3+i);ryz_j_i_3=*(ryz_j_3+i);rxz_j_i_3=*(rxz_j_3+i);
                            vxyyx_j_i=*(vxyyx_j+i);vyzzy_j_i=*(vyzzy_j+i);vxzzx_j_i=*(vxzzx_j+i);vxxyyzz_j_i=*(vxxyyzz_j+i);vyyzz_j_i=*(vyyzz_j+i);vxxzz_j_i=*(vxxzz_j+i);vxxyy_j_i=*(vxxyy_j+i);
                            vxyyx_j_i_2=*(vxyyx_j_2+i);vyzzy_j_i_2=*(vyzzy_j_2+i);vxzzx_j_i_2=*(vxzzx_j_2+i);vxxyyzz_j_i_2=*(vxxyyzz_j_2+i);vyyzz_j_i_2=*(vyyzz_j_2+i);vxxzz_j_i_2=*(vxxzz_j_2+i);vxxyy_j_i_2=*(vxxyy_j_2+i);
                            vxyyx_j_i_3=*(vxyyx_j_3+i);vyzzy_j_i_3=*(vyzzy_j_3+i);vxzzx_j_i_3=*(vxzzx_j_3+i); vxxyyzz_j_i_3=*(vxxyyzz_j_3+i);vyyzz_j_i_3=*(vyyzz_j_3+i);vxxzz_j_i_3=*(vxxzz_j_3+i); vxxyy_j_i_3=*(vxxyy_j_3+i);
                            
                            for (k=nz1;k<=nz2;k++){
                                
                                /* spatial derivatives of the components of the velocities
                                 are computed */
                                
                                vxx = (b1*(vx[j][i][k]-vx[j][i-1][k])+
                                       b2*(vx[j][i+1][k]-vx[j][i-2][k])+
                                       b3*(vx[j][i+2][k]-vx[j][i-3][k])+
                                       b4*(vx[j][i+3][k]-vx[j][i-4][k])+
                                       b5*(vx[j][i+4][k]-vx[j][i-5][k]))/DX;
                                
                                vxy = (b1*(vx[j+1][i][k]-vx[j][i][k])+
                                       b2*(vx[j+2][i][k]-vx[j-1][i][k])+
                                       b3*(vx[j+3][i][k]-vx[j-2][i][k])+
                                       b4*(vx[j+4][i][k]-vx[j-3][i][k])+
                                       b5*(vx[j+5][i][k]-vx[j-4][i][k]))/DY;
                                
                                vxz = (b1*(vx[j][i][k+1]-vx[j][i][k])+
                                       b2*(vx[j][i][k+2]-vx[j][i][k-1])+
                                       b3*(vx[j][i][k+3]-vx[j][i][k-2])+
                                       b4*(vx[j][i][k+4]-vx[j][i][k-3])+
                                       b5*(vx[j][i][k+5]-vx[j][i][k-4]))/DZ;
                                
                                vyx = (b1*(vy[j][i+1][k]-vy[j][i][k])+
                                       b2*(vy[j][i+2][k]-vy[j][i-1][k])+
                                       b3*(vy[j][i+3][k]-vy[j][i-2][k])+
                                       b4*(vy[j][i+4][k]-vy[j][i-3][k])+
                                       b5*(vy[j][i+5][k]-vy[j][i-4][k]))/DX;
                                
                                vyy = (b1*(vy[j][i][k]-vy[j-1][i][k])+
                                       b2*(vy[j+1][i][k]-vy[j-2][i][k])+
                                       b3*(vy[j+2][i][k]-vy[j-3][i][k])+
                                       b4*(vy[j+3][i][k]-vy[j-4][i][k])+
                                       b5*(vy[j+4][i][k]-vy[j-5][i][k]))/DY;
                                
                                vyz = (b1*(vy[j][i][k+1]-vy[j][i][k])+
                                       b2*(vy[j][i][k+2]-vy[j][i][k-1])+
                                       b3*(vy[j][i][k+3]-vy[j][i][k-2])+
                                       b4*(vy[j][i][k+4]-vy[j][i][k-3])+
                                       b5*(vy[j][i][k+5]-vy[j][i][k-4]))/DZ;
                                
                                vzx = (b1*(vz[j][i+1][k]-vz[j][i][k])+
                                       b2*(vz[j][i+2][k]-vz[j][i-1][k])+
                                       b3*(vz[j][i+3][k]-vz[j][i-2][k])+
                                       b4*(vz[j][i+4][k]-vz[j][i-3][k])+
                                       b5*(vz[j][i+5][k]-vz[j][i-4][k]))/DX;
                                
                                vzy = (b1*(vz[j+1][i][k]-vz[j][i][k])+
                                       b2*(vz[j+2][i][k]-vz[j-1][i][k])+
                                       b3*(vz[j+3][i][k]-vz[j-2][i][k])+
                                       b4*(vz[j+4][i][k]-vz[j-3][i][k])+
                                       b5*(vz[j+5][i][k]-vz[j-4][i][k]))/DY;
                                
                                vzz = (b1*(vz[j][i][k]-vz[j][i][k-1])+
                                       b2*(vz[j][i][k+1]-vz[j][i][k-2])+
                                       b3*(vz[j][i][k+2]-vz[j][i][k-3])+
                                       b4*(vz[j][i][k+3]-vz[j][i][k-4])+
                                       b5*(vz[j][i][k+4]-vz[j][i][k-5]))/DZ;
                                
                                
                                /* computing sums of the old memory variables in this version only one mechanism is possi
                                 ble */
                                
                                sumrxy=c1*(*(rxy_j_i+k))+c2*(*(rxy_j_i_2+k))+c3*(*(rxy_j_i_3+k));
                                sumryz=c1*(*(ryz_j_i+k))+c2*(*(ryz_j_i_2+k))+c3*(*(ryz_j_i_3+k));
                                sumrxz=c1*(*(rxz_j_i+k))+c2*(*(rxz_j_i_2+k))+c3*(*(rxz_j_i_3+k));
                                sumrxx=c1*(*(rxx_j_i+k))+c2*(*(rxx_j_i_2+k))+c3*(*(rxx_j_i_3+k));
                                sumryy=c1*(*(ryy_j_i+k))+c2*(*(ryy_j_i_2+k))+c3*(*(ryy_j_i_3+k));
                                sumrzz=c1*(*(rzz_j_i+k))+c2*(*(rzz_j_i_2+k))+c3*(*(rzz_j_i_3+k));
                                
                                /* updating components of the stress tensor, partially */
                                fipjp=uipjp[j][i][k]*DT*(1.0+L*tausipjp[j][i][k]);
                                fjpkp=ujpkp[j][i][k]*DT*(1.0+L*tausjpkp[j][i][k]);
                                fipkp=uipkp[j][i][k]*DT*(1.0+L*tausipkp[j][i][k]);
                                g=pi[j][i][k]*(1.0+L*taup[j][i][k]);
                                f=2.0*u[j][i][k]*(1.0+L*taus[j][i][k]);
                                
                                /* Save spatial derivations */
                                *(vxyyx_j_i+k)=vxy+vyx;
                                *(vyzzy_j_i+k)=vyz+vzy;
                                *(vxzzx_j_i+k)=vxz+vzx;
                                *(vxxyyzz_j_i+k)=vxx+vyy+vzz;
                                *(vyyzz_j_i+k)=vyy+vzz;
                                *(vxxzz_j_i+k)=vxx+vzz;
                                *(vxxyy_j_i+k)=vxx+vyy;
                                
                                /* Do the Adams Bashforth stuff and save it, we need it for the memory variables again */
                                vxyyx_T2=(c1*(*(vxyyx_j_i+k))+c2*(*(vxyyx_j_i_2+k))+c3*(*(vxyyx_j_i_3+k)));
                                vyzzy_T2=(c1*(*(vyzzy_j_i+k))+c2*(*(vyzzy_j_i_2+k))+c3*(*(vyzzy_j_i_3+k)));
                                vxzzx_T2=(c1*(*(vxzzx_j_i+k))+c2*(*(vxzzx_j_i_2+k))+c3*(*(vxzzx_j_i_3+k)));
                                vyyzz_T2=(c1*(*(vyyzz_j_i+k))+c2*(*(vyyzz_j_i_2+k))+c3*(*(vyyzz_j_i_3+k)));
                                vxxzz_T2=(c1*(*(vxxzz_j_i+k))+c2*(*(vxxzz_j_i_2+k))+c3*(*(vxxzz_j_i_3+k)));
                                vxxyy_T2=(c1*(*(vxxyy_j_i+k))+c2*(*(vxxyy_j_i_2+k))+c3*(*(vxxyy_j_i_3+k)));
                                vdiag=((c1*(*(vxxyyzz_j_i+k))+c2*(*(vxxyyzz_j_i_2+k))+c3*(*(vxxyyzz_j_i_3+k))));
                                
                                /* Update the stress tensor */
                                sxy[j][i][k]+=(fipjp*vxyyx_T2)+(dthalbe*sumrxy);
                                syz[j][i][k]+=(fjpkp*vyzzy_T2)+(dthalbe*sumryz);
                                sxz[j][i][k]+=(fipkp*vxzzx_T2)+(dthalbe*sumrxz);
                                sxx[j][i][k]+=DT*((g*vdiag)-(f*vyyzz_T2))+(dthalbe*sumrxx);
                                syy[j][i][k]+=DT*((g*vdiag)-(f*vxxzz_T2))+(dthalbe*sumryy);
                                szz[j][i][k]+=DT*((g*vdiag)-(f*vxxyy_T2))+(dthalbe*sumrzz);
                                
                                /* Update the memory variables */
                                
                                /* in the case of L=1 update the memory-variables without suming*/
                                l=1;
                                n1=1.0/(1.0+(c1*eta[l]*0.5)); /* 1/(1+DT*c1/(2*tau_sigl)) */
                                n2=1.0-(c1*eta[l]*0.5);
                                n3=c2*eta[l]/2.0;
                                n4=c3*eta[l]/2.0;
                                dipjp=uipjp[j][i][k]*eta[l]*tausipjp[j][i][k];
                                djpkp=ujpkp[j][i][k]*eta[l]*tausjpkp[j][i][k];
                                dipkp=uipkp[j][i][k]*eta[l]*tausipkp[j][i][k];
                                d=2.0*u[j][i][k]*eta[l]*taus[j][i][k];
                                e=pi[j][i][k]*eta[l]*taup[j][i][k];
                                
                                
                                *(rxy_j_i_3+k)=n1*((*(rxy_j_i+k))*n2-n3*((*(rxy_j_i+k))+(*(rxy_j_i_2+k)))-n4*((*(rxy_j_i_2+k))+(*(rxy_j_i_3+k)))-(dipjp*vxyyx_T2));
                                *(ryz_j_i_3+k)=n1*((*(ryz_j_i+k))*n2-n3*((*(ryz_j_i+k))+(*(ryz_j_i_2+k)))-n4*((*(ryz_j_i_2+k))+(*(ryz_j_i_3+k)))-(djpkp*vyzzy_T2));
                                *(rxz_j_i_3+k)=n1*((*(rxz_j_i+k))*n2-n3*((*(rxz_j_i+k))+(*(rxz_j_i_2+k)))-n4*((*(rxz_j_i_2+k))+(*(rxz_j_i_3+k)))-(dipkp*vxzzx_T2));
                                *(rxx_j_i_3+k)=n1*((*(rxx_j_i+k))*n2-n3*((*(rxx_j_i+k))+(*(rxx_j_i_2+k)))-n4*((*(rxx_j_i_2+k))+(*(rxx_j_i_3+k)))-(e*vdiag)+(d*vyyzz_T2));
                                *(ryy_j_i_3+k)=n1*((*(ryy_j_i+k))*n2-n3*((*(ryy_j_i+k))+(*(ryy_j_i_2+k)))-n4*((*(ryy_j_i_2+k))+(*(ryy_j_i_3+k)))-(e*vdiag)+(d*vxxzz_T2));
                                *(rzz_j_i_3+k)=n1*((*(rzz_j_i+k))*n2-n3*((*(rzz_j_i+k))+(*(rzz_j_i_2+k)))-n4*((*(rzz_j_i_2+k))+(*(rzz_j_i_3+k)))-(e*vdiag)+(d*vxxyy_T2));
                                
                                
                                sumrxy=c1*(*(rxy_j_i_3+k))+c2*(*(rxy_j_i+k))+c3*(*(rxy_j_i_2+k));
                                sumryz=c1*(*(ryz_j_i_3+k))+c2*(*(ryz_j_i+k))+c3*(*(ryz_j_i_2+k));
                                sumrxz=c1*(*(rxz_j_i_3+k))+c2*(*(rxz_j_i+k))+c3*(*(rxz_j_i_2+k));
                                sumrxx=c1*(*(rxx_j_i_3+k))+c2*(*(rxx_j_i+k))+c3*(*(rxx_j_i_2+k));
                                sumryy=c1*(*(ryy_j_i_3+k))+c2*(*(ryy_j_i+k))+c3*(*(ryy_j_i_2+k));
                                sumrzz=c1*(*(rzz_j_i_3+k))+c2*(*(rzz_j_i+k))+c3*(*(rzz_j_i_2+k));
                                
                                
                                /* and now the components of the stress tensor are
                                 completely updated */
                                sxy[j][i][k]+=(dthalbe*sumrxy);
                                syz[j][i][k]+=(dthalbe*sumryz);
                                sxz[j][i][k]+=(dthalbe*sumrxz);
                                sxx[j][i][k]+=(dthalbe*sumrxx);
                                syy[j][i][k]+=(dthalbe*sumryy);
                                szz[j][i][k]+=(dthalbe*sumrzz);
                                
                                
                                
                            }
                        }
                    }
                    break;
                    
                case 12 :
                    
                    /* Taylor coefficients */
                    b1=160083.0/131072.0; b2=-12705.0/131072.0; b3=22869.0/1310720.0;
                    b4=-5445.0/1835008.0; b5=847.0/2359296.0; b6=-63.0/2883584;
                    
                    /* Holberg coefficients E=0.1 %*/
                    if(FDCOEFF==2){
                        b1=1.2508; b2=-0.12034; b3=0.032131; b4=-0.010142; b5=0.0029857; b6=-0.00066667;}
                    
                    
                    for (j=ny1;j<=ny2;j++){
                        
                        rxx_j=*(rxx+j);ryy_j=*(ryy+j);rzz_j=*(rzz+j);rxy_j=*(rxy+j);ryz_j=*(ryz+j);rxz_j=*(rxz+j);
                        rxx_j_2=*(rxx_2+j);ryy_j_2=*(ryy_2+j);rzz_j_2=*(rzz_2+j);rxy_j_2=*(rxy_2+j);ryz_j_2=*(ryz_2+j); rxz_j_2=*(rxz_2+j);
                        rxx_j_3=*(rxx_3+j);ryy_j_3=*(ryy_3+j);rzz_j_3=*(rzz_3+j);rxy_j_3=*(rxy_3+j);ryz_j_3=*(ryz_3+j);rxz_j_3=*(rxz_3+j);
                        vxyyx_j=*(vxyyx+j);vyzzy_j=*(vyzzy+j);vxzzx_j=*(vxzzx+j);vxxyyzz_j=*(vxxyyzz+j);vyyzz_j=*(vyyzz+j);vxxzz_j=*(vxxzz+j);vxxyy_j=*(vxxyy+j);
                        vxyyx_j_2=*(vxyyx_2+j);vyzzy_j_2=*(vyzzy_2+j); vxzzx_j_2=*(vxzzx_2+j);vxxyyzz_j_2=*(vxxyyzz_2+j);vyyzz_j_2=*(vyyzz_2+j);vxxzz_j_2=*(vxxzz_2+j);vxxyy_j_2=*(vxxyy_2+j);
                        vxyyx_j_3=*(vxyyx_3+j);vyzzy_j_3=*(vyzzy_3+j);vxzzx_j_3=*(vxzzx_3+j);vxxyyzz_j_3=*(vxxyyzz_3+j);vyyzz_j_3=*(vyyzz_3+j);vxxzz_j_3=*(vxxzz_3+j);vxxyy_j_3=*(vxxyy_3+j);
                        for (i=nx1;i<=nx2;i++){
                            
                            rxx_j_i=*(rxx_j+i);ryy_j_i=*(ryy_j+i);rzz_j_i=*(rzz_j+i);rxy_j_i=*(rxy_j+i);ryz_j_i=*(ryz_j+i);rxz_j_i=*(rxz_j+i);
                            rxx_j_i_2=*(rxx_j_2+i);ryy_j_i_2=*(ryy_j_2+i);rzz_j_i_2=*(rzz_j_2+i);rxy_j_i_2=*(rxy_j_2+i);ryz_j_i_2=*(ryz_j_2+i);rxz_j_i_2=*(rxz_j_2+i);
                            rxx_j_i_3=*(rxx_j_3+i);ryy_j_i_3=*(ryy_j_3+i);rzz_j_i_3=*(rzz_j_3+i);rxy_j_i_3=*(rxy_j_3+i);ryz_j_i_3=*(ryz_j_3+i);rxz_j_i_3=*(rxz_j_3+i);
                            vxyyx_j_i=*(vxyyx_j+i);vyzzy_j_i=*(vyzzy_j+i);vxzzx_j_i=*(vxzzx_j+i);vxxyyzz_j_i=*(vxxyyzz_j+i);vyyzz_j_i=*(vyyzz_j+i);vxxzz_j_i=*(vxxzz_j+i);vxxyy_j_i=*(vxxyy_j+i);
                            vxyyx_j_i_2=*(vxyyx_j_2+i);vyzzy_j_i_2=*(vyzzy_j_2+i);vxzzx_j_i_2=*(vxzzx_j_2+i);vxxyyzz_j_i_2=*(vxxyyzz_j_2+i);vyyzz_j_i_2=*(vyyzz_j_2+i);vxxzz_j_i_2=*(vxxzz_j_2+i);vxxyy_j_i_2=*(vxxyy_j_2+i);
                            vxyyx_j_i_3=*(vxyyx_j_3+i);vyzzy_j_i_3=*(vyzzy_j_3+i);vxzzx_j_i_3=*(vxzzx_j_3+i); vxxyyzz_j_i_3=*(vxxyyzz_j_3+i);vyyzz_j_i_3=*(vyyzz_j_3+i);vxxzz_j_i_3=*(vxxzz_j_3+i); vxxyy_j_i_3=*(vxxyy_j_3+i);
                            
                            for (k=nz1;k<=nz2;k++){
                                
                                /* spatial derivatives of the components of the velocities
                                 are computed */
                                
                                vxx = (b1*(vx[j][i][k]-vx[j][i-1][k])+
                                       b2*(vx[j][i+1][k]-vx[j][i-2][k])+
                                       b3*(vx[j][i+2][k]-vx[j][i-3][k])+
                                       b4*(vx[j][i+3][k]-vx[j][i-4][k])+
                                       b5*(vx[j][i+4][k]-vx[j][i-5][k])+
                                       b6*(vx[j][i+5][k]-vx[j][i-6][k]))/DX;
                                
                                vxy = (b1*(vx[j+1][i][k]-vx[j][i][k])+
                                       b2*(vx[j+2][i][k]-vx[j-1][i][k])+
                                       b3*(vx[j+3][i][k]-vx[j-2][i][k])+
                                       b4*(vx[j+4][i][k]-vx[j-3][i][k])+
                                       b5*(vx[j+5][i][k]-vx[j-4][i][k])+
                                       b6*(vx[j+6][i][k]-vx[j-5][i][k]))/DY;
                                
                                vxz = (b1*(vx[j][i][k+1]-vx[j][i][k])+
                                       b2*(vx[j][i][k+2]-vx[j][i][k-1])+
                                       b3*(vx[j][i][k+3]-vx[j][i][k-2])+
                                       b4*(vx[j][i][k+4]-vx[j][i][k-3])+
                                       b5*(vx[j][i][k+5]-vx[j][i][k-4])+
                                       b6*(vx[j][i][k+6]-vx[j][i][k-5]))/DZ;
                                
                                vyx = (b1*(vy[j][i+1][k]-vy[j][i][k])+
                                       b2*(vy[j][i+2][k]-vy[j][i-1][k])+
                                       b3*(vy[j][i+3][k]-vy[j][i-2][k])+
                                       b4*(vy[j][i+4][k]-vy[j][i-3][k])+
                                       b5*(vy[j][i+5][k]-vy[j][i-4][k])+
                                       b6*(vy[j][i+6][k]-vy[j][i-5][k]))/DX;
                                
                                vyy = (b1*(vy[j][i][k]-vy[j-1][i][k])+
                                       b2*(vy[j+1][i][k]-vy[j-2][i][k])+
                                       b3*(vy[j+2][i][k]-vy[j-3][i][k])+
                                       b4*(vy[j+3][i][k]-vy[j-4][i][k])+
                                       b5*(vy[j+4][i][k]-vy[j-5][i][k])+
                                       b6*(vy[j+5][i][k]-vy[j-6][i][k]))/DY;
                                
                                vyz = (b1*(vy[j][i][k+1]-vy[j][i][k])+
                                       b2*(vy[j][i][k+2]-vy[j][i][k-1])+
                                       b3*(vy[j][i][k+3]-vy[j][i][k-2])+
                                       b4*(vy[j][i][k+4]-vy[j][i][k-3])+
                                       b5*(vy[j][i][k+5]-vy[j][i][k-4])+
                                       b6*(vy[j][i][k+6]-vy[j][i][k-5]))/DZ;
                                
                                vzx = (b1*(vz[j][i+1][k]-vz[j][i][k])+
                                       b2*(vz[j][i+2][k]-vz[j][i-1][k])+
                                       b3*(vz[j][i+3][k]-vz[j][i-2][k])+
                                       b4*(vz[j][i+4][k]-vz[j][i-3][k])+
                                       b5*(vz[j][i+5][k]-vz[j][i-4][k])+
                                       b6*(vz[j][i+6][k]-vz[j][i-5][k]))/DX;
                                
                                vzy = (b1*(vz[j+1][i][k]-vz[j][i][k])+
                                       b2*(vz[j+2][i][k]-vz[j-1][i][k])+
                                       b3*(vz[j+3][i][k]-vz[j-2][i][k])+
                                       b4*(vz[j+4][i][k]-vz[j-3][i][k])+
                                       b5*(vz[j+5][i][k]-vz[j-4][i][k])+
                                       b6*(vz[j+6][i][k]-vz[j-5][i][k]))/DY;
                                
                                vzz = (b1*(vz[j][i][k]-vz[j][i][k-1])+
                                       b2*(vz[j][i][k+1]-vz[j][i][k-2])+
                                       b3*(vz[j][i][k+2]-vz[j][i][k-3])+
                                       b4*(vz[j][i][k+3]-vz[j][i][k-4])+
                                       b5*(vz[j][i][k+4]-vz[j][i][k-5])+
                                       b6*(vz[j][i][k+5]-vz[j][i][k-6]))/DZ;
                                
                                
                                /* computing sums of the old memory variables in this version only one mechanism is possi
                                 ble */
                                
                                sumrxy=c1*(*(rxy_j_i+k))+c2*(*(rxy_j_i_2+k))+c3*(*(rxy_j_i_3+k));
                                sumryz=c1*(*(ryz_j_i+k))+c2*(*(ryz_j_i_2+k))+c3*(*(ryz_j_i_3+k));
                                sumrxz=c1*(*(rxz_j_i+k))+c2*(*(rxz_j_i_2+k))+c3*(*(rxz_j_i_3+k));
                                sumrxx=c1*(*(rxx_j_i+k))+c2*(*(rxx_j_i_2+k))+c3*(*(rxx_j_i_3+k));
                                sumryy=c1*(*(ryy_j_i+k))+c2*(*(ryy_j_i_2+k))+c3*(*(ryy_j_i_3+k));
                                sumrzz=c1*(*(rzz_j_i+k))+c2*(*(rzz_j_i_2+k))+c3*(*(rzz_j_i_3+k));
                                
                                /* updating components of the stress tensor, partially */
                                fipjp=uipjp[j][i][k]*DT*(1.0+L*tausipjp[j][i][k]);
                                fjpkp=ujpkp[j][i][k]*DT*(1.0+L*tausjpkp[j][i][k]);
                                fipkp=uipkp[j][i][k]*DT*(1.0+L*tausipkp[j][i][k]);
                                g=pi[j][i][k]*(1.0+L*taup[j][i][k]);
                                f=2.0*u[j][i][k]*(1.0+L*taus[j][i][k]);
                                
                                /* Save spatial derivations */
                                *(vxyyx_j_i+k)=vxy+vyx;
                                *(vyzzy_j_i+k)=vyz+vzy;
                                *(vxzzx_j_i+k)=vxz+vzx;
                                *(vxxyyzz_j_i+k)=vxx+vyy+vzz;
                                *(vyyzz_j_i+k)=vyy+vzz;
                                *(vxxzz_j_i+k)=vxx+vzz;
                                *(vxxyy_j_i+k)=vxx+vyy;
                                
                                /* Do the Adams Bashforth stuff and save it, we need it for the memory variables again */
                                vxyyx_T2=(c1*(*(vxyyx_j_i+k))+c2*(*(vxyyx_j_i_2+k))+c3*(*(vxyyx_j_i_3+k)));
                                vyzzy_T2=(c1*(*(vyzzy_j_i+k))+c2*(*(vyzzy_j_i_2+k))+c3*(*(vyzzy_j_i_3+k)));
                                vxzzx_T2=(c1*(*(vxzzx_j_i+k))+c2*(*(vxzzx_j_i_2+k))+c3*(*(vxzzx_j_i_3+k)));
                                vyyzz_T2=(c1*(*(vyyzz_j_i+k))+c2*(*(vyyzz_j_i_2+k))+c3*(*(vyyzz_j_i_3+k)));
                                vxxzz_T2=(c1*(*(vxxzz_j_i+k))+c2*(*(vxxzz_j_i_2+k))+c3*(*(vxxzz_j_i_3+k)));
                                vxxyy_T2=(c1*(*(vxxyy_j_i+k))+c2*(*(vxxyy_j_i_2+k))+c3*(*(vxxyy_j_i_3+k)));
                                vdiag=((c1*(*(vxxyyzz_j_i+k))+c2*(*(vxxyyzz_j_i_2+k))+c3*(*(vxxyyzz_j_i_3+k))));
                                
                                /* Update the stress tensor */
                                sxy[j][i][k]+=(fipjp*vxyyx_T2)+(dthalbe*sumrxy);
                                syz[j][i][k]+=(fjpkp*vyzzy_T2)+(dthalbe*sumryz);
                                sxz[j][i][k]+=(fipkp*vxzzx_T2)+(dthalbe*sumrxz);
                                sxx[j][i][k]+=DT*((g*vdiag)-(f*vyyzz_T2))+(dthalbe*sumrxx);
                                syy[j][i][k]+=DT*((g*vdiag)-(f*vxxzz_T2))+(dthalbe*sumryy);
                                szz[j][i][k]+=DT*((g*vdiag)-(f*vxxyy_T2))+(dthalbe*sumrzz);
                                
                                /* Update the memory variables */
                                
                                /* in the case of L=1 update the memory-variables without suming*/
                                l=1;
                                n1=1.0/(1.0+(c1*eta[l]*0.5)); /* 1/(1+DT*c1/(2*tau_sigl)) */
                                n2=1.0-(c1*eta[l]*0.5);
                                n3=c2*eta[l]/2.0;
                                n4=c3*eta[l]/2.0;
                                dipjp=uipjp[j][i][k]*eta[l]*tausipjp[j][i][k];
                                djpkp=ujpkp[j][i][k]*eta[l]*tausjpkp[j][i][k];
                                dipkp=uipkp[j][i][k]*eta[l]*tausipkp[j][i][k];
                                d=2.0*u[j][i][k]*eta[l]*taus[j][i][k];
                                e=pi[j][i][k]*eta[l]*taup[j][i][k];
                                
                                
                                *(rxy_j_i_3+k)=n1*((*(rxy_j_i+k))*n2-n3*((*(rxy_j_i+k))+(*(rxy_j_i_2+k)))-n4*((*(rxy_j_i_2+k))+(*(rxy_j_i_3+k)))-(dipjp*vxyyx_T2));
                                *(ryz_j_i_3+k)=n1*((*(ryz_j_i+k))*n2-n3*((*(ryz_j_i+k))+(*(ryz_j_i_2+k)))-n4*((*(ryz_j_i_2+k))+(*(ryz_j_i_3+k)))-(djpkp*vyzzy_T2));
                                *(rxz_j_i_3+k)=n1*((*(rxz_j_i+k))*n2-n3*((*(rxz_j_i+k))+(*(rxz_j_i_2+k)))-n4*((*(rxz_j_i_2+k))+(*(rxz_j_i_3+k)))-(dipkp*vxzzx_T2));
                                *(rxx_j_i_3+k)=n1*((*(rxx_j_i+k))*n2-n3*((*(rxx_j_i+k))+(*(rxx_j_i_2+k)))-n4*((*(rxx_j_i_2+k))+(*(rxx_j_i_3+k)))-(e*vdiag)+(d*vyyzz_T2));
                                *(ryy_j_i_3+k)=n1*((*(ryy_j_i+k))*n2-n3*((*(ryy_j_i+k))+(*(ryy_j_i_2+k)))-n4*((*(ryy_j_i_2+k))+(*(ryy_j_i_3+k)))-(e*vdiag)+(d*vxxzz_T2));
                                *(rzz_j_i_3+k)=n1*((*(rzz_j_i+k))*n2-n3*((*(rzz_j_i+k))+(*(rzz_j_i_2+k)))-n4*((*(rzz_j_i_2+k))+(*(rzz_j_i_3+k)))-(e*vdiag)+(d*vxxyy_T2));
                                
                                
                                sumrxy=c1*(*(rxy_j_i_3+k))+c2*(*(rxy_j_i+k))+c3*(*(rxy_j_i_2+k));
                                sumryz=c1*(*(ryz_j_i_3+k))+c2*(*(ryz_j_i+k))+c3*(*(ryz_j_i_2+k));
                                sumrxz=c1*(*(rxz_j_i_3+k))+c2*(*(rxz_j_i+k))+c3*(*(rxz_j_i_2+k));
                                sumrxx=c1*(*(rxx_j_i_3+k))+c2*(*(rxx_j_i+k))+c3*(*(rxx_j_i_2+k));
                                sumryy=c1*(*(ryy_j_i_3+k))+c2*(*(ryy_j_i+k))+c3*(*(ryy_j_i_2+k));
                                sumrzz=c1*(*(rzz_j_i_3+k))+c2*(*(rzz_j_i+k))+c3*(*(rzz_j_i_2+k));
                                
                                
                                /* and now the components of the stress tensor are
                                 completely updated */
                                sxy[j][i][k]+=(dthalbe*sumrxy);
                                syz[j][i][k]+=(dthalbe*sumryz);
                                sxz[j][i][k]+=(dthalbe*sumrxz);
                                sxx[j][i][k]+=(dthalbe*sumrxx);
                                syy[j][i][k]+=(dthalbe*sumryy);
                                szz[j][i][k]+=(dthalbe*sumrzz);
                                
                                
                                
                            }
                        }
                    }
                    break;
                    
            }
            break; /* Break FDORDER_TIME=3 */
            
        case 4:
            c1=13.0/12.0; c2=-5.0/24.0; c3=1.0/6.0; c4=-1.0/24.0;
            switch (FDORDER){
                    
                case 2 :
                    
                    for (j=ny1;j<=ny2;j++){
                        
                        rxx_j=*(rxx+j);ryy_j=*(ryy+j);rzz_j=*(rzz+j);rxy_j=*(rxy+j);ryz_j=*(ryz+j);rxz_j=*(rxz+j);
                        rxx_j_2=*(rxx_2+j);ryy_j_2=*(ryy_2+j);rzz_j_2=*(rzz_2+j);rxy_j_2=*(rxy_2+j);ryz_j_2=*(ryz_2+j); rxz_j_2=*(rxz_2+j);
                        rxx_j_3=*(rxx_3+j);ryy_j_3=*(ryy_3+j);rzz_j_3=*(rzz_3+j);rxy_j_3=*(rxy_3+j);ryz_j_3=*(ryz_3+j);rxz_j_3=*(rxz_3+j);
                        rxx_j_4=*(rxx_4+j);ryy_j_4=*(ryy_4+j);rzz_j_4=*(rzz_4+j);rxy_j_4=*(rxy_4+j);ryz_j_4=*(ryz_4+j);rxz_j_4=*(rxz_4+j);
                        vxyyx_j=*(vxyyx+j);vyzzy_j=*(vyzzy+j);vxzzx_j=*(vxzzx+j);vxxyyzz_j=*(vxxyyzz+j);vyyzz_j=*(vyyzz+j);vxxzz_j=*(vxxzz+j);vxxyy_j=*(vxxyy+j);
                        vxyyx_j_2=*(vxyyx_2+j);vyzzy_j_2=*(vyzzy_2+j); vxzzx_j_2=*(vxzzx_2+j);vxxyyzz_j_2=*(vxxyyzz_2+j);vyyzz_j_2=*(vyyzz_2+j);vxxzz_j_2=*(vxxzz_2+j);vxxyy_j_2=*(vxxyy_2+j);
                        vxyyx_j_3=*(vxyyx_3+j);vyzzy_j_3=*(vyzzy_3+j);vxzzx_j_3=*(vxzzx_3+j);vxxyyzz_j_3=*(vxxyyzz_3+j);vyyzz_j_3=*(vyyzz_3+j);vxxzz_j_3=*(vxxzz_3+j);vxxyy_j_3=*(vxxyy_3+j);
                        vxyyx_j_4=*(vxyyx_4+j);vyzzy_j_4=*(vyzzy_4+j);vxzzx_j_4=*(vxzzx_4+j);vxxyyzz_j_4=*(vxxyyzz_4+j);vyyzz_j_4=*(vyyzz_4+j);vxxzz_j_4=*(vxxzz_4+j);vxxyy_j_4=*(vxxyy_4+j);

                        for (i=nx1;i<=nx2;i++){
                            
                            rxx_j_i=*(rxx_j+i);ryy_j_i=*(ryy_j+i);rzz_j_i=*(rzz_j+i);rxy_j_i=*(rxy_j+i);ryz_j_i=*(ryz_j+i);rxz_j_i=*(rxz_j+i);
                            rxx_j_i_2=*(rxx_j_2+i);ryy_j_i_2=*(ryy_j_2+i);rzz_j_i_2=*(rzz_j_2+i);rxy_j_i_2=*(rxy_j_2+i);ryz_j_i_2=*(ryz_j_2+i);rxz_j_i_2=*(rxz_j_2+i);
                            rxx_j_i_3=*(rxx_j_3+i);ryy_j_i_3=*(ryy_j_3+i);rzz_j_i_3=*(rzz_j_3+i);rxy_j_i_3=*(rxy_j_3+i);ryz_j_i_3=*(ryz_j_3+i);rxz_j_i_3=*(rxz_j_3+i);
                            rxx_j_i_4=*(rxx_j_4+i);ryy_j_i_4=*(ryy_j_4+i);rzz_j_i_4=*(rzz_j_4+i);rxy_j_i_4=*(rxy_j_4+i);ryz_j_i_4=*(ryz_j_4+i);rxz_j_i_4=*(rxz_j_4+i);
                            vxyyx_j_i=*(vxyyx_j+i);vyzzy_j_i=*(vyzzy_j+i);vxzzx_j_i=*(vxzzx_j+i);vxxyyzz_j_i=*(vxxyyzz_j+i);vyyzz_j_i=*(vyyzz_j+i);vxxzz_j_i=*(vxxzz_j+i);vxxyy_j_i=*(vxxyy_j+i);
                            vxyyx_j_i_2=*(vxyyx_j_2+i);vyzzy_j_i_2=*(vyzzy_j_2+i);vxzzx_j_i_2=*(vxzzx_j_2+i);vxxyyzz_j_i_2=*(vxxyyzz_j_2+i);vyyzz_j_i_2=*(vyyzz_j_2+i);vxxzz_j_i_2=*(vxxzz_j_2+i);vxxyy_j_i_2=*(vxxyy_j_2+i);
                            vxyyx_j_i_3=*(vxyyx_j_3+i);vyzzy_j_i_3=*(vyzzy_j_3+i);vxzzx_j_i_3=*(vxzzx_j_3+i); vxxyyzz_j_i_3=*(vxxyyzz_j_3+i);vyyzz_j_i_3=*(vyyzz_j_3+i);vxxzz_j_i_3=*(vxxzz_j_3+i); vxxyy_j_i_3=*(vxxyy_j_3+i);
                            vxyyx_j_i_4=*(vxyyx_j_4+i);vyzzy_j_i_4=*(vyzzy_j_4+i);vxzzx_j_i_4=*(vxzzx_j_4+i); vxxyyzz_j_i_4=*(vxxyyzz_j_4+i);vyyzz_j_i_4=*(vyyzz_j_4+i);vxxzz_j_i_4=*(vxxzz_j_4+i);vxxyy_j_i_4=*(vxxyy_j_4+i);
                            
                            for (k=nz1;k<=nz2;k++){
                                
                                /* spatial derivatives of the components of the velocities
                                 are computed */
                                
                                vxx = (vx[j][i][k]-vx[j][i-1][k])/DX;
                                vxy = (vx[j+1][i][k]-vx[j][i][k])/DY;
                                vxz = (vx[j][i][k+1]-vx[j][i][k])/DZ;
                                vyx = (vy[j][i+1][k]-vy[j][i][k])/DX;
                                vyy = (vy[j][i][k]-vy[j-1][i][k])/DY;
                                vyz = (vy[j][i][k+1]-vy[j][i][k])/DZ;
                                vzx = (vz[j][i+1][k]-vz[j][i][k])/DX;
                                vzy = (vz[j+1][i][k]-vz[j][i][k])/DY;
                                vzz = (vz[j][i][k]-vz[j][i][k-1])/DZ;
                                
                                /* computing sums of the old memory variables in this version only one mechanism is possi
                                 ble */
                                
                                sumrxy=c1*(*(rxy_j_i+k))+c2*(*(rxy_j_i_2+k))+c3*(*(rxy_j_i_3+k))+c4*(*(rxy_j_i_4+k));
                                sumryz=c1*(*(ryz_j_i+k))+c2*(*(ryz_j_i_2+k))+c3*(*(ryz_j_i_3+k))+c4*(*(ryz_j_i_4+k));
                                sumrxz=c1*(*(rxz_j_i+k))+c2*(*(rxz_j_i_2+k))+c3*(*(rxz_j_i_3+k))+c4*(*(rxz_j_i_4+k));
                                sumrxx=c1*(*(rxx_j_i+k))+c2*(*(rxx_j_i_2+k))+c3*(*(rxx_j_i_3+k))+c4*(*(rxx_j_i_4+k));
                                sumryy=c1*(*(ryy_j_i+k))+c2*(*(ryy_j_i_2+k))+c3*(*(ryy_j_i_3+k))+c4*(*(ryy_j_i_4+k));
                                sumrzz=c1*(*(rzz_j_i+k))+c2*(*(rzz_j_i_2+k))+c3*(*(rzz_j_i_3+k))+c4*(*(rzz_j_i_4+k));
                                
                                /* updating components of the stress tensor, partially */
                                fipjp=uipjp[j][i][k]*DT*(1.0+L*tausipjp[j][i][k]);
                                fjpkp=ujpkp[j][i][k]*DT*(1.0+L*tausjpkp[j][i][k]);
                                fipkp=uipkp[j][i][k]*DT*(1.0+L*tausipkp[j][i][k]);
                                g=pi[j][i][k]*(1.0+L*taup[j][i][k]);
                                f=2.0*u[j][i][k]*(1.0+L*taus[j][i][k]);
                                
                                /* Save spatial derivations */
                                *(vxyyx_j_i+k)=vxy+vyx;
                                *(vyzzy_j_i+k)=vyz+vzy;
                                *(vxzzx_j_i+k)=vxz+vzx;
                                *(vxxyyzz_j_i+k)=vxx+vyy+vzz;
                                *(vyyzz_j_i+k)=vyy+vzz;
                                *(vxxzz_j_i+k)=vxx+vzz;
                                *(vxxyy_j_i+k)=vxx+vyy;
                                
                                /* Do the Adams Bashforth stuff and save it, we need it for the memory variables again */
                                vxyyx_T2=(c1*(*(vxyyx_j_i+k))+c2*(*(vxyyx_j_i_2+k))+c3*(*(vxyyx_j_i_3+k))+c4*(*(vxyyx_j_i_4+k)));
                                vyzzy_T2=(c1*(*(vyzzy_j_i+k))+c2*(*(vyzzy_j_i_2+k))+c3*(*(vyzzy_j_i_3+k))+c4*(*(vyzzy_j_i_4+k)));
                                vxzzx_T2=(c1*(*(vxzzx_j_i+k))+c2*(*(vxzzx_j_i_2+k))+c3*(*(vxzzx_j_i_3+k))+c4*(*(vxzzx_j_i_4+k)));
                                vyyzz_T2=(c1*(*(vyyzz_j_i+k))+c2*(*(vyyzz_j_i_2+k))+c3*(*(vyyzz_j_i_3+k))+c4*(*(vyyzz_j_i_4+k)));
                                vxxzz_T2=(c1*(*(vxxzz_j_i+k))+c2*(*(vxxzz_j_i_2+k))+c3*(*(vxxzz_j_i_3+k))+c4*(*(vxxzz_j_i_4+k)));
                                vxxyy_T2=(c1*(*(vxxyy_j_i+k))+c2*(*(vxxyy_j_i_2+k))+c3*(*(vxxyy_j_i_3+k))+c4*(*(vxxyy_j_i_4+k)));
                                vdiag=((c1*(*(vxxyyzz_j_i+k))+c2*(*(vxxyyzz_j_i_2+k))+c3*(*(vxxyyzz_j_i_3+k)))+c4*(*(vxxyyzz_j_i_4+k)));
                                
                                /* Update the stress tensor */
                                sxy[j][i][k]+=(fipjp*vxyyx_T2)+(dthalbe*sumrxy);
                                syz[j][i][k]+=(fjpkp*vyzzy_T2)+(dthalbe*sumryz);
                                sxz[j][i][k]+=(fipkp*vxzzx_T2)+(dthalbe*sumrxz);
                                sxx[j][i][k]+=DT*((g*vdiag)-(f*vyyzz_T2))+(dthalbe*sumrxx);
                                syy[j][i][k]+=DT*((g*vdiag)-(f*vxxzz_T2))+(dthalbe*sumryy);
                                szz[j][i][k]+=DT*((g*vdiag)-(f*vxxyy_T2))+(dthalbe*sumrzz);
                                
                                /* Update the memory variables */
                                
                                /* in the case of L=1 update the memory-variables without suming*/
                                l=1;
                                n1=1.0/(1.0+(c1*eta[l]*0.5)); /* 1/(1+DT*c1/(2*tau_sigl)) */
                                n2=1.0-(c1*eta[l]*0.5);
                                n3=c2*eta[l]/2.0;
                                n4=c3*eta[l]/2.0;
                                n5=c4*eta[l]/2.0;
                                dipjp=uipjp[j][i][k]*eta[l]*tausipjp[j][i][k];
                                djpkp=ujpkp[j][i][k]*eta[l]*tausjpkp[j][i][k];
                                dipkp=uipkp[j][i][k]*eta[l]*tausipkp[j][i][k];
                                d=2.0*u[j][i][k]*eta[l]*taus[j][i][k];
                                e=pi[j][i][k]*eta[l]*taup[j][i][k];
                                
                                
                                *(rxy_j_i_4+k)=n1*((*(rxy_j_i+k))*n2-n3*((*(rxy_j_i+k))+(*(rxy_j_i_2+k)))-n4*((*(rxy_j_i_2+k))+(*(rxy_j_i_3+k)))-n5*((*(rxy_j_i_3+k))+(*(rxy_j_i_4+k)))-(dipjp*vxyyx_T2));
                                *(ryz_j_i_4+k)=n1*((*(ryz_j_i+k))*n2-n3*((*(ryz_j_i+k))+(*(ryz_j_i_2+k)))-n4*((*(ryz_j_i_2+k))+(*(ryz_j_i_3+k)))-n5*((*(ryz_j_i_3+k))+(*(ryz_j_i_4+k)))-(djpkp*vyzzy_T2));
                                *(rxz_j_i_4+k)=n1*((*(rxz_j_i+k))*n2-n3*((*(rxz_j_i+k))+(*(rxz_j_i_2+k)))-n4*((*(rxz_j_i_2+k))+(*(rxz_j_i_3+k)))-n5*((*(rxz_j_i_3+k))+(*(rxz_j_i_4+k)))-(dipkp*vxzzx_T2));
                                *(rxx_j_i_4+k)=n1*((*(rxx_j_i+k))*n2-n3*((*(rxx_j_i+k))+(*(rxx_j_i_2+k)))-n4*((*(rxx_j_i_2+k))+(*(rxx_j_i_3+k)))-n5*((*(rxx_j_i_3+k))+(*(rxx_j_i_4+k)))-(e*vdiag)+(d*vyyzz_T2));
                                *(ryy_j_i_4+k)=n1*((*(ryy_j_i+k))*n2-n3*((*(ryy_j_i+k))+(*(ryy_j_i_2+k)))-n4*((*(ryy_j_i_2+k))+(*(ryy_j_i_3+k)))-n5*((*(ryy_j_i_3+k))+(*(ryy_j_i_4+k)))-(e*vdiag)+(d*vxxzz_T2));
                                *(rzz_j_i_4+k)=n1*((*(rzz_j_i+k))*n2-n3*((*(rzz_j_i+k))+(*(rzz_j_i_2+k)))-n4*((*(rzz_j_i_2+k))+(*(rzz_j_i_3+k)))-n5*((*(rzz_j_i_3+k))+(*(rzz_j_i_4+k)))-(e*vdiag)+(d*vxxyy_T2));
                                
                                
                                sumrxy=c1*(*(rxy_j_i_4+k))+c2*(*(rxy_j_i+k))+c3*(*(rxy_j_i_2+k))+c4*(*(rxy_j_i_3+k));
                                sumryz=c1*(*(ryz_j_i_4+k))+c2*(*(ryz_j_i+k))+c3*(*(ryz_j_i_2+k))+c4*(*(ryz_j_i_3+k));
                                sumrxz=c1*(*(rxz_j_i_4+k))+c2*(*(rxz_j_i+k))+c3*(*(rxz_j_i_2+k))+c4*(*(rxz_j_i_3+k));
                                sumrxx=c1*(*(rxx_j_i_4+k))+c2*(*(rxx_j_i+k))+c3*(*(rxx_j_i_2+k))+c4*(*(rxx_j_i_3+k));
                                sumryy=c1*(*(ryy_j_i_4+k))+c2*(*(ryy_j_i+k))+c3*(*(ryy_j_i_2+k))+c4*(*(ryy_j_i_3+k));
                                sumrzz=c1*(*(rzz_j_i_4+k))+c2*(*(rzz_j_i+k))+c3*(*(rzz_j_i_2+k))+c4*(*(rzz_j_i_3+k));
                                
                                
                                /* and now the components of the stress tensor are
                                 completely updated */
                                sxy[j][i][k]+=(dthalbe*sumrxy);
                                syz[j][i][k]+=(dthalbe*sumryz);
                                sxz[j][i][k]+=(dthalbe*sumrxz);
                                sxx[j][i][k]+=(dthalbe*sumrxx);
                                syy[j][i][k]+=(dthalbe*sumryy);
                                szz[j][i][k]+=(dthalbe*sumrzz);
                                
                                
                            }
                        }
                    }
                    break;
                    
                case 4 :
                    
                    b1=9.0/8.0; b2=-1.0/24.0; /* Taylor coefficients */
                    if(FDCOEFF==2){
                        b1=1.1382; b2=-0.046414;} /* Holberg coefficients E=0.1 %*/
                    
                    for (j=ny1;j<=ny2;j++){
                        
                        rxx_j=*(rxx+j);ryy_j=*(ryy+j);rzz_j=*(rzz+j);rxy_j=*(rxy+j);ryz_j=*(ryz+j);rxz_j=*(rxz+j);
                        rxx_j_2=*(rxx_2+j);ryy_j_2=*(ryy_2+j);rzz_j_2=*(rzz_2+j);rxy_j_2=*(rxy_2+j);ryz_j_2=*(ryz_2+j); rxz_j_2=*(rxz_2+j);
                        rxx_j_3=*(rxx_3+j);ryy_j_3=*(ryy_3+j);rzz_j_3=*(rzz_3+j);rxy_j_3=*(rxy_3+j);ryz_j_3=*(ryz_3+j);rxz_j_3=*(rxz_3+j);
                        rxx_j_4=*(rxx_4+j);ryy_j_4=*(ryy_4+j);rzz_j_4=*(rzz_4+j);rxy_j_4=*(rxy_4+j);ryz_j_4=*(ryz_4+j);rxz_j_4=*(rxz_4+j);
                        vxyyx_j=*(vxyyx+j);vyzzy_j=*(vyzzy+j);vxzzx_j=*(vxzzx+j);vxxyyzz_j=*(vxxyyzz+j);vyyzz_j=*(vyyzz+j);vxxzz_j=*(vxxzz+j);vxxyy_j=*(vxxyy+j);
                        vxyyx_j_2=*(vxyyx_2+j);vyzzy_j_2=*(vyzzy_2+j); vxzzx_j_2=*(vxzzx_2+j);vxxyyzz_j_2=*(vxxyyzz_2+j);vyyzz_j_2=*(vyyzz_2+j);vxxzz_j_2=*(vxxzz_2+j);vxxyy_j_2=*(vxxyy_2+j);
                        vxyyx_j_3=*(vxyyx_3+j);vyzzy_j_3=*(vyzzy_3+j);vxzzx_j_3=*(vxzzx_3+j);vxxyyzz_j_3=*(vxxyyzz_3+j);vyyzz_j_3=*(vyyzz_3+j);vxxzz_j_3=*(vxxzz_3+j);vxxyy_j_3=*(vxxyy_3+j);
                        vxyyx_j_4=*(vxyyx_4+j);vyzzy_j_4=*(vyzzy_4+j);vxzzx_j_4=*(vxzzx_4+j);vxxyyzz_j_4=*(vxxyyzz_4+j);vyyzz_j_4=*(vyyzz_4+j);vxxzz_j_4=*(vxxzz_4+j);vxxyy_j_4=*(vxxyy_4+j);
                        
                        for (i=nx1;i<=nx2;i++){
                            
                            rxx_j_i=*(rxx_j+i);ryy_j_i=*(ryy_j+i);rzz_j_i=*(rzz_j+i);rxy_j_i=*(rxy_j+i);ryz_j_i=*(ryz_j+i);rxz_j_i=*(rxz_j+i);
                            rxx_j_i_2=*(rxx_j_2+i);ryy_j_i_2=*(ryy_j_2+i);rzz_j_i_2=*(rzz_j_2+i);rxy_j_i_2=*(rxy_j_2+i);ryz_j_i_2=*(ryz_j_2+i);rxz_j_i_2=*(rxz_j_2+i);
                            rxx_j_i_3=*(rxx_j_3+i);ryy_j_i_3=*(ryy_j_3+i);rzz_j_i_3=*(rzz_j_3+i);rxy_j_i_3=*(rxy_j_3+i);ryz_j_i_3=*(ryz_j_3+i);rxz_j_i_3=*(rxz_j_3+i);
                            rxx_j_i_4=*(rxx_j_4+i);ryy_j_i_4=*(ryy_j_4+i);rzz_j_i_4=*(rzz_j_4+i);rxy_j_i_4=*(rxy_j_4+i);ryz_j_i_4=*(ryz_j_4+i);rxz_j_i_4=*(rxz_j_4+i);
                            vxyyx_j_i=*(vxyyx_j+i);vyzzy_j_i=*(vyzzy_j+i);vxzzx_j_i=*(vxzzx_j+i);vxxyyzz_j_i=*(vxxyyzz_j+i);vyyzz_j_i=*(vyyzz_j+i);vxxzz_j_i=*(vxxzz_j+i);vxxyy_j_i=*(vxxyy_j+i);
                            vxyyx_j_i_2=*(vxyyx_j_2+i);vyzzy_j_i_2=*(vyzzy_j_2+i);vxzzx_j_i_2=*(vxzzx_j_2+i);vxxyyzz_j_i_2=*(vxxyyzz_j_2+i);vyyzz_j_i_2=*(vyyzz_j_2+i);vxxzz_j_i_2=*(vxxzz_j_2+i);vxxyy_j_i_2=*(vxxyy_j_2+i);
                            vxyyx_j_i_3=*(vxyyx_j_3+i);vyzzy_j_i_3=*(vyzzy_j_3+i);vxzzx_j_i_3=*(vxzzx_j_3+i); vxxyyzz_j_i_3=*(vxxyyzz_j_3+i);vyyzz_j_i_3=*(vyyzz_j_3+i);vxxzz_j_i_3=*(vxxzz_j_3+i); vxxyy_j_i_3=*(vxxyy_j_3+i);
                            vxyyx_j_i_4=*(vxyyx_j_4+i);vyzzy_j_i_4=*(vyzzy_j_4+i);vxzzx_j_i_4=*(vxzzx_j_4+i); vxxyyzz_j_i_4=*(vxxyyzz_j_4+i);vyyzz_j_i_4=*(vyyzz_j_4+i);vxxzz_j_i_4=*(vxxzz_j_4+i);vxxyy_j_i_4=*(vxxyy_j_4+i);
                            
                            for (k=nz1;k<=nz2;k++){
                                /* spatial derivatives of the components of the velocities
                                 are computed */
                                
                                vxx = (b1*(vx[j][i][k]-vx[j][i-1][k])+b2*(vx[j][i+1][k]-vx[j][i-2][k]))/DX;
                                vxy = (b1*(vx[j+1][i][k]-vx[j][i][k])+b2*(vx[j+2][i][k]-vx[j-1][i][k]))/DY;
                                vxz = (b1*(vx[j][i][k+1]-vx[j][i][k])+b2*(vx[j][i][k+2]-vx[j][i][k-1]))/DZ;
                                vyx = (b1*(vy[j][i+1][k]-vy[j][i][k])+b2*(vy[j][i+2][k]-vy[j][i-1][k]))/DX;
                                vyy = (b1*(vy[j][i][k]-vy[j-1][i][k])+b2*(vy[j+1][i][k]-vy[j-2][i][k]))/DY;
                                vyz = (b1*(vy[j][i][k+1]-vy[j][i][k])+b2*(vy[j][i][k+2]-vy[j][i][k-1]))/DZ;
                                vzx = (b1*(vz[j][i+1][k]-vz[j][i][k])+b2*(vz[j][i+2][k]-vz[j][i-1][k]))/DX;
                                vzy = (b1*(vz[j+1][i][k]-vz[j][i][k])+b2*(vz[j+2][i][k]-vz[j-1][i][k]))/DY;
                                vzz = (b1*(vz[j][i][k]-vz[j][i][k-1])+b2*(vz[j][i][k+1]-vz[j][i][k-2]))/DZ;
                                
                                
                                /* computing sums of the old memory variables in this version only one mechanism is possi
                                 ble */
                                
                                sumrxy=c1*(*(rxy_j_i+k))+c2*(*(rxy_j_i_2+k))+c3*(*(rxy_j_i_3+k))+c4*(*(rxy_j_i_4+k));
                                sumryz=c1*(*(ryz_j_i+k))+c2*(*(ryz_j_i_2+k))+c3*(*(ryz_j_i_3+k))+c4*(*(ryz_j_i_4+k));
                                sumrxz=c1*(*(rxz_j_i+k))+c2*(*(rxz_j_i_2+k))+c3*(*(rxz_j_i_3+k))+c4*(*(rxz_j_i_4+k));
                                sumrxx=c1*(*(rxx_j_i+k))+c2*(*(rxx_j_i_2+k))+c3*(*(rxx_j_i_3+k))+c4*(*(rxx_j_i_4+k));
                                sumryy=c1*(*(ryy_j_i+k))+c2*(*(ryy_j_i_2+k))+c3*(*(ryy_j_i_3+k))+c4*(*(ryy_j_i_4+k));
                                sumrzz=c1*(*(rzz_j_i+k))+c2*(*(rzz_j_i_2+k))+c3*(*(rzz_j_i_3+k))+c4*(*(rzz_j_i_4+k));
                                
                                /* updating components of the stress tensor, partially */
                                fipjp=uipjp[j][i][k]*DT*(1.0+L*tausipjp[j][i][k]);
                                fjpkp=ujpkp[j][i][k]*DT*(1.0+L*tausjpkp[j][i][k]);
                                fipkp=uipkp[j][i][k]*DT*(1.0+L*tausipkp[j][i][k]);
                                g=pi[j][i][k]*(1.0+L*taup[j][i][k]);
                                f=2.0*u[j][i][k]*(1.0+L*taus[j][i][k]);
                                
                                /* Save spatial derivations */
                                *(vxyyx_j_i+k)=vxy+vyx;
                                *(vyzzy_j_i+k)=vyz+vzy;
                                *(vxzzx_j_i+k)=vxz+vzx;
                                *(vxxyyzz_j_i+k)=vxx+vyy+vzz;
                                *(vyyzz_j_i+k)=vyy+vzz;
                                *(vxxzz_j_i+k)=vxx+vzz;
                                *(vxxyy_j_i+k)=vxx+vyy;
                                
                                /* Do the Adams Bashforth stuff and save it, we need it for the memory variables again */
                                vxyyx_T2=(c1*(*(vxyyx_j_i+k))+c2*(*(vxyyx_j_i_2+k))+c3*(*(vxyyx_j_i_3+k))+c4*(*(vxyyx_j_i_4+k)));
                                vyzzy_T2=(c1*(*(vyzzy_j_i+k))+c2*(*(vyzzy_j_i_2+k))+c3*(*(vyzzy_j_i_3+k))+c4*(*(vyzzy_j_i_4+k)));
                                vxzzx_T2=(c1*(*(vxzzx_j_i+k))+c2*(*(vxzzx_j_i_2+k))+c3*(*(vxzzx_j_i_3+k))+c4*(*(vxzzx_j_i_4+k)));
                                vyyzz_T2=(c1*(*(vyyzz_j_i+k))+c2*(*(vyyzz_j_i_2+k))+c3*(*(vyyzz_j_i_3+k))+c4*(*(vyyzz_j_i_4+k)));
                                vxxzz_T2=(c1*(*(vxxzz_j_i+k))+c2*(*(vxxzz_j_i_2+k))+c3*(*(vxxzz_j_i_3+k))+c4*(*(vxxzz_j_i_4+k)));
                                vxxyy_T2=(c1*(*(vxxyy_j_i+k))+c2*(*(vxxyy_j_i_2+k))+c3*(*(vxxyy_j_i_3+k))+c4*(*(vxxyy_j_i_4+k)));
                                vdiag=((c1*(*(vxxyyzz_j_i+k))+c2*(*(vxxyyzz_j_i_2+k))+c3*(*(vxxyyzz_j_i_3+k)))+c4*(*(vxxyyzz_j_i_4+k)));
                                
                                /* Update the stress tensor */
                                sxy[j][i][k]+=(fipjp*vxyyx_T2)+(dthalbe*sumrxy);
                                syz[j][i][k]+=(fjpkp*vyzzy_T2)+(dthalbe*sumryz);
                                sxz[j][i][k]+=(fipkp*vxzzx_T2)+(dthalbe*sumrxz);
                                sxx[j][i][k]+=DT*((g*vdiag)-(f*vyyzz_T2))+(dthalbe*sumrxx);
                                syy[j][i][k]+=DT*((g*vdiag)-(f*vxxzz_T2))+(dthalbe*sumryy);
                                szz[j][i][k]+=DT*((g*vdiag)-(f*vxxyy_T2))+(dthalbe*sumrzz);
                                
                                /* Update the memory variables */
                                
                                /* in the case of L=1 update the memory-variables without suming*/
                                l=1;
                                n1=1.0/(1.0+(c1*eta[l]*0.5)); /* 1/(1+DT*c1/(2*tau_sigl)) */
                                n2=1.0-(c1*eta[l]*0.5);
                                n3=c2*eta[l]/2.0;
                                n4=c3*eta[l]/2.0;
                                n5=c4*eta[l]/2.0;
                                dipjp=uipjp[j][i][k]*eta[l]*tausipjp[j][i][k];
                                djpkp=ujpkp[j][i][k]*eta[l]*tausjpkp[j][i][k];
                                dipkp=uipkp[j][i][k]*eta[l]*tausipkp[j][i][k];
                                d=2.0*u[j][i][k]*eta[l]*taus[j][i][k];
                                e=pi[j][i][k]*eta[l]*taup[j][i][k];
                                
                                
                                *(rxy_j_i_4+k)=n1*((*(rxy_j_i+k))*n2-n3*((*(rxy_j_i+k))+(*(rxy_j_i_2+k)))-n4*((*(rxy_j_i_2+k))+(*(rxy_j_i_3+k)))-n5*((*(rxy_j_i_3+k))+(*(rxy_j_i_4+k)))-(dipjp*vxyyx_T2));
                                *(ryz_j_i_4+k)=n1*((*(ryz_j_i+k))*n2-n3*((*(ryz_j_i+k))+(*(ryz_j_i_2+k)))-n4*((*(ryz_j_i_2+k))+(*(ryz_j_i_3+k)))-n5*((*(ryz_j_i_3+k))+(*(ryz_j_i_4+k)))-(djpkp*vyzzy_T2));
                                *(rxz_j_i_4+k)=n1*((*(rxz_j_i+k))*n2-n3*((*(rxz_j_i+k))+(*(rxz_j_i_2+k)))-n4*((*(rxz_j_i_2+k))+(*(rxz_j_i_3+k)))-n5*((*(rxz_j_i_3+k))+(*(rxz_j_i_4+k)))-(dipkp*vxzzx_T2));
                                *(rxx_j_i_4+k)=n1*((*(rxx_j_i+k))*n2-n3*((*(rxx_j_i+k))+(*(rxx_j_i_2+k)))-n4*((*(rxx_j_i_2+k))+(*(rxx_j_i_3+k)))-n5*((*(rxx_j_i_3+k))+(*(rxx_j_i_4+k)))-(e*vdiag)+(d*vyyzz_T2));
                                *(ryy_j_i_4+k)=n1*((*(ryy_j_i+k))*n2-n3*((*(ryy_j_i+k))+(*(ryy_j_i_2+k)))-n4*((*(ryy_j_i_2+k))+(*(ryy_j_i_3+k)))-n5*((*(ryy_j_i_3+k))+(*(ryy_j_i_4+k)))-(e*vdiag)+(d*vxxzz_T2));
                                *(rzz_j_i_4+k)=n1*((*(rzz_j_i+k))*n2-n3*((*(rzz_j_i+k))+(*(rzz_j_i_2+k)))-n4*((*(rzz_j_i_2+k))+(*(rzz_j_i_3+k)))-n5*((*(rzz_j_i_3+k))+(*(rzz_j_i_4+k)))-(e*vdiag)+(d*vxxyy_T2));
                                
                                
                                sumrxy=c1*(*(rxy_j_i_4+k))+c2*(*(rxy_j_i+k))+c3*(*(rxy_j_i_2+k))+c4*(*(rxy_j_i_3+k));
                                sumryz=c1*(*(ryz_j_i_4+k))+c2*(*(ryz_j_i+k))+c3*(*(ryz_j_i_2+k))+c4*(*(ryz_j_i_3+k));
                                sumrxz=c1*(*(rxz_j_i_4+k))+c2*(*(rxz_j_i+k))+c3*(*(rxz_j_i_2+k))+c4*(*(rxz_j_i_3+k));
                                sumrxx=c1*(*(rxx_j_i_4+k))+c2*(*(rxx_j_i+k))+c3*(*(rxx_j_i_2+k))+c4*(*(rxx_j_i_3+k));
                                sumryy=c1*(*(ryy_j_i_4+k))+c2*(*(ryy_j_i+k))+c3*(*(ryy_j_i_2+k))+c4*(*(ryy_j_i_3+k));
                                sumrzz=c1*(*(rzz_j_i_4+k))+c2*(*(rzz_j_i+k))+c3*(*(rzz_j_i_2+k))+c4*(*(rzz_j_i_3+k));
                                
                                
                                /* and now the components of the stress tensor are
                                 completely updated */
                                sxy[j][i][k]+=(dthalbe*sumrxy);
                                syz[j][i][k]+=(dthalbe*sumryz);
                                sxz[j][i][k]+=(dthalbe*sumrxz);
                                sxx[j][i][k]+=(dthalbe*sumrxx);
                                syy[j][i][k]+=(dthalbe*sumryy);
                                szz[j][i][k]+=(dthalbe*sumrzz);
                                
                                

                                
                            }
                        }
                    }
                    break;
                    
                case 6 :
                    
                    b1=75.0/64.0; b2=-25.0/384.0; b3=3.0/640.0; /* Taylor coefficients */
                    if(FDCOEFF==2){
                        b1=1.1965; b2=-0.078804; b3=0.0081781;}   /* Holberg coefficients E=0.1 %*/
                    
                    
                    for (j=ny1;j<=ny2;j++){
                        
                        rxx_j=*(rxx+j);ryy_j=*(ryy+j);rzz_j=*(rzz+j);rxy_j=*(rxy+j);ryz_j=*(ryz+j);rxz_j=*(rxz+j);
                        rxx_j_2=*(rxx_2+j);ryy_j_2=*(ryy_2+j);rzz_j_2=*(rzz_2+j);rxy_j_2=*(rxy_2+j);ryz_j_2=*(ryz_2+j); rxz_j_2=*(rxz_2+j);
                        rxx_j_3=*(rxx_3+j);ryy_j_3=*(ryy_3+j);rzz_j_3=*(rzz_3+j);rxy_j_3=*(rxy_3+j);ryz_j_3=*(ryz_3+j);rxz_j_3=*(rxz_3+j);
                        rxx_j_4=*(rxx_4+j);ryy_j_4=*(ryy_4+j);rzz_j_4=*(rzz_4+j);rxy_j_4=*(rxy_4+j);ryz_j_4=*(ryz_4+j);rxz_j_4=*(rxz_4+j);
                        vxyyx_j=*(vxyyx+j);vyzzy_j=*(vyzzy+j);vxzzx_j=*(vxzzx+j);vxxyyzz_j=*(vxxyyzz+j);vyyzz_j=*(vyyzz+j);vxxzz_j=*(vxxzz+j);vxxyy_j=*(vxxyy+j);
                        vxyyx_j_2=*(vxyyx_2+j);vyzzy_j_2=*(vyzzy_2+j); vxzzx_j_2=*(vxzzx_2+j);vxxyyzz_j_2=*(vxxyyzz_2+j);vyyzz_j_2=*(vyyzz_2+j);vxxzz_j_2=*(vxxzz_2+j);vxxyy_j_2=*(vxxyy_2+j);
                        vxyyx_j_3=*(vxyyx_3+j);vyzzy_j_3=*(vyzzy_3+j);vxzzx_j_3=*(vxzzx_3+j);vxxyyzz_j_3=*(vxxyyzz_3+j);vyyzz_j_3=*(vyyzz_3+j);vxxzz_j_3=*(vxxzz_3+j);vxxyy_j_3=*(vxxyy_3+j);
                        vxyyx_j_4=*(vxyyx_4+j);vyzzy_j_4=*(vyzzy_4+j);vxzzx_j_4=*(vxzzx_4+j);vxxyyzz_j_4=*(vxxyyzz_4+j);vyyzz_j_4=*(vyyzz_4+j);vxxzz_j_4=*(vxxzz_4+j);vxxyy_j_4=*(vxxyy_4+j);
                        
                        for (i=nx1;i<=nx2;i++){
                            
                            rxx_j_i=*(rxx_j+i);ryy_j_i=*(ryy_j+i);rzz_j_i=*(rzz_j+i);rxy_j_i=*(rxy_j+i);ryz_j_i=*(ryz_j+i);rxz_j_i=*(rxz_j+i);
                            rxx_j_i_2=*(rxx_j_2+i);ryy_j_i_2=*(ryy_j_2+i);rzz_j_i_2=*(rzz_j_2+i);rxy_j_i_2=*(rxy_j_2+i);ryz_j_i_2=*(ryz_j_2+i);rxz_j_i_2=*(rxz_j_2+i);
                            rxx_j_i_3=*(rxx_j_3+i);ryy_j_i_3=*(ryy_j_3+i);rzz_j_i_3=*(rzz_j_3+i);rxy_j_i_3=*(rxy_j_3+i);ryz_j_i_3=*(ryz_j_3+i);rxz_j_i_3=*(rxz_j_3+i);
                            rxx_j_i_4=*(rxx_j_4+i);ryy_j_i_4=*(ryy_j_4+i);rzz_j_i_4=*(rzz_j_4+i);rxy_j_i_4=*(rxy_j_4+i);ryz_j_i_4=*(ryz_j_4+i);rxz_j_i_4=*(rxz_j_4+i);
                            vxyyx_j_i=*(vxyyx_j+i);vyzzy_j_i=*(vyzzy_j+i);vxzzx_j_i=*(vxzzx_j+i);vxxyyzz_j_i=*(vxxyyzz_j+i);vyyzz_j_i=*(vyyzz_j+i);vxxzz_j_i=*(vxxzz_j+i);vxxyy_j_i=*(vxxyy_j+i);
                            vxyyx_j_i_2=*(vxyyx_j_2+i);vyzzy_j_i_2=*(vyzzy_j_2+i);vxzzx_j_i_2=*(vxzzx_j_2+i);vxxyyzz_j_i_2=*(vxxyyzz_j_2+i);vyyzz_j_i_2=*(vyyzz_j_2+i);vxxzz_j_i_2=*(vxxzz_j_2+i);vxxyy_j_i_2=*(vxxyy_j_2+i);
                            vxyyx_j_i_3=*(vxyyx_j_3+i);vyzzy_j_i_3=*(vyzzy_j_3+i);vxzzx_j_i_3=*(vxzzx_j_3+i); vxxyyzz_j_i_3=*(vxxyyzz_j_3+i);vyyzz_j_i_3=*(vyyzz_j_3+i);vxxzz_j_i_3=*(vxxzz_j_3+i); vxxyy_j_i_3=*(vxxyy_j_3+i);
                            vxyyx_j_i_4=*(vxyyx_j_4+i);vyzzy_j_i_4=*(vyzzy_j_4+i);vxzzx_j_i_4=*(vxzzx_j_4+i); vxxyyzz_j_i_4=*(vxxyyzz_j_4+i);vyyzz_j_i_4=*(vyyzz_j_4+i);vxxzz_j_i_4=*(vxxzz_j_4+i);vxxyy_j_i_4=*(vxxyy_j_4+i);
                            
                            for (k=nz1;k<=nz2;k++){                                /* spatial derivatives of the components of the velocities
                                 are computed */
                                
                                vxx = (b1*(vx[j][i][k]-vx[j][i-1][k])+
                                       b2*(vx[j][i+1][k]-vx[j][i-2][k])+
                                       b3*(vx[j][i+2][k]-vx[j][i-3][k]))/DX;
                                
                                vxy = (b1*(vx[j+1][i][k]-vx[j][i][k])+
                                       b2*(vx[j+2][i][k]-vx[j-1][i][k])+
                                       b3*(vx[j+3][i][k]-vx[j-2][i][k]))/DY;
                                
                                vxz = (b1*(vx[j][i][k+1]-vx[j][i][k])+
                                       b2*(vx[j][i][k+2]-vx[j][i][k-1])+
                                       b3*(vx[j][i][k+3]-vx[j][i][k-2]))/DZ;
                                
                                vyx = (b1*(vy[j][i+1][k]-vy[j][i][k])+
                                       b2*(vy[j][i+2][k]-vy[j][i-1][k])+
                                       b3*(vy[j][i+3][k]-vy[j][i-2][k]))/DX;
                                
                                vyy = (b1*(vy[j][i][k]-vy[j-1][i][k])+
                                       b2*(vy[j+1][i][k]-vy[j-2][i][k])+
                                       b3*(vy[j+2][i][k]-vy[j-3][i][k]))/DY;
                                
                                vyz = (b1*(vy[j][i][k+1]-vy[j][i][k])+
                                       b2*(vy[j][i][k+2]-vy[j][i][k-1])+
                                       b3*(vy[j][i][k+3]-vy[j][i][k-2]))/DZ;
                                
                                vzx = (b1*(vz[j][i+1][k]-vz[j][i][k])+
                                       b2*(vz[j][i+2][k]-vz[j][i-1][k])+
                                       b3*(vz[j][i+3][k]-vz[j][i-2][k]))/DX;
                                
                                vzy = (b1*(vz[j+1][i][k]-vz[j][i][k])+
                                       b2*(vz[j+2][i][k]-vz[j-1][i][k])+
                                       b3*(vz[j+3][i][k]-vz[j-2][i][k]))/DY;
                                
                                vzz = (b1*(vz[j][i][k]-vz[j][i][k-1])+
                                       b2*(vz[j][i][k+1]-vz[j][i][k-2])+
                                       b3*(vz[j][i][k+2]-vz[j][i][k-3]))/DZ;
                                
                                
                                /* computing sums of the old memory variables in this version only one mechanism is possi
                                 ble */
                                
                                sumrxy=c1*(*(rxy_j_i+k))+c2*(*(rxy_j_i_2+k))+c3*(*(rxy_j_i_3+k))+c4*(*(rxy_j_i_4+k));
                                sumryz=c1*(*(ryz_j_i+k))+c2*(*(ryz_j_i_2+k))+c3*(*(ryz_j_i_3+k))+c4*(*(ryz_j_i_4+k));
                                sumrxz=c1*(*(rxz_j_i+k))+c2*(*(rxz_j_i_2+k))+c3*(*(rxz_j_i_3+k))+c4*(*(rxz_j_i_4+k));
                                sumrxx=c1*(*(rxx_j_i+k))+c2*(*(rxx_j_i_2+k))+c3*(*(rxx_j_i_3+k))+c4*(*(rxx_j_i_4+k));
                                sumryy=c1*(*(ryy_j_i+k))+c2*(*(ryy_j_i_2+k))+c3*(*(ryy_j_i_3+k))+c4*(*(ryy_j_i_4+k));
                                sumrzz=c1*(*(rzz_j_i+k))+c2*(*(rzz_j_i_2+k))+c3*(*(rzz_j_i_3+k))+c4*(*(rzz_j_i_4+k));
                                
                                /* updating components of the stress tensor, partially */
                                fipjp=uipjp[j][i][k]*DT*(1.0+L*tausipjp[j][i][k]);
                                fjpkp=ujpkp[j][i][k]*DT*(1.0+L*tausjpkp[j][i][k]);
                                fipkp=uipkp[j][i][k]*DT*(1.0+L*tausipkp[j][i][k]);
                                g=pi[j][i][k]*(1.0+L*taup[j][i][k]);
                                f=2.0*u[j][i][k]*(1.0+L*taus[j][i][k]);
                                
                                /* Save spatial derivations */
                                *(vxyyx_j_i+k)=vxy+vyx;
                                *(vyzzy_j_i+k)=vyz+vzy;
                                *(vxzzx_j_i+k)=vxz+vzx;
                                *(vxxyyzz_j_i+k)=vxx+vyy+vzz;
                                *(vyyzz_j_i+k)=vyy+vzz;
                                *(vxxzz_j_i+k)=vxx+vzz;
                                *(vxxyy_j_i+k)=vxx+vyy;
                                
                                /* Do the Adams Bashforth stuff and save it, we need it for the memory variables again */
                                vxyyx_T2=(c1*(*(vxyyx_j_i+k))+c2*(*(vxyyx_j_i_2+k))+c3*(*(vxyyx_j_i_3+k))+c4*(*(vxyyx_j_i_4+k)));
                                vyzzy_T2=(c1*(*(vyzzy_j_i+k))+c2*(*(vyzzy_j_i_2+k))+c3*(*(vyzzy_j_i_3+k))+c4*(*(vyzzy_j_i_4+k)));
                                vxzzx_T2=(c1*(*(vxzzx_j_i+k))+c2*(*(vxzzx_j_i_2+k))+c3*(*(vxzzx_j_i_3+k))+c4*(*(vxzzx_j_i_4+k)));
                                vyyzz_T2=(c1*(*(vyyzz_j_i+k))+c2*(*(vyyzz_j_i_2+k))+c3*(*(vyyzz_j_i_3+k))+c4*(*(vyyzz_j_i_4+k)));
                                vxxzz_T2=(c1*(*(vxxzz_j_i+k))+c2*(*(vxxzz_j_i_2+k))+c3*(*(vxxzz_j_i_3+k))+c4*(*(vxxzz_j_i_4+k)));
                                vxxyy_T2=(c1*(*(vxxyy_j_i+k))+c2*(*(vxxyy_j_i_2+k))+c3*(*(vxxyy_j_i_3+k))+c4*(*(vxxyy_j_i_4+k)));
                                vdiag=((c1*(*(vxxyyzz_j_i+k))+c2*(*(vxxyyzz_j_i_2+k))+c3*(*(vxxyyzz_j_i_3+k)))+c4*(*(vxxyyzz_j_i_4+k)));
                                
                                /* Update the stress tensor */
                                sxy[j][i][k]+=(fipjp*vxyyx_T2)+(dthalbe*sumrxy);
                                syz[j][i][k]+=(fjpkp*vyzzy_T2)+(dthalbe*sumryz);
                                sxz[j][i][k]+=(fipkp*vxzzx_T2)+(dthalbe*sumrxz);
                                sxx[j][i][k]+=DT*((g*vdiag)-(f*vyyzz_T2))+(dthalbe*sumrxx);
                                syy[j][i][k]+=DT*((g*vdiag)-(f*vxxzz_T2))+(dthalbe*sumryy);
                                szz[j][i][k]+=DT*((g*vdiag)-(f*vxxyy_T2))+(dthalbe*sumrzz);
                                
                                /* Update the memory variables */
                                
                                /* in the case of L=1 update the memory-variables without suming*/
                                l=1;
                                n1=1.0/(1.0+(c1*eta[l]*0.5)); /* 1/(1+DT*c1/(2*tau_sigl)) */
                                n2=1.0-(c1*eta[l]*0.5);
                                n3=c2*eta[l]/2.0;
                                n4=c3*eta[l]/2.0;
                                n5=c4*eta[l]/2.0;
                                dipjp=uipjp[j][i][k]*eta[l]*tausipjp[j][i][k];
                                djpkp=ujpkp[j][i][k]*eta[l]*tausjpkp[j][i][k];
                                dipkp=uipkp[j][i][k]*eta[l]*tausipkp[j][i][k];
                                d=2.0*u[j][i][k]*eta[l]*taus[j][i][k];
                                e=pi[j][i][k]*eta[l]*taup[j][i][k];
                                
                                
                                *(rxy_j_i_4+k)=n1*((*(rxy_j_i+k))*n2-n3*((*(rxy_j_i+k))+(*(rxy_j_i_2+k)))-n4*((*(rxy_j_i_2+k))+(*(rxy_j_i_3+k)))-n5*((*(rxy_j_i_3+k))+(*(rxy_j_i_4+k)))-(dipjp*vxyyx_T2));
                                *(ryz_j_i_4+k)=n1*((*(ryz_j_i+k))*n2-n3*((*(ryz_j_i+k))+(*(ryz_j_i_2+k)))-n4*((*(ryz_j_i_2+k))+(*(ryz_j_i_3+k)))-n5*((*(ryz_j_i_3+k))+(*(ryz_j_i_4+k)))-(djpkp*vyzzy_T2));
                                *(rxz_j_i_4+k)=n1*((*(rxz_j_i+k))*n2-n3*((*(rxz_j_i+k))+(*(rxz_j_i_2+k)))-n4*((*(rxz_j_i_2+k))+(*(rxz_j_i_3+k)))-n5*((*(rxz_j_i_3+k))+(*(rxz_j_i_4+k)))-(dipkp*vxzzx_T2));
                                *(rxx_j_i_4+k)=n1*((*(rxx_j_i+k))*n2-n3*((*(rxx_j_i+k))+(*(rxx_j_i_2+k)))-n4*((*(rxx_j_i_2+k))+(*(rxx_j_i_3+k)))-n5*((*(rxx_j_i_3+k))+(*(rxx_j_i_4+k)))-(e*vdiag)+(d*vyyzz_T2));
                                *(ryy_j_i_4+k)=n1*((*(ryy_j_i+k))*n2-n3*((*(ryy_j_i+k))+(*(ryy_j_i_2+k)))-n4*((*(ryy_j_i_2+k))+(*(ryy_j_i_3+k)))-n5*((*(ryy_j_i_3+k))+(*(ryy_j_i_4+k)))-(e*vdiag)+(d*vxxzz_T2));
                                *(rzz_j_i_4+k)=n1*((*(rzz_j_i+k))*n2-n3*((*(rzz_j_i+k))+(*(rzz_j_i_2+k)))-n4*((*(rzz_j_i_2+k))+(*(rzz_j_i_3+k)))-n5*((*(rzz_j_i_3+k))+(*(rzz_j_i_4+k)))-(e*vdiag)+(d*vxxyy_T2));
                                
                                
                                sumrxy=c1*(*(rxy_j_i_4+k))+c2*(*(rxy_j_i+k))+c3*(*(rxy_j_i_2+k))+c4*(*(rxy_j_i_3+k));
                                sumryz=c1*(*(ryz_j_i_4+k))+c2*(*(ryz_j_i+k))+c3*(*(ryz_j_i_2+k))+c4*(*(ryz_j_i_3+k));
                                sumrxz=c1*(*(rxz_j_i_4+k))+c2*(*(rxz_j_i+k))+c3*(*(rxz_j_i_2+k))+c4*(*(rxz_j_i_3+k));
                                sumrxx=c1*(*(rxx_j_i_4+k))+c2*(*(rxx_j_i+k))+c3*(*(rxx_j_i_2+k))+c4*(*(rxx_j_i_3+k));
                                sumryy=c1*(*(ryy_j_i_4+k))+c2*(*(ryy_j_i+k))+c3*(*(ryy_j_i_2+k))+c4*(*(ryy_j_i_3+k));
                                sumrzz=c1*(*(rzz_j_i_4+k))+c2*(*(rzz_j_i+k))+c3*(*(rzz_j_i_2+k))+c4*(*(rzz_j_i_3+k));
                                
                                
                                /* and now the components of the stress tensor are
                                 completely updated */
                                sxy[j][i][k]+=(dthalbe*sumrxy);
                                syz[j][i][k]+=(dthalbe*sumryz);
                                sxz[j][i][k]+=(dthalbe*sumrxz);
                                sxx[j][i][k]+=(dthalbe*sumrxx);
                                syy[j][i][k]+=(dthalbe*sumryy);
                                szz[j][i][k]+=(dthalbe*sumrzz);
                                
                                

                            }
                        }
                    }
                    break;
                    
                case 8 :
                    
                    b1=1225.0/1024.0; b2=-245.0/3072.0; b3=49.0/5120.0; b4=-5.0/7168.0; /* Taylor coefficients */
                    if(FDCOEFF==2){
                        b1=1.2257; b2=-0.099537; b3=0.018063; b4=-0.0026274;} /* Holberg coefficients E=0.1 %*/
                    
                    for (j=ny1;j<=ny2;j++){
                        
                        rxx_j=*(rxx+j);ryy_j=*(ryy+j);rzz_j=*(rzz+j);rxy_j=*(rxy+j);ryz_j=*(ryz+j);rxz_j=*(rxz+j);
                        rxx_j_2=*(rxx_2+j);ryy_j_2=*(ryy_2+j);rzz_j_2=*(rzz_2+j);rxy_j_2=*(rxy_2+j);ryz_j_2=*(ryz_2+j); rxz_j_2=*(rxz_2+j);
                        rxx_j_3=*(rxx_3+j);ryy_j_3=*(ryy_3+j);rzz_j_3=*(rzz_3+j);rxy_j_3=*(rxy_3+j);ryz_j_3=*(ryz_3+j);rxz_j_3=*(rxz_3+j);
                        rxx_j_4=*(rxx_4+j);ryy_j_4=*(ryy_4+j);rzz_j_4=*(rzz_4+j);rxy_j_4=*(rxy_4+j);ryz_j_4=*(ryz_4+j);rxz_j_4=*(rxz_4+j);
                        vxyyx_j=*(vxyyx+j);vyzzy_j=*(vyzzy+j);vxzzx_j=*(vxzzx+j);vxxyyzz_j=*(vxxyyzz+j);vyyzz_j=*(vyyzz+j);vxxzz_j=*(vxxzz+j);vxxyy_j=*(vxxyy+j);
                        vxyyx_j_2=*(vxyyx_2+j);vyzzy_j_2=*(vyzzy_2+j); vxzzx_j_2=*(vxzzx_2+j);vxxyyzz_j_2=*(vxxyyzz_2+j);vyyzz_j_2=*(vyyzz_2+j);vxxzz_j_2=*(vxxzz_2+j);vxxyy_j_2=*(vxxyy_2+j);
                        vxyyx_j_3=*(vxyyx_3+j);vyzzy_j_3=*(vyzzy_3+j);vxzzx_j_3=*(vxzzx_3+j);vxxyyzz_j_3=*(vxxyyzz_3+j);vyyzz_j_3=*(vyyzz_3+j);vxxzz_j_3=*(vxxzz_3+j);vxxyy_j_3=*(vxxyy_3+j);
                        vxyyx_j_4=*(vxyyx_4+j);vyzzy_j_4=*(vyzzy_4+j);vxzzx_j_4=*(vxzzx_4+j);vxxyyzz_j_4=*(vxxyyzz_4+j);vyyzz_j_4=*(vyyzz_4+j);vxxzz_j_4=*(vxxzz_4+j);vxxyy_j_4=*(vxxyy_4+j);
                        
                        for (i=nx1;i<=nx2;i++){
                            
                            rxx_j_i=*(rxx_j+i);ryy_j_i=*(ryy_j+i);rzz_j_i=*(rzz_j+i);rxy_j_i=*(rxy_j+i);ryz_j_i=*(ryz_j+i);rxz_j_i=*(rxz_j+i);
                            rxx_j_i_2=*(rxx_j_2+i);ryy_j_i_2=*(ryy_j_2+i);rzz_j_i_2=*(rzz_j_2+i);rxy_j_i_2=*(rxy_j_2+i);ryz_j_i_2=*(ryz_j_2+i);rxz_j_i_2=*(rxz_j_2+i);
                            rxx_j_i_3=*(rxx_j_3+i);ryy_j_i_3=*(ryy_j_3+i);rzz_j_i_3=*(rzz_j_3+i);rxy_j_i_3=*(rxy_j_3+i);ryz_j_i_3=*(ryz_j_3+i);rxz_j_i_3=*(rxz_j_3+i);
                            rxx_j_i_4=*(rxx_j_4+i);ryy_j_i_4=*(ryy_j_4+i);rzz_j_i_4=*(rzz_j_4+i);rxy_j_i_4=*(rxy_j_4+i);ryz_j_i_4=*(ryz_j_4+i);rxz_j_i_4=*(rxz_j_4+i);
                            vxyyx_j_i=*(vxyyx_j+i);vyzzy_j_i=*(vyzzy_j+i);vxzzx_j_i=*(vxzzx_j+i);vxxyyzz_j_i=*(vxxyyzz_j+i);vyyzz_j_i=*(vyyzz_j+i);vxxzz_j_i=*(vxxzz_j+i);vxxyy_j_i=*(vxxyy_j+i);
                            vxyyx_j_i_2=*(vxyyx_j_2+i);vyzzy_j_i_2=*(vyzzy_j_2+i);vxzzx_j_i_2=*(vxzzx_j_2+i);vxxyyzz_j_i_2=*(vxxyyzz_j_2+i);vyyzz_j_i_2=*(vyyzz_j_2+i);vxxzz_j_i_2=*(vxxzz_j_2+i);vxxyy_j_i_2=*(vxxyy_j_2+i);
                            vxyyx_j_i_3=*(vxyyx_j_3+i);vyzzy_j_i_3=*(vyzzy_j_3+i);vxzzx_j_i_3=*(vxzzx_j_3+i); vxxyyzz_j_i_3=*(vxxyyzz_j_3+i);vyyzz_j_i_3=*(vyyzz_j_3+i);vxxzz_j_i_3=*(vxxzz_j_3+i); vxxyy_j_i_3=*(vxxyy_j_3+i);
                            vxyyx_j_i_4=*(vxyyx_j_4+i);vyzzy_j_i_4=*(vyzzy_j_4+i);vxzzx_j_i_4=*(vxzzx_j_4+i); vxxyyzz_j_i_4=*(vxxyyzz_j_4+i);vyyzz_j_i_4=*(vyyzz_j_4+i);vxxzz_j_i_4=*(vxxzz_j_4+i);vxxyy_j_i_4=*(vxxyy_j_4+i);
                            
                            for (k=nz1;k<=nz2;k++){
                                /* spatial derivatives of the components of the velocities
                                 are computed */
                                
                                vxx = (b1*(vx[j][i][k]-vx[j][i-1][k])+
                                       b2*(vx[j][i+1][k]-vx[j][i-2][k])+
                                       b3*(vx[j][i+2][k]-vx[j][i-3][k])+
                                       b4*(vx[j][i+3][k]-vx[j][i-4][k]))/DX;
                                
                                vxy = (b1*(vx[j+1][i][k]-vx[j][i][k])+
                                       b2*(vx[j+2][i][k]-vx[j-1][i][k])+
                                       b3*(vx[j+3][i][k]-vx[j-2][i][k])+
                                       b4*(vx[j+4][i][k]-vx[j-3][i][k]))/DY;
                                
                                vxz = (b1*(vx[j][i][k+1]-vx[j][i][k])+
                                       b2*(vx[j][i][k+2]-vx[j][i][k-1])+
                                       b3*(vx[j][i][k+3]-vx[j][i][k-2])+
                                       b4*(vx[j][i][k+4]-vx[j][i][k-3]))/DZ;
                                
                                vyx = (b1*(vy[j][i+1][k]-vy[j][i][k])+
                                       b2*(vy[j][i+2][k]-vy[j][i-1][k])+
                                       b3*(vy[j][i+3][k]-vy[j][i-2][k])+
                                       b4*(vy[j][i+4][k]-vy[j][i-3][k]))/DX;
                                
                                vyy = (b1*(vy[j][i][k]-vy[j-1][i][k])+
                                       b2*(vy[j+1][i][k]-vy[j-2][i][k])+
                                       b3*(vy[j+2][i][k]-vy[j-3][i][k])+
                                       b4*(vy[j+3][i][k]-vy[j-4][i][k]))/DY;
                                
                                vyz = (b1*(vy[j][i][k+1]-vy[j][i][k])+
                                       b2*(vy[j][i][k+2]-vy[j][i][k-1])+
                                       b3*(vy[j][i][k+3]-vy[j][i][k-2])+
                                       b4*(vy[j][i][k+4]-vy[j][i][k-3]))/DZ;
                                
                                vzx = (b1*(vz[j][i+1][k]-vz[j][i][k])+
                                       b2*(vz[j][i+2][k]-vz[j][i-1][k])+
                                       b3*(vz[j][i+3][k]-vz[j][i-2][k])+
                                       b4*(vz[j][i+4][k]-vz[j][i-3][k]))/DX;
                                
                                vzy = (b1*(vz[j+1][i][k]-vz[j][i][k])+
                                       b2*(vz[j+2][i][k]-vz[j-1][i][k])+
                                       b3*(vz[j+3][i][k]-vz[j-2][i][k])+
                                       b4*(vz[j+4][i][k]-vz[j-3][i][k]))/DY;
                                
                                vzz = (b1*(vz[j][i][k]-vz[j][i][k-1])+
                                       b2*(vz[j][i][k+1]-vz[j][i][k-2])+
                                       b3*(vz[j][i][k+2]-vz[j][i][k-3])+
                                       b4*(vz[j][i][k+3]-vz[j][i][k-4]))/DZ;
                                
                                
                                /* computing sums of the old memory variables in this version only one mechanism is possi
                                 ble */
                                
                                sumrxy=c1*(*(rxy_j_i+k))+c2*(*(rxy_j_i_2+k))+c3*(*(rxy_j_i_3+k))+c4*(*(rxy_j_i_4+k));
                                sumryz=c1*(*(ryz_j_i+k))+c2*(*(ryz_j_i_2+k))+c3*(*(ryz_j_i_3+k))+c4*(*(ryz_j_i_4+k));
                                sumrxz=c1*(*(rxz_j_i+k))+c2*(*(rxz_j_i_2+k))+c3*(*(rxz_j_i_3+k))+c4*(*(rxz_j_i_4+k));
                                sumrxx=c1*(*(rxx_j_i+k))+c2*(*(rxx_j_i_2+k))+c3*(*(rxx_j_i_3+k))+c4*(*(rxx_j_i_4+k));
                                sumryy=c1*(*(ryy_j_i+k))+c2*(*(ryy_j_i_2+k))+c3*(*(ryy_j_i_3+k))+c4*(*(ryy_j_i_4+k));
                                sumrzz=c1*(*(rzz_j_i+k))+c2*(*(rzz_j_i_2+k))+c3*(*(rzz_j_i_3+k))+c4*(*(rzz_j_i_4+k));
                                
                                /* updating components of the stress tensor, partially */
                                fipjp=uipjp[j][i][k]*DT*(1.0+L*tausipjp[j][i][k]);
                                fjpkp=ujpkp[j][i][k]*DT*(1.0+L*tausjpkp[j][i][k]);
                                fipkp=uipkp[j][i][k]*DT*(1.0+L*tausipkp[j][i][k]);
                                g=pi[j][i][k]*(1.0+L*taup[j][i][k]);
                                f=2.0*u[j][i][k]*(1.0+L*taus[j][i][k]);
                                
                                /* Save spatial derivations */
                                *(vxyyx_j_i+k)=vxy+vyx;
                                *(vyzzy_j_i+k)=vyz+vzy;
                                *(vxzzx_j_i+k)=vxz+vzx;
                                *(vxxyyzz_j_i+k)=vxx+vyy+vzz;
                                *(vyyzz_j_i+k)=vyy+vzz;
                                *(vxxzz_j_i+k)=vxx+vzz;
                                *(vxxyy_j_i+k)=vxx+vyy;
                                
                                /* Do the Adams Bashforth stuff and save it, we need it for the memory variables again */
                                vxyyx_T2=(c1*(*(vxyyx_j_i+k))+c2*(*(vxyyx_j_i_2+k))+c3*(*(vxyyx_j_i_3+k))+c4*(*(vxyyx_j_i_4+k)));
                                vyzzy_T2=(c1*(*(vyzzy_j_i+k))+c2*(*(vyzzy_j_i_2+k))+c3*(*(vyzzy_j_i_3+k))+c4*(*(vyzzy_j_i_4+k)));
                                vxzzx_T2=(c1*(*(vxzzx_j_i+k))+c2*(*(vxzzx_j_i_2+k))+c3*(*(vxzzx_j_i_3+k))+c4*(*(vxzzx_j_i_4+k)));
                                vyyzz_T2=(c1*(*(vyyzz_j_i+k))+c2*(*(vyyzz_j_i_2+k))+c3*(*(vyyzz_j_i_3+k))+c4*(*(vyyzz_j_i_4+k)));
                                vxxzz_T2=(c1*(*(vxxzz_j_i+k))+c2*(*(vxxzz_j_i_2+k))+c3*(*(vxxzz_j_i_3+k))+c4*(*(vxxzz_j_i_4+k)));
                                vxxyy_T2=(c1*(*(vxxyy_j_i+k))+c2*(*(vxxyy_j_i_2+k))+c3*(*(vxxyy_j_i_3+k))+c4*(*(vxxyy_j_i_4+k)));
                                vdiag=((c1*(*(vxxyyzz_j_i+k))+c2*(*(vxxyyzz_j_i_2+k))+c3*(*(vxxyyzz_j_i_3+k)))+c4*(*(vxxyyzz_j_i_4+k)));
                                
                                /* Update the stress tensor */
                                sxy[j][i][k]+=(fipjp*vxyyx_T2)+(dthalbe*sumrxy);
                                syz[j][i][k]+=(fjpkp*vyzzy_T2)+(dthalbe*sumryz);
                                sxz[j][i][k]+=(fipkp*vxzzx_T2)+(dthalbe*sumrxz);
                                sxx[j][i][k]+=DT*((g*vdiag)-(f*vyyzz_T2))+(dthalbe*sumrxx);
                                syy[j][i][k]+=DT*((g*vdiag)-(f*vxxzz_T2))+(dthalbe*sumryy);
                                szz[j][i][k]+=DT*((g*vdiag)-(f*vxxyy_T2))+(dthalbe*sumrzz);
                                
                                /* Update the memory variables */
                                
                                /* in the case of L=1 update the memory-variables without suming*/
                                l=1;
                                n1=1.0/(1.0+(c1*eta[l]*0.5)); /* 1/(1+DT*c1/(2*tau_sigl)) */
                                n2=1.0-(c1*eta[l]*0.5);
                                n3=c2*eta[l]/2.0;
                                n4=c3*eta[l]/2.0;
                                n5=c4*eta[l]/2.0;
                                dipjp=uipjp[j][i][k]*eta[l]*tausipjp[j][i][k];
                                djpkp=ujpkp[j][i][k]*eta[l]*tausjpkp[j][i][k];
                                dipkp=uipkp[j][i][k]*eta[l]*tausipkp[j][i][k];
                                d=2.0*u[j][i][k]*eta[l]*taus[j][i][k];
                                e=pi[j][i][k]*eta[l]*taup[j][i][k];
                                
                                
                                *(rxy_j_i_4+k)=n1*((*(rxy_j_i+k))*n2-n3*((*(rxy_j_i+k))+(*(rxy_j_i_2+k)))-n4*((*(rxy_j_i_2+k))+(*(rxy_j_i_3+k)))-n5*((*(rxy_j_i_3+k))+(*(rxy_j_i_4+k)))-(dipjp*vxyyx_T2));
                                *(ryz_j_i_4+k)=n1*((*(ryz_j_i+k))*n2-n3*((*(ryz_j_i+k))+(*(ryz_j_i_2+k)))-n4*((*(ryz_j_i_2+k))+(*(ryz_j_i_3+k)))-n5*((*(ryz_j_i_3+k))+(*(ryz_j_i_4+k)))-(djpkp*vyzzy_T2));
                                *(rxz_j_i_4+k)=n1*((*(rxz_j_i+k))*n2-n3*((*(rxz_j_i+k))+(*(rxz_j_i_2+k)))-n4*((*(rxz_j_i_2+k))+(*(rxz_j_i_3+k)))-n5*((*(rxz_j_i_3+k))+(*(rxz_j_i_4+k)))-(dipkp*vxzzx_T2));
                                *(rxx_j_i_4+k)=n1*((*(rxx_j_i+k))*n2-n3*((*(rxx_j_i+k))+(*(rxx_j_i_2+k)))-n4*((*(rxx_j_i_2+k))+(*(rxx_j_i_3+k)))-n5*((*(rxx_j_i_3+k))+(*(rxx_j_i_4+k)))-(e*vdiag)+(d*vyyzz_T2));
                                *(ryy_j_i_4+k)=n1*((*(ryy_j_i+k))*n2-n3*((*(ryy_j_i+k))+(*(ryy_j_i_2+k)))-n4*((*(ryy_j_i_2+k))+(*(ryy_j_i_3+k)))-n5*((*(ryy_j_i_3+k))+(*(ryy_j_i_4+k)))-(e*vdiag)+(d*vxxzz_T2));
                                *(rzz_j_i_4+k)=n1*((*(rzz_j_i+k))*n2-n3*((*(rzz_j_i+k))+(*(rzz_j_i_2+k)))-n4*((*(rzz_j_i_2+k))+(*(rzz_j_i_3+k)))-n5*((*(rzz_j_i_3+k))+(*(rzz_j_i_4+k)))-(e*vdiag)+(d*vxxyy_T2));
                                
                                
                                sumrxy=c1*(*(rxy_j_i_4+k))+c2*(*(rxy_j_i+k))+c3*(*(rxy_j_i_2+k))+c4*(*(rxy_j_i_3+k));
                                sumryz=c1*(*(ryz_j_i_4+k))+c2*(*(ryz_j_i+k))+c3*(*(ryz_j_i_2+k))+c4*(*(ryz_j_i_3+k));
                                sumrxz=c1*(*(rxz_j_i_4+k))+c2*(*(rxz_j_i+k))+c3*(*(rxz_j_i_2+k))+c4*(*(rxz_j_i_3+k));
                                sumrxx=c1*(*(rxx_j_i_4+k))+c2*(*(rxx_j_i+k))+c3*(*(rxx_j_i_2+k))+c4*(*(rxx_j_i_3+k));
                                sumryy=c1*(*(ryy_j_i_4+k))+c2*(*(ryy_j_i+k))+c3*(*(ryy_j_i_2+k))+c4*(*(ryy_j_i_3+k));
                                sumrzz=c1*(*(rzz_j_i_4+k))+c2*(*(rzz_j_i+k))+c3*(*(rzz_j_i_2+k))+c4*(*(rzz_j_i_3+k));
                                
                                
                                /* and now the components of the stress tensor are
                                 completely updated */
                                sxy[j][i][k]+=(dthalbe*sumrxy);
                                syz[j][i][k]+=(dthalbe*sumryz);
                                sxz[j][i][k]+=(dthalbe*sumrxz);
                                sxx[j][i][k]+=(dthalbe*sumrxx);
                                syy[j][i][k]+=(dthalbe*sumryy);
                                szz[j][i][k]+=(dthalbe*sumrzz);
                                
                                

                                
                            }
                        }
                    }
                    break;
                    
                case 10 :
                    
                    b1=19845.0/16384.0; b2=-735.0/8192.0; b3=567.0/40960.0; b4=-405.0/229376.0; b5=35.0/294912.0; /* Taylor coefficients */
                    if(FDCOEFF==2){
                        b1=1.2415; b2=-0.11231; b3=0.026191; b4=-0.0064682; b5=0.001191;} /* Holberg coefficients E=0.1 %*/
                    
                    
                    for (j=ny1;j<=ny2;j++){
                        
                        rxx_j=*(rxx+j);ryy_j=*(ryy+j);rzz_j=*(rzz+j);rxy_j=*(rxy+j);ryz_j=*(ryz+j);rxz_j=*(rxz+j);
                        rxx_j_2=*(rxx_2+j);ryy_j_2=*(ryy_2+j);rzz_j_2=*(rzz_2+j);rxy_j_2=*(rxy_2+j);ryz_j_2=*(ryz_2+j); rxz_j_2=*(rxz_2+j);
                        rxx_j_3=*(rxx_3+j);ryy_j_3=*(ryy_3+j);rzz_j_3=*(rzz_3+j);rxy_j_3=*(rxy_3+j);ryz_j_3=*(ryz_3+j);rxz_j_3=*(rxz_3+j);
                        rxx_j_4=*(rxx_4+j);ryy_j_4=*(ryy_4+j);rzz_j_4=*(rzz_4+j);rxy_j_4=*(rxy_4+j);ryz_j_4=*(ryz_4+j);rxz_j_4=*(rxz_4+j);
                        vxyyx_j=*(vxyyx+j);vyzzy_j=*(vyzzy+j);vxzzx_j=*(vxzzx+j);vxxyyzz_j=*(vxxyyzz+j);vyyzz_j=*(vyyzz+j);vxxzz_j=*(vxxzz+j);vxxyy_j=*(vxxyy+j);
                        vxyyx_j_2=*(vxyyx_2+j);vyzzy_j_2=*(vyzzy_2+j); vxzzx_j_2=*(vxzzx_2+j);vxxyyzz_j_2=*(vxxyyzz_2+j);vyyzz_j_2=*(vyyzz_2+j);vxxzz_j_2=*(vxxzz_2+j);vxxyy_j_2=*(vxxyy_2+j);
                        vxyyx_j_3=*(vxyyx_3+j);vyzzy_j_3=*(vyzzy_3+j);vxzzx_j_3=*(vxzzx_3+j);vxxyyzz_j_3=*(vxxyyzz_3+j);vyyzz_j_3=*(vyyzz_3+j);vxxzz_j_3=*(vxxzz_3+j);vxxyy_j_3=*(vxxyy_3+j);
                        vxyyx_j_4=*(vxyyx_4+j);vyzzy_j_4=*(vyzzy_4+j);vxzzx_j_4=*(vxzzx_4+j);vxxyyzz_j_4=*(vxxyyzz_4+j);vyyzz_j_4=*(vyyzz_4+j);vxxzz_j_4=*(vxxzz_4+j);vxxyy_j_4=*(vxxyy_4+j);
                        
                        for (i=nx1;i<=nx2;i++){
                            
                            rxx_j_i=*(rxx_j+i);ryy_j_i=*(ryy_j+i);rzz_j_i=*(rzz_j+i);rxy_j_i=*(rxy_j+i);ryz_j_i=*(ryz_j+i);rxz_j_i=*(rxz_j+i);
                            rxx_j_i_2=*(rxx_j_2+i);ryy_j_i_2=*(ryy_j_2+i);rzz_j_i_2=*(rzz_j_2+i);rxy_j_i_2=*(rxy_j_2+i);ryz_j_i_2=*(ryz_j_2+i);rxz_j_i_2=*(rxz_j_2+i);
                            rxx_j_i_3=*(rxx_j_3+i);ryy_j_i_3=*(ryy_j_3+i);rzz_j_i_3=*(rzz_j_3+i);rxy_j_i_3=*(rxy_j_3+i);ryz_j_i_3=*(ryz_j_3+i);rxz_j_i_3=*(rxz_j_3+i);
                            rxx_j_i_4=*(rxx_j_4+i);ryy_j_i_4=*(ryy_j_4+i);rzz_j_i_4=*(rzz_j_4+i);rxy_j_i_4=*(rxy_j_4+i);ryz_j_i_4=*(ryz_j_4+i);rxz_j_i_4=*(rxz_j_4+i);
                            vxyyx_j_i=*(vxyyx_j+i);vyzzy_j_i=*(vyzzy_j+i);vxzzx_j_i=*(vxzzx_j+i);vxxyyzz_j_i=*(vxxyyzz_j+i);vyyzz_j_i=*(vyyzz_j+i);vxxzz_j_i=*(vxxzz_j+i);vxxyy_j_i=*(vxxyy_j+i);
                            vxyyx_j_i_2=*(vxyyx_j_2+i);vyzzy_j_i_2=*(vyzzy_j_2+i);vxzzx_j_i_2=*(vxzzx_j_2+i);vxxyyzz_j_i_2=*(vxxyyzz_j_2+i);vyyzz_j_i_2=*(vyyzz_j_2+i);vxxzz_j_i_2=*(vxxzz_j_2+i);vxxyy_j_i_2=*(vxxyy_j_2+i);
                            vxyyx_j_i_3=*(vxyyx_j_3+i);vyzzy_j_i_3=*(vyzzy_j_3+i);vxzzx_j_i_3=*(vxzzx_j_3+i); vxxyyzz_j_i_3=*(vxxyyzz_j_3+i);vyyzz_j_i_3=*(vyyzz_j_3+i);vxxzz_j_i_3=*(vxxzz_j_3+i); vxxyy_j_i_3=*(vxxyy_j_3+i);
                            vxyyx_j_i_4=*(vxyyx_j_4+i);vyzzy_j_i_4=*(vyzzy_j_4+i);vxzzx_j_i_4=*(vxzzx_j_4+i); vxxyyzz_j_i_4=*(vxxyyzz_j_4+i);vyyzz_j_i_4=*(vyyzz_j_4+i);vxxzz_j_i_4=*(vxxzz_j_4+i);vxxyy_j_i_4=*(vxxyy_j_4+i);
                            
                            for (k=nz1;k<=nz2;k++){
                                /* spatial derivatives of the components of the velocities
                                 are computed */
                                
                                vxx = (b1*(vx[j][i][k]-vx[j][i-1][k])+
                                       b2*(vx[j][i+1][k]-vx[j][i-2][k])+
                                       b3*(vx[j][i+2][k]-vx[j][i-3][k])+
                                       b4*(vx[j][i+3][k]-vx[j][i-4][k])+
                                       b5*(vx[j][i+4][k]-vx[j][i-5][k]))/DX;
                                
                                vxy = (b1*(vx[j+1][i][k]-vx[j][i][k])+
                                       b2*(vx[j+2][i][k]-vx[j-1][i][k])+
                                       b3*(vx[j+3][i][k]-vx[j-2][i][k])+
                                       b4*(vx[j+4][i][k]-vx[j-3][i][k])+
                                       b5*(vx[j+5][i][k]-vx[j-4][i][k]))/DY;
                                
                                vxz = (b1*(vx[j][i][k+1]-vx[j][i][k])+
                                       b2*(vx[j][i][k+2]-vx[j][i][k-1])+
                                       b3*(vx[j][i][k+3]-vx[j][i][k-2])+
                                       b4*(vx[j][i][k+4]-vx[j][i][k-3])+
                                       b5*(vx[j][i][k+5]-vx[j][i][k-4]))/DZ;
                                
                                vyx = (b1*(vy[j][i+1][k]-vy[j][i][k])+
                                       b2*(vy[j][i+2][k]-vy[j][i-1][k])+
                                       b3*(vy[j][i+3][k]-vy[j][i-2][k])+
                                       b4*(vy[j][i+4][k]-vy[j][i-3][k])+
                                       b5*(vy[j][i+5][k]-vy[j][i-4][k]))/DX;
                                
                                vyy = (b1*(vy[j][i][k]-vy[j-1][i][k])+
                                       b2*(vy[j+1][i][k]-vy[j-2][i][k])+
                                       b3*(vy[j+2][i][k]-vy[j-3][i][k])+
                                       b4*(vy[j+3][i][k]-vy[j-4][i][k])+
                                       b5*(vy[j+4][i][k]-vy[j-5][i][k]))/DY;
                                
                                vyz = (b1*(vy[j][i][k+1]-vy[j][i][k])+
                                       b2*(vy[j][i][k+2]-vy[j][i][k-1])+
                                       b3*(vy[j][i][k+3]-vy[j][i][k-2])+
                                       b4*(vy[j][i][k+4]-vy[j][i][k-3])+
                                       b5*(vy[j][i][k+5]-vy[j][i][k-4]))/DZ;
                                
                                vzx = (b1*(vz[j][i+1][k]-vz[j][i][k])+
                                       b2*(vz[j][i+2][k]-vz[j][i-1][k])+
                                       b3*(vz[j][i+3][k]-vz[j][i-2][k])+
                                       b4*(vz[j][i+4][k]-vz[j][i-3][k])+
                                       b5*(vz[j][i+5][k]-vz[j][i-4][k]))/DX;
                                
                                vzy = (b1*(vz[j+1][i][k]-vz[j][i][k])+
                                       b2*(vz[j+2][i][k]-vz[j-1][i][k])+
                                       b3*(vz[j+3][i][k]-vz[j-2][i][k])+
                                       b4*(vz[j+4][i][k]-vz[j-3][i][k])+
                                       b5*(vz[j+5][i][k]-vz[j-4][i][k]))/DY;
                                
                                vzz = (b1*(vz[j][i][k]-vz[j][i][k-1])+
                                       b2*(vz[j][i][k+1]-vz[j][i][k-2])+
                                       b3*(vz[j][i][k+2]-vz[j][i][k-3])+
                                       b4*(vz[j][i][k+3]-vz[j][i][k-4])+
                                       b5*(vz[j][i][k+4]-vz[j][i][k-5]))/DZ;
                                
                                
                                /* computing sums of the old memory variables in this version only one mechanism is possi
                                 ble */
                                
                                sumrxy=c1*(*(rxy_j_i+k))+c2*(*(rxy_j_i_2+k))+c3*(*(rxy_j_i_3+k))+c4*(*(rxy_j_i_4+k));
                                sumryz=c1*(*(ryz_j_i+k))+c2*(*(ryz_j_i_2+k))+c3*(*(ryz_j_i_3+k))+c4*(*(ryz_j_i_4+k));
                                sumrxz=c1*(*(rxz_j_i+k))+c2*(*(rxz_j_i_2+k))+c3*(*(rxz_j_i_3+k))+c4*(*(rxz_j_i_4+k));
                                sumrxx=c1*(*(rxx_j_i+k))+c2*(*(rxx_j_i_2+k))+c3*(*(rxx_j_i_3+k))+c4*(*(rxx_j_i_4+k));
                                sumryy=c1*(*(ryy_j_i+k))+c2*(*(ryy_j_i_2+k))+c3*(*(ryy_j_i_3+k))+c4*(*(ryy_j_i_4+k));
                                sumrzz=c1*(*(rzz_j_i+k))+c2*(*(rzz_j_i_2+k))+c3*(*(rzz_j_i_3+k))+c4*(*(rzz_j_i_4+k));
                                
                                /* updating components of the stress tensor, partially */
                                fipjp=uipjp[j][i][k]*DT*(1.0+L*tausipjp[j][i][k]);
                                fjpkp=ujpkp[j][i][k]*DT*(1.0+L*tausjpkp[j][i][k]);
                                fipkp=uipkp[j][i][k]*DT*(1.0+L*tausipkp[j][i][k]);
                                g=pi[j][i][k]*(1.0+L*taup[j][i][k]);
                                f=2.0*u[j][i][k]*(1.0+L*taus[j][i][k]);
                                
                                /* Save spatial derivations */
                                *(vxyyx_j_i+k)=vxy+vyx;
                                *(vyzzy_j_i+k)=vyz+vzy;
                                *(vxzzx_j_i+k)=vxz+vzx;
                                *(vxxyyzz_j_i+k)=vxx+vyy+vzz;
                                *(vyyzz_j_i+k)=vyy+vzz;
                                *(vxxzz_j_i+k)=vxx+vzz;
                                *(vxxyy_j_i+k)=vxx+vyy;
                                
                                /* Do the Adams Bashforth stuff and save it, we need it for the memory variables again */
                                vxyyx_T2=(c1*(*(vxyyx_j_i+k))+c2*(*(vxyyx_j_i_2+k))+c3*(*(vxyyx_j_i_3+k))+c4*(*(vxyyx_j_i_4+k)));
                                vyzzy_T2=(c1*(*(vyzzy_j_i+k))+c2*(*(vyzzy_j_i_2+k))+c3*(*(vyzzy_j_i_3+k))+c4*(*(vyzzy_j_i_4+k)));
                                vxzzx_T2=(c1*(*(vxzzx_j_i+k))+c2*(*(vxzzx_j_i_2+k))+c3*(*(vxzzx_j_i_3+k))+c4*(*(vxzzx_j_i_4+k)));
                                vyyzz_T2=(c1*(*(vyyzz_j_i+k))+c2*(*(vyyzz_j_i_2+k))+c3*(*(vyyzz_j_i_3+k))+c4*(*(vyyzz_j_i_4+k)));
                                vxxzz_T2=(c1*(*(vxxzz_j_i+k))+c2*(*(vxxzz_j_i_2+k))+c3*(*(vxxzz_j_i_3+k))+c4*(*(vxxzz_j_i_4+k)));
                                vxxyy_T2=(c1*(*(vxxyy_j_i+k))+c2*(*(vxxyy_j_i_2+k))+c3*(*(vxxyy_j_i_3+k))+c4*(*(vxxyy_j_i_4+k)));
                                vdiag=((c1*(*(vxxyyzz_j_i+k))+c2*(*(vxxyyzz_j_i_2+k))+c3*(*(vxxyyzz_j_i_3+k)))+c4*(*(vxxyyzz_j_i_4+k)));
                                
                                /* Update the stress tensor */
                                sxy[j][i][k]+=(fipjp*vxyyx_T2)+(dthalbe*sumrxy);
                                syz[j][i][k]+=(fjpkp*vyzzy_T2)+(dthalbe*sumryz);
                                sxz[j][i][k]+=(fipkp*vxzzx_T2)+(dthalbe*sumrxz);
                                sxx[j][i][k]+=DT*((g*vdiag)-(f*vyyzz_T2))+(dthalbe*sumrxx);
                                syy[j][i][k]+=DT*((g*vdiag)-(f*vxxzz_T2))+(dthalbe*sumryy);
                                szz[j][i][k]+=DT*((g*vdiag)-(f*vxxyy_T2))+(dthalbe*sumrzz);
                                
                                /* Update the memory variables */
                                
                                /* in the case of L=1 update the memory-variables without suming*/
                                l=1;
                                n1=1.0/(1.0+(c1*eta[l]*0.5)); /* 1/(1+DT*c1/(2*tau_sigl)) */
                                n2=1.0-(c1*eta[l]*0.5);
                                n3=c2*eta[l]/2.0;
                                n4=c3*eta[l]/2.0;
                                n5=c4*eta[l]/2.0;
                                dipjp=uipjp[j][i][k]*eta[l]*tausipjp[j][i][k];
                                djpkp=ujpkp[j][i][k]*eta[l]*tausjpkp[j][i][k];
                                dipkp=uipkp[j][i][k]*eta[l]*tausipkp[j][i][k];
                                d=2.0*u[j][i][k]*eta[l]*taus[j][i][k];
                                e=pi[j][i][k]*eta[l]*taup[j][i][k];
                                
                                
                                *(rxy_j_i_4+k)=n1*((*(rxy_j_i+k))*n2-n3*((*(rxy_j_i+k))+(*(rxy_j_i_2+k)))-n4*((*(rxy_j_i_2+k))+(*(rxy_j_i_3+k)))-n5*((*(rxy_j_i_3+k))+(*(rxy_j_i_4+k)))-(dipjp*vxyyx_T2));
                                *(ryz_j_i_4+k)=n1*((*(ryz_j_i+k))*n2-n3*((*(ryz_j_i+k))+(*(ryz_j_i_2+k)))-n4*((*(ryz_j_i_2+k))+(*(ryz_j_i_3+k)))-n5*((*(ryz_j_i_3+k))+(*(ryz_j_i_4+k)))-(djpkp*vyzzy_T2));
                                *(rxz_j_i_4+k)=n1*((*(rxz_j_i+k))*n2-n3*((*(rxz_j_i+k))+(*(rxz_j_i_2+k)))-n4*((*(rxz_j_i_2+k))+(*(rxz_j_i_3+k)))-n5*((*(rxz_j_i_3+k))+(*(rxz_j_i_4+k)))-(dipkp*vxzzx_T2));
                                *(rxx_j_i_4+k)=n1*((*(rxx_j_i+k))*n2-n3*((*(rxx_j_i+k))+(*(rxx_j_i_2+k)))-n4*((*(rxx_j_i_2+k))+(*(rxx_j_i_3+k)))-n5*((*(rxx_j_i_3+k))+(*(rxx_j_i_4+k)))-(e*vdiag)+(d*vyyzz_T2));
                                *(ryy_j_i_4+k)=n1*((*(ryy_j_i+k))*n2-n3*((*(ryy_j_i+k))+(*(ryy_j_i_2+k)))-n4*((*(ryy_j_i_2+k))+(*(ryy_j_i_3+k)))-n5*((*(ryy_j_i_3+k))+(*(ryy_j_i_4+k)))-(e*vdiag)+(d*vxxzz_T2));
                                *(rzz_j_i_4+k)=n1*((*(rzz_j_i+k))*n2-n3*((*(rzz_j_i+k))+(*(rzz_j_i_2+k)))-n4*((*(rzz_j_i_2+k))+(*(rzz_j_i_3+k)))-n5*((*(rzz_j_i_3+k))+(*(rzz_j_i_4+k)))-(e*vdiag)+(d*vxxyy_T2));
                                
                                
                                sumrxy=c1*(*(rxy_j_i_4+k))+c2*(*(rxy_j_i+k))+c3*(*(rxy_j_i_2+k))+c4*(*(rxy_j_i_3+k));
                                sumryz=c1*(*(ryz_j_i_4+k))+c2*(*(ryz_j_i+k))+c3*(*(ryz_j_i_2+k))+c4*(*(ryz_j_i_3+k));
                                sumrxz=c1*(*(rxz_j_i_4+k))+c2*(*(rxz_j_i+k))+c3*(*(rxz_j_i_2+k))+c4*(*(rxz_j_i_3+k));
                                sumrxx=c1*(*(rxx_j_i_4+k))+c2*(*(rxx_j_i+k))+c3*(*(rxx_j_i_2+k))+c4*(*(rxx_j_i_3+k));
                                sumryy=c1*(*(ryy_j_i_4+k))+c2*(*(ryy_j_i+k))+c3*(*(ryy_j_i_2+k))+c4*(*(ryy_j_i_3+k));
                                sumrzz=c1*(*(rzz_j_i_4+k))+c2*(*(rzz_j_i+k))+c3*(*(rzz_j_i_2+k))+c4*(*(rzz_j_i_3+k));
                                
                                
                                /* and now the components of the stress tensor are
                                 completely updated */
                                sxy[j][i][k]+=(dthalbe*sumrxy);
                                syz[j][i][k]+=(dthalbe*sumryz);
                                sxz[j][i][k]+=(dthalbe*sumrxz);
                                sxx[j][i][k]+=(dthalbe*sumrxx);
                                syy[j][i][k]+=(dthalbe*sumryy);
                                szz[j][i][k]+=(dthalbe*sumrzz);
                                
                                
                                
                            }
                        }
                    }
                    break;
                    
                case 12 :
                    
                    /* Taylor coefficients */
                    b1=160083.0/131072.0; b2=-12705.0/131072.0; b3=22869.0/1310720.0;
                    b4=-5445.0/1835008.0; b5=847.0/2359296.0; b6=-63.0/2883584;
                    
                    /* Holberg coefficients E=0.1 %*/
                    if(FDCOEFF==2){
                        b1=1.2508; b2=-0.12034; b3=0.032131; b4=-0.010142; b5=0.0029857; b6=-0.00066667;}
                    
                    
                    for (j=ny1;j<=ny2;j++){
                        
                        rxx_j=*(rxx+j);ryy_j=*(ryy+j);rzz_j=*(rzz+j);rxy_j=*(rxy+j);ryz_j=*(ryz+j);rxz_j=*(rxz+j);
                        rxx_j_2=*(rxx_2+j);ryy_j_2=*(ryy_2+j);rzz_j_2=*(rzz_2+j);rxy_j_2=*(rxy_2+j);ryz_j_2=*(ryz_2+j); rxz_j_2=*(rxz_2+j);
                        rxx_j_3=*(rxx_3+j);ryy_j_3=*(ryy_3+j);rzz_j_3=*(rzz_3+j);rxy_j_3=*(rxy_3+j);ryz_j_3=*(ryz_3+j);rxz_j_3=*(rxz_3+j);
                        rxx_j_4=*(rxx_4+j);ryy_j_4=*(ryy_4+j);rzz_j_4=*(rzz_4+j);rxy_j_4=*(rxy_4+j);ryz_j_4=*(ryz_4+j);rxz_j_4=*(rxz_4+j);
                        vxyyx_j=*(vxyyx+j);vyzzy_j=*(vyzzy+j);vxzzx_j=*(vxzzx+j);vxxyyzz_j=*(vxxyyzz+j);vyyzz_j=*(vyyzz+j);vxxzz_j=*(vxxzz+j);vxxyy_j=*(vxxyy+j);
                        vxyyx_j_2=*(vxyyx_2+j);vyzzy_j_2=*(vyzzy_2+j); vxzzx_j_2=*(vxzzx_2+j);vxxyyzz_j_2=*(vxxyyzz_2+j);vyyzz_j_2=*(vyyzz_2+j);vxxzz_j_2=*(vxxzz_2+j);vxxyy_j_2=*(vxxyy_2+j);
                        vxyyx_j_3=*(vxyyx_3+j);vyzzy_j_3=*(vyzzy_3+j);vxzzx_j_3=*(vxzzx_3+j);vxxyyzz_j_3=*(vxxyyzz_3+j);vyyzz_j_3=*(vyyzz_3+j);vxxzz_j_3=*(vxxzz_3+j);vxxyy_j_3=*(vxxyy_3+j);
                        vxyyx_j_4=*(vxyyx_4+j);vyzzy_j_4=*(vyzzy_4+j);vxzzx_j_4=*(vxzzx_4+j);vxxyyzz_j_4=*(vxxyyzz_4+j);vyyzz_j_4=*(vyyzz_4+j);vxxzz_j_4=*(vxxzz_4+j);vxxyy_j_4=*(vxxyy_4+j);
                        
                        for (i=nx1;i<=nx2;i++){
                            
                            rxx_j_i=*(rxx_j+i);ryy_j_i=*(ryy_j+i);rzz_j_i=*(rzz_j+i);rxy_j_i=*(rxy_j+i);ryz_j_i=*(ryz_j+i);rxz_j_i=*(rxz_j+i);
                            rxx_j_i_2=*(rxx_j_2+i);ryy_j_i_2=*(ryy_j_2+i);rzz_j_i_2=*(rzz_j_2+i);rxy_j_i_2=*(rxy_j_2+i);ryz_j_i_2=*(ryz_j_2+i);rxz_j_i_2=*(rxz_j_2+i);
                            rxx_j_i_3=*(rxx_j_3+i);ryy_j_i_3=*(ryy_j_3+i);rzz_j_i_3=*(rzz_j_3+i);rxy_j_i_3=*(rxy_j_3+i);ryz_j_i_3=*(ryz_j_3+i);rxz_j_i_3=*(rxz_j_3+i);
                            rxx_j_i_4=*(rxx_j_4+i);ryy_j_i_4=*(ryy_j_4+i);rzz_j_i_4=*(rzz_j_4+i);rxy_j_i_4=*(rxy_j_4+i);ryz_j_i_4=*(ryz_j_4+i);rxz_j_i_4=*(rxz_j_4+i);
                            vxyyx_j_i=*(vxyyx_j+i);vyzzy_j_i=*(vyzzy_j+i);vxzzx_j_i=*(vxzzx_j+i);vxxyyzz_j_i=*(vxxyyzz_j+i);vyyzz_j_i=*(vyyzz_j+i);vxxzz_j_i=*(vxxzz_j+i);vxxyy_j_i=*(vxxyy_j+i);
                            vxyyx_j_i_2=*(vxyyx_j_2+i);vyzzy_j_i_2=*(vyzzy_j_2+i);vxzzx_j_i_2=*(vxzzx_j_2+i);vxxyyzz_j_i_2=*(vxxyyzz_j_2+i);vyyzz_j_i_2=*(vyyzz_j_2+i);vxxzz_j_i_2=*(vxxzz_j_2+i);vxxyy_j_i_2=*(vxxyy_j_2+i);
                            vxyyx_j_i_3=*(vxyyx_j_3+i);vyzzy_j_i_3=*(vyzzy_j_3+i);vxzzx_j_i_3=*(vxzzx_j_3+i); vxxyyzz_j_i_3=*(vxxyyzz_j_3+i);vyyzz_j_i_3=*(vyyzz_j_3+i);vxxzz_j_i_3=*(vxxzz_j_3+i); vxxyy_j_i_3=*(vxxyy_j_3+i);
                            vxyyx_j_i_4=*(vxyyx_j_4+i);vyzzy_j_i_4=*(vyzzy_j_4+i);vxzzx_j_i_4=*(vxzzx_j_4+i); vxxyyzz_j_i_4=*(vxxyyzz_j_4+i);vyyzz_j_i_4=*(vyyzz_j_4+i);vxxzz_j_i_4=*(vxxzz_j_4+i);vxxyy_j_i_4=*(vxxyy_j_4+i);
                            
                            for (k=nz1;k<=nz2;k++){
                                /* spatial derivatives of the components of the velocities
                                 are computed */
                                
                                vxx = (b1*(vx[j][i][k]-vx[j][i-1][k])+
                                       b2*(vx[j][i+1][k]-vx[j][i-2][k])+
                                       b3*(vx[j][i+2][k]-vx[j][i-3][k])+
                                       b4*(vx[j][i+3][k]-vx[j][i-4][k])+
                                       b5*(vx[j][i+4][k]-vx[j][i-5][k])+
                                       b6*(vx[j][i+5][k]-vx[j][i-6][k]))/DX;
                                
                                vxy = (b1*(vx[j+1][i][k]-vx[j][i][k])+
                                       b2*(vx[j+2][i][k]-vx[j-1][i][k])+
                                       b3*(vx[j+3][i][k]-vx[j-2][i][k])+
                                       b4*(vx[j+4][i][k]-vx[j-3][i][k])+
                                       b5*(vx[j+5][i][k]-vx[j-4][i][k])+
                                       b6*(vx[j+6][i][k]-vx[j-5][i][k]))/DY;
                                
                                vxz = (b1*(vx[j][i][k+1]-vx[j][i][k])+
                                       b2*(vx[j][i][k+2]-vx[j][i][k-1])+
                                       b3*(vx[j][i][k+3]-vx[j][i][k-2])+
                                       b4*(vx[j][i][k+4]-vx[j][i][k-3])+
                                       b5*(vx[j][i][k+5]-vx[j][i][k-4])+
                                       b6*(vx[j][i][k+6]-vx[j][i][k-5]))/DZ;
                                
                                vyx = (b1*(vy[j][i+1][k]-vy[j][i][k])+
                                       b2*(vy[j][i+2][k]-vy[j][i-1][k])+
                                       b3*(vy[j][i+3][k]-vy[j][i-2][k])+
                                       b4*(vy[j][i+4][k]-vy[j][i-3][k])+
                                       b5*(vy[j][i+5][k]-vy[j][i-4][k])+
                                       b6*(vy[j][i+6][k]-vy[j][i-5][k]))/DX;
                                
                                vyy = (b1*(vy[j][i][k]-vy[j-1][i][k])+
                                       b2*(vy[j+1][i][k]-vy[j-2][i][k])+
                                       b3*(vy[j+2][i][k]-vy[j-3][i][k])+
                                       b4*(vy[j+3][i][k]-vy[j-4][i][k])+
                                       b5*(vy[j+4][i][k]-vy[j-5][i][k])+
                                       b6*(vy[j+5][i][k]-vy[j-6][i][k]))/DY;
                                
                                vyz = (b1*(vy[j][i][k+1]-vy[j][i][k])+
                                       b2*(vy[j][i][k+2]-vy[j][i][k-1])+
                                       b3*(vy[j][i][k+3]-vy[j][i][k-2])+
                                       b4*(vy[j][i][k+4]-vy[j][i][k-3])+
                                       b5*(vy[j][i][k+5]-vy[j][i][k-4])+
                                       b6*(vy[j][i][k+6]-vy[j][i][k-5]))/DZ;
                                
                                vzx = (b1*(vz[j][i+1][k]-vz[j][i][k])+
                                       b2*(vz[j][i+2][k]-vz[j][i-1][k])+
                                       b3*(vz[j][i+3][k]-vz[j][i-2][k])+
                                       b4*(vz[j][i+4][k]-vz[j][i-3][k])+
                                       b5*(vz[j][i+5][k]-vz[j][i-4][k])+
                                       b6*(vz[j][i+6][k]-vz[j][i-5][k]))/DX;
                                
                                vzy = (b1*(vz[j+1][i][k]-vz[j][i][k])+
                                       b2*(vz[j+2][i][k]-vz[j-1][i][k])+
                                       b3*(vz[j+3][i][k]-vz[j-2][i][k])+
                                       b4*(vz[j+4][i][k]-vz[j-3][i][k])+
                                       b5*(vz[j+5][i][k]-vz[j-4][i][k])+
                                       b6*(vz[j+6][i][k]-vz[j-5][i][k]))/DY;
                                
                                vzz = (b1*(vz[j][i][k]-vz[j][i][k-1])+
                                       b2*(vz[j][i][k+1]-vz[j][i][k-2])+
                                       b3*(vz[j][i][k+2]-vz[j][i][k-3])+
                                       b4*(vz[j][i][k+3]-vz[j][i][k-4])+
                                       b5*(vz[j][i][k+4]-vz[j][i][k-5])+
                                       b6*(vz[j][i][k+5]-vz[j][i][k-6]))/DZ;
                                
                                
                                /* computing sums of the old memory variables in this version only one mechanism is possi
                                 ble */
                                
                                sumrxy=c1*(*(rxy_j_i+k))+c2*(*(rxy_j_i_2+k))+c3*(*(rxy_j_i_3+k))+c4*(*(rxy_j_i_4+k));
                                sumryz=c1*(*(ryz_j_i+k))+c2*(*(ryz_j_i_2+k))+c3*(*(ryz_j_i_3+k))+c4*(*(ryz_j_i_4+k));
                                sumrxz=c1*(*(rxz_j_i+k))+c2*(*(rxz_j_i_2+k))+c3*(*(rxz_j_i_3+k))+c4*(*(rxz_j_i_4+k));
                                sumrxx=c1*(*(rxx_j_i+k))+c2*(*(rxx_j_i_2+k))+c3*(*(rxx_j_i_3+k))+c4*(*(rxx_j_i_4+k));
                                sumryy=c1*(*(ryy_j_i+k))+c2*(*(ryy_j_i_2+k))+c3*(*(ryy_j_i_3+k))+c4*(*(ryy_j_i_4+k));
                                sumrzz=c1*(*(rzz_j_i+k))+c2*(*(rzz_j_i_2+k))+c3*(*(rzz_j_i_3+k))+c4*(*(rzz_j_i_4+k));
                                
                                /* updating components of the stress tensor, partially */
                                fipjp=uipjp[j][i][k]*DT*(1.0+L*tausipjp[j][i][k]);
                                fjpkp=ujpkp[j][i][k]*DT*(1.0+L*tausjpkp[j][i][k]);
                                fipkp=uipkp[j][i][k]*DT*(1.0+L*tausipkp[j][i][k]);
                                g=pi[j][i][k]*(1.0+L*taup[j][i][k]);
                                f=2.0*u[j][i][k]*(1.0+L*taus[j][i][k]);
                                
                                /* Save spatial derivations */
                                *(vxyyx_j_i+k)=vxy+vyx;
                                *(vyzzy_j_i+k)=vyz+vzy;
                                *(vxzzx_j_i+k)=vxz+vzx;
                                *(vxxyyzz_j_i+k)=vxx+vyy+vzz;
                                *(vyyzz_j_i+k)=vyy+vzz;
                                *(vxxzz_j_i+k)=vxx+vzz;
                                *(vxxyy_j_i+k)=vxx+vyy;
                                
                                /* Do the Adams Bashforth stuff and save it, we need it for the memory variables again */
                                vxyyx_T2=(c1*(*(vxyyx_j_i+k))+c2*(*(vxyyx_j_i_2+k))+c3*(*(vxyyx_j_i_3+k))+c4*(*(vxyyx_j_i_4+k)));
                                vyzzy_T2=(c1*(*(vyzzy_j_i+k))+c2*(*(vyzzy_j_i_2+k))+c3*(*(vyzzy_j_i_3+k))+c4*(*(vyzzy_j_i_4+k)));
                                vxzzx_T2=(c1*(*(vxzzx_j_i+k))+c2*(*(vxzzx_j_i_2+k))+c3*(*(vxzzx_j_i_3+k))+c4*(*(vxzzx_j_i_4+k)));
                                vyyzz_T2=(c1*(*(vyyzz_j_i+k))+c2*(*(vyyzz_j_i_2+k))+c3*(*(vyyzz_j_i_3+k))+c4*(*(vyyzz_j_i_4+k)));
                                vxxzz_T2=(c1*(*(vxxzz_j_i+k))+c2*(*(vxxzz_j_i_2+k))+c3*(*(vxxzz_j_i_3+k))+c4*(*(vxxzz_j_i_4+k)));
                                vxxyy_T2=(c1*(*(vxxyy_j_i+k))+c2*(*(vxxyy_j_i_2+k))+c3*(*(vxxyy_j_i_3+k))+c4*(*(vxxyy_j_i_4+k)));
                                vdiag=((c1*(*(vxxyyzz_j_i+k))+c2*(*(vxxyyzz_j_i_2+k))+c3*(*(vxxyyzz_j_i_3+k)))+c4*(*(vxxyyzz_j_i_4+k)));
                                
                                /* Update the stress tensor */
                                sxy[j][i][k]+=(fipjp*vxyyx_T2)+(dthalbe*sumrxy);
                                syz[j][i][k]+=(fjpkp*vyzzy_T2)+(dthalbe*sumryz);
                                sxz[j][i][k]+=(fipkp*vxzzx_T2)+(dthalbe*sumrxz);
                                sxx[j][i][k]+=DT*((g*vdiag)-(f*vyyzz_T2))+(dthalbe*sumrxx);
                                syy[j][i][k]+=DT*((g*vdiag)-(f*vxxzz_T2))+(dthalbe*sumryy);
                                szz[j][i][k]+=DT*((g*vdiag)-(f*vxxyy_T2))+(dthalbe*sumrzz);
                                
                                /* Update the memory variables */
                                
                                /* in the case of L=1 update the memory-variables without suming*/
                                l=1;
                                n1=1.0/(1.0+(c1*eta[l]*0.5)); /* 1/(1+DT*c1/(2*tau_sigl)) */
                                n2=1.0-(c1*eta[l]*0.5);
                                n3=c2*eta[l]/2.0;
                                n4=c3*eta[l]/2.0;
                                n5=c4*eta[l]/2.0;
                                dipjp=uipjp[j][i][k]*eta[l]*tausipjp[j][i][k];
                                djpkp=ujpkp[j][i][k]*eta[l]*tausjpkp[j][i][k];
                                dipkp=uipkp[j][i][k]*eta[l]*tausipkp[j][i][k];
                                d=2.0*u[j][i][k]*eta[l]*taus[j][i][k];
                                e=pi[j][i][k]*eta[l]*taup[j][i][k];
                                
                                
                                *(rxy_j_i_4+k)=n1*((*(rxy_j_i+k))*n2-n3*((*(rxy_j_i+k))+(*(rxy_j_i_2+k)))-n4*((*(rxy_j_i_2+k))+(*(rxy_j_i_3+k)))-n5*((*(rxy_j_i_3+k))+(*(rxy_j_i_4+k)))-(dipjp*vxyyx_T2));
                                *(ryz_j_i_4+k)=n1*((*(ryz_j_i+k))*n2-n3*((*(ryz_j_i+k))+(*(ryz_j_i_2+k)))-n4*((*(ryz_j_i_2+k))+(*(ryz_j_i_3+k)))-n5*((*(ryz_j_i_3+k))+(*(ryz_j_i_4+k)))-(djpkp*vyzzy_T2));
                                *(rxz_j_i_4+k)=n1*((*(rxz_j_i+k))*n2-n3*((*(rxz_j_i+k))+(*(rxz_j_i_2+k)))-n4*((*(rxz_j_i_2+k))+(*(rxz_j_i_3+k)))-n5*((*(rxz_j_i_3+k))+(*(rxz_j_i_4+k)))-(dipkp*vxzzx_T2));
                                *(rxx_j_i_4+k)=n1*((*(rxx_j_i+k))*n2-n3*((*(rxx_j_i+k))+(*(rxx_j_i_2+k)))-n4*((*(rxx_j_i_2+k))+(*(rxx_j_i_3+k)))-n5*((*(rxx_j_i_3+k))+(*(rxx_j_i_4+k)))-(e*vdiag)+(d*vyyzz_T2));
                                *(ryy_j_i_4+k)=n1*((*(ryy_j_i+k))*n2-n3*((*(ryy_j_i+k))+(*(ryy_j_i_2+k)))-n4*((*(ryy_j_i_2+k))+(*(ryy_j_i_3+k)))-n5*((*(ryy_j_i_3+k))+(*(ryy_j_i_4+k)))-(e*vdiag)+(d*vxxzz_T2));
                                *(rzz_j_i_4+k)=n1*((*(rzz_j_i+k))*n2-n3*((*(rzz_j_i+k))+(*(rzz_j_i_2+k)))-n4*((*(rzz_j_i_2+k))+(*(rzz_j_i_3+k)))-n5*((*(rzz_j_i_3+k))+(*(rzz_j_i_4+k)))-(e*vdiag)+(d*vxxyy_T2));
                                
                                
                                sumrxy=c1*(*(rxy_j_i_4+k))+c2*(*(rxy_j_i+k))+c3*(*(rxy_j_i_2+k))+c4*(*(rxy_j_i_3+k));
                                sumryz=c1*(*(ryz_j_i_4+k))+c2*(*(ryz_j_i+k))+c3*(*(ryz_j_i_2+k))+c4*(*(ryz_j_i_3+k));
                                sumrxz=c1*(*(rxz_j_i_4+k))+c2*(*(rxz_j_i+k))+c3*(*(rxz_j_i_2+k))+c4*(*(rxz_j_i_3+k));
                                sumrxx=c1*(*(rxx_j_i_4+k))+c2*(*(rxx_j_i+k))+c3*(*(rxx_j_i_2+k))+c4*(*(rxx_j_i_3+k));
                                sumryy=c1*(*(ryy_j_i_4+k))+c2*(*(ryy_j_i+k))+c3*(*(ryy_j_i_2+k))+c4*(*(ryy_j_i_3+k));
                                sumrzz=c1*(*(rzz_j_i_4+k))+c2*(*(rzz_j_i+k))+c3*(*(rzz_j_i_2+k))+c4*(*(rzz_j_i_3+k));
                                
                                
                                /* and now the components of the stress tensor are
                                 completely updated */
                                sxy[j][i][k]+=(dthalbe*sumrxy);
                                syz[j][i][k]+=(dthalbe*sumryz);
                                sxz[j][i][k]+=(dthalbe*sumrxz);
                                sxx[j][i][k]+=(dthalbe*sumrxx);
                                syy[j][i][k]+=(dthalbe*sumryy);
                                szz[j][i][k]+=(dthalbe*sumrzz);
                                

                                
                            }
                        }
                    }
                    break;
                    
            }
            break; /* Break FDORDER_TIME=4 */
    }
    
    

	if (LOG)
		if ((MYID==0) && ((nt+(OUTNTIMESTEPINFO-1))%OUTNTIMESTEPINFO)==0){
			time2=MPI_Wtime();
			time=time2-time1;
			fprintf(FP," Real time for stress tensor update: \t\t %4.2f s.\n",time);
		}
	return time;

}
