/*-----------------------------------------------------------------------------------------
 * Copyright (C) 2015  For the list of authors, see file AUTHORS.
 *
 * This file is part of SOFI2D.
 *
 * SOFI2D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 *
 * SOFI2D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SOFI2D. See file COPYING and/or <http://www.gnu.org/licenses/gpl-2.0.html>.
 -----------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
 *   updating stress components at interior gridpoints (excluding boundarys) [gx2+1...gx3][gy2+1...gy3]
 *   by a staggered grid finite difference scheme of arbitrary (FDORDER) order accuracy in space
 *   and second order accuracy in time
 *   T. Bohlen
 *
 *   gx and gy are arrays with the locations of the boundary specified in subgrid_bounds.c
 *   for each subgrid
 *  ----------------------------------------------------------------------*/


#include "fd.h"

void update_s_elastic_interior_4 ( int nx1, int nx2, int ny1, int ny2, int * gx, int * gy, int nt,
                                float **  vx, float **   vy, float **   sxx, float **   syy,
                                float **   sxy, float ** pi, float ** u, float ** uipjp, float *hc,float ** vxx_1,float ** vxx_2,float ** vxx_3,float ** vxx_4,float ** vyy_1,float ** vyy_2,float ** vyy_3,float ** vyy_4,float ** vxy_1,float ** vxy_2,float ** vxy_3,float ** vxy_4,float ** vyx_1,float ** vyx_2,float ** vyx_3,float ** vyx_4)
{
    
    
    int i,j;
    float fipjp, f, g;
    float  dhi;
    extern float DT, DH;
    extern int MYID, FDORDER;
    extern FILE *FP;
    extern int OUTNTIMESTEPINFO;
    double time1=0.0, time2=0.0;
    float c1, c2, c3, c4; /* Coefficients for Adam Bashforth */
    c1=13.0/12.0; c2=-5.0/24.0; c3=1.0/6.0; c4=-1.0/24.0;
    
    float * vxx_1_j, *vyy_1_j,* vxy_1_j, *vyx_1_j;
    float * vxx_2_j, *vyy_2_j,* vxy_2_j, *vyx_2_j;
    float * vxx_3_j, *vyy_3_j,* vxy_3_j, *vyx_3_j;
    float * vxx_4_j, *vyy_4_j,* vxy_4_j, *vyx_4_j;
    
    float sumxx=0.0,sumyy=0.0,sumxy=0.0,sumyx=0.0;
    
    dhi = 1.0/DH;
    
    if ( ( MYID==0 ) && ( ( nt+ ( OUTNTIMESTEPINFO-1 ) ) %OUTNTIMESTEPINFO ) ==0 ) {
        time1=MPI_Wtime();
        fprintf ( FP,"\n **Message from update_s_interior_4 (printed by PE %d):\n",MYID );
        fprintf ( FP," Updating stress components ..." );
    }
    
    
    
    switch ( FDORDER ) {
            
        case 2:
            for ( j=gy[2]+1; j<=gy[3]; j++ ) {
                vxx_1_j=*(vxx_1+j),vyy_1_j=*(vyy_1+j),vxy_1_j=*(vxy_1+j),vyx_1_j=*(vyx_1+j);
                vxx_2_j=*(vxx_2+j),vyy_2_j=*(vyy_2+j),vxy_2_j=*(vxy_2+j),vyx_2_j=*(vyx_2+j);
                vxx_3_j=*(vxx_3+j),vyy_3_j=*(vyy_3+j),vxy_3_j=*(vxy_3+j),vyx_3_j=*(vyx_3+j);
                vxx_4_j=*(vxx_4+j),vyy_4_j=*(vyy_4+j),vxy_4_j=*(vxy_4+j),vyx_4_j=*(vyx_4+j);
                for ( i=gx[2]+1; i<=gx[3]; i++ ) {
    
                    
                    /* spatial derivatives of the components of the velocities */
                    /* using Holberg coefficients */
                    *(vxx_1_j+i) = hc[1]* ( vx[j][i]  -vx[j][i-1] );
                    *(vyy_1_j+i) = hc[1]* ( vy[j][i]  -vy[j-1][i] );
                    *(vyx_1_j+i) = hc[1]* ( vy[j][i+1]-vy[j][i] );
                    *(vxy_1_j+i) = hc[1]* ( vx[j+1][i]-vx[j][i] );
                    
                    /* updating components of the stress tensor, partially */
                    fipjp=uipjp[j][i]*DT;
                    f=u[j][i]*DT;
                    g=pi[j][i]*DT;
                    
                    // Calculate Adams-Bashforth stuff
                    sumxx=c1*(*(vxx_1_j+i))+c2*(*(vxx_2_j+i))+c3*(*(vxx_3_j+i))+c4*(*(vxx_4_j+i));
                    sumyy=c1*(*(vyy_1_j+i))+c2*(*(vyy_2_j+i))+c3*(*(vyy_3_j+i))+c4*(*(vyy_4_j+i));
                    sumxy=c1*(*(vxy_1_j+i))+c2*(*(vxy_2_j+i))+c3*(*(vxy_3_j+i))+c4*(*(vxy_4_j+i));
                    sumyx=c1*(*(vyx_1_j+i))+c2*(*(vyx_2_j+i))+c3*(*(vyx_3_j+i))+c4*(*(vyx_4_j+i));

                    // Update stress
                    sxy[j][i]+= ( fipjp* ( sumxy+sumyx ) )*dhi;
                    sxx[j][i]+= (( g* ( sumxx+sumyy ) )- ( 2.0*f*sumyy ))*dhi;
                    syy[j][i]+= (( g* ( sumxx+sumyy ) )- ( 2.0*f*sumxx ))*dhi;
                }
            }
            break;
            
        case 4:
            for ( j=gy[2]+1; j<=gy[3]; j++ ) {
                vxx_1_j=*(vxx_1+j),vyy_1_j=*(vyy_1+j),vxy_1_j=*(vxy_1+j),vyx_1_j=*(vyx_1+j);
                vxx_2_j=*(vxx_2+j),vyy_2_j=*(vyy_2+j),vxy_2_j=*(vxy_2+j),vyx_2_j=*(vyx_2+j);
                vxx_3_j=*(vxx_3+j),vyy_3_j=*(vyy_3+j),vxy_3_j=*(vxy_3+j),vyx_3_j=*(vyx_3+j);
                vxx_4_j=*(vxx_4+j),vyy_4_j=*(vyy_4+j),vxy_4_j=*(vxy_4+j),vyx_4_j=*(vyx_4+j);
                for ( i=gx[2]+1; i<=gx[3]; i++ ) {
                    *(vxx_1_j+i)= ( hc[1]* ( vx[j][i]  -vx[j][i-1] )
                           + hc[2]* ( vx[j][i+1]-vx[j][i-2] )
                                   );
                    *(vyy_1_j+i) = ( hc[1]* ( vy[j][i]  -vy[j-1][i] )
                           + hc[2]* ( vy[j+1][i]-vy[j-2][i] )
                                    );
                    *(vyx_1_j+i)= ( hc[1]* ( vy[j][i+1]-vy[j][i] )
                           + hc[2]* ( vy[j][i+2]-vy[j][i-1] )
                                   );
                    *(vxy_1_j+i) = ( hc[1]* ( vx[j+1][i]-vx[j][i] )
                           + hc[2]* ( vx[j+2][i]-vx[j-1][i] )
                                    );
                    
                    fipjp=uipjp[j][i]*DT;
                    f=u[j][i]*DT;
                    g=pi[j][i]*DT;
                    
                    // Calculate Adams-Bashforth stuff
                    sumxx=c1*(*(vxx_1_j+i))+c2*(*(vxx_2_j+i))+c3*(*(vxx_3_j+i))+c4*(*(vxx_4_j+i));
                    sumyy=c1*(*(vyy_1_j+i))+c2*(*(vyy_2_j+i))+c3*(*(vyy_3_j+i))+c4*(*(vyy_4_j+i));
                    sumxy=c1*(*(vxy_1_j+i))+c2*(*(vxy_2_j+i))+c3*(*(vxy_3_j+i))+c4*(*(vxy_4_j+i));
                    sumyx=c1*(*(vyx_1_j+i))+c2*(*(vyx_2_j+i))+c3*(*(vyx_3_j+i))+c4*(*(vyx_4_j+i));
                    
                    // Update stress
                    sxy[j][i]+= ( fipjp* ( sumxy+sumyx ) )*dhi;
                    sxx[j][i]+= (( g* ( sumxx+sumyy ) )- ( 2.0*f*sumyy ))*dhi;
                    syy[j][i]+= (( g* ( sumxx+sumyy ) )- ( 2.0*f*sumxx ))*dhi;

                }
            }
            break;
            
        case 6:
            for ( j=gy[2]+1; j<=gy[3]; j++ ) {
                vxx_1_j=*(vxx_1+j),vyy_1_j=*(vyy_1+j),vxy_1_j=*(vxy_1+j),vyx_1_j=*(vyx_1+j);
                vxx_2_j=*(vxx_2+j),vyy_2_j=*(vyy_2+j),vxy_2_j=*(vxy_2+j),vyx_2_j=*(vyx_2+j);
                vxx_3_j=*(vxx_3+j),vyy_3_j=*(vyy_3+j),vxy_3_j=*(vxy_3+j),vyx_3_j=*(vyx_3+j);
                vxx_4_j=*(vxx_4+j),vyy_4_j=*(vyy_4+j),vxy_4_j=*(vxy_4+j),vyx_4_j=*(vyx_4+j);
                for ( i=gx[2]+1; i<=gx[3]; i++ ) {
                    *(vxx_1_j+i) = ( hc[1]* ( vx[j][i]  -vx[j][i-1] )
                           + hc[2]* ( vx[j][i+1]-vx[j][i-2] )
                           + hc[3]* ( vx[j][i+2]-vx[j][i-3] )
                           );
                    *(vyy_1_j+i) = ( hc[1]* ( vy[j][i]  -vy[j-1][i] )
                           + hc[2]* ( vy[j+1][i]-vy[j-2][i] )
                           + hc[3]* ( vy[j+2][i]-vy[j-3][i] )
                           );
                    *(vyx_1_j+i) = ( hc[1]* ( vy[j][i+1]-vy[j][i] )
                           + hc[2]* ( vy[j][i+2]-vy[j][i-1] )
                           + hc[3]* ( vy[j][i+3]-vy[j][i-2] )
                           );
                    *(vxy_1_j+i) = ( hc[1]* ( vx[j+1][i]-vx[j][i] )
                           + hc[2]* ( vx[j+2][i]-vx[j-1][i] )
                           + hc[3]* ( vx[j+3][i]-vx[j-2][i] )
                           );
                    
                    fipjp=uipjp[j][i]*DT;
                    f=u[j][i]*DT;
                    g=pi[j][i]*DT;
                    
                    // Calculate Adams-Bashforth stuff
                    sumxx=c1*(*(vxx_1_j+i))+c2*(*(vxx_2_j+i))+c3*(*(vxx_3_j+i))+c4*(*(vxx_4_j+i));
                    sumyy=c1*(*(vyy_1_j+i))+c2*(*(vyy_2_j+i))+c3*(*(vyy_3_j+i))+c4*(*(vyy_4_j+i));
                    sumxy=c1*(*(vxy_1_j+i))+c2*(*(vxy_2_j+i))+c3*(*(vxy_3_j+i))+c4*(*(vxy_4_j+i));
                    sumyx=c1*(*(vyx_1_j+i))+c2*(*(vyx_2_j+i))+c3*(*(vyx_3_j+i))+c4*(*(vyx_4_j+i));
                    
                    // Update stress
                    sxy[j][i]+= ( fipjp* ( sumxy+sumyx ) )*dhi;
                    sxx[j][i]+= (( g* ( sumxx+sumyy ) )- ( 2.0*f*sumyy ))*dhi;
                    syy[j][i]+= (( g* ( sumxx+sumyy ) )- ( 2.0*f*sumxx ))*dhi;

                }
            }
            break;
            
        case 8:
            for ( j=gy[2]+1; j<=gy[3]; j++ ) {
                vxx_1_j=*(vxx_1+j),vyy_1_j=*(vyy_1+j),vxy_1_j=*(vxy_1+j),vyx_1_j=*(vyx_1+j);
                vxx_2_j=*(vxx_2+j),vyy_2_j=*(vyy_2+j),vxy_2_j=*(vxy_2+j),vyx_2_j=*(vyx_2+j);
                vxx_3_j=*(vxx_3+j),vyy_3_j=*(vyy_3+j),vxy_3_j=*(vxy_3+j),vyx_3_j=*(vyx_3+j);
                vxx_4_j=*(vxx_4+j),vyy_4_j=*(vyy_4+j),vxy_4_j=*(vxy_4+j),vyx_4_j=*(vyx_4+j);
                for ( i=gx[2]+1; i<=gx[3]; i++ ) {
                    *(vxx_1_j+i) = ( hc[1]* ( vx[j][i]  -vx[j][i-1] )
                           + hc[2]* ( vx[j][i+1]-vx[j][i-2] )
                           + hc[3]* ( vx[j][i+2]-vx[j][i-3] )
                           + hc[4]* ( vx[j][i+3]-vx[j][i-4] )
                           );
                    *(vyy_1_j+i) = ( hc[1]* ( vy[j][i]  -vy[j-1][i] )
                           + hc[2]* ( vy[j+1][i]-vy[j-2][i] )
                           + hc[3]* ( vy[j+2][i]-vy[j-3][i] )
                           + hc[4]* ( vy[j+3][i]-vy[j-4][i] )
                           );
                    *(vyx_1_j+i) = ( hc[1]* ( vy[j][i+1]-vy[j][i] )
                           + hc[2]* ( vy[j][i+2]-vy[j][i-1] )
                           + hc[3]* ( vy[j][i+3]-vy[j][i-2] )
                           + hc[4]* ( vy[j][i+4]-vy[j][i-3] )
                           );
                    *(vxy_1_j+i) = ( hc[1]* ( vx[j+1][i]-vx[j][i] )
                           + hc[2]* ( vx[j+2][i]-vx[j-1][i] )
                           + hc[3]* ( vx[j+3][i]-vx[j-2][i] )
                           + hc[4]* ( vx[j+4][i]-vx[j-3][i] )
                           );
                    
                    fipjp=uipjp[j][i]*DT;
                    f=u[j][i]*DT;
                    g=pi[j][i]*DT;
                    
                    // Calculate Adams-Bashforth stuff
                    sumxx=c1*(*(vxx_1_j+i))+c2*(*(vxx_2_j+i))+c3*(*(vxx_3_j+i))+c4*(*(vxx_4_j+i));
                    sumyy=c1*(*(vyy_1_j+i))+c2*(*(vyy_2_j+i))+c3*(*(vyy_3_j+i))+c4*(*(vyy_4_j+i));
                    sumxy=c1*(*(vxy_1_j+i))+c2*(*(vxy_2_j+i))+c3*(*(vxy_3_j+i))+c4*(*(vxy_4_j+i));
                    sumyx=c1*(*(vyx_1_j+i))+c2*(*(vyx_2_j+i))+c3*(*(vyx_3_j+i))+c4*(*(vyx_4_j+i));
                    
                    // Update stress
                    sxy[j][i]+= ( fipjp* ( sumxy+sumyx ) )*dhi;
                    sxx[j][i]+= (( g* ( sumxx+sumyy ) )- ( 2.0*f*sumyy ))*dhi;
                    syy[j][i]+= (( g* ( sumxx+sumyy ) )- ( 2.0*f*sumxx ))*dhi;

                }
            }
            break;
            
        case 10:
            for ( j=gy[2]+1; j<=gy[3]; j++ ) {
                vxx_1_j=*(vxx_1+j),vyy_1_j=*(vyy_1+j),vxy_1_j=*(vxy_1+j),vyx_1_j=*(vyx_1+j);
                vxx_2_j=*(vxx_2+j),vyy_2_j=*(vyy_2+j),vxy_2_j=*(vxy_2+j),vyx_2_j=*(vyx_2+j);
                vxx_3_j=*(vxx_3+j),vyy_3_j=*(vyy_3+j),vxy_3_j=*(vxy_3+j),vyx_3_j=*(vyx_3+j);
                vxx_4_j=*(vxx_4+j),vyy_4_j=*(vyy_4+j),vxy_4_j=*(vxy_4+j),vyx_4_j=*(vyx_4+j);
                for ( i=gx[2]+1; i<=gx[3]; i++ ) {
                    *(vxx_1_j+i) = ( hc[1]* ( vx[j][i]  -vx[j][i-1] )
                           + hc[2]* ( vx[j][i+1]-vx[j][i-2] )
                           + hc[3]* ( vx[j][i+2]-vx[j][i-3] )
                           + hc[4]* ( vx[j][i+3]-vx[j][i-4] )
                           + hc[5]* ( vx[j][i+4]-vx[j][i-5] )
                           );
                    *(vyy_1_j+i)= ( hc[1]* ( vy[j][i]  -vy[j-1][i] )
                           + hc[2]* ( vy[j+1][i]-vy[j-2][i] )
                           + hc[3]* ( vy[j+2][i]-vy[j-3][i] )
                           + hc[4]* ( vy[j+3][i]-vy[j-4][i] )
                           + hc[5]* ( vy[j+4][i]-vy[j-5][i] )
                           );
                    *(vyx_1_j+i) = ( hc[1]* ( vy[j][i+1]-vy[j][i] )
                           + hc[2]* ( vy[j][i+2]-vy[j][i-1] )
                           + hc[3]* ( vy[j][i+3]-vy[j][i-2] )
                           + hc[4]* ( vy[j][i+4]-vy[j][i-3] )
                           + hc[5]* ( vy[j][i+5]-vy[j][i-4] )
                           );
                    *(vxy_1_j+i) = ( hc[1]* ( vx[j+1][i]-vx[j][i] )
                           + hc[2]* ( vx[j+2][i]-vx[j-1][i] )
                           + hc[3]* ( vx[j+3][i]-vx[j-2][i] )
                           + hc[4]* ( vx[j+4][i]-vx[j-3][i] )
                           + hc[5]* ( vx[j+5][i]-vx[j-4][i] )
                           );
                    
                    fipjp=uipjp[j][i]*DT;
                    f=u[j][i]*DT;
                    g=pi[j][i]*DT;
                    
                    // Calculate Adams-Bashforth stuff
                    sumxx=c1*(*(vxx_1_j+i))+c2*(*(vxx_2_j+i))+c3*(*(vxx_3_j+i))+c4*(*(vxx_4_j+i));
                    sumyy=c1*(*(vyy_1_j+i))+c2*(*(vyy_2_j+i))+c3*(*(vyy_3_j+i))+c4*(*(vyy_4_j+i));
                    sumxy=c1*(*(vxy_1_j+i))+c2*(*(vxy_2_j+i))+c3*(*(vxy_3_j+i))+c4*(*(vxy_4_j+i));
                    sumyx=c1*(*(vyx_1_j+i))+c2*(*(vyx_2_j+i))+c3*(*(vyx_3_j+i))+c4*(*(vyx_4_j+i));
                    
                    // Update stress
                    sxy[j][i]+= ( fipjp* ( sumxy+sumyx ) )*dhi;
                    sxx[j][i]+= (( g* ( sumxx+sumyy ) )- ( 2.0*f*sumyy ))*dhi;
                    syy[j][i]+= (( g* ( sumxx+sumyy ) )- ( 2.0*f*sumxx ))*dhi;

                    
                }
            }
            break;
            
        case 12:
            for ( j=gy[2]+1; j<=gy[3]; j++ ) {
                vxx_1_j=*(vxx_1+j),vyy_1_j=*(vyy_1+j),vxy_1_j=*(vxy_1+j),vyx_1_j=*(vyx_1+j);
                vxx_2_j=*(vxx_2+j),vyy_2_j=*(vyy_2+j),vxy_2_j=*(vxy_2+j),vyx_2_j=*(vyx_2+j);
                vxx_3_j=*(vxx_3+j),vyy_3_j=*(vyy_3+j),vxy_3_j=*(vxy_3+j),vyx_3_j=*(vyx_3+j);
                vxx_4_j=*(vxx_4+j),vyy_4_j=*(vyy_4+j),vxy_4_j=*(vxy_4+j),vyx_4_j=*(vyx_4+j);
                for ( i=gx[2]+1; i<=gx[3]; i++ ) {
                    *(vxx_1_j+i) = ( hc[1]* ( vx[j][i]  -vx[j][i-1] )
                           + hc[2]* ( vx[j][i+1]-vx[j][i-2] )
                           + hc[3]* ( vx[j][i+2]-vx[j][i-3] )
                           + hc[4]* ( vx[j][i+3]-vx[j][i-4] )
                           + hc[5]* ( vx[j][i+4]-vx[j][i-5] )
                           + hc[6]* ( vx[j][i+5]-vx[j][i-6] )
                           );
                    *(vyy_1_j+i) = ( hc[1]* ( vy[j][i]  -vy[j-1][i] )
                           + hc[2]* ( vy[j+1][i]-vy[j-2][i] )
                           + hc[3]* ( vy[j+2][i]-vy[j-3][i] )
                           + hc[4]* ( vy[j+3][i]-vy[j-4][i] )
                           + hc[5]* ( vy[j+4][i]-vy[j-5][i] )
                           + hc[6]* ( vy[j+5][i]-vy[j-6][i] )
                           );
                    *(vyx_1_j+i)= ( hc[1]* ( vy[j][i+1]-vy[j][i] )
                           + hc[2]* ( vy[j][i+2]-vy[j][i-1] )
                           + hc[3]* ( vy[j][i+3]-vy[j][i-2] )
                           + hc[4]* ( vy[j][i+4]-vy[j][i-3] )
                           + hc[5]* ( vy[j][i+5]-vy[j][i-4] )
                           + hc[6]* ( vy[j][i+6]-vy[j][i-5] )
                           );
                    *(vxy_1_j+i)= ( hc[1]* ( vx[j+1][i]-vx[j][i] )
                           + hc[2]* ( vx[j+2][i]-vx[j-1][i] )
                           + hc[3]* ( vx[j+3][i]-vx[j-2][i] )
                           + hc[4]* ( vx[j+4][i]-vx[j-3][i] )
                           + hc[5]* ( vx[j+5][i]-vx[j-4][i] )
                           + hc[6]* ( vx[j+6][i]-vx[j-5][i] )
                           );
                    
                    fipjp=uipjp[j][i]*DT;
                    f=u[j][i]*DT;
                    g=pi[j][i]*DT;
                    
                    // Calculate Adams-Bashforth stuff
                    sumxx=c1*(*(vxx_1_j+i))+c2*(*(vxx_2_j+i))+c3*(*(vxx_3_j+i))+c4*(*(vxx_4_j+i));
                    sumyy=c1*(*(vyy_1_j+i))+c2*(*(vyy_2_j+i))+c3*(*(vyy_3_j+i))+c4*(*(vyy_4_j+i));
                    sumxy=c1*(*(vxy_1_j+i))+c2*(*(vxy_2_j+i))+c3*(*(vxy_3_j+i))+c4*(*(vxy_4_j+i));
                    sumyx=c1*(*(vyx_1_j+i))+c2*(*(vyx_2_j+i))+c3*(*(vyx_3_j+i))+c4*(*(vyx_4_j+i));
                    
                    // Update stress
                    sxy[j][i]+= ( fipjp* ( sumxy+sumyx ) )*dhi;
                    sxx[j][i]+= (( g* ( sumxx+sumyy ) )- ( 2.0*f*sumyy ))*dhi;
                    syy[j][i]+= (( g* ( sumxx+sumyy ) )- ( 2.0*f*sumxx ))*dhi;

                }
            }
            break;
            
        default: /* CASE 2 */
            for ( j=gy[2]+1; j<=gy[3]; j++ ) {
                vxx_1_j=*(vxx_1+j),vyy_1_j=*(vyy_1+j),vxy_1_j=*(vxy_1+j),vyx_1_j=*(vyx_1+j);
                vxx_2_j=*(vxx_2+j),vyy_2_j=*(vyy_2+j),vxy_2_j=*(vxy_2+j),vyx_2_j=*(vyx_2+j);
                vxx_3_j=*(vxx_3+j),vyy_3_j=*(vyy_3+j),vxy_3_j=*(vxy_3+j),vyx_3_j=*(vyx_3+j);
                vxx_4_j=*(vxx_4+j),vyy_4_j=*(vyy_4+j),vxy_4_j=*(vxy_4+j),vyx_4_j=*(vyx_4+j);
                for ( i=gx[2]+1; i<=gx[3]; i++ ) {
                    
                    
                    /* spatial derivatives of the components of the velocities */
                    /* using Holberg coefficients */
                    *(vxx_1_j+i) = hc[1]* ( vx[j][i]  -vx[j][i-1] );
                    *(vyy_1_j+i) = hc[1]* ( vy[j][i]  -vy[j-1][i] );
                    *(vyx_1_j+i) = hc[1]* ( vy[j][i+1]-vy[j][i] );
                    *(vxy_1_j+i) = hc[1]* ( vx[j+1][i]-vx[j][i] );
                    
                    /* updating components of the stress tensor, partially */
                    fipjp=uipjp[j][i]*DT;
                    f=u[j][i]*DT;
                    g=pi[j][i]*DT;
                    
                    // Calculate Adams-Bashforth stuff
                    sumxx=c1*(*(vxx_1_j+i))+c2*(*(vxx_2_j+i))+c3*(*(vxx_3_j+i))+c4*(*(vxx_4_j+i));
                    sumyy=c1*(*(vyy_1_j+i))+c2*(*(vyy_2_j+i))+c3*(*(vyy_3_j+i))+c4*(*(vyy_4_j+i));
                    sumxy=c1*(*(vxy_1_j+i))+c2*(*(vxy_2_j+i))+c3*(*(vxy_3_j+i))+c4*(*(vxy_4_j+i));
                    sumyx=c1*(*(vyx_1_j+i))+c2*(*(vyx_2_j+i))+c3*(*(vyx_3_j+i))+c4*(*(vyx_4_j+i));
                    
                    // Update stress
                    sxy[j][i]+= ( fipjp* ( sumxy+sumyx ) )*dhi;
                    sxx[j][i]+= (( g* ( sumxx+sumyy ) )- ( 2.0*f*sumyy ))*dhi;
                    syy[j][i]+= (( g* ( sumxx+sumyy ) )- ( 2.0*f*sumxx ))*dhi;

                }
            }
            break;
            
    } /* end of switch(FDORDER) */
    
    
    
    if ( ( MYID==0 ) && ( ( nt+ ( OUTNTIMESTEPINFO-1 ) ) %OUTNTIMESTEPINFO ) ==0 ) {
        time2=MPI_Wtime();
        fprintf ( FP," finished (real time: %4.3f s).\n",time2-time1 );
    }
}


