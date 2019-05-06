/*-----------------------------------------------------------------------------------------
 * Copyright (C) 2016  For the list of authors, see file AUTHORS.
 *
 * This file is part of IFOS.
 *
 * IFOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 *
 * IFOS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with IFOS. See file COPYING and/or <http://www.gnu.org/licenses/gpl-2.0.html>.
 -----------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
 *   stress free surface condition
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void surface_elastic_PML(int ndepth, float ** vx, float ** vy, float ** sxx, float ** syy,
                         float ** sxy, float ** syz, float  **  pi, float  **  u, float ** rho, float * hc, float * K_x, float * a_x, float * b_x, float ** psi_vxx, float ** ux, float ** uy, float ** uxy, float ** uyz,float ** sxz,float **uxz){
    
    
    int i,j,m,h,h1;
    int fdoh;
    float fjm, g;
    float  vxx, vyy;
    float  dh24, dthalbe;
    extern float DT, DH;
    extern int NX, PARAMETERIZATION;
    extern int FW, BOUNDARY;
    extern int NPROCX, NPROCY, POS[3], MYID;
    extern int FDORDER, WAVETYPE;
    
    fdoh = FDORDER/2;
    dthalbe=DT/2.0;
    dh24=1.0/DH;
    
    j=ndepth;     /* The free surface is located exactly in y=1/2*dh !! */
    if (WAVETYPE==1||WAVETYPE==3){
        for (i=1;i<=NX;i++){
            
            /*Mirroring the components of the stress tensor to make
             a stress free surface (method of imaging)*/
            syy[j][i]=0.0;
            uy[j][i]=0.0;
            
            /*syy[j-1][i]=-syy[j+1][i];
             sxy[j-1][i]=-sxy[j][i];*/
            
            vxx = 0.0;
            vyy = 0.0;
            for (m=1; m<=fdoh; m++) {
                
                /*Mirroring the components of the stress tensor to make
                 a stress free surface (method of imaging)*/
                syy[j-m][i]=-syy[j+m][i];
                sxy[j-m][i]=-sxy[j+m-1][i];
                
                uy[j-m][i]=-uy[j+m][i];
                uxy[j-m][i]=-uxy[j+m-1][i];
                
                vxx += hc[m]*(vx[j][i+m-1] -vx[j][i-m]);
                vyy += hc[m]*(vy[j+m-1][i] -vy[j-m][i]);
            }
            vxx *= dh24;
            vyy *= dh24;
            
            /* apply PML boundary */
            /* left boundary */
            if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
                
                psi_vxx[j][i] = b_x[i] * psi_vxx[j][i] + a_x[i] * vxx;
                vxx = vxx / K_x[i] + psi_vxx[j][i];
            }
            
            /* right boundary */
            if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=NX-FW+1)){
                
                h1 = (i-NX+2*FW);
                h = i;
                
                psi_vxx[j][h1] = b_x[h1] * psi_vxx[j][h1] + a_x[h1] * vxx;
                vxx = vxx / K_x[h1] + psi_vxx[j][h1];
            }
            
            
            if(PARAMETERIZATION==3){
                fjm=u[j][i]*2.0;
                g=pi[j][i];}
            
            if(PARAMETERIZATION==1){
                fjm=rho[j][i] * u[j][i] * u[j][i] * 2.0;
                g=rho[j][i] * ((pi[j][i] * pi[j][i]) - 2 * u[j][i] * u[j][i]);}
            
            
            /*sxx[j][i]+= DT*((4.0*((g*fjm)+(fjm*fjm))/(g+2*fjm))*vxx);*/
            sxx[j][i]+= -DT*((g*g)/(g+fjm)*vxx+g*vyy);
            ux[j][i]=((2.0*g*fjm+fjm*fjm)/(g+fjm))*vxx;
            
            
        }
    }
    if (WAVETYPE==2||WAVETYPE==3){
        for (i=1;i<=NX;i++){
        
            syz[j][i]=0.0;
            uyz[j][i]=0.0;
            for (m=1; m<=fdoh; m++) {
                syz[j-m][i]=-syz[j+m][i];
                sxz[j-m][i]=-sxz[j+m-1][i];
                
                uyz[j-m][i]=-uyz[j+m][i];
                uxz[j-m][i]=-uxz[j+m-1][i];
            }
        }
    }
}
