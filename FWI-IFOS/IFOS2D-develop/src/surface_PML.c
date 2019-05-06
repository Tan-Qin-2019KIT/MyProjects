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
 *   stress free surface condition *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void surface_PML(int ndepth, float ** vx, float ** vy, float ** sxx, float ** syy,
                 float ** sxy, float ** syz, float ***p, float ***q, float  **  ppi, float  **  pu, float **prho, float **ptaup, float **ptaus, float *etajm, float *peta, float * hc, float * K_x, float * a_x, float * b_x, float ** psi_vxx,
                 float ** ux, float ** uy, float ** uxy, float ** uyz,float ** sxz,float **uxz){
    
    
    int i,j,m,h,h1,l;
    int fdoh;
    float bjm, djm, e, fjm, g;
    float  vxx, vyy, sump=0.0;
    float  dh24, dthalbe;
    float *pts = NULL, ws, sumu = 0.0, sumpi = 0.0, mu = 0.0, pi = 0.0;
    extern float DT, DH, *FL;
    extern int NX, PARAMETERIZATION, L;
    extern int FW, BOUNDARY;
    extern int NPROCX, NPROCY, POS[3], MYID;
    extern int FDORDER,WAVETYPE;
    extern float F_REF;
    
    fdoh = FDORDER/2;
    dthalbe=DT/2.0;
    dh24=1.0/DH;
    
    if (WAVETYPE==1||WAVETYPE==3){
        /* vector for maxwellbodies */
        pts=vector(1,L);
        for (l=1;l<=L;l++) {
            pts[l]=1.0/(2.0*PI*FL[l]);
        }
        
        
        /*ws=2.0*PI*FL[1];*/
        ws=2.0*PI*F_REF;
        
        sumu=0.0;
        sumpi=0.0;
        for (l=1;l<=L;l++){
            sumu=sumu+((ws*ws*pts[l]*pts[l])/(1.0+ws*ws*pts[l]*pts[l]));
            sumpi=sumpi+((ws*ws*pts[l]*pts[l])/(1.0+ws*ws*pts[l]*pts[l]));
        }
    }
    
    
    j=ndepth;     /* The free surface is located exactly in y=1/2*dh !! */
    if (WAVETYPE==1||WAVETYPE==3){
        for (i=1;i<=NX;i++){
            for (l=1;l<=L;l++){
                etajm[l]=peta[l];
            }
            
            /*Mirroring the components of the stress tensor to make
             a stress free surface (method of imaging)*/
            syy[j][i]=0.0;
            uy[j][i]=0.0;
            
            /* since syy is zero on the free surface also the
             corresponding memory-variables must set to zero */
            for (l=1;l<=L;l++) q[j][i][l]=0.0;
            
            /* now updating the stress component sxx and the memory-
             variables p[j][i][l] at the free surface */
            
            /* first calculate spatial derivatives of components
             of particle velocities */
            
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
            
            if (PARAMETERIZATION==1){
                mu=(pu[j][i]*pu[j][i]*prho[j][i])/(1.0+sumu*ptaus[j][i]);
                pi=(ppi[j][i]*ppi[j][i]*prho[j][i])/(1.0+sumpi*ptaup[j][i]);
            }
            if (PARAMETERIZATION==3){
                mu=pu[j][i]/(1.0+sumu*ptaus[j][i]);
                pi=(ppi[j][i]+2*pu[j][i])/(1.0+sumpi*ptaup[j][i]);
            }
            
            /* sums used in updating sxx */
            sump=0.0;
            for (l=1;l<=L;l++) sump+=p[j][i][l];
            
            
            fjm=mu*2.0*(1.0+L*ptaus[j][i]);
            g=pi*(1.0+L*ptaup[j][i]);
            
            /* partially updating sxx */
            sxx[j][i]+= -(DT*(g-fjm)*(g-fjm)*vxx/g)-(DT*(g-fjm)*vyy)-(dthalbe*sump);
            ux[j][i]+= -((g-fjm)*(g-fjm)*vxx/g)-((g-fjm)*vyy)-(0.5*sump);
            
            /* updating the memory-variable p[j][i][l] at the free surface */
            sump=0.0;
            for (l=1;l<=L;l++){
                bjm=etajm[l]/(1.0+(etajm[l]*0.5));
                djm=2.0*mu*ptaus[j][i];
                e=pi*ptaup[j][i];
                p[j][i][l]+=bjm*(((djm-e)*((fjm/g)-1.0)*vxx)-((djm-e)*vyy));
                sump+=p[j][i][l];
            }
            /*completely updating the stress sxx */
            sxx[j][i]+=(dthalbe*sump);
            ux[j][i]+=(0.5*sump);
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
    if (WAVETYPE==1||WAVETYPE==3){
        free_vector(pts,1,L);
    }
}
