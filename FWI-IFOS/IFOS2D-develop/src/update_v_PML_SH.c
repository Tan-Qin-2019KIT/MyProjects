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
 *   updating particle velocities at gridpoints [nx1...nx2][ny1...ny2]
 *   by a staggered grid finite difference scheme of FDORDER accuracy in space
 *   and second order accuracy in time
 *  ----------------------------------------------------------------------*/

#include "fd.h"



void update_v_PML_SH(int nx1, int nx2, int ny1, int ny2, int nt,
                     float **  vz, float **  vzp1, float **  vzm1, float ** sxz, float ** syz,float **rjp, float **  srcpos_loc, float ** signals, float ** signals1, int nsrc, float ** absorb_coeff,
                     float *hc, int infoout,int sw, float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
                     float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,float ** psi_sxz_x,float ** psi_syz_y){
    
    int i, j,l,fdoh, h, h1;
    float dtdh;
    float sxz_x, syz_y;
    extern float DT, DH;
    double time1 = 0.0, time2;
    extern int MYID, SOURCE_TYPE, ADJOINT_TYPE, FDORDER;
    extern int FDORDER, PARAMETERIZATION;
    extern int FREE_SURF, BOUNDARY, FW;
    extern int NPROCX, NPROCY, POS[3];
    extern FILE *FP;
    extern int VELOCITY;
    
    fdoh = FDORDER/2;
    dtdh = DT*DT/DH;
    
    if (infoout && (MYID==0)){
        time1=MPI_Wtime();
        fprintf(FP,"\n **Message from update_v_SH (printed by PE %d):\n",MYID);
        fprintf(FP," Updating particle velocities ...");
    }
    
    /* ------------------------------------------------------------
     * Important!
     * rip and rjp are reciprocal values of averaged densities
     * ------------------------------------------------------------ */
    
    switch (FDORDER){
        case 2:
            for (j=ny1;j<=ny2;j++){
                for (i=nx1;i<=nx2;i++){
                    
                    sxz_x=  hc[1]*(sxz[j][i]-sxz[j][i-1]);
                    
                    syz_y=  hc[1]*(syz[j+1][i]-syz[j][i]);
                    
                    /* left boundary */
                    if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
                        psi_sxz_x[j][i] = b_x[i] * psi_sxz_x[j][i] + a_x[i] * sxz_x;
                        sxz_x = sxz_x / K_x[i] + psi_sxz_x[j][i];
                        
                    }
                    
                    /* right boundary */
                    if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){
                        h1 = (i-nx2+2*FW);
                        h=i;
                        psi_sxz_x[j][h1] = b_x[h1] * psi_sxz_x[j][h1] + a_x[h1] * sxz_x;
                        sxz_x = sxz_x / K_x[h1] + psi_sxz_x[j][h1];
                    }
                    
                    /* top boundary */
                    if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
                        psi_syz_y[j][i] = b_y_half[j] * psi_syz_y[j][i] + a_y_half[j] * syz_y;
                        syz_y = syz_y / K_y_half[j] + psi_syz_y[j][i];
                    }
                    
                    /* bottom boundary */
                    if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){
                        h1 = (j-ny2+2*FW);
                        h = j;
                        psi_syz_y[h1][i] = b_y_half[h1] * psi_syz_y[h1][i] + a_y_half[h1] * syz_y;
                        syz_y = syz_y / K_y_half[h1] + psi_syz_y[h1][i];
                    }
                    
                    vzp1[j][i] = vz[j][i];
                    
                    vz[j][i] += rjp[j][i]*DT*(sxz_x+syz_y)/(DH);
                    
                }
            }
            break;
            
        case 4:
            for (j=ny1;j<=ny2;j++){
                for (i=nx1;i<=nx2;i++){
                    
                    sxz_x=  hc[1]*(sxz[j][i]-sxz[j][i-1])+hc[2]*(sxz[j][i+1]-sxz[j][i-2]);
                    syz_y=  hc[1]*(syz[j+1][i]-syz[j][i])+hc[2]*(syz[j+2][i]-syz[j-1][i]);
                    
                    /* left boundary */
                    if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
                        psi_sxz_x[j][i] = b_x[i] * psi_sxz_x[j][i] + a_x[i] * sxz_x;
                        sxz_x = sxz_x / K_x[i] + psi_sxz_x[j][i];
                        
                    }
                    
                    /* right boundary */
                    if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){
                        h1 = (i-nx2+2*FW);
                        h=i;
                        psi_sxz_x[j][h1] = b_x[h1] * psi_sxz_x[j][h1] + a_x[h1] * sxz_x;
                        sxz_x = sxz_x / K_x[h1] + psi_sxz_x[j][h1];
                    }
                    
                    /* top boundary */
                    if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
                        psi_syz_y[j][i] = b_y_half[j] * psi_syz_y[j][i] + a_y_half[j] * syz_y;
                        syz_y = syz_y / K_y_half[j] + psi_syz_y[j][i];
                    }
                    
                    /* bottom boundary */
                    if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){
                        h1 = (j-ny2+2*FW);
                        h = j;
                        psi_syz_y[h1][i] = b_y_half[h1] * psi_syz_y[h1][i] + a_y_half[h1] * syz_y;
                        syz_y = syz_y / K_y_half[h1] + psi_syz_y[h1][i];
                    }
                    
                    vzp1[j][i] = vz[j][i];
                    
                    vz[j][i] += rjp[j][i]*DT*(sxz_x+syz_y)/(DH);
                    
                }
            }
            
            break;
            
        case 6:
            for (j=ny1;j<=ny2;j++){
                for (i=nx1;i<=nx2;i++){
                    sxz_x=  hc[1]*(sxz[j][i]-sxz[j][i-1])+hc[2]*(sxz[j][i+1]-sxz[j][i-2])+hc[3]*(sxz[j][i+2]-sxz[j][i-3]);
                    syz_y=  hc[1]*(syz[j+1][i]-syz[j][i])+hc[2]*(syz[j+2][i]-syz[j-1][i])+hc[3]*(syz[j+3][i]-syz[j-2][i]);
                    
                    /* left boundary */
                    if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
                        psi_sxz_x[j][i] = b_x[i] * psi_sxz_x[j][i] + a_x[i] * sxz_x;
                        sxz_x = sxz_x / K_x[i] + psi_sxz_x[j][i];
                        
                    }
                    
                    /* right boundary */
                    if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){
                        h1 = (i-nx2+2*FW);
                        h=i;
                        psi_sxz_x[j][h1] = b_x[h1] * psi_sxz_x[j][h1] + a_x[h1] * sxz_x;
                        sxz_x = sxz_x / K_x[h1] + psi_sxz_x[j][h1];
                    }
                    
                    /* top boundary */
                    if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
                        psi_syz_y[j][i] = b_y_half[j] * psi_syz_y[j][i] + a_y_half[j] * syz_y;
                        syz_y = syz_y / K_y_half[j] + psi_syz_y[j][i];
                    }
                    
                    /* bottom boundary */
                    if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){
                        h1 = (j-ny2+2*FW);
                        h = j;
                        psi_syz_y[h1][i] = b_y_half[h1] * psi_syz_y[h1][i] + a_y_half[h1] * syz_y;
                        syz_y = syz_y / K_y_half[h1] + psi_syz_y[h1][i];
                    }
                    
                    vzp1[j][i] = vz[j][i];
                    
                    vz[j][i] += rjp[j][i]*DT*(sxz_x+syz_y)/(DH);
                    
                }
            }
            break;
            
        case 8:
            for (j=ny1;j<=ny2;j++){
                for (i=nx1;i<=nx2;i++){
                    sxz_x=  hc[1]*(sxz[j][i]-sxz[j][i-1])+hc[2]*(sxz[j][i+1]-sxz[j][i-2])+hc[3]*(sxz[j][i+2]-sxz[j][i-3])+hc[4]*(sxz[j][i+3]-sxz[j][i-4]);
                    syz_y=  hc[1]*(syz[j+1][i]-syz[j][i])+hc[2]*(syz[j+2][i]-syz[j-1][i])+hc[3]*(syz[j+3][i]-syz[j-2][i])+hc[4]*(syz[j+4][i]-syz[j-3][i]);
                    
                    /* left boundary */
                    if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
                        psi_sxz_x[j][i] = b_x[i] * psi_sxz_x[j][i] + a_x[i] * sxz_x;
                        sxz_x = sxz_x / K_x[i] + psi_sxz_x[j][i];
                        
                    }
                    
                    /* right boundary */
                    if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){
                        h1 = (i-nx2+2*FW);
                        h=i;
                        psi_sxz_x[j][h1] = b_x[h1] * psi_sxz_x[j][h1] + a_x[h1] * sxz_x;
                        sxz_x = sxz_x / K_x[h1] + psi_sxz_x[j][h1];
                    }
                    
                    /* top boundary */
                    if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
                        psi_syz_y[j][i] = b_y_half[j] * psi_syz_y[j][i] + a_y_half[j] * syz_y;
                        syz_y = syz_y / K_y_half[j] + psi_syz_y[j][i];
                    }
                    
                    /* bottom boundary */
                    if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){
                        h1 = (j-ny2+2*FW);
                        h = j;
                        psi_syz_y[h1][i] = b_y_half[h1] * psi_syz_y[h1][i] + a_y_half[h1] * syz_y;
                        syz_y = syz_y / K_y_half[h1] + psi_syz_y[h1][i];
                    }
                    
                    vzp1[j][i] = vz[j][i];
                    
                    vz[j][i] += rjp[j][i]*DT*(sxz_x+syz_y)/(DH);
                    
                }
            }
            break;
            
        case 10:
            for (j=ny1;j<=ny2;j++){
                for (i=nx1;i<=nx2;i++){
                    sxz_x=  hc[1]*(sxz[j][i]-sxz[j][i-1])+hc[2]*(sxz[j][i+1]-sxz[j][i-2])+hc[3]*(sxz[j][i+2]-sxz[j][i-3])+hc[4]*(sxz[j][i+3]-sxz[j][i-4])+hc[5]*(sxz[j][i+4]-sxz[j][i-5]);
                    syz_y=  hc[1]*(syz[j+1][i]-syz[j][i])+hc[2]*(syz[j+2][i]-syz[j-1][i])+hc[3]*(syz[j+3][i]-syz[j-2][i])+hc[4]*(syz[j+4][i]-syz[j-3][i])+hc[5]*(syz[j+5][i]-syz[j-4][i]);
                    
                    /* left boundary */
                    if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
                        psi_sxz_x[j][i] = b_x[i] * psi_sxz_x[j][i] + a_x[i] * sxz_x;
                        sxz_x = sxz_x / K_x[i] + psi_sxz_x[j][i];
                        
                    }
                    
                    /* right boundary */
                    if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){
                        h1 = (i-nx2+2*FW);
                        h=i;
                        psi_sxz_x[j][h1] = b_x[h1] * psi_sxz_x[j][h1] + a_x[h1] * sxz_x;
                        sxz_x = sxz_x / K_x[h1] + psi_sxz_x[j][h1];
                    }
                    
                    /* top boundary */
                    if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
                        psi_syz_y[j][i] = b_y_half[j] * psi_syz_y[j][i] + a_y_half[j] * syz_y;
                        syz_y = syz_y / K_y_half[j] + psi_syz_y[j][i];
                    }
                    
                    /* bottom boundary */
                    if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){
                        h1 = (j-ny2+2*FW);
                        h = j;
                        psi_syz_y[h1][i] = b_y_half[h1] * psi_syz_y[h1][i] + a_y_half[h1] * syz_y;
                        syz_y = syz_y / K_y_half[h1] + psi_syz_y[h1][i];
                    }
                    
                    vzp1[j][i] = vz[j][i];
                    
                    vz[j][i] += rjp[j][i]*DT*(sxz_x+syz_y)/(DH);
                    
                }
            }
            break;
            
        case 12:
            for (j=ny1;j<=ny2;j++){
                for (i=nx1;i<=nx2;i++){
                    sxz_x=  hc[1]*(sxz[j][i]-sxz[j][i-1])+hc[2]*(sxz[j][i+1]-sxz[j][i-2])+hc[3]*(sxz[j][i+2]-sxz[j][i-3])+hc[4]*(sxz[j][i+3]-sxz[j][i-4])+hc[5]*(sxz[j][i+4]-sxz[j][i-5])+hc[6]*(sxz[j][i+5]-sxz[j][i-6]);
                    syz_y=  hc[1]*(syz[j+1][i]-syz[j][i])+hc[2]*(syz[j+2][i]-syz[j-1][i])+hc[3]*(syz[j+3][i]-syz[j-2][i])+hc[4]*(syz[j+4][i]-syz[j-3][i])+hc[5]*(syz[j+5][i]-syz[j-4][i])+hc[6]*(syz[j+6][i]-syz[j-5][i]);
                    
                    /* left boundary */
                    if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
                        psi_sxz_x[j][i] = b_x[i] * psi_sxz_x[j][i] + a_x[i] * sxz_x;
                        sxz_x = sxz_x / K_x[i] + psi_sxz_x[j][i];
                        
                    }
                    
                    /* right boundary */
                    if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){
                        h1 = (i-nx2+2*FW);
                        h=i;
                        psi_sxz_x[j][h1] = b_x[h1] * psi_sxz_x[j][h1] + a_x[h1] * sxz_x;
                        sxz_x = sxz_x / K_x[h1] + psi_sxz_x[j][h1];
                    }
                    
                    /* top boundary */
                    if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
                        psi_syz_y[j][i] = b_y_half[j] * psi_syz_y[j][i] + a_y_half[j] * syz_y;
                        syz_y = syz_y / K_y_half[j] + psi_syz_y[j][i];
                    }
                    
                    /* bottom boundary */
                    if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){
                        h1 = (j-ny2+2*FW);
                        h = j;
                        psi_syz_y[h1][i] = b_y_half[h1] * psi_syz_y[h1][i] + a_y_half[h1] * syz_y;
                        syz_y = syz_y / K_y_half[h1] + psi_syz_y[h1][i];
                    }
                    
                    vzp1[j][i] = vz[j][i];
                    
                    vz[j][i] += rjp[j][i]*DT*(sxz_x+syz_y)/(DH);
                    
                }
            }
            break;
    } /* end of switch(FDORDER) */
    
    
    
    /* Forward Modelling (sw==0) */
    if(sw==0){
        for (l=1;l<=nsrc;l++) {
            i=(int)srcpos_loc[1][l];
            j=(int)srcpos_loc[2][l];
            vz[j][i] += (rjp[j][i]*DT*1.0*signals1[l][nt])/(DH*DH);  /* single force in z */
            
        }
    }
    
    /* Backpropagation (sw==1) */
    if(sw==1){
        for (l=1;l<=nsrc;l++) {
            i=(int)srcpos_loc[1][l];
            j=(int)srcpos_loc[2][l];
            
            vz[j][i] += signals[l][nt];  /* + single force in z */
            
        }
    }
    
    
    if (infoout && (MYID==0)){
        time2=MPI_Wtime();
        fprintf(FP," finished (real time: %4.2f s).\n",time2-time1);
    }
}