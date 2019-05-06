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

/*
 * Smoothing gradient / model with a median filter
 */

#include "fd.h"

void smooth(float ** mat, int sws, int filter, float Vs_avg, float F_LOW_PASS)
{
    
    /* extern variables */
    
    extern float DH, A;
    extern int FREE_SURF, NX, NY, NXG, NYG;
    extern int NPROCX, NPROCY, MYID, POS[3],WAVETYPE,VERBOSE;
    extern char JACOBIAN[STRING_SIZE], INV_MODELFILE[STRING_SIZE];
    extern FILE *FP;
    extern int FILT_SIZE_GRAD, FILT_SIZE, MODEL_FILTER, GRAD_FILTER, TIME_FILT,GRAD_FILT_WAVELENGTH;
    extern int VERBOSE;
    /* local variables */
    int i, j, ii, jj;
    int i1, j1, filtsize, hfs;
    
    float **model_tmp, **model_med, **filterpart, grad, normgauss, smooth_meter;
    
    char jac_tmp[STRING_SIZE];
    
    FILE *model;
    
    char modfile[STRING_SIZE];
    
    float ** global_matrix;
    
    global_matrix=get_global_from_local_matrix(mat);
    
    switch (filter){
        case 1:
            if((GRAD_FILT_WAVELENGTH==1)&&(TIME_FILT==1)){
                if(VERBOSE) printf("\n -------------------------------------------------------------------------- \n");
                if(VERBOSE) printf("\n Calculating a wavelength dependent filter size for smoothing the gradient: \n");
                FILT_SIZE_GRAD = (int)(Vs_avg/F_LOW_PASS*A/DH);
                if(VERBOSE) printf("\n FILT_SIZE_GRAD = Vs_avg = %4.2f m/s / F_LOW_PASS = %4.2f Hz * weighting factor A = %4.2f / grid spacing DH = %4.2f m  \n",Vs_avg,F_LOW_PASS,A,DH);
                if(VERBOSE) printf("\n New FILT_SIZE_GRAD = %d (grid points) is used (-> %4.2f m).                \n",FILT_SIZE_GRAD,FILT_SIZE_GRAD*DH);
            }
            if (FILT_SIZE_GRAD==0)  return;
            if (!(FILT_SIZE_GRAD % 2)) {
                if (FILT_SIZE_GRAD > 0) FILT_SIZE_GRAD += 1;
                else            FILT_SIZE_GRAD -= 1;
            }
            hfs = abs(FILT_SIZE_GRAD)/2;
            if(VERBOSE) printf("\n ----------------------------------------------------------------\n");
            if(VERBOSE) printf("\n Filter size is %d gridpoints, half filter size is %d gridpoints.\n",FILT_SIZE_GRAD,hfs);
            filterpart=matrix(1,abs(FILT_SIZE_GRAD),1,abs(FILT_SIZE_GRAD));
            model_tmp = matrix(-hfs+1,NYG+hfs,-hfs+1,NXG+hfs);
            break;
            
        case 2:
            if (FILT_SIZE==0)   return;
            if (!(FILT_SIZE % 2)) {
                if (FILT_SIZE > 0)  FILT_SIZE += 1;
                else            FILT_SIZE -= 1;
            }
            hfs = abs(FILT_SIZE)/2;
            if(VERBOSE) printf("\n ----------------------------------------------------------------\n");
            if(VERBOSE) printf("\n Filter size is %d gridpoints, half filter size is %d gridpoints.\n",FILT_SIZE,hfs);
            filterpart=matrix(1,abs(FILT_SIZE),1,abs(FILT_SIZE));
            model_tmp = matrix(-hfs+1,NYG+hfs,-hfs+1,NXG+hfs);
            break;
    }
    
    /* load merged model */
    for (i=1;i<=NXG;i++){
        for (j=1;j<=NYG;j++){
            model_tmp[j][i]=global_matrix[j][i];
        }
    }
    
    
    /* apply 2D-Gaussian filter on vp and vs model */
    /* extrapolate array */
    /* left/right boundary */
    for (j=1;j<=NYG;j++){
        for (i=-hfs+1;i<=0;i++){
            model_tmp[j][i] = model_tmp[j][1];}
        for (i=NXG+1;i<=NXG+hfs;i++){
            model_tmp[j][i] = model_tmp[j][NXG];}
    }
    /* top/bottom boundary incl. corners */
    for (j=-hfs+1;j<=0;j++){
        for (i=-hfs+1;i<=NXG+hfs;i++){
            model_tmp[j][i] = model_tmp[1][i];}
    }
    for (j=NYG+1;j<=NYG+hfs;j++){
        for (i=-hfs+1;i<=NXG+hfs;i++){
            model_tmp[j][i] = model_tmp[NYG][i];}
    }
    
    /* filter */
    for (j=1;j<=NY;j++){
        for (i=1;i<=NX;i++){
            /* create a filtersize x filtersize matrix */
            for (ii=-hfs;ii<=hfs;ii++){
                for (jj=-hfs;jj<=hfs;jj++){
                    
                    filterpart[ii+hfs+1][jj+hfs+1] = model_tmp[j+POS[2]*NY+jj][i+POS[1]*NX+ii];
                }
            }
            /* filter */
            switch (filter){
                case 1:
                    mat[j][i] = median2d(filterpart,abs(FILT_SIZE_GRAD),abs(FILT_SIZE_GRAD));
                    break;
                case 2:
                    mat[j][i] = median2d(filterpart,abs(FILT_SIZE),abs(FILT_SIZE));
                    break;
            }            
        }
    }
        
    switch (filter){
        case 1:
            free_matrix(model_tmp,-hfs+1,NYG+hfs,-hfs+1,NXG+hfs);
            free_matrix(filterpart,1,abs(FILT_SIZE_GRAD),1,abs(FILT_SIZE_GRAD));
            break;
        case 2:
            free_matrix(model_tmp,-hfs+1,NYG+hfs,-hfs+1,NXG+hfs);
            free_matrix(filterpart,1,abs(FILT_SIZE),1,abs(FILT_SIZE));
            break;
    }    
    
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    switch (filter){
        case 1:
            smooth_meter=FILT_SIZE_GRAD*DH;
            if(VERBOSE) fprintf(FP,"\n Gradient %s is smoothed with filter length of %4.2f meter.\n",jac_tmp,smooth_meter);
            break;
        case 2:
            smooth_meter=FILT_SIZE*DH;
            if(VERBOSE) fprintf(FP,"\n Model %s is smoothed with filter length of %4.2f meter.\n",jac_tmp,smooth_meter);
            break ;
    }
    
}/* end of smoothing */
