/*------------------------------------------------------------------------
 * Copyright (C) 2016 For the list of authors, see file AUTHORS.
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
 * along with 3D-AWAIT. See file COPYING and/or
 * <http://www.gnu.org/licenses/gpl-2.0.html>.
 --------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------
 * Calculatipon of L-BFGS update
 --------------------------------------------------------------------------*/

#include "fd.h"
void lbfgs_reset(int iter, int N_LBFGS, int NPAR_LBFGS,float ** s_LBFGS1, float ** y_LBFGS1, float * rho_LBFGS1);
void lbfgs_core(int iteration, int N_LBFGS, int NPAR_LBFGS,float ** s_LBFGS, float ** y_LBFGS, float * rho_LBFGS,float *q_LBFGS,float *alpha_LBFGS,float *r_LBFGS);

void lbfgs(float **grad_vs, float **grad_rho, float **grad_vp,float Vs_avg,float rho_avg,float Vp_avg, float *rho_LBFGS, float **s_LBFGS, float **y_LBFGS,int N_LBFGS,int NPAR_LBFGS, int iteration, int *LBFGS_iter_start){
    
    /* global */
    extern int NX,NY,MYID;
    extern FILE *FP;
    extern int POS[3];
    extern char JACOBIAN[STRING_SIZE];
    extern int WAVETYPE;
    extern int ACOUSTIC;
    extern float LBFGS_SCALE_GRADIENTS;
    
    /* local */
    int w=0;
    int i,j,l;
    float *q_LBFGS,*alpha_LBFGS,*r_LBFGS;
    char jac[225];
    FILE *FP_JAC = NULL;
    
    
    /*---------------------*/
    /*      DEBUGGING      */
    /*---------------------*/
    
    if(!ACOUSTIC) {
        sprintf(jac,"%s_grad1_vs_it%d",JACOBIAN,iteration);
        write_matrix_disk(grad_vs, jac);
    }
    
    sprintf(jac,"%s_grad1_rho_it%d",JACOBIAN,iteration);
    write_matrix_disk(grad_rho, jac);
    
    if(WAVETYPE==1||WAVETYPE==3) {
        sprintf(jac,"%s_grad1_vp_it%d",JACOBIAN,iteration);
        write_matrix_disk(grad_vp, jac);
    }
    
    /*---------------------*/
    /*      Experimental   */
    /*---------------------*/
    /* Scale the gradients with a constant factor */
    /* This is to avoid numerical instabilities if the gradients are to small (absolute value) */
    /* Do not use this feature, unless you have to */
    if(LBFGS_SCALE_GRADIENTS!=1) {
        if(MYID==0) printf("\n\n Scaling the gradients to ensure L-BFGS stability");
        if(MYID==0) printf("\n Scaling with %f",LBFGS_SCALE_GRADIENTS);
        if(MYID==0) printf("\n This is an experimental feature.");

        for (i=1;i<=NX;i++){
            for (j=1;j<=NY;j++){
                if(!ACOUSTIC) grad_vs[j][i]=grad_vs[j][i]*LBFGS_SCALE_GRADIENTS;
                if(NPAR_LBFGS>1) grad_rho[j][i]=grad_rho[j][i]*LBFGS_SCALE_GRADIENTS;
                if(NPAR_LBFGS>2) grad_vp[j][i]=grad_vp[j][i]*LBFGS_SCALE_GRADIENTS;
            }
        }
    }
    
    /*-------------------------------------------------*/
    /*      Init L-BFGS at iter==LBFGS_iter_start      */
    /*-------------------------------------------------*/
    if(iteration==*LBFGS_iter_start) {
        w=iteration%N_LBFGS;
        if(w==0) w=N_LBFGS;
        
        l=0;
        for (i=1;i<=NX;i++){
            for (j=1;j<=NY;j++){
                l++;
                if(!ACOUSTIC) y_LBFGS[w][l]=-grad_vs[j][i]*Vs_avg; /* VS */
                if(NPAR_LBFGS>1) y_LBFGS[w][l+NX*NY]=-grad_rho[j][i]*rho_avg; /* RHO */
                if(NPAR_LBFGS>2) y_LBFGS[w][l+2*NX*NY]=-grad_vp[j][i]*Vp_avg; /* VP */
            }
        }
    }
    
    /*------------------------*/
    /*      Start L-BFGS      */
    /*------------------------*/
    if(iteration>*LBFGS_iter_start) {
        
        alpha_LBFGS = vector(1,N_LBFGS);
        q_LBFGS = vector(1,NPAR_LBFGS*NX*NY);
        r_LBFGS = vector(1,NPAR_LBFGS*NX*NY);
        
        w=(iteration-1)%N_LBFGS;
        if(w==0) w=N_LBFGS;
        
        /* Debugging */
        if(!ACOUSTIC) {
            sprintf(jac,"%s_y_LBFGS_vs_it%d.bin.%i.%i",JACOBIAN,iteration,POS[1],POS[2]);
            FP_JAC=fopen(jac,"wb");
        }
        
        if(MYID==0) printf("\n\n ------------ L-BFGS ---------------");
        if(MYID==0) printf("\n Start calculation L-BFGS update");
        if(MYID==0) printf("\n At Iteration %i in L-BFGS vector %i\n",iteration,w);
        
        l=0;
        for (i=1;i<=NX;i++){
            for (j=1;j<=NY;j++){
                
                l++;
                
                /* VS */
                if(!ACOUSTIC){
                    y_LBFGS[w][l]+=grad_vs[j][i]*Vs_avg; /* add grad(i) to build grad(i)-grad(i-1) */
                    q_LBFGS[l]=grad_vs[j][i]*Vs_avg; /* Normalisation */
                }
                
                /* RHO */
                if(NPAR_LBFGS>1) {
                    y_LBFGS[w][l+NY*NX]+=grad_rho[j][i]*rho_avg; /* add grad(i) to build grad(i)-grad(i-1) */
                    q_LBFGS[l+NY*NX]=grad_rho[j][i]*rho_avg; /* Normalisation */
                }
                
                /* VP */
                if(NPAR_LBFGS>2) {
                    y_LBFGS[w][l+2*NY*NX]+=grad_vp[j][i]*Vp_avg; /* add grad(i) to build grad(i)-grad(i-1) */
                    q_LBFGS[l+2*NY*NX]=grad_vp[j][i]*Vp_avg; /* Normalisation */
                }
                
                /* Debugging */
                if(!ACOUSTIC) fwrite(&y_LBFGS[w][l],sizeof(float),1,FP_JAC);
            }
        }
        
        /*---------------------*/
        /*      DEBUGGING      */
        /*---------------------*/
        if(!ACOUSTIC) {
            fclose(FP_JAC);
            MPI_Barrier(MPI_COMM_WORLD);
            sprintf(jac,"%s_y_LBFGS_vs_it%d.bin",JACOBIAN,iteration);
            if (MYID==0) mergemod(jac,3);
            MPI_Barrier(MPI_COMM_WORLD);
            sprintf(jac,"%s_y_LBFGS_vs_it%d.bin.%i.%i",JACOBIAN,iteration,POS[1],POS[2]);
            remove(jac);
        }
        
        /*----------------------------------*/
        /*      call L-BFGS Algorithm       */
        /*----------------------------------*/
        lbfgs_core(iteration,N_LBFGS, NX*NY*NPAR_LBFGS,s_LBFGS,y_LBFGS,rho_LBFGS,q_LBFGS,alpha_LBFGS,r_LBFGS);
        
        /*-------------------------------------------------------------*/
        /* Save model pertubation and save gradient for next iteration */
        /*-------------------------------------------------------------*/
        w=iteration%N_LBFGS;
        if(w==0) w=N_LBFGS;
        
        l=0;
        for (i=1;i<=NX;i++){
            for (j=1;j<=NY;j++){
                
                l++;
                
                /* VS */
                if(!ACOUSTIC) {
                    y_LBFGS[w][l]=-grad_vs[j][i]*Vs_avg; /* add -grad(i-1) to build grad(i)-grad(i-1) */
                    grad_vs[j][i]=r_LBFGS[l]*Vs_avg; /* Denormalization */
                }
                
                /* RHO */
                if(NPAR_LBFGS>1) {
                    y_LBFGS[w][l+NY*NX]=-grad_rho[j][i]*rho_avg; /* add -grad(i-1) to build grad(i)-grad(i-1) */
                    grad_rho[j][i]=r_LBFGS[l+NY*NX]*rho_avg; /* Denormalization */
                }
                
                /* VP */
                if(NPAR_LBFGS>2) {
                    y_LBFGS[w][l+2*NY*NX]=-grad_vp[j][i]*Vp_avg; /* add -grad(i-1) to build grad(i)-grad(i-1) */
                    grad_vp[j][i]=r_LBFGS[l+2*NY*NX]*Vp_avg; /* Denormalization */
                }
            }
        }
        
        free_vector(r_LBFGS,1,NPAR_LBFGS*NX*NY);
        free_vector(q_LBFGS,1,NPAR_LBFGS*NX*NY);
        free_vector(alpha_LBFGS, 1, N_LBFGS);
        
    }
    
    /*---------------------*/
    /*      DEBUGGING      */
    /*---------------------*/
    if(!ACOUSTIC){
        sprintf(jac,"%s_grad2_vs_it%d",JACOBIAN,iteration);
        write_matrix_disk(grad_vs, jac);
    }
    
    sprintf(jac,"%s_grad2_rho_it%d",JACOBIAN,iteration);
    write_matrix_disk(grad_rho, jac);
    
    if(WAVETYPE==1||WAVETYPE==3) {
        sprintf(jac,"%s_grad2_vp_it%d",JACOBIAN,iteration);
        write_matrix_disk(grad_vp, jac);
    }
}

void lbfgs_core(int iteration, int N_LBFGS, int NPAR_LBFGS,float ** s_LBFGS, float ** y_LBFGS, float * rho_LBFGS,float *q_LBFGS,float *alpha_LBFGS,float *r_LBFGS) {
    
    extern FILE * FP;
    extern int VERBOSE;
    
    float beta_LBFGS=0.0;
    float dum1=0.0, dum2=0.0, buf1=0.0, buf2=0.0;
    float h0;
    int m=0,v=0,w=0,l=0;
    m=iteration-N_LBFGS;
    if(m<1) m=1;
    
    /*----------------------------------*/
    /* calculate H0 and rho_LBFGS       */
    /*----------------------------------*/
    w=(iteration-1)%N_LBFGS; if(w==0) w=N_LBFGS;
    dum1=0.0; dum2=0.0;
    
    for (l=1;l<=NPAR_LBFGS;l++){
        dum1+=y_LBFGS[w][l]*s_LBFGS[w][l];
        dum2+=y_LBFGS[w][l]*y_LBFGS[w][l];
    }
    
    buf1=0.0; buf2=0.0;
    MPI_Allreduce(&dum1,&buf1,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&dum2,&buf2,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
    
    rho_LBFGS[w]=1/buf1;
    
    h0=buf1/buf2;
    
    
    /* give output so stdout */
    if(VERBOSE || 1) {
        for(w=1;w<=N_LBFGS;w++) {
            fprintf(FP,"\n rho_LBFGS(%2d)=%e",w,rho_LBFGS[w]);
        }
        fprintf(FP,"\n h0=%e\n",h0);
    }
    
    /*----------------------------------*/
    /*       L-BFGS loop 1              */
    /*----------------------------------*/
		  
    for(v=iteration-1; v>=m;v--){
        
        w=v%N_LBFGS;
        if(w==0) w=N_LBFGS;
        
        alpha_LBFGS[w]=0.0;
        for (l=1;l<=NPAR_LBFGS;l++){
            alpha_LBFGS[w]+=rho_LBFGS[w]*s_LBFGS[w][l]*q_LBFGS[l];
            
        }
        
        buf1=0.0;
        dum2=alpha_LBFGS[w];
        MPI_Allreduce(&dum2,&buf1,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
        alpha_LBFGS[w]=buf1;
        
        for (l=1;l<=NPAR_LBFGS;l++){
            q_LBFGS[l]=q_LBFGS[l]-alpha_LBFGS[w]*y_LBFGS[w][l];
            
        }
    }
    
    /*----------------------------------*/
    /*       Apply H0^-1                */
    /*----------------------------------*/
    for (l=1;l<=NPAR_LBFGS;l++){
        r_LBFGS[l]=h0*q_LBFGS[l];
    }
    
    /*----------------------------------*/
    /*       L-BFGS loop 2              */
    /*----------------------------------*/
    for(v=m; v<=iteration-1;v++){
        
        w=v%N_LBFGS;
        if(w==0) w=N_LBFGS;
        
        beta_LBFGS=0.0;
        
        for (l=1;l<=NPAR_LBFGS;l++){
            beta_LBFGS+=rho_LBFGS[w]*y_LBFGS[w][l]*r_LBFGS[l];
        }
        
        buf1=0.0;
        buf2=beta_LBFGS;
        MPI_Allreduce(&buf2,&buf1,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
        beta_LBFGS=buf1;
        
        for (l=1;l<=NPAR_LBFGS;l++){
            r_LBFGS[l]=r_LBFGS[l]+s_LBFGS[w][l]*(alpha_LBFGS[w]-beta_LBFGS);
        }
    }
    
}

void lbfgs_reset(int iter, int N_LBFGS, int NPAR_LBFGS,float ** s_LBFGS1, float ** y_LBFGS1, float * rho_LBFGS1) {
    
    /* local variables */
    int l,m;
    
    /* global variables */
    extern int NX,NY,MYID;
    
    if(MYID==0) printf("\n\n ------------ L-BFGS ---------------");
    if(MYID==0&&iter>1) printf("\n Reset L-BFGS at iteration %d",iter);
    if(MYID==0&&iter==1) printf("\n L-BFGS will be used from iteration %d on",iter+1);
    
    for(l=1;l<=N_LBFGS;l++){
        for(m=1;m<=NPAR_LBFGS*NX*NY;m++){
            s_LBFGS1[l][m]=0.0;
            y_LBFGS1[l][m]=0.0;
        }
        rho_LBFGS1[l]=0.0;
    }
}
