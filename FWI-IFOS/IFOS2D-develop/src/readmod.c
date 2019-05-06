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
 *   Read viscoelastic model properties (vp,vs,density,Qp,Qs) from files
 *
 *  ----------------------------------------------------------------------*/


/* This file contains function readmod, which has the purpose
 to read data from model-files for viscoelastic simulation */

#include "fd.h"

void readmod(float  **  rho, float **  pi, float **  u, float ** taus, float ** taup, float * eta){
    
    extern float DT, *FL, TAU;
    extern int L,WAVETYPE, VERBOSE;
    extern int NX, NY, NXG, NYG,  POS[3], MYID, PARAMETERIZATION;
    extern char  MFILE[STRING_SIZE];
    extern FILE *FP;
    
    
    /* local variables */
    float rhov, muv, piv, vp, vs, qp, qs, *pts;
    int i, j, ii, jj, l, sw_Qp=1, sw_Qs=1;
    FILE *fp_vs, *fp_vp, *fp_rho, *fp_qp, *fp_qs;
    char filename[STRING_SIZE];
    
    
    /* vector for maxwellbodies */
    pts=vector(1,L);
    for (l=1;l<=L;l++) {
        pts[l]=1.0/(2.0*PI*FL[l]);
        eta[l]=DT/pts[l];
    }
    
	   fprintf(FP,"\n...reading model information from model-files...\n");
    
	   /* read density and seismic velocities */
	   /* ----------------------------------- */
	   if(PARAMETERIZATION==1){
           
           if(WAVETYPE==1||WAVETYPE==3){
               fprintf(FP,"\t Vp:\n\t %s.vp\n\n",MFILE);
               sprintf(filename,"%s.vp",MFILE);
               fp_vp=fopen(filename,"r");
               if (fp_vp==NULL) declare_error(" Could not open model file for Vp ! ");
               
               
               fprintf(FP,"\t Qp:\n\t %s.qp\n\n",MFILE);
               sprintf(filename,"%s.qp",MFILE);
               fp_qp=fopen(filename,"r");
               if (fp_qp==NULL){
                   if (MYID==0){
                       printf(" Could not open model file for Qp-values ! \n");
                       printf(" Uses input file value for tau! \n");}
                   sw_Qp=0;
               }
               
           }
           
           
           fprintf(FP,"\t Vs:\n\t %s.vs\n\n",MFILE);
           sprintf(filename,"%s.vs",MFILE);
           fp_vs=fopen(filename,"r");
           if (fp_vs==NULL) declare_error(" Could not open model file for Vs ! ");
           
           fprintf(FP,"\t Density:\n\t %s.rho\n\n",MFILE);
           sprintf(filename,"%s.rho",MFILE);
           fp_rho=fopen(filename,"r");
           if (fp_rho==NULL) declare_error(" Could not open model file for densities ! ");
           
           fprintf(FP,"\t Qs:\n\t %s.qs\n\n",MFILE);
           sprintf(filename,"%s.qs",MFILE);
           fp_qs=fopen(filename,"r");
           if (fp_qs==NULL){
               if (MYID==0){
                   printf(" Could not open model file for Qs-values ! ");
                   printf(" Uses input file value for tau! \n");}
               sw_Qs=0;
           }
           
       }
	   
	   /* read density and Lame parameters */
	   /* ----------------------------------- */
	   if(PARAMETERIZATION==3){
           fprintf(FP,"\t Lame parameter lambda:\n\t %s.lam\n\n",MFILE);
           sprintf(filename,"%s.lam",MFILE);
           fp_vp=fopen(filename,"r");
           if (fp_vp==NULL) declare_error(" Could not open model file for Lame parameter lambda ! ");
           
           
           fprintf(FP,"\t Lame parameter mu:\n\t %s.mu\n\n",MFILE);
           sprintf(filename,"%s.mu",MFILE);
           fp_vs=fopen(filename,"r");
           if (fp_vs==NULL) declare_error(" Could not open model file for Lame parameter mu ! ");
           
           fprintf(FP,"\t Density:\n\t %s.rho\n\n",MFILE);
           sprintf(filename,"%s.rho",MFILE);
           fp_rho=fopen(filename,"r");
           if (fp_rho==NULL) declare_error(" Could not open model file for densities ! ");
           
           fprintf(FP,"\t Qp:\n\t %s.qp\n\n",MFILE);
           sprintf(filename,"%s.qp",MFILE);
           fp_qp=fopen(filename,"r");
           if (fp_qp==NULL){
               if (MYID==0){
                   printf(" Could not open model file for Qp-values ! \n");
                   printf(" Uses input file value for tau! \n");}
               sw_Qp=0;
           }
           
           fprintf(FP,"\t Qs:\n\t %s.qs\n\n",MFILE);
           sprintf(filename,"%s.qs",MFILE);
           fp_qs=fopen(filename,"r");
           if (fp_qs==NULL){
               if (MYID==0){
                   printf(" Could not open model file for Qs-values ! ");
                   printf(" Uses input file value for tau! \n");}
               sw_Qs=0;
           }
       }
	   
    
    /* loop over global grid */
    for (i=1;i<=NXG;i++){
        for (j=1;j<=NYG;j++){
            
            if(WAVETYPE==1||WAVETYPE==3){
                
                if(feof(fp_vp)){
                    declare_error("Model file VP is to small. Check dimensions NX*NY of file.");
                }
                fread(&vp, sizeof(float), 1, fp_vp);
                
                if (sw_Qp){
                    
                    if(feof(fp_qp)){
                        declare_error("Model file QP is to small. Check dimensions NX*NY of file.");
                    }
                    
                    fread(&qp, sizeof(float), 1, fp_qp);
                }
            }
            
            if(feof(fp_vs)){
                declare_error("Model file VS is to small. Check dimensions NX*NY of file.");
            }
            
            fread(&vs, sizeof(float), 1, fp_vs);
            
            if(feof(fp_rho)){
                declare_error("Model file RHO is to small. Check dimensions NX*NY of file.");
            }
            
            fread(&rhov, sizeof(float), 1, fp_rho);
            
            if (sw_Qs){
                
                if(feof(fp_vs)){
                    declare_error("Model file QS is to small. Check dimensions NX*NY of file.");
                }
                
                fread(&qs, sizeof(float), 1, fp_qs);
            }
            
            /* only the PE which belongs to the current global gridpoint
             is saving model parameters in his local arrays */
            if ((POS[1]==((i-1)/NX)) &&
                (POS[2]==((j-1)/NY))){
                ii=i-POS[1]*NX;
                jj=j-POS[2]*NY;
                
                u[jj][ii]=vs;
                rho[jj][ii]=rhov;
                if (sw_Qs){
                    taus[jj][ii]=2.0/qs;}
                else taus[jj][ii]=TAU;
                if(WAVETYPE==1||WAVETYPE==3){
                    pi[jj][ii]=vp;
                    if (sw_Qp){
                        taup[jj][ii]=2.0/qp;}
                    else taup[jj][ii]=TAU;
                }
                
                
                
            }
        }
    }
    

    
    if(WAVETYPE==1||WAVETYPE==3){
        
        fread(&vp, sizeof(float), 1, fp_vp);
        if(!feof(fp_vp)){
            declare_error("Model file VP is to big. Check dimensions NX*NY of file.");
        }
        fclose(fp_vp);
        
        if (sw_Qp) {
            fread(&qp, sizeof(float), 1, fp_qp);
            if(!feof(fp_qp)){
                declare_error("Model file QP is to big. Check dimensions NX*NY of file.");
            }
            fclose(fp_qp);
        }
    }
    
    fread(&vs, sizeof(float), 1, fp_vs);
    if(!feof(fp_vs)){
        declare_error("Model file VS is to big. Check dimensions NX*NY of file.");
    }
    fclose(fp_vs);
    
    fread(&rho, sizeof(float), 1, fp_rho);
    if(!feof(fp_rho)){
        declare_error("Model file RHO is to big. Check dimensions NX*NY of file.");
    }
    fclose(fp_rho);
    
    if (sw_Qs){
        fread(&qs, sizeof(float), 1, fp_qs);
        if(!feof(fp_qs)){
            declare_error("Model file QS is to big. Check dimensions NX*NY of file.");
        }
        fclose(fp_qs);
    }
    
    
    /* each PE writes his model to disk */
    if(WAVETYPE==1||WAVETYPE==3){
        if(PARAMETERIZATION==1) sprintf(filename,"%s.out.vp",MFILE);
        if(PARAMETERIZATION==3) sprintf(filename,"%s.out.pi",MFILE);
        write_matrix_disk(pi, filename);
        
        sprintf(filename,"%s.out.qp",MFILE);
        write_matrix_disk(taup, filename);
    }
    
    if(PARAMETERIZATION==1) sprintf(filename,"%s.out.vs",MFILE);
    if(PARAMETERIZATION==3) sprintf(filename,"%s.out.mu",MFILE);
    write_matrix_disk(u, filename);
    
    sprintf(filename,"%s.out.rho",MFILE);
    write_matrix_disk(rho, filename);
    
    sprintf(filename,"%s.out.qs",MFILE);
    write_matrix_disk(taus, filename);
    
    free_vector(pts,1,L);
}




