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
 *   Read elastic model properties (vp,vs,density) from files
 *
 *  ----------------------------------------------------------------------*/


/* This file contains function readmod, which has the purpose
 to read data from model-files for viscoelastic simulation */

#include "fd.h"

void readmod_elastic(float  **  rho, float **  pi, float **  u){
    
    extern int NX, NY, NXG, NYG,  POS[3], MYID, PARAMETERIZATION,WAVETYPE;
    extern char  MFILE[STRING_SIZE];
    extern FILE *FP;
    
    
    /* local variables */
    float rhov, vp = 0.0, vs;
    int i, j, ii, jj;
    FILE *fp_vs, *fp_vp = NULL, *fp_rho;
    char filename[STRING_SIZE];
    
    
    
    
	   fprintf(FP,"\n...reading model information from model-files...\n");
       /* ----------------------------------- */
	   /* read density and seismic velocities */
	   /* ----------------------------------- */
	   if(PARAMETERIZATION==1){
           
           if(WAVETYPE==1||WAVETYPE==3){
               fprintf(FP,"\t Vp:\n\t %s.vp\n\n",MFILE);
               sprintf(filename,"%s.vp",MFILE);
               fp_vp=fopen(filename,"r");
               if (fp_vp==NULL) declare_error(" Could not open model file for Vp ! ");
           }
           
           
           fprintf(FP,"\t Vs:\n\t %s.vs\n\n",MFILE);
           sprintf(filename,"%s.vs",MFILE);
           fp_vs=fopen(filename,"r");
           if (fp_vs==NULL) declare_error(" Could not open model file for Vs ! ");
           
           fprintf(FP,"\t Density:\n\t %s.rho\n\n",MFILE);
           sprintf(filename,"%s.rho",MFILE);
           fp_rho=fopen(filename,"r");
           if (fp_rho==NULL) declare_error(" Could not open model file for densities ! ");
           
           
       } else {
	   
	   /* read density and Lame parameters */
	   /* ----------------------------------- */
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
       }
	   
    
    /* loop over global grid */
    for (i=1;i<=NXG;i++){
        for (j=1;j<=NYG;j++){
            
            if(WAVETYPE==1||WAVETYPE==3){
                
                if(feof(fp_vp)){
                    declare_error("Model file VP is to small. Check dimensions NX*NY of file.");
                }
                
                fread(&vp, sizeof(float), 1, fp_vp);
                
            }
            
            if(feof(fp_vs) || feof(fp_rho)){
                declare_error("Model file VS or RHO is to small. Check dimensions NX*NY of file.");
            }
            
            fread(&vs, sizeof(float), 1, fp_vs);
            fread(&rhov, sizeof(float), 1, fp_rho);
            
            
            /* only the PE which belongs to the current global gridpoint
             is saving model parameters in his local arrays */
            if ((POS[1]==((i-1)/NX)) &&
                (POS[2]==((j-1)/NY))){
                ii=i-POS[1]*NX;
                jj=j-POS[2]*NY;
                
                u[jj][ii]=vs;
                rho[jj][ii]=rhov;
                if(WAVETYPE==1||WAVETYPE==3){
                    pi[jj][ii]=vp;
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
    }
    
    fread(&vs, sizeof(float), 1, fp_vs);
    fread(&rho, sizeof(float), 1, fp_rho);
    if(!feof(fp_vs) || !feof(fp_rho)){
        declare_error("Model file VS or RHO is to big. Check dimensions NX*NY of file.");
    }
    fclose(fp_vs);
    fclose(fp_rho);
    
    
    
    /* each PE writes his model to disk */
    if(WAVETYPE==1||WAVETYPE==3){
        if(PARAMETERIZATION==1) sprintf(filename,"%s.out.vp",MFILE);
        if(PARAMETERIZATION==3) sprintf(filename,"%s.out.pi",MFILE);
        write_matrix_disk(pi, filename);
    }
    
    if(PARAMETERIZATION==1) sprintf(filename,"%s.out.vs",MFILE);
    if(PARAMETERIZATION==3) sprintf(filename,"%s.out.mu",MFILE);
    write_matrix_disk(u, filename);
    
    sprintf(filename,"%s.out.rho",MFILE);
    write_matrix_disk(rho, filename);
    
}




