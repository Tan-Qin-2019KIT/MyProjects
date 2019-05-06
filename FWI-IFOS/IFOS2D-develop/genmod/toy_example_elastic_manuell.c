/*-----------------------------------------------------------------------------------------
 * Copyright (C) 2016
 *
 * This file is part of IFOS2D.
 *
 * IFOS2D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 *
 * IFOS2D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with IFOS2D. See file COPYING and/or <http://www.gnu.org/licenses/gpl-2.0.html>.
 -----------------------------------------------------------------------------------------*/


#include "fd.h"

void model_elastic(float  **  rho, float **  pi, float **  u){
    
    /*--------------------------------------------------------------------------*/
    /* extern variables */
    
    extern int NX, NY, NXG, NYG,  POS[3], MYID, L;
    extern float DH, DT, *FL, TAU;
    extern char  MFILE[STRING_SIZE];
    extern char INV_MODELFILE[STRING_SIZE];
    
    /* local variables */
    float muv, piv, vp, vs, rhov, ts, tp, *pts;
    int i, j, ii, jj, l;
    char modfile[STRING_SIZE];
    
    FILE *flfile;
    int nodes;
    char cline[256];
    
    float *fldepth, *flrho, *flvp, *flvs;
    
    
    
    
    /**************************************************/
    /* creation of shear wave velocity and TAU models */
    /**************************************************/
    
    /*read FL nodes from File*/
    nodes=3;
    fldepth=vector(1,nodes);
    flvs=vector(1,nodes);
    flrho=vector(1,nodes);
    flvp=vector(1,nodes);
    flfile=fopen("model_true/flnodes.toy_example.start","r");
    if (flfile==NULL) declare_error(" FL-file could not be opened !");
    
    
    
    for (l=1;l<=nodes;l++){
        fgets(cline,255,flfile);
        if (cline[0]!='#'){
            sscanf(cline,"%f%f%f%f",&fldepth[l], &flrho[l], &flvp[l], &flvs[l]);
        }
        else l=l-1;
        
    }
    
    if(MYID==0){
        printf(" ------------------------------------------------------------------ \n\n");
        printf(" Information of FL nodes: \n\n");
        printf(" \t depth \t vs \n\n");
        
        for (l=1;l<=nodes;l++){
            printf(" \t %f \t %f\n\n",fldepth[l],flvs[l]);
        }
        printf(" ------------------------------------------------------------------ \n\n");
    }
    
    /* loop over global grid */
    for (i=1;i<=NXG;i++){
        for (l=1;l<nodes;l++){
            if (fldepth[l]==fldepth[l+1]){
                if ((i==1) && (MYID==0)){
                    printf("depth: %f m: double node\n",fldepth[l]);}}
            else{
                for (j=(int)(fldepth[l]/DH)+1;j<=(int)(fldepth[l+1]/DH);j++){
                    
                    vs=0.0;
                    
                    vs=(DH*(j-1)-fldepth[l])*(flvs[l+1]-flvs[l])/(fldepth[l+1]-fldepth[l])+flvs[l];
                    vs=vs*1000.0;
                    
                    muv=vs;
                    ts=TAU;
                    tp=TAU;
                    
                    /* only the PE which belongs to the current global gridpoint
                     is saving model parameters in his local arrays */
                    if ((POS[1]==((i-1)/NX)) &&
                        (POS[2]==((j-1)/NY))){
                        ii=i-POS[1]*NX;
                        jj=j-POS[2]*NY;
                        
                        u[jj][ii]=muv;
                    }
                }
            }
        }
        
        for (j=(int)(fldepth[nodes]/DH)+1;j<=NYG;j++){
            
            vs=0.0;
            vs=flvs[nodes]*1000.0;
            
            muv=vs;
            
            /* only the PE which belongs to the current global gridpoint
             is saving model parameters in his local arrays */
            if ((POS[1]==((i-1)/NX)) &&
                (POS[2]==((j-1)/NY))){
                ii=i-POS[1]*NX;
                jj=j-POS[2]*NY;
                
                u[jj][ii]=muv;
            }
        }
    }
    
    free_vector(fldepth,1,nodes);
    free_vector(flrho,1,nodes);
    free_vector(flvp,1,nodes);
    free_vector(flvs,1,nodes);
    
    
    
    /**************************************************/
    /* creation of P wave velocity and density models */
    /**************************************************/
    
    /*read FL nodes from File*/
    nodes=7;
    fldepth=vector(1,nodes);
    flvs=vector(1,nodes);
    flrho=vector(1,nodes);
    flvp=vector(1,nodes);
    flfile=fopen("model_true/flnodes.toy_example","r");
    flfile=fopen("model_true/flnodes.toy_example.start","r");
    if (flfile==NULL) declare_error(" FL-file could not be opened !");
    
    for (l=1;l<=nodes;l++){
        fgets(cline,255,flfile);
        if (cline[0]!='#'){
            sscanf(cline,"%f%f%f%f",&fldepth[l], &flrho[l], &flvp[l], &flvs[l]);
        }
        else l=l-1;
    }
    
    
    if(MYID==0){
        printf(" ------------------------------------------------------------------ \n\n");
        printf(" Information of FL nodes: \n\n");
        printf(" \t depth \t vp \t rho \n\n");
        
        for (l=1;l<=nodes;l++){
            printf(" \t %f \t %f \t %f \n\n",fldepth[l],flvp[l],flrho[l]);
        }
        printf(" ------------------------------------------------------------------ \n\n");
    }
    
    /* loop over global grid */
    for (i=1;i<=NXG;i++){
        for (l=1;l<nodes;l++){
            if (fldepth[l]==fldepth[l+1]){
                if ((i==1) && (MYID==0)){
                    printf("depth: %f m: double node\n",fldepth[l]);}}
            else{
                for (j=(int)(fldepth[l]/DH)+1;j<=(int)(fldepth[l+1]/DH);j++){
                    
                    vp=0.0; rhov=0.0;
                    
                    vp=(DH*(j-1)-fldepth[l])*(flvp[l+1]-flvp[l])/(fldepth[l+1]-fldepth[l])+flvp[l];
                    vp=vp*1000.0;
                    rhov=(DH*(j-1)-fldepth[l])*(flrho[l+1]-flrho[l])/(fldepth[l+1]-fldepth[l])+flrho[l];
                    rhov=rhov*1000.0;
                    
                    piv=vp;
                    
                    /* only the PE which belongs to the current global gridpoint
                     is saving model parameters in his local arrays */
                    if ((POS[1]==((i-1)/NX)) &&
                        (POS[2]==((j-1)/NY))){
                        ii=i-POS[1]*NX;
                        jj=j-POS[2]*NY;
                        
                        rho[jj][ii]=rhov;
                        pi[jj][ii]=piv;
                    }
                }
            }
        }
        
        for (j=(int)(fldepth[nodes]/DH)+1;j<=NYG;j++){
            
            vp=0.0; rhov=0.0;
            vp=flvp[nodes]*1000.0; rhov=flrho[nodes]*1000.0;
            
            piv=vp;
            
            /* only the PE which belongs to the current global gridpoint 
             is saving model parameters in his local arrays */
            if ((POS[1]==((i-1)/NX)) && 
                (POS[2]==((j-1)/NY))){
                ii=i-POS[1]*NX;
                jj=j-POS[2]*NY;
                
                rho[jj][ii]=rhov;
                pi[jj][ii]=piv;
            }
        }
    }
    
//    free_vector(fldepth,1,nodes);
//    free_vector(flrho,1,nodes);
//    free_vector(flvp,1,nodes);
//    free_vector(flvs,1,nodes);
//    
    
    
    
    sprintf(modfile,"%s_rho_it0.bin",INV_MODELFILE);
    writemod(modfile,rho,3);
    MPI_Barrier(MPI_COMM_WORLD);
    if (MYID==0) mergemod(modfile,3);
    MPI_Barrier(MPI_COMM_WORLD); 
    sprintf(modfile,"%s_rho_it0.bin.%i%i",INV_MODELFILE,POS[1],POS[2]);
    remove(modfile);
    
    sprintf(modfile,"%s_vs_it0.bin",INV_MODELFILE);
    writemod(modfile,u,3);
    MPI_Barrier(MPI_COMM_WORLD);
    if (MYID==0) mergemod(modfile,3);
    MPI_Barrier(MPI_COMM_WORLD); 
    sprintf(modfile,"%s_vs_it0.bin.%i%i",INV_MODELFILE,POS[1],POS[2]);
    remove(modfile);
    
    sprintf(modfile,"%s_vp_it0.bin",INV_MODELFILE);
    writemod(modfile,pi,3);
    MPI_Barrier(MPI_COMM_WORLD);
    if (MYID==0) mergemod(modfile,3);
    MPI_Barrier(MPI_COMM_WORLD); 
    sprintf(modfile,"%s_vp_it0.bin.%i%i",INV_MODELFILE,POS[1],POS[2]);
    remove(modfile);
    
    
    
    
//    free_vector(pts,1,L);
    
}



