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
 *   Calculate energy of measured data
 *  ----------------------------------------------------------------------*/
#include "fd.h"

double calc_energy(float **sectiondata, int ntr, int ns, float energy, int ntr_glob, int **recpos_loc, int nsrc_glob, int ishot, int iter, float ** srcpos, int ** recpos){
    
    /* declaration of variables */
    extern float DT;
    extern int MYID, USE_WORKFLOW, WORKFLOW_STAGE;
    extern int TRKILL, NORMALIZE, F_LOW_PASS, TIMEWIN;
    extern char TRKILL_FILE[STRING_SIZE];
    extern int VELOCITY;
    int i, j;
    float energy_dummy;
    int h;
    float intseis_sd;
    float **intseis_sectiondata=NULL;
    float *picked_times=NULL;
    
    intseis_sectiondata = matrix(1,ntr,1,ns); /* declaration of variables for integration */
    if(TIMEWIN) picked_times = vector(1,ntr); /* declaration of variables for TIMEWIN */
    
    /*extern int MYID;*/
    
    /*printf("MYID: %d; ns=%d\n",MYID,ns);
     printf("MYID: %d; ntr=%d\n",MYID,ntr);*/
    
    /* TRACE KILLING */
    int ** kill_tmp, *kill_vector;	/* declaration of variables for trace killing */
    char trace_kill_file[STRING_SIZE];
    FILE *ftracekill;
    extern int TRKILL_OFFSET;
    extern float TRKILL_OFFSET_LOWER;
    extern float TRKILL_OFFSET_UPPER;
    
    /*-----------------*/
    /* Start Tracekill */
    /*-----------------*/
    if(TRKILL){
        kill_tmp = imatrix(1,ntr_glob,1,nsrc_glob);
        kill_vector = ivector(1,ntr);
        
        /*------------------*/
        /* clear kill table */
        /*------------------*/
        for(i=1;i<=nsrc_glob;i++){
            for (j=1; j<=ntr_glob; j++) {
                kill_tmp[j][i]=0;
            }
        }
        
        /* Use ONLY offset based TraceKill */
        if(TRKILL_OFFSET==1) {
            /* Generate TraceKill file on the fly based on the offset from the source */
            create_trkill_table(kill_tmp,ntr_glob,recpos,nsrc_glob,srcpos,ishot,TRKILL_OFFSET_LOWER,TRKILL_OFFSET_UPPER);
        } else {
            
            /* READ TraceKill file from disk */
            
            if(USE_WORKFLOW){
                sprintf(trace_kill_file,"%s_%i.dat",TRKILL_FILE,WORKFLOW_STAGE);
                ftracekill=fopen(trace_kill_file,"r");
                if (ftracekill==NULL){
                    /* If Workflow TraceKill file not found use File without workflow extensions */
                    sprintf(trace_kill_file,"%s.dat",TRKILL_FILE);
                    ftracekill=fopen(trace_kill_file,"r");
                    if (ftracekill==NULL){
                        declare_error(" Trace kill file could not be opened!");
                    }
                }
            }else{
                sprintf(trace_kill_file,"%s.dat",TRKILL_FILE);
                ftracekill=fopen(trace_kill_file,"r");
                if (ftracekill==NULL){
                    declare_error(" Trace kill file could not be opened!");
                }
            }
            
            for(i=1;i<=ntr_glob;i++){
                for(j=1;j<=nsrc_glob;j++){
                    fscanf(ftracekill,"%d",&kill_tmp[i][j]);
                    if(feof(ftracekill)){
                        declare_error(" Error while reading TraceKill file. Check dimensions!");
                    }
                }
            }
            
            fclose(ftracekill);
            
            /* Use Tracekill FILE and add the offset based TraceKill */
            if(TRKILL_OFFSET==2) {
                create_trkill_table(kill_tmp,ntr_glob,recpos,nsrc_glob,srcpos,ishot,TRKILL_OFFSET_LOWER,TRKILL_OFFSET_UPPER);
            }
        }
        
        /* Generate local kill vector */
        h=1;
        for(i=1;i<=ntr;i++){
            kill_vector[h] = kill_tmp[recpos_loc[3][i]][ishot];
            h++;
        }
    }
    /*---------------*/
    /* End Tracekill */
    /*---------------*/
    
    /* Integration of measured data  */
    for(i=1;i<=ntr;i++){
        
        intseis_sd=0;
        
        if (VELOCITY==0){	/* only integtration if displacement seismograms are inverted */
            for(j=1;j<=ns;j++){
                intseis_sd += sectiondata[i][j];
                intseis_sectiondata[i][j]=intseis_sd*DT;
            }
        }
        else{
            for(j=1;j<=ns;j++){
                intseis_sectiondata[i][j]=sectiondata[i][j];
            }
        }
    }
    
    /* TIME WINDOWING */
    if(TIMEWIN==1) time_window(intseis_sectiondata, iter, ntr_glob,recpos_loc, ntr, ns, ishot);
    
    /* NORMALIZE TRACES */
    if(NORMALIZE==1) normalize_data(intseis_sectiondata,ntr,ns);
    
    
    /* calculate energy of measured data */
    for(i=1;i<=ntr;i++){
        if((TRKILL==1)&&(kill_vector[i]==1))
            continue;
        
        for(j=1;j<=ns;j++){
            energy += intseis_sectiondata[i][j]*intseis_sectiondata[i][j];
        }
    }
    
    energy_dummy=energy;
    
    /* printf("\n MYID = %i   IN CALC_ENERGY: energy = %10.12f \n",MYID,energy); */
    
    /* free memory for integrated seismogram */
    free_matrix(intseis_sectiondata,1,ntr,1,ns);
    
    /* free memory for time windowing and trace killing */
    if(TIMEWIN==1) free_vector(picked_times,1,ntr);
    if(TRKILL==1){
        free_imatrix(kill_tmp,1,ntr_glob,1,nsrc_glob);
        free_ivector(kill_vector,1,ntr);}
    
    return energy_dummy;
} /* end of function */
