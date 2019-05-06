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
 *   Calculate number of killed traces
 *  ----------------------------------------------------------------------*/
#include "fd.h"

void count_killed_traces(int ntr, int swstestshot, int ntr_glob, int **recpos_loc, int nsrc_glob, int ishot, int* ptr_killed_traces, int* ptr_killed_traces_testshots,float ** srcpos, int ** recpos){
    
    /* declaration of variables */
    extern int USE_WORKFLOW, WORKFLOW_STAGE;
    extern char TRKILL_FILE[STRING_SIZE];
    extern int TRKILL_OFFSET;
    extern float TRKILL_OFFSET_LOWER;
    extern float TRKILL_OFFSET_UPPER;
    extern int MYID;
    int i,j,h;
    
    
    /* declaration of variables for trace killing */
    int ** kill_tmp, *kill_vector;
    char trace_kill_file[STRING_SIZE];
    FILE *ftracekill;
    
    
    
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
        if(MYID==0) {
            printf("\n\n ----- Offset based Tracekill ------");
            printf("\n Kill offsets between %.1f m and %.1f m",TRKILL_OFFSET_LOWER,TRKILL_OFFSET_UPPER);
        }
        
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
            create_trkill_table(kill_tmp,ntr_glob,recpos,nsrc_glob,srcpos,-100,TRKILL_OFFSET_LOWER,TRKILL_OFFSET_UPPER);
            if(MYID==0) {
                printf("\n\n ----- Offset based Tracekill ------");
                printf("\n Kill offsets between %.1f m and %.1f m",TRKILL_OFFSET_LOWER,TRKILL_OFFSET_UPPER);
                printf("\n In addition the TraceKill File is used");
            }
            
        }
    }
    
    
    h=1;
    for(i=1;i<=ntr;i++){
        kill_vector[h] = kill_tmp[recpos_loc[3][i]][ishot];
        if (kill_vector[h]==1) *ptr_killed_traces=*ptr_killed_traces+1;
        if ((kill_vector[h]==1) && (swstestshot==1)) *ptr_killed_traces_testshots=*ptr_killed_traces_testshots+1;
        h++;
    }
}

