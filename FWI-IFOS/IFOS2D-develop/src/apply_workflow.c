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
 *  Apply workflow
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void apply_workflow(float ** workflow,int workflow_lines,char workflow_header[STRING_SIZE],int *iter,float *F_LOW_PASS,int wavetype_start, int * change_wavetype_iter, int * LBFGS_iter_start){
    
    /* local variables */
    int x;
    
    /* extern variables */
    extern int INV_RHO_ITER,INV_VS_ITER,INV_VP_ITER;
    extern int TIME_FILT,MYID;
    extern float F_HIGH_PASS;
    extern float PRO;
    extern int WAVETYPE;
    extern float JOINT_INVERSION_PSV_SH_ALPHA_VS;
    extern float JOINT_INVERSION_PSV_SH_ALPHA_RHO;
    extern int EPRECOND;
    extern float EPSILON_WE;
    extern float GAMMA;
    extern int GRAD_METHOD;
    extern int WORKFLOW_STAGE;
    
    /******************/
    /* Apply Workflow */
    /******************/
    
    /* Print current workflow */
    if(MYID==0){
        printf("\n ---------- Applying Workflow -----------\n");
        printf(" %s ",workflow_header);
        for(x=1;x<=WORKFLOW_MAX_VAR;x++){
            printf("%.2f\t",workflow[WORKFLOW_STAGE][x]);
        }
    }
    
    /* Inversion of material parameter */
    if(workflow[WORKFLOW_STAGE][2]!=-1) {
        if(workflow[WORKFLOW_STAGE][2]==1) {
            if(INV_VS_ITER>*iter) INV_VS_ITER=*iter;
        } else {
            /* detect change and reset LBFGS */
            if(INV_VS_ITER<*iter) *LBFGS_iter_start=*iter;
            INV_VS_ITER=*iter+10;
        }
    }
    
    if(workflow[WORKFLOW_STAGE][3]!=-1){
        if(workflow[WORKFLOW_STAGE][3]==1) {
            if(INV_VP_ITER>*iter) INV_VP_ITER=*iter;
        } else {
            /* detect change and reset LBFGS */
            if(INV_VP_ITER<*iter) *LBFGS_iter_start=*iter;
            INV_VP_ITER=*iter+10;
        }
    }
    
    if(workflow[WORKFLOW_STAGE][4]!=-1){
        if(workflow[WORKFLOW_STAGE][4]==1) {
            if(INV_RHO_ITER>*iter) INV_RHO_ITER=*iter;
        } else {
            /* detect change and reset LBFGS */
            if(INV_RHO_ITER<*iter) *LBFGS_iter_start=*iter;
            INV_RHO_ITER=*iter+10;
        }
    }
    
    /* Abort criterium */
    PRO=workflow[WORKFLOW_STAGE][5];
    
    /* Frequency filtering  */
    if( TIME_FILT == 1 ) {
        
        TIME_FILT=workflow[WORKFLOW_STAGE][6];
        
        if( TIME_FILT > 0 ) {
            if( F_HIGH_PASS != workflow[WORKFLOW_STAGE][7] ) *LBFGS_iter_start=*iter;
            F_HIGH_PASS=workflow[WORKFLOW_STAGE][7];
            
            if( *F_LOW_PASS != workflow[WORKFLOW_STAGE][8] ) *LBFGS_iter_start=*iter;
            *F_LOW_PASS=workflow[WORKFLOW_STAGE][8];
        }
        
    } else {
        if(MYID==0&&(workflow[WORKFLOW_STAGE][6]>0))printf("\n TIME_FILT cannot be activated due to it is not activated in the JSON File \n");
    }
    
    /* Change of wavetype  */
    if(wavetype_start!=3&&(WAVETYPE!=workflow[WORKFLOW_STAGE][9])){
        if(MYID==0)printf("\n Sorry, change of WAVETYPE with workflow only possible if WAVETYPE==3 in *.json");
        if(MYID==0)printf("\n WAVETYPE will remain unchanged %i",WAVETYPE);
    } else {
        /* detect change and reset some things */
        if(WAVETYPE!=workflow[WORKFLOW_STAGE][9]) {
            *change_wavetype_iter=*iter;
            *LBFGS_iter_start=*iter;
        }
        WAVETYPE=workflow[WORKFLOW_STAGE][9];
    }
    
    /* Joint inversion PSV and SH  */
    JOINT_INVERSION_PSV_SH_ALPHA_VS=workflow[WORKFLOW_STAGE][10];
    JOINT_INVERSION_PSV_SH_ALPHA_RHO=workflow[WORKFLOW_STAGE][11];
    
    /* Approx. Hessian  */
    if(EPRECOND==0 && workflow[WORKFLOW_STAGE][12]!=0){
        if(MYID==0) printf(" WARNING: EPRECOND have to be set >0 in JSON (if so, ignore this message)");
    }
    
    EPRECOND=workflow[WORKFLOW_STAGE][12];
    EPSILON_WE=workflow[WORKFLOW_STAGE][13];
    GAMMA=workflow[WORKFLOW_STAGE][14];
    
    if(*LBFGS_iter_start==*iter && GRAD_METHOD==2){
        if(MYID==0)printf("\n L-BFGS will be used from iteration %d on.",*LBFGS_iter_start+1);
    }
}
