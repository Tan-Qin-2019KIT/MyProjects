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
/* Return values:
 * 1: condition 1 not fullfilled
 * 2: condition 2 not fullfilled
 * 0: Wolfe Condition fullfilled
 */
#include "fd.h"
int check_wolfe(float steplength, float misfit_old, float misfit_new, float ** grad_old_vs, float ** grad_new_vs, float ** update_vs, float ** grad_old_rho, float ** grad_new_rho, float ** update_rho, float ** grad_old_vp, float ** grad_new_vp, float ** update_vp, float c1, float c2, int NPAR_LBFGS){
    
    /* global */
    extern int NX,NY,MYID;
    extern FILE *FP;
    extern int ACOUSTIC;
    
    /* local */
    float left=0.0, right=0.0;
    float grad_dot_p_old=0.0,grad_dot_p_old_norm=0.0, grad_dot_p_new=0.0;
    float temp;
    float temp2;
    
    if(!ACOUSTIC){
        temp=matrix_product(grad_old_vs, update_vs);
        grad_dot_p_old+=temp;
        temp2=global_maximum(grad_old_vs);
        if(fabs(temp2)>0) grad_dot_p_old_norm+=temp/temp2;
        grad_dot_p_new+=matrix_product(grad_new_vs, update_vs);
    }
    
    if(NPAR_LBFGS>1) {
        temp=matrix_product(grad_old_rho, update_rho);
        grad_dot_p_old+=temp;
        temp2=global_maximum(grad_old_rho);
        if(fabs(temp2)>0) grad_dot_p_old_norm+=temp/temp2;
        grad_dot_p_new+=matrix_product(grad_new_rho, update_rho);
    }
    
    if(NPAR_LBFGS>2) {
        temp=matrix_product(grad_old_vp, update_vp);
        grad_dot_p_old+=temp;
        temp2=global_maximum(grad_old_vp);
        if(fabs(temp2)>0) grad_dot_p_old_norm+=temp/temp2;
        grad_dot_p_new+=matrix_product(grad_new_vp, update_vp);
    }
    
    fprintf(FP,"\n\n ------- Wolfe condition -------- ");
    
    /*---------------------------------*/
    /*       Check condition 1         */
    /* f(x+a*p) <= f(x)+c1*a*grad(x)*p */
    /*---------------------------------*/
    
    left  = misfit_new;
    right = (misfit_old-c1*steplength*grad_dot_p_old_norm*1e-3);
    
    if(left<=right) {
        fprintf(FP,"\n Condition 1 OK: (%f <= %f) ",left,right);
    } else {
        fprintf(FP,"\n Condition 1 NOT OK: (%f <= %f) ",left,right);
        return 1;
    }
    
    /*-------------------------------*/
    /*       Check condition 2       */
    /* grad(x+a*p)*p >= c2*grad(x)*p */
    /*-------------------------------*/
    
    left  = -grad_dot_p_new;
    right = -c2*grad_dot_p_old;
    
    if(left>=right) {
        fprintf(FP,"\n Condition 2 OK: (%.1f>= %.1f)", left, right);
    } else {
        fprintf(FP,"\n Condition 2 NOT OK: (%.1f>= %.1f)", left, right);
        return 2;
    }
    
    fprintf(FP,"\n Steplength %.3f satisfy wolfe condition",steplength);
    return 0;
}


void wolfe_linesearch(int wolfe_status, float *alpha_SL_min, float *alpha_SL_max, float *alpha_SL){
    
    switch (wolfe_status) {
        case 1:
            *alpha_SL_max=*alpha_SL;
            *alpha_SL=0.5*(*alpha_SL_max+*alpha_SL_min);
            break;
            
        case 2:
            *alpha_SL_min=*alpha_SL;
            if(*alpha_SL_max==0) *alpha_SL=10*(*alpha_SL);
            if(*alpha_SL_max!=0) *alpha_SL=0.5*(*alpha_SL_max+*alpha_SL_min);
            break;
            
        default:
            declare_error("\n wolfe_status not known to wolfe_linesearch");
            break;
    }
    
}


