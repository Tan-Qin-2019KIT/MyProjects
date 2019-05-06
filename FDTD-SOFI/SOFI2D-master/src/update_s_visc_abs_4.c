/*------------------------------------------------------------------------
 * Copyright (C) 2015 For the list of authors, see file AUTHORS.
 *
 * This file is part of SOFI2D.
 *
 * SOFI2D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 *
 * SOFI2D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SOFI2D. See file COPYING and/or
 * <http://www.gnu.org/licenses/gpl-2.0.html>.
 --------------------------------------------------------------------------*/

/* $Id: update_s_visc_abs.c 819 2015-04-17 11:07:06Z tmetz $*/
/*------------------------------------------------------------------------
 *   updating stress components at gridpoints of the absorbing frame (ABS=2 in the json file)
 *   by a staggered grid finite difference scheme of FDORDER accuracy in space
 *   and second order accuracy in time
 *
 *   gx and gy are arrays with the locations of the boundary specified in subgrid_bounds.c
 *   for each subgrid
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void update_s_visc_abs_4(int nx1, int nx2, int ny1, int ny2, int *gx, int *gy, int nt,
                       float **vx, float **vy, float **sxx, float **syy,
                       float **sxy, float ***r, float *** p, float ***q,
                       float **pi,
                       float ** fipjp, float **f, float **g, float *bip, float *bjm, float *cip,
                       float *cjm, float ***d, float ***e, float ***dip,
                       float ** absorb_coeff,float *hc,  float ** vxx_1,float ** vxx_2,float ** vxx_3,float ** vxx_4,float ** vyy_1,float ** vyy_2,float ** vyy_3,float ** vyy_4,float ** vxy_1,float ** vxy_2,float ** vxy_3,float ** vxy_4,float ** vyx_1,float ** vyx_2,float ** vyx_3,float ** vyx_4,float ** svx_1,float ** svx_2,float ** svx_3,float ** svx_4,float ** svy_1,float ** svy_2,float ** svy_3,float ** svy_4,float ***pr_2,float ***pr_3,float ***pr_4, float ***pp_2, float ***pp_3, float ***pp_4, float ***pq_2, float ***pq_3, float ***pq_4){
    
    
    int i,j,fdoh;
    float  vxx, vyy, vxy, vyx;
    extern int MYID, FDORDER;
    extern FILE *FP;
    extern int OUTNTIMESTEPINFO;
    double time1=0.0, time2=0.0;
    
    
    
    fdoh = FDORDER/2;
    
    /*Pointer array to the locations of the fd-operator functions*/
    void ( *FD_op_s[7] ) ();
    FD_op_s[1] = &operator_s_fd2;
    FD_op_s[2] = &operator_s_fd4;
    FD_op_s[3] = &operator_s_fd6;
    FD_op_s[4] = &operator_s_fd8;
    FD_op_s[5] = &operator_s_fd10;
    FD_op_s[6] = &operator_s_fd12;
    
    if ((MYID==0) && ((nt+(OUTNTIMESTEPINFO-1))%OUTNTIMESTEPINFO)==0) {
        time1=MPI_Wtime();
        fprintf(FP,"\n **Message from update_s_4 (printed by PE %d):\n",MYID);
        fprintf(FP," Updating stress components ...");
    }
    
    
    for ( j=gy[2]+1; j<=gy[3]; j++ ) {
        for ( i=gx[1]; i<=gx[2]; i++ ) {
            FD_op_s[fdoh] ( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy,hc );
            
            wavefield_update_s_visc_4 ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,r,p,q,
                                       fipjp,f,g,bip,bjm,cip,cjm,d,e,dip,vxx_1,vxx_2,vxx_3,vxx_4,vyy_1,vyy_2,vyy_3,vyy_4,vxy_1,vxy_2,vxy_3,vxy_4,vyx_1,vyx_2,vyx_3,vyx_4,svx_1,svx_2,svx_3,svx_4,svy_1,svy_2,svy_3,svy_4,pr_2,pr_3,pr_4,pp_2,pp_3,pp_4,pq_2,pq_3,pq_4);
            abs_update_s (i,j,sxx,sxy,syy,absorb_coeff );
        }
    }
    
    /* right boundary */
    for ( j=gy[2]+1; j<=gy[3]; j++ ) {
        for ( i=gx[3]+1; i<=gx[4]; i++ ) {
            FD_op_s[fdoh] ( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy,hc );
            
            wavefield_update_s_visc_4 ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,r,p,q,
                                       fipjp,f,g,bip,bjm,cip,cjm,d,e,dip,vxx_1,vxx_2,vxx_3,vxx_4,vyy_1,vyy_2,vyy_3,vyy_4,vxy_1,vxy_2,vxy_3,vxy_4,vyx_1,vyx_2,vyx_3,vyx_4,svx_1,svx_2,svx_3,svx_4,svy_1,svy_2,svy_3,svy_4,pr_2,pr_3,pr_4,pp_2,pp_3,pp_4,pq_2,pq_3,pq_4);
            abs_update_s (i,j,sxx,sxy,syy,absorb_coeff );
            
        }
    }
    
    /* top boundary */
    for ( j=gy[1]; j<=gy[2]; j++ ) {
        for ( i=gx[2]+1; i<=gx[3]; i++ ) {
            FD_op_s[fdoh] ( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy,hc );
            wavefield_update_s_visc_4 ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,r,p,q,
                                       fipjp,f,g,bip,bjm,cip,cjm,d,e,dip,vxx_1,vxx_2,vxx_3,vxx_4,vyy_1,vyy_2,vyy_3,vyy_4,vxy_1,vxy_2,vxy_3,vxy_4,vyx_1,vyx_2,vyx_3,vyx_4,svx_1,svx_2,svx_3,svx_4,svy_1,svy_2,svy_3,svy_4,pr_2,pr_3,pr_4,pp_2,pp_3,pp_4,pq_2,pq_3,pq_4);
            abs_update_s (i,j,sxx,sxy,syy,absorb_coeff );
        }
    }
    
    /* bottom boundary */
    for ( j=gy[3]+1; j<=gy[4]; j++ ) {
        for ( i=gx[2]+1; i<=gx[3]; i++ ) {
            FD_op_s[fdoh] ( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy,hc );
            
            wavefield_update_s_visc_4 ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,r,p,q,
                                       fipjp,f,g,bip,bjm,cip,cjm,d,e,dip,vxx_1,vxx_2,vxx_3,vxx_4,vyy_1,vyy_2,vyy_3,vyy_4,vxy_1,vxy_2,vxy_3,vxy_4,vyx_1,vyx_2,vyx_3,vyx_4,svx_1,svx_2,svx_3,svx_4,svy_1,svy_2,svy_3,svy_4,pr_2,pr_3,pr_4,pp_2,pp_3,pp_4,pq_2,pq_3,pq_4);
            
            abs_update_s (i,j,sxx,sxy,syy,absorb_coeff );
        }
    }
    
    /* corners */
    
    /*left-top*/
    for ( j=gy[1]; j<=gy[2]; j++ ) {
        for ( i=gx[1]; i<=gx[2]; i++ ) {
            FD_op_s[fdoh] ( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy,hc );
            
            wavefield_update_s_visc_4 ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,r,p,q,
                                       fipjp,f,g,bip,bjm,cip,cjm,d,e,dip,vxx_1,vxx_2,vxx_3,vxx_4,vyy_1,vyy_2,vyy_3,vyy_4,vxy_1,vxy_2,vxy_3,vxy_4,vyx_1,vyx_2,vyx_3,vyx_4,svx_1,svx_2,svx_3,svx_4,svy_1,svy_2,svy_3,svy_4,pr_2,pr_3,pr_4,pp_2,pp_3,pp_4,pq_2,pq_3,pq_4);
            
            abs_update_s (i,j,sxx,sxy,syy,absorb_coeff );
            
        }
    }
    
    /*left-bottom*/
    for ( j=gy[3]+1; j<=gy[4]; j++ ) {
        for ( i=gx[1]; i<=gx[2]; i++ ) {
            FD_op_s[fdoh] ( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy,hc );
            
            wavefield_update_s_visc_4 ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,r,p,q,
                                       fipjp,f,g,bip,bjm,cip,cjm,d,e,dip,vxx_1,vxx_2,vxx_3,vxx_4,vyy_1,vyy_2,vyy_3,vyy_4,vxy_1,vxy_2,vxy_3,vxy_4,vyx_1,vyx_2,vyx_3,vyx_4,svx_1,svx_2,svx_3,svx_4,svy_1,svy_2,svy_3,svy_4,pr_2,pr_3,pr_4,pp_2,pp_3,pp_4,pq_2,pq_3,pq_4);
            
            abs_update_s (i,j,sxx,sxy,syy,absorb_coeff );
        }
    }
    
    /* right-top */
    for ( j=gy[1]; j<=gy[2]; j++ ) {
        for ( i=gx[3]+1; i<=gx[4]; i++ ) {
            FD_op_s[fdoh] ( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy,hc );
            
            wavefield_update_s_visc_4 ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,r,p,q,
                                       fipjp,f,g,bip,bjm,cip,cjm,d,e,dip,vxx_1,vxx_2,vxx_3,vxx_4,vyy_1,vyy_2,vyy_3,vyy_4,vxy_1,vxy_2,vxy_3,vxy_4,vyx_1,vyx_2,vyx_3,vyx_4,svx_1,svx_2,svx_3,svx_4,svy_1,svy_2,svy_3,svy_4,pr_2,pr_3,pr_4,pp_2,pp_3,pp_4,pq_2,pq_3,pq_4);
            
            abs_update_s (i,j,sxx,sxy,syy,absorb_coeff );
        }
    }
    
    /* right-bottom */
    for ( j=gy[3]+1; j<=gy[4]; j++ ) {
        for ( i=gx[3]+1; i<=gx[4]; i++ ) {
            FD_op_s[fdoh] ( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy,hc );
            
            wavefield_update_s_visc_4 ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,r,p,q,
                                       fipjp,f,g,bip,bjm,cip,cjm,d,e,dip,vxx_1,vxx_2,vxx_3,vxx_4,vyy_1,vyy_2,vyy_3,vyy_4,vxy_1,vxy_2,vxy_3,vxy_4,vyx_1,vyx_2,vyx_3,vyx_4,svx_1,svx_2,svx_3,svx_4,svy_1,svy_2,svy_3,svy_4,pr_2,pr_3,pr_4,pp_2,pp_3,pp_4,pq_2,pq_3,pq_4);
            
            abs_update_s (i,j,sxx,sxy,syy,absorb_coeff );
        }
    }
    
    
    if ((MYID==0) && ((nt+(OUTNTIMESTEPINFO-1))%OUTNTIMESTEPINFO)==0) {
        time2=MPI_Wtime();
        fprintf(FP," finished (real time: %4.3f s).\n",time2-time1);
    }
}
