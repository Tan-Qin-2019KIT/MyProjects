/*-----------------------------------------------------------------------------------------
 * Copyright (C) 2015  For the list of authors, see file AUTHORS.
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
 * along with SOFI2D See file COPYING and/or <http://www.gnu.org/licenses/gpl-2.0.html>.
 -----------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
 *   updating particle velocities at gridpoints of the CPML-frame (ABS=1 in the json file)
 *   by a staggered grid finite difference scheme of FDORDER accuracy in space
 *   and second order accuracy in time
 *   T. Bohlen
 *
 *   gx and gy are arrays with the locations of the boundary specified in subgrid_bounds.c
 *   for each subgrid
 *  ----------------------------------------------------------------------*/

#include "fd.h"




void update_v_PML_4 ( int nx1, int nx2, int ny1, int ny2, int * gx, int * gy, int nt, float **  vx, float ** vy,
                   float ** sxx, float ** syy, float ** sxy,  float  **rip, float **rjp,
                   float *hc, float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half,
                   float * b_x_half, float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,
                   float ** psi_sxx_x, float ** psi_syy_y, float ** psi_sxy_y, float ** psi_sxy_x,float ** svx_1,float ** svx_2,float ** svx_3,float ** svx_4,float ** svy_1,float ** svy_2,float ** svy_3,float ** svy_4)
{
    
    int i, j, h1,fdoh;
    float sxx_x, syy_y, sxy_y, sxy_x;
    
    extern int OUTNTIMESTEPINFO;
    double time1=0.0, time2=0.0;
    
    extern int MYID, FDORDER;
    extern int  FW;
    extern FILE *FP;
    
    fdoh=FDORDER/2;
    
    /*Pointer array to the locations of the fd-operator functions*/
    void ( *FD_op_v[7] ) ();
    FD_op_v[1] = &operator_v_fd2;
    FD_op_v[2] = &operator_v_fd4;
    FD_op_v[3] = &operator_v_fd6;
    FD_op_v[4] = &operator_v_fd8;
    FD_op_v[5] = &operator_v_fd10;
    FD_op_v[6] = &operator_v_fd12;
    
    if ( ( MYID==0 ) && ( ( nt+ ( OUTNTIMESTEPINFO-1 ) ) %OUTNTIMESTEPINFO ) ==0 ) {
        time1=MPI_Wtime();
        fprintf ( FP,"\n **Message from update_v_PML_4 (printed by PE %d):\n",MYID );
        fprintf ( FP," Updating particle velocities ..." );
    }
    
    /* ------------------------------------------------------------
     * Important!
     * rip and rjp are reciprocal values of averaged densities
     * ------------------------------------------------------------ */
    
    /*-------------------------update-CPML-Boundarys------------------------------------*/
    
    /* interior */
    /*for ( j=gy[2]+1; j<=gy[3]; j++ ) {
     for ( i=gx[2]+1; i<=gx[3]; i++ ) {
     FD_op_v[fdoh] ( i,j,&sxx_x, &sxy_x, &sxy_y,&syy_y, sxx,syy,sxy,hc );
     wavefield_update_v ( i,j,sxx_x,sxy_x,sxy_y,syy_y,vx,vy, rip, rjp );
     
     }
     }*/
    
    /* left boundary */
    
    for ( j=gy[2]+1; j<=gy[3]; j++ ) {
        for ( i=gx[1]; i<=gx[2]; i++ ) {
            
            FD_op_v[fdoh] ( i,j,&sxx_x, &sxy_x, &sxy_y,&syy_y, sxx,syy,sxy,hc );
            
            cpml_update_v_x ( i,j,&sxx_x, &sxy_x,K_x,a_x, b_x,K_x_half,
                             a_x_half,b_x_half ,psi_sxx_x,psi_sxy_x );
            
            wavefield_update_v_4 ( i,j,sxx_x,sxy_x,sxy_y,syy_y,vx,vy, rip, rjp,svx_1,svx_2,svx_3,svx_4,svy_1,svy_2,svy_3,svy_4);
            
        }
    }
    
    /* right boundary */
    
    
    for ( j=gy[2]+1; j<=gy[3]; j++ ) {
        for ( i=gx[3]+1; i<=gx[4]; i++ ) {
            
            
            FD_op_v[fdoh] ( i,j,&sxx_x, &sxy_x, &sxy_y,&syy_y, sxx,syy,sxy,hc );
            
            h1 = ( i-nx2+2*FW );
            cpml_update_v_x ( h1,j,&sxx_x, &sxy_x,K_x,a_x, b_x,K_x_half,
                             a_x_half,b_x_half ,psi_sxx_x,psi_sxy_x );
            
            wavefield_update_v_4 ( i,j,sxx_x,sxy_x,sxy_y,syy_y,vx,vy, rip, rjp,svx_1,svx_2,svx_3,svx_4,svy_1,svy_2,svy_3,svy_4);
            
        }
    }
    
    
    
    /* top boundary */
    
    for ( j=gy[1]; j<=gy[2]; j++ ) {
        for ( i=gx[2]+1; i<=gx[3]; i++ ) {
            
            FD_op_v[fdoh] ( i,j,&sxx_x, &sxy_x, &sxy_y,&syy_y, sxx,syy,sxy,hc );
            
            cpml_update_v_y (i,j,&sxy_y,&syy_y,K_y,a_y,b_y,K_y_half,a_y_half,b_y_half ,
                             psi_syy_y,psi_sxy_y);
            
            wavefield_update_v_4 ( i,j,sxx_x,sxy_x,sxy_y,syy_y,vx,vy, rip, rjp,svx_1,svx_2,svx_3,svx_4,svy_1,svy_2,svy_3,svy_4);
            
            
        }
    }
    
    
    
    /* bottom boundary */
    
    
    for ( j=gy[3]+1; j<=gy[4]; j++ ) {
        for ( i=gx[2]+1; i<=gx[3]; i++ ) {
            
            FD_op_v[fdoh] ( i,j,&sxx_x, &sxy_x, &sxy_y,&syy_y, sxx,syy,sxy,hc );
            
            h1 = ( j-ny2+2*FW );
            
            cpml_update_v_y (i,h1,&sxy_y,&syy_y,K_y,a_y,b_y,K_y_half,a_y_half,b_y_half ,
                             psi_syy_y,psi_sxy_y);
            
            wavefield_update_v_4 ( i,j,sxx_x,sxy_x,sxy_y,syy_y,vx,vy, rip, rjp,svx_1,svx_2,svx_3,svx_4,svy_1,svy_2,svy_3,svy_4);
        }
    }
    
    /* corners */
    
    /*left-top*/
    
    for ( j=gy[1]; j<=gy[2]; j++ ) {
        for ( i=gx[1]; i<=gx[2]; i++ ) {
            
            FD_op_v[fdoh] ( i,j,&sxx_x, &sxy_x, &sxy_y,&syy_y, sxx,syy,sxy,hc );
            
            cpml_update_v_x ( i,j,&sxx_x, &sxy_x,K_x,a_x, b_x,K_x_half,
                             a_x_half,b_x_half ,psi_sxx_x,psi_sxy_x );
            
            cpml_update_v_y (i,j,&sxy_y,&syy_y,K_y,a_y,b_y,K_y_half,a_y_half,b_y_half ,
                             psi_syy_y,psi_sxy_y);
            
            wavefield_update_v_4 ( i,j,sxx_x,sxy_x,sxy_y,syy_y,vx,vy, rip, rjp,svx_1,svx_2,svx_3,svx_4,svy_1,svy_2,svy_3,svy_4);
            
        }
    }
    
    
    /*left-bottom*/
    
    for ( j=gy[3]+1; j<=gy[4]; j++ ) {
        for ( i=gx[1]; i<=gx[2]; i++ ) {
            
            FD_op_v[fdoh] ( i,j,&sxx_x, &sxy_x, &sxy_y,&syy_y, sxx,syy,sxy,hc );
            
            cpml_update_v_x ( i,j,&sxx_x, &sxy_x,K_x,a_x, b_x,K_x_half,
                             a_x_half,b_x_half ,psi_sxx_x,psi_sxy_x );
            
            h1 = ( j-ny2+2*FW );
            
            cpml_update_v_y (i,h1,&sxy_y,&syy_y,K_y,a_y,b_y,K_y_half,a_y_half,b_y_half ,
                             psi_syy_y,psi_sxy_y);
            
            wavefield_update_v_4 ( i,j,sxx_x,sxy_x,sxy_y,syy_y,vx,vy, rip, rjp,svx_1,svx_2,svx_3,svx_4,svy_1,svy_2,svy_3,svy_4);
            
        }
    }
    
    
    
    /* right-top */
    
    
    for ( j=gy[1]; j<=gy[2]; j++ ) {
        for ( i=gx[3]+1; i<=gx[4]; i++ ) {
            
            
            FD_op_v[fdoh] ( i,j,&sxx_x, &sxy_x, &sxy_y,&syy_y, sxx,syy,sxy,hc );
            
            h1 = ( i-nx2+2*FW );
            
            cpml_update_v_x ( h1,j,&sxx_x, &sxy_x,K_x,a_x, b_x,K_x_half,
                             a_x_half,b_x_half ,psi_sxx_x,psi_sxy_x );
            
            cpml_update_v_y (i,j,&sxy_y,&syy_y,K_y,a_y,b_y,K_y_half,a_y_half,b_y_half ,
                             psi_syy_y,psi_sxy_y);
            
            wavefield_update_v_4 ( i,j,sxx_x,sxy_x,sxy_y,syy_y,vx,vy, rip, rjp,svx_1,svx_2,svx_3,svx_4,svy_1,svy_2,svy_3,svy_4);
            
        }
    }
    
    
    /* right-bottom */
    
    
    for ( j=gy[3]+1; j<=gy[4]; j++ ) {
        for ( i=gx[3]+1; i<=gx[4]; i++ ) {
            
            
            FD_op_v[fdoh] ( i,j,&sxx_x, &sxy_x, &sxy_y,&syy_y, sxx,syy,sxy,hc );
            
            h1 = ( i-nx2+2*FW );
            
            cpml_update_v_x ( h1,j,&sxx_x, &sxy_x,K_x,a_x, b_x,K_x_half,
                             a_x_half,b_x_half ,psi_sxx_x,psi_sxy_x );
            
            h1 = ( j-ny2+2*FW );
            
            cpml_update_v_y (i,h1,&sxy_y,&syy_y,K_y,a_y,b_y,K_y_half,a_y_half,b_y_half ,
                             psi_syy_y,psi_sxy_y);
            
            wavefield_update_v_4 ( i,j,sxx_x,sxy_x,sxy_y,syy_y,vx,vy, rip, rjp,svx_1,svx_2,svx_3,svx_4,svy_1,svy_2,svy_3,svy_4);
            
        }
    }
    
    if ( ( MYID==0 ) && ( ( nt+ ( OUTNTIMESTEPINFO-1 ) ) %OUTNTIMESTEPINFO ) ==0 ) {
        time2=MPI_Wtime();
        fprintf ( FP," finished (real time: %4.3f s).\n",time2-time1 );
    }
}