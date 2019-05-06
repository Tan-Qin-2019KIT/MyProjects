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

/* $Id: update_s_elastic_PML.c 819 2015-04-17 11:07:06Z tmetz $*/
/*------------------------------------------------------------------------
 *   updating stress components at gridpoints of the CPML-frame (ABS=1 in the json file)
 *   by a staggered grid finite difference scheme of arbitrary (FDORDER) order accuracy in space
 *   and second order accuracy in time
 *   T. Bohlen
 *
 *   gx and gy are arrays with the locations of the boundary specified in subgrid_bounds.c
 *   for each subgrid
 *  ----------------------------------------------------------------------*/

#include "fd.h"


void update_s_elastic_PML ( int nx1, int nx2, int ny1, int ny2, int * gx, int * gy, int nt,
                            float **  vx, float **   vy, float **   sxx, float **   syy,
                            float **   sxy, float ** pi, float ** u, float ** uipjp, float *hc,
                            float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
                            float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,
                            float ** psi_vxx, float ** psi_vyy, float ** psi_vxy, float ** psi_vyx )
{
	int i,j, h1,fdoh;
	float  vxx, vyy, vxy, vyx;
	extern int MYID, FDORDER, FW;
	extern FILE *FP;
	extern int OUTNTIMESTEPINFO;
	double time1=0.0, time2=0.0;

	fdoh=FDORDER/2;
	
	/*Pointer array to the locations of the fd-operator functions*/
	void ( *FD_op_s[7] ) ();
	FD_op_s[1] = &operator_s_fd2;
	FD_op_s[2] = &operator_s_fd4;
	FD_op_s[3] = &operator_s_fd6;
	FD_op_s[4] = &operator_s_fd8;
	FD_op_s[5] = &operator_s_fd10;
	FD_op_s[6] = &operator_s_fd12;

	if ( ( MYID==0 ) && ( ( nt+ ( OUTNTIMESTEPINFO-1 ) ) %OUTNTIMESTEPINFO ) ==0 ) {
		time1=MPI_Wtime();
		fprintf ( FP,"\n **Message from update_s_PML (printed by PE %d):\n",MYID );
		fprintf ( FP," Updating stress components ..." );
	}


	/* interior */
	/*for ( j=gy[2]+1; j<=gy[3]; j++ ) {
		for ( i=gx[2]+1; i<=gx[3]; i++ ) {
			FD_op_s[fdoh] ( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy,hc );

			wavefield_update_s_el ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,pi,u,uipjp );
		}
	}*/
	/* left boundary */
	for ( j=gy[2]+1; j<=gy[3]; j++ ) {
		for ( i=gx[1]; i<=gx[2]; i++ ) {
			FD_op_s[fdoh] ( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy,hc );

			cpml_update_s_x ( i,j,&vxx,&vyx,K_x,a_x,
			               b_x, K_x_half, a_x_half, b_x_half ,psi_vxx,psi_vyx );

			wavefield_update_s_el ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,pi,u,uipjp );
		}
	}

	/* right boundary */
	for ( j=gy[2]+1; j<=gy[3]; j++ ) {
		for ( i=gx[3]+1; i<=gx[4]; i++ ) {
			FD_op_s[fdoh] ( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy,hc );
			h1 = ( i-nx2+2*FW );

			cpml_update_s_x ( h1,j,&vxx,&vyx,K_x,a_x,
			               b_x, K_x_half, a_x_half, b_x_half ,psi_vxx,psi_vyx );

			wavefield_update_s_el ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,pi,u,uipjp );

		}
	}

	/* top boundary */
	for ( j=gy[1]; j<=gy[2]; j++ ) {
		for ( i=gx[2]+1; i<=gx[3]; i++ ) {
			FD_op_s[fdoh] ( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy,hc );

			cpml_update_s_y ( i,j,&vxy,&vyy,K_y,a_y,
			               b_y, K_y_half, a_y_half, b_y_half ,psi_vyy,psi_vxy );

			wavefield_update_s_el ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,pi,u,uipjp );
		}
	}

	/* bottom boundary */
	for ( j=gy[3]+1; j<=gy[4]; j++ ) {
		for ( i=gx[2]+1; i<=gx[3]; i++ ) {
			FD_op_s[fdoh] ( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy,hc );

			h1 = ( j-ny2+2*FW );
			cpml_update_s_y ( i,h1,&vxy,&vyy,K_y,a_y,
			               b_y, K_y_half, a_y_half, b_y_half ,psi_vyy,psi_vxy );

			wavefield_update_s_el ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,pi,u,uipjp );
		}
	}

	/* corners */

	/*left-top*/
	for ( j=gy[1]; j<=gy[2]; j++ ) {
		for ( i=gx[1]; i<=gx[2]; i++ ) {
			FD_op_s[fdoh] ( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy,hc );

			cpml_update_s_x ( i,j,&vxx,&vyx,K_x,a_x,
			               b_x, K_x_half, a_x_half, b_x_half ,psi_vxx,psi_vyx );

			cpml_update_s_y ( i,j,&vxy,&vyy,K_y,a_y,
			               b_y, K_y_half, a_y_half, b_y_half ,psi_vyy,psi_vxy );

			wavefield_update_s_el ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,pi,u,uipjp );

		}
	}

	/*left-bottom*/
	for ( j=gy[3]+1; j<=gy[4]; j++ ) {
		for ( i=gx[1]; i<=gx[2]; i++ ) {
			FD_op_s[fdoh] ( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy,hc );

			cpml_update_s_x ( i,j,&vxx,&vyx,K_x,a_x,
			               b_x, K_x_half, a_x_half, b_x_half ,psi_vxx,psi_vyx );

			h1 = ( j-ny2+2*FW );

			cpml_update_s_y ( i,h1,&vxy,&vyy,K_y,a_y,
			               b_y, K_y_half, a_y_half, b_y_half ,psi_vyy,psi_vxy );

			wavefield_update_s_el ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,pi,u,uipjp );
		}
	}

	/* right-top */
	for ( j=gy[1]; j<=gy[2]; j++ ) {
		for ( i=gx[3]+1; i<=gx[4]; i++ ) {
			FD_op_s[fdoh] ( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy,hc );

			h1 = ( i-nx2+2*FW );
			cpml_update_s_x ( h1,j,&vxx,&vyx,K_x,a_x,
			               b_x, K_x_half, a_x_half, b_x_half ,psi_vxx,psi_vyx );

			cpml_update_s_y ( i,j,&vxy,&vyy,K_y,a_y,
			               b_y, K_y_half, a_y_half, b_y_half ,psi_vyy,psi_vxy );

			wavefield_update_s_el ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,pi,u,uipjp );
		}
	}

	/* right-bottom */
	for ( j=gy[3]+1; j<=gy[4]; j++ ) {
		for ( i=gx[3]+1; i<=gx[4]; i++ ) {
			FD_op_s[fdoh] ( i,j,&vxx,&vyx,&vxy,&vyy,vx,vy,hc );

			h1 = ( i-nx2+2*FW );
			cpml_update_s_x ( h1,j,&vxx,&vyx,K_x,a_x,
			               b_x, K_x_half, a_x_half, b_x_half ,psi_vxx,psi_vyx );

			h1 = ( j-ny2+2*FW );
			cpml_update_s_y ( i,h1,&vxy,&vyy,K_y,a_y,
			               b_y, K_y_half, a_y_half, b_y_half ,psi_vyy,psi_vxy );

			wavefield_update_s_el ( i,j,vxx,vyx,vxy,vyy,sxy,sxx,syy,pi,u,uipjp );
		}
	}


	if ( ( MYID==0 ) && ( ( nt+ ( OUTNTIMESTEPINFO-1 ) ) %OUTNTIMESTEPINFO ) ==0 ) {
		time2=MPI_Wtime();
		fprintf ( FP," finished (real time: %4.3f s).\n",time2-time1 );
	}
}
