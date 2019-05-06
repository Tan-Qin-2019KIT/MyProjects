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

/* $Id: update_v_abs.c 819 2015-04-17 11:07:06Z tmetz $*/
/*---------------------------------------------------------------------------------------
 *   updating particle velocities at gridpoints of the absorbing frame (ABS=2 in the json file)
 *   by a staggered grid finite difference scheme of FDORDER accuracy in space
 *   and second order accuracy in time
 *   T. Bohlen
 *   
 *   gx and gy are arrays with the locations of the boundary specified in subgrid_bounds.c 
 *   for each subgrid
 *  --------------------------------------------------------------------------------------*/

#include "fd.h"


void update_v_abs ( int nx1, int nx2, int ny1, int ny2, int * gx, int * gy, int nt, 
		    float **  vx, float ** vy, float ** sxx, float ** syy, float ** sxy,  
		    float  **rip, float **rjp, float ** absorb_coeff,float *hc)
{

	int i, j, fdoh;
	float sxx_x, syy_y, sxy_y, sxy_x;

	extern int OUTNTIMESTEPINFO;
	double time1=0.0, time2=0.0;
	extern int MYID, FDORDER; 
	extern FILE *FP;

	fdoh=FDORDER/2;
	
	/* Pointer to function */
	void ( *FD_op_v[7] ) ();
	
	/* array with locations of the fd-op-functions */
	FD_op_v[1] = &operator_v_fd2;
	FD_op_v[2] = &operator_v_fd4;
	FD_op_v[3] = &operator_v_fd6;
	FD_op_v[4] = &operator_v_fd8;
	FD_op_v[5] = &operator_v_fd10;
	FD_op_v[6] = &operator_v_fd12;

	if ( ( MYID==0 ) && ( ( nt+ ( OUTNTIMESTEPINFO-1 ) ) %OUTNTIMESTEPINFO ) ==0 ) {
		time1=MPI_Wtime();
		fprintf ( FP,"\n **Message from update_v_abs2 (printed by PE %d):\n",MYID );
		fprintf ( FP," Updating particle velocities ..." );
	}

	/* ------------------------------------------------------------
	 * Important!
	 * rip and rjp are reciprocal values of averaged densities
	 * ------------------------------------------------------------ */
	


/* left boundary */

	for ( j=gy[2]+1; j<=gy[3]; j++ ) {
		for ( i=gx[1]; i<=gx[2]; i++ ) {
			
			
			FD_op_v[fdoh] ( i,j,&sxx_x, &sxy_x, &sxy_y,&syy_y, sxx,syy,sxy,hc );

			wavefield_update_v ( i,j,sxx_x,sxy_x,sxy_y,syy_y,vx,vy, rip, rjp );
			
			/* Damping the wavfield */
			abs_update_v (i, j, vx,vy, absorb_coeff);

		}
	}

	/* right boundary */


	for ( j=gy[2]+1; j<=gy[3]; j++ ) {
		for ( i=gx[3]+1; i<=gx[4]; i++ ) {


			FD_op_v[fdoh] ( i,j,&sxx_x, &sxy_x, &sxy_y,&syy_y, sxx,syy,sxy,hc );
			
			wavefield_update_v ( i,j,sxx_x,sxy_x,sxy_y,syy_y,vx,vy, rip, rjp );
			
			abs_update_v (i, j, vx,vy, absorb_coeff);

		}
	}



	/* top boundary */

	for ( j=gy[1]; j<=gy[2]; j++ ) {
		for ( i=gx[2]+1; i<=gx[3]; i++ ) {

			FD_op_v[fdoh] ( i,j,&sxx_x, &sxy_x, &sxy_y,&syy_y, sxx,syy,sxy,hc );

			wavefield_update_v ( i,j,sxx_x,sxy_x,sxy_y,syy_y,vx,vy, rip, rjp );
			
			abs_update_v (i, j, vx,vy, absorb_coeff);
		

		}
	}



	/* bottom boundary */


	for ( j=gy[3]+1; j<=gy[4]; j++ ) {
		for ( i=gx[2]+1; i<=gx[3]; i++ ) {

			FD_op_v[fdoh] ( i,j,&sxx_x, &sxy_x, &sxy_y,&syy_y, sxx,syy,sxy,hc );

			wavefield_update_v ( i,j,sxx_x,sxy_x,sxy_y,syy_y,vx,vy, rip, rjp );
			
			abs_update_v (i, j, vx,vy, absorb_coeff);
		}
	}

	/* corners */

	/*left-top*/

	for ( j=gy[1]; j<=gy[2]; j++ ) {
		for ( i=gx[1]; i<=gx[2]; i++ ) {

			FD_op_v[fdoh] ( i,j,&sxx_x, &sxy_x, &sxy_y,&syy_y, sxx,syy,sxy,hc );

			wavefield_update_v ( i,j,sxx_x,sxy_x,sxy_y,syy_y,vx,vy, rip, rjp );
			
			abs_update_v (i, j, vx,vy, absorb_coeff);

		}
	}


	/*left-bottom*/

	for ( j=gy[3]+1; j<=gy[4]; j++ ) {
		for ( i=gx[1]; i<=gx[2]; i++ ) {

			FD_op_v[fdoh] ( i,j,&sxx_x, &sxy_x, &sxy_y,&syy_y, sxx,syy,sxy,hc );

			wavefield_update_v ( i,j,sxx_x,sxy_x,sxy_y,syy_y,vx,vy, rip, rjp );
			
			abs_update_v (i, j, vx,vy, absorb_coeff);

		}
	}



	/* right-top */


	for ( j=gy[1]; j<=gy[2]; j++ ) {
		for ( i=gx[3]+1; i<=gx[4]; i++ ) {

			FD_op_v[fdoh] ( i,j,&sxx_x, &sxy_x, &sxy_y,&syy_y, sxx,syy,sxy,hc );

			wavefield_update_v ( i,j,sxx_x,sxy_x,sxy_y,syy_y,vx,vy, rip, rjp );
			
			abs_update_v (i, j, vx,vy, absorb_coeff);

		}
	}


	/* right-bottom */


	for ( j=gy[3]+1; j<=gy[4]; j++ ) {
		for ( i=gx[3]+1; i<=gx[4]; i++ ) {


			FD_op_v[fdoh] ( i,j,&sxx_x, &sxy_x, &sxy_y,&syy_y, sxx,syy,sxy,hc );

			wavefield_update_v ( i,j,sxx_x,sxy_x,sxy_y,syy_y,vx,vy, rip, rjp );
			
			abs_update_v (i, j, vx,vy, absorb_coeff);

		}
	}
	
	if ( ( MYID==0 ) && ( ( nt+ ( OUTNTIMESTEPINFO-1 ) ) %OUTNTIMESTEPINFO ) ==0 ) {
		time2=MPI_Wtime();
		fprintf ( FP," finished (real time: %4.3f s).\n",time2-time1 );
	}
}
