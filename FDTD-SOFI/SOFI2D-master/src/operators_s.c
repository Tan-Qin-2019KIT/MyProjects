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

/* $Id: operators_s.c 819 2015-04-17 11:07:06Z tmetz $ */

/*------------------------------------------------------------------------
 * FD Operators (Order 2-12) for the particle velocity components (used in 
 * the update functions of the boundary frame).
 *    
 * ----------------------------------------------------------------------*/

#include "fd.h"

/*FD Operators Order 2-12  ----------------------------------------------------------------------------------------------------*/
void operator_s_fd2 ( int i, int j,float  * vxx, float * vyx,float * vxy,float * vyy, float **vx, float **vy,float * hc )
{
	float  dhi;
	extern float DH;

	dhi = 1.0/DH;

	*vxx = hc[1]* ( vx[j][i]  -vx[j][i-1] ) *dhi;
	*vyx = hc[1]* ( vy[j][i+1]-vy[j][i] ) *dhi;
	*vxy = hc[1]* ( vx[j+1][i]-vx[j][i] ) *dhi;
	*vyy = hc[1]* ( vy[j][i]  -vy[j-1][i] ) *dhi;
}

void operator_s_fd4 ( int i, int j,float  * vxx, float * vyx,float * vxy,float * vyy, float **vx, float **vy,float * hc )
{
	float  dhi;
	extern float DH;

	dhi = 1.0/DH;

	*vxx = ( hc[1]* ( vx[j][i]  -vx[j][i-1] )
	         + hc[2]* ( vx[j][i+1]-vx[j][i-2] ) ) *dhi;

	*vyx = ( hc[1]* ( vy[j][i+1]-vy[j][i] )
	         + hc[2]* ( vy[j][i+2]-vy[j][i-1] ) ) *dhi;

	*vxy = ( hc[1]* ( vx[j+1][i]-vx[j][i] )
	         + hc[2]* ( vx[j+2][i]-vx[j-1][i] ) ) *dhi;

	*vyy = ( hc[1]* ( vy[j][i]  -vy[j-1][i] )
	         + hc[2]* ( vy[j+1][i]-vy[j-2][i] ) ) *dhi;
}

void operator_s_fd6 ( int i, int j,float  * vxx, float * vyx,float * vxy,float * vyy, float **vx, float **vy,float * hc )
{
	float  dhi;
	extern float DH;

	dhi = 1.0/DH;
	*vxx = ( hc[1]* ( vx[j][i]  -vx[j][i-1] )
	         + hc[2]* ( vx[j][i+1]-vx[j][i-2] )
	         + hc[3]* ( vx[j][i+2]-vx[j][i-3] ) ) *dhi;

	*vyx = ( hc[1]* ( vy[j][i+1]-vy[j][i] )
	         + hc[2]* ( vy[j][i+2]-vy[j][i-1] )
	         + hc[3]* ( vy[j][i+3]-vy[j][i-2] ) ) *dhi;

	*vxy = ( hc[1]* ( vx[j+1][i]-vx[j][i] )
	         + hc[2]* ( vx[j+2][i]-vx[j-1][i] )
	         + hc[3]* ( vx[j+3][i]-vx[j-2][i] ) ) *dhi;

	*vyy = ( hc[1]* ( vy[j][i]  -vy[j-1][i] )
	         + hc[2]* ( vy[j+1][i]-vy[j-2][i] )
	         + hc[3]* ( vy[j+2][i]-vy[j-3][i] ) ) *dhi;
}

void operator_s_fd8 ( int i, int j,float  * vxx, float * vyx,float * vxy,float * vyy, float **vx, float **vy,float * hc )
{
	float  dhi;
	extern float DH;

	dhi = 1.0/DH;

	*vxx = ( hc[1]* ( vx[j][i]  -vx[j][i-1] )
	         + hc[2]* ( vx[j][i+1]-vx[j][i-2] )
	         + hc[3]* ( vx[j][i+2]-vx[j][i-3] )
	         + hc[4]* ( vx[j][i+3]-vx[j][i-4] ) ) *dhi;

	*vyx = ( hc[1]* ( vy[j][i+1]-vy[j][i] )
	         + hc[2]* ( vy[j][i+2]-vy[j][i-1] )
	         + hc[3]* ( vy[j][i+3]-vy[j][i-2] )
	         + hc[4]* ( vy[j][i+4]-vy[j][i-3] ) ) *dhi;

	*vxy = ( hc[1]* ( vx[j+1][i]-vx[j][i] )
	         + hc[2]* ( vx[j+2][i]-vx[j-1][i] )
	         + hc[3]* ( vx[j+3][i]-vx[j-2][i] )
	         + hc[4]* ( vx[j+4][i]-vx[j-3][i] ) ) *dhi;

	*vyy = ( hc[1]* ( vy[j][i]  -vy[j-1][i] )
	         + hc[2]* ( vy[j+1][i]-vy[j-2][i] )
	         + hc[3]* ( vy[j+2][i]-vy[j-3][i] )
	         + hc[4]* ( vy[j+3][i]-vy[j-4][i] ) ) *dhi;
}

void operator_s_fd10 ( int i, int j,float  * vxx, float * vyx,float * vxy,float * vyy, float **vx, float **vy,float * hc )
{
	float  dhi;
	extern float DH;

	dhi = 1.0/DH;
	*vxx = ( hc[1]* ( vx[j][i]  -vx[j][i-1] )
	         + hc[2]* ( vx[j][i+1]-vx[j][i-2] )
	         + hc[3]* ( vx[j][i+2]-vx[j][i-3] )
	         + hc[4]* ( vx[j][i+3]-vx[j][i-4] )
	         + hc[5]* ( vx[j][i+4]-vx[j][i-5] )
	       ) *dhi;
	*vyy = ( hc[1]* ( vy[j][i]  -vy[j-1][i] )
	         + hc[2]* ( vy[j+1][i]-vy[j-2][i] )
	         + hc[3]* ( vy[j+2][i]-vy[j-3][i] )
	         + hc[4]* ( vy[j+3][i]-vy[j-4][i] )
	         + hc[5]* ( vy[j+4][i]-vy[j-5][i] )
	       ) *dhi;
	*vyx = ( hc[1]* ( vy[j][i+1]-vy[j][i] )
	         + hc[2]* ( vy[j][i+2]-vy[j][i-1] )
	         + hc[3]* ( vy[j][i+3]-vy[j][i-2] )
	         + hc[4]* ( vy[j][i+4]-vy[j][i-3] )
	         + hc[5]* ( vy[j][i+5]-vy[j][i-4] )
	       ) *dhi;
	*vxy = ( hc[1]* ( vx[j+1][i]-vx[j][i] )
	         + hc[2]* ( vx[j+2][i]-vx[j-1][i] )
	         + hc[3]* ( vx[j+3][i]-vx[j-2][i] )
	         + hc[4]* ( vx[j+4][i]-vx[j-3][i] )
	         + hc[5]* ( vx[j+5][i]-vx[j-4][i] )
	       ) *dhi;
}

void operator_s_fd12 ( int i, int j,float  * vxx, float * vyx,float * vxy,float * vyy, float **vx, float **vy,float * hc )
{
	float  dhi;
	extern float DH;

	dhi = 1.0/DH;
	*vxx = ( hc[1]* ( vx[j][i]  -vx[j][i-1] )
	         + hc[2]* ( vx[j][i+1]-vx[j][i-2] )
	         + hc[3]* ( vx[j][i+2]-vx[j][i-3] )
	         + hc[4]* ( vx[j][i+3]-vx[j][i-4] )
	         + hc[5]* ( vx[j][i+4]-vx[j][i-5] )
	         + hc[6]* ( vx[j][i+5]-vx[j][i-6] )
	       ) *dhi;
	*vyy = ( hc[1]* ( vy[j][i]  -vy[j-1][i] )
	         + hc[2]* ( vy[j+1][i]-vy[j-2][i] )
	         + hc[3]* ( vy[j+2][i]-vy[j-3][i] )
	         + hc[4]* ( vy[j+3][i]-vy[j-4][i] )
	         + hc[5]* ( vy[j+4][i]-vy[j-5][i] )
	         + hc[6]* ( vy[j+5][i]-vy[j-6][i] )
	       ) *dhi;
	*vyx = ( hc[1]* ( vy[j][i+1]-vy[j][i] )
	         + hc[2]* ( vy[j][i+2]-vy[j][i-1] )
	         + hc[3]* ( vy[j][i+3]-vy[j][i-2] )
	         + hc[4]* ( vy[j][i+4]-vy[j][i-3] )
	         + hc[5]* ( vy[j][i+5]-vy[j][i-4] )
	         + hc[6]* ( vy[j][i+6]-vy[j][i-5] )
	       ) *dhi;
	*vxy = ( hc[1]* ( vx[j+1][i]-vx[j][i] )
	         + hc[2]* ( vx[j+2][i]-vx[j-1][i] )
	         + hc[3]* ( vx[j+3][i]-vx[j-2][i] )
	         + hc[4]* ( vx[j+4][i]-vx[j-3][i] )
	         + hc[5]* ( vx[j+5][i]-vx[j-4][i] )
	         + hc[6]* ( vx[j+6][i]-vx[j-5][i] )
	       ) *dhi;
}