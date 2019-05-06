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

/* $Id: operators_v.c 819 2015-04-17 11:07:06Z tmetz $ */

/*------------------------------------------------------------------------
 * FD Operators (Order 2-12) for the stress components (used in 
 * the update functions of the boundary frame).
 *    
 * ----------------------------------------------------------------------*/

#include "fd.h"


/*FD Operators Order 2-12  ----------------------------------------------------------------------------------------------------*/
void operator_v_fd2 ( int i, int j,float  * sxx_x, float * sxy_x,float * sxy_y,float * syy_y, float **sxx, float **syy,float **sxy, float * hc )
{

	*sxx_x =  hc[1]* ( sxx[j][i+1]-sxx[j][i] );
	*sxy_x =  hc[1]* ( sxy[j][i]-sxy[j][i-1] );
	*sxy_y =  hc[1]* ( sxy[j][i]-sxy[j-1][i] );
	*syy_y =  hc[1]* ( syy[j+1][i]-syy[j][i] );
}

void operator_v_fd4 ( int i, int j,float  * sxx_x, float * sxy_x,float * sxy_y,float * syy_y, float **sxx, float **syy,float **sxy, float * hc )
{
	*sxx_x =  hc[1]* ( sxx[j][i+1]-sxx[j][i] )
	          + hc[2]* ( sxx[j][i+2]-sxx[j][i-1] );

	*sxy_x = hc[1]* ( sxy[j][i]-sxy[j][i-1] )
	         + hc[2]* ( sxy[j][i+1]-sxy[j][i-2] );

	*sxy_y =  hc[1]* ( sxy[j][i]-sxy[j-1][i] )
	          + hc[2]* ( sxy[j+1][i]-sxy[j-2][i] );

	*syy_y = hc[1]* ( syy[j+1][i]-syy[j][i] )
	         + hc[2]* ( syy[j+2][i]-syy[j-1][i] );

}

void operator_v_fd6 ( int i, int j,float  * sxx_x, float * sxy_x,float * sxy_y,float * syy_y, float **sxx, float **syy,float **sxy, float * hc )
{
	*sxx_x =  hc[1]* ( sxx[j][i+1]-sxx[j][i] )
	          + hc[2]* ( sxx[j][i+2]-sxx[j][i-1] )
	          + hc[3]* ( sxx[j][i+3]-sxx[j][i-2] );

	*sxy_x = hc[1]* ( sxy[j][i]-sxy[j][i-1] )
	         + hc[2]* ( sxy[j][i+1]-sxy[j][i-2] )
	         + hc[3]* ( sxy[j][i+2]-sxy[j][i-3] );

	*sxy_y =  hc[1]* ( sxy[j][i]-sxy[j-1][i] )
	          + hc[2]* ( sxy[j+1][i]-sxy[j-2][i] )
	          + hc[3]* ( sxy[j+2][i]-sxy[j-3][i] );

	*syy_y = hc[1]* ( syy[j+1][i]-syy[j][i] )
	         + hc[2]* ( syy[j+2][i]-syy[j-1][i] )
	         + hc[3]* ( syy[j+3][i]-syy[j-2][i] );

}

void operator_v_fd8 ( int i, int j,float  * sxx_x, float * sxy_x,float * sxy_y,float * syy_y, float **sxx, float **syy,float **sxy, float * hc )
{
	*sxx_x =  hc[1]* ( sxx[j][i+1]-sxx[j][i] )
	          + hc[2]* ( sxx[j][i+2]-sxx[j][i-1] )
	          + hc[3]* ( sxx[j][i+3]-sxx[j][i-2] )
	          + hc[4]* ( sxx[j][i+4]-sxx[j][i-3] );

	*sxy_x = hc[1]* ( sxy[j][i]-sxy[j][i-1] )
	         + hc[2]* ( sxy[j][i+1]-sxy[j][i-2] )
	         + hc[3]* ( sxy[j][i+2]-sxy[j][i-3] )
	         + hc[4]* ( sxy[j][i+3]-sxy[j][i-4] );

	*sxy_y =  hc[1]* ( sxy[j][i]-sxy[j-1][i] )
	          + hc[2]* ( sxy[j+1][i]-sxy[j-2][i] )
	          + hc[3]* ( sxy[j+2][i]-sxy[j-3][i] )
	          + hc[4]* ( sxy[j+3][i]-sxy[j-4][i] );

	*syy_y = hc[1]* ( syy[j+1][i]-syy[j][i] )
	         + hc[2]* ( syy[j+2][i]-syy[j-1][i] )
	         + hc[3]* ( syy[j+3][i]-syy[j-2][i] )
	         + hc[4]* ( syy[j+4][i]-syy[j-3][i] );


}

void operator_v_fd10 ( int i, int j,float  * sxx_x, float * sxy_x,float * sxy_y,float * syy_y, float **sxx, float **syy,float **sxy, float * hc )
{
	*sxx_x= hc[1]* ( sxx[j][i+1]-sxx[j][i] )
	        + hc[2]* ( sxx[j][i+2]-sxx[j][i-1] )
	        + hc[3]* ( sxx[j][i+3]-sxx[j][i-2] )
	        + hc[4]* ( sxx[j][i+4]-sxx[j][i-3] )
	        + hc[5]* ( sxx[j][i+5]-sxx[j][i-4] );

	*sxy_x= hc[1]* ( sxy[j][i]  -sxy[j][i-1] )
	        + hc[2]* ( sxy[j][i+1]-sxy[j][i-2] )
	        + hc[3]* ( sxy[j][i+2]-sxy[j][i-3] )
	        + hc[4]* ( sxy[j][i+3]-sxy[j][i-4] )
	        + hc[5]* ( sxy[j][i+4]-sxy[j][i-5] );


	*sxy_y= hc[1]* ( sxy[j][i]  -sxy[j-1][i] )
	        + hc[2]* ( sxy[j+1][i]-sxy[j-2][i] )
	        + hc[3]* ( sxy[j+2][i]-sxy[j-3][i] )
	        + hc[4]* ( sxy[j+3][i]-sxy[j-4][i] )
	        + hc[5]* ( sxy[j+4][i]-sxy[j-5][i] );


	*syy_y= hc[1]* ( syy[j+1][i]-syy[j][i] )
	        + hc[2]* ( syy[j+2][i]-syy[j-1][i] )
	        + hc[3]* ( syy[j+3][i]-syy[j-2][i] )
	        + hc[4]* ( syy[j+4][i]-syy[j-3][i] )
	        + hc[5]* ( syy[j+5][i]-syy[j-4][i] );

}

void operator_v_fd12 ( int i, int j,float  * sxx_x, float * sxy_x,float * sxy_y,float * syy_y, float **sxx, float **syy,float **sxy, float * hc )
{
	*sxx_x = hc[1]* ( sxx[j][i+1]-sxx[j][i] )
	         + hc[2]* ( sxx[j][i+2]-sxx[j][i-1] )
	         + hc[3]* ( sxx[j][i+3]-sxx[j][i-2] )
	         + hc[4]* ( sxx[j][i+4]-sxx[j][i-3] )
	         + hc[5]* ( sxx[j][i+5]-sxx[j][i-4] )
	         + hc[6]* ( sxx[j][i+6]-sxx[j][i-5] );

	*sxy_x = hc[1]* ( sxy[j][i]  -sxy[j][i-1] )
	         + hc[2]* ( sxy[j][i+1]-sxy[j][i-2] )
	         + hc[3]* ( sxy[j][i+2]-sxy[j][i-3] )
	         + hc[4]* ( sxy[j][i+3]-sxy[j][i-4] )
	         + hc[5]* ( sxy[j][i+4]-sxy[j][i-5] )
	         + hc[6]* ( sxy[j][i+5]-sxy[j][i-6] );

	*sxy_y = hc[1]* ( sxy[j][i]  -sxy[j-1][i] )
	         + hc[2]* ( sxy[j+1][i]-sxy[j-2][i] )
	         + hc[3]* ( sxy[j+2][i]-sxy[j-3][i] )
	         + hc[4]* ( sxy[j+3][i]-sxy[j-4][i] )
	         + hc[5]* ( sxy[j+4][i]-sxy[j-5][i] )
	         + hc[6]* ( sxy[j+5][i]-sxy[j-6][i] );


	*syy_y = hc[1]* ( syy[j+1][i]-syy[j][i] )
	         + hc[2]* ( syy[j+2][i]-syy[j-1][i] )
	         + hc[3]* ( syy[j+3][i]-syy[j-2][i] )
	         + hc[4]* ( syy[j+4][i]-syy[j-3][i] )
	         + hc[5]* ( syy[j+5][i]-syy[j-4][i] )
	         + hc[6]* ( syy[j+6][i]-syy[j-5][i] );

}



