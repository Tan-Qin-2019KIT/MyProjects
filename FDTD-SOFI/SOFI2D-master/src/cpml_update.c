/*---------------------------------------------------------------------------------
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
---------------------------------------------------------------------------------*/

/* $Id: cpml_update.c 819 2015-04-17 11:07:06Z tmetz $ */

/*--------------------------------------------------------------------------------
 * CPML update functions (ABS=1) for the stress and particle velocitiy components used
 * in the update functions of the CPML-boundary area.
 *    
 * -------------------------------------------------------------------------------*/


#include "fd.h"


/* CPML Functions for update_s ---------------------------------------------------*/

void cpml_update_s_x ( int i, int j,float  * vxx, float * vyx,float * K_x, 
		       float * a_x, float * b_x, float * K_x_half, float * a_x_half, 
		       float * b_x_half ,float ** psi_vxx,float ** psi_vyx )
{

	psi_vxx[j][i] = b_x[i] * psi_vxx[j][i] + a_x[i] * ( *vxx );
	*vxx = ( *vxx ) / K_x[i] + psi_vxx[j][i];

	psi_vyx[j][i] = b_x_half[i] * psi_vyx[j][i] + a_x_half[i] * ( *vyx );
	*vyx = ( *vyx ) / K_x_half[i] + psi_vyx[j][i];
}

void cpml_update_s_y ( int i, int j,float * vxy,float * vyy,float * K_y, float * a_y,
                    float * b_y, float * K_y_half, float * a_y_half, float * b_y_half ,
		    float ** psi_vyy,float ** psi_vxy )
{

	psi_vyy[j][i] = b_y[j] * psi_vyy[j][i] + a_y[j] * ( *vyy );
	*vyy = ( *vyy ) / K_y[j] + psi_vyy[j][i];
	
	psi_vxy[j][i] = b_y_half[j] * psi_vxy[j][i] + a_y_half[j] * ( *vxy );
        *vxy = ( *vxy ) / K_y_half[j] + psi_vxy[j][i];
}


/* CPML Functions for update_v ---------------------------------------------------*/

void cpml_update_v_x ( int i, int j,float  * sxx_x, float * sxy_x,float * K_x, 
		       float * a_x,float * b_x, float * K_x_half, float * a_x_half, 
		       float * b_x_half ,float ** psi_sxx_x,float ** psi_sxy_x )
{

	psi_sxx_x[j][i] = b_x_half[i] * psi_sxx_x[j][i] + a_x_half[i] * ( *sxx_x );
	*sxx_x = ( *sxx_x ) / K_x_half[i] + psi_sxx_x[j][i];

	psi_sxy_x[j][i] = b_x[i] * psi_sxy_x[j][i] + a_x[i] * ( *sxy_x );
	*sxy_x = ( *sxy_x ) / K_x[i] + psi_sxy_x[j][i];
}

void cpml_update_v_y ( int i, int j,float * sxy_y,float * syy_y,float * K_y, float * a_y,
                       float * b_y, float * K_y_half, float * a_y_half, float * b_y_half ,
		       float ** psi_syy_y,float ** psi_sxy_y )
{

	psi_syy_y[j][i] = b_y_half[j] * psi_syy_y[j][i] + a_y_half[j] * ( *syy_y );
	*syy_y = ( *syy_y ) / K_y_half[j] + psi_syy_y[j][i];

	psi_sxy_y[j][i] = b_y[j] * psi_sxy_y[j][i] + a_y[j] * ( *sxy_y );
	*sxy_y = ( *sxy_y ) / K_y[j] + psi_sxy_y[j][i];
}

