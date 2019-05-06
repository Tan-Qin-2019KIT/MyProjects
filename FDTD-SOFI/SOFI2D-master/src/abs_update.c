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

/* $Id: abs_update.c 819 2015-04-17 11:07:06Z tmetz $ */

/*------------------------------------------------------------------------
 * Damping functions for the stress fields and the particle velocity fields
 * used at the absorbing boundaries (ABS=2).
 *    
 * ----------------------------------------------------------------------*/


#include "fd.h"

void abs_update_s ( int i, int j, float **sxx,float **sxy,float **syy, float ** absorb_coeff ){
                               sxy[j][i]*=absorb_coeff[j][i];
				sxx[j][i]*=absorb_coeff[j][i];
				syy[j][i]*=absorb_coeff[j][i]; 
}

void abs_update_v ( int i, int j, float **vx,float **vy, float ** absorb_coeff ){
                                vx[j][i]*=absorb_coeff[j][i];
				vy[j][i]*=absorb_coeff[j][i];
}