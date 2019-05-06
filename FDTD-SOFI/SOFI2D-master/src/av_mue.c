/*------------------------------------------------------------------------
 * Copyright (C) 2011 For the list of authors, see file AUTHORS.
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
/* -------------------------------------------------------------
 * Averaging of material parameters (mue)
 * -------------------------------------------------------------*/

#include "fd.h"

void av_mue(float ** u, float ** uipjp){

	extern int NX, NY;
	int i, j;
	for (j=1;j<=NY;j++){
		for (i=1;i<=NX;i++){
	       
		       uipjp[j][i]=4.0/((1.0/u[j][i])+
			  	(1.0/u[j][i+1])+(1.0/u[j+1][i])+(1.0/u[j+1][i+1])); 
	
		}
	}
}