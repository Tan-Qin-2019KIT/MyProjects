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

/* $Id: wavefield_update_s_el.c 819 2015-04-17 11:07:06Z tmetz $ */

/*Update Function of the stress-Wavefields in the elastic case*/

#include "fd.h"

void wavefield_update_s_el ( int i, int j,float   vxx, float  vyx,float vxy,float  vyy, float **sxy,
                   float **sxx, float ** syy, float ** pi, float ** u, float ** uipjp )
{

	extern float DT;
	float fipjp, f, g;


	fipjp=uipjp[j][i]*DT;
	f = u[j][i]*DT;
	g = pi[j][i]*DT;

	/*Update  */
	sxy[j][i] += fipjp* ( vyx+vxy );
	sxx[j][i] += ( g* ( vxx+vyy ) )- ( 2.0*f*vyy );
	syy[j][i] += ( g* ( vxx+vyy ) )- ( 2.0*f*vxx );

}