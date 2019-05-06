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

/* $Id: wavefield_update_s_visc.c 819 2015-04-17 11:07:06Z tmetz $ */

/*Update Function of the stress-Wavefields in the viscoelastic case*/

#include "fd.h"

void wavefield_update_s_visc ( int i, int j,float   vxx, float  vyx,float vxy,float  vyy, float **sxy,
                               float **sxx, float ** syy, float ***r, float ***p,
                               float ***q,float **fipjp, float **f, float **g, float *bip, 
			       float *bjm,float *cip, float *cjm, float ***d, float ***e, float ***dip )
{
	int l;
	float  dthalbe;
	extern float DT;
	extern int L;
	float sumr=0.0, sump=0.0, sumq=0.0;
	/* computing sums of the old memory variables */
	
	dthalbe = DT/2.0;
	
	sumr=sump=sumq=0.0;
	for ( l=1; l<=L; l++ ) {
		sumr+=r[j][i][l];
		sump+=p[j][i][l];
		sumq+=q[j][i][l];
	}


	/* updating components of the stress tensor, partially */
	sxy[j][i] += ( fipjp[j][i]* ( vxy+vyx ) ) + ( dthalbe*sumr );
	sxx[j][i] += ( g[j][i]* ( vxx+vyy ) )- ( 2.0*f[j][i]*vyy ) + ( dthalbe*sump );
	syy[j][i] += ( g[j][i]* ( vxx+vyy ) )- ( 2.0*f[j][i]*vxx ) + ( dthalbe*sumq );



	/* now updating the memory-variables and sum them up*/
	sumr=sump=sumq=0.0;
	for ( l=1; l<=L; l++ ) {
		r[j][i][l] = bip[l]* ( r[j][i][l]*cip[l]- ( dip[j][i][l]* ( vxy+vyx ) ) );
		p[j][i][l] = bjm[l]* ( p[j][i][l]*cjm[l]- ( e[j][i][l]* ( vxx+vyy ) ) + ( 2.0*d[j][i][l]*vyy ) );
		q[j][i][l] = bjm[l]* ( q[j][i][l]*cjm[l]- ( e[j][i][l]* ( vxx+vyy ) ) + ( 2.0*d[j][i][l]*vxx ) );
		sumr += r[j][i][l];
		sump += p[j][i][l];
		sumq += q[j][i][l];
	}


	/* and now the components of the stress tensor are
	   completely updated */
	sxy[j][i]+= ( dthalbe*sumr );
	sxx[j][i]+= ( dthalbe*sump );
	syy[j][i]+= ( dthalbe*sumq );

}