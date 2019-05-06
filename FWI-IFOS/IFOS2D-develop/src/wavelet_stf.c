/*-----------------------------------------------------------------------------------------
 * Copyright (C) 2016  For the list of authors, see file AUTHORS.
 *
 * This file is part of IFOS.
 * 
 * IFOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 * 
 * IFOS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with IFOS. See file COPYING and/or <http://www.gnu.org/licenses/gpl-2.0.html>.
-----------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
*   Passing the source time function of stf.c to signals
*  ----------------------------------------------------------------------*/

#include "fd.h"


float ** wavelet_stf(int nsrc, int ishot, float ** signals_stf){


	/* extern variables */
	extern int SOURCE_SHAPE, NT, MYID, INV_STF;
	extern float  DT;
	extern FILE *FP;

	/*local variables */
	int nt, k, z=1;
	float ** signals;
	
	signals=fmatrix(1,nsrc,1,NT);
		
		for (k=1;k<=nsrc;k++){
			for (nt=1;nt<=NT;nt++) {
						signals[z][nt]=signals_stf[k][nt];	
						}
					++z;	
					}
	
	fprintf(FP," Message from function wavelet_stf written by PE %d \n",MYID);
	fprintf(FP," %d source positions located in subdomain of PE %d \n",nsrc,MYID);
	fprintf(FP," have been assigned with a source signal out of source time function. \n");
			
	return signals;	

}
