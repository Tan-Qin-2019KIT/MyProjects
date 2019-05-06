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

/*
 * Normalize Gradient
 *
 */

#include "fd.h"

float norm(float ** waveconv, int iter, int sws)
{

	/* extern variables */

    extern float DH;
	extern int FREE_SURF, NX, NY, NXG, NYG;
	extern int NPROCX, NPROCY, MYID, POS[3];
	extern char JACOBIAN[STRING_SIZE];
	extern FILE *FP;
	
	/* local variables */
    int i, j;

	float min2, max2, norm;
	float tmp_max, tmp_min;
	
	min2 = waveconv[1][1];
	max2 = waveconv[1][1];
	
	/* estimate min and max of the global gradient */
	 for (i=1;i<=NX;i++){
	    for (j=1;j<=NY;j++){
	            
	       if(waveconv[j][i]<min2){min2=waveconv[j][i];}
	       if(waveconv[j][i]>max2){max2=waveconv[j][i];}
	            
	    }
	 }
	
   /* Compute maximum value on all CPUs */
   tmp_max = 0.0;
   MPI_Allreduce(&max2,&tmp_max,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
   max2 = tmp_max;

   /* Compute minimum value on all CPUs */
   tmp_min = 0.0;
   MPI_Allreduce(&min2,&tmp_min,1,MPI_FLOAT,MPI_MIN,MPI_COMM_WORLD);
   min2 = tmp_min;

	/* Normalize gradient */
	norm = max2;
	
	/*if(sqrt(min2*min2)>sqrt(max2*max2)){norm=sqrt(min2*min2);}*/
	
	for (i=1;i<=NX;i++){
	    for (j=1;j<=NY;j++){	    	    
	       waveconv[j][i] = waveconv[j][i]/norm;	    
	    }
	}                
	
    printf("MYID=%d \t norm_fac=%e \n",MYID,norm);
    return norm;
}



