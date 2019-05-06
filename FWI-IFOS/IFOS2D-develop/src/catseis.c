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
 *-------------------------------------------------------------
 *
 * Cat seismograms
 *  
 *
 *-------------------------------------------------------------
*/

#include "fd.h"

void	catseis(float **data, float **fulldata, int *recswitch, int ntr_glob, MPI_Comm newcomm_nodentr) {

	extern int	NT, MYID;
	
	int		i, j, k;
	float		**fulldata2;

	/* temporary global data array for MPI-exchange */
	fulldata2 = matrix(1,ntr_glob,1,NT);
	
	k = 0;	/* trace counter for local data array */

	/* loop over global traces: copy traces of local array	*/
	/* to appropriate locations in the global array		*/
	for(i=1;i<=ntr_glob;i++)
{	
	if (recswitch[i]) {
			k++;
			for(j=1;j<=NT;j++)	fulldata2[i][j] = data[k][j];
		}}

	MPI_Allreduce(&fulldata2[1][1], &fulldata[1][1], ntr_glob*NT, MPI_FLOAT, MPI_SUM, newcomm_nodentr);
	
	free_matrix(fulldata2, 1,ntr_glob,1,NT);
}
