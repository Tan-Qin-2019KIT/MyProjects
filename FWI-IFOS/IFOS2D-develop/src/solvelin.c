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
 *------------------------------------------------------------------------
 *
 *   Solve linear equation systems
 *
 *------------------------------------------------------------------------
*/


#include "fd.h"

void solvelin(float  **AA, float *bb, float *x, int e, int method)
{

	/* local variables */
	int k, m, n, rows, columns;
	float a, c, **A, *b;

	
	rows = e;
	columns = e;

	A = matrix(1,rows,1,columns);
	b = vector(1,rows);

	/* temporary variables */
	for (k=1;k<=rows;k++) {
		for (n=1;n<=columns;n++)  A[k][n] = AA[k][n];
		b[k] = bb[k];
	}


	switch (method)
	{
	case 1:	/* Gauﬂ algorithm */
	{
		for (k=1;k<=rows-1;k++)
			for (n=k;n<=rows-1;n++)
			{
				a = A[n+1][k]/A[k][k];
				for (m=1;m<=columns;m++) A[n+1][m] = A[n+1][m] - a*A[k][m];
				b[n+1] = b[n+1] - a*b[k];
			}
		
		for (k=rows;k>=1;k--)
		{
			c = b[k];
			for (m=columns;m>=k+1;m--) c = c - A[k][m]*x[m];
			x[k] = c/A[k][k];
		}
		break;
	} /* END of case Gauﬂ */
		
	} /* END of switch (method) */
	
	
	free_matrix(A,1,rows,1,columns);
	free_vector(b,1,rows);

	return;
}




