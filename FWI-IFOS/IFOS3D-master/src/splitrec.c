/*------------------------------------------------------------------------
 * Copyright (C) 2015 For the list of authors, see file AUTHORS.
 *
 * This file is part of IFOS3D.
 * 
 * IFOS3D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 * 
 * IFOS3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with IFOS3D. See file COPYING and/or 
 * <http://www.gnu.org/licenses/gpl-2.0.html>.
--------------------------------------------------------------------------*/

/*-------------------------------------------------------------------
 *  Computation of local receiver coordinates (within each subgrid)
  --------------------------------------------------------------------*/
 
 #include "fd.h"
int **splitrec(int **recpos,int *ntr_loc, int ntr, int *rnum_loc)
{

	extern int IENDX, IENDY, IENDZ, MYID, POS[4], VERBOSE;
	extern FILE *FP;

	int a,b,c,i=0,j,k;
	int ** recpos_dummy, **recpos_local=NULL;
	recpos_dummy = imatrix(1,4,1,ntr);

	for (j=1;j<=ntr;j++) {
		a=(recpos[1][j]-1)/IENDX;
		b=(recpos[2][j]-1)/IENDY;
		c=(recpos[3][j]-1)/IENDZ;


		if ((POS[1]==a)&&(POS[2]==b)&&(POS[3]==c)) {
			i++;
			recpos_dummy[1][i] = ((recpos[1][j]-1)%IENDX)+1;
			recpos_dummy[2][i] = ((recpos[2][j]-1)%IENDY)+1;
			recpos_dummy[3][i] = ((recpos[3][j]-1)%IENDZ)+1;
			recpos_dummy[4][i] = j;
			rnum_loc[j]        =1;
		}
		else 	rnum_loc[j]	   =0;
	}
   
	if (i>0) recpos_local = imatrix(1,4,1,i);
	for (k=1;k<=i;k++){
		recpos_local[1][k] = recpos_dummy[1][k];
		recpos_local[2][k] = recpos_dummy[2][k];
		recpos_local[3][k] = recpos_dummy[3][k];
		recpos_local[4][k] = recpos_dummy[4][k];
	}
	free_imatrix(recpos_dummy,1,4,1,ntr);
	if(VERBOSE){
	fprintf(FP,"\n **Message from split_rec:\n");
	fprintf(FP," Splitting of receivers from global to local grids finished.\n");
	fprintf(FP," MYID= %d \t \t no. of receivers= %d\n",MYID,i);
	}
/*
	fprintf(FP,"\n **Message from split_rec:\n");
	fprintf(FP," Table of local receiver positions (in gridpoints):\n");
	fprintf(FP," MYID \t x \t y \t z\n");

	for (j=1;j<=i;j++)
		fprintf(FP," %d \t %d \t %d \t %d \n",
		    MYID,recpos_local[1][j],recpos_local[2][j],recpos_local[3][j]);
*/
   *ntr_loc=i;
	return recpos_local;

}
