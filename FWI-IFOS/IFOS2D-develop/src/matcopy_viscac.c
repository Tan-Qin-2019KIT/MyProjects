 /*-----------------------------------------------------------------------------------------
 * Copyright (C) 2013  For the list of authors, see file AUTHORS.
 *
 * This file is part of DENISE.
 * 
 * DENISE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 * 
 * DENISE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with DENISE. See file COPYING and/or <http://www.gnu.org/licenses/gpl-2.0.html>.
-----------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
 * For the averaging of material properties each process requires values
 * at the indices 0 and NX+1 etc. These lie on the neighbouring processes.
 * Thus, they have to be copied which is done by this function.
 *
 *   last update 12/02/02, T. Bohlen
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void matcopy_viscac(float ** rho, float ** pi, float ** taup){

	extern int MYID, NX, NY, INDEX[5];
	extern const int TAG1,TAG2,TAG5,TAG6;
	extern FILE *FP;


	MPI_Status status;	
	double time1, time2;	
	int i, j;
	float ** bufferlef_to_rig, ** bufferrig_to_lef;
	float ** buffertop_to_bot, ** bufferbot_to_top;

	bufferlef_to_rig = matrix(0,NY+1,1,5);
	bufferrig_to_lef = matrix(0,NY+1,1,5);
	buffertop_to_bot = matrix(0,NX+1,1,5);
	bufferbot_to_top = matrix(0,NX+1,1,5);
	
	
	fprintf(FP,"\n\n **Message from matcopy (written by PE %d):",MYID);
	fprintf(FP,"\n Copy material properties at inner boundaries ... \n");
	time1=MPI_Wtime();




//	if (POS[2]!=0)	/* no boundary exchange at top of global grid */
	for (i=0;i<=NX+1;i++){
			/* storage of top of local volume into buffer */
			buffertop_to_bot[i][1]  =  rho[1][i];
			buffertop_to_bot[i][2]  =  pi[1][i];
			buffertop_to_bot[i][5]  = taup[1][i];
	}


//	if (POS[2]!=NPROCY-1)	/* no boundary exchange at bottom of global grid */
	for (i=0;i<=NX+1;i++){
			
			/* storage of bottom of local volume into buffer */
			bufferbot_to_top[i][1]  =  rho[NY][i];
			bufferbot_to_top[i][2]  =  pi[NY][i];
			bufferbot_to_top[i][5]  = taup[NY][i];
	}


 	/*=========sending and receiving of the boundaries==========*/

	MPI_Bsend(&buffertop_to_bot[1][1],NX*5,MPI_FLOAT,INDEX[3],TAG5,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Recv(&buffertop_to_bot[1][1],NX*5,MPI_FLOAT,INDEX[4],TAG5,MPI_COMM_WORLD,&status);
	MPI_Bsend(&bufferbot_to_top[1][1],NX*5,MPI_FLOAT,INDEX[4],TAG6,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Recv(&bufferbot_to_top[1][1],NX*5,MPI_FLOAT,INDEX[3],TAG6,MPI_COMM_WORLD,&status);   


//	if (POS[2]!=NPROCY-1)	/* no boundary exchange at bottom of global grid */
	for (i=0;i<=NX+1;i++){
			rho[NY+1][i] = 	buffertop_to_bot[i][1];
			pi[NY+1][i] = 	buffertop_to_bot[i][2];
			taup[NY+1][i] = buffertop_to_bot[i][5];
	}

//	if (POS[2]!=0)	/* no boundary exchange at top of global grid */
	for (i=0;i<=NX+1;i++){
			rho[0][i] = 	bufferbot_to_top[i][1];
			pi[0][i] = 	bufferbot_to_top[i][2];
			taup[0][i] = 	bufferbot_to_top[i][5];
	}




//	if (POS[1]!=0)	/* no boundary exchange at left edge of global grid */
		for (j=0;j<=NY+1;j++)
		{
			/* storage of left edge of local volume into buffer */
			bufferlef_to_rig[j][1] =  rho[j][1];
			bufferlef_to_rig[j][2] =  pi[j][1];
			bufferlef_to_rig[j][5] =  taup[j][1];
		}


//	if (POS[1]!=NPROCX-1)	/* no boundary exchange at right edge of global grid */
	for (j=0;j<=NY+1;j++){
			/* storage of right edge of local volume into buffer */
			bufferrig_to_lef[j][1] =  rho[j][NX];
			bufferrig_to_lef[j][2] =  pi[j][NX];
			bufferrig_to_lef[j][5] =  taup[j][NX];
	}



 	MPI_Bsend(&bufferlef_to_rig[0][1],(NY+2)*5,MPI_FLOAT,INDEX[1],TAG1,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Recv(&bufferlef_to_rig[0][1],(NY+2)*5,MPI_FLOAT,INDEX[2],TAG1,MPI_COMM_WORLD,&status);
	MPI_Bsend(&bufferrig_to_lef[0][1],(NY+2)*5,MPI_FLOAT,INDEX[2],TAG2,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Recv(&bufferrig_to_lef[0][1],(NY+2)*5,MPI_FLOAT,INDEX[1],TAG2,MPI_COMM_WORLD,&status);


//	if (POS[1]!=NPROCX-1)	/* no boundary exchange at right edge of global grid */
	for (j=0;j<=NY+1;j++){
			rho[j][NX+1] = 	bufferlef_to_rig[j][1];
			pi[j][NX+1] = 	bufferlef_to_rig[j][2];
			taup[j][NX+1] = bufferlef_to_rig[j][5];
	}

//	if (POS[1]!=0)	/* no boundary exchange at left edge of global grid */
	for (j=0;j<=NY+1;j++){
			rho[j][0] = 	bufferrig_to_lef[j][1];
			pi[j][0] = 	bufferrig_to_lef[j][2];
			taup[j][0] = 	bufferrig_to_lef[j][5];
	}


	if (MYID==0){
		time2=MPI_Wtime();
		fprintf(FP," finished (real time: %4.2f s).\n",time2-time1);
	}

	free_matrix(bufferlef_to_rig,0,NY+1,1,5);
	free_matrix(bufferrig_to_lef,0,NY+1,1,5);
	free_matrix(buffertop_to_bot,0,NX+1,1,5);
	free_matrix(bufferbot_to_top,0,NX+1,1,5);
}
