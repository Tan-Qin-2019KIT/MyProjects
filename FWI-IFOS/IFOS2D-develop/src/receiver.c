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
 *  compute receiver positions or read them from external file        
 *
 *  See COPYING file for copying and redistribution conditions.
 *  ----------------------------------------------------------------------*/

#include "fd.h"

int **receiver(int *ntr, float **srcpos, int shotno){

	/* declaration of extern variables */
	extern  char REC_FILE[STRING_SIZE];
	extern float DH, REFREC[4], REC_ARRAY_DEPTH, REC_ARRAY_DIST, FW;
	extern int READREC, NGEOPH, DRX, REC_ARRAY, NXG, NYG;
	extern int MYID, VERBOSE,TRKILL;
	extern FILE *FP;

	int **recpos1, **recpos = NULL, nxrec=0, nyrec=0, nzrec=0, recdist;
	int   itr=1, itr1=0, itr2=0, recflag=0, c, ifw, n, i, j, ntr1=0;
	int nxrec1, nxrec2, nxrec3, nyrec1, nyrec2, nzrec1=1, nzrec2=1, offset, shotdist;
	extern float XREC1, YREC1, XREC2, YREC2;
	float xrec, yrec;
	FILE *fpr,*f;
	char filename[STRING_SIZE];

	if (MYID==0)
	{
     	if (READREC==1){ /* read receiver positions from file */
     		fprintf(FP,"\n Reading receiver positions from file: \n\t%s\n",REC_FILE);
		    fpr=fopen(REC_FILE,"r");
     		if (fpr==NULL) declare_error(" Receiver file could not be opened !");
     		*ntr=0;
     		while ((c=fgetc(fpr)) != EOF)
     			if (c=='\n') ++(ntr1);
     		rewind(fpr);
     
     		recpos1=imatrix(1,3,1,ntr1);
     		for (itr=1;itr<=ntr1;itr++){
     			fscanf(fpr,"%f%f\n",&xrec, &yrec);
     			recpos1[1][itr]=iround((xrec+REFREC[1])/DH);
     			recpos1[2][itr]=iround((yrec+REFREC[2])/DH);
     			recpos1[3][itr]=iround((0.0+REFREC[3])/DH);
     		}
     		fclose(fpr);
     		fprintf(FP," Message from function receiver (written by PE %d):\n",MYID);/***/
     		fprintf(FP," Number of receiver positions found: %i\n",ntr1);
            
            /* recpos1[1][itr] X in m
             * recpos1[2][itr] Y in m
             *
             */
            
     		/* check if more than one receiver is located
     				         at the same gridpoint */
     		for (itr=1;itr<=(ntr1-1);itr++)
     			for (itr1=itr+1;itr1<=ntr1;itr1++)
     				if ((recpos1[1][itr]==recpos1[1][itr1])
     				    && (recpos1[2][itr]==recpos1[2][itr1])
     				    && (recpos1[3][itr]==recpos1[3][itr1]))
     					recpos1[1][itr1]=-(++recflag);
     
     		recpos=imatrix(1,3,1,ntr1-recflag);
     		for (itr=1;itr<=ntr1;itr++)
     			if (recpos1[1][itr]>0){
     				recpos[1][++itr2]=recpos1[1][itr];
     				recpos[2][itr2]=recpos1[2][itr];
     				recpos[3][itr2]=recpos1[3][itr];
     			}
     
     		*ntr=itr2;
     		if (recflag>0){
     			fprintf(FP,"\n\n");
     			fprintf(FP," Warning:\n");
     			fprintf(FP," Several receivers located at the same gridpoint !\n");
     			fprintf(FP," Number of receivers reduced to %i\n", *ntr);
     			fprintf(FP,"\n\n");
     		}
     
     		free_imatrix(recpos1,1,3,1,ntr1);
  
     	}
     
     	else if (READREC==2){
		fprintf(FP,"\n Moving streamer acquisition geometry choosen. \n");
		/* decide if XREC1 or XREC2 is the nearest receiver to shot 1 */
		if (abs(XREC1-srcpos[1][1])<abs(XREC2-srcpos[1][1])){
			nxrec1=iround(XREC1/DH);   /* (nxrec1,nyrec1) and (nxrec2,nyrec2) are */
			nyrec1=iround(YREC1/DH);   /* the positions of the first and last receiver*/
			nxrec2=iround(XREC2/DH);   /* in gridpoints */
			nyrec2=iround(YREC2/DH);
			offset=iround((XREC1-srcpos[1][1])/DH);
		}
		else {
			nxrec1=iround(XREC2/DH);   /* (nxrec1,nyrec1) and (nxrec2,nyrec2) are */
			nyrec1=iround(YREC2/DH);   /* the positions of the first and last receiver*/
			nxrec2=iround(XREC1/DH);   /* in gridpoints */
			nyrec2=iround(YREC1/DH);
			offset=iround((XREC2-srcpos[1][1])/DH);
		}
		if (nyrec1 != nyrec2){
			fprintf(FP,"\n\n");
			fprintf(FP," Warning:\n");
			fprintf(FP," Receiver positions are not in the same depth !\n");
			fprintf(FP," Depth of all receivers is seth to %i m\n",YREC1);
			fprintf(FP,"\n\n");
		}
		if (offset<0) recdist=-NGEOPH;
		else recdist=NGEOPH;
		*ntr=iround((nxrec2-nxrec1)/recdist)+1;
		if (shotno>0) nxrec1=nxrec1+((srcpos[1][shotno]-srcpos[1][1])/DH);
		recpos=imatrix(1,3,1,*ntr);
		n=0;
		
		for (itr=1;itr<=*ntr;itr++){
			nxrec=nxrec1+(itr-1)*recdist;
			if (nxrec>FW && nxrec<(NXG-FW) && nyrec1>0 && nyrec1<(NYG-FW)){
				recpos[1][itr]=nxrec;
				recpos[2][itr]=nyrec1;
			} else declare_error(" Receiver positions of current shot are out of model boundaries !");
		}

		if (shotno>0){
			/* write receiver positions to file */
			sprintf(filename,"%s.shot%i",REC_FILE,shotno);
			f=fopen(filename,"wb");
			for (i=1;i<=*ntr;i++) {
				fprintf(f,"%f \t %f \n",recpos[1][i]*DH,recpos[2][i]*DH);
			}
			fclose(f);
		}
			
		if (VERBOSE==1) {
		  fprintf(FP," Message from function receiver (written by PE %d):\n",MYID);
		  fprintf(FP," Number of receiver positions found: %i\n",*ntr);
		  fprintf(FP," First receiver position for shot %i: %f (x), %f (y)\n",shotno,recpos[1][1]*DH,recpos[2][1]*DH);
		  fprintf(FP," Last receiver position for shot %i: %f (x), %f (y)\n",shotno,recpos[1][*ntr]*DH,recpos[2][*ntr]*DH);  
		}
     	}
     
     	else{         /* straight horizontal or vertical
     		                     line of receivers */
			nxrec1=iround(XREC1/DH);   /* (nxrec1,nyrec1) and (nxrec2,nyrec2) are */
			nyrec1=iround(YREC1/DH);   /* the positions of the first and last receiver*/
			nxrec2=iround(XREC2/DH);	 /* in gridpoints */
			nyrec2=iround(YREC2/DH);

     		if ((abs(nyrec2-nyrec1)<=abs(nxrec2-nxrec1))||
     		    (abs(nyrec2-nyrec1)<=abs(nzrec2-nzrec1))){
     			if (abs(nzrec2-nzrec1)<=abs(nxrec2-nxrec1)){
     				/* geophone-array horizontal x-dirextion */
     				*ntr=iround((nxrec2-nxrec1)/NGEOPH)+1;
     				recpos=imatrix(1,3,1,*ntr);
     				for (nxrec=nxrec1;nxrec<=nxrec2;nxrec+=NGEOPH){
     					nyrec=nyrec1+((nyrec2-nyrec1)/(nxrec2-nxrec1)*(nxrec-nxrec1));
     					nzrec=nzrec1+((nzrec2-nzrec1)/(nxrec2-nxrec1)*(nxrec-nxrec1));
     					itr=iround((nxrec-nxrec1)/NGEOPH)+1;
     					recpos[1][itr]=nxrec;
     					recpos[2][itr]=nyrec;
     					recpos[3][itr]=nzrec;
     				}
     			}

     		}
     		else{         /* receiver-line vertical */
     			*ntr=iround((nyrec2-nyrec1)/NGEOPH)+1;
     			recpos=imatrix(1,3,1,*ntr);
     			for (nyrec=nyrec1;nyrec<=nyrec2;nyrec+=NGEOPH){
     				nxrec=nxrec1+((nxrec2-nxrec1)/(nyrec2-nyrec1)*(nyrec-nyrec1));
     				/**wird nzrec noch gebraucht?**/
     				nzrec=nzrec1+((nzrec2-nzrec1)/(nyrec2-nyrec1)*(nyrec-nyrec1));
     				itr=iround((nyrec-nyrec1)/NGEOPH)+1;
     				recpos[1][itr]=nxrec;
     				recpos[2][itr]=nyrec;
     				recpos[3][itr]=nzrec;
     			}
     		}
     
     	}
     	/*   fprintf(fp,"Gridpoints of receiver positions (x,y,z):\n");
     		for (itr=1;itr<=*ntr;itr++)
     				fprintf(fp,"%i\t%i\t%i\n",recpos[1][itr],recpos[2][itr],recpos[3][itr]);*/
     
    

	}

      
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(ntr,1,MPI_INT,0,MPI_COMM_WORLD);
	if (MYID!=0) recpos=imatrix(1,3,1,*ntr);
	MPI_Bcast(&recpos[1][1],(*ntr)*3,MPI_INT,0,MPI_COMM_WORLD);

/*	if (MYID==0)
	{
		fprintf(fp,"\n **Message from function receiver (written by PE %d):\n",MYID);
		fprintf(fp," Number of receiver positions found: %i\n",*ntr);
		fprintf(fp," Receiver positions (in gridpoints) in the global model-system:\n");
		fprintf(fp," x  \ty \n");
		fprintf(fp," -  \t- \n");
		for (l=1;l<=*ntr;l++)
			fprintf(fp," %i\t%i\n",recpos[1][l],recpos[2][l]);
		fprintf(fp,"\n\n");
	}
*/

	return recpos;
}
