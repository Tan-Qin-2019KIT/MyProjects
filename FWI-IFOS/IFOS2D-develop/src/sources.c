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
   Reading (distributed) source positions, timeshift, centre frequency 
   and amplitude from SOURCE_FILE.
   
*/

#include "fd.h"

float **sources(int *nsrc){

	/* declaration of extern variables */
	extern float PLANE_WAVE_DEPTH, PHI, TS, DH, F_REF;
	extern  char SOURCE_FILE[STRING_SIZE];
	extern int MYID, NXG, NYG, SRCREC, RUN_MULTIPLE_SHOTS, SOURCE_TYPE;
	extern FILE *FP;

	float **srcpos = NULL;
	int   i, l, isrc=0, current_source=0, current_azimuth=0, nvarin=0;
	float xsrc, ysrc, zsrc, tshift, fc = 0.0, tan_phi, dz, x;
	char  cline[256];
	FILE *fpsrc;


	if (MYID==0){
		if (SRCREC){ /* read source positions from file */
			fprintf(FP,"\n Reading source positions, time-shift, centre frequency \n");
			fprintf(FP," and amplitude from file: %s\n",SOURCE_FILE);
			fpsrc=fopen(SOURCE_FILE,"r");
	
			if (fpsrc==NULL) declare_error(" Source file could not be opened !");
			*nsrc=0;

			
			/* read number of source positions */	
			fscanf(fpsrc,"%d",nsrc);
						
			fprintf(FP," Number of source positions specified in %s : %d \n",SOURCE_FILE,*nsrc);
			
			srcpos=matrix(1,8,1,*nsrc);

			rewind(fpsrc);		
			fgets(cline,255,fpsrc);		/* Dummy fgets for ignoring first line */

			for (l=1;l<=*nsrc;l++){
				fgets(cline,255,fpsrc);

				nvarin=sscanf(cline,"%f%f%f%f%f%f%f%f",&xsrc, &zsrc, &ysrc, &tshift, &srcpos[5][l], &srcpos[6][l], &srcpos[7][l], &srcpos[8][l]);
				switch(nvarin){
					case 0: xsrc=0.0;
					case 1: zsrc=0.0;
					case 2: ysrc=0.0;
					case 3: if (MYID==0) fprintf(FP," No time shift defined for source %i in %s!\n",l, SOURCE_FILE);
						declare_error("Missing parameter in SOURCE_FILE!");
					case 4: if (MYID==0) fprintf(FP," No frequency defined for source %i in %s!\n",l, SOURCE_FILE);
						declare_error("Missing parameter in SOURCE_FILE!");
					case 5: if (MYID==0) fprintf(FP," No amplitude defined for source %i in %s!\n",l, SOURCE_FILE);
						declare_error("Missing parameter in SOURCE_FILE!");
					case 6: srcpos[7][l]=0.0;
					case 7: srcpos[8][l]=SOURCE_TYPE;
                        
				}
				if ((srcpos[8][l]!=4) && (nvarin>6)) {
				current_source=(int)srcpos[8][l];
				if (MYID==0) fprintf(FP," SOURCE_TYPE of source #%i is specified as %i, SOURCE_AZIMUTH is ignored.\n", l, current_source);
				}
				if ((srcpos[7][l]==0) && (srcpos[8][l]==4)) {
				current_azimuth=(int)srcpos[7][l];
				if (MYID==0) fprintf(FP,"\n WARNING: SOURCE_TYPE of source #%i is specified as rotated force, but no SOURCE_AZIMUTH is specified in the source file! The SOURCE_AZIMUTH is set to %i degrees.\n",l,current_azimuth);
				}
				/* fscanf(fpsrc,"%f%f%f%f%f",&xsrc, &ysrc, &tshift, &fc, &amp); */ 
				srcpos[1][l]=xsrc;
				srcpos[2][l]=ysrc;
				srcpos[3][l]=0.0;
				srcpos[4][l]=tshift;
				fc=srcpos[5][l];
			}

			fclose(fpsrc);

			/* Compute maximum frequency */
			for (l=1;l<=*nsrc;l++)
				if (srcpos[5][l]>fc) fc=srcpos[5][l];
			fprintf(FP," Maximum frequency defined in %s: %6.2e Hz\n",SOURCE_FILE,fc);
			TS=1.0/fc;

			/* outputs all sources per each subdomain / node*/
		
			if (MYID==0){
				fprintf(FP," number\t    x\t\t    y\t\t  tshift\t    fc\t\t   amp\t	source_azimuth\tsource_type\n");
	   			for (l=1;l<=*nsrc;l++)
		      			fprintf(FP,"    %i \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f   \t %6.2f  \t   %1.0f\n\n",
					l, srcpos[1][l],srcpos[2][l],srcpos[4][l],srcpos[5][l],srcpos[6][l],srcpos[7][l],srcpos[8][l]);
				if (RUN_MULTIPLE_SHOTS) fprintf(FP," All sources will be modelled individually because of RUN_MULTIPLE_SHOTS=1!\n");
				else fprintf(FP," All sources will be modelled simultaneously because of RUN_MULTIPLE_SHOTS=0!\n");
	      		}

            /*
             * srcpos[1][l] X-Position in m
             * srcpos[2][l] Y-Position in m
             * srcpos[3][l] tshift
             * srcpos[4][l] fv
             * srcpos[5][l] amp
             * srcpos[6][l] azimuth
             * srcpos[7][l] type
             */
		} 
    else if (PLANE_WAVE_DEPTH > 0) {  /* plane wave excitation */
				fprintf(FP," Computing source nodes for plane wave excitation.\n");
				fprintf(FP," depth= %5.2f meter, incidence angle= %5.2f degrees.\n",PLANE_WAVE_DEPTH, PHI);

	
				tan_phi=tan(PHI*PI/180.0);
				
				dz=(float)NXG*DH*tan_phi;
				fprintf(FP," Message from function sources (written by PE %d):\n",MYID);				
				fprintf(FP," Maximum depth of plane wave: %5.2f meter \n",PLANE_WAVE_DEPTH+dz);				
				if ((PLANE_WAVE_DEPTH+dz)<=NYG*DH){			
					*nsrc=NXG;
					srcpos=matrix(1,6,1,*nsrc);
					isrc=0;
					for (i=1;i<=NXG;i++){
						isrc++;
						x=(float)i*DH;
						srcpos[1][isrc]=x;
						srcpos[2][isrc]=PLANE_WAVE_DEPTH+(tan_phi*x);
						srcpos[3][isrc]=0.0;
						srcpos[4][isrc]=0.0;
						srcpos[5][isrc]=1.0/TS;
						srcpos[6][isrc]=1.0;
					}
				}
				else declare_error(" Maximum depth of plane wave exceeds model depth. ");
			}
				
			fprintf(FP," Message from function sources (written by PE %d):\n",MYID);
			
			if (F_REF<0)	F_REF=1.0/TS;			

		}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(nsrc,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&TS,1,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Bcast(&F_REF,1,MPI_FLOAT,0,MPI_COMM_WORLD);
	if (MYID!=0) srcpos=matrix(1,8,1,*nsrc);
	MPI_Bcast(&srcpos[1][1],(*nsrc)*8,MPI_FLOAT,0,MPI_COMM_WORLD);

/*	if (MYID==0){
		fprintf(FP,"\n **Message from function source (written by PE %d):\n",MYID);
		fprintf(FP," Number of global source positions found: %i\n",*nsrc);
		fprintf(FP," x\t\ty\t\tz\t\ttshift\t\tfc\t\tamp\n");
		for (l=1;l<=*nsrc;l++)
			fprintf(FP," %6.2e\t%6.2e\t%6.2e\t%6.2e\t%6.2e\t%6.2e\n",
					srcpos[1][l],srcpos[2][l],srcpos[3][l],srcpos[4][l],srcpos[5][l],srcpos[6][l]);
		fprintf(FP,"\n\n");
	}
*/

	return srcpos;
}
