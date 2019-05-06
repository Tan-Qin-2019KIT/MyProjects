/*------------------------------------------------------------------------
 * Copyright (C) 2011 For the list of authors, see file AUTHORS.
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
--------------------------------------------------------------------------*/
/* ----------------------------------------------------------------------
 * Reading (distributed) source positions, timeshift, centre frequency
 * and amplitude from SOURCE_FILE.
 *
 * ---------------------------------------------------------------------- */

#include "fd.h"

float **sources(int *nsrc){

	/* declaration of extern variables */
	extern float PLANE_WAVE_DEPTH, PLANE_WAVE_ANGLE, TS, DH, SRCPOSXYZ[3];
	extern char SOURCE_FILE[STRING_SIZE];
	extern int MYID, NXG, NYG, SRCREC, RUNMODE, RUN_MULTIPLE_SHOTS, SOURCE_TYPE;
	extern FILE *FP;

	float **srcpos=NULL;
	int   i, l, isrc=0, current_source=0,nvarin=0;
	float xsrc, ysrc, tshift, tan_phi, dz, x, fc=0.0;
	char buffer[STRING_SIZE], bufferstring[10], cline[256];
	FILE *fpsrc;


	if (MYID==0){
		fprintf(FP," Message from function sources (written by PE %d):\n",MYID);
		if (SRCREC==1){ /* read source positions from file */
			*nsrc=0;
			fprintf(FP,"\n ------------------ READING SOURCE PARAMETERS ------------------- \n");
			fprintf(FP,"\n Reading source parameters from file: %s (SOFI2D source format)\n",SOURCE_FILE);
			if ((fpsrc=fopen(SOURCE_FILE,"r"))==NULL) declare_error(" Source file could not be opened !");
			while(fgets(buffer, STRING_SIZE, fpsrc))
			{
				sscanf(buffer,"%s",bufferstring);
				/* checks if the line contains a '%'character which indicates a comment line,
					and if the reading of a string was successful, which is not the case for an empty line*/
				/*if (sscanf(buffer,"%s",bufferstring)==1) printf("string %s \n",bufferstring);*/
				if ((strchr(buffer,'#')==0) && (sscanf(buffer,"%s",bufferstring)==1)) ++(*nsrc);
			}

			rewind(fpsrc);

			if ((nsrc)==0) fprintf(FP,"\n WARNING: Could not determine number of sources parameter sets in input file. Assuming %i.\n",(*nsrc=0));
			else fprintf(FP," Number of source positions specified in %s : %d \n",SOURCE_FILE,*nsrc);

			/* AUTO MODE: there is only 1 source => read first line, ignore the others */
			if (RUNMODE == 1)
				*nsrc = 1;

			/* memory for source position definition */
			srcpos=matrix(1,8,1,*nsrc);

			for (l=1;l<=*nsrc;l++){

				fgets(cline,255,fpsrc);
				nvarin=sscanf(cline,"%f%f%f%f%f%f%f",&xsrc, &ysrc, &tshift, &srcpos[5][l], &srcpos[6][l], &srcpos[7][l], &srcpos[8][l]);
				switch(nvarin){
				case 0: xsrc=0.0;
				case 1: ysrc=0.0;
				case 2: if (MYID==0) fprintf(FP," No time shift defined for source %i in %s!\n",l, SOURCE_FILE);
				declare_error("Missing parameter in SOURCE_FILE!");
				case 3: if (MYID==0) fprintf(FP," No frequency defined for source %i in %s!\n",l, SOURCE_FILE);
				declare_error("Missing parameter in SOURCE_FILE!");
				case 4: if (MYID==0) fprintf(FP," No amplitude defined for source %i in %s!\n",l, SOURCE_FILE);
				declare_error("Missing parameter in SOURCE_FILE!");
				case 5: srcpos[7][l]=0.0;
				case 6: srcpos[8][l]=SOURCE_TYPE;
				}
				if ((srcpos[8][l]!=4) && (nvarin>5)) {
					current_source=(int)srcpos[8][l];
					if (MYID==0) fprintf(FP," SOURCE_TYPE of source #%i is specified as %i, SOURCE_AZIMUTH is ignored.\n", l, current_source);
				}
				/* fscanf(fpsrc,"%f%f%f%f%f",&xsrc, &ysrc, &tshift, &fc, &amp); */

				/* AUTO MODE: use source positions from auto input file */
				if (RUNMODE == 1) {
					srcpos[1][l] = SRCPOSXYZ[0];
					srcpos[2][l] = SRCPOSXYZ[1];
					srcpos[3][l] = SRCPOSXYZ[2];
					/* otherwise: from source file */
				} 
				else {
					srcpos[1][l]=xsrc;
					srcpos[2][l]=ysrc;
					srcpos[3][l]=0.0;
				}
				srcpos[4][l]=tshift;
				fc=srcpos[5][l];
			}

			fclose(fpsrc);

			/* Compute maximum frequency */
			for (l=1;l<=*nsrc;l++)
				if (srcpos[5][l]>fc) fc=srcpos[5][l];
			fprintf(FP," Maximum frequency defined in %s: %6.2f Hz\n",SOURCE_FILE,fc);
			TS=1.0/fc;

			/* outputs all sources per each subdomain / node*/

			if (MYID==0){
				/*fprintf(FP," number\t    x\t\t    y\t\t  tshift\t    fc\t\t   amp\t	source_azimuth\tsource_type\n");

				for (l=1;l<=*nsrc;l++)
					fprintf(FP,"    %i \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f   \t %6.2f  \t   %1.0f\n\n",
							l, srcpos[1][l],srcpos[2][l],srcpos[4][l],srcpos[5][l],srcpos[6][l],srcpos[7][l],srcpos[8][l]);*/
				if (RUN_MULTIPLE_SHOTS) fprintf(FP," All sources will be modelled individually because of RUN_MULTIPLE_SHOTS=1!\n\n");
				else fprintf(FP," All sources will be modelled simultaneously because of RUN_MULTIPLE_SHOTS=0!\n\n");

			}

		} 
		else if (SRCREC==2) {
			if (PLANE_WAVE_DEPTH > 0) {  /* plane wave excitation */

				fprintf(FP," Computing source nodes for plane wave excitation.\n");
				fprintf(FP," depth= %5.2f meter, incidence angle= %5.2f degrees.\n",PLANE_WAVE_DEPTH, PLANE_WAVE_ANGLE);


				tan_phi=tan(PLANE_WAVE_ANGLE*PI/180.0);

				dz=(float)NXG*DH*tan_phi;
				fprintf(FP," Message from function sources (written by PE %d):\n",MYID);
				fprintf(FP," Maximum depth of plane wave: %5.2f meter \n",PLANE_WAVE_DEPTH+dz);
				if ((PLANE_WAVE_DEPTH+dz)<=NYG*DH){
					*nsrc=NXG;
					srcpos=matrix(1,8,1,*nsrc);
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
						srcpos[7][isrc]=0.0;
						srcpos[8][isrc]=SOURCE_TYPE;
					}
				}
				else declare_error(" Maximum depth of plane wave exceeds model depth. ");
			}
			else declare_error("SRCREC parameter specifies PLANE_WAVE excitation, but PLANE_WAVE_DEPTH<=0!");
		}
		else declare_error("SRCREC parameter is invalid (SRCREC!=1 or SRCREC!=2)! No source parameters specified!");
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(nsrc,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&TS,1,MPI_FLOAT,0,MPI_COMM_WORLD);

	if (MYID!=0) srcpos=matrix(1,8,1,*nsrc);
	MPI_Bcast(&srcpos[1][1],(*nsrc)*8,MPI_FLOAT,0,MPI_COMM_WORLD);

	if (MYID==0){
		if (*nsrc>50) fprintf(FP," The following table is quite large (%i lines) and will, thus, be truncated to the first 50 entries! \n\n",*nsrc);
		fprintf(FP," number\t    x\t\t    y\t\t  tshift\t    fc\t\t   amp\t	source_azimuth\tsource_type\n");

		if (*nsrc>50) { for (l=1;l<=50;l++)
			fprintf(FP,"    %i \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f   \t %6.2f  \t   %1.0f\n",
					l, srcpos[1][l],srcpos[2][l],srcpos[4][l],srcpos[5][l],srcpos[6][l],srcpos[7][l],srcpos[8][l]);
		}
		else for (l=1;l<=*nsrc;l++)
					fprintf(FP,"    %i \t %6.2f \t %6.2f \t %6.2f \t %6.2f \t %6.2f   \t %6.2f  \t   %1.0f\n",
							l, srcpos[1][l],srcpos[2][l],srcpos[4][l],srcpos[5][l],srcpos[6][l],srcpos[7][l],srcpos[8][l]);
	}

	return srcpos;
}
