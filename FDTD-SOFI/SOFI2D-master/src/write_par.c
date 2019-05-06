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
/*------------------------------------------------------------------------
 *   Write FD-Parameters to stdout                           
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

/* printing all important parameters on stdout */
void write_par(FILE *fp){

	/* declaration of extern variables */
	extern int   NX, NY, NT, SOURCE_SHAPE, SOURCE_TYPE;
	extern int  SNAP, SNAP_FORMAT, L, SRCREC;
	extern float DH, TIME, DT, TS, *FL, TAU, DAMPING, PLANE_WAVE_DEPTH, PLANE_WAVE_ANGLE;
	extern float REC_ARRAY_DEPTH, REC_ARRAY_DIST;
	extern float FPML, VPPML, NPOWER, K_MAX_CPML;
	extern float XREC1, XREC2, YREC1, YREC2;
	/*extern int NGEOPH;*/
	extern float NGEOPH;
	extern int SEISMO, NDT, SEIS_FORMAT, FREE_SURF, ABS_TYPE, FW;
	extern int  READMOD, READREC, BOUNDARY, REC_ARRAY, DRX, FDORDER, RSG;
	extern float TSNAP1, TSNAP2, TSNAPINC, REFREC[4];
	extern char SNAP_FILE[STRING_SIZE], SOURCE_FILE[STRING_SIZE], REC_FILE[STRING_SIZE];
	extern char SEIS_FILE[STRING_SIZE];
	extern char SIGNAL_FILE[STRING_SIZE];
	extern char  MFILE[STRING_SIZE];
	extern int NP, NPROCX, NPROCY, MYID,FDORDER_TIME;

	/* definition of local variables */
	int l;
	char file_ext[8];

	fprintf(fp,"\n **********************************************************");
	fprintf(fp,"\n ********* PARAMETERS AS SPECIFIED IN INPUT FILE **********");
	fprintf(fp,"\n **********************************************************\n\n");

	fprintf(fp,"\n **Message from write_par (printed by PE %d):\n",MYID);
	fprintf(fp,"\n");
	fprintf(fp,"------------------------- Processors ------------------------\n");
	fprintf(fp," Number of PEs in x-direction (NPROCX): %d\n",NPROCX);
	fprintf(fp," Number of PEs in vertical direction (NPROCY): %d\n",NPROCY);
	fprintf(fp," Total number of PEs in use: %d\n",NP);
	fprintf(fp,"\n");

	fprintf(fp,"------------------------- FD Algorithm ------------------------\n");
	if (RSG) {
		fprintf(fp," Rotated Staggered Grid (RSG) is used. \n");
		fprintf(fp," Order of spatial FD operators (FDORDER) is set to 2.");

	} else 	{
		fprintf(fp," Standard Staggered Grid (SSG) (Virieux-grid) is used. \n");
		fprintf(fp," Order of spatial FD operators (FDORDER) is %d\n", FDORDER);
        fprintf(fp," Order of temporal FD operator (FDORDER_TIME) is %d\n", FDORDER_TIME);
	}
	fprintf(fp,"\n");


	fprintf(fp," ----------------------- Discretization  ---------------------\n");
	fprintf(fp," Number of gridpoints in x-direction (NX): %i\n", NX);
	fprintf(fp," Number of gridpoints in y-direction (NY): %i\n", NY);
	fprintf(fp," Grid-spacing (DH): %e meter\n", DH);
	fprintf(fp," Model size in x-direction: %.5g m\n", NX*DH);
	fprintf(fp," Model size in y-direction: %.5g m\n", NY*DH);
	fprintf(fp," Time of wave propagation (T): %e seconds\n",TIME);
	fprintf(fp," Timestep (DT): %e seconds\n", DT);
	fprintf(fp," Number of timesteps: %i \n",NT);
	fprintf(fp,"\n");
	fprintf(fp," ------------------------- SOURCE -----------------------------\n");

	if (SRCREC){
		fprintf(fp," reading source positions, time delay, centre frequency \n");
		fprintf(fp," and initial amplitude from ASCII-file: %s\n\n",SOURCE_FILE);
	} else {
		fprintf(fp," plane wave excitation: depth= %5.2f meter \n",PLANE_WAVE_DEPTH);
		fprintf(fp," incidence angle of plane P-wave (from vertical) PLANE_WAVE_ANGLE= %5.2f degrees \n",PLANE_WAVE_ANGLE);
 		fprintf(fp," duration of source signal: %e seconds\n",TS);
 		fprintf(fp," (centre frequency is approximately %e Hz)\n",1.0/TS);
	}


	fprintf(fp," wavelet of source:");

	switch (SOURCE_SHAPE){
	case 1 :
		fprintf(fp," Ricker\n");
		break;
	case 2 :
		fprintf(fp," Fuchs-Mueller\n");
		break;
	case 3 :
		fprintf(fp," reading from \n\t %s\n",SIGNAL_FILE);
		break;
	case 4 :
		fprintf(fp," sinus raised to the power of 3.0 \n");
		break;
	default :
		declare_error(" Sorry, incorrect specification of source wavelet ! ");
		break;
	}

	fprintf(fp," Type of source:");
	switch (SOURCE_TYPE){
	case 1 :
		fprintf(fp," explosive source \n");
		break;
	case 2 :
		fprintf(fp," point source with directive force in x-direction\n");
		break;
	case 3 :
		fprintf(fp," point source with directive force in (vertical) y-direction\n");
		break;
	case 4 :
		fprintf(fp," point source with directive force in  z-direction\n");
		break;
	default :
		declare_error(" Sorry, wrong source type specification ! ");
		break;
	}
	
	fprintf(fp,"\n");

	if (SEISMO){
		fprintf(fp," ------------------------- RECEIVER  ------- -------------------\n");
		if (READREC){
			fprintf(fp," reading receiver positions from file \n");
			fprintf(fp,"\t%s\n\n",REC_FILE);
			fprintf(fp," reference point for receiver coordinate system:\n");
			fprintf(fp," x=%f \ty=%f\t z=%f\n",REFREC[1], REFREC[2], REFREC[3]);
		} else if (REC_ARRAY>0){
				fprintf(fp," horizontal lines of receivers.\n");
				fprintf(fp," number of lines: %d \n",REC_ARRAY);
				fprintf(fp," depth of upper line: %e m \n",REC_ARRAY_DEPTH);
				fprintf(fp," vertical increment between lines: %e m \n",REC_ARRAY_DIST);
				fprintf(fp," distance between receivers in x-direction within line: %i \n", DRX);		
		}else{

			fprintf(fp," first receiver position (XREC1,YREC1) = (%e, %e) m\n",
			    XREC1,YREC1);
			fprintf(fp," last receiver position (XREC2,YREC2) = (%e, %e) m\n",
			    XREC2,YREC2);
		}
		fprintf(fp,"\n");
	}

	fprintf(fp," ------------------------- FREE SURFACE ------------------------\n");
	if (FREE_SURF) fprintf(fp," free surface at the top of the model ! \n");
	else fprintf(fp," no free surface at the top of the model ! \n");
	fprintf(fp,"\n");

	fprintf(fp," ------------------------- ABSORBING FRAME ---------------------\n");
	if (FW>0){
		fprintf(fp," width of absorbing frame is %i grid points (%5.2f m).\n",FW,(float)FW*DH);
		if (ABS_TYPE==1) {
                	fprintf(fp," CPML damping applied. \n");
                	fprintf(fp," Damping velocity in the PML frame in m/s: %f .\n",VPPML);
                	fprintf(fp," Frequency within the PML frame in Hz: %f \n",FPML);
                	fprintf(fp," NPOWER: %f \n",NPOWER);
                	fprintf(fp," K_MAX: %f \n",K_MAX_CPML);
        	}

		
		if (ABS_TYPE==2) {
			fprintf(fp," Exponential damping applied. \n");
			fprintf(fp," Percentage of amplitude decay: %f .\n",DAMPING);
		}
	}
	else fprintf(fp," absorbing frame not installed ! \n");


	switch (BOUNDARY){
	case 0 :
		fprintf(fp," No periodic boundary condition.\n");
		break;
	case 1 :
		fprintf(fp," Periodic boundary condition at left and right edges.\n");
		break;
	default :
		warning(" Wrong integer value for BOUNDARY specified ! ");
		warning(" No periodic boundary condition will be applied ");
		BOUNDARY=0;
		break;
	}

	if (READMOD){
		fprintf(fp," ------------------------- MODEL-FILES -------------------------\n");
		fprintf(fp," names of model-files: \n");
		fprintf(fp,"\t shear wave velocities:\n\t %s.vs\n",MFILE);
		fprintf(fp,"\t tau for shear waves:\n\t %s.ts\n",MFILE);
		fprintf(fp,"\t density:\n\t %s.rho\n",MFILE);
		fprintf(fp,"\t compressional wave velocities:\n\t %s.vp\n",MFILE);
		fprintf(fp,"\t tau for P-waves:\n\t %s.tp\n",MFILE);
		for (l=1;l<=L;l++) fprintf(fp,"\t %1i. relaxation frequencies: %s.f%1i\n",l,MFILE,l);
	}

	fprintf(fp,"\n");
	fprintf(fp," ------------------------- Q-APROXIMATION --------------------\n");
	fprintf(fp," Number of relaxation mechanisms (L): %i\n",L);
	fprintf(fp," The L relaxation frequencies are at:  \n");
	for (l=1;l<=L;l++) fprintf(fp,"\t%f",FL[l]);
	fprintf(fp," Hz\n");
	fprintf(fp," Value for tau is : %f\n",TAU);


	if (SNAP){
		fprintf(fp,"\n");
		fprintf(fp," -----------------------  SNAPSHOTS  -----------------------\n");
		fprintf(fp," Snapshots of");
		switch(SNAP){
		case 1:
			fprintf(fp," x- and y-component");
			fprintf(fp," of particle velocity.\n");
			break;
		case 2:
			fprintf(fp," pressure field.\n");
			break;
		case 3:
			fprintf(fp," curl and divergence energy of the wavefield.\n");
			break;
		case 4:
			fprintf(fp," curl and divergence energy of the wavefield.\n");
			fprintf(fp," x- and y-component of particle velocity.\n");
			break;
		default:
			declare_error(" sorry, incorrect value for SNAP ! \n");
			break;
		}

		fprintf(fp," \t first (TSNAP1)= %8.5f s\n", TSNAP1);
		fprintf(fp," \t last (TSNAP2)=%8.5f s\n",TSNAP2);
		fprintf(fp," \t increment (TSNAPINC) =%8.5f s\n\n",TSNAPINC);
		fprintf(fp," \t first_and_last_horizontal(x)_gridpoint = %i, %i \n",1,NX);
		fprintf(fp," \t first_and_last_vertical_gridpoint = %i, %i \n",1,NY);
		fprintf(fp," \n name of output-file (SNAP_FILE):\n\t %s\n",SNAP_FILE);
		switch (SNAP_FORMAT){
		case 1 :
			declare_error(" SU-Format not yet available !!");
			break;
		case 2 :
			fprintf(fp," The data is written in ASCII. \n");
			break;
		case 3 :
			fprintf(fp," The data is written binary (IEEE) (4 byte per float)");
			break;
		default:
			declare_error(" Don't know the format for the Snapshot-data ! \n");
			break;
		}
	
		fprintf(fp,"\n\n");
	}
	if (SEISMO){
		fprintf(fp,"\n");
		fprintf(fp," -----------------------  SEISMOGRAMS  ----------------------\n");
		switch (SEIS_FORMAT){
			case 1: sprintf(file_ext,"su");  break;
			case 2: sprintf(file_ext,"txt"); break;
			case 3: sprintf(file_ext,"bin"); break;
		}
		if ((SEISMO==1) || (SEISMO==4)){
			fprintf(fp," Seismograms of ");
			fprintf(fp," x-, y-, and z-component");
			fprintf(fp," of particle velocity.\n");
			fprintf(fp," output-files: \n ");
			fprintf(fp,"\t%s_vx.%s\n\t%s_vy.%s\n",SEIS_FILE,file_ext,SEIS_FILE,file_ext);
		}
		if ((SEISMO==2) || (SEISMO==4)){
			fprintf(fp," Seismograms of pressure field (hydrophones).\n");
			fprintf(fp," output-file: \n ");
			fprintf(fp,"\t%s_p.%s\n",SEIS_FILE,file_ext);
		}
		if ((SEISMO==3) || (SEISMO==4)){
			fprintf(fp," Seismograms of curl (S-wave component) and div (P-wave component of wavefield).\n");
			fprintf(fp," output-files: \n ");
			fprintf(fp,"\t%s_rot.%s \n\t%s_div.%s\n",SEIS_FILE,file_ext,SEIS_FILE,file_ext);
			
		}		
			
		switch (SEIS_FORMAT){
		case 1 :
			fprintf(fp," The data is written in IEEE SU-format . \n");
			break;
		case 2 :
			fprintf(fp," The data is written in ASCII. \n");
			break;
		case 3 :
			fprintf(fp," The data is written binary IEEE (4 byte per float)");
			break;
		default:
			declare_error(" Sorry. I don't know the format for the seismic data ! \n");
			break;
		}
		fprintf(fp," samplingrate of seismic data: %f s\n",NDT*DT);
		if (!READREC) fprintf(fp," Trace-spacing: %5.2f m\n", NGEOPH*DH);
		fprintf(fp," Number of samples per trace: %i \n", iround(NT/NDT));
		fprintf(fp," ----------------------------------------------------------\n");
		fprintf(fp,"\n");
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n **********************************************************");
	fprintf(fp,"\n ******* PARAMETERS READ or PROCESSED within SOFI2D ********");
	fprintf(fp,"\n **********************************************************\n\n");


}
