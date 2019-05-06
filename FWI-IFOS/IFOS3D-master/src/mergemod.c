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

/*------------------------------------------------------------------------
 *   merge model files written by the different processes to 
 *   a single file                                 
 *   last update 10/10/00   T. Bohlen
 *  ----------------------------------------------------------------------*/

#include "fd.h"


void mergemod(char modfile[STRING_SIZE], int format){



	extern int MYID, NPROCX, NPROCY, NPROCZ;
	extern int NX, NY, NZ, NPROC, IDX, IDY, IDZ, VERBOSE;
	extern FILE *FP;


	char file[STRING_SIZE];
	FILE *fp[NPROCY_MAX][NPROCX_MAX][NPROCZ_MAX], *fpout;
	int i, j, k, ip, jp, kp;
	float a;
	

	if ((NPROCX>NPROCX_MAX)||(NPROCY>NPROCY_MAX)||(NPROCZ>NPROCZ_MAX))
		err(" mergemod.c: constant expression NPROC?_MAX < NPROC? ");

	
	if (VERBOSE) printf(" PE %d starts merge of %d model files \n",MYID,NPROC);

	fprintf(FP," writing merged model file to  %s \n",modfile);
	fpout=fopen(modfile,"w");


	for (kp=0;kp<=NPROCZ-1; kp++)
	    for (ip=0;ip<=NPROCX-1; ip++)
   		for (jp=0;jp<=NPROCY-1; jp++){
      		sprintf(file,"%s.%i.%i.%i",modfile,ip,jp,kp);
      		fp[jp][ip][kp]=fopen(file,"r");
      		if (fp[jp][ip][kp]==NULL) err("mergemod.c: can't read modfile !"); 
     	 }

	//fprintf(FP," ... finished. \n");



	fprintf(FP," Copying...");

	for (kp=0;kp<=NPROCZ-1; kp++)
	for (k=1;k<=NZ;k+=IDZ)
	for (ip=0;ip<=NPROCX-1; ip++)
	for (i=1;i<=NX;i+=IDX)
	for (jp=0;jp<=NPROCY-1; jp++)
	for (j=1;j<=NY;j+=IDY){
	
	      	a=readdsk(fp[jp][ip][kp],format);
	      	writedsk(fpout,a,format);
	}

	fprintf(FP," ... finished. \n");

	for (kp=0;kp<=NPROCZ-1; kp++)
	for (ip=0;ip<=NPROCX-1; ip++)
   	for (jp=0;jp<=NPROCY-1; jp++){
      		fclose(fp[jp][ip][kp]);
      	}
      	
	fclose(fpout);
		
	/*fprintf(FP," Use \n");
	fprintf(FP," xmovie n1=%d n2=%d < %s loop=1 label1=Y label2=X title=%%g \n",
		((NYG-1)/IDY)+1,((NXG-1)/IDY)+1,modfile);
	fprintf(FP," to visualize model. \n");*/


}


