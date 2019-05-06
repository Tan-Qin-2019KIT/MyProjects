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
/* ------------------------------------------------------------------------
 * reads checkpoint file to continue a previously started modeling
 * (to allow for splitting of modeling jobs)
 *
 * ------------------------------------------------------------------------*/

#include "fd.h"

void read_checkpoint(int nx1, int nx2, int ny1, int ny2,
float **  vx, float ** vy, float ** sxx, float ** syy, float ** sxy){

	int i,j;
	char myid[5];
	FILE *fp;
	char checkptfile[STRING_SIZE];
	extern int MYID;
	extern char  CHECKPTFILE[STRING_SIZE];



	sprintf(checkptfile,"%s",CHECKPTFILE);
	sprintf(myid,".%d",MYID);
	strcat(checkptfile,myid);



	fp=fopen(checkptfile,"rb");
	if (fp==NULL) declare_error("CHECKPTFILE can't be opened !");
	
	for (j=ny1;j<=ny2;j++){
		for (i=nx1;i<=nx2;i++){

		fread( &vx[j][i],sizeof(float),1,fp);
		fread( &vy[j][i],sizeof(float),1,fp);
		fread(&sxx[j][i],sizeof(float),1,fp);
		fread(&syy[j][i],sizeof(float),1,fp);
		fread(&sxy[j][i],sizeof(float),1,fp);
		}
	}

	fclose(fp);

}
