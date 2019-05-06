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
 *   output source signal e.g. for cross-correlation of comparison with analytical solutions                                  
 *  ----------------------------------------------------------------------*/
#include "fd.h"
#include "segy.h"

void  output_source_signal(FILE *fp, float **signals, int ns, int seis_form){

	/* declaration of extern variables */
	extern float DT;
	extern int NDT;
	
	/* declaration of local variables */
	int i,j, ntr=1;
	segy tr;
	int tracl ;
	
	
	
	switch(seis_form){
	case 1 :
		for(tracl=1;tracl<=ntr;tracl++){        /*SEGY (without file-header)*/
         		
			tr.tracl=tracl;     
			tr.trid=(short)1;           /* trace identification code: 1=seismic*/
			
			tr.ns=(unsigned short)iround(ns/NDT); /* number of samples in this trace */
			tr.dt=(unsigned short)iround(((float)NDT*DT)*1.0e6); /* sample interval in micro-seconds */
			tr.d1=(float)NDT*DT;        /* sample spacing for non-seismic data */
			
			tr.scalel=(signed short)1;
			tr.scalco=(signed short)1;
			tr.sx=(signed int)0;
			tr.sy=(signed int)0;
			tr.gx=(signed int)0;
			tr.gy=(signed int)0;
			tr.selev=0;
			tr.gelev=0;
			tr.sdepth=0;
			tr.gdel=0;
			tr.sdel=0;
			tr.gwdep=0;
			tr.offset=(signed int)0;
			tr.swdep=0;
			
			tr.ep=0;
			tr.cdp=(int)0;
			tr.tracr=0;
			tr.fldr=0;
			tr.tracf=0;
			tr.ep=0;
			tr.cdpt=0;
			tr.nvs=0;
			tr.nhs=0;
			tr.duse=0;
			tr.counit=0;
			tr.wevel=0;
			tr.swevel=0;
			tr.sut=0;
			tr.gut=0;
			tr.sstat=0;
			tr.gstat=0;
			tr.tstat=0;
			tr.laga=0;
			tr.lagb=0;
			tr.delrt=0;
			tr.muts=0;
			tr.mute=0;
			tr.gain=0;
			tr.igc=0;
			tr.igi=0;
			tr.corr=0;
			tr.sfs=0;
			tr.sfe=0;
			tr.slen=0;
			tr.styp=0;
			tr.stas=0;
			tr.stae=0;
			tr.tatyp=0;
			tr.afilf=0;
			tr.afils=0;
			tr.nofilf=0;
			tr.nofils=0;
			tr.lcf=0;
			tr.hcf=0;
			tr.lcs=0;
			tr.hcs=0;
			tr.year=0;
			tr.day=0;
			tr.hour=0;
			tr.minute=0;
			tr.sec=0;
			tr.timbas=0;
			tr.trwf=0;
			tr.grnors=0;
			tr.grnofr=0;
			tr.grnlof=0;
			tr.gaps=0;
			tr.otrav=0;
			tr.f1=0.0;
			tr.d2=0.0;
			tr.f2=0.0;
			tr.ungpow=0.0;
			tr.unscale=0.0;
			tr.ntr=0;
			tr.mark=0;

			for(j=1;j<=tr.ns;j++) tr.data[j]=signals[tracl][j*NDT];

			fwrite(&tr,240,1,fp);
			fwrite(&tr.data[1],4,tr.ns,fp);
		}
		break;


	case 2 :
		for(i=1;i<=ntr;i++){         /*ASCII ONE COLUMN*/
			for(j=1;j<=ns;j+=NDT) fprintf(fp,"%e\n", signals[i][j]);
		}
		break;

	case 3 :                             /*BINARY */

		for(i=1;i<=ntr;i++)
			for(j=1;j<=ns;j+=NDT){
				fwrite(&signals[i][j],sizeof(float),1,fp); }
		break;

	default :
		fprintf(stdout," Message from output_source_signal: Don't know data format for seismograms !\n");
		fprintf(stdout," No output written. ");
	}

	fclose(fp);
}
