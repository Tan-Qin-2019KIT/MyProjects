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
 *   Read source wavelet in su format                                  
 *  ----------------------------------------------------------------------*/
#include "fd.h"
#include "segy.h"

void  inseis_source_wavelet(float *section, int ns, int ishot, int SH, int STF){

	/* declaration of extern variables */
	extern int MYID;
	extern char SIGNAL_FILE[STRING_SIZE];
    extern char SIGNAL_FILE_SH[STRING_SIZE];
    
    extern int USE_WORKFLOW;
    extern int WORKFLOW_STAGE;
	/* declaration of local variables */
	int j;
	float dump;
	segy tr;
	char data[STRING_SIZE];
    FILE *fpdata;
    
    if (STF==0){ /* reading inverted signals */
        if(SH==0) {
            if(USE_WORKFLOW){
                sprintf(data,"%s.stage%d.shot%d.su",SIGNAL_FILE,WORKFLOW_STAGE,ishot);
            } else {
                sprintf(data,"%s.shot%d.su",SIGNAL_FILE,ishot);
            }
        } else {
            if(USE_WORKFLOW){
                sprintf(data,"%s.stage%d.shot%d.su",SIGNAL_FILE_SH,WORKFLOW_STAGE,ishot);
            } else {
                sprintf(data,"%s.shot%d.su",SIGNAL_FILE_SH,ishot);
            }
        }
    } else { /* reading signals for STF inversion */
        if(SH==0) {
            sprintf(data,"%s.shot%d_start.su",SIGNAL_FILE,ishot);
        } else {
            sprintf(data,"%s.shot%d_start.su",SIGNAL_FILE_SH,ishot);
        }
    }
    
	fpdata = fopen(data,"r");
	if (fpdata==NULL) declare_error(" Source wavelet not found ");

	/* SEGY (without file-header) */
	fread(&tr,240,1,fpdata);
	fread(&tr.data,4,ns,fpdata);
			
	for(j=0;j<ns;j++){
	    dump=tr.data[j];
	    section[j+1]=dump;
	}
			  
	fclose(fpdata);
}
