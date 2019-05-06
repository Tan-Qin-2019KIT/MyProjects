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
 *   Write 2D snapshot for current timestep  to file
 *
 *  See COPYING file for copying and redistribution conditions.
 *  ----------------------------------------------------------------------*/

#include "fd.h"


void snap_SH(FILE *fp,int nt, int nsnap, float ** vz, float **u, float **pi, float *hc, int ishot){
    
    /*
     different data formats of output:
     SNAP_FORMAT=1  :  SU (IEEE)
     SNAP_FORMAT=2  :  ASCII
     SNAP_FORMAT=3  :  BINARY (IEEE)
     
     different types:
     SNAP=1 : values in vx and vy
     SNAP=2 : -(vx+vy) (pressure field)
     SNAP=3 : divergence of vx and vy (energy of compressional waves)
     and curl of vx and vy (energy of shear waves)
     SNAP=4 : both particle velocities (type=1) and energy (type=3)
     */
    
    
    int i,j, fdoh;
    float  dhi;
    char snapfile_z[STRING_SIZE];
    char  ext[8], wm[2];
    FILE *fpx1;
    
    extern float DH, DT;
    extern char SNAP_FILE[STRING_SIZE];
    extern int NX, NY,  SNAP_FORMAT, SNAP, FDORDER, ACOUSTIC;
    extern int MYID, POS[3], IDX, IDY,VERBOSE;
    
    /* Check if snapshots should be writen for this shot*/
    extern int SNAPSHOT_START, SNAPSHOT_END, SNAPSHOT_INCR;
    int check=0;
    for(i=SNAPSHOT_START;i<=SNAPSHOT_END;i+=SNAPSHOT_INCR) {
        if(ishot==i) check=1;
    }
    if(check==0) return;
    
    
    dhi = 1.0/DH;
    fdoh = FDORDER/2;
    
    switch(SNAP_FORMAT){
        case 1:
            sprintf(ext,".su");
            break;
        case 2:
            sprintf(ext,".asc");
            break;
        case 3:
            sprintf(ext,".bin");
            break;
    }
    
    sprintf(snapfile_z,"%s%s.z.shot%i.%i.%i",SNAP_FILE,ext,ishot,POS[1],POS[2]);

    
    if(VERBOSE) fprintf(fp,"\n\n PE %d is writing snapshot-data SH at T=%fs to \n",MYID,nt*DT);
    
    if (nsnap==1)
        sprintf(wm,"w");
    else
        sprintf(wm,"a");
    
    switch(SNAP){
        case 1 :
            if(VERBOSE) fprintf(fp,"%s\n", snapfile_z);
            fpx1=fopen(snapfile_z,wm);
            for (i=1;i<=NX;i+=IDX)
                for (j=1;j<=NY;j+=IDY){
                    writedsk(fpx1,vz[j][i],SNAP_FORMAT);
                }
            fclose(fpx1);
            break;
            
            
        case 2 :
            break;
            
        case 4 :
           if(VERBOSE)  fprintf(fp,"%s\n", snapfile_z);
            fpx1=fopen(snapfile_z,wm);
            for (i=1;i<=NX;i+=IDX)
                for (j=1;j<=NY;j+=IDY){
                    writedsk(fpx1,vz[j][i],SNAP_FORMAT);
                }
            fclose(fpx1);
            break;
        case 3 :
            break;
    }
    
    
}


