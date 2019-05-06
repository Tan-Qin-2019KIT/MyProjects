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
 *   merge snapshots files written by the different processes to
 *   a single file
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"


void merge(int nsnap, int type){
    
    
    extern char SNAP_FILE[STRING_SIZE];
    extern int NXG, NYG, SNAP_FORMAT, NPROCX, NPROCY;
    extern int NX, NY, IDX, IDY;
    extern FILE *FP;
    extern float DH;
    extern int SNAPSHOT_START, SNAPSHOT_END, SNAPSHOT_INCR;
    char file[STRING_SIZE], mfile[STRING_SIZE], outfile[STRING_SIZE], ext[10];
    FILE *fp[NPROCY][NPROCX], *fpout;
    int i, j, ip, jp, n;
    float a, max = 0.0;
    int shot;
    
    
    
    
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
    
    switch(type){
        case 1:
            fprintf(FP," x-component of particle velocity");
            strcat(ext,".x");
            break;
        case 2:
            fprintf(FP," y-component of particle velocity");
            strcat(ext,".y");
            break;
        case 4:
            fprintf(FP," P-wave energyfield");
            strcat(ext,".div");
            break;
        case 5:
            fprintf(FP," S-wave energyfield");
            strcat(ext,".rot");
            break;
        case 6:
            fprintf(FP," pressure");
            strcat(ext,".p");
            break;
        case 7:
            fprintf(FP," z-component of particle velocity");
            strcat(ext,".z");
            break;
        default:
            declare_error(" merge: cannot find snapfiles! ");
            break;
    }
    
    sprintf(mfile,"%s%s",SNAP_FILE,ext);
    fprintf(FP," (files: %s.??? ).\n",mfile);
    
    
    fprintf(FP," Opening snapshot files: %s.??? ",mfile);
    

    for(shot=SNAPSHOT_START;shot<=SNAPSHOT_END;shot+=SNAPSHOT_INCR) {
        sprintf(outfile,"%s.shot%i%s",SNAP_FILE,shot,ext);
        fprintf(FP,"\n writing merged snapshot file to  %s \n",outfile);
        fpout=fopen(outfile,"w");
        
        for (ip=0;ip<=NPROCX-1; ip++)
            for (jp=0;jp<=NPROCY-1; jp++){
                sprintf(file,"%s.shot%i.%i.%i",mfile,shot,ip,jp);
                fp[jp][ip]=fopen(file,"r");
                if (fp[jp][ip]==NULL) declare_error("merge: can't read snapfile !");
            }
        
        
        fprintf(FP," ... finished. \n");
        
        
        fprintf(FP," Copying...");
        
        max=0.0;
        for (n=0;n<=nsnap-2; n++){
            for (ip=0;ip<=NPROCX-1; ip++){
                for (i=1;i<=NX;i+=IDX){
                    for (jp=0;jp<=NPROCY-1; jp++){
                        for (j=1;j<=NY;j+=IDY){
                            a=readdsk(fp[jp][ip],SNAP_FORMAT);
                            if (a>max) max=a;
                            writedsk(fpout,a,SNAP_FORMAT);
                        }
                    }
                }
            }
        }
        fprintf(FP," ... finished. \n");
        
        for (ip=0;ip<=NPROCX-1; ip++)
            for (jp=0;jp<=NPROCY-1; jp++)
            {
                fclose(fp[jp][ip]);
            }
        for (ip=0;ip<=NPROCX-1; ip++){
            for (jp=0;jp<=NPROCY-1; jp++){
                sprintf(file,"%s.shot%i.%i.%i",mfile,shot,ip,jp);
                remove(file);
            }
        }
    }
    if (SNAP_FORMAT==3){
        fprintf(FP," Use \n");
        fprintf(FP," xmovie n1=%d n2=%d  < %s loop=1 label1=Y label2=X title=%%g d1=%f d2=%f f1=%f f2=%f clip=%e \n",
                ((NYG-1)/IDY)+1,((NXG-1)/IDY)+1,outfile,DH,DH,DH,DH,max/10.0);
        fprintf(FP," to play the movie. \n");
    }
    
    
}


