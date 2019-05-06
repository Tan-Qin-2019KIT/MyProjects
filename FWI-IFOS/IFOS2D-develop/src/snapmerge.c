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
 *   loop over snapshotfiles which have to be merged.
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "globvar.h"      /* definition of global variables  */

int main(int argc, char **argv){
    
    int nsnap;
    char *fileinp;
    
    /* read parameters from parameter-file (stdin) */
    /*read_par(stdin);*/
    
    
    /* read parameters from parameter-file (stdin) */
    fileinp=argv[1];
    FP=fopen(fileinp,"r");
    if(FP==NULL) {
        if (MYID == 0){
            printf("\n==================================================================\n");
            printf(" Cannot open IFOS input file %s \n",fileinp);
            printf("\n==================================================================\n\n");
            declare_error(" --- ");
        }
    }
    
    /* read json formatted input file */
    read_par_json(stdout,fileinp);
    
    
    
    NXG=NX;
    NYG=NY;
    NX = NXG/NPROCX;
    NY = NYG/NPROCY;
    
    nsnap=1+iround((TSNAP2-TSNAP1)/TSNAPINC);
    
    FP=stdout;
    
    
    switch(SNAP){
        case 1 : /*particle velocity*/
            if(WAVETYPE==1||WAVETYPE==3) {
                merge(nsnap,1);
                merge(nsnap,2);
            }
            if(WAVETYPE==2||WAVETYPE==3) {
                merge(nsnap,7);
            }
            break;
        case 2 : /*pressure */
            merge(nsnap,6);
            break;
        case 4 : /*particle velocity*/
            if(WAVETYPE==1||WAVETYPE==3) {
                merge(nsnap,1);
                merge(nsnap,2);
            }
            if(WAVETYPE==2||WAVETYPE==3) {
                merge(nsnap,7);
            }
        case 3 :
            merge(nsnap,4);
            merge(nsnap,5);
            break;
        default :
            warning(" snapmerge: cannot identify content of snapshot !");
            break;
            
    }	
    return 0;	
    
}
