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

#include "fd.h"

void create_trkill_table(int ** killtable, int ntr_glob, int **recpos, int nsrc_glob, float **srcpos, int ishot, float kill_offset_lower, float kill_offset_upper) {
    
    /* Local variables */
    int s, r;
    float offset=0.0;
    
    /* extern variables */
    extern float DH;
    extern int MYID;
    extern char TRKILL_FILE[STRING_SIZE];
    
    /*---------------------*/
    /* Generate Killtable  */
    /*---------------------*/
    for(s=1;s<=nsrc_glob;s++){
        for (r=1; r<=ntr_glob; r++) {
            
            /* calculate Offset for current reciever and source */
            offset=fabs(sqrt((srcpos[1][s]-recpos[1][r]*DH)*(srcpos[1][s]-recpos[1][r]*DH)+(srcpos[2][s]-recpos[2][r]*DH)*(srcpos[2][s]-recpos[2][r]*DH)));
            
            /* If offset lies inside the bound kill */
            if( offset>kill_offset_lower && offset<kill_offset_upper ) {
                killtable[r][s]=1;
            }
        }
    }
    
    /*---------------------*/
    /*       Debug         */
    /*---------------------*/
    /* Debug is declared with ishot=-100 */
    if((ishot==-100) && (MYID==0)){
        FILE *FP_TK_OUT;
        char filename[225];
        
        sprintf(filename,"%stracekill.out.%.0f_%.0f.txt",TRKILL_FILE,kill_offset_lower,kill_offset_upper);
        FP_TK_OUT=fopen(filename,"w");
        
        for (r=1; r<=ntr_glob; r++) {
            for(s=1;s<=nsrc_glob;s++){
                fprintf(FP_TK_OUT,"%i\t",killtable[r][s]);
            }
            fprintf(FP_TK_OUT,"\n");
        }
        fclose(FP_TK_OUT);
    }
    
}
