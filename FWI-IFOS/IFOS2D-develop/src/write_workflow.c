/*-----------------------------------------------------------------------------------------
 * Copyright (C) 2015  For the list of authors, see file AUTHORS.
 *
 * This file is part of DENISE.
 *
 * DENISE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 *
 * DENISE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with DENISE. See file COPYING and/or <http://www.gnu.org/licenses/gpl-2.0.html>.
 -----------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
 *  Write workflow after each iteration
 *
 *  Written by Wittkamp Oct 2015
 *
 *  If you want to adjust the workflow, I think it is not necessary to
 *  modify this file. Have a look at apply_workflow.c and adjust
 *  WORKFLOW_MAX_VAR in fd.h.
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void write_workflow(float ** workflow, int iter, int workflow_lines){
    
    /* local variables */
    int x,y;
    int  workflow_line_next;
    FILE *fwork;
    char workflow_out[STRING_SIZE];
    
    /* extern variables */
    extern char FILE_WORKFLOW[STRING_SIZE];
    extern FILE *FP;
    extern int VERBOSE;
    /* WORKFLOW_MAX_VAR is set in fd.h */
    
    sprintf(workflow_out,"%s.restart",FILE_WORKFLOW);
    
    if(VERBOSE) fprintf(FP,"\n Writing workflow to file %s\n",workflow_out);
    
    /* Open Workflow file */
    fwork=fopen(workflow_out,"w+");
    if (fwork==NULL) declare_error(" Workflow output file could no be opened !");
    
    workflow_line_next=1;
    for(y=1;y<workflow_lines;y++){
        if(iter>workflow[y][1]){
            workflow_line_next=y;
        }
    }
    
    fprintf(fwork,"%i\n",iter);
    /* write workflow */
    for(y=workflow_line_next;y<=workflow_lines;y++){
        for(x=1;x<=WORKFLOW_MAX_VAR;x++){
            if(x>1)fprintf(fwork,"\t");
            fprintf(fwork,"%f",workflow[y][x]);
        }
        fprintf(fwork,"\n");
    }
    fclose(fwork);
}
