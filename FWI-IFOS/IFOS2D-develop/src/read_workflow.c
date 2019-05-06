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
 *  Read extern workflow                          
 *
 *
 *  If you want to adjust the workflow, I think it is not necessary to
 *  modify this file. Have a look at apply_workflow.c and adjust
 *  WORKFLOW_MAX_VAR in fd.h.
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void read_workflow(char file_in[STRING_SIZE],float *** workflow, int *workflow_lines, char header[STRING_SIZE]){
    
    /* workflow is a pointer to a pointer, keep care... */
    
    /* intern variables */
    int i, c;
    int nw=0,x,y;
    
    FILE *fwork;
    /* extern variables */
    extern FILE *FP;
    /* WORKFLOW_MAX_VAR is set in fd.h */
    
    fprintf(FP,"\n Reading Workflow from file: %s\n",file_in);
    
    /* Open Workflow file */
    fwork=fopen(file_in,"r");
    if (fwork==NULL) declare_error(" Workflow file could no be opened !");
    
    /* Count how many lines the work flow file has */
    fgets(header, 200, fwork); /* Read header */
    while(!feof(fwork)) {
        c = fgetc(fwork);
        if(c == '\n'){nw++;}
    }
    fseek(fwork, SEEK_SET, 0); /* Reset */
    fgets(header, 200, fwork); /* Read header */
    
    fprintf(FP," Number of lines in workflow file: %i\n",nw);
    
    /* Allocate memory for workflow */
    (*workflow) = matrix(1,nw,1,WORKFLOW_MAX_VAR);
    if ((*workflow)==NULL) declare_error(" Was not able to allocate memory for workflow file !");
    
    /* Read workflow */
    fprintf(FP,"\n %s",header);
    fprintf(FP," ");
    for(y=1;y<=nw;y++){
        for(x=1;x<=WORKFLOW_MAX_VAR;x++){
            (*workflow)[y][x]=0;
            fscanf(fwork,"%f",&(*workflow)[y][x]);
            fprintf(FP,"%.3f\t",(*workflow)[y][x]);
        }
    fprintf(FP,"\n ");
    }
    fclose(fwork);
    *workflow_lines=nw;
}
