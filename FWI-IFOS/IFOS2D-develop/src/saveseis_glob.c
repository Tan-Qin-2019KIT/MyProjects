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
 *   write seismograms to files
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void saveseis_glob(FILE *fp, float **sectionvx, float **sectionvy,float **sectionvz,float **sectionp,float **sectioncurl, float **sectiondiv, int  **recpos, int  **recpos_loc,int ntr, float ** srcpos, int ishot, int ns, int iter, int type_switch){
    
    /* type_switch:
     *  1== synthetic data
     *  2== measured - synthetic data (residuals)
     *  3== filtered measured data
     */
    
    extern int SEISMO, SEIS_FORMAT, RUN_MULTIPLE_SHOTS, WAVETYPE, VERBOSE,FORWARD_ONLY;
    extern char SEIS_FILE[STRING_SIZE];
    extern int VELOCITY, WRITE_FILTERED_DATA, ADJOINT_TYPE;;
    
    char vxf[STRING_SIZE], vyf[STRING_SIZE],vzf[STRING_SIZE], curlf[STRING_SIZE], divf[STRING_SIZE], pf[STRING_SIZE];
    int nsrc=1;
    
    switch (type_switch) {
        case 1:
            if(ADJOINT_TYPE==1) {
                sprintf(vxf,"%s_vx.su.syn.shot%d.it%d",SEIS_FILE,ishot,iter);
                sprintf(vyf,"%s_vy.su.syn.shot%d.it%d",SEIS_FILE,ishot,iter);
            }
            if(ADJOINT_TYPE==2){
                sprintf(vyf,"%s_vy.su.syn.shot%d.it%d",SEIS_FILE,ishot,iter);
            }
            if(ADJOINT_TYPE==3){
                sprintf(vxf,"%s_vx.su.syn.shot%d.it%d",SEIS_FILE,ishot,iter);
            }

            if(WAVETYPE==2 || WAVETYPE==3) {
                sprintf(vzf,"%s_vz.su.syn.shot%d.it%d",SEIS_FILE,ishot,iter);
            }
            sprintf(pf,"%s_p.su.syn.shot%d.it%d",SEIS_FILE,ishot,iter);
            sprintf(divf,"%s_div.su.syn.shot%d.it%d",SEIS_FILE,ishot,iter);
            sprintf(curlf,"%s_curl.su.syn.shot%d.it%d",SEIS_FILE,ishot,iter);
            break;
            
        case 2:
            if(ADJOINT_TYPE==1) {
                sprintf(vxf,"%s_vx.su.shot%d_adjoint_src",SEIS_FILE,ishot);
                sprintf(vyf,"%s_vy.su.shot%d_adjoint_src",SEIS_FILE,ishot);
            }
            if(ADJOINT_TYPE==2){
                sprintf(vyf,"%s_vy.su.shot%d_adjoint_src",SEIS_FILE,ishot);
            }
            if(ADJOINT_TYPE==3){
                sprintf(vxf,"%s_vx.su.shot%d_adjoint_src",SEIS_FILE,ishot);
            }
            if(WAVETYPE==2 || WAVETYPE==3) {
                sprintf(vzf,"%s_vz.su.shot%d_adjoint_src",SEIS_FILE,ishot);
            }
            sprintf(pf,"%s_p.su.shot%d_adjoint_src",SEIS_FILE,ishot);
            sprintf(divf,"%s_div.su.shot%d_adjoint_src",SEIS_FILE,ishot);
            sprintf(curlf,"%s_curl.su.shot%d_adjoint_src",SEIS_FILE,ishot);
            break;
            
        case 3:
            if(WRITE_FILTERED_DATA==1){
                if(ADJOINT_TYPE==1) {
                    sprintf(vxf,"%s_vx.su.obs.shot%d.it%d",SEIS_FILE,ishot,iter);
                    sprintf(vyf,"%s_vy.su.obs.shot%d.it%d",SEIS_FILE,ishot,iter);
                }
                if(ADJOINT_TYPE==2){
                    sprintf(vyf,"%s_vy.su.obs.shot%d.it%d",SEIS_FILE,ishot,iter);
                }
                if(ADJOINT_TYPE==3){
                    sprintf(vxf,"%s_vx.su.obs.shot%d.it%d",SEIS_FILE,ishot,iter);
                }

                if(WAVETYPE==2 || WAVETYPE==3) {
                    sprintf(vzf,"%s_vz.su.obs.shot%d.it%d",SEIS_FILE,ishot,iter);
                }
                sprintf(pf,"%s_p.su.obs.shot%d.it%d",SEIS_FILE,ishot,iter);
                sprintf(divf,"%s_div.su.obs.shot%d.it%d",SEIS_FILE,ishot,iter);
                sprintf(curlf,"%s_curl.su.obs.shot%d.it%d",SEIS_FILE,ishot,iter);
            }else if(WRITE_FILTERED_DATA==2){
                if(ADJOINT_TYPE==1) {
                    sprintf(vxf,"%s_vx.su.obs.adj.shot%d.it%d",SEIS_FILE,ishot,iter);
                    sprintf(vyf,"%s_vy.su.obs.adj.shot%d.it%d",SEIS_FILE,ishot,iter);
                }
                if(ADJOINT_TYPE==2){
                    sprintf(vyf,"%s_vy.su.obs.adj.shot%d.it%d",SEIS_FILE,ishot,iter);
                }
                if(ADJOINT_TYPE==3){
                    sprintf(vxf,"%s_vx.su.obs.adj.shot%d.it%d",SEIS_FILE,ishot,iter);
                }

                if(WAVETYPE==2 || WAVETYPE==3) {
                    sprintf(vzf,"%s_vz.su.obs.adj.shot%d.it%d",SEIS_FILE,ishot,iter);
                }
                sprintf(pf,"%s_p.su.obs.adj.shot%d.it%d",SEIS_FILE,ishot,iter);
                sprintf(divf,"%s_div.su.obs.adj.shot%d.it%d",SEIS_FILE,ishot,iter);
                sprintf(curlf,"%s_curl.su.obs.adj.shot%d.it%d",SEIS_FILE,ishot,iter);
            }
            break;
            
        case 4:
            if(ADJOINT_TYPE==1) {
                sprintf(vxf,"%s_vx.su.syn.adj.shot%d.it%d",SEIS_FILE,ishot,iter);
                sprintf(vyf,"%s_vy.su.syn.adj.shot%d.it%d",SEIS_FILE,ishot,iter);
            }
            if(ADJOINT_TYPE==2){
                sprintf(vyf,"%s_vy.su.syn.adj.shot%d.it%d",SEIS_FILE,ishot,iter);
            }
            if(ADJOINT_TYPE==3){
                sprintf(vxf,"%s_vx.su.syn.adj.shot%d.it%d",SEIS_FILE,ishot,iter);
            }

            if(WAVETYPE==2 || WAVETYPE==3) {
                sprintf(vzf,"%s_vz.su.syn.adj.shot%d.it%d",SEIS_FILE,ishot,iter);
            }
            sprintf(pf,"%s_p.su.syn.adj.shot%d.it%d",SEIS_FILE,ishot,iter);
            sprintf(divf,"%s_div.su.syn.adj.shot%d.it%d",SEIS_FILE,ishot,iter);
            sprintf(curlf,"%s_curl.su.syn.adj.shot%d.it%d",SEIS_FILE,ishot,iter);
            break;
            
        default:
            declare_error("saveseis_glob: Unkown type_switch");
            break;
    }
    
    if(FORWARD_ONLY==1){
        sprintf(vxf,"%s_vx.su.shot%d",SEIS_FILE,ishot);
        sprintf(vyf,"%s_vy.su.shot%d",SEIS_FILE,ishot);
        if(WAVETYPE==2 || WAVETYPE==3) {
            sprintf(vzf,"%s_vz.su.shot%d",SEIS_FILE,ishot);
        }
        sprintf(pf,"%s_p.su.shot%d",SEIS_FILE,ishot);
        sprintf(divf,"%s_div.su.shot%d",SEIS_FILE,ishot);
        sprintf(curlf,"%s_curl.su.shot%d",SEIS_FILE,ishot);
    }
    
    switch (SEISMO){
        case 1 : /* particle velocities only */
            if (WAVETYPE==1 || WAVETYPE==3) {
                if(ADJOINT_TYPE==1 || ADJOINT_TYPE==0) {
                    
                    if(VERBOSE==1) fprintf(fp," PE %d is writing %d seismograms (vx) to\n\t %s \n",0,ntr,vxf);
                    outseis_glob(fp,fopen(vxf,"w"),1,sectionvx,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
                    if(VERBOSE==1) fprintf(fp," PE %d is writing %d seismograms (vy) to\n\t %s \n",0,ntr,vyf);
                    outseis_glob(fp,fopen(vyf,"w"),2,sectionvy,recpos,recpos_loc, ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
                }
                if(ADJOINT_TYPE==2){
                    if(VERBOSE==1) fprintf(fp," PE %d is writing %d seismograms (vy) to\n\t %s \n",0,ntr,vyf);
                    outseis_glob(fp,fopen(vyf,"w"),2,sectionvy,recpos,recpos_loc, ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
                }
                if(ADJOINT_TYPE==3){
                    if(VERBOSE==1) fprintf(fp," PE %d is writing %d seismograms (vx) to\n\t %s \n",0,ntr,vxf);
                    outseis_glob(fp,fopen(vxf,"w"),1,sectionvx,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
                }
            }
            if(WAVETYPE==2 || WAVETYPE==3) {
                if(VERBOSE==1) fprintf(fp," PE %d is writing %d seismograms (vz) to\n\t %s \n",0,ntr,vzf);
                outseis_glob(fp,fopen(vzf,"w"),2,sectionvz,recpos,recpos_loc, ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
            }
            break;
            
        case 2 : /* pressure only */
            if(VERBOSE==1) fprintf(fp," PE %d is writing %d seismograms of pressure to\n\t %s \n",0,ntr,pf);
            outseis_glob(fp,fopen(pf,"w"), 0, sectionp,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
            break;
            
        case 3 : /* curl and div only */
            if(VERBOSE==1) fprintf(fp," PE %d is writing %d seismograms of divergence to\n\t %s \n",0,ntr,divf);
            outseis_glob(fp,fopen(divf,"w"), 0, sectiondiv,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
            if(VERBOSE==1) fprintf(fp," PE %d is writing %d seismograms of curl to\n\t %s \n",0,ntr,curlf);
            outseis_glob(fp,fopen(curlf,"w"), 0, sectioncurl,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
            break;
            
        case 4 : /* everything */
            if (WAVETYPE==1 || WAVETYPE==3) {
                if(ADJOINT_TYPE==1 || ADJOINT_TYPE==0) {
                    if(VERBOSE==1) fprintf(fp," PE %d is writing %d seismograms (vx) to\n\t %s \n",0,ntr,vxf);
                    outseis_glob(fp,fopen(vxf,"w"),1,sectionvx,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
                    if(VERBOSE==1) fprintf(fp," PE %d is writing %d seismograms (vy) to\n\t %s \n",0,ntr,vyf);
                    outseis_glob(fp,fopen(vyf,"w"),2,sectionvy,recpos,recpos_loc, ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
                }
                if(ADJOINT_TYPE==2){
                    if(VERBOSE==1) fprintf(fp," PE %d is writing %d seismograms (vy) to\n\t %s \n",0,ntr,vyf);
                    outseis_glob(fp,fopen(vyf,"w"),2,sectionvy,recpos,recpos_loc, ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
                }
                if(ADJOINT_TYPE==3){
                    if(VERBOSE==1) fprintf(fp," PE %d is writing %d seismograms (vx) to\n\t %s \n",0,ntr,vxf);
                    outseis_glob(fp,fopen(vxf,"w"),1,sectionvx,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
                }

                if(VERBOSE==1) fprintf(fp," PE %d is writing %d seismograms of pressure to\n\t %s \n",0,ntr,pf);
                outseis_glob(fp,fopen(pf,"w"), 0, sectionp,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
                
                if(VERBOSE==1) fprintf(fp," PE %d is writing %d seismograms of divergence to\n\t %s \n",0,ntr,divf);
                outseis_glob(fp,fopen(divf,"w"),0, sectiondiv,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
                if(VERBOSE==1) fprintf(fp," PE %d is writing %d seismograms of curl to\n\t %s \n",0,ntr,curlf);
                outseis_glob(fp,fopen(curlf,"w"),0, sectioncurl,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
            }
            if(WAVETYPE==2 || WAVETYPE==3) {
                if(VERBOSE==1) fprintf(fp," PE %d is writing %d seismograms (vz) to\n\t %s \n",0,ntr,vzf);
                outseis_glob(fp,fopen(vzf,"w"),2,sectionvz,recpos,recpos_loc, ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
            }
            break;
            
        case 5 : /* everything except curl and div */
            if (WAVETYPE==1 || WAVETYPE==3) {
                if(ADJOINT_TYPE==1 || ADJOINT_TYPE==0) {
                    if(VERBOSE==1) fprintf(fp," PE %d is writing %d seismograms (vx) to\n\t %s \n",0,ntr,vxf);
                    outseis_glob(fp,fopen(vxf,"w"),1,sectionvx,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
                    if(VERBOSE==1) fprintf(fp," PE %d is writing %d seismograms (vy) to\n\t %s \n",0,ntr,vyf);
                    outseis_glob(fp,fopen(vyf,"w"),2,sectionvy,recpos,recpos_loc, ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
                }
                if(ADJOINT_TYPE==2){
                    if(VERBOSE==1) fprintf(fp," PE %d is writing %d seismograms (vy) to\n\t %s \n",0,ntr,vyf);
                    outseis_glob(fp,fopen(vyf,"w"),2,sectionvy,recpos,recpos_loc, ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
                }
                if(ADJOINT_TYPE==3){
                    if(VERBOSE==1) fprintf(fp," PE %d is writing %d seismograms (vx) to\n\t %s \n",0,ntr,vxf);
                    outseis_glob(fp,fopen(vxf,"w"),1,sectionvx,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
                }

                if(VERBOSE==1) fprintf(fp," PE %d is writing %d seismograms of pressure to\n\t %s \n",0,ntr,pf);
                outseis_glob(fp,fopen(pf,"w"), 0, sectionp,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
            }
            if(WAVETYPE==2 || WAVETYPE==3) {
                if(VERBOSE==1) fprintf(fp," PE %d is writing %d seismograms (vz) to\n\t %s \n",0,ntr,vzf);
                outseis_glob(fp,fopen(vzf,"w"),2,sectionvz,recpos,recpos_loc, ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
            }
            break;
            
    }
}
