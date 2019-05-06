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

/*  Computation of local source coordinates
 *
 *
 */
#include "fd.h"
 
void alloc_sections(int ntr,int ns,float ***sectionvx,float ***sectionvy,float ***sectionvz,float ***sectionp,float ***sectionpnp1,float ***sectionpn,float ***sectioncurl,float ***sectiondiv,
                    float ***sectionpdata,float ***sectionpdiff,float ***sectionpdiffold,float ***sectionvxdata,float ***sectionvxdiff,float ***sectionvxdiffold,float ***sectionvydata,
                    float ***sectionvydiff,float ***sectionvydiffold,float ***sectionvzdata,float ***sectionvzdiff,float ***sectionvzdiffold)
{
	extern int SEISMO, WAVETYPE;
	extern FILE *FP;
	
	switch (SEISMO){
            case 1 : /* particle velocities only */
                switch (WAVETYPE) {
                    case 1:
                        *sectionvx=matrix(1,ntr,1,ns);
                        *sectionvy=matrix(1,ntr,1,ns);
                        break;
                    case 2:
                        *sectionvz=matrix(1,ntr,1,ns);
                        break;
                    case 3:
                        *sectionvx=matrix(1,ntr,1,ns);
                        *sectionvy=matrix(1,ntr,1,ns);
                        *sectionvz=matrix(1,ntr,1,ns);
                        break;
                }
                break;
            case 2 : /* pressure only */
                *sectionp=matrix(1,ntr,1,ns);
                *sectionpnp1=matrix(1,ntr,1,ns);
                *sectionpn=matrix(1,ntr,1,ns);
                break;
            case 3 : /* curl and div only */
                *sectioncurl=matrix(1,ntr,1,ns);
                *sectiondiv=matrix(1,ntr,1,ns);
                break;
            case 4 : /* everything */
                switch (WAVETYPE) {
                    case 1:
                        *sectionvx=matrix(1,ntr,1,ns);
                        *sectionvy=matrix(1,ntr,1,ns);
                        break;
                    case 2:
                        *sectionvz=matrix(1,ntr,1,ns);
                        break;
                    case 3:
                        *sectionvx=matrix(1,ntr,1,ns);
                        *sectionvy=matrix(1,ntr,1,ns);
                        *sectionvz=matrix(1,ntr,1,ns);
                        break;
                }
                *sectioncurl=matrix(1,ntr,1,ns);
                *sectiondiv=matrix(1,ntr,1,ns);
                *sectionp=matrix(1,ntr,1,ns);
                break;
            case 5 : /* everything except curl and div*/
                switch (WAVETYPE) {
                    case 1:
                        *sectionvx=matrix(1,ntr,1,ns);
                        *sectionvy=matrix(1,ntr,1,ns);
                        break;
                    case 2:
                        *sectionvz=matrix(1,ntr,1,ns);
                        break;
                    case 3:
                        *sectionvx=matrix(1,ntr,1,ns);
                        *sectionvy=matrix(1,ntr,1,ns);
                        *sectionvz=matrix(1,ntr,1,ns);
                        break;
                }
                *sectionp=matrix(1,ntr,1,ns);
                break;
        }
        
        *sectionpdata=matrix(1,ntr,1,ns);
    *sectionpdiff=matrix(1,ntr,1,ns);
    *sectionpdiffold=matrix(1,ntr,1,ns);
    switch (WAVETYPE) {
        case 1:
            *sectionvxdata=matrix(1,ntr,1,ns);
            *sectionvxdiff=matrix(1,ntr,1,ns);
            *sectionvxdiffold=matrix(1,ntr,1,ns);
            *sectionvydata=matrix(1,ntr,1,ns);
            *sectionvydiff=matrix(1,ntr,1,ns);
            *sectionvydiffold=matrix(1,ntr,1,ns);
            break;
            
        case 2:
            *sectionvzdata=matrix(1,ntr,1,ns);
            *sectionvzdiff=matrix(1,ntr,1,ns);
            *sectionvzdiffold=matrix(1,ntr,1,ns);
            break;
            
        case 3:
            *sectionvxdata=matrix(1,ntr,1,ns);
            *sectionvxdiff=matrix(1,ntr,1,ns);
            *sectionvxdiffold=matrix(1,ntr,1,ns);
            *sectionvydata=matrix(1,ntr,1,ns);
            *sectionvydiff=matrix(1,ntr,1,ns);
            *sectionvydiffold=matrix(1,ntr,1,ns);
            *sectionvzdata=matrix(1,ntr,1,ns);
            *sectionvzdiff=matrix(1,ntr,1,ns);
            *sectionvzdiffold=matrix(1,ntr,1,ns);
            break;
    }
 }
 
 
 
 void dealloc_sections(int ntr,int ns,int **recpos_loc,float **sectionvx,float **sectionvy,float **sectionvz,float **sectionp,float **sectionpnp1,float **sectionpn,float **sectioncurl,float **sectiondiv,
                    float **sectionpdata,float **sectionpdiff,float **sectionpdiffold,float **sectionvxdata,float **sectionvxdiff,float **sectionvxdiffold,float **sectionvydata,
                    float **sectionvydiff,float **sectionvydiffold,float **sectionvzdata,float **sectionvzdiff,float **sectionvzdiffold)
{
	extern int SEISMO, WAVETYPE;
	extern FILE *FP;
	
	free_imatrix(recpos_loc,1,3,1,ntr);
        switch (SEISMO){
            case 1 : /* particle velocities only */
                if (WAVETYPE==1 || WAVETYPE==3) {
                    free_matrix(sectionvx,1,ntr,1,ns);
                    free_matrix(sectionvy,1,ntr,1,ns);
                }
                if (WAVETYPE==2 || WAVETYPE==3) {
                    free_matrix(sectionvz,1,ntr,1,ns);
                }
                break;
            case 2 : /* pressure only */
                if (WAVETYPE==1 || WAVETYPE==3) {
                    free_matrix(sectionp,1,ntr,1,ns);
                    free_matrix(sectionpn,1,ntr,1,ns);
                    free_matrix(sectionpnp1,1,ntr,1,ns);
                }
                break;
            case 3 : /* curl and div only */
                if (WAVETYPE==1 || WAVETYPE==3) {
                    free_matrix(sectioncurl,1,ntr,1,ns);
                    free_matrix(sectiondiv,1,ntr,1,ns);
                }
                break;
            case 4 : /* everything */
                if (WAVETYPE==1 || WAVETYPE==3) {
                    free_matrix(sectionvx,1,ntr,1,ns);
                    free_matrix(sectionvy,1,ntr,1,ns);
                    free_matrix(sectionp,1,ntr,1,ns);
                    free_matrix(sectioncurl,1,ntr,1,ns);
                    free_matrix(sectiondiv,1,ntr,1,ns);
                }
                if (WAVETYPE==2 || WAVETYPE==3) {
                    free_matrix(sectionvz,1,ntr,1,ns);
                }
                break;
            case 5 : /* everything except curl and div */
                if (WAVETYPE==1 || WAVETYPE==3) {
                    free_matrix(sectionvx,1,ntr,1,ns);
                    free_matrix(sectionvy,1,ntr,1,ns);
                    free_matrix(sectionp,1,ntr,1,ns);
                }
                if (WAVETYPE==2 || WAVETYPE==3) {
                    free_matrix(sectionvz,1,ntr,1,ns);
                }
                break;
        }
    
    
    if (WAVETYPE==1 || WAVETYPE==3) {
        free_matrix(sectionvxdata,1,ntr,1,ns);
        free_matrix(sectionvxdiff,1,ntr,1,ns);
        free_matrix(sectionvydata,1,ntr,1,ns);
        free_matrix(sectionvydiff,1,ntr,1,ns);
        free_matrix(sectionvydiffold,1,ntr,1,ns);
        free_matrix(sectionvxdiffold,1,ntr,1,ns);
        free_matrix(sectionpdata,1,ntr,1,ns);
        free_matrix(sectionpdiff,1,ntr,1,ns);
        free_matrix(sectionpdiffold,1,ntr,1,ns);
    }
    if (WAVETYPE==2 || WAVETYPE==3) {
        free_matrix(sectionvzdata,1,ntr,1,ns);
        free_matrix(sectionvzdiff,1,ntr,1,ns);
        free_matrix(sectionvzdiffold,1,ntr,1,ns);
    }
 }