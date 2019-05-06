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
 *   write values of dynamic field variables at the edges of the
 *   local grid into buffer arrays and  exchanged between
 *   processes.
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void exchange_v(float ** vx, float ** vy, float ** vz,
                float ** bufferlef_to_rig, float ** bufferrig_to_lef,
                float ** buffertop_to_bot, float ** bufferbot_to_top,
                MPI_Request * req_send, MPI_Request * req_rec, int wavetyp_start){
    
    
    extern int NX, NY, POS[3], NPROCX, NPROCY, BOUNDARY, FDORDER, WAVETYPE;
    extern int INDEX[5];
    extern const int TAG1,TAG2,TAG5,TAG6;
    MPI_Status  status;
    int i, j, fdo, fdo3, n, l;
    
    fdo = FDORDER/2 + 1;

    switch (wavetyp_start) {
        case 1:
            fdo3 = 2*fdo;
            break;
        case 2:
            fdo3 = 1*fdo;
            break;
        case 3:
            fdo3 = 3*fdo;
            break;
        default:
            fdo3 = 2*fdo;
            break;
    }
    
    /* ************ */
    /* top - bottom */
    /* ************ */
    
    if (POS[2]!=0)	/* no boundary exchange at top of global grid */
        for (i=1;i<=NX;i++){
            n = 1;
            /* storage of top of local volume into buffer */
            if (WAVETYPE==1 || WAVETYPE==3) {
                for (l=1;l<=fdo-1;l++) {
                    buffertop_to_bot[i][n++]  =  vy[l][i];
                }
                for (l=1;l<=fdo;l++) {
                    buffertop_to_bot[i][n++]  =  vx[l][i];
                }
            }
            if(WAVETYPE==2 || WAVETYPE==3) {
                for (l=1;l<=fdo;l++) {
                    buffertop_to_bot[i][n++]  =  vz[l][i];
                }
            }
        }
    
    if (POS[2]!=NPROCY-1)	/* no boundary exchange at bottom of global grid */
        for (i=1;i<=NX;i++){
            /* storage of bottom of local volume into buffer */
            n = 1;
            if (WAVETYPE==1 || WAVETYPE==3) {
                for (l=1;l<=fdo;l++) {
                    bufferbot_to_top[i][n++]  =  vy[NY-l+1][i];
                }
                for (l=1;l<=fdo-1;l++) {
                    bufferbot_to_top[i][n++]  =  vx[NY-l+1][i];
                }
            }
            if(WAVETYPE==2 || WAVETYPE==3) {
                for (l=1;l<=fdo-1;l++) {
                    bufferbot_to_top[i][n++]  =  vz[NY-l+1][i];
                }
            }
        }
    
    /* still blocking communication */
    MPI_Sendrecv_replace(&buffertop_to_bot[1][1],NX*fdo3,MPI_FLOAT,INDEX[3],TAG5,INDEX[4],TAG5,MPI_COMM_WORLD,&status);
    MPI_Sendrecv_replace(&bufferbot_to_top[1][1],NX*fdo3,MPI_FLOAT,INDEX[4],TAG6,INDEX[3],TAG6,MPI_COMM_WORLD,&status);
    
    if (POS[2]!=NPROCY-1)	/* no boundary exchange at bottom of global grid */
        for (i=1;i<=NX;i++){
            n = 1;
            if (WAVETYPE==1 || WAVETYPE==3) {
                for (l=1;l<=fdo-1;l++) {
                    vy[NY+l][i] = buffertop_to_bot[i][n++];
                }
                for (l=1;l<=fdo;l++) {
                    vx[NY+l][i] = buffertop_to_bot[i][n++];
                }
            }
            if(WAVETYPE==2 || WAVETYPE==3) {
                for (l=1;l<=fdo;l++) {
                    vz[NY+l][i] = buffertop_to_bot[i][n++];
                }
            }
        }
    
    if (POS[2]!=0)	/* no boundary exchange at top of global grid */
        for (i=1;i<=NX;i++){
            n = 1;
            if (WAVETYPE==1 || WAVETYPE==3) {
                for (l=1;l<=fdo;l++) {
                    vy[1-l][i] = bufferbot_to_top[i][n++];
                }
                for (l=1;l<=fdo-1;l++) {
                    vx[1-l][i] = bufferbot_to_top[i][n++];
                }
            }
            if(WAVETYPE==2 || WAVETYPE==3) {
                for (l=1;l<=fdo-1;l++) {
                    vz[1-l][i] = bufferbot_to_top[i][n++];
                }
            }
        }
    
    
    /* ************ */
    /* left - right */
    /* ************ */
    
    /* exchange if periodic boundary condition is applied */
    if ((BOUNDARY) || (POS[1]!=0))
        for (j=1;j<=NY;j++){
            /* storage of left edge of local volume into buffer */
            n = 1;
            if (WAVETYPE==1 || WAVETYPE==3) {
                for (l=1;l<fdo;l++) {
                    bufferlef_to_rig[j][n++] =  vy[j][l];
                }
                for (l=1;l<fdo-1;l++) {
                    bufferlef_to_rig[j][n++] =  vx[j][l];
                }
            }
            if(WAVETYPE==2 || WAVETYPE==3) {
                for (l=1;l<fdo;l++) {
                    bufferlef_to_rig[j][n++] =  vz[j][l];
                }
            }
        }
    /* no exchange if periodic boundary condition is applied */
    if ((BOUNDARY) || (POS[1]!=NPROCX-1))	/* no boundary exchange at right edge of global grid */
        for (j=1;j<=NY;j++){
            /* storage of right edge of local volume into buffer */
            n = 1;
            if (WAVETYPE==1 || WAVETYPE==3) {
                for (l=1;l<fdo-1;l++) {
                    bufferrig_to_lef[j][n++] =  vy[j][NX-l+1];
                }
                for (l=1;l<fdo;l++) {
                    bufferrig_to_lef[j][n++] =  vx[j][NX-l+1];
                }
            }
            if(WAVETYPE==2 || WAVETYPE==3) {
                for (l=1;l<fdo-1;l++) {
                    bufferrig_to_lef[j][n++] =  vz[j][NX-l+1];
                }
            }
        }
    
    /* send and reveive values for points at inner boundaries */
    
    /* alternative communication */
    /* still blocking communication */
    MPI_Sendrecv_replace(&bufferlef_to_rig[1][1],NY*fdo3,MPI_FLOAT,INDEX[1],TAG1,INDEX[2],TAG1,MPI_COMM_WORLD,&status);
    MPI_Sendrecv_replace(&bufferrig_to_lef[1][1],NY*fdo3,MPI_FLOAT,INDEX[2],TAG2,INDEX[1],TAG2,MPI_COMM_WORLD,&status);
    
    
    /* no exchange if periodic boundary condition is applied */
    if ((BOUNDARY) || (POS[1]!=NPROCX-1))	/* no boundary exchange at right edge of global grid */
        for (j=1;j<=NY;j++){
            n = 1;
            if (WAVETYPE==1 || WAVETYPE==3) {
                for (l=1;l<fdo;l++) {
                    vy[j][NX+l] = bufferlef_to_rig[j][n++];
                }
                for (l=1;l<fdo-1;l++) {
                    vx[j][NX+l] = bufferlef_to_rig[j][n++];
                }
            }
            if(WAVETYPE==2 || WAVETYPE==3) {
                for (l=1;l<fdo;l++) {
                    vz[j][NX+l] = bufferlef_to_rig[j][n++];
                }
            }
        }
    
    /* no exchange if periodic boundary condition is applied */
    if ((BOUNDARY) || (POS[1]!=0))	/* no boundary exchange at left edge of global grid */
        for (j=1;j<=NY;j++){
            n = 1;
            if (WAVETYPE==1 || WAVETYPE==3) {
                for (l=1;l<fdo-1;l++) {
                    vy[j][1-l] = bufferrig_to_lef[j][n++];
                }
                for (l=1;l<fdo;l++) {
                    vx[j][1-l] = bufferrig_to_lef[j][n++];
                }
            }
            if(WAVETYPE==2 || WAVETYPE==3) {
                for (l=1;l<fdo-1;l++) {
                    vz[j][1-l] = bufferrig_to_lef[j][n++];
                }
            }
        }
    
    

}

