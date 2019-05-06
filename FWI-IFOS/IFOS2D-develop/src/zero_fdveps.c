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
 *   zero wavefield
 *
 *
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void zero_fdveps(int ny1, int ny2, int nx1, int nx2, float ** vx, float ** vy, float ** vz, float ** sxx, float ** syy, float ** sxy,float ** sxz,float ** syz,float ** vxm1, float ** vym1, float ** uxy, float ** vxp1, float ** vyp1,float ** psi_sxx_x, float ** psi_sxy_x, float ** psi_sxz_x, float ** psi_vxx, float ** psi_vyx, float ** psi_vzx, float ** psi_syy_y, float ** psi_sxy_y, float ** psi_syz_y, float ** psi_vyy, float ** psi_vxy,float ** psi_vzy,float ** psi_vxxs){
    
    
    
    register int i, j;
    extern int FW, NX, NY,WAVETYPE;
    
    
    for (j=ny1;j<=ny2;j++){
        for (i=nx1;i<=nx2;i++){
            if (WAVETYPE==1 || WAVETYPE==3) {
                vx[j][i]=0.0;
                vy[j][i]=0.0;
                sxx[j][i]=0.0;
                syy[j][i]=0.0;
                sxy[j][i]=0.0;
                vxm1[j][i]=0.0;
                vym1[j][i]=0.0;
                uxy[j][i]=0.0;
                vxp1[j][i]=0.0;
                vyp1[j][i]=0.0;
            }
            if (WAVETYPE==2 || WAVETYPE==3) {
                vz[j][i]=0.0;
                sxz[j][i]=0.0;
                syz[j][i]=0.0;
            }
        }
    }
    
    for (j=1;j<=NY;j++){
        for (i=1;i<=2*FW;i++){
            if (WAVETYPE==1 || WAVETYPE==3) {
                psi_sxx_x[j][i] = 0.0;
                psi_sxy_x[j][i] = 0.0;
                psi_vxx[j][i] = 0.0;
                psi_vxxs[j][i] = 0.0;
                psi_vyx[j][i] = 0.0;
            }
            if (WAVETYPE==2 || WAVETYPE==3) {
                psi_sxz_x[j][i] = 0.0;
                psi_vzx[j][i] = 0.0;
            }
        }
    }
    
    for (j=1;j<=2*FW;j++){
        for (i=1;i<=NX;i++){
            if (WAVETYPE==1 || WAVETYPE==3) {
                psi_syy_y[j][i] = 0.0;
                psi_sxy_y[j][i] = 0.0;
                psi_vyy[j][i] = 0.0;
                psi_vxy[j][i] = 0.0;
            }
            if (WAVETYPE==2 || WAVETYPE==3) {
                psi_syz_y[j][i] = 0.0;
                psi_vzy[j][i] = 0.0;
            }
        }
    }
    
}
