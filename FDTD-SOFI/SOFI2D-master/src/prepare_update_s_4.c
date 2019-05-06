/*------------------------------------------------------------------------
 * Copyright (C) 2011 For the list of authors, see file AUTHORS.
 *
 * This file is part of SOFI2D.
 *
 * SOFI2D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 *
 * SOFI2D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SOFI2D. See file COPYING and/or
 * <http://www.gnu.org/licenses/gpl-2.0.html>.
 --------------------------------------------------------------------------*/
/* ------------------------------------------------------------------------
 * prepare of update of the stress tensor
 * ------------------------------------------------------------------------*/


/* ------------------------------------------------------------------------
 * ATTENTION: The parameters below will be scaled by factor c* due to
 * Adams-Bashforth method, so be aware and only call this function when
 * FDORDER_TIME is set to 4
 * ------------------------------------------------------------------------*/

#include "fd.h"

void prepare_update_s_4(float *etajm, float *etaip, float *peta, float **fipjp, float **pu,
                      float **puipjp, float **ppi, float **ptaus, float **ptaup,
                      float **ptausipjp, float **f, float **g, float *bip, float *bjm,
                      float *cip, float *cjm, float ***dip, float ***d, float ***e) {
    
    extern int NX, NY, L;
    extern float DT;
    int i, j, l;
    float c1; /* Coefficients for Adam Bashforth */
    c1=13.0/12.0;
    
    for (l=1;l<=L;l++){
        etajm[l] = peta[l];
        etaip[l] = peta[l];
    }
    for (j=1;j<=NY;j++){
        for (i=1;i<=NX;i++){
            fipjp[j][i] = puipjp[j][i]*DT*(1.0+L*ptausipjp[j][i]);
            f[j][i] = pu[j][i]*DT*(1.0+L*ptaus[j][i]);
            g[j][i] = ppi[j][i]*DT*(1.0+L*ptaup[j][i]);
            for (l=1;l<=L;l++){
                bip[l] = 1.0/(1.0+(c1*etaip[l]*0.5));
                bjm[l] = 1.0/(1.0+(c1*etajm[l]*0.5));
                cip[l] = 1.0-(c1*etaip[l]*0.5);
                cjm[l] = 1.0-(c1*etajm[l]*0.5);
                dip[j][i][l] = puipjp[j][i]*etaip[l]*ptausipjp[j][i];
                d[j][i][l] = pu[j][i]*etajm[l]*ptaus[j][i];
                e[j][i][l] = ppi[j][i]*etajm[l]*ptaup[j][i];
            }
        }
    }
    
}