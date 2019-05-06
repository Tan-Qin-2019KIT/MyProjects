/*---------------------------------------------------------------------------------
 * Copyright (C) 2015 For the list of authors, see file AUTHORS.
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
 ---------------------------------------------------------------------------------*/

/* $Id: wavefield_update_s_visc.c 819 2015-04-17 11:07:06Z tmetz $ */

/*Update Function of the stress-Wavefields in the viscoelastic case*/

#include "fd.h"

void wavefield_update_s_visc_4 ( int i, int j,float   vxx, float  vyx,float vxy,float  vyy, float **sxy,
                              float **sxx, float ** syy, float ***r, float ***p,
                              float ***q,float **fipjp, float **f, float **g, float *bip,
                              float *bjm,float *cip, float *cjm, float ***d, float ***e, float ***dip ,  float ** vxx_1,float ** vxx_2,float ** vxx_3,float ** vxx_4,float ** vyy_1,float ** vyy_2,float ** vyy_3,float ** vyy_4,float ** vxy_1,float ** vxy_2,float ** vxy_3,float ** vxy_4,float ** vyx_1,float ** vyx_2,float ** vyx_3,float ** vyx_4,float ** svx_1,float ** svx_2,float ** svx_3,float ** svx_4,float ** svy_1,float ** svy_2,float ** svy_3,float ** svy_4,float ***r_2,float ***r_3,float ***r_4, float ***p_2, float ***p_3, float ***p_4, float ***q_2, float ***q_3, float ***q_4)
{
    int l;
    float  dthalbe;
    extern float DT;
    extern int L;
    extern float DH;
    float sumr=0.0, sump=0.0, sumq=0.0;
    float c1, c2, c3, c4; /* Coefficients for Adam Bashforth */
    c1=13.0/12.0; c2=-5.0/24.0; c3=1.0/6.0; c4=-1.0/24.0;
    float sumxx=0.0,sumyy=0.0,sumxy=0.0,sumyx=0.0;
    float  dhi;
    float ctemp;
    
    dhi = 1.0/DH;
    dthalbe = DT/2.0;
    
    // Save derviations
    vxx_1[j][i]=vxx*DH;
    vyy_1[j][i]=vyy*DH;
    vxy_1[j][i]=vxy*DH;
    vyx_1[j][i]=vyx*DH;
    
    sumr=sump=sumq=0.0;
    for ( l=1; l<=L; l++ ) {
        sumr+=c1*r[j][i][l]+c2*r_2[j][i][l]+c3*r_3[j][i][l]+c4*r_4[j][i][l];
        sump+=c1*p[j][i][l]+c2*p_2[j][i][l]+c3*p_3[j][i][l]+c4*p_4[j][i][l];
        sumq+=c1*q[j][i][l]+c2*q_2[j][i][l]+c3*q_3[j][i][l]+c4*q_4[j][i][l];
    }
    
    // Calculate Adams-Bashforth stuff
    sumxx=c1*vxx_1[j][i]+c2*vxx_2[j][i]+c3*vxx_3[j][i]+c4*vxx_4[j][i];
    sumyy=c1*vyy_1[j][i]+c2*vyy_2[j][i]+c3*vyy_3[j][i]+c4*vyy_4[j][i];
    sumxy=c1*vxy_1[j][i]+c2*vxy_2[j][i]+c3*vxy_3[j][i]+c4*vxy_4[j][i];
    sumyx=c1*vyx_1[j][i]+c2*vyx_2[j][i]+c3*vyx_3[j][i]+c4*vyx_4[j][i];
    
    /* updating components of the stress tensor, partially */
    sxy[j][i] += ( fipjp[j][i]* ( sumxy+sumyx ) )*dhi + ( dthalbe*sumr );
    sxx[j][i] += ( g[j][i]* ( sumxx+sumyy )*dhi )- ( 2.0*f[j][i]*sumyy*dhi ) + ( dthalbe*sump );
    syy[j][i] += ( g[j][i]* (sumxx+sumyy)*dhi )- ( 2.0*f[j][i]*sumxx*dhi ) + ( dthalbe*sumq );
    
    sumr=sump=sumq=0.0;
    for ( l=1; l<=L; l++ ) {
        ctemp=2*(1-cip[l])/c1;
        r_4[j][i][l] = bip[l]* ( r[j][i][l]*cip[l]-ctemp*c2*(r[j][i][l]+r_2[j][i][l])-ctemp*c3*(r_2[j][i][l]+r_3[j][i][l])-ctemp*c4*(r_3[j][i][l]+r_4[j][i][l])- ( dip[j][i][l]* ( sumxy+sumyx)*dhi ) );
        p_4[j][i][l] = bjm[l]* ( p[j][i][l]*cjm[l]-ctemp*c2*(p[j][i][l]+p_2[j][i][l])-ctemp*c3*(p_2[j][i][l]+p_3[j][i][l])-ctemp*c4*(p_3[j][i][l]+p_4[j][i][l])- ( e[j][i][l]* ( sumxx+sumyy )*dhi ) + ( 2.0*d[j][i][l]*sumyy*dhi ) );
        q_4[j][i][l] = bjm[l]* ( q[j][i][l]*cjm[l]-ctemp*c2*(q[j][i][l]+q_2[j][i][l])-ctemp*c3*(q_2[j][i][l]+q_3[j][i][l])-ctemp*c4*(q_3[j][i][l]+q_4[j][i][l])- ( e[j][i][l]* ( sumxx+sumyy )*dhi ) + ( 2.0*d[j][i][l]*sumxx*dhi ) );
        sumr+=c1*r_4[j][i][l]+c2*r[j][i][l]+c3*r_2[j][i][l]+c4*r_3[j][i][l];
        sump+=c1*p_4[j][i][l]+c2*p[j][i][l]+c3*p_2[j][i][l]+c4*p_3[j][i][l];
        sumq+=c1*q_4[j][i][l]+c2*q[j][i][l]+c3*q_2[j][i][l]+c4*q_3[j][i][l];
    }
    
    sxy[j][i]+= ( dthalbe*sumr );
    sxx[j][i]+= ( dthalbe*sump );
    syy[j][i]+= ( dthalbe*sumq );
    
    
}