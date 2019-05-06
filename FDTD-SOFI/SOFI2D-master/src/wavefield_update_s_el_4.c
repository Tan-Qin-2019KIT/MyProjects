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

/*Update Function of the stress-Wavefields in the elastic case*/

#include "fd.h"

void wavefield_update_s_el_4 ( int i, int j,float   vxx, float  vyx,float vxy,float  vyy, float **sxy,
                            float **sxx, float ** syy, float ** pi, float ** u, float ** uipjp,float ** vxx_1,float ** vxx_2,float ** vxx_3,float ** vxx_4,float ** vyy_1,float ** vyy_2,float ** vyy_3,float ** vyy_4,float ** vxy_1,float ** vxy_2,float ** vxy_3,float ** vxy_4,float ** vyx_1,float ** vyx_2,float ** vyx_3,float ** vyx_4)
{
    
    extern float DT;
    float fipjp, f, g;
    extern float DH;
    float  dhi;
    float c1, c2, c3, c4; /* Coefficients for Adam Bashforth */
    c1=13.0/12.0; c2=-5.0/24.0; c3=1.0/6.0; c4=-1.0/24.0;
    float sumxx=0.0,sumyy=0.0,sumxy=0.0,sumyx=0.0;
    
    dhi = 1.0/DH;
    fipjp=uipjp[j][i]*DT;
    f = u[j][i]*DT;
    g = pi[j][i]*DT;

    // Save derviations
    vxx_1[j][i]=vxx*DH;
    vyy_1[j][i]=vyy*DH;
    vxy_1[j][i]=vxy*DH;
    vyx_1[j][i]=vyx*DH;

    // Calculate Adams-Bashforth stuff
    sumxx=c1*vxx_1[j][i]+c2*vxx_2[j][i]+c3*vxx_3[j][i]+c4*vxx_4[j][i];
    sumyy=c1*vyy_1[j][i]+c2*vyy_2[j][i]+c3*vyy_3[j][i]+c4*vyy_4[j][i];
    sumxy=c1*vxy_1[j][i]+c2*vxy_2[j][i]+c3*vxy_3[j][i]+c4*vxy_4[j][i];
    sumyx=c1*vyx_1[j][i]+c2*vyx_2[j][i]+c3*vyx_3[j][i]+c4*vyx_4[j][i];
    
    // Update stress
    sxy[j][i]+= (( fipjp* ( sumxy+sumyx ) ))*dhi;
    sxx[j][i]+= (( g* ( sumxx+sumyy ) )- ( 2.0*f*sumyy ))*dhi;
    syy[j][i]+= (( g* ( sumxx+sumyy ) )- ( 2.0*f*sumxx ))*dhi;
    
}