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

/*
   Holberg coefficients for a certain FD order and a margin of error E
   (MAXRELERROR)

   MAXRELERROR = 0 -> Taylor-coeff.
   MAXRELERROR = 1 -> Holberg-coeff.: E = 0.1 %
   MAXRELERROR = 2 ->                 E = 0.5 %
   MAXRELERROR = 3 ->                 E = 1.0 %
   MAXRELERROR = 4 ->                 E = 3.0 %
 
  hc: column 0 = minimum number of grid points per shortest wavelength
      columns 1-6 = Holberg coefficients
     
*/


#include "fd.h"           /* general include file for viscoelastic FD programs */

float *holbergcoeff() {

	int nl,i;
	extern int FDORDER, MAXRELERROR, MYID;
	float *hc;
	
	float hcall[5][6][7] =
	{
	    {	
		{ 23.0, 1.0,                0.0,              0.0,                0.0,              0.0,              0.0            }, 
		{  8.0, 9.0/8.0,           -1.0/24.0,         0.0,                0.0,              0.0,              0.0            },
		{  6.0, 75.0/64.0,         -25.0/384.0,       3.0/640.0,          0.0,              0.0,              0.0            },
		{  5.0, 1225.0/1024.0,     -245.0/3072.0,     49.0/5120.0,       -5.0/7168.0,       0.0,              0.0            },
		{  5.0, 19845.0/16384.0,   -735.0/8192.0,     567.0/40960.0,     -405.0/229376.0,   35.0/294912.0,    0.0            },
		{  4.0, 160083.0/131072.0, -12705.0/131072.0, 22869.0/1310720.0, -5445.0/1835008.0, 847.0/2359296.0, -63.0/2883584.0 }
	    },
	    {	
		{ 49.7 , 1.0010,  0.0,      0.0,        0.0,       0.0,        0.0       }, 
		{  8.32, 1.1382, -0.046414, 0.0,        0.0,       0.0,        0.0       },
		{  4.77, 1.1965, -0.078804, 0.0081781,  0.0,       0.0,        0.0       },
		{  3.69, 1.2257, -0.099537, 0.018063,  -0.0026274, 0.0,        0.0       },
		{  3.19, 1.2415, -0.11231,  0.026191,  -0.0064682, 0.001191,   0.0       },
		{  2.91, 1.2508, -0.12034,  0.032131,  -0.010142,  0.0029857, -0.00066667}
	    },
	    {	
		{ 22.2 , 1.0050,  0.0,      0.0,        0.0,       0.0,        0.0       }, 
		{  5.65, 1.1534, -0.052806, 0.0,        0.0,       0.0,        0.0       },
		{  3.74, 1.2111, -0.088313, 0.011768,   0.0,       0.0,        0.0       },
		{  3.11, 1.2367, -0.10815,  0.023113,  -0.0046905, 0.0,        0.0       },
		{  2.80, 1.2496, -0.11921,  0.031130,  -0.0093272, 0.0025161,  0.0       },
		{  2.62, 1.2568, -0.12573,  0.036423,  -0.013132,  0.0047484, -0.0015979 }
	    },
	    {
		{ 15.8,  1.0100,  0.0,      0.0,        0.0,       0.0,        0.0       }, 
		{  4.80, 1.1640, -0.057991, 0.0,        0.0,       0.0,        0.0       },
		{  3.39, 1.2192, -0.094070, 0.014608,   0.0,       0.0,        0.0       },
		{  2.90, 1.2422, -0.11269,  0.026140,  -0.0064054, 0.0,        0.0       },
		{  2.65, 1.2534, -0.12257,  0.033755,  -0.011081,  0.0036784,  0.0       },
		{  2.51, 1.2596, -0.12825,  0.038550,  -0.014763,  0.0058619, -0.0024538 }
	    },
	    {
		{  9.16, 1.0300,  0.0,      0.0,        0.0,       0.0,        0.0       }, 
		{  3.47, 1.1876, -0.072518, 0.0,        0.0,       0.0,        0.0       },
		{  2.91, 1.2341, -0.10569,  0.022589,   0.0,       0.0,        0.0       },
		{  2.61, 1.2516, -0.12085,  0.032236,  -0.011459,  0.0,        0.0       },
		{  2.45, 1.2596, -0.12829,  0.038533,  -0.014681,  0.0072580,  0.0       },
		{  2.36, 1.2640, -0.13239,  0.042217,  -0.017803,  0.0081959, -0.0051848 }
	    },
	};
	
	if (MYID == 0) {
		if ((FDORDER!=2) && (FDORDER!=4) && (FDORDER!=6) && (FDORDER!=8) && (FDORDER!=10) && (FDORDER!=12)) {
			declare_error(" Error in selection of FD coefficients: wrong FDORDER! ");
		}

		if ((MAXRELERROR<0) || (MAXRELERROR>4)) {
			declare_error(" Error in selection of FD coefficients: wrong choice of maximum relative error! ");
		}
	}

	nl = (FDORDER/2) - 1;
	
	hc = vector(0,6);

	for (i=0; i<=6; i++) {
		hc[i] = hcall[MAXRELERROR][nl][i];
	}
	
	
	return hc;

}
