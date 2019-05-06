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
 *   local grid into buffer arrays (which will be exchanged between
 *   processes)
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void writebufv(float ** vx, float ** vy,
float ** bufferlef_to_rig, float ** bufferrig_to_lef, 
float ** buffertop_to_bot, float ** bufferbot_to_top){


	extern int NX, NY, POS[3], NPROCX, NPROCY, BOUNDARY;
	int i, j;


	/* exchange if periodic boundary condition is applied */
	if ((BOUNDARY) || (POS[1]!=0))	
	for (j=1;j<=NY;j++){


			/* storage of left edge of local volume into buffer */
			bufferlef_to_rig[j][1] =  vx[j][1];
			bufferlef_to_rig[j][2] =  vy[j][1];
			bufferlef_to_rig[j][3] =  vy[j][2];


	}


														/* no exchange if periodic boundary condition is applied */
	if ((BOUNDARY) || (POS[1]!=NPROCX-1))	/* no boundary exchange at right edge of global grid */
	for (j=1;j<=NY;j++){
			/* storage of right edge of local volume into buffer */
			bufferrig_to_lef[j][1] =  vx[j][NX];
			bufferrig_to_lef[j][2] =  vy[j][NX];
			bufferrig_to_lef[j][3] =  vx[j][NX-1];


	}

	if (POS[2]!=0)	/* no boundary exchange at top of global grid */
	for (i=1;i<=NX;i++){

			/* storage of top of local volume into buffer */
			buffertop_to_bot[i][1]  =  vx[1][i];
			buffertop_to_bot[i][2]  =  vy[1][i];
			buffertop_to_bot[i][3]  =  vx[2][i];
			
	}


	if (POS[2]!=NPROCY-1)	/* no boundary exchange at bottom of global grid */
	for (i=1;i<=NX;i++){
		
			/* storage of bottom of local volume into buffer */
			bufferbot_to_top[i][1]  =  vx[NY][i];
			bufferbot_to_top[i][2]  =  vy[NY][i];
			bufferbot_to_top[i][3]  =  vy[NY-1][i];

	}

}
