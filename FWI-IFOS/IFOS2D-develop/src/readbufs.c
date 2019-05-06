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
 *   read values of components of stress tensor at the edges of the
 *   local grid from buffer arrays (which have been exchanged between
 *   processes)
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void readbufs(float ** sxx, float ** syy, 
float ** sxy, float ** bufferlef_to_rig, float ** bufferrig_to_lef, 
float ** buffertop_to_bot, float ** bufferbot_to_top) {

	extern int NX, NY, POS[3], NPROCX, NPROCY, BOUNDARY;
	int i, j;


	if ((BOUNDARY) || (POS[1]!=NPROCX-1))	/* no boundary exchange at right edge of global grid */
	for (j=1;j<=NY;j++){

			sxx[j][NX+1] = bufferlef_to_rig[j][1];
			sxy[j][NX+1] = bufferlef_to_rig[j][2];
			sxx[j][NX+2] = bufferlef_to_rig[j][3];
		
	}

	if ((BOUNDARY) || (POS[1]!=0))	/* no boundary exchange at left edge of global grid */
	for (j=1;j<=NY;j++){
			sxx[j][0] = bufferrig_to_lef[j][1];
			sxy[j][0] = bufferrig_to_lef[j][2];
			sxy[j][-1] = bufferrig_to_lef[j][3];
	}

	if (POS[2]!=NPROCY-1)	/* no boundary exchange at bottom of global grid */
	for (i=1;i<=NX;i++){
			syy[NY+1][i] = buffertop_to_bot[i][1];
			sxy[NY+1][i] = buffertop_to_bot[i][2];
			syy[NY+2][i] = buffertop_to_bot[i][3];
	}

	if (POS[2]!=0)	/* no boundary exchange at top of global grid */
	for (i=1;i<=NX;i++){
			syy[0][i] = bufferbot_to_top[i][1];
			sxy[0][i] = bufferbot_to_top[i][2];
			sxy[-1][i] = bufferbot_to_top[i][3];
	}

}
