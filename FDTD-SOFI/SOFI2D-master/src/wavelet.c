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
/*------------------------------------------------------------------------
*   Calculating source signal at different source positions with different
*   time-shift, centre frequency and amplitude (as specified in SOURCE_FILE).
*   Source signals are written to array signals 
*
*  ----------------------------------------------------------------------*/

#include "fd.h"


float ** wavelet(float ** srcpos_loc, int nsrc){


	/* extern variables */
	extern int SOURCE_SHAPE, NT, MYID;
	extern float  DT;
	extern char SIGNAL_FILE[STRING_SIZE];
	extern FILE *FP;

	/*local variables */
	int nts=0, nt, k;
	float *psource=NULL, tshift, amp=0.0, a, fc, tau, t, ts;
	/*float a2;*/ /*different Ricker signal used for comparison with analytical solution*/
	float ** signals;


	
	if (SOURCE_SHAPE==3) psource=rd_sour(&nts,fopen(SIGNAL_FILE,"r"));
		
	signals=matrix(1,nsrc,1,NT);
	
	
	for (nt=1;nt<=NT;nt++){
			t=(float)nt*DT;
			
			for (k=1;k<=nsrc;k++) {
				tshift=srcpos_loc[4][k];
				fc=srcpos_loc[5][k];
				a=srcpos_loc[6][k];
				ts=1.0/fc;

				switch (SOURCE_SHAPE){
					case 1 : 
						/* standard Ricker signal */
						tau=PI*(t-1.5*ts-tshift)/(ts);
						amp=(((1.0-2.0*tau*tau)*exp(-tau*tau)));
						
						/* different Ricker signal used for comparison with analytical
						solution*/
						/*tau=t-ts-tshift;
						a2=2.0*(PI*fc*PI*fc);
						amp=tau*(a2*tau*tau-1.5)*exp(-a2*tau*tau);*/
						
					break;
					case 2 : 
						if ((t<tshift) || (t>(tshift+ts))) amp=0.0;
						else amp=((sin(2.0*PI*(t-tshift)*fc) 
			    				-0.5*sin(4.0*PI*(t-tshift)*fc)));

/*						amp=((sin(2.0*PI*(t+tshift)*fc) 
			    				-0.5*sin(4.0*PI*(t+tshift)*fc)));
*/
					break;
					case 3 : 
						if (nt<=nts) amp=psource[nt]; 
						else amp=0.0;
					break;  /* source wavelet from file SOURCE_FILE */
					case 4 : 
						if ((t<tshift) || (t>(tshift+ts))) amp=0.0;
						else amp=pow(sin(PI*(t-tshift)/ts),3.0);
						break; /* sinus raised to the power of three */
					default : 
						declare_error("No valid source-wavelet (SOURCE_SHAPE) specified! ");
						break;
					}
					
					
					signals[k][nt]=amp*a;
		}
	}
	
	fprintf(FP," Message from function wavelet written by PE %d \n",MYID);
	fprintf(FP," %d source positions located in subdomain of PE %d \n",nsrc,MYID);
	fprintf(FP," have been assigned with a source signal. \n");
			
		
	if (SOURCE_SHAPE==3) free_vector(psource,1,nts);
	
	return signals;	

}
