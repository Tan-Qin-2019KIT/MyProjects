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
 *   updating particle velocities at gridpoints [nx1...nx2][ny1...ny2]
 *   by a staggered grid finite difference scheme of FDORDER accuracy in space
 *   and second order accuracy in time
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"



void update_v_acoustic_PML(int nx1, int nx2, int ny1, int ny2, int nt,
	float **  vx, float **  vxp1, float **  vxm1, float ** vy, float **  vyp1, float **  vym1, float ** sp,
	float  **rip, float **rjp, float **  srcpos_loc, float ** signals, float ** signals1, int nsrc, float ** absorb_coeff,
	float *hc, int infoout,int sw, float * K_x_half, float * a_x_half, float * b_x_half,
        float * K_y_half, float * a_y_half, float * b_y_half,
        float ** psi_sxx_x, float ** psi_syy_y){

	int i, j,l,fdoh,m, h, h1;
	float amp, dtdh, azi_rad;
	float vxtmp, vytmp;
        float sp_x, sp_y;
	
	extern float DT, DH;
	double time1, time2;
	extern int MYID, SOURCE_TYPE, ADJOINT_TYPE, FDORDER;
        extern int FDORDER, PARAMETERIZATION;
        extern int FREE_SURF, BOUNDARY, FW;
        extern int NPROCX, NPROCY, POS[3];
	extern FILE *FP;
	extern int VELOCITY;

	
	fdoh = FDORDER/2;
	dtdh = DT*DT/DH;
        
     /* drad = PI/180.0; 
        angle = 135.0; */
         
	if (infoout && (MYID==0)){
		time1=MPI_Wtime();
		fprintf(FP,"\n **Message from update_v (printed by PE %d):\n",MYID);
		fprintf(FP," Updating particle velocities ...");
	}

	
			
	/* ------------------------------------------------------------
	 * Important!
	 * rip and rjp are reciprocal values of averaged densities
	 * ------------------------------------------------------------ */

	switch (FDORDER){
	
	case 2:
		for (j=ny1;j<=ny2;j++){
		for (i=nx1;i<=nx2;i++){
			sp_x = hc[1]*(sp[j][i+1]-sp[j][i]);
			sp_y = hc[1]*(sp[j+1][i]-sp[j][i]);
			

			/* left boundary */                                         
			if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
				
				psi_sxx_x[j][i] = b_x_half[i] * psi_sxx_x[j][i] + a_x_half[i] * sp_x;                                                
				sp_x = sp_x / K_x_half[i] + psi_sxx_x[j][i];
				
			}

			/* right boundary */                                         
			if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){

				h1 = (i-nx2+2*FW);
				h=i; 
				
				/*psi_sxx_x[j][i] = b_x_half[i] * psi_sxx_x[j][i] + a_x_half[i] * sxx_x;                                                
				sxx_x = sxx_x / K_x_half[i] + psi_sxx_x[j][i];*/
				
				psi_sxx_x[j][h1] = b_x_half[h1] * psi_sxx_x[j][h1] + a_x_half[h1] * sp_x;                                                
				sp_x = sp_x / K_x_half[h1] + psi_sxx_x[j][h1];
			}
	
			/* top boundary */                                         
			if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
				
				psi_syy_y[j][i] = b_y_half[j] * psi_syy_y[j][i] + a_y_half[j] * sp_y;                                                
				sp_y = sp_y / K_y_half[j] + psi_syy_y[j][i];
			}
	
			/* bottom boundary */                                         
			if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){
			
				h1 = (j-ny2+2*FW);
				h = j;
							
				/*psi_syy_y[j][i] = b_y_half[j] * psi_syy_y[j][i] + a_y_half[j] * syy_y;                                                
				syy_y = syy_y / K_y_half[j] + psi_syy_y[j][i];*/
				
				psi_syy_y[h1][i] = b_y_half[h1] * psi_syy_y[h1][i] + a_y_half[h1] * sp_y;                                                
				sp_y = sp_y / K_y_half[h1] + psi_syy_y[h1][i];
			}                       
			
			if(sw==0){
				if (VELOCITY==0){
                    vxp1[j][i] = rip[j][i]*(sp_x)/DH;
                    vyp1[j][i] = rjp[j][i]*(sp_y)/DH;
                }
				else{
                    vxp1[j][i] = rip[j][i]*(sp_x)/DH;
                    vyp1[j][i] = rjp[j][i]*(sp_y)/DH;                    
                }
			}
			
			if(sw==1){
				if (VELOCITY==0){
                    vxp1[j][i] += vx[j][i]*DT;
                    vyp1[j][i] += vy[j][i]*DT;
                }
				else{
                    vxp1[j][i] = vx[j][i];
                    vyp1[j][i] = vy[j][i];
                }
			}
			
			vx[j][i] += DT*rip[j][i]*(sp_x)/DH;
			vy[j][i] += DT*rjp[j][i]*(sp_y)/DH;
      
		}
		}
	break;
	
		
	case 4:
		for (j=ny1;j<=ny2;j++){
		for (i=nx1;i<=nx2;i++){
			sp_x =  hc[1]*(sp[j][i+1]-sp[j][i])
					+ hc[2]*(sp[j][i+2]-sp[j][i-1]);

			sp_y = hc[1]*(sp[j+1][i]-sp[j][i])
					+ hc[2]*(sp[j+2][i]-sp[j-1][i]);
			
			/* left boundary */                                         
			if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
				
				psi_sxx_x[j][i] = b_x_half[i] * psi_sxx_x[j][i] + a_x_half[i] * sp_x;                                                
				sp_x = sp_x / K_x_half[i] + psi_sxx_x[j][i];

			}

			/* right boundary */                                         
			if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){

				h1 = (i-nx2+2*FW);
				h=i; 
					
				/*psi_sxx_x[j][i] = b_x_half[i] * psi_sxx_x[j][i] + a_x_half[i] * sxx_x;                                                
				sxx_x = sxx_x / K_x_half[i] + psi_sxx_x[j][i];*/
				
				psi_sxx_x[j][h1] = b_x_half[h1] * psi_sxx_x[j][h1] + a_x_half[h1] * sp_x;                                                
				sp_x = sp_x / K_x_half[h1] + psi_sxx_x[j][h1];

			}

	
			/* top boundary */                                         
			if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){

				psi_syy_y[j][i] = b_y_half[j] * psi_syy_y[j][i] + a_y_half[j] * sp_y;                                                
				sp_y = sp_y / K_y_half[j] + psi_syy_y[j][i]; 
				
			}
		

			/* bottom boundary */                                         
			if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){
				
				h1 = (j-ny2+2*FW);
				h = j;
							
				/*psi_syy_y[j][i] = b_y_half[j] * psi_syy_y[j][i] + a_y_half[j] * syy_y;                                                
				syy_y = syy_y / K_y_half[j] + psi_syy_y[j][i];*/
				
				psi_syy_y[h1][i] = b_y_half[h1] * psi_syy_y[h1][i] + a_y_half[h1] * sp_y;                                                
				sp_y = sp_y / K_y_half[h1] + psi_syy_y[h1][i];
				
			}                   
			
			if(sw==0){
				if (VELOCITY==0){
					vxp1[j][i] = rip[j][i]*(sp_x)/DH;
					vyp1[j][i] = rjp[j][i]*(sp_y)/DH;}
				else{
					vxp1[j][i] = rip[j][i]*(sp_x)/DH;
					vyp1[j][i] = rjp[j][i]*(sp_y)/DH;}
			}
			
			if(sw==1){
				if (VELOCITY==0){
					vxp1[j][i] += vx[j][i]*DT;
					vyp1[j][i] += vy[j][i]*DT;}
				else{
					vxp1[j][i] = vx[j][i];
					vyp1[j][i] = vy[j][i];}
			}
			
			vx[j][i] += DT*rip[j][i]*(sp_x)/DH;
			vy[j][i] += DT*rjp[j][i]*(sp_y)/DH;

		}
		}
	
	break;
		
	case 6:
		for (j=ny1;j<=ny2;j++){
		for (i=nx1;i<=nx2;i++){
			sp_x =  hc[1]*(sp[j][i+1]-sp[j][i])
					+ hc[2]*(sp[j][i+2]-sp[j][i-1])
					+ hc[3]*(sp[j][i+3]-sp[j][i-2]);
			
			sp_y = hc[1]*(sp[j+1][i]-sp[j][i])
					+ hc[2]*(sp[j+2][i]-sp[j-1][i])
					+ hc[3]*(sp[j+3][i]-sp[j-2][i]);
			

			/* left boundary */                                         
			if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
					
				psi_sxx_x[j][i] = b_x_half[i] * psi_sxx_x[j][i] + a_x_half[i] * sp_x;                                                
				sp_x = sp_x / K_x_half[i] + psi_sxx_x[j][i];
			}

			/* right boundary */                                         
			if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){

				h1 = (i-nx2+2*FW);
				h=i; 
					
				/*psi_sxx_x[j][i] = b_x_half[i] * psi_sxx_x[j][i] + a_x_half[i] * sxx_x;                                                
				sxx_x = sxx_x / K_x_half[i] + psi_sxx_x[j][i];*/
				
				psi_sxx_x[j][h1] = b_x_half[h1] * psi_sxx_x[j][h1] + a_x_half[h1] * sp_x;                                                
				sp_x = sp_x / K_x_half[h1] + psi_sxx_x[j][h1];
			}

		
			/* top boundary */                                         
			if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
						
				psi_syy_y[j][i] = b_y_half[j] * psi_syy_y[j][i] + a_y_half[j] * sp_y;                                                
				sp_y = sp_y / K_y_half[j] + psi_syy_y[j][i];
			}
			

			/* bottom boundary */                                         
			if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){
				
				h1 = (j-ny2+2*FW);
				h = j;
							
				/*psi_syy_y[j][i] = b_y_half[j] * psi_syy_y[j][i] + a_y_half[j] * syy_y;                                                
				syy_y = syy_y / K_y_half[j] + psi_syy_y[j][i];*/
				
				psi_syy_y[h1][i] = b_y_half[h1] * psi_syy_y[h1][i] + a_y_half[h1] * sp_y;                                                
				sp_y = sp_y / K_y_half[h1] + psi_syy_y[h1][i];
			}                       
			
			if(sw==0){
				if (VELOCITY==0){
					vxp1[j][i] = rip[j][i]*(sp_x)/DH;
					vyp1[j][i] = rjp[j][i]*(sp_y)/DH;}
				else{
					vxp1[j][i] = rip[j][i]*(sp_x)/DH;
					vyp1[j][i] = rjp[j][i]*(sp_y)/DH;}
			}
			
			if(sw==1){
				if (VELOCITY==0){
					vxp1[j][i] += vx[j][i]*DT;
					vyp1[j][i] += vy[j][i]*DT;}
				else{
					vxp1[j][i] = vx[j][i];
					vyp1[j][i] = vy[j][i];}
			}
			
			vx[j][i] += DT*rip[j][i]*(sp_x)/DH;
			vy[j][i] += DT*rjp[j][i]*(sp_y)/DH; 		         

		}
		}

	
	break;
		
	case 8:

		for (j=ny1;j<=ny2;j++){
		for (i=nx1;i<=nx2;i++){

			sp_x =  hc[1]*(sp[j][i+1]-sp[j][i])
					+ hc[2]*(sp[j][i+2]-sp[j][i-1])
					+ hc[3]*(sp[j][i+3]-sp[j][i-2])
					+ hc[4]*(sp[j][i+4]-sp[j][i-3]);

			sp_y = hc[1]*(sp[j+1][i]-sp[j][i])
					+ hc[2]*(sp[j+2][i]-sp[j-1][i])
					+ hc[3]*(sp[j+3][i]-sp[j-2][i])
					+ hc[4]*(sp[j+4][i]-sp[j-3][i]);
			

			/* left boundary */                                         
			if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
					
				psi_sxx_x[j][i] = b_x_half[i] * psi_sxx_x[j][i] + a_x_half[i] * sp_x;                                                
				sp_x = sp_x / K_x_half[i] + psi_sxx_x[j][i];
			}

			/* right boundary */                                         
			if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){

				h1 = (i-nx2+2*FW);
				h=i; 
					
				/*psi_sxx_x[j][i] = b_x_half[i] * psi_sxx_x[j][i] + a_x_half[i] * sxx_x;                                                
				sxx_x = sxx_x / K_x_half[i] + psi_sxx_x[j][i];*/
				
				psi_sxx_x[j][h1] = b_x_half[h1] * psi_sxx_x[j][h1] + a_x_half[h1] * sp_x;                                                
				sp_x = sp_x / K_x_half[h1] + psi_sxx_x[j][h1];
			}

		
			/* top boundary */                                         
			if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
						
				psi_syy_y[j][i] = b_y_half[j] * psi_syy_y[j][i] + a_y_half[j] * sp_y;                                                
				sp_y = sp_y / K_y_half[j] + psi_syy_y[j][i];
			}
			

			/* bottom boundary */                                         
			if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){
				
				h1 = (j-ny2+2*FW);
				h = j;
							
				/*psi_syy_y[j][i] = b_y_half[j] * psi_syy_y[j][i] + a_y_half[j] * syy_y;                                                
				syy_y = syy_y / K_y_half[j] + psi_syy_y[j][i];*/
				
				psi_syy_y[h1][i] = b_y_half[h1] * psi_syy_y[h1][i] + a_y_half[h1] * sp_y;                                                
				sp_y = sp_y / K_y_half[h1] + psi_syy_y[h1][i];
				
			}                       
			
			if(sw==0){
				if (VELOCITY==0){
					vxp1[j][i] = rip[j][i]*(sp_x)/DH;
					vyp1[j][i] = rjp[j][i]*(sp_y)/DH;}
				else{
					vxp1[j][i] = rip[j][i]*(sp_x)/DH;
					vyp1[j][i] = rjp[j][i]*(sp_y)/DH;}
			}
			
			if(sw==1){
				if (VELOCITY==0){
					vxp1[j][i] += vx[j][i]*DT;
					vyp1[j][i] += vy[j][i]*DT;}
				else{
					vxp1[j][i] = vx[j][i];
					vyp1[j][i] = vy[j][i];}
			}
			
			vx[j][i] += DT*rip[j][i]*(sp_x)/DH;
			vy[j][i] += DT*rjp[j][i]*(sp_y)/DH;
		}
		}
	
	break;
		
	default:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				vxtmp = 0;
				vytmp = 0;
				for (m=1; m<=fdoh; m++) {
					vxtmp +=   hc[m]*( sp[j][i+m]   - sp[j][i-m+1] );
							
					vytmp +=   hc[m]*( sp[j+m][i]   - sp[j-m+1][i] );
				}
					
				vx[j][i] += vxtmp*dtdh*rip[j][i];
				vy[j][i] += vytmp*dtdh*rjp[j][i];
			}
		}
		break;

	} /* end of switch(FDORDER) */


 		/* Forward Modelling (sw==0) */
	        if(sw==0){
	        for (l=1;l<=nsrc;l++) {
		    i=(int)srcpos_loc[1][l];
		    j=(int)srcpos_loc[2][l];
		    azi_rad=srcpos_loc[7][l]*PI/180;
		    SOURCE_TYPE=(int)srcpos_loc[8][l];
		
		    if(SOURCE_TYPE==2){vx[j][i] += (DT*rip[j][i]*signals1[l][nt])/(DH*DH);}  /* single force in x */
		    if(SOURCE_TYPE==3){vy[j][i] += (DT*rjp[j][i]*signals1[l][nt])/(DH*DH);}  /* single force in y */
		    if(SOURCE_TYPE==4){vx[j][i] += (DT*rip[j][i]*sin(azi_rad) * signals[l][nt])/(DH*DH);    /* rotated force in x */
		                    vy[j][i] += (DT*rjp[j][i]*cos(azi_rad) * signals[l][nt])/(DH*DH);}  /* rotated force in y */          
		              
		}}
		
		/* Backpropagation (sw==1) */
		if(sw==1){
		for (l=1;l<=nsrc;l++) {
			i=(int)srcpos_loc[1][l];
			j=(int)srcpos_loc[2][l];
		    
			if(ADJOINT_TYPE==1){
				vx[j][i] += signals[l][nt];   /* single force in x */
				vy[j][i] += signals1[l][nt];  /* + single force in y */
			}
			
			if(ADJOINT_TYPE==2)
				vy[j][i] += signals1[l][nt];  /* single force in y */
			if(ADJOINT_TYPE==3)
				vx[j][i] += signals[l][nt];   /* single force in x */
			
		}}                         
	

	if (infoout && (MYID==0)){
		time2=MPI_Wtime();
		fprintf(FP," finished (real time: %4.2f s).\n",time2-time1);
	}
}
