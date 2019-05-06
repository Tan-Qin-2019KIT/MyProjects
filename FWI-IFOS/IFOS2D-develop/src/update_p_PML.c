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
 *   updating stress components at gridpoints [nx1...nx2][ny1...ny2]
 *   by a staggered grid finite difference scheme of arbitrary (FDORDER) order accuracy in space
 *   and second order accuracy in time
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void update_p_PML(int nx1, int nx2, int ny1, int ny2, 	float **  vx, float ** vy, float ** sp, float ** u, float ** pi, float ** absorb_coeff, float **rho, float *hc, int infoout,
		  float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
		  float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half, 
		  float ** psi_vxx, float ** psi_vyy, float ** psi_vxy, float ** psi_vyx){
	
	int i,j, m, fdoh, h, h1;
	float g;
	float  vxx, vyy;
	float  dhi;
	extern float DT, DH;
	extern int MYID, FDORDER, PARAMETERIZATION, FW;
        extern int FREE_SURF, BOUNDARY;
	extern int NPROCX, NPROCY, POS[3];
	extern FILE *FP;
	double time1, time2;
	
	
	dhi = DT/DH;
	fdoh = FDORDER/2;
	
	if (infoout && (MYID==0)){
		time1=MPI_Wtime();
		fprintf(FP,"\n **Message from update_p (printed by PE %d):\n",MYID);
		fprintf(FP," Updating stress components ...");
	}
	
	switch (FDORDER){
	
	case 2:
		for (j=ny1;j<=ny2;j++){
		for (i=nx1;i<=nx2;i++){
			vxx = (  hc[1]*(vx[j][i]  -vx[j][i-1]))*dhi;
			vyy = (  hc[1]*(vy[j][i]  -vy[j-1][i]))*dhi; 
			
			/* left boundary */                                         
			if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
				
				psi_vxx[j][i] = b_x[i] * psi_vxx[j][i] + a_x[i] * vxx;
				vxx = vxx / K_x[i] + psi_vxx[j][i];
			}
			
			/* right boundary */                                         
			if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){
				
				h1 = (i-nx2+2*FW);
				h = i;
				
				psi_vxx[j][h1] = b_x[h1] * psi_vxx[j][h1] + a_x[h1] * vxx;
				vxx = vxx / K_x[h1] + psi_vxx[j][h1];
			}
			
			/* top boundary */                                         
			if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
				
				psi_vyy[j][i] = b_y[j] * psi_vyy[j][i] + a_y[j] * vyy;                                            
				vyy = vyy / K_y[j] + psi_vyy[j][i];
			}
			
			/* bottom boundary */                                         
			if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){
				
				h1 = (j-ny2+2*FW);                                        
				h = j;
				
				psi_vyy[h1][i] = b_y[h1] * psi_vyy[h1][i] + a_y[h1] * vyy;                                            
				vyy = vyy / K_y[h1] + psi_vyy[h1][i];
			}
			
			/* lambda - mu relationship*/
			if (PARAMETERIZATION==1){
				g = rho[j][i] * (pi[j][i] * pi[j][i]);
			}
			
			sp[j][i] += g*(vxx+vyy);
			
            u[j][i] = (g/DT)*(vxx+vyy);
		}
		}
	break;

	case 4:
		for (j=ny1;j<=ny2;j++){
		for (i=nx1;i<=nx2;i++){
			vxx = (  hc[1]*(vx[j][i]  -vx[j][i-1])
				       + hc[2]*(vx[j][i+1]-vx[j][i-2]))*dhi;
			
                        vyy = (  hc[1]*(vy[j][i]  -vy[j-1][i])
				       + hc[2]*(vy[j+1][i]-vy[j-2][i]))*dhi; 

			/* left boundary */                                         
			if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
				
				psi_vxx[j][i] = b_x[i] * psi_vxx[j][i] + a_x[i] * vxx;
				vxx = vxx / K_x[i] + psi_vxx[j][i];             
			}

			/* right boundary */                                         
			if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){
			
				h1 = (i-nx2+2*FW);
				h = i;
				
				psi_vxx[j][h1] = b_x[h1] * psi_vxx[j][h1] + a_x[h1] * vxx;
				vxx = vxx / K_x[h1] + psi_vxx[j][h1];
			}

			/* top boundary */                                         
			if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
				
				psi_vyy[j][i] = b_y[j] * psi_vyy[j][i] + a_y[j] * vyy;
				vyy = vyy / K_y[j] + psi_vyy[j][i];
			}
		
			/* bottom boundary */                                         
			if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){

				h1 = (j-ny2+2*FW);                                        
				h = j;
				
				psi_vyy[h1][i] = b_y[h1] * psi_vyy[h1][i] + a_y[h1] * vyy;                                            
				vyy = vyy / K_y[h1] + psi_vyy[h1][i];
			}
			
			if (PARAMETERIZATION==1){
				g = rho[j][i] * (pi[j][i] * pi[j][i]);}
			
			sp[j][i] += g*(vxx+vyy);
            
            u[j][i] = (g/DT)*(vxx+vyy);
		}
		}
	break;

	case 6:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
			vxx = (  hc[1]*(vx[j][i]  -vx[j][i-1])
				       + hc[2]*(vx[j][i+1]-vx[j][i-2])
				       + hc[3]*(vx[j][i+2]-vx[j][i-3]))*dhi;
			
                        vyy = (  hc[1]*(vy[j][i]  -vy[j-1][i])
				       + hc[2]*(vy[j+1][i]-vy[j-2][i])
				       + hc[3]*(vy[j+2][i]-vy[j-3][i]))*dhi; 

			/* left boundary */                                         
			if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
				
				psi_vxx[j][i] = b_x[i] * psi_vxx[j][i] + a_x[i] * vxx;
				vxx = vxx / K_x[i] + psi_vxx[j][i];
			}

			/* right boundary */                                         
			if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){
			
				h1 = (i-nx2+2*FW);
				h = i;
				
				psi_vxx[j][h1] = b_x[h1] * psi_vxx[j][h1] + a_x[h1] * vxx;
				vxx = vxx / K_x[h1] + psi_vxx[j][h1]; 
			}

			/* top boundary */                                         
			if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
				
				psi_vyy[j][i] = b_y[j] * psi_vyy[j][i] + a_y[j] * vyy;
				vyy = vyy / K_y[j] + psi_vyy[j][i];

			}
			
			/* bottom boundary */                                         
			if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){

				h1 = (j-ny2+2*FW);                                        
				h = j;
				
				psi_vyy[h1][i] = b_y[h1] * psi_vyy[h1][i] + a_y[h1] * vyy;                                            
				vyy = vyy / K_y[h1] + psi_vyy[h1][i];
			}
			
			if (PARAMETERIZATION==1){
				g = rho[j][i] * (pi[j][i] * pi[j][i]);}
			
			sp[j][i] += g*(vxx+vyy);
            
            u[j][i] = (g/DT)*(vxx+vyy);
		}
		}
	break;

	case 8:

		for (j=ny1;j<=ny2;j++){
		for (i=nx1;i<=nx2;i++){

			vxx = (  hc[1]*(vx[j][i]  -vx[j][i-1])
				       + hc[2]*(vx[j][i+1]-vx[j][i-2])
				       + hc[3]*(vx[j][i+2]-vx[j][i-3])
				       + hc[4]*(vx[j][i+3]-vx[j][i-4]))*dhi;
			
			vyy = (  hc[1]*(vy[j][i]  -vy[j-1][i])
				       + hc[2]*(vy[j+1][i]-vy[j-2][i])
				       + hc[3]*(vy[j+2][i]-vy[j-3][i])
				       + hc[4]*(vy[j+3][i]-vy[j-4][i]))*dhi; 

			/* left boundary */                                         
			if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
				
				psi_vxx[j][i] = b_x[i] * psi_vxx[j][i] + a_x[i] * vxx;
				vxx = vxx / K_x[i] + psi_vxx[j][i];        
			}

			/* right boundary */                                         
			if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=nx2-FW+1)){
			
				h1 = (i-nx2+2*FW);
				h = i;
				
				psi_vxx[j][h1] = b_x[h1] * psi_vxx[j][h1] + a_x[h1] * vxx;
				vxx = vxx / K_x[h1] + psi_vxx[j][h1]; 
			}

			/* top boundary */                                         
			if((POS[2]==0) && (!(FREE_SURF)) && (j<=FW)){
				
				psi_vyy[j][i] = b_y[j] * psi_vyy[j][i] + a_y[j] * vyy;
				vyy = vyy / K_y[j] + psi_vyy[j][i];
			}
			
			/* bottom boundary */                                         
			if((POS[2]==NPROCY-1) && (j>=ny2-FW+1)){

				h1 = (j-ny2+2*FW);                                        
				h = j;
				
				psi_vyy[h1][i] = b_y[h1] * psi_vyy[h1][i] + a_y[h1] * vyy;                                            
				vyy = vyy / K_y[h1] + psi_vyy[h1][i];
			}
			
			if (PARAMETERIZATION==1){
				g = rho[j][i] * (pi[j][i] * pi[j][i]);}
			
			sp[j][i] += g*(vxx+vyy);
            
            u[j][i] = (g/DT)*(vxx+vyy);
		}
		}
	break;
	
	default:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				vxx = 0.0;
				vyy = 0.0;
				for (m=1; m<=fdoh; m++) {
					vxx += hc[m]*(vx[j][i+m-1] -vx[j][i-m]  );
					vyy += hc[m]*(vy[j+m-1][i] -vy[j-m][i]  );
				}
				
				g=pi[j][i]*DT;
				
				sp[j][i]+=(g*(vxx+vyy))*dhi;
			}
		}
	break;
		
	} /* end of switch(FDORDER) */
	
	if (infoout && (MYID==0)){
		time2=MPI_Wtime();
		fprintf(FP," finished (real time: %4.2f s).\n",time2-time1);
	}
}
