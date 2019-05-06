/*------------------------------------------------------------------------
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
--------------------------------------------------------------------------*/

/* $Id: update_s_visc_interior.c 819 2015-04-17 11:07:06Z tmetz $*/
/*------------------------------------------------------------------------
 *   updating stress components at interior gridpoints (excluding boundarys) [gx2+1...gx3][gy2+1...gy3]
 *   by a staggered grid finite difference scheme of FDORDER accuracy in space
 *   and second order accuracy in time
 *
 *   gx and gy are arrays with the locations of the boundary specified in subgrid_bounds.c
 *   for each subgrid
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void update_s_visc_interior ( int nx1, int nx2, int ny1, int ny2, int *gx, int *gy, int nt,
                              float **vx, float **vy, float **sxx, float **syy,
                              float **sxy, float ***r, float *** p, float ***q,
                              float ** fipjp, float **f, float **g, float *bip, float *bjm, float *cip,
                              float *cjm, float ***d, float ***e, float ***dip,
                              float *hc )
{


	int i,j,l;
	float  vxx, vyy, vxy, vyx, sumr=0.0, sump=0.0, sumq=0.0;
	float  dthalbe, dhi;// dhi2;
	extern float DT, DH;
	extern int L, MYID, FDORDER;
	extern FILE *FP;
	extern int OUTNTIMESTEPINFO;
	double time1=0.0, time2=0.0;

	dthalbe=DT/2.0;
	dhi = 1.0/DH;


	if ( ( MYID==0 ) && ( ( nt+ ( OUTNTIMESTEPINFO-1 ) ) %OUTNTIMESTEPINFO ) ==0 ) {
		time1=MPI_Wtime();
		fprintf ( FP,"\n **Message from update_s_visc_interior (printed by PE %d):\n",MYID );
		fprintf ( FP," Updating stress components ..." );
	}


	switch ( FDORDER ) { /* standard staggered grid (SSG) */
	case 2:
		for ( j=gy[2]+1; j<=gy[3]; j++ ) {
			for ( i=gx[2]+1; i<=gx[3]; i++ ) {
				/* Compute values for shearmodulus u[j][i],
				P-wave modulus pi[j][i],
				tau for S-waves and P-waves taus[j][i],
				taup[j][i] at staggered grid points: */

				/* spatial derivatives of the components of the velocities */
				/* using Holberg coefficients */
				vxx = hc[1]* ( vx[j][i]  -vx[j][i-1] ) *dhi;
				vyy = hc[1]* ( vy[j][i]  -vy[j-1][i] ) *dhi;
				vyx = hc[1]* ( vy[j][i+1]-vy[j][i] ) *dhi;
				vxy = hc[1]* ( vx[j+1][i]-vx[j][i] ) *dhi;

				/* computing sums of the old memory variables */
				sumr=sump=sumq=0.0;
				for ( l=1; l<=L; l++ ) {
					sumr+=r[j][i][l];
					sump+=p[j][i][l];
					sumq+=q[j][i][l];
				}

				/* updating components of the stress tensor, partially */
				sxy[j][i] += ( fipjp[j][i]* ( vxy+vyx ) ) + ( dthalbe*sumr );
				sxx[j][i] += ( g[j][i]* ( vxx+vyy ) )- ( 2.0*f[j][i]*vyy ) + ( dthalbe*sump );
				syy[j][i] += ( g[j][i]* ( vxx+vyy ) )- ( 2.0*f[j][i]*vxx ) + ( dthalbe*sumq );

				/* now updating the memory-variables and sum them up*/
				sumr=sump=sumq=0.0;
				for ( l=1; l<=L; l++ ) {
					r[j][i][l] = bip[l]* ( r[j][i][l]*cip[l]- ( dip[j][i][l]* ( vxy+vyx ) ) );
					p[j][i][l] = bjm[l]* ( p[j][i][l]*cjm[l]- ( e[j][i][l]* ( vxx+vyy ) ) + ( 2.0*d[j][i][l]*vyy ) );
					q[j][i][l] = bjm[l]* ( q[j][i][l]*cjm[l]- ( e[j][i][l]* ( vxx+vyy ) ) + ( 2.0*d[j][i][l]*vxx ) );
					sumr += r[j][i][l];
					sump += p[j][i][l];
					sumq += q[j][i][l];
				}

				/* and now the components of the stress tensor are
				completely updated */
				sxy[j][i]+= ( dthalbe*sumr );
				sxx[j][i]+= ( dthalbe*sump );
				syy[j][i]+= ( dthalbe*sumq );
			}
		}
		break;

	case 4:
		for ( j=gy[2]+1; j<=gy[3]; j++ ) {
			for ( i=gx[2]+1; i<=gx[3]; i++ ) {
				vxx = ( hc[1]* ( vx[j][i]  -vx[j][i-1] )
				        + hc[2]* ( vx[j][i+1]-vx[j][i-2] )
				      ) *dhi;
				vyy = ( hc[1]* ( vy[j][i]  -vy[j-1][i] )
				        + hc[2]* ( vy[j+1][i]-vy[j-2][i] )
				      ) *dhi;
				vyx = ( hc[1]* ( vy[j][i+1]-vy[j][i] )
				        + hc[2]* ( vy[j][i+2]-vy[j][i-1] )
				      ) *dhi;
				vxy = ( hc[1]* ( vx[j+1][i]-vx[j][i] )
				        + hc[2]* ( vx[j+2][i]-vx[j-1][i] )
				      ) *dhi;

				sumr=sump=sumq=0.0;
				for ( l=1; l<=L; l++ ) {
					sumr+=r[j][i][l];
					sump+=p[j][i][l];
					sumq+=q[j][i][l];
				}

				sxy[j][i] += ( fipjp[j][i]* ( vxy+vyx ) ) + ( dthalbe*sumr );
				sxx[j][i] += ( g[j][i]* ( vxx+vyy ) )- ( 2.0*f[j][i]*vyy ) + ( dthalbe*sump );
				syy[j][i] += ( g[j][i]* ( vxx+vyy ) )- ( 2.0*f[j][i]*vxx ) + ( dthalbe*sumq );

				sumr=sump=sumq=0.0;
				for ( l=1; l<=L; l++ ) {
					r[j][i][l] = bip[l]* ( r[j][i][l]*cip[l]- ( dip[j][i][l]* ( vxy+vyx ) ) );
					p[j][i][l] = bjm[l]* ( p[j][i][l]*cjm[l]- ( e[j][i][l]* ( vxx+vyy ) ) + ( 2.0*d[j][i][l]*vyy ) );
					q[j][i][l] = bjm[l]* ( q[j][i][l]*cjm[l]- ( e[j][i][l]* ( vxx+vyy ) ) + ( 2.0*d[j][i][l]*vxx ) );
					sumr += r[j][i][l];
					sump += p[j][i][l];
					sumq += q[j][i][l];
				}

				sxy[j][i]+= ( dthalbe*sumr );
				sxx[j][i]+= ( dthalbe*sump );
				syy[j][i]+= ( dthalbe*sumq );
			}
		}
		break;

	case 6:
		for ( j=gy[2]+1; j<=gy[3]; j++ ) {
			for ( i=gx[2]+1; i<=gx[3]; i++ ) {
				vxx = ( hc[1]* ( vx[j][i]  -vx[j][i-1] )
				        + hc[2]* ( vx[j][i+1]-vx[j][i-2] )
				        + hc[3]* ( vx[j][i+2]-vx[j][i-3] )
				      ) *dhi;
				vyy = ( hc[1]* ( vy[j][i]  -vy[j-1][i] )
				        + hc[2]* ( vy[j+1][i]-vy[j-2][i] )
				        + hc[3]* ( vy[j+2][i]-vy[j-3][i] )
				      ) *dhi;
				vyx = ( hc[1]* ( vy[j][i+1]-vy[j][i] )
				        + hc[2]* ( vy[j][i+2]-vy[j][i-1] )
				        + hc[3]* ( vy[j][i+3]-vy[j][i-2] )
				      ) *dhi;
				vxy = ( hc[1]* ( vx[j+1][i]-vx[j][i] )
				        + hc[2]* ( vx[j+2][i]-vx[j-1][i] )
				        + hc[3]* ( vx[j+3][i]-vx[j-2][i] )
				      ) *dhi;

				sumr=sump=sumq=0.0;
				for ( l=1; l<=L; l++ ) {
					sumr+=r[j][i][l];
					sump+=p[j][i][l];
					sumq+=q[j][i][l];
				}

				sxy[j][i] += ( fipjp[j][i]* ( vxy+vyx ) ) + ( dthalbe*sumr );
				sxx[j][i] += ( g[j][i]* ( vxx+vyy ) )- ( 2.0*f[j][i]*vyy ) + ( dthalbe*sump );
				syy[j][i] += ( g[j][i]* ( vxx+vyy ) )- ( 2.0*f[j][i]*vxx ) + ( dthalbe*sumq );

				sumr=sump=sumq=0.0;
				for ( l=1; l<=L; l++ ) {
					r[j][i][l] = bip[l]* ( r[j][i][l]*cip[l]- ( dip[j][i][l]* ( vxy+vyx ) ) );
					p[j][i][l] = bjm[l]* ( p[j][i][l]*cjm[l]- ( e[j][i][l]* ( vxx+vyy ) ) + ( 2.0*d[j][i][l]*vyy ) );
					q[j][i][l] = bjm[l]* ( q[j][i][l]*cjm[l]- ( e[j][i][l]* ( vxx+vyy ) ) + ( 2.0*d[j][i][l]*vxx ) );
					sumr += r[j][i][l];
					sump += p[j][i][l];
					sumq += q[j][i][l];
				}

				sxy[j][i]+= ( dthalbe*sumr );
				sxx[j][i]+= ( dthalbe*sump );
				syy[j][i]+= ( dthalbe*sumq );
			}
		}
		break;

	case 8:
		for ( j=gy[2]+1; j<=gy[3]; j++ ) {
			for ( i=gx[2]+1; i<=gx[3]; i++ ) {
				vxx = ( hc[1]* ( vx[j][i]  -vx[j][i-1] )
				        + hc[2]* ( vx[j][i+1]-vx[j][i-2] )
				        + hc[3]* ( vx[j][i+2]-vx[j][i-3] )
				        + hc[4]* ( vx[j][i+3]-vx[j][i-4] )
				      ) *dhi;
				vyy = ( hc[1]* ( vy[j][i]  -vy[j-1][i] )
				        + hc[2]* ( vy[j+1][i]-vy[j-2][i] )
				        + hc[3]* ( vy[j+2][i]-vy[j-3][i] )
				        + hc[4]* ( vy[j+3][i]-vy[j-4][i] )
				      ) *dhi;
				vyx = ( hc[1]* ( vy[j][i+1]-vy[j][i] )
				        + hc[2]* ( vy[j][i+2]-vy[j][i-1] )
				        + hc[3]* ( vy[j][i+3]-vy[j][i-2] )
				        + hc[4]* ( vy[j][i+4]-vy[j][i-3] )
				      ) *dhi;
				vxy = ( hc[1]* ( vx[j+1][i]-vx[j][i] )
				        + hc[2]* ( vx[j+2][i]-vx[j-1][i] )
				        + hc[3]* ( vx[j+3][i]-vx[j-2][i] )
				        + hc[4]* ( vx[j+4][i]-vx[j-3][i] )
				      ) *dhi;

				sumr=sump=sumq=0.0;
				for ( l=1; l<=L; l++ ) {
					sumr+=r[j][i][l];
					sump+=p[j][i][l];
					sumq+=q[j][i][l];
				}

				sxy[j][i] += ( fipjp[j][i]* ( vxy+vyx ) ) + ( dthalbe*sumr );
				sxx[j][i] += ( g[j][i]* ( vxx+vyy ) )- ( 2.0*f[j][i]*vyy ) + ( dthalbe*sump );
				syy[j][i] += ( g[j][i]* ( vxx+vyy ) )- ( 2.0*f[j][i]*vxx ) + ( dthalbe*sumq );

				sumr=sump=sumq=0.0;
				for ( l=1; l<=L; l++ ) {
					r[j][i][l] = bip[l]* ( r[j][i][l]*cip[l]- ( dip[j][i][l]* ( vxy+vyx ) ) );
					p[j][i][l] = bjm[l]* ( p[j][i][l]*cjm[l]- ( e[j][i][l]* ( vxx+vyy ) ) + ( 2.0*d[j][i][l]*vyy ) );
					q[j][i][l] = bjm[l]* ( q[j][i][l]*cjm[l]- ( e[j][i][l]* ( vxx+vyy ) ) + ( 2.0*d[j][i][l]*vxx ) );
					sumr += r[j][i][l];
					sump += p[j][i][l];
					sumq += q[j][i][l];
				}

				sxy[j][i]+= ( dthalbe*sumr );
				sxx[j][i]+= ( dthalbe*sump );
				syy[j][i]+= ( dthalbe*sumq );
			}
		}
		break;

	case 10:
		for ( j=gy[2]+1; j<=gy[3]; j++ ) {
			for ( i=gx[2]+1; i<=gx[3]; i++ ) {
				vxx = ( hc[1]* ( vx[j][i]  -vx[j][i-1] )
				        + hc[2]* ( vx[j][i+1]-vx[j][i-2] )
				        + hc[3]* ( vx[j][i+2]-vx[j][i-3] )
				        + hc[4]* ( vx[j][i+3]-vx[j][i-4] )
				        + hc[5]* ( vx[j][i+4]-vx[j][i-5] )
				      ) *dhi;
				vyy = ( hc[1]* ( vy[j][i]  -vy[j-1][i] )
				        + hc[2]* ( vy[j+1][i]-vy[j-2][i] )
				        + hc[3]* ( vy[j+2][i]-vy[j-3][i] )
				        + hc[4]* ( vy[j+3][i]-vy[j-4][i] )
				        + hc[5]* ( vy[j+4][i]-vy[j-5][i] )
				      ) *dhi;
				vyx = ( hc[1]* ( vy[j][i+1]-vy[j][i] )
				        + hc[2]* ( vy[j][i+2]-vy[j][i-1] )
				        + hc[3]* ( vy[j][i+3]-vy[j][i-2] )
				        + hc[4]* ( vy[j][i+4]-vy[j][i-3] )
				        + hc[5]* ( vy[j][i+5]-vy[j][i-4] )
				      ) *dhi;
				vxy = ( hc[1]* ( vx[j+1][i]-vx[j][i] )
				        + hc[2]* ( vx[j+2][i]-vx[j-1][i] )
				        + hc[3]* ( vx[j+3][i]-vx[j-2][i] )
				        + hc[4]* ( vx[j+4][i]-vx[j-3][i] )
				        + hc[5]* ( vx[j+5][i]-vx[j-4][i] )
				      ) *dhi;

				sumr = sump = sumq = 0.0;
				for ( l=1; l<=L; l++ ) {
					sumr += r[j][i][l];
					sump += p[j][i][l];
					sumq += q[j][i][l];
				}

				sxy[j][i] += ( fipjp[j][i]* ( vxy+vyx ) ) + ( dthalbe*sumr );
				sxx[j][i] += ( g[j][i]* ( vxx+vyy ) )- ( 2.0*f[j][i]*vyy ) + ( dthalbe*sump );
				syy[j][i] += ( g[j][i]* ( vxx+vyy ) )- ( 2.0*f[j][i]*vxx ) + ( dthalbe*sumq );

				sumr=sump=sumq=0.0;
				for ( l=1; l<=L; l++ ) {
					r[j][i][l] = bip[l]* ( r[j][i][l]*cip[l]- ( dip[j][i][l]* ( vxy+vyx ) ) );
					p[j][i][l] = bjm[l]* ( p[j][i][l]*cjm[l]- ( e[j][i][l]* ( vxx+vyy ) ) + ( 2.0*d[j][i][l]*vyy ) );
					q[j][i][l] = bjm[l]* ( q[j][i][l]*cjm[l]- ( e[j][i][l]* ( vxx+vyy ) ) + ( 2.0*d[j][i][l]*vxx ) );
					sumr += r[j][i][l];
					sump += p[j][i][l];
					sumq += q[j][i][l];
				}

				sxy[j][i] += ( dthalbe*sumr );
				sxx[j][i] += ( dthalbe*sump );
				syy[j][i] += ( dthalbe*sumq );
			}
		}
		break;

	case 12:
		for ( j=gy[2]+1; j<=gy[3]; j++ ) {
			for ( i=gx[2]+1; i<=gx[3]; i++ ) {
				vxx = ( hc[1]* ( vx[j][i]  -vx[j][i-1] )
				        + hc[2]* ( vx[j][i+1]-vx[j][i-2] )
				        + hc[3]* ( vx[j][i+2]-vx[j][i-3] )
				        + hc[4]* ( vx[j][i+3]-vx[j][i-4] )
				        + hc[5]* ( vx[j][i+4]-vx[j][i-5] )
				        + hc[6]* ( vx[j][i+5]-vx[j][i-6] )
				      ) *dhi;
				vyy = ( hc[1]* ( vy[j][i]  -vy[j-1][i] )
				        + hc[2]* ( vy[j+1][i]-vy[j-2][i] )
				        + hc[3]* ( vy[j+2][i]-vy[j-3][i] )
				        + hc[4]* ( vy[j+3][i]-vy[j-4][i] )
				        + hc[5]* ( vy[j+4][i]-vy[j-5][i] )
				        + hc[6]* ( vy[j+5][i]-vy[j-6][i] )
				      ) *dhi;
				vyx = ( hc[1]* ( vy[j][i+1]-vy[j][i] )
				        + hc[2]* ( vy[j][i+2]-vy[j][i-1] )
				        + hc[3]* ( vy[j][i+3]-vy[j][i-2] )
				        + hc[4]* ( vy[j][i+4]-vy[j][i-3] )
				        + hc[5]* ( vy[j][i+5]-vy[j][i-4] )
				        + hc[6]* ( vy[j][i+6]-vy[j][i-5] )
				      ) *dhi;
				vxy = ( hc[1]* ( vx[j+1][i]-vx[j][i] )
				        + hc[2]* ( vx[j+2][i]-vx[j-1][i] )
				        + hc[3]* ( vx[j+3][i]-vx[j-2][i] )
				        + hc[4]* ( vx[j+4][i]-vx[j-3][i] )
				        + hc[5]* ( vx[j+5][i]-vx[j-4][i] )
				        + hc[6]* ( vx[j+6][i]-vx[j-5][i] )
				      ) *dhi;

				sumr=sump=sumq=0.0;
				for ( l=1; l<=L; l++ ) {
					sumr+=r[j][i][l];
					sump+=p[j][i][l];
					sumq+=q[j][i][l];
				}

				sxy[j][i] += ( fipjp[j][i]* ( vxy+vyx ) ) + ( dthalbe*sumr );
				sxx[j][i] += ( g[j][i]* ( vxx+vyy ) )- ( 2.0*f[j][i]*vyy ) + ( dthalbe*sump );
				syy[j][i] += ( g[j][i]* ( vxx+vyy ) )- ( 2.0*f[j][i]*vxx ) + ( dthalbe*sumq );

				sumr=sump=sumq=0.0;
				for ( l=1; l<=L; l++ ) {
					r[j][i][l] = bip[l]* ( r[j][i][l]*cip[l]- ( dip[j][i][l]* ( vxy+vyx ) ) );
					p[j][i][l] = bjm[l]* ( p[j][i][l]*cjm[l]- ( e[j][i][l]* ( vxx+vyy ) ) + ( 2.0*d[j][i][l]*vyy ) );
					q[j][i][l] = bjm[l]* ( q[j][i][l]*cjm[l]- ( e[j][i][l]* ( vxx+vyy ) ) + ( 2.0*d[j][i][l]*vxx ) );
					sumr += r[j][i][l];
					sump += p[j][i][l];
					sumq += q[j][i][l];
				}

				sxy[j][i]+= ( dthalbe*sumr );
				sxx[j][i]+= ( dthalbe*sump );
				syy[j][i]+= ( dthalbe*sumq );
			}
		}
		break;


	default: /* Case 2 */
		for ( j=gy[2]+1; j<=gy[3]; j++ ) {
			for ( i=gx[2]+1; i<=gx[3]; i++ ) {
				/* Compute values for shearmodulus u[j][i],
				P-wave modulus pi[j][i],
				tau for S-waves and P-waves taus[j][i],
				taup[j][i] at staggered grid points: */

				/* spatial derivatives of the components of the velocities */
				/* using Holberg coefficients */
				vxx = hc[1]* ( vx[j][i]  -vx[j][i-1] ) *dhi;
				vyy = hc[1]* ( vy[j][i]  -vy[j-1][i] ) *dhi;
				vyx = hc[1]* ( vy[j][i+1]-vy[j][i] ) *dhi;
				vxy = hc[1]* ( vx[j+1][i]-vx[j][i] ) *dhi;

				/* computing sums of the old memory variables */
				sumr=sump=sumq=0.0;
				for ( l=1; l<=L; l++ ) {
					sumr+=r[j][i][l];
					sump+=p[j][i][l];
					sumq+=q[j][i][l];
				}

				/* updating components of the stress tensor, partially */
				sxy[j][i] += ( fipjp[j][i]* ( vxy+vyx ) ) + ( dthalbe*sumr );
				sxx[j][i] += ( g[j][i]* ( vxx+vyy ) )- ( 2.0*f[j][i]*vyy ) + ( dthalbe*sump );
				syy[j][i] += ( g[j][i]* ( vxx+vyy ) )- ( 2.0*f[j][i]*vxx ) + ( dthalbe*sumq );

				/* now updating the memory-variables and sum them up*/
				sumr=sump=sumq=0.0;
				for ( l=1; l<=L; l++ ) {
					r[j][i][l] = bip[l]* ( r[j][i][l]*cip[l]- ( dip[j][i][l]* ( vxy+vyx ) ) );
					p[j][i][l] = bjm[l]* ( p[j][i][l]*cjm[l]- ( e[j][i][l]* ( vxx+vyy ) ) + ( 2.0*d[j][i][l]*vyy ) );
					q[j][i][l] = bjm[l]* ( q[j][i][l]*cjm[l]- ( e[j][i][l]* ( vxx+vyy ) ) + ( 2.0*d[j][i][l]*vxx ) );
					sumr += r[j][i][l];
					sump += p[j][i][l];
					sumq += q[j][i][l];
				}

				/* and now the components of the stress tensor are
				completely updated */
				sxy[j][i]+= ( dthalbe*sumr );
				sxx[j][i]+= ( dthalbe*sump );
				syy[j][i]+= ( dthalbe*sumq );
			}
		}
		break;

	} /* end of switch(FDORDER) */



	if ( ( MYID==0 ) && ( ( nt+ ( OUTNTIMESTEPINFO-1 ) ) %OUTNTIMESTEPINFO ) ==0 ) {
		time2=MPI_Wtime();
		fprintf ( FP," finished (real time: %4.3f s).\n",time2-time1 );
	}
}
