/*                          Acknowledgement                         */
/* This function is copied from DENISE Black Edition from D. Koehn  */
/* Licence: GNU GENERAL PUBLIC LICENSE Version 2, June 1991         */
/* https://github.com/daniel-koehn/DENISE-Black-Edition             */

#include "fd.h"
void eprecond1(float ** We, float ** Ws, float ** Wr, float epsilon){
    
    extern int NX, NY, IDX, IDY, DTINV, EPRECOND, VERBOSE;
    extern int POS[3], NXG;
    extern float DH;
    int i, j, ii, jj;
    float maxWetmp, maxWe, x, y, xmin, xmax;
    xmin = DH;
    xmax = NXG*DH;
    
    
    maxWetmp=0.0;
    
    if(EPRECOND==1){
        /* calculate energy weighting */
        for (j=1;j<=NY;j=j+IDY){
            for (i=1;i<=NX;i=i+IDX){
                
                We[j][i]=sqrt(Ws[j][i]*Wr[j][i]); /* energy weighted source and receiver contribution */
                
                /* estimate maximum We on this CPU*/
                if(We[j][i]>maxWetmp){
                    maxWetmp = We[j][i];
                }
                
            }
        }
    }
    
    if(EPRECOND==3){
        /* Forward wavefield + approximation of the receiver Greens function (Plessix and Mulder, 2004) */
        for (j=1;j<=NY;j=j+IDY){
            for (i=1;i<=NX;i=i+IDX){
                
                /* calculate global coordinates */
                ii=i+POS[1]*NX;
                jj=j+POS[2]*NY;
                
                x = ii*DH;
                y = jj*DH;
                
                We[j][i]= sqrt(Ws[j][i]) * (asinh((xmax-x)/y)-asinh((xmin-x)/y));
                
                /* estimate maximum We on this CPU*/
                if(We[j][i]>maxWetmp){
                    maxWetmp = We[j][i];
                }
                
            }
        }
    }
    
    /* estimate maximum of We */
    MPI_Allreduce(&maxWetmp,&maxWe,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
    
    /* regularize energy weighting to avoid divison by zero */
    for (j=1;j<=NY;j=j+IDY){
        for (i=1;i<=NX;i=i+IDX){
            
            We[j][i] = We[j][i] + (epsilon*maxWe);
            
        }
    }
    
    
}
