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
 *   Write FD-Parameters to stdout
 *  ----------------------------------------------------------------------*/

#include "fd.h"

/* printing all important parameters on stdout */
void write_par(FILE *fp){
    
    /* declaration of extern variables */
    extern int   NX, NY, NT, SOURCE_SHAPE, SOURCE_TYPE, FDORDER, MAXRELERROR;
    extern int  SNAP, SNAP_FORMAT, ACOUSTIC, L, SRCREC, TAPER;
    extern float DH, TIME, DT, TS, *FL, TAU, VPPML, PLANE_WAVE_DEPTH, PHI, FPML, npower, k_max_PML, F_REF;
    extern float REC_ARRAY_DEPTH, REC_ARRAY_DIST;
    extern float XREC1, XREC2, YREC1, YREC2;
    extern int SEISMO, NDT, NGEOPH, SEIS_FORMAT, FREE_SURF, FW;
    extern int  READMOD, READREC, BOUNDARY, REC_ARRAY, DRX, INVTYPE;
    extern float TSNAP1, TSNAP2, TSNAPINC, REFREC[4];
    extern char SNAP_FILE[STRING_SIZE], SOURCE_FILE[STRING_SIZE], REC_FILE[STRING_SIZE];
    extern char SEIS_FILE[STRING_SIZE];
    extern char SIGNAL_FILE[STRING_SIZE];
    extern char  MFILE[STRING_SIZE], JACOBIAN[STRING_SIZE], DATA_DIR[STRING_SIZE],FREQ_FILE[STRING_SIZE];
    extern int NP, NPROCX, NPROCY, MYID;
    
    extern int GRADT1, GRADT2, GRADT3, GRADT4, ITERMAX, PARAMETERIZATION, FORWARD_ONLY, ADJOINT_TYPE;
    extern int  GRAD_METHOD;
    extern float TSHIFT_back;
    extern int FILT_SIZE, MODEL_FILTER;
    extern int FILT_SIZE_GRAD, GRAD_FILTER;
    
    extern int TESTSHOT_START, TESTSHOT_END, TESTSHOT_INCR, NO_OF_TESTSHOTS;
    extern int SWS_TAPER_GRAD_VERT, SWS_TAPER_GRAD_HOR, SWS_TAPER_GRAD_SOURCES, SWS_TAPER_CIRCULAR_PER_SHOT, SRTSHAPE, FILTSIZE;
    extern int SWS_TAPER_FILE, SWS_TAPER_FILE_PER_SHOT;
    extern float SRTRADIUS;
    extern char TAPER_FILE_NAME[STRING_SIZE];
    extern int SPATFILTER, SPAT_FILT_SIZE, SPAT_FILT_1, SPAT_FILT_ITER;
    extern int INV_RHO_ITER, INV_VP_ITER, INV_VS_ITER;
    extern int MIN_ITER;;
    extern char INV_MODELFILE[STRING_SIZE];
    extern int nfstart, nf;
    extern int nfstart_jac, nf_jac;
    extern float VPUPPERLIM, VPLOWERLIM, VSUPPERLIM, VSLOWERLIM, RHOUPPERLIM, RHOLOWERLIM;
    
    extern int INV_STF, N_STF, N_STF_START;
    extern char PARA[STRING_SIZE];
    
    extern int TIME_FILT, ORDER;
    extern float F_LOW_PASS_START, F_LOW_PASS_END, F_LOW_PASS_INCR;
    extern int LNORM, DTINV;
    extern int STEPMAX, TRKILL, TRKILL_STF;
    
    extern int TRKILL_OFFSET;
    extern float TRKILL_OFFSET_LOWER;
    extern float TRKILL_OFFSET_UPPER;
    
    extern float EPS_SCALE, SCALEFAC;
    extern float PRO;
    extern char  TRKILL_FILE[STRING_SIZE], TRKILL_FILE_STF[STRING_SIZE];
    
    extern int TIMEWIN, NORMALIZE;
    extern float TWLENGTH_PLUS, TWLENGTH_MINUS, GAMMA;
    extern char PICKS_FILE[STRING_SIZE];
    
    extern char MISFIT_LOG_FILE[STRING_SIZE];
    
    extern int VELOCITY;
    
    extern float WATERLEVEL_LNORM8;
    
    extern float VP_VS_RATIO;
    
    extern int S;
    
    extern float S_VP, S_VS, S_RHO;
    
    extern int GRAD_FILT_WAVELENGTH;
    extern float A;
    
    /* definition of local variables */
    int l;
    
    
    fprintf(fp,"\n **Message from write_par (printed by PE %d):\n",MYID);
    fprintf(fp,"\n");
    fprintf(fp,"------------------------- Processors ------------------------\n");
    fprintf(fp," Number of PEs in x-direction (NPROCX): %d\n",NPROCX);
    fprintf(fp," Number of PEs in vertical direction (NPROCY): %d\n",NPROCY);
    fprintf(fp," Total number of PEs in use: %d\n",NP);
    fprintf(fp,"\n");
    fprintf(fp," ----------------------- Discretization  ---------------------\n");
    fprintf(fp," Number of gridpoints in x-direction (NX): %i\n", NX);
    fprintf(fp," Number of gridpoints in y-direction (NY): %i\n", NY);
    fprintf(fp," Grid-spacing (DH): %e meter\n", DH);
    fprintf(fp," Time of wave propagation (T): %e seconds\n",TIME);
    fprintf(fp," Timestep (DT): %e seconds\n", DT);
    fprintf(fp," Number of timesteps: %i \n",NT);
    fprintf(fp,"\n");
    fprintf(fp," ------------------------- FD ORDER -----------------------------\n");
    fprintf(fp," FDORDER = %d\n",FDORDER);
    fprintf(fp," MAXRELERROR = %d\n",MAXRELERROR);
    fprintf(fp,"\n");
    fprintf(fp," ------------------------- SOURCE -----------------------------\n");
    
    if (SRCREC){
        fprintf(fp," reading source positions, time delay, centre frequency \n");
        fprintf(fp," and initial amplitude from ASCII-file \n");
        fprintf(fp,"\t%s\n\n",SOURCE_FILE);
    } else {
        fprintf(fp," plane wave excitation: depth= %5.2f meter \n",PLANE_WAVE_DEPTH);
        fprintf(fp," incidence angle of plane P-wave (from vertical) PHI= %5.2f degrees \n",PHI);
        fprintf(fp," duration of source signal: %e seconds\n",TS);
        fprintf(fp," (centre frequency is approximately %e Hz)\n",1.0/TS);
    }
    
    
    fprintf(fp," wavelet of source:");
    
    switch (SOURCE_SHAPE){
        case 1 :
            fprintf(fp," Ricker\n");
            break;
        case 2 :
            fprintf(fp," Fuchs-Mueller\n");
            break;
        case 3 :
            fprintf(fp," reading from \n\t %s\n",SIGNAL_FILE);
            break;
        case 4 :
            fprintf(fp," sinus raised to the power of 3.0 \n");
            break;
        case 5 :
            fprintf(fp," 1st derivative of Gaussian \n");
            break;
        case 6 :
            fprintf(fp," Bandlimited Spike \n");
            break;
        case 7 :
            fprintf(fp," reading from \n\t %s.shot*.su in su format (one file for each shot)\n",SIGNAL_FILE);
            break;
        case 8 :
            fprintf(fp," Integral of sin^3 function\n");
            break;
        default :
            declare_error(" Sorry, incorrect specification of source wavelet ! ");
    }
    
    fprintf(fp," Type of source:");
    switch (SOURCE_TYPE){
        case 1 :
            fprintf(fp," explosive source \n");
            break;
        case 2 :
            fprintf(fp," point source with directive force in x-direction\n");
            break;
        case 3 :
            fprintf(fp," point source with directive force in (vertical) y-direction\n");
            break;
        case 4 :
            fprintf(fp," rotated point source with directive force in x- and y-direction\n");
            break;
        default :
            declare_error(" Sorry, wrong source type specification ! ");
    }
    
    fprintf(fp,"\n");
    
    if (SEISMO){
        fprintf(fp," ------------------------- RECEIVER  --------------------------\n");
        if (READREC==1){
            fprintf(fp," reading receiver positions from file \n");
            fprintf(fp,"\t%s\n\n",REC_FILE);
            fprintf(fp," reference_point_for_receiver_coordinate_system:\n");
            fprintf(fp," x=%f \ty=%f\t z=%f\n",REFREC[1], REFREC[2], REFREC[3]);
        } else if (REC_ARRAY>0){
            fprintf(fp," Horizontal lines of receivers.\n");
            fprintf(fp," number of lines: %d \n",REC_ARRAY);
            fprintf(fp," depth of upper line: %e m \n",REC_ARRAY_DEPTH);
            fprintf(fp," vertical increment between lines: %e m \n",REC_ARRAY_DIST);
            fprintf(fp," distance between receivers in x-direction within line: %i \n", DRX);
        }else{
            
            fprintf(fp," first receiver position (XREC1,YREC1) = (%e, %e) m\n",
                    XREC1,YREC1);
            fprintf(fp," last receiver position (XREC1,YREC1) = (%e, %e) m\n",
                    XREC2,YREC2);
            fprintf(fp,"\n");
        }
    }
    
    fprintf(fp," ------------------------- FREE SURFACE ------------------------\n");
    if (FREE_SURF) fprintf(fp," free surface at the top of the model ! \n");
    else fprintf(fp," no free surface at the top of the model ! \n");
    fprintf(fp,"\n");
    
    fprintf(fp," ------------------------- CPML ---------------------\n");
    if (FW>0.0){
        fprintf(fp," width of absorbing frame is %i gridpoints.\n",FW);
        fprintf(fp," CPML VPPML applied. \n");
        fprintf(fp," VPPML velocity in the PML frame in m/s: %f .\n",VPPML);
        fprintf(fp," Frequency within the PML frame in Hz: %f \n",FPML);
        fprintf(fp," npower: %f \n",npower);
        fprintf(fp," k_max: %f \n",k_max_PML);
    }
    else fprintf(fp," absorbing frame not installed ! \n");
    
    
    switch (BOUNDARY){
        case 0 :
            fprintf(fp," No periodic boundary condition.\n");
            break;
        case 1 :
            fprintf(fp," Periodic boundary condition at left and right edges.\n");
            break;
        default :
            warning(" Wrong integer value for BOUNDARY specified ! ");
            warning(" No periodic boundary condition will be applied ");
            BOUNDARY=0;
            break;
    }
    
    if (READMOD){
        fprintf(fp," ------------------------- MODEL-FILES -------------------------\n");
        fprintf(fp," names of model-files: \n");
        if(!ACOUSTIC){
            fprintf(fp,"\t shear wave velocities:\n\t %s.vs\n",MFILE);
            fprintf(fp,"\t tau for shear waves:\n\t %s.ts\n",MFILE);}
        if(ACOUSTIC) fprintf(fp,"\t acoustic modelling!\n");
        fprintf(fp,"\t density:\n\t %s.rho\n",MFILE);
        fprintf(fp,"\t compressional wave velocities:\n\t %s.vp\n",MFILE);
        if(!ACOUSTIC) fprintf(fp,"\t tau for P-waves:\n\t %s.tp\n",MFILE);
        for (l=1;l<=L;l++) fprintf(fp,"\t %1i. relaxation frequencies: %s.f%1i\n",l,MFILE,l);
    }
    
    if(L){
        fprintf(fp,"\n");
        fprintf(fp," ------------------------- Q-APROXIMATION --------------------\n");
        fprintf(fp," Number of relaxation mechanisms (L): %i\n",L);
        fprintf(fp," The L relaxation frequencies are at:  \n");
        for (l=1;l<=L;l++) fprintf(fp,"\t%f",FL[l]);
        fprintf(fp," Hz\n");
        fprintf(fp," Value for tau is : %f\n",TAU);
        if (F_REF<0) fprintf(fp," Center frequency of source wavelet is used as reference frequency.\n");
        else fprintf(fp," Reference frequency: F_REF = %f\n",F_REF);
    }
    
    if (SNAP){
        fprintf(fp,"\n");
        fprintf(fp," -----------------------  SNAPSHOTS  -----------------------\n");
        fprintf(fp," Snapshots of");
        switch(SNAP){
            case 1:
                fprintf(fp," x- and y-component");
                fprintf(fp," of particle velocity.\n");
                break;
            case 2:
                fprintf(fp," pressure field.\n");
                break;
            case 3:
                fprintf(fp," curl and divergence energy of the wavefield.\n");
                break;
            case 4:
                fprintf(fp," curl and divergence energy of the wavefield.\n");
                fprintf(fp," x- and y-component of particle velocity.\n");
                break;
            default:
                declare_error(" sorry, incorrect value for SNAP ! \n");
        }
        
        fprintf(fp," \t first (TSNAP1)= %8.5f s\n", TSNAP1);
        fprintf(fp," \t last (TSNAP2)=%8.5f s\n",TSNAP2);
        fprintf(fp," \t increment (TSNAPINC) =%8.5f s\n\n",TSNAPINC);
        fprintf(fp," \t first_and_last_horizontal(x)_gridpoint = %i, %i \n",1,NX);
        fprintf(fp," \t first_and_last_vertical_gridpoint = %i, %i \n",1,NY);
        fprintf(fp," \n name of output-file (SNAP_FILE):\n\t %s\n",SNAP_FILE);
        switch (SNAP_FORMAT){
            case 1 :
                declare_error(" SU-Format not yet available !!");
                break;
            case 2 :
                fprintf(fp," The data is written in ASCII. \n");
                break;
            case 3 :
                fprintf(fp," The data is written binary (IEEE) (4 byte per float)");
                break;
            default:
                declare_error(" Don't know the format for the Snapshot-data ! \n");
        }
        
        fprintf(fp,"\n\n");
    }
    if (SEISMO){
        fprintf(fp,"\n");
        fprintf(fp," -----------------------  SEISMOGRAMS  ----------------------\n");
        if ((SEISMO==1) || (SEISMO==4) || (SEISMO==5)){
            fprintf(fp," seismograms of ");
            fprintf(fp," x- and y-component");
            fprintf(fp," of particle velocity.\n");
            fprintf(fp," output-files: \n ");
            fprintf(fp,"\t%s\n\t%s\n",SEIS_FILE,SEIS_FILE);
        }
        if ((SEISMO==2) || (SEISMO==4) || (SEISMO==5)){
            fprintf(fp," seismograms of pressure field (hydrophones).\n");
            fprintf(fp," output-file: \n ");
            fprintf(fp,"\t%s\n",SEIS_FILE);
        }
        if ((SEISMO==3) || (SEISMO==4)){
            fprintf(fp," seismograms of curl and div.\n");
            fprintf(fp," output-files: \n ");
            fprintf(fp,"\t%s\n\t%s\n",SEIS_FILE,SEIS_FILE);
            
        }
        
        switch (SEIS_FORMAT){
            case 1 :
                fprintf(fp," The data is written in IEEE SU-format . \n");
                break;
            case 2 :
                fprintf(fp," The data is written in ASCII. \n");
                break;
            case 3 :
                fprintf(fp," The data is written binary IEEE (4 byte per float)");
                break;
            default:
                declare_error(" Sorry. I don't know the format for the seismic data ! \n");
        }
        fprintf(fp," samplingrate of seismic data: %f s\n",NDT*DT);
        if (!READREC) fprintf(fp," Trace-spacing: %5.2f m\n", NGEOPH*DH);
        fprintf(fp," Number of samples per trace: %i \n", iround(NT/NDT));
        fprintf(fp," ----------------------------------------------------------\n");
        fprintf(fp,"\n");
        fprintf(fp,"\n");
    }
    fprintf(fp,"\n");
    fprintf(fp," -----------------------  IFOS elastic specific parameters  ----------------------\n");
    if (FORWARD_ONLY==1){
        fprintf(fp," FORWARD_ONLY=%d: Only forward modeling is applied.\n",FORWARD_ONLY);}
    else {
        if (FORWARD_ONLY==0){
            fprintf(fp," FORWARD_ONLY=%d: FWI is applied.\n",FORWARD_ONLY);}
        
        fprintf(fp,"\n Maximum number of iterations: %d\n",ITERMAX);
        fprintf(fp," location of the measured seismograms : \n ");
        fprintf(fp,"\t%s\n\n",DATA_DIR);
        if (PARAMETERIZATION==1){
            fprintf(fp," PARAMETERIZATION=%d: Inversion parameters are vp, vs and rho.\n",PARAMETERIZATION);}
        if (PARAMETERIZATION==2){
            fprintf(fp," PARAMETERIZATION=%d: Inversion parameters are Zp, Zs and rho.\n",PARAMETERIZATION);}
        if (PARAMETERIZATION==3){
            fprintf(fp," PARAMETERIZATION=%d: Inversion parameters are lambda, mu and rho.\n",PARAMETERIZATION);}
        fprintf(fp,"\n INVTYPE = %d\n\n",INVTYPE);
        if (ADJOINT_TYPE==1){
            fprintf(fp," ADJOINT_TYPE=%d: Inversion of x and y component.\n\n",ADJOINT_TYPE);}
        if (ADJOINT_TYPE==2){
            fprintf(fp," ADJOINT_TYPE=%d: Inversion of y component.\n\n",ADJOINT_TYPE);}
        if (ADJOINT_TYPE==3){
            fprintf(fp," ADJOINT_TYPE=%d: Inversion of x component.\n\n",ADJOINT_TYPE);}
        if (ADJOINT_TYPE==4){
            fprintf(fp," ADJOINT_TYPE=%d: Inversion of pressure component.\n\n",ADJOINT_TYPE);}
        
        if (VELOCITY==1){
            fprintf(fp," VELOCITY=%d: Minimization of misfit in velocity seismograms.\n\n",VELOCITY);}
        
        fprintf(fp," Shots used for step length estimation (in total %d testshots are used):\n",NO_OF_TESTSHOTS);
        fprintf(fp,"\t TESTSHOT_START = %d \n",TESTSHOT_START);
        fprintf(fp,"\t TESTSHOT_END = %d \n",TESTSHOT_END);
        fprintf(fp,"\t TESTSHOT_INCR = %d \n\n",TESTSHOT_INCR);
        
        fprintf(fp," Cosine Taper used : \n ");
        fprintf(fp,"\t%d\n",TAPER);
        
        fprintf(fp," Log file for misfit in each iteration step: \n ");
        fprintf(fp,"\t%s \n\n",MISFIT_LOG_FILE);
        
        fprintf(fp," Output of inverted models: \n ");
        fprintf(fp,"\t%s (nfstart=%d, nf=%d)\n\n",INV_MODELFILE,nfstart,nf);
        
        fprintf(fp," Output of jacobian: \n ");
        fprintf(fp,"\t%s (nfstart_jac=%d, nf_jac=%d)\n\n\n",JACOBIAN,nfstart_jac,nf_jac);
        
        fprintf(fp,"\n");
        fprintf(fp," Density is inverted from iteration step %d on.\n",INV_RHO_ITER);
        
        fprintf(fp,"\n");
        fprintf(fp," Vp is inverted from iteration step %d on.\n",INV_VP_ITER);
        
        if(!ACOUSTIC){
            fprintf(fp,"\n");
            fprintf(fp," Vs is inverted from iteration step %d on.\n\n\n",INV_VS_ITER);
        }
        
        fprintf(fp,"\n");
        fprintf(fp," Minimum Vp/Vs-ratio is set to %4.2f.\n",VP_VS_RATIO);
        if(VP_VS_RATIO<1)fprintf(fp," which means that it is disregarded. \n\n\n");
        
        if(S==1){
            fprintf(fp,"\n\n");
            fprintf(fp," Limited update of Vs in reference to the starting model is set to %4.2f %%.\n",S_VS);
            fprintf(fp," Limited update of Vp in reference to the starting model is set to %4.2f %%.\n",S_VP);
            fprintf(fp," Limited update of Rho in reference to the starting model is set to %4.2f %%.\n\n\n",S_RHO);
        }
        
        fprintf(fp," --------------- Gradient tapering -------------------\n");
        if (SWS_TAPER_GRAD_VERT==1){
            fprintf(fp," SWS_TAPER_GRAD_VERT=%d: Vertical taper applied.\n",SWS_TAPER_GRAD_VERT);
            fprintf(fp," (GRADT1=%d, GRADT2=%d, GRADT3=%d, GRADT4=%d)\n\n",GRADT1,GRADT2,GRADT3,GRADT4);}
        else	fprintf(fp," SWS_TAPER_GRAD_VERT=%d: No vertical taper applied.\n\n",SWS_TAPER_GRAD_VERT);
        
        if (SWS_TAPER_GRAD_HOR==1){
            fprintf(fp," SWS_TAPER_GRAD_HOR=%d: Horizontal taper applied.\n",SWS_TAPER_GRAD_HOR);
            fprintf(fp," (GRADT1=%d, GRADT2=%d, GRADT3=%d, GRADT4=%d)\n\n",GRADT1,GRADT2,GRADT3,GRADT4);}
        else	fprintf(fp," SWS_TAPER_GRAD_HOR=%d: No horizontal taper applied.\n\n",SWS_TAPER_GRAD_HOR);
        
        if (SWS_TAPER_GRAD_SOURCES==1){
            fprintf(fp," SWS_TAPER_GRAD_SOURCES=%d: Taper around the sources.\n",SWS_TAPER_GRAD_SOURCES);
            fprintf(fp," (SRTSHAPE=%d, SRTRADIUS=%f, FILTSIZE=%d)\n\n",SRTSHAPE,SRTRADIUS,FILTSIZE);}
        else	fprintf(fp," SWS_TAPER_GRAD_SOURCES=%d: No taper around the sources applied.\n\n",SWS_TAPER_GRAD_SOURCES);
        
        if (SWS_TAPER_CIRCULAR_PER_SHOT==1){
            fprintf(fp," SWS_TAPER_CIRCULAR_PER_SHOT=%d: Taper around the source for each shot.\n",SWS_TAPER_CIRCULAR_PER_SHOT);
            fprintf(fp," (SRTSHAPE=%d, SRTRADIUS=%f, FILTSIZE=%d)\n\n",SRTSHAPE,SRTRADIUS,FILTSIZE);}
        else	fprintf(fp," SWS_TAPER_CIRCULAR_PER_SHOT=%d: No taper around the sources applied.\n\n",SWS_TAPER_CIRCULAR_PER_SHOT);
        
        if (SWS_TAPER_FILE==1){
            fprintf(fp," SWS_TAPER_FILE=%d: Taper files taper.bin, taper_u.bin taper_rho.bin are read in and applied to the gradients.\n",SWS_TAPER_FILE);}
        else	fprintf(fp," SWS_TAPER_FILE=%d: No taper files are applied to the summed gradients.\n\n",SWS_TAPER_FILE);
        
        if (SWS_TAPER_FILE_PER_SHOT==1){
            fprintf(fp," SWS_TAPER_FILE_PER_SHOT=%d: Taper files for single shots are read in and applied to the gradients.\n",SWS_TAPER_FILE_PER_SHOT);
            fprintf(fp,"     File for vp or lambda gradients: %s.vp\n",TAPER_FILE_NAME);
            fprintf(fp,"     File for vs or mu gradients: %s.vs\n",TAPER_FILE_NAME);
            fprintf(fp,"     File for rho gradients: %s.rho\n",TAPER_FILE_NAME);}
        else	fprintf(fp," SWS_TAPER_FILE_PER_SHOT=%d: No taper files are applied to the gradients before summation.\n\n",SWS_TAPER_FILE_PER_SHOT);
        
        fprintf(fp,"\n");
        fprintf(fp," Smoothing (spatial filtering) of the gradients: \n ");
        if(SPATFILTER==1){
            fprintf(fp," \tSPATFILTER=%d: Gradients are smoothed.\n",SPATFILTER);
            fprintf(fp," \t(SPAT_FILT_SIZE=%d, SPAT_FILT_1=%d, SPAT_FILT_ITER=%d)\n",SPAT_FILT_SIZE,SPAT_FILT_1,SPAT_FILT_ITER);}
        else 	fprintf(fp," \tSPATFILTER=%d: Gradients are not smoothed.\n",SPATFILTER);
        
        fprintf(fp,"\n\n");
        fprintf(fp," --------------- Gradient smoothing with 2D-Gaussian filter -------------------\n");
        if(GRAD_FILTER==1){
            if(GRAD_FILT_WAVELENGTH==0)fprintf(fp," GRAD_FILTER=%d: Gradients are filtered.(FILT_SIZE_GRAD=%d)\n",GRAD_FILTER,FILT_SIZE_GRAD);
            if(GRAD_FILT_WAVELENGTH==1)fprintf(fp," GRAD_FILTER=%d: FILT_SIZE_GRAD is ignored. Gradients are filtered with a wavelength dependent filter size. Weighting factor A = %4.2f \n",GRAD_FILTER,A);}
        else 	fprintf(fp," GRAD_FILTER=%d: Jacobians are not filtered.\n",GRAD_FILTER);
        
        fprintf(fp,"\n\n");
        fprintf(fp," --------------- Limits of model parameters -------------------\n");
        fprintf(fp," VPLOWERLIM = %f \t\t VPUPPERLIM = %f \n",VPLOWERLIM,VPUPPERLIM);
        fprintf(fp," VSLOWERLIM = %f \t\t VSUPPERLIM = %f \n",VSLOWERLIM,VSUPPERLIM);
        fprintf(fp," RHOLOWERLIM = %f \t RHOUPPERLIM = %f \n",RHOLOWERLIM,RHOUPPERLIM);
        
        
        fprintf(fp,"\n\n");
        fprintf(fp," --------------- Calculation of the diagonal elements of the approximate Hessian matrix -------------------\n");
        switch(GRAD_METHOD){
            case 1:
                fprintf(fp," GRAD_METHOD=%d: PCG\n",GRAD_METHOD);
                break;
            case 2:
                fprintf(fp," GRAD_METHOD=%d: LBFGS\n",GRAD_METHOD);
                break;
            case 0: break;	/* only forward modeling is applied */
            default:
                declare_error(" Sorry, incorrect value for GRAD_METHOD ! \n");
        }
        
        
        fprintf(fp,"\n\n");
        fprintf(fp," --------------- Model smoothing -------------------\n");
        if(MODEL_FILTER==1){
            fprintf(fp," MODEL_FILTER=%d: vp and vs models are filtered after each iteration step.\n",MODEL_FILTER);
            fprintf(fp," (FILT_SIZE=%d)\n",FILT_SIZE);}
        else 	fprintf(fp," MODEL_FILTER=%d: vp and vs models are not filtered after each iteration step.\n",MODEL_FILTER);
        
        
        fprintf(fp,"\n\n");
        fprintf(fp," --------------- Inversion of the source time function -------------------\n");
        if(INV_STF==1){
            fprintf(fp," INV_STF=%d: Source time function will be inverted.\n",INV_STF);
            fprintf(fp," (PARA=%s, N_STF=%d, N_STF_START=%d)\n",PARA,N_STF,N_STF_START);}
        else 	fprintf(fp," INV_STF=%d: No inversion of the source time function.\n",INV_STF);
        
        
        
        fprintf(fp,"\n\n");
        
        fprintf(fp," --------------- Trace kill STF -------------------\n");
        if (TRKILL_STF){
            fprintf(fp," TRKILL_STF=%d: Trace kill STF is applied \n",TRKILL_STF);
            fprintf(fp," Reading trace kill STF matrix from file: %s \n\n",TRKILL_FILE_STF);}
        else fprintf(fp," TRKILL_STF=%d: No trace kill STF is applied \n",TRKILL_STF);
        
        
        
        fprintf(fp,"\n\n");
        fprintf(fp," --------------- Frequency filtering -------------------\n");
        if (TIME_FILT){
            if(TIME_FILT==1){
                fprintf(fp," TIME_FILT=%d: Time domain filtering is applied \n",TIME_FILT);
                fprintf(fp," Starting at frequencies of %.2f Hz\n",F_LOW_PASS_START);
                fprintf(fp," Increasing the bandwidth up to %.2f Hz in steps of %.2f Hz\n",F_LOW_PASS_END,F_LOW_PASS_INCR);
            }
            if(TIME_FILT==2){
                fprintf(fp," TIME_FILT=%d: Time domain filtering is applied \n Frequencies will be read from file: %s\n",TIME_FILT,FREQ_FILE);}
            fprintf(fp," Order of lowpass filter is:\t%d\n",ORDER);
            if ((ORDER%2)!=0){
                declare_error(" Order of time domain filter must be an even number! \n");}
        } else {
            fprintf(fp," TIME_FILT=%d: No time domain filtering is applied.\n",TIME_FILT);
        }
        
        
        fprintf(fp,"\n\n");
        fprintf(fp," --------------- Trace kill -------------------\n");
        if (TRKILL){
            fprintf(fp," TRKILL=%d: Trace kill is applied \n",TRKILL);
            if(TRKILL_OFFSET) {
                fprintf(fp," Traces with offset between %f m and %f m are killed.\n\n",TRKILL_OFFSET_LOWER,TRKILL_OFFSET_UPPER);
            } else {
                fprintf(fp," Reading trace kill matrix from file: %s \n\n",TRKILL_FILE);
            }
        } else {
            fprintf(fp," TRKILL=%d: No trace kill is applied \n",TRKILL);
        }
        
        
        
        fprintf(fp,"\n\n");
        fprintf(fp," --------------- Time windowing and damping -------------------\n");
        if (TIMEWIN){
            fprintf(fp," TIMEWIN=%d: Time windowing and damping is applied \n",TIMEWIN);
            fprintf(fp," Reading picked times from files: %s \n",PICKS_FILE);
            fprintf(fp," length of window after pick in s is: %f \n",TWLENGTH_PLUS);
            fprintf(fp," length of window befor pick in s is: %f \n",TWLENGTH_MINUS);
            fprintf(fp," gamma is : %f \n\n",GAMMA);}
        else fprintf(fp," TIMEWIN=%d: No time windowing and damping is applied \n",TIMEWIN);
        
        
        
        fprintf(fp,"\n\n");
        fprintf(fp," --------------- Trace normalization -------------------\n");
        if (NORMALIZE){
            fprintf(fp," NORMALIZE=%d: The measured and synthetic seismograms will be normalized.\n",NORMALIZE);
            fprintf(fp," before calculating the residuals. \n\n");}
        else fprintf(fp," NORMALIZE=%d: No normalization of measured and synthetic seismograms.\n",NORMALIZE);
        
        
        fprintf(fp,"\n\n");
        fprintf(fp," --------------- Gradient calculation -------------------\n");
        fprintf(fp," Used normed:\n");
        fprintf(fp,"   LNORM==1 corresponds to L1 Norm\n");
        fprintf(fp,"   LNORM==2 corresponds to L2 Norm\n");
        fprintf(fp,"   LNORM==3 corresponds to Cauchy\n");
        fprintf(fp,"   LNORM==4 corresponds to SECH\n");
        fprintf(fp,"   LNORM==5 corresponds to global correlation\n");
        fprintf(fp,"   LNORM==7 corresponds to normalized L2 Norm (each trace is normalized by its RMS value)\n");
        fprintf(fp,"   LNORM==8 corresponds to enveloped-based Norm\n\n");
        fprintf(fp," Switched LNORM=%d\n\n",LNORM);
        fprintf(fp,"   WATERLEVEL_LNORM8 is %e\n\n",WATERLEVEL_LNORM8);
        
        fprintf(fp," Every %d time sample is used for the calculation of the gradients.\n\n",DTINV);
        
        
        fprintf(fp,"\n\n");
        fprintf(fp," --------------- Step length estimation -------------------\n");
        fprintf(fp," EPS_SCALE = %f\n",EPS_SCALE);
        fprintf(fp," STEPMAX = %d\n",STEPMAX);
        fprintf(fp," SCALEFAC = %f\n",SCALEFAC);
        
        
        fprintf(fp,"\n\n");
        fprintf(fp," --------------- termination of the program -------------------\n");
        fprintf(fp," Misfit change during the last two iterations is smaller than %f percent.\n\n",(PRO*100.0));
        
        fprintf(fp,"\n\n");
        fprintf(fp," --------------- minimum number of iteration per frequency -------------------\n");
        fprintf(fp," MIN_ITER = %d \n\n",MIN_ITER);
    }
    
    fprintf(fp,"\n");
    fprintf(fp," **************************************************************\n");
    fprintf(fp,"\n");
}
