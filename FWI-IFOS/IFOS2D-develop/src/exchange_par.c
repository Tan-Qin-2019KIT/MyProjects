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
 *   Exchange FD-Parameters between PEs
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void exchange_par(void){
    
    /* declaration of extern variables */
    extern int   NX, NY, FDORDER, MAXRELERROR, SOURCE_SHAPE, SOURCE_TYPE, SNAP, SNAP_FORMAT, L;
    extern float DH, TIME, DT, TS, *FL, TAU, VPPML, PLANE_WAVE_DEPTH, PHI, F_REF;
    extern float XREC1, XREC2, YREC1, YREC2, FPML;
    extern float REC_ARRAY_DEPTH, REC_ARRAY_DIST, MUN, EPSILON, EPSILON_u, EPSILON_rho;
    extern int SEISMO, NDT, NGEOPH, SEIS_FORMAT, FREE_SURF, READMOD, READREC, SRCREC;
    extern int BOUNDARY, REC_ARRAY, DRX, FW, STF_FULL;
    extern int SNAPSHOT_START,SNAPSHOT_END,SNAPSHOT_INCR;
    extern float TSNAP1, TSNAP2, TSNAPINC, REFREC[4];
    extern char  MFILE[STRING_SIZE], SIGNAL_FILE[STRING_SIZE],SIGNAL_FILE_SH[STRING_SIZE], LOG_FILE[STRING_SIZE];
    extern char SNAP_FILE[STRING_SIZE], SOURCE_FILE[STRING_SIZE], REC_FILE[STRING_SIZE];
    extern char SEIS_FILE[STRING_SIZE];
    extern char JACOBIAN[STRING_SIZE], DATA_DIR[STRING_SIZE], INV_MODELFILE[STRING_SIZE], FREQ_FILE[STRING_SIZE];
    extern int RUN_MULTIPLE_SHOTS, TAPERLENGTH, INVTYPE;
    extern int NPROC, NPROCX, NPROCY, MYID, IDX, IDY;
    extern int GRADT1, GRADT2, GRADT3, GRADT4, ITERMAX, PARAMETERIZATION, FORWARD_ONLY, ADJOINT_TYPE;
    extern int GRAD_METHOD;
    extern float TSHIFT_back;
    extern int MODEL_FILTER, FILT_SIZE;
    extern int FILT_SIZE_GRAD, GRAD_FILTER;
    
    extern int TESTSHOT_START, TESTSHOT_END, TESTSHOT_INCR, NO_OF_TESTSHOTS;
    extern int SWS_TAPER_GRAD_VERT, SWS_TAPER_GRAD_HOR, SWS_TAPER_GRAD_SOURCES, SWS_TAPER_CIRCULAR_PER_SHOT, SRTSHAPE, FILTSIZE;
    extern int SWS_TAPER_FILE, SWS_TAPER_FILE_PER_SHOT;
    extern float SRTRADIUS;
    extern char TAPER_FILE_NAME[STRING_SIZE];
    extern int SPATFILTER, SPAT_FILT_SIZE, SPAT_FILT_1, SPAT_FILT_ITER;
    extern int INV_RHO_ITER, INV_VP_ITER, INV_VS_ITER;
    extern int MIN_ITER;
    extern int nfstart, nf;
    extern int nfstart_jac, nf_jac;
    extern float VPUPPERLIM, VPLOWERLIM, VSUPPERLIM, VSLOWERLIM, RHOUPPERLIM, RHOLOWERLIM;
    extern float npower, k_max_PML;
    extern int INV_STF, N_STF, N_STF_START;
    extern char PARA[STRING_SIZE];
    extern int TIME_FILT, ORDER,WRITE_FILTERED_DATA;
    extern float F_LOW_PASS_START, F_LOW_PASS_END, F_LOW_PASS_INCR, F_HIGH_PASS;
    extern int LNORM, DTINV;
    extern int STEPMAX;
    extern float EPS_SCALE, SCALEFAC;
    extern float PRO;
    extern int TRKILL, TRKILL_STF;
    extern char TRKILL_FILE[STRING_SIZE], TRKILL_FILE_STF[STRING_SIZE];
    extern int TAPER_STF;
    extern int TIMEWIN, NORMALIZE, TW_IND;
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
    extern int ACOUSTIC;
    extern int VERBOSE;
    
    extern int TRKILL_STF_OFFSET;
    extern int TRKILL_STF_OFFSET_INVERT;
    extern float TRKILL_STF_OFFSET_LOWER;
    extern float TRKILL_STF_OFFSET_UPPER;
    
    extern int TRKILL_OFFSET;
    extern float TRKILL_OFFSET_LOWER;
    extern float TRKILL_OFFSET_UPPER;
    
    // Parameter for inversion of SH waves and joint inversion
    extern int WAVETYPE;
    extern int SOURCE_SHAPE_SH;
    extern int JOINT_INVERSION_PSV_SH_TYPE;
    extern int JOINT_EQUAL_WEIGHTING;
    /* Workflow  */
    extern char FILE_WORKFLOW[STRING_SIZE];
    extern int USE_WORKFLOW;
    
    extern float JOINT_INVERSION_PSV_SH_ALPHA_VS;
    extern float JOINT_INVERSION_PSV_SH_ALPHA_RHO;
    
    extern int EPRECOND;
    extern int EPRECOND_ITER;
    extern float EPSILON_WE,EPSILON_WE_SH;
    extern int EPRECOND_PER_SHOT;
    extern int EPRECOND_PER_SHOT_SH;
    
    extern float LBFGS_SCALE_GRADIENTS;
    extern int LBFGS_STEP_LENGTH;
    
    extern int N_LBFGS;
    
    extern int WOLFE_CONDITION;
    extern int WOLFE_NUM_TEST;
    extern int WOLFE_TRY_OLD_STEPLENGTH;
    extern float WOLFE_C1_SL;
    extern float WOLFE_C2_SL;
    
    
    /* definition of local variables */
    /* NPAR is set in fd.h and must be set to the max. number of elements in idum/fdum +1 (because vector index starts at 0) */
    int idum[NPAR];
    float fdum[NPAR];
    
    
    if (MYID == 0){
        
        /***********/
        /*  Float  */
        /***********/
        fdum[1]  = DH;
        fdum[2]  = TIME;
        fdum[3]  = DT;
        fdum[4]  = TS;
        fdum[5]  = 0.0;
        fdum[6]  = 0.0;
        
        fdum[8]  = TAU;
        fdum[10]  = TSNAP1;
        fdum[11]  = TSNAP2;
        fdum[12]  = TSNAPINC;
        fdum[13]  = REFREC[1];
        fdum[14]  = REFREC[2];
        fdum[15]  = PHI;
        
        fdum[16]  = XREC1;
        fdum[17]  = YREC1;
        
        fdum[19]  = XREC2;
        fdum[20]  = YREC2;
        
        fdum[22]  = VPPML;
        fdum[23]  = REC_ARRAY_DEPTH;
        fdum[24]  = REC_ARRAY_DIST;
        fdum[25]  = PLANE_WAVE_DEPTH;
        
        fdum[26]  = MUN;
        fdum[27]  = EPSILON;
        fdum[28]  = EPSILON_u;
        fdum[29]  = EPSILON_rho;
        fdum[31]  = FPML;
        
        fdum[32]  = SRTRADIUS;
        
        fdum[33]  = VPUPPERLIM;
        fdum[34]  = VPLOWERLIM;
        fdum[35]  = VSUPPERLIM;
        fdum[36]  = VSLOWERLIM;
        fdum[37]  = RHOUPPERLIM;
        fdum[38]  = RHOLOWERLIM;
        
        fdum[39]  = npower;
        fdum[40]  = k_max_PML;

        fdum[42]  = F_LOW_PASS_START;
        fdum[43]  = F_LOW_PASS_END;
        fdum[44]  = F_LOW_PASS_INCR;
        
        fdum[45]  = EPS_SCALE;
        fdum[46]  = SCALEFAC;
        fdum[47]  = PRO;
        
        fdum[48]  = TWLENGTH_PLUS;
        fdum[49]  = TWLENGTH_MINUS;
        fdum[50]  = GAMMA;
        
        fdum[51]  = F_REF;
        
        fdum[52]  = TSHIFT_back;
        
        fdum[53]  = WATERLEVEL_LNORM8;
        
        fdum[54]  = VP_VS_RATIO;
        
        fdum[55]  = S_VS;
        
        fdum[56]  = S_VP;
        
        fdum[57]  = S_RHO;
        
        fdum[58]  = A;
        
        fdum[59]  = F_HIGH_PASS;
        
        fdum[60] = JOINT_INVERSION_PSV_SH_ALPHA_VS;
        fdum[61] = JOINT_INVERSION_PSV_SH_ALPHA_RHO;
        
        fdum[62]=EPSILON_WE;
        fdum[63]=EPSILON_WE_SH;
        
        fdum[64]=WOLFE_C1_SL;
        fdum[65]=WOLFE_C2_SL;
        
        fdum[66]=TRKILL_STF_OFFSET_LOWER;
        fdum[67]=TRKILL_STF_OFFSET_UPPER;
        fdum[68]=TRKILL_OFFSET_LOWER;
        fdum[69]=TRKILL_OFFSET_UPPER;
        
        fdum[70]=LBFGS_SCALE_GRADIENTS;

        /***********/
        /* Integer */
        /***********/
        idum[1]  = NPROCX;
        idum[2]  = NPROCY;
        
        idum[4]  = NPROCX*NPROCY;
        idum[5]  = NX;
        idum[6]  = NY;
        
        idum[8]  = SOURCE_SHAPE;
        idum[9]  = SOURCE_TYPE;
        idum[10]  = READMOD;
        idum[11]  = L;
        idum[12]  = FREE_SURF;
        idum[13]  = SNAP;
        idum[14]  = DRX;
        
        idum[16]  = BOUNDARY;
        idum[17]  = REC_ARRAY;
        idum[18]  = SRCREC;
        idum[19]  = IDX;
        idum[20]  = IDY;
        idum[21]  = TRKILL_STF;
        idum[22]  = 0;
        idum[23]  = SNAP_FORMAT;
        idum[24]  = SEISMO;
        idum[25]  = READREC;
        idum[26]  = NGEOPH;
        idum[27]  = NDT;
        idum[28]  = SEIS_FORMAT;

        
        idum[31]  = FDORDER;
        idum[32]  = MAXRELERROR;
        idum[33]  = RUN_MULTIPLE_SHOTS;
        idum[34]  = TAPERLENGTH;
        idum[35]  = INVTYPE;
        idum[36]  = GRADT1;
        idum[37]  = GRADT2;
        idum[38]  = GRADT3;
        idum[39]  = GRADT4;
        idum[40]  = ITERMAX;
        idum[41]  = PARAMETERIZATION;
        idum[42]  = FW;
        idum[43]  = FORWARD_ONLY;
        idum[44]  = ADJOINT_TYPE;
        
        idum[45]  = TESTSHOT_START;
        idum[46]  = TESTSHOT_END;
        idum[47]  = TESTSHOT_INCR;
        
        idum[48]  = SWS_TAPER_GRAD_VERT;
        idum[49]  = SWS_TAPER_GRAD_HOR;
        idum[50]  = SWS_TAPER_GRAD_SOURCES;
        idum[51]  = SWS_TAPER_CIRCULAR_PER_SHOT;
        idum[52]  = SRTSHAPE;
        idum[53]  = FILTSIZE;
        
        idum[54]  = SPATFILTER;
        idum[55]  = SPAT_FILT_SIZE;
        idum[56]  = SPAT_FILT_1;
        idum[57]  = SPAT_FILT_ITER;
        
        idum[58]  = INV_RHO_ITER;
        idum[59]  = nfstart;
        idum[60]  = nf;
        
        idum[61]  = nfstart_jac;
        idum[62]  = nf_jac;
        idum[63]  = SWS_TAPER_FILE;
        idum[65]  = GRAD_METHOD;
        
        idum[66]  = MODEL_FILTER;
        idum[67]  = FILT_SIZE;
        
        idum[69]  = INV_STF;
        idum[70]  = N_STF;
        idum[71]  = N_STF_START;
        
        idum[72]  = TIME_FILT;
        idum[73]  = ORDER;
        
        idum[74]  = LNORM;
        idum[75]  = DTINV;
        
        idum[76]  = STEPMAX;
        
        idum[77]  = TRKILL;
        
        idum[78]  = TIMEWIN;
        
        idum[79]  = NORMALIZE;
        
        idum[80]  = INV_VP_ITER;
        idum[81]  = INV_VS_ITER;
        
        idum[82]  = MIN_ITER;
        
        idum[83]  = GRAD_FILTER;
        idum[84]  = FILT_SIZE_GRAD;
        
        idum[85]  = NO_OF_TESTSHOTS;
        
        // idum[86]  = EMPTY;
        
        idum[87]  = VELOCITY;
        
        idum[88]  = SWS_TAPER_FILE_PER_SHOT;
        
        idum[89]  = S;
        
        idum[90]  = GRAD_FILT_WAVELENGTH;
        
        idum[91]  = ACOUSTIC;
        
        
        idum[92]  = VERBOSE;
        
        idum[93]  = WAVETYPE;
        
        idum[94]  = SOURCE_SHAPE_SH;
        
        idum[95] = JOINT_INVERSION_PSV_SH_TYPE;
        
        idum[96]  = SNAPSHOT_START;
        idum[97]  = SNAPSHOT_END;
        idum[98]  = SNAPSHOT_INCR;
        
        idum[100]=USE_WORKFLOW;
        
        idum[101]=EPRECOND;
        idum[102]=EPRECOND_ITER;
        
        idum[104]=LBFGS_STEP_LENGTH;
        
        idum[105]=EPRECOND_PER_SHOT;
        idum[106]=EPRECOND_PER_SHOT_SH;
        
        idum[107]=N_LBFGS;

        idum[108]=WOLFE_CONDITION;
        idum[109]=WOLFE_NUM_TEST;
        idum[110]=WOLFE_TRY_OLD_STEPLENGTH;
        
        idum[111]=WRITE_FILTERED_DATA;
        
        idum[112]=TAPER_STF;
        idum[113]=TW_IND;
        
        idum[114]=TRKILL_OFFSET;
        idum[115]=TRKILL_STF_OFFSET;
        idum[116]=TRKILL_STF_OFFSET_INVERT;
        
        idum[117]=JOINT_EQUAL_WEIGHTING;
        idum[118]=STF_FULL;
    } /** if (MYID == 0) **/
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    MPI_Bcast(&idum,NPAR,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&fdum,NPAR,MPI_FLOAT,0,MPI_COMM_WORLD);
    
    
    MPI_Bcast(&FILE_WORKFLOW,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
    MPI_Bcast(&SOURCE_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
    MPI_Bcast(&MFILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
    MPI_Bcast(&SNAP_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
    MPI_Bcast(&REC_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
    MPI_Bcast(&SEIS_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
    MPI_Bcast(&LOG_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
    MPI_Bcast(&SIGNAL_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
    MPI_Bcast(&SIGNAL_FILE_SH,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
    MPI_Bcast(&JACOBIAN,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
    MPI_Bcast(&DATA_DIR,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
    MPI_Bcast(&INV_MODELFILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
    MPI_Bcast(&FREQ_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
    
    
    MPI_Bcast(&PARA,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
    
    MPI_Bcast(&MISFIT_LOG_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
    
    MPI_Bcast(&TRKILL_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
    MPI_Bcast(&TRKILL_FILE_STF,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
    
    MPI_Bcast(&PICKS_FILE,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
    
    MPI_Bcast(&TAPER_FILE_NAME,STRING_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    /***********/
    /*  Float  */
    /***********/
    DH=fdum[1];
    TIME=fdum[2];
    DT=fdum[3];
    TS=fdum[4];
    
    TAU=fdum[8];
    TSNAP1=fdum[10];
    TSNAP2=fdum[11];
    TSNAPINC=fdum[12];
    REFREC[1]=fdum[13];
    REFREC[2]=fdum[14];
    PHI=fdum[15];
    XREC1=fdum[16];
    YREC1=fdum[17];
    
    XREC2=fdum[19];
    YREC2=fdum[20];
    
    VPPML=fdum[22];
    REC_ARRAY_DEPTH=fdum[23];
    REC_ARRAY_DIST=fdum[24];
    PLANE_WAVE_DEPTH=fdum[25];
    
    MUN = fdum[26];
    EPSILON = fdum[27];
    EPSILON_u = fdum[28];
    EPSILON_rho = fdum[29];
    FPML = fdum[31];
    
    SRTRADIUS = fdum[32];
    
    VPUPPERLIM = fdum[33];
    VPLOWERLIM = fdum[34];
    VSUPPERLIM = fdum[35];
    VSLOWERLIM = fdum[36];
    RHOUPPERLIM = fdum[37];
    RHOLOWERLIM = fdum[38];
    
    npower = fdum[39];
    k_max_PML = fdum[40];
    
    
    F_LOW_PASS_START = fdum[42];
    F_LOW_PASS_END = fdum[43];
    F_LOW_PASS_INCR = fdum[44];
    
    EPS_SCALE = fdum[45];
    SCALEFAC = fdum[46];
    
    PRO = fdum[47];
    
    TWLENGTH_PLUS = fdum[48];
    TWLENGTH_MINUS = fdum[49];
    GAMMA = fdum[50];
    
    F_REF = fdum[51];
    
    TSHIFT_back = fdum[52];
    
    WATERLEVEL_LNORM8 = fdum[53];
    
    VP_VS_RATIO = fdum[54];
    
    S_VS = fdum[55];
    
    S_VP = fdum[56];
    
    S_RHO = fdum[57];
    
    A = fdum[58];
    
    F_HIGH_PASS = fdum[59];
    
    JOINT_INVERSION_PSV_SH_ALPHA_VS = fdum[60];
    JOINT_INVERSION_PSV_SH_ALPHA_RHO = fdum[61];
    
    EPSILON_WE=fdum[62];
    EPSILON_WE_SH=fdum[63];
    
    WOLFE_C1_SL=fdum[64];
    WOLFE_C2_SL=fdum[65];
    
    TRKILL_STF_OFFSET_LOWER=fdum[66];
    TRKILL_STF_OFFSET_UPPER=fdum[67];
    TRKILL_OFFSET_LOWER=fdum[68];
    TRKILL_OFFSET_UPPER=fdum[69];
    
    LBFGS_SCALE_GRADIENTS=fdum[70];
    
    /***********/
    /* Integer */
    /***********/
    
    NPROCX = idum[1];
    NPROCY = idum[2];
    NPROC  = idum[4];
    NX = idum[5];
    NY = idum[6];
    
    SOURCE_SHAPE = idum[8];
    SOURCE_TYPE = idum[9];
    READMOD = idum[10];
    L = idum[11];
    FREE_SURF = idum[12];
    SNAP = idum[13];
    DRX = idum[14];
    
    BOUNDARY = idum[16];
    REC_ARRAY = idum[17];
    SRCREC = idum[18];
    IDX = idum[19];
    IDY = idum[20];
    TRKILL_STF=idum[21];
    
    SNAP_FORMAT = idum[23];
    SEISMO = idum[24];
    READREC = idum[25];
    NGEOPH = idum[26];
    NDT = idum[27];
    SEIS_FORMAT = idum[28];

    
    FDORDER = idum[31];
    MAXRELERROR = idum[32];
    RUN_MULTIPLE_SHOTS = idum[33];
    TAPERLENGTH = idum[34];
    INVTYPE = idum[35];
    GRADT1 = idum[36];
    GRADT2 = idum[37];
    GRADT3 = idum[38];
    GRADT4 = idum[39];
    ITERMAX = idum[40];
    PARAMETERIZATION = idum[41];
    FW = idum[42];
    FORWARD_ONLY  = idum[43];  
    ADJOINT_TYPE = idum[44];
    
    TESTSHOT_START = idum[45];
    TESTSHOT_END = idum[46];
    TESTSHOT_INCR = idum[47];
    
    SWS_TAPER_GRAD_VERT = idum[48];
    SWS_TAPER_GRAD_HOR = idum[49];
    SWS_TAPER_GRAD_SOURCES = idum[50];
    SWS_TAPER_CIRCULAR_PER_SHOT = idum[51];
    SRTSHAPE = idum[52];
    FILTSIZE = idum[53];
    
    SPATFILTER = idum[54];
    SPAT_FILT_SIZE = idum[55];
    SPAT_FILT_1 = idum[56];
    SPAT_FILT_ITER = idum[57];
    
    INV_RHO_ITER = idum[58];
    nfstart = idum[59];
    nf = idum[60];
    
    nfstart_jac = idum[61];
    nf_jac = idum[62];
    SWS_TAPER_FILE = idum[63];
    GRAD_METHOD = idum[65];
    
    
    MODEL_FILTER = idum[66];
    FILT_SIZE = idum[67];
    
    
    INV_STF = idum[69];
    N_STF = idum[70];
    N_STF_START = idum[71];
    
    TIME_FILT = idum[72];
    ORDER = idum[73];
    
    LNORM = idum[74];
    DTINV = idum[75];
    
    STEPMAX = idum[76];
    
    TRKILL = idum[77];
    
    TIMEWIN = idum[78];
    
    NORMALIZE = idum[79];
    
    INV_VP_ITER = idum[80];
    INV_VS_ITER = idum[81];
    
    MIN_ITER = idum[82];
    
    GRAD_FILTER = idum[83];
    FILT_SIZE_GRAD = idum[84];
    
    NO_OF_TESTSHOTS = idum[85];
    
    // EMPTY = idum[86];
    
    VELOCITY = idum[87];
    
    SWS_TAPER_FILE_PER_SHOT = idum[88];
    
    S = idum[89];
    
    GRAD_FILT_WAVELENGTH = idum[90];
    
    ACOUSTIC = idum[91];
    
    VERBOSE = idum[92];
    
    WAVETYPE = idum[93];
    
    SOURCE_SHAPE_SH = idum[94];
    
    JOINT_INVERSION_PSV_SH_TYPE = idum[95];
    
    SNAPSHOT_START=idum[96];
    SNAPSHOT_END=idum[97];
    SNAPSHOT_INCR=idum[98];
    
    USE_WORKFLOW=idum[100];
    
    EPRECOND=idum[101];
    EPRECOND_ITER=idum[102];
    
    LBFGS_STEP_LENGTH=idum[104];
    
    EPRECOND_PER_SHOT= idum[105];
    EPRECOND_PER_SHOT_SH= idum[106];
    
    N_LBFGS=idum[107];
    
    WOLFE_CONDITION=idum[108];
    
    WOLFE_NUM_TEST=idum[109];
    
    WOLFE_TRY_OLD_STEPLENGTH=idum[110];
    
    WRITE_FILTERED_DATA=idum[111];
    
    TAPER_STF=idum[112];
    TW_IND=idum[113];
    
    TRKILL_OFFSET=idum[114];
    TRKILL_STF_OFFSET=idum[115];
    TRKILL_STF_OFFSET_INVERT=idum[116];
    
    JOINT_EQUAL_WEIGHTING=idum[117];
    STF_FULL=idum[118];
    if ( MYID!=0 && L>0 ) {
        FL=vector(1,L);
    }
    
    if ( L>0 ) {
        MPI_Bcast(&FL[1],L,MPI_FLOAT,0,MPI_COMM_WORLD);
    }

}
