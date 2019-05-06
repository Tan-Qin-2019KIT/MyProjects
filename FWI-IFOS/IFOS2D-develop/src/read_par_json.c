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
 *   program IFOS, reading input-parameters from input-file or stdin
 *  ----------------------------------------------------------------------*/

#include <unistd.h>
#include "fd.h"

char ** varname_list,** value_list;

void read_par_json(FILE *fp, char *fileinp){
    
    /* declaration of extern variables */
    extern int   NX, NY, FDORDER, MAXRELERROR, SOURCE_SHAPE,SOURCE_SHAPE_SH, SOURCE_TYPE, SNAP, SNAP_FORMAT, ACOUSTIC, L, VERBOSE, WAVETYPE,JOINT_INVERSION_PSV_SH_TYPE,JOINT_EQUAL_WEIGHTING;
    extern float DH, TIME, DT, TS, *FL, TAU, VPPML, PLANE_WAVE_DEPTH, PHI, F_REF,JOINT_INVERSION_PSV_SH_ALPHA_VS,JOINT_INVERSION_PSV_SH_ALPHA_RHO;
    extern float XREC1, XREC2, YREC1, YREC2, FPML;
    extern float REC_ARRAY_DEPTH, REC_ARRAY_DIST;
    extern int SEISMO, NDT, NGEOPH, SEIS_FORMAT, FREE_SURF, READMOD, READREC, SRCREC, RUN_MULTIPLE_SHOTS;
    extern int BOUNDARY, REC_ARRAY, DRX, TAPER, TAPERLENGTH, INVTYPE, FW;
    extern float TSNAP1, TSNAP2, TSNAPINC, REFREC[4];
    extern int SNAPSHOT_START,SNAPSHOT_END,SNAPSHOT_INCR;
    extern char  MFILE[STRING_SIZE], SIGNAL_FILE[STRING_SIZE], SIGNAL_FILE_SH[STRING_SIZE];
    extern char SNAP_FILE[STRING_SIZE], SOURCE_FILE[STRING_SIZE], REC_FILE[STRING_SIZE];
    extern char SEIS_FILE[STRING_SIZE];
    extern char JACOBIAN[STRING_SIZE],DATA_DIR[STRING_SIZE],FREQ_FILE[STRING_SIZE];
    extern int  NPROCX, NPROCY, MYID, IDX, IDY;
    extern int GRADT1, GRADT2, GRADT3, GRADT4, ITERMAX, PARAMETERIZATION, FORWARD_ONLY, ADJOINT_TYPE;
    extern int  GRAD_METHOD;
    extern int FILT_SIZE, MODEL_FILTER;
    extern int FILT_SIZE_GRAD, GRAD_FILTER;
    
    extern int TESTSHOT_START, TESTSHOT_END, TESTSHOT_INCR, NO_OF_TESTSHOTS;
    extern int SWS_TAPER_GRAD_VERT, SWS_TAPER_GRAD_HOR, SWS_TAPER_GRAD_SOURCES, SWS_TAPER_CIRCULAR_PER_SHOT, SRTSHAPE, FILTSIZE;
    extern int SWS_TAPER_FILE, SWS_TAPER_FILE_PER_SHOT;
    extern float SRTRADIUS;
    extern char TAPER_FILE_NAME[STRING_SIZE];
    extern int SPATFILTER, SPAT_FILT_SIZE, SPAT_FILT_1, SPAT_FILT_ITER;
    extern int INV_RHO_ITER, INV_VS_ITER, INV_VP_ITER;
    extern char INV_MODELFILE[STRING_SIZE];
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
    extern int MIN_ITER;
    
    extern int TRKILL;
    extern char TRKILL_FILE[STRING_SIZE];
    
    extern int TRKILL_OFFSET;
    extern float TRKILL_OFFSET_LOWER;
    extern float TRKILL_OFFSET_UPPER;
    
    extern int TRKILL_STF;
    extern char TRKILL_FILE_STF[STRING_SIZE];

    extern int TRKILL_STF_OFFSET;
    extern int TRKILL_STF_OFFSET_INVERT;
    extern float TRKILL_STF_OFFSET_LOWER;
    extern float TRKILL_STF_OFFSET_UPPER;
    
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
    
    extern char FILE_WORKFLOW[STRING_SIZE];
    extern int USE_WORKFLOW;
    
    extern int EPRECOND;
    extern int EPRECOND_ITER;
    extern float EPSILON_WE,EPSILON_WE_SH;
    extern int EPRECOND_PER_SHOT;
    extern int EPRECOND_PER_SHOT_SH;
    
    extern int LBFGS_STEP_LENGTH;
    extern int N_LBFGS;
    
    extern float LBFGS_SCALE_GRADIENTS;
    extern int WOLFE_CONDITION;
    extern int WOLFE_NUM_TEST;
    extern int WOLFE_TRY_OLD_STEPLENGTH;
    extern float WOLFE_C1_SL;
    extern float WOLFE_C2_SL;
    
    extern int STF_FULL;
    /* definition of local variables */
    
    int number_readobjects=0,fserr=0;
    
    
    /* Support naming of variables  */
    float FC_START, FC_END, FC_INCR, F_HP;
    
    
    char errormessage[STRING_SIZE2];
    
    char ** varname_list, ** value_list;
    
    if (MYID == 0){
        
        /* allocate first object in list */
        varname_list = malloc(STRING_SIZE2*sizeof(char*));
        value_list = malloc(STRING_SIZE2*sizeof(char*));
        
        /* read in objects from file */
        number_readobjects=read_objects_from_intputfile(fp, fileinp, varname_list, value_list);
        fprintf(fp,"\nFrom input file %s, %i objects have been read in. \n",fileinp, number_readobjects);
        
        /* print objects to screen */
        fprintf(fp, "\n===========================================================");
        fprintf(fp, "\n=   List of Parameters read by the built in Json Parser   =");
        print_objectlist_screen(fp, number_readobjects, varname_list, value_list);
        
        /* extract variables form object list */
        
        /*=================================
         section general grid and discretization parameters
         =================================*/
        if (get_int_from_objectlist("NPROCX",number_readobjects,&NPROCX,varname_list, value_list))
            declare_error("Variable NPROCX could not be retrieved from the json input file!");
        if (get_int_from_objectlist("NPROCY",number_readobjects,&NPROCY,varname_list, value_list))
            declare_error("Variable NPROCY could not be retrieved from the json input file!");
        if (get_int_from_objectlist("FDORDER",number_readobjects,&FDORDER,varname_list, value_list))
            declare_error("Variable FDORDER could not be retrieved from the json input file!");
        if (get_int_from_objectlist("MAXRELERROR",number_readobjects,&MAXRELERROR,varname_list, value_list))
            declare_error("Variable MAXRELERROR could not be retrieved from the json input file!");
        if (get_int_from_objectlist("NX",number_readobjects,&NX,varname_list, value_list))
            declare_error("Variable NX could not be retrieved from the json input file!");
        if (get_int_from_objectlist("NY",number_readobjects,&NY,varname_list, value_list))
            declare_error("Variable NY could not be retrieved from the json input file!");
        if (get_float_from_objectlist("DH",number_readobjects,&DH,varname_list, value_list))
            declare_error("Variable DH could not be retrieved from the json input file!");
        if (get_float_from_objectlist("TIME",number_readobjects,&TIME,varname_list, value_list))
            declare_error("Variable TIME could not be retrieved from the json input file!");
        if (get_float_from_objectlist("DT",number_readobjects,&DT,varname_list, value_list))
            declare_error("Variable DT could not be retrieved from the json input file!");
        
        /*=================================
         section source parameters
         =================================*/
        fprintf(fp,"The following default values are set:\n");
        fprintf(fp,"=====================================\n\n");
        
        if (get_int_from_objectlist("SOURCE_TYPE",number_readobjects,&SOURCE_TYPE,varname_list, value_list))
            declare_error("Variable SOURCE_TYPE could not be retrieved from the json input file!");
        
        /* Definition of inversion for source time function */
        if (get_int_from_objectlist("INV_STF",number_readobjects,&INV_STF,varname_list, value_list)){
            INV_STF=0;
            fprintf(fp,"Variable INV_STF is set to default value %d.\n",INV_STF);}
        else {
            if (INV_STF==1) {
                if (get_string_from_objectlist("PARA",number_readobjects,PARA,varname_list, value_list))
                    declare_error("Variable PARA could not be retrieved from the json input file!");
                if (get_int_from_objectlist("N_STF",number_readobjects,&N_STF,varname_list, value_list))
                    declare_error("Variable N_STF could not be retrieved from the json input file!");
                if (get_int_from_objectlist("N_STF_START",number_readobjects,&N_STF_START,varname_list, value_list))
                    declare_error("Variable N_STF_START could not be retrieved from the json input file!");
            }
        }
        
        if (get_int_from_objectlist("ACOUSTIC",number_readobjects,&ACOUSTIC,varname_list, value_list)){
            ACOUSTIC=0;
            fprintf(fp,"Variable ACOUSTIC is set to default value %d.\n",ACOUSTIC);}
        if (get_int_from_objectlist("WAVETYPE",number_readobjects,&WAVETYPE,varname_list, value_list)){
            WAVETYPE=1;
            fprintf(fp,"Variable WAVETYPE is set to default value %d.\n",WAVETYPE);
        } else {
            if (ACOUSTIC && WAVETYPE!=1) {
                WAVETYPE=1;
                fprintf(fp,"For acoustic modelling WAVETYPE is set to %d.\n",WAVETYPE);
            }
            if(WAVETYPE==3) {
                
                if (get_int_from_objectlist("JOINT_EQUAL_WEIGHTING",number_readobjects,&JOINT_EQUAL_WEIGHTING,varname_list, value_list)){
                    JOINT_EQUAL_WEIGHTING=0;
                    fprintf(fp,"Variable JOINT_EQUAL_WEIGHTING is set to default value %d.\n",JOINT_EQUAL_WEIGHTING);
                }
                
                if (get_int_from_objectlist("JOINT_INVERSION_PSV_SH_TYPE",number_readobjects,&JOINT_INVERSION_PSV_SH_TYPE,varname_list, value_list)){
                    JOINT_INVERSION_PSV_SH_TYPE=1;
                    fprintf(fp,"Variable JOINT_INVERSION_PSV_SH_TYPE is set to default value %d.\n",JOINT_INVERSION_PSV_SH_TYPE);
                } else {
                    /* Check herer possible dependencies */
                }
                if (get_float_from_objectlist("JOINT_INVERSION_PSV_SH_ALPHA_VS",number_readobjects,&JOINT_INVERSION_PSV_SH_ALPHA_VS,varname_list, value_list)){
                    JOINT_INVERSION_PSV_SH_ALPHA_VS=0.5;
                    fprintf(fp,"Variable JOINT_INVERSION_PSV_SH_ALPHA_VS is set to default value %f.\n",JOINT_INVERSION_PSV_SH_ALPHA_VS);
                } else {
                    /* Check herer possible dependencies */
                }
                if (get_float_from_objectlist("JOINT_INVERSION_PSV_SH_ALPHA_RHO",number_readobjects,&JOINT_INVERSION_PSV_SH_ALPHA_RHO,varname_list, value_list)){
                    JOINT_INVERSION_PSV_SH_ALPHA_RHO=0.5;
                    fprintf(fp,"Variable JOINT_INVERSION_PSV_SH_ALPHA_RHO is set to default value %f.\n",JOINT_INVERSION_PSV_SH_ALPHA_RHO);
                } else {
                    /* Check herer possible dependencies */
                }
            }
        }
        
        if (get_int_from_objectlist("SOURCE_SHAPE",number_readobjects,&SOURCE_SHAPE,varname_list, value_list))
            declare_error("Variable SOURCE_SHAPE could not be retrieved from the json input file!");
        else {
            if ((SOURCE_SHAPE==3)||(SOURCE_SHAPE==7)||(INV_STF==1)) {
                if (WAVETYPE==1 || WAVETYPE==3){
                    if (get_string_from_objectlist("SIGNAL_FILE",number_readobjects,SIGNAL_FILE,varname_list, value_list))
                        declare_error("Variable SIGNAL_FILE could not be retrieved from the json input file!");
                }
                if (WAVETYPE==2 || WAVETYPE==3){
                    if (get_string_from_objectlist("SIGNAL_FILE_SH",number_readobjects,SIGNAL_FILE_SH,varname_list, value_list))
                        declare_error("Variable SIGNAL_FILE_SH could not be retrieved from the json input file!");
                }
            }
        }
        if (get_int_from_objectlist("SOURCE_SHAPE_SH",number_readobjects,&SOURCE_SHAPE_SH,varname_list, value_list)){
            if (WAVETYPE==2 || WAVETYPE==3){
                declare_error("Variable SOURCE_SHAPE_SH could not be retrieved from the json input file!");
            }
        }
        if (get_int_from_objectlist("SRCREC",number_readobjects,&SRCREC,varname_list, value_list)){
            SRCREC=1;
            fprintf(fp,"Variable SRCREC is set to default value %d\n",SRCREC);}
        
        if (SRCREC==1) {
            if (get_string_from_objectlist("SOURCE_FILE",number_readobjects,SOURCE_FILE,varname_list, value_list))
                declare_error("Variable SOURCE_FILE could not be retrieved from the json input file!");
            if (get_int_from_objectlist("RUN_MULTIPLE_SHOTS",number_readobjects,&RUN_MULTIPLE_SHOTS,varname_list, value_list))
                declare_error("Variable RUN_MULTIPLE_SHOTS could not be retrieved from the json input file!");
        }
        if (SRCREC==2) {
            if (get_float_from_objectlist("PLANE_WAVE_DEPTH",number_readobjects,&PLANE_WAVE_DEPTH,varname_list, value_list))
                declare_error("Variable PLANE_WAVE_DEPTH could not be retrieved from the json input file!");
            else {
                if (PLANE_WAVE_DEPTH>0.0) {
                    if (get_float_from_objectlist("PHI",number_readobjects,&PHI,varname_list, value_list))
                        declare_error("Variable PHI could not be retrieved from the json input file!");
                    if (get_float_from_objectlist("TS",number_readobjects,&TS,varname_list, value_list))
                        declare_error("Variable TS could not be retrieved from the json input file!");
                }
            }
        }
        
        
        
        /*=================================
         section boundary parameters
         =================================*/
        
        if (get_int_from_objectlist("FREE_SURF",number_readobjects,&FREE_SURF,varname_list, value_list))
            declare_error("Variable FREE_SURF could not be retrieved from the json input file!");
        if (get_int_from_objectlist("BOUNDARY",number_readobjects,&BOUNDARY,varname_list, value_list))
            declare_error("Variable BOUNDARY could not be retrieved from the json input file!");
        if (get_int_from_objectlist("FW",number_readobjects,&FW,varname_list, value_list))
            declare_error("Variable FW could not be retrieved from the json input file!");
        if (get_float_from_objectlist("VPPML",number_readobjects,&VPPML,varname_list, value_list))
            declare_error("Variable VPPML could not be retrieved from the json input file!");
        if (get_float_from_objectlist("FPML",number_readobjects,&FPML,varname_list, value_list))
            declare_error("Variable FPML could not be retrieved from the json input file!");
        if (get_float_from_objectlist("npower",number_readobjects,&npower,varname_list, value_list))
            declare_error("Variable npower could not be retrieved from the json input file!");
        if (get_float_from_objectlist("k_max_PML",number_readobjects,&k_max_PML,varname_list, value_list))
            declare_error("Variable k_max_PML could not be retrieved from the json input file!");
        
        
        /*=================================
         section snapshot parameters
         =================================*/
        if (get_int_from_objectlist("SNAP",number_readobjects,&SNAP,varname_list, value_list)){
            SNAP=0;
            fprintf(fp,"Variable SNAP is set to default value %d.\n",SNAP);}
        else {
            if (SNAP>0) {
                if (get_int_from_objectlist("SNAP_FORMAT",number_readobjects,&SNAP_FORMAT,varname_list, value_list))
                    declare_error("Variable SNAP_FORMAT could not be retrieved from the json input file!");
                if (get_float_from_objectlist("TSNAP1",number_readobjects,&TSNAP1,varname_list, value_list))
                    declare_error("Variable TSNAP1 could not be retrieved from the json input file!");
                if (get_float_from_objectlist("TSNAP2",number_readobjects,&TSNAP2,varname_list, value_list))
                    declare_error("Variable TSNAP2 could not be retrieved from the json input file!");
                if (get_float_from_objectlist("TSNAPINC",number_readobjects,&TSNAPINC,varname_list, value_list))
                    declare_error("Variable TSNAPINC could not be retrieved from the json input file!");
                if (get_string_from_objectlist("SNAP_FILE",number_readobjects,SNAP_FILE,varname_list, value_list))
                    declare_error("Variable SNAP_FILE could not be retrieved from the json input file!");
                if (get_int_from_objectlist("SNAPSHOT_START",number_readobjects,&SNAPSHOT_START,varname_list, value_list))
                    declare_error("Variable SNAPSHOT_START could not be retrieved from the json input file!");
                if (get_int_from_objectlist("SNAPSHOT_END",number_readobjects,&SNAPSHOT_END,varname_list, value_list))
                    declare_error("Variable SNAPSHOT_START could not be retrieved from the json input file!");
                if (get_int_from_objectlist("SNAPSHOT_INCR",number_readobjects,&SNAPSHOT_INCR,varname_list, value_list))
                    declare_error("Variable SNAPSHOT_INCR could not be retrieved from the json input file!");
            }
        }
        /* increments are read in any case, because they will be also used as increment for model output */
        if (get_int_from_objectlist("IDX",number_readobjects,&IDX,varname_list, value_list)){
            IDX=1;
            fprintf(fp,"Variable IDX is set to default value %d.\n",IDX);}
        if (get_int_from_objectlist("IDY",number_readobjects,&IDY,varname_list, value_list)){
            IDY=1;
            fprintf(fp,"Variable IDY is set to default value %d.\n",IDY);}
        
        /*=================================
         section seismogramm parameters
         =================================*/
        if (get_int_from_objectlist("SEISMO",number_readobjects,&SEISMO,varname_list, value_list))
            declare_error("Variable SEISMO could not be retrieved from the json input file!");
        else {
            if (SEISMO>0){
                if (get_string_from_objectlist("SEIS_FILE",number_readobjects,SEIS_FILE,varname_list, value_list))
                    declare_error("Variable SEIS_FILE could not be retrieved from the json input file!");
                
                if (get_int_from_objectlist("READREC",number_readobjects,&READREC,varname_list, value_list))
                    declare_error("Variable READREC could not be retrieved from the json input file!");
                else {
                    if (READREC==0 || READREC==2) {
                        if (get_float_from_objectlist("XREC1",number_readobjects,&XREC1,varname_list, value_list))
                            declare_error("Variable XREC1 could not be retrieved from the json input file!");
                        if (get_float_from_objectlist("XREC2",number_readobjects,&XREC2,varname_list, value_list))
                            declare_error("Variable XREC2T could not be retrieved from the json input file!");
                        if (get_float_from_objectlist("YREC1",number_readobjects,&YREC1,varname_list, value_list))
                            declare_error("Variable YREC1 could not be retrieved from the json input file!");
                        if (get_float_from_objectlist("YREC2",number_readobjects,&YREC2,varname_list, value_list))
                            declare_error("Variable YREC2 could not be retrieved from the json input file!");
                        if (get_int_from_objectlist("NGEOPH",number_readobjects,&NGEOPH,varname_list, value_list))
                            declare_error("Variable NGEOPH could not be retrieved from the json input file!");
                    }
                    if (READREC>0) {
                        if (get_string_from_objectlist("REC_FILE",number_readobjects,REC_FILE,varname_list, value_list))
                            declare_error("Variable REC_FILE could not be retrieved from the json input file!");
                    }
                }
                if (get_int_from_objectlist("NDT",number_readobjects,&NDT,varname_list, value_list)){
                    NDT=1;
                    fprintf(fp,"Variable NDT is set to default value %d.\n",NDT);}
                if (get_int_from_objectlist("SEIS_FORMAT",number_readobjects,&SEIS_FORMAT,varname_list, value_list))
                    declare_error("Variable SEIS_FORMAT could not be retrieved from the json input file!");
                
                if (get_int_from_objectlist("REC_ARRAY",number_readobjects,&REC_ARRAY,varname_list, value_list)){
                    REC_ARRAY=0;
                    fprintf(fp,"Variable REC_ARRAY is set to default value %d.\n",REC_ARRAY);}
                else {
                    if (REC_ARRAY>0) {
                        if (get_int_from_objectlist("DRX",number_readobjects,&DRX,varname_list, value_list))
                            declare_error("Variable DRX could not be retrieved from the json input file!");
                        if (get_float_from_objectlist("REC_ARRAY_DEPTH",number_readobjects,&REC_ARRAY_DEPTH,varname_list, value_list))
                            declare_error("Variable REC_ARRAY_DEPTH could not be retrieved from the json input file!");
                        if (get_float_from_objectlist("REC_ARRAY_DIST",number_readobjects,&REC_ARRAY_DIST,varname_list, value_list))
                            declare_error("Variable REC_ARRAY_DIST could not be retrieved from the json input file!");
                    }
                }
                if (get_float_from_objectlist("REFRECX",number_readobjects,&REFREC[1],varname_list, value_list))
                    declare_error("Variable REFRECX could not be retrieved from the json input file!");
                if (get_float_from_objectlist("REFRECY",number_readobjects,&REFREC[2],varname_list, value_list))
                    declare_error("Variable REFRECY could not be retrieved from the json input file!");
            }
        }
        
        /*=================================
         section general model and log parameters
         =================================*/
        if (get_int_from_objectlist("VERBOSE",number_readobjects,&VERBOSE,varname_list, value_list)){
            VERBOSE=1;
            fprintf(fp,"Variable VERBOSE is set to value %d.\n",VERBOSE);
        } else {
        }
        if (get_string_from_objectlist("MFILE",number_readobjects,MFILE,varname_list, value_list))
            declare_error("Variable MFILE could not be retrieved from the json input file!");
        if (get_int_from_objectlist("READMOD",number_readobjects,&READMOD,varname_list, value_list))
            declare_error("Variable READMOD could not be retrieved from the json input file!");
        
        if (get_int_from_objectlist("L",number_readobjects,&L,varname_list, value_list)){
            L=0;
            F_REF=0;
            fprintf(fp,"Variable L is set to default value %d.\n",L);}
        else {
            FL=vector(1,L);
            switch(L) {
                case 0:
                    break;
                case 5:
                    if (get_float_from_objectlist("FL5",number_readobjects,&FL[5],varname_list, value_list))
                        declare_error("Variable FL5 could not be retrieved from the json input file!");
                case 4:
                    if (get_float_from_objectlist("FL4",number_readobjects,&FL[4],varname_list, value_list))
                        declare_error("Variable FL4 could not be retrieved from the json input file!");
                case 3:
                    if (get_float_from_objectlist("FL3",number_readobjects,&FL[3],varname_list, value_list))
                        declare_error("Variable FL3 could not be retrieved from the json input file!");
                case 2:
                    if (get_float_from_objectlist("FL2",number_readobjects,&FL[2],varname_list, value_list))
                        declare_error("Variable FL2 could not be retrieved from the json input file!");
                case 1:
                    if (get_float_from_objectlist("FL1",number_readobjects,&FL[1],varname_list, value_list))
                        declare_error("Variable FL1 could not be retrieved from the json input file!");
                    break;
                default:
                    declare_error("More than four relaxation Parameter (L>5) are not implemented yet!");
                    break;
            }
            if (get_float_from_objectlist("TAU",number_readobjects,&TAU,varname_list, value_list))
                declare_error("Variable TAU could not be retrieved from the json input file!");
            if (get_float_from_objectlist("F_REF",number_readobjects,&F_REF,varname_list, value_list)){
                F_REF=-1.0;
                fprintf(fp,"Reference frequency for viscoelastic modeling is set to center frequency of the source wavelet.\n");}
        }
        
        
        /*=================================
         section inversion parameters
         =================================*/
        
        if (get_int_from_objectlist("PARAMETERIZATION",number_readobjects,&PARAMETERIZATION,varname_list, value_list))
            declare_error("Variable PARAMETERIZATION could not be retrieved from the json input file!");
        else
            if(ACOUSTIC){
                PARAMETERIZATION=1;
                fprintf(fp,"For acoustic modelling only PARAMETERIZATION=%d possible.\n",PARAMETERIZATION);}
        
        if (get_int_from_objectlist("FORWARD_ONLY",number_readobjects,&FORWARD_ONLY,varname_list, value_list))
            declare_error("Variable FORWARD_ONLY could not be retrieved from the json input file!");
        else {
            if (FORWARD_ONLY==0) {	/* FWI is calculated */
                /* Overwrite IDX/IDY option from forward modeling (used for snapshots), interpolation for FWI not yet implemented*/
                IDX=1;
                IDY=1;
                /* General inversion parameters */
                if (get_int_from_objectlist("ITERMAX",number_readobjects,&ITERMAX,varname_list, value_list))
                    declare_error("Variable ITERMAX could not be retrieved from the json input file!");
                if (get_string_from_objectlist("DATA_DIR",number_readobjects,DATA_DIR,varname_list, value_list))
                    declare_error("Variable DATA_DIR could not be retrieved from the json input file!");
                if (get_int_from_objectlist("INVTYPE",number_readobjects,&INVTYPE,varname_list, value_list)){
                    INVTYPE=2;
                    fprintf(fp,"\nVariable INVTYPE is set to default value %d.\n",INVTYPE);}
                if (get_int_from_objectlist("ADJOINT_TYPE",number_readobjects,&ADJOINT_TYPE,varname_list, value_list))
                    declare_error("Variable ADJOINT_TYPE could not be retrieved from the json input file!");
                
                if (get_string_from_objectlist("MISFIT_LOG_FILE",number_readobjects,MISFIT_LOG_FILE,varname_list, value_list)){
                    strcpy(MISFIT_LOG_FILE,"L2_LOG.dat");
                    fprintf(fp,"\nVariable MISFIT_LOG_FILE is set to default value %s.\n",MISFIT_LOG_FILE);}
                
                if (get_int_from_objectlist("VELOCITY",number_readobjects,&VELOCITY,varname_list, value_list)){
                    VELOCITY=0;
                    fprintf(fp,"Variable VELOCITY is set to default value %d.\n",VELOCITY);
                }
                
                if (get_int_from_objectlist("USE_WORKFLOW",number_readobjects,&USE_WORKFLOW,varname_list, value_list)){
                    USE_WORKFLOW=0;
                } else {
                    
                    if (get_string_from_objectlist("FILE_WORKFLOW",number_readobjects,FILE_WORKFLOW,varname_list, value_list))
                        declare_error("Variable FILE_WORKFLOW could not be retrieved from the json input file!");
                }
                
                if (get_int_from_objectlist("EPRECOND",number_readobjects,&EPRECOND,varname_list, value_list)){
                    EPRECOND=0;
                    fprintf(fp,"Variable EPRECOND is set to default value %d.\n",EPRECOND);
                } else {
                    if (get_int_from_objectlist("EPRECOND_ITER",number_readobjects,&EPRECOND_ITER,varname_list, value_list)){
                        EPRECOND_ITER=0;
                        fprintf(fp,"Variable EPRECOND_ITER is set to default value %d.\n",EPRECOND_ITER);
                    }
                    if (get_int_from_objectlist("EPRECOND_PER_SHOT",number_readobjects,&EPRECOND_PER_SHOT,varname_list, value_list)){
                        EPRECOND_PER_SHOT=0;
                        fprintf(fp,"Variable EPRECOND_PER_SHOT is set to default value %d.\n",EPRECOND_PER_SHOT);
                        if(EPRECOND_ITER!=0) {
                            EPRECOND_ITER=0;
                            fprintf(fp," EPRECOND_PER_SHOT and EPRECOND_ITER>0 not supported.\n");
                            fprintf(fp," EPRECOND_ITER is set to EPRECOND_ITER=%d.\n",EPRECOND_ITER);
                        }
                    }
                    if (get_int_from_objectlist("EPRECOND_PER_SHOT_SH",number_readobjects,&EPRECOND_PER_SHOT_SH,varname_list, value_list)){
                        EPRECOND_PER_SHOT_SH=0;
                        fprintf(fp,"Variable EPRECOND_PER_SHOT_SH is set to default value %d.\n",EPRECOND_PER_SHOT_SH);
                    }
                    if (get_float_from_objectlist("EPSILON_WE",number_readobjects,&EPSILON_WE,varname_list, value_list))
                        declare_error("Variable EPSILON_WE could not be retrieved from the json input file!");
                    
                    if (get_float_from_objectlist("EPSILON_WE_SH",number_readobjects,&EPSILON_WE_SH,varname_list, value_list)) {
                        EPSILON_WE_SH=EPSILON_WE;
                        fprintf(fp,"Variable EPSILON_WE_SH is set to EPSILON_WE=%f.\n",EPSILON_WE_SH);
                    }
                    
                }
                
                if (get_int_from_objectlist("TESTSHOT_START",number_readobjects,&TESTSHOT_START,varname_list, value_list))
                    declare_error("Variable TESTSHOT_START could not be retrieved from the json input file!");
                if (get_int_from_objectlist("TESTSHOT_END",number_readobjects,&TESTSHOT_END,varname_list, value_list))
                    declare_error("Variable TESTSHOT_START could not be retrieved from the json input file!");
                if (get_int_from_objectlist("TESTSHOT_INCR",number_readobjects,&TESTSHOT_INCR,varname_list, value_list))
                    declare_error("Variable TESTSHOT_INCR could not be retrieved from the json input file!");
                NO_OF_TESTSHOTS=(TESTSHOT_END-TESTSHOT_START)/TESTSHOT_INCR+1;	/* calculation of number of testsshots */
                
                
                /* Cosine taper */
                if (get_int_from_objectlist("TAPER",number_readobjects,&TAPER,varname_list, value_list)){
                    TAPER=0;
                    fprintf(fp,"Variable TAPER is set to default value %d.\n",TAPER);}
                if (get_int_from_objectlist("TAPERLENGTH",number_readobjects,&TAPERLENGTH,varname_list, value_list)){
                    TAPERLENGTH=1000;
                    fprintf(fp,"Variable TAPERLENGTH is set to default value %d.\n",TAPERLENGTH);}
                
                
                /* Definition of gradient taper geometry */
                if (get_int_from_objectlist("SWS_TAPER_GRAD_VERT",number_readobjects,&SWS_TAPER_GRAD_VERT,varname_list, value_list)){
                    SWS_TAPER_GRAD_VERT=0;
                    fprintf(fp,"Variable SWS_TAPER_GRAD_VERT is set to default value %d.\n",SWS_TAPER_GRAD_VERT);}
                if (get_int_from_objectlist("SWS_TAPER_GRAD_HOR",number_readobjects,&SWS_TAPER_GRAD_HOR,varname_list, value_list)){
                    SWS_TAPER_GRAD_HOR=0;
                    fprintf(fp,"Variable SWS_TAPER_GRAD_HOR is set to default value %d.\n",SWS_TAPER_GRAD_HOR);}
                if ((SWS_TAPER_GRAD_VERT==1) || (SWS_TAPER_GRAD_HOR==1)){
                    if (get_int_from_objectlist("GRADT1",number_readobjects,&GRADT1,varname_list, value_list))
                        declare_error("Variable GRADT1 could not be retrieved from the json input file!");
                    if (get_int_from_objectlist("GRADT2",number_readobjects,&GRADT2,varname_list, value_list))
                        declare_error("Variable GRADT2 could not be retrieved from the json input file!");
                    if (get_int_from_objectlist("GRADT3",number_readobjects,&GRADT3,varname_list, value_list))
                        declare_error("Variable GRADT3 could not be retrieved from the json input file!");
                    if (get_int_from_objectlist("GRADT4",number_readobjects,&GRADT4,varname_list, value_list))
                        declare_error("Variable GRADT4 could not be retrieved from the json input file!");
                }
                if (get_int_from_objectlist("SWS_TAPER_GRAD_SOURCES",number_readobjects,&SWS_TAPER_GRAD_SOURCES,varname_list, value_list)){
                    SWS_TAPER_GRAD_SOURCES=0;
                    fprintf(fp,"Variable SWS_TAPER_GRAD_SOURCES is set to default value %d.\n",SWS_TAPER_GRAD_SOURCES);}
                if (get_int_from_objectlist("SWS_TAPER_CIRCULAR_PER_SHOT",number_readobjects,&SWS_TAPER_CIRCULAR_PER_SHOT,varname_list, value_list)){
                    SWS_TAPER_CIRCULAR_PER_SHOT=0;
                    fprintf(fp,"Variable SWS_TAPER_CIRCULAR_PER_SHOT is set to default value %d.\n",SWS_TAPER_CIRCULAR_PER_SHOT);}
                if ((SWS_TAPER_GRAD_SOURCES==1) || (SWS_TAPER_CIRCULAR_PER_SHOT==1)){
                    if (get_int_from_objectlist("SRTSHAPE",number_readobjects,&SRTSHAPE,varname_list, value_list))
                        declare_error("Variable SRTSHAPE could not be retrieved from the json input file!");
                    if (get_float_from_objectlist("SRTRADIUS",number_readobjects,&SRTRADIUS,varname_list, value_list))
                        declare_error("Variable SRTRADIUS could not be retrieved from the json input file!");
                    if (get_int_from_objectlist("FILTSIZE",number_readobjects,&FILTSIZE,varname_list, value_list))
                        declare_error("Variable FILTSIZE could not be retrieved from the json input file!");
                }
                if (get_int_from_objectlist("SWS_TAPER_FILE",number_readobjects,&SWS_TAPER_FILE,varname_list, value_list)){
                    SWS_TAPER_FILE=0;
                    fprintf(fp,"Variable SWS_TAPER_FILE is set to default value %d.\n",SWS_TAPER_FILE);
                }
                if (SWS_TAPER_FILE==1){
                    if (get_string_from_objectlist("TAPER_FILE_NAME",number_readobjects,TAPER_FILE_NAME,varname_list, value_list))
                        declare_error("Variable TAPER_FILE_NAME could not be retrieved from the json input file!");
                }
                if (get_int_from_objectlist("SWS_TAPER_FILE_PER_SHOT",number_readobjects,&SWS_TAPER_FILE_PER_SHOT,varname_list, value_list)){
                    SWS_TAPER_FILE_PER_SHOT=0;
                    fprintf(fp,"Variable SWS_TAPER_FILE_PER_SHOT is set to default value %d.\n",SWS_TAPER_FILE_PER_SHOT);
                }
                if (SWS_TAPER_FILE_PER_SHOT==1){
                    if (get_string_from_objectlist("TAPER_FILE_NAME",number_readobjects,TAPER_FILE_NAME,varname_list, value_list))
                        declare_error("Variable TAPER_FILE_NAME could not be retrieved from the json input file!");
                }
                
                
                
                /* Definition of smoothing (spatial filtering) of the gradients */
                if (get_int_from_objectlist("SPATFILTER",number_readobjects,&SPATFILTER,varname_list, value_list)){
                    SPATFILTER=0;
                    fprintf(fp,"Variable SPATFILTER is set to default value %d.\n",SPATFILTER);}
                
                else {
                    if (SPATFILTER==1) {
                        if (get_int_from_objectlist("SPAT_FILT_SIZE",number_readobjects,&SPAT_FILT_SIZE,varname_list, value_list))
                            declare_error("Variable SPAT_FILT_SIZE could not be retrieved from the json input file!");
                        if (get_int_from_objectlist("SPAT_FILT_1",number_readobjects,&SPAT_FILT_1,varname_list, value_list))
                            declare_error("Variable SPAT_FILT_1 could not be retrieved from the json input file!");
                        if (get_int_from_objectlist("SPAT_FILT_ITER",number_readobjects,&SPAT_FILT_ITER,varname_list, value_list))
                            declare_error("Variable SPAT_FILT_ITER could not be retrieved from the json input file!");
                    }
                }
                
                /* Definition of 2D-Gaussian filter of the gradients */
                if (get_int_from_objectlist("GRAD_FILTER",number_readobjects,&GRAD_FILTER,varname_list, value_list)){
                    GRAD_FILTER=0;
                    fprintf(fp,"Variable GRAD_FILTER is set to default value %d.\n",GRAD_FILTER);}
                
                if (get_int_from_objectlist("FILT_SIZE_GRAD",number_readobjects,&FILT_SIZE_GRAD,varname_list, value_list)){
                    FILT_SIZE_GRAD=0;
                    fprintf(fp,"Variable FILT_SIZE_GRAD is set to default value %d.\n",FILT_SIZE_GRAD);}
                
                if (get_int_from_objectlist("GRAD_FILT_WAVELENGTH",number_readobjects,&GRAD_FILT_WAVELENGTH,varname_list, value_list)){
                    GRAD_FILT_WAVELENGTH=0;
                    fprintf(fp,"Variable GRAD_FILT_WAVELENGTH is set to default value %d.\n",GRAD_FILT_WAVELENGTH);}
                else{
                    if (GRAD_FILT_WAVELENGTH==1) {
                        if (get_float_from_objectlist("A",number_readobjects,&A,varname_list, value_list))
                            declare_error("Variable A could not be retrieved from the json input file!");
                    }
                }
                
                
                /* Output of inverted models */
                if (get_string_from_objectlist("INV_MODELFILE",number_readobjects,INV_MODELFILE,varname_list, value_list))
                    declare_error("Variable INV_MODELFILE could not be retrieved from the json input file!");
                if (get_int_from_objectlist("nfstart",number_readobjects,&nfstart,varname_list, value_list))
                    declare_error("Variable nfstart could not be retrieved from the json input file!");
                if (get_int_from_objectlist("nf",number_readobjects,&nf,varname_list, value_list))
                    declare_error("Variable nf could not be retrieved from the json input file!");
                
                
                /* Output of gradients */
                if (get_string_from_objectlist("JACOBIAN",number_readobjects,JACOBIAN,varname_list, value_list))
                    declare_error("Variable JACOBIAN could not be retrieved from the json input file!");
                if (get_int_from_objectlist("nfstart_jac",number_readobjects,&nfstart_jac,varname_list, value_list))
                    declare_error("Variable nfstart_jac could not be retrieved from the json input file!");
                if (get_int_from_objectlist("nf_jac",number_readobjects,&nf_jac,varname_list, value_list))
                    declare_error("Variable nf_jac could not be retrieved from the json input file!");
                
                
                /* Inversion for density */
                if (get_int_from_objectlist("INV_RHO_ITER",number_readobjects,&INV_RHO_ITER,varname_list, value_list)){
                    INV_RHO_ITER=0;
                    fprintf(fp,"Variable INV_RHO_ITER is set to default value %d.\n",INV_RHO_ITER);}
                
                /* Inversion for Vp */
                if (get_int_from_objectlist("INV_VP_ITER",number_readobjects,&INV_VP_ITER,varname_list, value_list)){
                    INV_VP_ITER=0;
                    fprintf(fp,"Variable INV_VP_ITER is set to default value %d.\n",INV_VP_ITER);}
                
                /* Inversion for Vs */
                if (get_int_from_objectlist("INV_VS_ITER",number_readobjects,&INV_VS_ITER,varname_list, value_list)){
                    INV_VS_ITER=0;
                    fprintf(fp,"Variable INV_VS_ITER is set to default value %d.\n",INV_VS_ITER);}
                
                /* Vp/Vs-Ratio */
                if (get_float_from_objectlist("VP_VS_RATIO",number_readobjects,&VP_VS_RATIO,varname_list, value_list)){
                    VP_VS_RATIO=0.0;
                    fprintf(fp,"Variable VP_VS_RATIO is set to default value %4.2f which means that it is disregarded.\n",VP_VS_RATIO);}
                
                /* Limited update of model parameters in reference to the starting model */
                if (get_int_from_objectlist("S",number_readobjects,&S,varname_list, value_list)){
                    S=0;
                    fprintf(fp,"Variable S is set to default value %d.\n",S);}
                else {
                    if (S==1) {
                        if (get_float_from_objectlist("S_VS",number_readobjects,&S_VS,varname_list, value_list))
                            declare_error("Variable S_VS could not be retrieved from the json input file!");
                        if (get_float_from_objectlist("S_VP",number_readobjects,&S_VP,varname_list, value_list))
                            declare_error("Variable S_VP could not be retrieved from the json input file!");
                        if (get_float_from_objectlist("S_RHO",number_readobjects,&S_RHO,varname_list, value_list))
                            declare_error("Variable S_RHO could not be retrieved from the json input file!");
                    }
                }
                
                
                /* Upper and lower limits for model parameters */
                if (get_float_from_objectlist("VPUPPERLIM",number_readobjects,&VPUPPERLIM,varname_list, value_list)){
                    VPUPPERLIM=25000.0;
                    fprintf(fp,"Variable VPUPPERLIM is set to default value %f.\n",VPUPPERLIM);}
                if (get_float_from_objectlist("VPLOWERLIM",number_readobjects,&VPLOWERLIM,varname_list, value_list)){
                    VPLOWERLIM=0.0;
                    fprintf(fp,"Variable VPLOWERLIM is set to default value %f.\n",VPLOWERLIM);}
                if (get_float_from_objectlist("VSUPPERLIM",number_readobjects,&VSUPPERLIM,varname_list, value_list)){
                    VSUPPERLIM=25000.0;
                    fprintf(fp,"Variable VSUPPERLIM is set to default value %f.\n",VSUPPERLIM);}
                if (get_float_from_objectlist("VSLOWERLIM",number_readobjects,&VSLOWERLIM,varname_list, value_list)){
                    VSLOWERLIM=0.0;
                    fprintf(fp,"Variable VSLOWERLIM is set to default value %f.\n",VSLOWERLIM);}
                if (get_float_from_objectlist("RHOUPPERLIM",number_readobjects,&RHOUPPERLIM,varname_list, value_list)){
                    RHOUPPERLIM=25000.0;
                    fprintf(fp,"Variable RHOUPPERLIM is set to default value %f.\n",RHOUPPERLIM);}
                if (get_float_from_objectlist("RHOLOWERLIM",number_readobjects,&RHOLOWERLIM,varname_list, value_list)){
                    RHOLOWERLIM=0.0;
                    fprintf(fp,"Variable RHOLOWERLIM is set to default value %f.\n",RHOLOWERLIM);}
                
                
                /* Hessian and Gradient-Method */
                if (get_int_from_objectlist("GRAD_METHOD",number_readobjects,&GRAD_METHOD,varname_list, value_list))
                    declare_error("Variable GRAD_METHOD could not be retrieved from the json input file!");
                else {
                    if(GRAD_METHOD==1) {
                        WOLFE_CONDITION=0;
                    }
                    if(GRAD_METHOD==2) {
                        if (get_int_from_objectlist("LBFGS_STEP_LENGTH",number_readobjects,&LBFGS_STEP_LENGTH,varname_list, value_list)){
                            LBFGS_STEP_LENGTH=1;
                            fprintf(fp,"Variable LBFGS_STEP_LENGTH is set to default value %d.\n",LBFGS_STEP_LENGTH);
                        }
                        if (get_int_from_objectlist("N_LBFGS",number_readobjects,&N_LBFGS,varname_list, value_list)){
                            N_LBFGS=5;
                            fprintf(fp,"Variable N_LBFGS is set to default value %d.\n",N_LBFGS);
                        }
                        
                        if (get_float_from_objectlist("LBFGS_SCALE_GRADIENTS",number_readobjects,&LBFGS_SCALE_GRADIENTS,varname_list, value_list)){
                            LBFGS_SCALE_GRADIENTS=1;
                        }
                        
                        if (get_int_from_objectlist("WOLFE_CONDITION",number_readobjects,&WOLFE_CONDITION,varname_list, value_list)){
                            WOLFE_CONDITION=1;
                            fprintf(fp,"Variable WOLFE_CONDITION is set to default value %d.\n",WOLFE_CONDITION);
                        } else {
                            if (get_int_from_objectlist("WOLFE_NUM_TEST",number_readobjects,&WOLFE_NUM_TEST,varname_list, value_list)){
                                WOLFE_NUM_TEST=10;
                                fprintf(fp,"Variable WOLFE_NUM_TEST is set to default value %d.\n",WOLFE_NUM_TEST);
                            }
                            if (get_int_from_objectlist("WOLFE_TRY_OLD_STEPLENGTH",number_readobjects,&WOLFE_TRY_OLD_STEPLENGTH,varname_list, value_list)){
                                WOLFE_TRY_OLD_STEPLENGTH=0;
                                fprintf(fp,"Variable WOLFE_TRY_OLD_STEPLENGTH is set to default value %d.\n",WOLFE_TRY_OLD_STEPLENGTH);
                            }
                            if (get_float_from_objectlist("WOLFE_C1_SL",number_readobjects,&WOLFE_C1_SL,varname_list, value_list)){
                                WOLFE_C1_SL=1e-4;
                                fprintf(fp,"Variable WOLFE_C1_SL is set to default value %f.\n",WOLFE_C1_SL);
                            }
                            if (get_float_from_objectlist("WOLFE_C2_SL",number_readobjects,&WOLFE_C2_SL,varname_list, value_list)){
                                WOLFE_C2_SL=0.9;
                                fprintf(fp,"Variable WOLFE_C2_SL is set to default value %f.\n",WOLFE_C2_SL);
                            }
                        }
                    }
                }
                
                /* Definition of smoothing the models vp and vs */
                if (get_int_from_objectlist("MODEL_FILTER",number_readobjects,&MODEL_FILTER,varname_list, value_list)){
                    MODEL_FILTER=0;
                    fprintf(fp,"Variable MODEL_FILTER is set to default value %d.\n",MODEL_FILTER);}
                else {
                    if (MODEL_FILTER==1) {
                        if (get_int_from_objectlist("FILT_SIZE",number_readobjects,&FILT_SIZE,varname_list, value_list))
                            declare_error("Variable FILT_SIZE could not be retrieved from the json input file!");
                    }
                }
                
                if(INV_STF){
                    /* Trace killing STF */
                    if (get_int_from_objectlist("TRKILL_STF",number_readobjects,&TRKILL_STF,varname_list, value_list)){
                        TRKILL_STF=0;
                        fprintf(fp,"Variable TRKILL_STF is set to default value %d.\n",TRKILL_STF);}
                    else {
                        if (get_int_from_objectlist("STF_FULL",number_readobjects,&STF_FULL,varname_list, value_list)){
                            STF_FULL=0;
                            fprintf(fp,"Variable STF_FULL is set to default value %d.\n",STF_FULL);}
                        if (TRKILL_STF==1) {
                            if (get_int_from_objectlist("TRKILL_STF_OFFSET",number_readobjects,&TRKILL_STF_OFFSET,varname_list, value_list)){
                                TRKILL_STF_OFFSET=0;
                                if (get_string_from_objectlist("TRKILL_FILE_STF",number_readobjects,TRKILL_FILE_STF,varname_list, value_list))
                                    declare_error("Variable TRKILL_FILE_STF could not be retrieved from the json input file!");
                            } else {
                                if(TRKILL_STF_OFFSET==0) { /* Only TRKILL File */
                                    if (get_string_from_objectlist("TRKILL_FILE_STF",number_readobjects,TRKILL_FILE_STF,varname_list, value_list))
                                        declare_error("Variable TRKILL_FILE_STF could not be retrieved from the json input file!");
                                }
                                if(TRKILL_STF_OFFSET==1) { /* Only Offset based TRKill */
                                    if (get_int_from_objectlist("TRKILL_STF_OFFSET_INVERT",number_readobjects,&TRKILL_STF_OFFSET_INVERT,varname_list, value_list)){
                                        TRKILL_STF_OFFSET_INVERT=0;
                                    }
                                    
                                    if (get_float_from_objectlist("TRKILL_STF_OFFSET_LOWER",number_readobjects,&TRKILL_STF_OFFSET_LOWER,varname_list, value_list)){
                                        TRKILL_STF_OFFSET_LOWER=0.0;
                                    }
                                    if (get_float_from_objectlist("TRKILL_STF_OFFSET_UPPER",number_readobjects,&TRKILL_STF_OFFSET_UPPER,varname_list, value_list)){
                                        declare_error("Variable TRKILL_STF_OFFSET_UPPER could not be retrieved from the json input file!");
                                    }
                                }
                                if(TRKILL_STF_OFFSET==2){ /* Both Offset based TRKill & File */
                                    if (get_int_from_objectlist("TRKILL_STF_OFFSET_INVERT",number_readobjects,&TRKILL_STF_OFFSET_INVERT,varname_list, value_list)){
                                        TRKILL_STF_OFFSET_INVERT=0;
                                    } else {
                                        if(TRKILL_STF_OFFSET_INVERT==1) {
                                            declare_error("Variable TRKILL_STF_OFFSET_INVERT==1 and TRKILL_STF_OFFSET==2 not possible!");
                                        }
                                    }
                                    if (get_int_from_objectlist("TRKILL_STF_OFFSET_INVERT",number_readobjects,&TRKILL_STF_OFFSET_INVERT,varname_list, value_list)){
                                        TRKILL_STF_OFFSET_INVERT=0;
                                    }
                                    
                                    if (get_float_from_objectlist("TRKILL_STF_OFFSET_LOWER",number_readobjects,&TRKILL_STF_OFFSET_LOWER,varname_list, value_list)){
                                        TRKILL_STF_OFFSET_LOWER=0.0;
                                    }
                                    if (get_float_from_objectlist("TRKILL_STF_OFFSET_UPPER",number_readobjects,&TRKILL_STF_OFFSET_UPPER,varname_list, value_list)){
                                        declare_error("Variable TRKILL_STF_OFFSET_UPPER could not be retrieved from the json input file!");
                                    }
                                    if (get_string_from_objectlist("TRKILL_FILE_STF",number_readobjects,TRKILL_FILE_STF,varname_list, value_list))
                                        declare_error("Variable TRKILL_FILE_STF could not be retrieved from the json input file!");
                                }
                            }
                        }
                    }
                    /* Taper STF */
                    if (get_int_from_objectlist("TAPER_STF",number_readobjects,&TAPER_STF,varname_list, value_list)){
                        TAPER_STF=0;
                        fprintf(fp,"Variable TAPER_STF is set to default value %d.\n",TAPER_STF);}
                }
                
                /* Frequency filtering during inversion */
                if (get_int_from_objectlist("TIME_FILT",number_readobjects,&TIME_FILT,varname_list, value_list)){
                    TIME_FILT=0;
                    fprintf(fp,"Variable TIME_FILT is set to default value %d.\n",TIME_FILT);}
                else {
                    if (get_int_from_objectlist("WRITE_FILTERED_DATA",number_readobjects,&WRITE_FILTERED_DATA,varname_list, value_list)){
                        WRITE_FILTERED_DATA=0;
                    }
                    if (get_float_from_objectlist("F_HIGH_PASS",number_readobjects,&F_HIGH_PASS,varname_list, value_list)){
                        /* Support of old variable naming: Test if old variable naming is used */
                        if (get_float_from_objectlist("F_HP",number_readobjects,&F_HIGH_PASS,varname_list, value_list)){
                            F_HIGH_PASS=0.0;
                            fprintf(fp,"Variable F_HIGH_PASS is set to default value %f.\n",F_HIGH_PASS);
                        }
                    }
                    if (TIME_FILT==1) {
                        if (get_float_from_objectlist("F_LOW_PASS_START",number_readobjects,&F_LOW_PASS_START,varname_list, value_list)){
                            /* Support of old variable naming: Test if old variable naming is used */
                            if (get_float_from_objectlist("FC_START",number_readobjects,&F_LOW_PASS_START,varname_list, value_list)){
                                declare_error("Variable F_LOW_PASS_START could not be retrieved from the json input file!");
                            }
                        }
                        if (get_float_from_objectlist("F_LOW_PASS_END",number_readobjects,&F_LOW_PASS_END,varname_list, value_list)){
                            /* Support of old variable naming: Test if old variable naming is used */
                            if (get_float_from_objectlist("FC_END",number_readobjects,&F_LOW_PASS_END,varname_list, value_list)){
                                declare_error("Variable F_LOW_PASS_END could not be retrieved from the json input file!");
                            }
                        }
                        if (get_float_from_objectlist("F_LOW_PASS_INCR",number_readobjects,&F_LOW_PASS_INCR,varname_list, value_list)){
                            /* Support of old variable naming: Test if old variable naming is used */
                            if (get_float_from_objectlist("FC_INCR",number_readobjects,&F_LOW_PASS_INCR,varname_list, value_list)){
                                declare_error("Variable F_LOW_PASS_INCR could not be retrieved from the json input file!");
                            }
                        }
                        if (get_int_from_objectlist("ORDER",number_readobjects,&ORDER,varname_list, value_list)){
                            declare_error("Variable ORDER could not be retrieved from the json input file!");
                        }
                    }
                    if (TIME_FILT==2) {
                        if (get_float_from_objectlist("F_HIGH_PASS",number_readobjects,&F_HIGH_PASS,varname_list, value_list)){
                            F_HIGH_PASS=0.0;
                            fprintf(fp,"Variable F_HIGH_PASS is set to default value %f.\n",F_HIGH_PASS);}
                        if (get_string_from_objectlist("FREQ_FILE",number_readobjects,FREQ_FILE,varname_list, value_list))
                            declare_error("Variable FREQ_FILE could not be retrieved from the json input file!");
                        if (get_int_from_objectlist("ORDER",number_readobjects,&ORDER,varname_list, value_list))
                            declare_error("Variable ORDER could not be retrieved from the json input file!");
                    }
                }
                
                
                
                /* Gradient calculation */
                if (get_int_from_objectlist("LNORM",number_readobjects,&LNORM,varname_list, value_list))
                    declare_error("Variable LNORM could not be retrieved from the json input file!");
                if (get_int_from_objectlist("NORMALIZE",number_readobjects,&NORMALIZE,varname_list, value_list))
                    declare_error("Variable NORMALIZE could not be retrieved from the json input file!");
                if (get_int_from_objectlist("DTINV",number_readobjects,&DTINV,varname_list, value_list))
                    declare_error("Variable DTINV could not be retrieved from the json input file!");
                if (LNORM==8){
                    if (get_float_from_objectlist("WATERLEVEL_LNORM8",number_readobjects,&WATERLEVEL_LNORM8,varname_list, value_list))
                        declare_error("Variable WATERLEVEL_LNORM8 could not be retrieved from the json input file!");
                }
                
                
                /* Step length estimation */
                if (get_float_from_objectlist("EPS_SCALE",number_readobjects,&EPS_SCALE,varname_list, value_list))
                    declare_error("Variable EPS_SCALE could not be retrieved from the json input file!");
                if (get_int_from_objectlist("STEPMAX",number_readobjects,&STEPMAX,varname_list, value_list))
                    declare_error("Variable STEPMAX could not be retrieved from the json input file!");
                if (get_float_from_objectlist("SCALEFAC",number_readobjects,&SCALEFAC,varname_list, value_list))
                    declare_error("Variable SCALEFAC could not be retrieved from the json input file!");
                
                
                /* Termination of the program */
                if (get_float_from_objectlist("PRO",number_readobjects,&PRO,varname_list, value_list))
                    declare_error("Variable PRO could not be retrieved from the json input file!");
                if (get_int_from_objectlist("MIN_ITER",number_readobjects,&MIN_ITER,varname_list, value_list)){
                    MIN_ITER=0;
                    fprintf(fp,"Variable MIN_ITER is set to default value %d.\n",MIN_ITER);}
                
                
                /* Trace killing */
                if (get_int_from_objectlist("TRKILL",number_readobjects,&TRKILL,varname_list, value_list)){
                    TRKILL=0;
                    fprintf(fp,"Variable TRKILL is set to default value %d.\n",TRKILL);}
                else {
                    if (TRKILL==1) {
                        if (get_int_from_objectlist("TRKILL_OFFSET",number_readobjects,&TRKILL_OFFSET,varname_list, value_list)){
                            TRKILL_OFFSET=0;
                            if (get_string_from_objectlist("TRKILL_FILE",number_readobjects,TRKILL_FILE,varname_list, value_list))
                                declare_error("Variable TRKILL_FILE could not be retrieved from the json input file!");
                        } else {
                            if(TRKILL_OFFSET==0) { /* Only TRKILL File */
                                if (get_string_from_objectlist("TRKILL_FILE",number_readobjects,TRKILL_FILE,varname_list, value_list))
                                    declare_error("Variable TRKILL_FILE could not be retrieved from the json input file!");
                            }
                            if(TRKILL_OFFSET==1) { /* Only Offset based TRKill */
                                if (get_float_from_objectlist("TRKILL_OFFSET_LOWER",number_readobjects,&TRKILL_OFFSET_LOWER,varname_list, value_list)){
                                    TRKILL_OFFSET_LOWER=0.0;
                                }
                                if (get_float_from_objectlist("TRKILL_OFFSET_UPPER",number_readobjects,&TRKILL_OFFSET_UPPER,varname_list, value_list)){
                                    declare_error("Variable TRKILL_OFFSET_UPPER could not be retrieved from the json input file!");
                                }
                            }
                            if(TRKILL_OFFSET==2){ /* Both Offset based TRKill & File */
                                if (get_float_from_objectlist("TRKILL_OFFSET_LOWER",number_readobjects,&TRKILL_OFFSET_LOWER,varname_list, value_list)){
                                    TRKILL_OFFSET_LOWER=0.0;
                                }
                                if (get_float_from_objectlist("TRKILL_OFFSET_UPPER",number_readobjects,&TRKILL_OFFSET_UPPER,varname_list, value_list)){
                                    declare_error("Variable TRKILL_OFFSET_UPPER could not be retrieved from the json input file!");
                                }
                                if (get_string_from_objectlist("TRKILL_FILE",number_readobjects,TRKILL_FILE,varname_list, value_list))
                                    declare_error("Variable TRKILL_FILE could not be retrieved from the json input file!");
                            }
                        }
                    }
                }
                
                
                /* Time windowing and VPPML */
                if (get_int_from_objectlist("TIMEWIN",number_readobjects,&TIMEWIN,varname_list, value_list)){
                    TIMEWIN=0;
                    fprintf(fp,"Variable TIMEWIN is set to default value %d.\n",TIMEWIN);}
                else {
                    if (TIMEWIN==1){
                        if (get_int_from_objectlist("TW_IND",number_readobjects,&TW_IND,varname_list, value_list)){
                            TW_IND=0;
                            fprintf(fp,"Variable TW_IND is set to default value %d.\n",TW_IND);
                        } else {
                            if (TW_IND>2) declare_error("Only TW_IND=1 (one time window) and TW_IND=2 (two time windows) possible");
                        }
                        if (get_string_from_objectlist("PICKS_FILE",number_readobjects,PICKS_FILE,varname_list, value_list))
                            declare_error("Variable PICKS_FILE could not be retrieved from the json input file!");
                        if (get_float_from_objectlist("TWLENGTH_PLUS",number_readobjects,&TWLENGTH_PLUS,varname_list, value_list))
                            declare_error("Variable TWLENGTH_PLUS could not be retrieved from the json input file!");
                        if (get_float_from_objectlist("TWLENGTH_MINUS",number_readobjects,&TWLENGTH_MINUS,varname_list, value_list))
                            declare_error("Variable TWLENGTH_MINUS could not be retrieved from the json input file!");
                        if (get_float_from_objectlist("GAMMA",number_readobjects,&GAMMA,varname_list, value_list))
                            declare_error("Variable GAMMA could not be retrieved from the json input file!");
                    }
                }
            } /* end if (FORWARD_ONLY==0) */
            else {
                if (FORWARD_ONLY>0){	/* only forward modeling is applied */
                    FORWARD_ONLY=1;
                    ITERMAX=1;
		    INV_STF=0;
                    strcpy(INV_MODELFILE,MFILE);
                    DTINV=1;
                    INVTYPE=2;
                    fprintf(fp,"Setting some default values needed for modeling:\n");
                    fprintf(fp,"Variable ITERMAX is set to default value %d.\n",ITERMAX);
                    fprintf(fp,"Variable INV_MODELFILE is set to default value %s.\n",MFILE);
                    fprintf(fp,"Variable DTINV is set to default value %d.\n",DTINV);
                    fprintf(fp,"Variable INVTYPE is set to default value %d.\n",INVTYPE);
                }
            }
            
        }
        
        fprintf(fp,"\nEnd of setting default values\n");
        fprintf(fp,"=====================================\n\n");
        
        
        /********************************************/
        /* Check files and directories if necessary */
        /********************************************/
        
        /* signal file */
        if (SOURCE_SHAPE == 3)
        {
            if (access(SIGNAL_FILE,0) != 0)
            {
                fprintf(fp, "\n==================================================================\n");
                fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
                fprintf(fp, "        The signal file does not exist!\n");
                fprintf(fp, "        File name: <%s>", SIGNAL_FILE);
                fprintf(fp, "\n==================================================================\n");
                fserr = 1;
            }
            else if (access(SIGNAL_FILE,4) != 0)
            {
                fprintf(fp, "\n==================================================================\n");
                fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
                fprintf(fp, "        The signal file does not have read access!\n");
                fprintf(fp, "        File name: <%s>", SIGNAL_FILE);
                fprintf(fp, "\n==================================================================\n");
                fserr = 1;
            }
        }
        
        /* source file */
        if (SRCREC==1)
        {
            if (access(SOURCE_FILE,0) != 0)
            {
                fprintf(fp, "\n==================================================================\n");
                fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
                fprintf(fp, "        The source file does not exist!\n");
                fprintf(fp, "        File name: <%s>", SOURCE_FILE);
                fprintf(fp, "\n==================================================================\n");
                fserr = 1;
            }
            else if (access(SOURCE_FILE,4) != 0)
            {
                fprintf(fp, "\n==================================================================\n");
                fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
                fprintf(fp, "        The source file does not have read access!\n");
                fprintf(fp, "        File name: <%s>", SOURCE_FILE);
                fprintf(fp, "\n==================================================================\n");
                fserr = 1;
            }
        }
        
        /* model file -> this is checked in more detail in readmod.c!
         if (READMOD)
         {
         if (access(MFILE,0) != 0)
         {
         fprintf(fp, "\n==================================================================\n");
         fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
         fprintf(fp, "        The model file does not exist!\n");
         fprintf(fp, "        File name: <%s>", MFILE);
         fprintf(fp, "\n==================================================================\n");
         fserr = 1;
         }
         else if (access(MFILE,4) != 0)
         {
         fprintf(fp, "\n==================================================================\n");
         fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
         fprintf(fp, "        The model file does not have read access!\n");
         fprintf(fp, "        File name: <%s>", MFILE);
         fprintf(fp, "\n==================================================================\n");
         fserr = 1;
         }
         }*/
        
        /* receiver file */
        if (READREC==1)
        {
            if (access(REC_FILE,0) != 0)
            {
                fprintf(fp, "\n==================================================================\n");
                fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
                fprintf(fp, "        The receiver file does not exist!\n");
                fprintf(fp, "        File name: <%s>", REC_FILE);
                fprintf(fp, "\n==================================================================\n");
                fserr = 1;
            }
            else if (access(REC_FILE,4) != 0)
            {
                fprintf(fp, "\n==================================================================\n");
                fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
                fprintf(fp, "        The receiver file does not have read access!\n");
                fprintf(fp, "        File name: <%s>", REC_FILE);
                fprintf(fp, "\n==================================================================\n");
                fserr = 1;
            }
        }
        
        
//         /* trace kill file */
//         if (TRKILL == 1)
//         {
//             if (access(TRKILL_FILE,0) != 0)
//             {
//                 fprintf(fp, "\n==================================================================\n");
//                 fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
//                 fprintf(fp, "        The trace kill file does not exist!\n");
//                 fprintf(fp, "        File name: <%s>", TRKILL_FILE);
//                 fprintf(fp, "\n==================================================================\n");
//                 fserr = 1;
//             }
//             else if (access(TRKILL_FILE,4) != 0)
//             {
//                 fprintf(fp, "\n==================================================================\n");
//                 fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
//                 fprintf(fp, "        The trace kill file does not have read access!\n");
//                 fprintf(fp, "        File name: <%s>", TRKILL_FILE);
//                 fprintf(fp, "\n==================================================================\n");
//                 fserr = 1;
//             }
//         }
//         
//         
//         /* trace kill STF file */
//         if (TRKILL_STF == 1)
//         {
//             if (access(TRKILL_FILE_STF,0) != 0)
//             {
//                 fprintf(fp, "\n==================================================================\n");
//                 fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
//                 fprintf(fp, "        The trace kill file does not exist!\n");
//                 fprintf(fp, "        File name: <%s>", TRKILL_FILE_STF);
//                 fprintf(fp, "\n==================================================================\n");
//                 fserr = 1;
//             }
//             else if (access(TRKILL_FILE_STF,4) != 0)
//             {
//                 fprintf(fp, "\n==================================================================\n");
//                 fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
//                 fprintf(fp, "        The trace kill file does not have read access!\n");
//                 fprintf(fp, "        File name: <%s>", TRKILL_FILE_STF);
//                 fprintf(fp, "\n==================================================================\n");
//                 fserr = 1;
//             }
//         }
        
        
        /* frequency file */
        if (TIME_FILT==2)
        {
            if (access(FREQ_FILE,0) != 0)
            {
                fprintf(fp, "\n==================================================================\n");
                fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
                fprintf(fp, "        The frequency file does not exist!\n");
                fprintf(fp, "        File name: <%s>", FREQ_FILE);
                fprintf(fp, "\n==================================================================\n");
                fserr = 1;
            }
            else if (access(FREQ_FILE,4) != 0)
            {
                fprintf(fp, "\n==================================================================\n");
                fprintf(fp, "  ERROR parsing input file <%s>:\n", fileinp);
                fprintf(fp, "        The frequency file does not have read access!\n");
                
                fprintf(fp, "\n==================================================================\n");
                fserr = 1;
            }
        }
        
        /********************************************/
        /* ERROR                                    */
        /********************************************/
        if (fserr)
        {
            fprintf(fp, "\n");
            sprintf(errormessage, "\n  in: <read_par_json.c> \n");
            declare_error(errormessage);
        }
        
        
    } /* End of if(MYID==0) */
}
