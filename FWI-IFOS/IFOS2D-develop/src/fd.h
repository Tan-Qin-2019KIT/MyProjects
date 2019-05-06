/*------------------------------------------------------------------------
 *  fd.h - include file for viscoelastic FD programs          
 *  See COPYING file for copying and redistribution conditions.
 *  ---------------------------------------------------------------------*/

/* files to include */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#define iround(x) ((int)(floor)(x+0.5))
#define min(x,y) ((x<y)?x:y)    
#define max(x,y) ((x<y)?y:x)
#define fsign(x) ((x<0.0)?(-1):1)    

#define PI (3.141592653589793)
#define NPAR 120
#define STRING_SIZE 74
#define STRING_SIZE2 256
#define REQUEST_COUNT 4
#define WORKFLOW_MAX_VAR 14

/* declaration of functions */



void window_cos(float **win, int npad, int nsrc, float it1, float it2, float it3, float it4);

void catseis(float **data, float **fulldata, int *recswitch, int ntr_glob, MPI_Comm newcomm_nodentr);

void stf(FILE *fp, float **sectionvy, float ** sectionvy_obs, float ** sectionvy_conv, float * source_time_function, int  **recpos, int  **recpos_loc,
         int ntr_glob,int ntr, float ** srcpos, int ishot, int ns, int iter, int nshots, float F_LOW_PASS, int SH,int nsrc_glob);

int **splitrec(int **recpos,int *ntr_loc, int ntr, int *recswitch);

void absorb(float ** absorb_coeff);

void taper_grad(float ** waveconv, float ** taper_coeff, float **srcpos, int nshots, int **recpos, int ntr, int sws);

void taper_grad_shot(float ** waveconv,float ** taper_coeff, float **srcpos, int nshots, int **recpos, int ntr, int ishot, int sws);

void spat_filt(float ** waveconv, int iter, int sws);

float norm(float ** waveconv, int iter, int sws);

void alloc_sections(int ntr,int ns,float ***sectionvx,float ***sectionvy,float ***sectionvz,float ***sectionp,float ***sectionpnp1,float ***sectionpn,float ***sectioncurl,float ***sectiondiv,
                    float ***sectionpdata,float ***sectionpdiff,float ***sectionpdiffold,float ***sectionvxdata,float ***sectionvxdiff,float ***sectionvxdiffold,float ***sectionvydata,
                    float ***sectionvydiff,float ***sectionvydiffold,float ***sectionvzdata,float ***sectionvzdiff,float ***sectionvzdiffold);

void av_mat(float **  pi, float **  u,
            float **  ppijm, float **  puip, float ** pujm);

void av_mue(float ** u, float ** uipjp,float ** rho);

void av_rho(float **rho, float **rip, float **rjp);

void av_tau(float **taus, float **tausipjp);

float median2d(float **mat, int ny, int nx);

void calc_mat_change(float  **  waveconv, float ** waveconv_rho, float ** waveconv_u, float  **  rho, float  **  rhonp1, float **  pi, float **  pinp1, float **  u,
                     float **  unp1, int iter, int epstest, int calcneweps, float eps_scale_vp, float eps_scale_vs);

void calc_mat_change_test(float  **  waveconv, float  **  waveconv_rho, float  **  waveconv_u, float  **  rho, float  **  rhonp1, float **  pi, float **  pinp1, float **  u, float **  unp1, int iter,
                          int epstest, int FORWARD_ONLY, float eps_scale, int itest, int nfstart, float ** u_start, float ** pi_start, float ** rho_start,int wavetype_start,float **bfgsmod,int bfgsnum,int bfgspar,float Vs_avg,float Vp_avg,float rho_avg,int LBFGS_iter_start);

double calc_res(float **sectiondata, float **section, float **sectiondiff, float **sectiondiffold, int ntr, int ns, int LNORM, float L2, int itest, int sws, int swstestshot, int ntr_glob, int **recpos_loc, int nsrc_glob, int ishot, int iter, float ** srcpos, int ** recpos);

double calc_misfit(float **sectiondata, float **section, int ntr, int ns, int LNORM, float L2, int itest, int sws, int swstestshot,int ntr_glob, int **recpos_loc, int nsrc_glob, int ishot, int iter, float ** srcpos, int ** recpos);

float calc_opt_step(float *  L2t, float ** waveconv, float ** gradp, float * epst, int sws, float C_vp);

float calc_opt_step_test(float *  L2t, float ** waveconv, float ** gradp, float * epst, int sws, float C_vp);

double calc_energy(float **sectiondata, int ntr, int ns, float energy, int ntr_glob, int **recpos_loc, int nsrc_glob, int ishot, int iter, float ** srcpos, int ** recpos);

void checkfd(FILE *fp, float ** prho, float ** ppi, float ** pu, float ** ptaus, float ** ptaup, float * peta, float *hc, float **srcpos, int nsrc, int **recpos, int ntr);

void checkfd_hc(FILE *fp, float ** prho, float ** ppi, float ** pu,
                float ** ptaus, float ** ptaup, float *peta, float *hc);

void checkfd_ssg_elastic(FILE *fp, float ** prho, float ** ppi, float ** pu, float *hc);
void checkfd_ssg_visc(FILE *fp, float ** prho, float ** ppi, float ** pu, float ** ptaus, float ** ptaup, float * peta, float *hc);

void checkfd_rsg(FILE *fp, float ** prho, float ** ppi, float ** pu,
                 float ** ptaus, float ** ptaup, float *peta);

void comm_ini(float ** bufferlef_to_rig, float ** bufferrig_to_lef,
              float ** buffertop_to_bot, float ** bufferbot_to_top,
              MPI_Request *req_send, MPI_Request *req_rec);

void conv_FD(float * temp_TS, float * temp_TS1, float * temp_conv, int ns);

void count_killed_traces(int ntr, int swstestshot, int ntr_glob, int **recpos_loc, int nsrc_glob, int ishot, int* ptr_killed_traces, int* ptr_killed_traces_testshots,float ** srcpos, int ** recpos);

void create_trkill_table(int ** killtable, int ntr_glob, int **recpos, int nsrc_glob, float **srcpos, int ishot, float kill_offset_lower, float kill_offset_upper);

void dealloc_sections(int ntr,int ns,int **recpos_loc,float **sectionvx,float **sectionvy,float **sectionvz,float **sectionp,float **sectionpnp1,float **sectionpn,float **sectioncurl,float **sectiondiv,
                    float **sectionpdata,float **sectionpdiff,float **sectionpdiffold,float **sectionvxdata,float **sectionvxdiff,float **sectionvxdiffold,float **sectionvydata,
                    float **sectionvydiff,float **sectionvydiffold,float **sectionvzdata,float **sectionvzdiff,float **sectionvzdiffold);

float exchange_L2(float L2, int sw, int bcast_l2);

void exchange_rsg(float ** vx, float ** vy, float ** vz,
                  float ** bufferlef_to_rig, float ** bufferrig_to_lef,
                  float ** buffertop_to_bot, float ** bufferbot_to_top);

void exchange_rsg_4th(float ** vx, float ** vy, float ** vz,
                      float ** bufferlef_to_rig, float ** bufferrig_to_lef,
                      float ** buffertop_to_bot, float ** bufferbot_to_top);

void exchange_v(float ** vx, float ** vy, float ** vz,
                float ** bufferlef_to_rig, float ** bufferrig_to_lef,
                float ** buffertop_to_bot, float ** bufferbot_to_top,
                MPI_Request * req_send, MPI_Request * req_rec, int wavetyp_start);

void exchange_s(float ** sxx, float ** syy,
                float ** sxy,float ** sxz,float ** syz, float ** bufferlef_to_rig, float ** bufferrig_to_lef,
                float ** buffertop_to_bot, float ** bufferbot_to_top,
                MPI_Request * req_send, MPI_Request * req_rec, int wavetyp_start);

void exchange_par(void);

void exchange_mod_es(float ** matmod, int ncptot, int nparameter);

void  FFT_filt(float ** data, float freqshift, int ntr, int ns,int method);

/*short FFT(short int dir,long m,float *x,float *y);*/

void  FFT(int isign, unsigned long nlog, float *re, float *im); /* NR version*/

float *filter_frequencies(int *nfrq);

float *holbergcoeff(void);

void info(FILE *fp);

void initproc(void);

void interpol(int ni1, int ni2, float **  intvar, int cfgt_check);

double LU_decomp(double  **A, double *x, double *b,int n);

float minimum_m(float **mat, int nx, int ny);
float maximum_m(float **mat, int nx, int ny);


void model(float  **  rho, float **  pi, float **  u,
           float **  taus, float **  taup, float *  eta);

void model_elastic(float  **  rho, float **  pi, float **  u);

void model_ani(float  **  rho, float **  c11, float **  c15, float **  c13,
               float **  c35, float **  c33, float **  c55,
               float **  taus, float **  taup, float *  eta);

void matcopy(float ** prho, float ** ppi, float ** pu, float ** ptaup,
             float ** ptaus);

void matcopy_elastic(float ** prho, float ** ppi, float ** pu);

void matcopy_ani(float ** rho, float **  c11, float **  c15, float **  c13,
                 float **  c35, float **  c33, float **  c55, float ** taus,
                 float ** taup);

void max_grad(float  **  waveconv, float  **  waveconv_rho, float  **  waveconv_u, float  **  rho, float **  pi, float **  u);

void merge(int nsnap, int type);

void merge2(int nsnap, int type);

void mergemod(char modfile[STRING_SIZE], int format);

void note(FILE *fp);


void  outseis_glob(FILE *fp, FILE *fpdata, int comp, float **section,
                   int **recpos, int **recpos_loc, int ntr, float ** srcpos_loc,
                   int nsrc, int ns, int seis_form, int ishot, int sws);

void  outseis_vector(FILE *fp, FILE *fpdata, int comp, float *section,
                     int **recpos, int **recpos_loc, int ntr, float ** srcpos_loc,
                     int nsrc, int ns, int seis_form, int ishot, int sws);

void  inseis(FILE *fp, int comp, float **section, int ntr, int ns, int sws, int iter);

void  inseis_source_wavelet(float *section, int ns, int ishot, int SH, int STF);

void  taper(float *section, int ns, float fc);

void  output_source_signal(FILE *fp, float **signals, int ns, int seis_form);

void PCG(float ** waveconv, float ** taper_coeff, int nsrc, float ** srcpos, int ** recpos, int ntr_glob, int iter, float C_vp, float ** gradp, int nfstart_jac, float ** waveconv_u, float C_vs, float ** gradp_u, float ** waveconv_rho, float C_rho, float ** gradp_rho, float Vs_avg, float F_LOW_PASS, int PCG_iter_start);

void PCG_SH(float ** taper_coeff, int nsrc, float ** srcpos, int ** recpos, int ntr_glob, int iter, int nfstart_jac, float ** waveconv_u, float C_vs, float ** gradp_u, float ** waveconv_rho, float C_rho, float ** gradp_rho, float Vs_avg, float F_LOW_PASS, int PCG_iter_start);

void PML_pro(float * d_x, float * K_x, float * alpha_prime_x, float * a_x, float * b_x,
             float * d_x_half, float * K_x_half, float * alpha_prime_x_half, float * a_x_half, float * b_x_half,
             float * d_y, float * K_y, float * alpha_prime_y, float * a_y, float * b_y,
             float * d_y_half, float * K_y_half, float * alpha_prime_y_half, float * a_y_half, float * b_y_half);

void psource(int nt, float ** sxx, float ** syy, float ** sp,
             float **  srcpos_loc, float ** signals, int nsrc, int sw);

void psource_rsg(int nt, float ** sxx, float ** syy,
                 float **  srcpos_loc, float ** signals, int nsrc);


float *rd_sour(int *nts,FILE* fp_source);

float readdsk(FILE *fp_in, int format);

void readbufs(float ** sxx, float ** syy,
              float ** sxy, float ** bufferlef_to_rig, float ** bufferrig_to_lef,
              float ** buffertop_to_bot, float ** bufferbot_to_top);

void readbufv(float ** vx, float ** vy,
              float ** bufferlef_to_rig, float ** bufferrig_to_lef,
              float ** buffertop_to_bot, float ** bufferbot_to_top);

void read_checkpoint(int nx1, int nx2, int ny1, int ny2,
                     float **  vx, float ** vy, float ** sxx, float ** syy, float ** sxy);

void read_par_json(FILE *fp, char *fileinp);

void readmod(float  **  rho, float **  pi, float **  u,
             float **  taus, float **  taup, float *  eta);

void readmod_elastic(float  **  rho, float **  pi, float **  u);

void readmod_elastic_es(float  **  rho, float **  pi, float **  u, float ** matmod, int is);

int **receiver(int* ntr, float** srcpos, int shotno);

void save_checkpoint(int nx1, int nx2, int ny1, int ny2,
                     float **  vx, float ** vy, float ** sxx, float ** syy, float ** sxy);

void saveseis_glob(FILE *fp, float **sectionvx, float **sectionvy,float **sectionvz,float **sectionp,float **sectioncurl, float **sectiondiv, int  **recpos, int  **recpos_loc,int ntr, float ** srcpos, int ishot, int ns, int iter, int type_switch);

void snap(FILE *fp,int nt, int nsnap, float **vx, float **vy, float **sxx,
          float **syy, float **sp, float **u, float **pi, float *hc, int ishot);


void snap_rsg(FILE *fp,int nt, int nsnap, float **vx, float **vy, float **sxx, float **syy, float **u, float **pi);

void snapmerge(int nsnap);

void snapmerge2(int nsnap);

float **sources(int *nsrc);

void solvelin(float  **AA, float *bb, float *x, int e, int method);

void seismo(int lsamp, int ntr, int **recpos, float **sectionvx,
            float **sectionvy, float **sectionp, float **sectioncurl, float **sectiondiv,
            float **pvx, float **pvy, float **psxx, float **psyy, float **ppi, float **pu);

void seismo_ssg(int lsamp, int ntr, int **recpos, float **sectionvx,
                float **sectionvy, float **sectionvz,float **sectionp, float **sectioncurl, float **sectiondiv,
                float **pvx, float **pvy,float **pvz, float **psxx, float **psyy, float **psp, float **ppi, float **pu,
                float *hc);

void seismo_rsg(int lsamp, int ntr, int **recpos, float **sectionvx,
                float **sectionvy, float **sectionp, float **sectioncurl, float **sectiondiv,
                float **pvx, float **pvy, float **psxx, float **psyy, float **ppi, float **pu);

float **splitsrc(float **srcpos,int *nsrc_loc, int nsrc);

float **splitsrc_back(int **recpos,int *nsrc_loc, int nsrc);

void stalta(float **sectionp, int ntr, int nst, float *picked_times);

void surface(int ndepth, float ** pvx, float ** pvy,
             float ** psxx, float ** psyy,
             float ** psxy, float *** pp, float *** pq,
             float  **  ppi, float  **  pu, float ** ptaup,
             float ** ptaus, float * etajm, float * peta);

void surface_elastic(int ndepth, float ** vx, float ** vy, float ** sxx, float ** syy,
                     float ** sxy, float  **  pi, float  **  u, float ** rho, float * hc);

void surface_elastic_PML(int ndepth, float ** vx, float ** vy, float ** sxx, float ** syy, float ** sxy, float ** syz, float  **  pi, float  **  u, float ** rho, float * hc, float * K_x, float * a_x, float * b_x, float ** psi_vxx, float ** ux, float ** uy, float ** uxy, float ** uyz,float ** sxz,float **uxz);

void surface_PML(int ndepth, float ** vx, float ** vy, float ** sxx, float ** syy, float ** sxy, float ** syz, float ***p, float ***q, float  **  ppi, float  **  pu, float **prho, float **ptaup, float **ptaus, float *etajm, float *peta, float * hc, float * K_x, float * a_x, float * b_x, float ** psi_vxx, float ** ux, float ** uy, float ** uxy, float ** uyz,float ** sxz,float **uxz);

void  timedomain_filt(float ** data, float fc, int order, int ntr, int ns, int method);
void  timedomain_filt_vector(float * data, float fc, int order, int ns, int method);

void time_window(float **sectiondata, int iter, int ntr_glob, int **recpos_loc, int ntr, int ns, int ishot);
void time_window_glob(float **sectiondata, int iter, int ntr_glob, int ns, int ishot);

void prepare_update_s(float *etajm, float *etaip, float *peta, float **fipjp, float **pu,
                      float **puipjp, float **ppi, float **prho, float **ptaus, float **ptaup,
                      float **ptausipjp, float **f, float **g, float *bip, float *bjm,
                      float *cip, float *cjm, float ***dip, float ***d, float ***e);

void update_s(int nx1, int nx2, int ny1, int ny2,
              float **  vx, float **   vy, float **   sxx, float **   syy,
              float **   sxy, float *** r, float *** p, float *** q,
              float ** ppi, float ** pu, float ** taup, float ** taus,
              float *   etaip, float *   etajm, float * peta);

void update_s_visc_hc(int nx1, int nx2, int ny1, int ny2,
                      float **  vx, float **   vy, float **   sxx, float **   syy,
                      float **   sxy, float *** r, float *** p, float *** q,
                      float ** ppi, float ** pu, float **uipjp, float ** taup, float ** taus,
                      float **tausipjp, float *   etaip, float *   etajm, float * peta, float *hc);

void update_s_elastic_PML_SH(int nx1, int nx2, int ny1, int ny2, float **  vz, float **   sxz, float **   syz, float ** uxz, float ** uyz, float *hc,  int infoout,float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
                             float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,float ** psi_vzx, float ** psi_vzy,float ** uipjp,float ** u,float ** rho);

void update_s_visc_PML_SH(int nx1, int nx2, int ny1, int ny2, float **  vz, float **   sxz, float **   syz, float ***t, float ***o, float *bip, float *bjm, float *cip, float *cjm, float ***d, float ***dip, float **fipjp, float **f, float *hc,  int infoout,float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,float ** psi_vzx, float ** psi_vzy);

void update_s_rsg(int nx1, int nx2, int ny1, int ny2,
                  float ** pvx, float ** pvy, float ** psxx, float ** psyy,
                  float ** psxy, float *** pr, float *** pp, float ***pq,
                  float  **  ppi, float  **  pu, float ** ptaup,
                  float ** ptaus, float * etaip,
                  float * etajm, float * peta, float ** absorb_coeff);

void update_s_rsg_4th(int nx1, int nx2, int ny1, int ny2,
                      float ** pvx, float ** pvy, float ** psxx, float ** psyy,
                      float ** psxy, float *** pr, float *** pp, float ***pq,
                      float  **  ppi, float  **  pu, float ** ptaup,
                      float ** ptaus, float * etaip,
                      float * etajm, float * peta);

void update_s_ani(int nx1, int nx2, int ny1, int ny2,
                  float **  vx, float **   vy, float **   sxx, float **   syy,
                  float **   sxy, float *** r, float *** p, float *** q,
                  float  ** c11, float  **  c15, float  ** c13, float  **  c35,
                  float  ** c33, float  **  c55, float **   ptaup, float **   ptaus,
                  float *   etaip, float *   etajm, float * peta);

void update_s_elastic(int nx1, int nx2, int ny1, int ny2,
                      float **  vx, float **   vy, float **   sxx, float **   syy,
                      float **   sxy, float ** pi, float ** u, float ** uipjm, float ** absorb_coeff);

void update_s_elastic_rsg(int nx1, int nx2, int ny1, int ny2,
                          float **  vx, float **  vy, float **  sxx, float **  syy,
                          float **  sxy, float  **   pi, float  **   u, float ** absorb_coeff);

void update_s_elastic_hc(int nx1, int nx2, int ny1, int ny2,
                         float **  vx, float **   vy, float **  ux, float **   uy, float **  uxy, float **   uyx, float **   sxx, float **   syy,
                         float **   sxy, float ** pi, float ** u, float ** uipjp, float ** absorb_coeff, float ** rho,
                         float *hc, int infoout);

void update_s_elastic_PML(int nx1, int nx2, int ny1, int ny2,
                          float **  vx, float **   vy, float **  ux, float **   uy, float **  uxy, float **   uyx, float **   sxx, float **   syy,
                          float **   sxy, float ** pi, float ** u, float ** uipjp, float ** absorb_coeff, float **rho, float *hc, int infoout,
                          float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
                          float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,
                          float ** psi_vxx, float ** psi_vyy, float ** psi_vxy, float ** psi_vyx);

void update_s_elastic_hh(int nx1, int nx2, int ny1, int ny2,
                         float **  vx, float **   vy, float **   sxx, float **   syy,
                         float **   sxy, float ** pi, float ** u );

void update_s_visc_PML(int nx1, int nx2, int ny1, int ny2,
                       float **  vx, float **   vy, float **  ux, float **   uy, float **  uxy, float **   uyx, float **   sxx, float **   syy,
                       float **   sxy, float ** pi, float ** u, float ** uipjp, float **rho, float *hc, int infoout,
                       float ***r, float ***p, float ***q, float **fipjp, float **f, float **g, float *bip, float *bjm, float *cip, float *cjm, float ***d, float ***e, float ***dip,
                       float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
                       float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,
                       float ** psi_vxx, float ** psi_vyy, float ** psi_vxy, float ** psi_vyx);


void update_v(int nx1, int nx2, int ny1, int ny2, int nt,
              float **  pvx, float ** pvy, float ** psxx, float ** psyy,
              float ** psxy, float  ** prho,
              float **  srcpos_loc, float ** signals, int nsrc, float ** absorb_coeff);


void update_v_hc(int nx1, int nx2, int ny1, int ny2, int nt,
                 float **  vx, float **  vxp1, float **  vxm1, float ** vy, float **  vyp1, float **  vym1, float **  uttx, float **  utty, float ** sxx, float ** syy,
                 float ** sxy, float  **rip, float **rjp,
                 float **  srcpos_loc, float ** signals, float ** signals1, int nsrc, float ** absorb_coeff,
                 float *hc, int infoout, int sw);

void update_v_PML(int nx1, int nx2, int ny1, int ny2, int nt,
                  float **  vx, float **  vxp1, float **  vxm1, float ** vy, float **  vyp1, float **  vym1, float **  uttx, float **  utty,float ** sxx, float ** syy,
                  float ** sxy, float  **rip, float **rjp, float **  srcpos_loc, float ** signals, float ** signals1, int nsrc, float ** absorb_coeff,
                  float *hc, int infoout,int sw, float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
                  float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,
                  float ** psi_sxx_x, float ** psi_syy_y, float ** psi_sxy_y, float ** psi_syx_x);

void update_v_PML_SH(int nx1, int nx2, int ny1, int ny2, int nt,
                     float **  vz, float **  vzp1, float **  vzm1, float ** sxz, float ** syz,float **rho, float **  srcpos_loc, float ** signals, float ** signals1, int nsrc, float ** absorb_coeff,
                     float *hc, int infoout,int sw, float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
                     float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,float ** psi_sxz_x,float ** psi_syz_y);

void update_v_hh(int nx1, int nx2, int ny1, int ny2, int nt,
                 float **  pvx, float ** pvy, float ** psxx, float ** psyy,
                 float ** psxy, float  ** prho,
                 float **  srcpos_loc, float ** signals, int nsrc, float ** absorb_coeff);

void update_v_rsg(int nx1, int nx2, int ny1, int ny2, int nt,
                  float **  pvx, float ** pvy, float ** psxx, float ** psyy,
                  float ** psxy, float  ** prho,
                  float **  srcpos_loc, float ** signals, int nsrc, float ** absorb_coeff);

void update_v_rsg_4th(int nx1, int nx2, int ny1, int ny2, int nt,
                      float **  pvx, float ** pvy, float ** psxx, float ** psyy,
                      float ** psxy, float  ** prho,
                      float **  srcpos_loc, float ** signals, int nsrc, float ** absorb_coeff);

float ** wavelet(float ** srcpos_loc, int nsrc, int ishot, int SH, int STF);
float ** wavelet_stf(int nsrc, int ishot, float ** signals_stf);

void writebufs(float ** sxx, float ** syy,
               float ** sxy, float ** bufferlef_to_rig, float ** bufferrig_to_lef,
               float ** buffertop_to_bot, float ** bufferbot_to_top);

void writebufv(float ** vx, float ** vy,
               float ** bufferlef_to_rig, float ** bufferrig_to_lef,
               float ** buffertop_to_bot, float ** bufferbot_to_top);

void write_par(FILE *fp);

void writedsk(FILE *fp_out, float amp, int format);

void writemod(char modfile[STRING_SIZE], float ** array, int format);

void zero_fdveps(int ny1, int ny2, int nx1, int nx2, float ** vx, float ** vy, float ** vz, float ** sxx, float ** syy, float ** sxy,float ** sxz,float ** syz,float ** vxm1, float ** vym1, float ** uxy, float ** vxp1, float ** vyp1,float ** psi_sxx_x, float ** psi_sxy_x, float ** psi_sxz_x, float ** psi_vxx, float ** psi_vyx, float ** psi_vzx, float ** psi_syy_y, float ** psi_sxy_y, float ** psi_syz_y, float ** psi_vyy, float ** psi_vxy,float ** psi_vzy,float ** psi_vxxs);

void zero_fdveps_visc(int ny1, int ny2, int nx1, int nx2, float ** vx, float ** vy,float ** vz, float ** sxx,
                      float ** syy, float ** sxy,float ** sxz,float ** syz, float ** vxm1, float ** vym1, float ** uxy, float ** vxp1, float ** vyp1,
                      float ** psi_sxx_x, float ** psi_sxy_x, float ** psi_sxz_x, float ** psi_vxx, float ** psi_vyx, float ** psi_vzx, float ** psi_syy_y, float ** psi_sxy_y, float ** psi_syz_y, float ** psi_vyy, float ** psi_vxy,float ** psi_vzy,
                      float ** psi_vxxs, float ***pr, float ***pp, float ***pq, float ***pt, float ***po);

void FLnode(float  **  rho, float **  pi, float **  u, float **  taus, float **  taup, float *  eta);

void smooth(float ** mat, int sws, int filter, float Vs_avg, float F_LOW_PASS);

/* declaration of functions for parser*/

/* declaration of functions for json parser in json_parser.c*/
int read_objects_from_intputfile(FILE *fp, char input_file[STRING_SIZE],char ** varname_list,char ** value_list);

void print_objectlist_screen(FILE *fp, int number_readobject,char ** varname_list,char ** value_list);

int count_occure_charinstring(char stringline[STRING_SIZE], char teststring[]);

void copy_str2str_uptochar(char string_in[STRING_SIZE], char string_out[STRING_SIZE], char teststring[]);

int get_int_from_objectlist(char string_in[STRING_SIZE], int number_readobject, int * int_buffer,
                            char ** varname_list,char ** value_list);

int get_float_from_objectlist(char string_in[STRING_SIZE], int number_readobject, float * double_buffer,
                              char ** varname_list,char ** value_list);

int get_string_from_objectlist(char string_in[STRING_SIZE], int number_readobject, char string_buffer[STRING_SIZE],
                               char ** varname_list,char ** value_list);

int is_string_blankspace(char string_in[STRING_SIZE]);

void remove_blankspaces_around_string(char string_in[STRING_SIZE] );

void add_object_tolist(char string_name[STRING_SIZE],char string_value[STRING_SIZE], int * number_read_object,
                       char ** varname_list,char ** value_list );

void calc_hilbert(float ** datatrace, float ** envelope, int ns, int ntr);

void calc_envelope(float ** datatrace, float ** envelope, int ns, int ntr);

/* utility functions */
void declare_error(char err_text[]);
void warning(char warn_text[]);
double maximum(float **a, int nx, int ny);
float *vector(int nl, int nh);
int *ivector(int nl, int nh);
double *dvector(int nl, int nh);
float **fmatrix(int nrl, int nrh, int ncl, int nch);
int *ivector(int nl, int nh);
void quicksort(float *arr, int dummy, int elements);
float **matrix(int nrl, int nrh, int ncl, int nch);
int **imatrix(int nrl, int nrh, int ncl, int nch);
float ***f3tensor(int nrl, int nrh, int ncl, int nch,int ndl, int ndh);
void free_vector(float *v, int nl, int nh);
void free_dvector(double *v, int nl, int nh);
void free_ivector(int *v, int nl, int nh);
void free_matrix(float **m, int nrl, int nrh, int ncl, int nch);
void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch);
void free_f3tensor(float ***t, int nrl, int nrh, int ncl, int nch, int ndl,
                   int ndh);
void zero(float *A, int u_max);
void normalize_data(float **data, int ntr, int ns);

/* functions for acoustic modelling */

void model_acoustic(float  **  rho, float **  pi);
void readmod_acoustic(float  **  rho, float **  pi);
void matcopy_acoustic(float ** prho, float ** ppi);
void zero_fdveps_ac(int ny1, int ny2, int nx1, int nx2, float ** vx, float ** vy, float ** psp, float ** vxp1, float ** vyp1,
                    float ** psi_sxx_x, float ** psi_sxy_x, float ** psi_vxx, float ** psi_vyx,
                    float ** psi_syy_y, float ** psi_sxy_y, float ** psi_vyy, float ** psi_vxy, float ** psi_vxxs);

void update_v_acoustic_PML(int nx1, int nx2, int ny1, int ny2, int nt,
                           float **  vx, float **  vxp1, float **  vxm1, float ** vy, float **  vyp1, float **  vym1,
                           float ** sp, float  **rip, float **rjp, float **  srcpos_loc, float ** signals, float ** signals1, int nsrc, float ** absorb_coeff,
                           float *hc, int infoout,int sw, float * K_x_half, float * a_x_half, float * b_x_half,
                           float * K_y_half, float * a_y_half, float * b_y_half,
                           float ** psi_sxx_x, float ** psi_syy_y);

void update_p_PML(int nx1, int nx2, int ny1, int ny2,
                  float **  vx, float **   vy, float **  sp, float ** u, float ** pi, float ** absorb_coeff, float **rho, float *hc, int infoout,
                  float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
                  float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,
                  float ** psi_vxx, float ** psi_vyy, float ** psi_vxy, float ** psi_vyx);

void surface_acoustic_PML(int ndepth, float ** sp);

void exchange_p(float ** sp, float ** bufferlef_to_rig, float ** bufferrig_to_lef, 
                float ** buffertop_to_bot, float ** bufferbot_to_top,
                MPI_Request * req_send, MPI_Request * req_rec);


void read_workflow(char file_in[STRING_SIZE],float *** workflow, int *workflow_lines, char header[STRING_SIZE]);
float ** joint_inversion_grad ( float ** gradiant_1,float ** gradiant_2, float alpha, int joint_type);

void snap_SH(FILE *fp,int nt, int nsnap, float ** vz, float **u, float **pi, float *hc,int ishot);
void apply_workflow(float ** workflow,int workflow_lines,char workflow_header[STRING_SIZE],int *iter,float *F_LOW_PASS,int wavetype_start, int * change_wavetype_iter, int * LBFGS_iter_start);

void eprecond(float ** W, float ** vx, float ** vy);
void eprecond_SH(float ** W, float ** vz);
void eprecond1(float ** We, float ** Ws, float ** Wr, float epsilon);

/* Matrix Operations */
float average_matrix(float ** matrix);
float global_maximum(float ** gradiant_1);
void write_matrix_disk(float ** gradient,char path_name[STRING_SIZE]);
float matrix_product(float ** matrix1, float **matrix2);
void get_local_from_global_matrix(float ** global_matrix,float ** local_matrix);
float ** get_global_from_local_matrix(float ** local_matrix);

/* L-BFGS */
void lbfgs(float **grad1, float **grad2, float **grad3,float Vs_avg,float rho_avg,float Vp_avg, float *bfgsscale, float **bfgsmod, float **bfgsgrad,int bfgsnum,int bfgspar, int iteration, int * LBFGS_iter_start);
void lbfgs_reset(int iter, int bfgsnum, int bfgspar,float ** bfgsmod1, float ** bfgsgrad1, float * bfgsscale1);
void lbfgs_core(int iteration, int N_LBFGS, int NPAR_LBFGS,float ** s_LBFGS, float ** y_LBFGS, float * rho_LBFGS,float *q_LBFGS,float *alpha_LBFGS,float *r_LBFGS);

/* Wolfe condition */
int check_wolfe(float steplength, float misfit_old, float misfit_new, float ** grad_old_vs, float ** grad_new_vs, float ** update_vs, float ** grad_old_rho, float ** grad_new_rho, float ** update_rho, float ** grad_old_vp, float ** grad_new_vp, float ** update_vp, float c1, float c2, int NPAR_LBFGS);
void wolfe_linesearch(int wolfe_status, float *alpha_SL_min, float *alpha_SL_max, float *alpha_SL);

/* functions for viscoacoustic modelling */
void model_viscac(float  **  rho, float **  pi, float **  taup, float *  eta);
void readmod_viscac(float  **  rho, float **  pi, float **  taup, float *  eta);
void matcopy_viscac(float ** prho, float ** ppi, float ** taup);
void prepare_update_p(float *etajm, float *peta, float **ppi, float **prho, float **ptaup, float **g, float *bjm, float *cjm, float ***e);
void zero_fdveps_viscac(int ny1, int ny2, int nx1, int nx2, float ** vx, float ** vy, float ** sp, float ** vxp1, float ** vyp1,
                        float ** psi_sxx_x, float ** psi_sxy_x, float ** psi_vxx, float ** psi_vyx, float ** psi_syy_y, float ** psi_sxy_y, float ** psi_vyy, float ** psi_vxy, float ** psi_vxxs, float ***pp);
void update_p_visc_PML(int nx1, int nx2, int ny1, int ny2, float ** vx, float ** vy, float ** sp, float ** pi, float **rho, float *hc, int infoout,
                       float ***p, float **g, float *bjm, float *cjm, float ***e,
                       float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half,
                       float * K_y, float * a_y, float * b_y, float * K_y_half, float * a_y_half, float * b_y_half,
                       float ** psi_vxx, float ** psi_vyy, float ** psi_vxy, float ** psi_vyx);

