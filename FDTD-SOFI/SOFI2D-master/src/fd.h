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
/* $Id: fd.h 865 2015-09-22 12:57:11Z tmetz $ */
/*------------------------------------------------------------------------
 *  fd.h - include file for viscoelastic FD program sofi2D
 *
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

#define PI (3.141592653589793238462643383279502884197169)
#define NPAR 41
//#define STRING_SIZE 74 // previous value, sometimes not enough to handle longer file names
#define STRING_SIZE 256
#define REQUEST_COUNT 4


/* declaration of functions */
void abs_update_s(int i, int j, float **sxx,float **sxy,float **syy, float **absorb_coeff);

void abs_update_v(int i, int j, float **vx,float **vy, float **absorb_coeff);

void absorb(float **absorb_coeff);

void av_mat(float   **pi, float   **u,
            float   **ppijm, float   **puip, float **pujm);

void av_mue(float **u, float **uipjp);

void av_rho(float **rho, float **rip, float **rjp);

void av_tau(float **taus, float **tausipjp);

void check_fs(FILE *fp, int argc, char *fileinp);

void checkfd(FILE *fp, float **prho, float **ppi, float **pu,
             float **ptaus, float **ptaup, float *peta, float *hc, float **srcpos, int nsrc, int **recpos, int ntr);

void catseis(float **data, float **fulldata, int *recswitch, int ntr_glob, int ns);

void comm_ini(float **bufferlef_to_rig, float **bufferrig_to_lef,
              float **buffertop_to_bot, float **bufferbot_to_top,
              MPI_Request *req_send, MPI_Request *req_rec);

void cpml_update_s_x(int i, int j,float   *vxx, float *vyx,float *K_x, float *a_x,
                     float *b_x, float *K_x_half, float *a_x_half, float *b_x_half ,float **psi_vxx,float **psi_vyx);

void cpml_update_s_y(int i, int j,float *vxy,float *vyy,float *K_y, float *a_y,
                     float *b_y, float *K_y_half, float *a_y_half, float *b_y_half ,float **psi_vyy,float **psi_vxy);

void cpml_update_v_x(int i, int j,float   *sxx_x, float *sxy_x,float *K_x, float *a_x,float *b_x,
                     float *K_x_half, float *a_x_half, float *b_x_half ,float **psi_sxx_x,float **psi_sxy_x);

void cpml_update_v_y(int i, int j,float *sxy_y,float *syy_y,float *K_y, float *a_y,
                     float *b_y, float *K_y_half, float *a_y_half, float *b_y_half ,float **psi_syy_y,float **psi_sxy_y);

void exchange_s_rsg(float **vx, float **vy, float **vz,
                    float **bufferlef_to_rig, float **bufferrig_to_lef,
                    float **buffertop_to_bot, float **bufferbot_to_top);

void exchange_v(int nd, float **vx, float **vy,
                float **bufferlef_to_rig, float **bufferrig_to_lef,
                float **buffertop_to_bot, float **bufferbot_to_top,
                MPI_Request *req_send, MPI_Request *req_rec);

void exchange_s(int nd, float **sxx, float **syy,
                float **sxy, float **bufferlef_to_rig, float **bufferrig_to_lef,
                float **buffertop_to_bot, float **bufferbot_to_top,
                MPI_Request *req_send, MPI_Request *req_rec);


void exchange_par(void);

float *holbergcoeff(void);

void info(FILE *fp);

void initproc(void);

void interpol(int ni1, int ni2, float   **intvar, int cfgt_check);

void model_visco(float    **rho, float   **pi, float   **u,
                 float   **taus, float   **taup, float   *eta);

void model_elastic(float    **rho, float   **pi, float   **u);

void model_ani(float    **rho, float   **c11, float   **c15, float   **c13,
               float   **c35, float   **c33, float   **c55,
               float   **taus, float   **taup, float   *eta);

void matcopy(float **prho, float **ppi, float **pu, float **ptaup,
             float **ptaus);

void matcopy_elastic(float **prho, float **ppi, float **pu);

void matcopy_ani(float **rho, float   **c11, float   **c15, float   **c13,
                 float   **c35, float   **c33, float   **c55, float **taus,
                 float **taup);

void merge(int nsnap, int type);

void mergemod(char modfile[STRING_SIZE], int format);

void note(FILE *fp);

void  outseis(FILE *fp, FILE *fpdata, float **section,
              int **recpos, int **recpos_loc, int ntr, float **srcpos_loc,
              int nsrc, int ns, int seis_form, int ishot);

void  outseis_glob(FILE *fp, FILE *fpdata, float **section,
                   int **recpos, int **recpos_loc, int ntr, float **srcpos_loc,
                   int nsrc, int ns, int seis_form, int ishot, int comp);

void  output_source_signal(FILE *fp, float **signals, int ns, int seis_form);

void operator_s_fd2(int i, int j,float   *vxx, float *vyx,float *vxy,
                    float *vyy, float **vx, float **vy,float *hc);

void operator_s_fd4(int i, int j,float   *vxx, float *vyx,float *vxy,
                    float *vyy, float **vx, float **vy,float *hc);

void operator_s_fd6(int i, int j,float   *vxx, float *vyx,float *vxy,
                    float *vyy, float **vx, float **vy,float *hc);

void operator_s_fd8(int i, int j,float   *vxx, float *vyx,float *vxy,
                    float *vyy, float **vx, float **vy,float *hc);

void operator_s_fd10(int i, int j,float   *vxx, float *vyx,float *vxy,
                     float *vyy, float **vx, float **vy,float *hc);

void operator_s_fd12(int i, int j,float   *vxx, float *vyx,float *vxy,
                     float *vyy, float **vx, float **vy,float *hc);

void operator_v_fd2(int i, int j,float   *sxx_x, float *sxy_x,float *sxy_y,
                    float *syy_y, float **sxx, float **syy,float **sxy, float *hc);

void operator_v_fd4(int i, int j,float   *sxx_x, float *sxy_x,float *sxy_y,
                    float *syy_y, float **sxx, float **syy,float **sxy, float *hc);

void operator_v_fd6(int i, int j,float   *sxx_x, float *sxy_x,float *sxy_y,
                    float *syy_y, float **sxx, float **syy,float **sxy, float *hc);

void operator_v_fd8(int i, int j,float   *sxx_x, float *sxy_x,float *sxy_y,
                    float *syy_y, float **sxx, float **syy,float **sxy, float *hc);

void operator_v_fd10(int i, int j,float   *sxx_x, float *sxy_x,float *sxy_y,
                     float *syy_y, float **sxx, float **syy,float **sxy, float *hc);

void operator_v_fd12(int i, int j,float   *sxx_x, float *sxy_x,float *sxy_y,
                     float *syy_y, float **sxx, float **syy,float **sxy, float *hc);

void par_mult_dt(float **pi, float **u, float **uipjp);

void psource(int nt, float **sxx, float **syy,
             float   **srcpos_loc, float **signals, int nsrc);

void psource_rsg(int nt, float **sxx, float **syy,
                 float   **srcpos_loc, float **signals, int nsrc);


float *rd_sour(int *nts,FILE *fp_source);

float readdsk(FILE *fp_in, int format);

void readbufs(float **sxx, float **syy,
              float **sxy, float **bufferlef_to_rig, float **bufferrig_to_lef,
              float **buffertop_to_bot, float **bufferbot_to_top);

void readbufv(float **vx, float **vy,
              float **bufferlef_to_rig, float **bufferrig_to_lef,
              float **buffertop_to_bot, float **bufferbot_to_top);

void read_checkpoint(int nx1, int nx2, int ny1, int ny2,
                     float   **vx, float **vy, float **sxx, float **syy, float **sxy);

void read_par_json(FILE *fp, char *fileinp);


void readmod_visco(float    **rho, float   **pi, float   **u,
                   float   **taus, float   **taup, float   *eta);

void readmod_elastic(float    **rho, float   **pi, float   **u);

int **receiver(FILE *fp, int *ntr);

void save_checkpoint(int nx1, int nx2, int ny1, int ny2,
                     float   **vx, float **vy, float **sxx, float **syy, float **sxy);

void saveseis(FILE *fp, float **sectionvx, float **sectionvy,float **sectionp,
              float **sectioncurl, float **sectiondiv, int  **recpos, int  **recpos_loc,
              int ntr, float **srcpos_loc, int nsrc,int ns);

void saveseis_glob(FILE *fp, float **sectiondata, int  **recpos, int  **recpos_loc,
                   int ntr, float **srcpos, int ishot,int ns, int sectiondatatype);

void snap(FILE *fp,int nt, int nsnap, float **vx, float **vy, float **sxx,
          float **syy, float **u, float **pi, float *hc);


void snap_rsg(FILE *fp,int nt, int nsnap, float **vx, float **vy, float **sxx, float **syy, float **u, float **pi);

void snapmerge(int nsnap);

float **sources(int *nsrc);

void seismo(int lsamp, int ntr, int **recpos, float **sectionvx,
            float **sectionvy, float **sectionp, float **sectioncurl, float **sectiondiv,
            float **pvx, float **pvy, float **psxx, float **psyy, float **ppi, float **pu);

void seismo_ssg(int lsamp, int ntr, int **recpos, float **sectionvx,
                float **sectionvy, float **sectionp, float **sectioncurl, float **sectiondiv,
                float **pvx, float **pvy, float **psxx, float **psyy, float **ppi, float **pu,
                float *hc);

void seismo_rsg(int lsamp, int ntr, int **recpos, float **sectionvx,
                float **sectionvy, float **sectionp, float **sectioncurl, float **sectiondiv,
                float **pvx, float **pvy, float **psxx, float **psyy, float **ppi, float **pu);

int **splitrec(int **recpos,int *ntr_loc, int ntr, int *recswitch);

float **splitsrc(float **srcpos,int *nsrc_loc, int nsrc);

void subgrid_bounds(int nx1, int nx2, int ny1, int ny2, int *gx, int *gy);

void surface(int ndepth, float **pvx, float **pvy,
             float **psxx, float **psyy,
             float **psxy, float *** pp, float *** pq,
             float    **ppi, float    **pu, float **ptaup,
             float **ptaus, float *etajm, float *peta, float *hc, float *K_x, float *a_x, float *b_x, float **psi_vxx);

void surface_elastic(int ndepth, int *gx, float **vx, float **vy, float **sxx, float **syy,
                     float **sxy, float    **pi, float    **u, float *hc, float *K_x, float *a_x, float *b_x, float**
                     psi_vxx);


void prepare_update_s(float *etajm, float *etaip, float *peta, float **fipjp, float **pu,
                      float **puipjp, float **ppi, float **ptaus, float **ptaup,
                      float **ptausipjp, float **f, float **g, float *bip, float *bjm,
                      float *cip, float *cjm, float ***dip, float ***d, float ***e);

void prepare_update_s_4(float *etajm, float *etaip, float *peta, float **fipjp, float **pu,
                        float **puipjp, float **ppi, float **ptaus, float **ptaup,
                        float **ptausipjp, float **f, float **g, float *bip, float *bjm,
                        float *cip, float *cjm, float ***dip, float ***d, float ***e);

void update_s_visc_abs(int nx1, int nx2, int ny1, int ny2, int *gx, int *gy, int nt,
                       float **vx, float **vy, float **sxx, float **syy,
                       float **sxy, float ***r, float *** p, float ***q,
                       float **pi,
                       float **fipjp, float **f, float **g, float *bip, float *bjm, float *cip,
                       float *cjm, float ***d, float ***e, float ***dip,
                       float **absorb_coeff,float *hc);

void update_s_visc_abs_4(int nx1, int nx2, int ny1, int ny2, int *gx, int *gy, int nt,
                         float **vx, float **vy, float **sxx, float **syy,
                         float **sxy, float ***r, float *** p, float ***q,
                         float **pi,
                         float **fipjp, float **f, float **g, float *bip, float *bjm, float *cip,
                         float *cjm, float ***d, float ***e, float ***dip,
                         float **absorb_coeff,float *hc,  float **vxx_1,float **vxx_2,float **vxx_3,float **vxx_4,float **vyy_1,float **vyy_2,float **vyy_3,float **vyy_4,float **vxy_1,float **vxy_2,float **vxy_3,float **vxy_4,float **vyx_1,float **vyx_2,float **vyx_3,float **vyx_4,float **svx_1,float **svx_2,float **svx_3,float **svx_4,float **svy_1,float **svy_2,float **svy_3,float **svy_4,float ***pr_2,float ***pr_3,float ***pr_4, float ***pp_2, float ***pp_3, float ***pp_4, float ***pq_2, float ***pq_3, float ***pq_4);

void update_s_elastic(int nx1, int nx2, int ny1, int ny2, int nt,
                      float   **vx, float    **vy, float    **sxx, float    **syy,
                      float    **sxy, float **pi, float **u, float **uipjp, float **absorb_coeff,
                      float *hc);

void update_s_elastic_abs(int nx1, int nx2, int ny1, int ny2, int *gx, int *gy, int nt,
                          float   **vx, float    **vy, float    **sxx, float    **syy,
                          float    **sxy, float **pi, float **u, float **uipjp,
                          float **absorb_coeff, float *hc);

void update_s_elastic_abs_4(int nx1, int nx2, int ny1, int ny2, int *gx, int *gy, int nt,
                            float   **vx, float    **vy, float    **sxx, float    **syy,
                            float    **sxy, float **pi, float **u, float **uipjp,
                            float **absorb_coeff, float *hc ,float **vxx_1,float **vxx_2,float **vxx_3,float **vxx_4,float **vyy_1,float **vyy_2,float **vyy_3,float **vyy_4,float **vxy_1,float **vxy_2,float **vxy_3,float **vxy_4,float **vyx_1,float **vyx_2,float **vyx_3,float **vyx_4);

void update_s_elastic_interior(int nx1, int nx2, int ny1, int ny2, int *gx, int *gy, int nt,
                               float   **vx, float    **vy, float    **sxx, float    **syy,
                               float    **sxy, float **pi, float **u, float **uipjp,
                               float *hc);

void update_s_elastic_interior_4(int nx1, int nx2, int ny1, int ny2, int *gx, int *gy, int nt,
                                 float   **vx, float    **vy, float    **sxx, float    **syy,
                                 float    **sxy, float **pi, float **u, float **uipjp, float *hc,float **vxx_1,float **vxx_2,float **vxx_3,float **vxx_4,float **vyy_1,float **vyy_2,float **vyy_3,float **vyy_4,float **vxy_1,float **vxy_2,float **vxy_3,float **vxy_4,float **vyx_1,float **vyx_2,float **vyx_3,float **vyx_4);

void update_s_elastic_PML(int nx1, int nx2, int ny1, int ny2, int *gx, int *gy, int nt,
                          float   **vx, float    **vy, float    **sxx, float    **syy,
                          float    **sxy, float **pi, float **u, float **uipjp, float *hc,
                          float *K_x, float *a_x, float *b_x, float *K_x_half, float *a_x_half, float *b_x_half,
                          float *K_y, float *a_y, float *b_y, float *K_y_half, float *a_y_half, float *b_y_half,
                          float **psi_vxx, float **psi_vyy, float **psi_vxy, float **psi_vyx);

void update_s_elastic_PML_4(int nx1, int nx2, int ny1, int ny2, int *gx, int *gy, int nt,
                            float   **vx, float    **vy, float    **sxx, float    **syy,
                            float    **sxy, float **pi, float **u, float **uipjp, float *hc,
                            float *K_x, float *a_x, float *b_x, float *K_x_half, float *a_x_half, float *b_x_half,
                            float *K_y, float *a_y, float *b_y, float *K_y_half, float *a_y_half, float *b_y_half,
                            float **psi_vxx, float **psi_vyy, float **psi_vxy, float **psi_vyx,float **vxx_1,float **vxx_2,float **vxx_3,float **vxx_4,float **vyy_1,float **vyy_2,float **vyy_3,float **vyy_4,float **vxy_1,float **vxy_2,float **vxy_3,float **vxy_4,float **vyx_1,float **vyx_2,float **vyx_3,float **vyx_4);

void update_s_visc_interior(int nx1, int nx2, int ny1, int ny2, int *gx, int *gy, int nt,
                            float **vx, float **vy, float **sxx, float **syy,
                            float **sxy, float ***r, float *** p, float ***q,
                            float **fipjp, float **f, float **g, float *bip, float *bjm, float *cip,
                            float *cjm, float ***d, float ***e, float ***dip,
                            float *hc);

void update_s_visc_PML(int nx1, int nx2, int ny1, int ny2, int *gx, int *gy,int nt,
                       float   **vx, float    **vy, float    **sxx, float    **syy,
                       float    **sxy, float *hc, float ***r, float ***p, float ***q, float **fipjp, float **f, float **g,
                       float *bip, float *bjm, float *cip, float *cjm, float ***d, float ***e, float ***dip,
                       float *K_x, float *a_x, float *b_x, float *K_x_half, float *a_x_half, float *b_x_half,
                       float *K_y, float *a_y, float *b_y, float *K_y_half, float *a_y_half, float *b_y_half,
                       float **psi_vxx, float **psi_vyy, float **psi_vxy, float **psi_vyx);

void update_s_visc_interior_4(int nx1, int nx2, int ny1, int ny2, int *gx, int *gy, int nt,
                              float **vx, float **vy, float **sxx, float **syy,
                              float **sxy, float ***r, float *** p, float ***q,
                              float **fipjp, float **f, float **g, float *bip, float *bjm, float *cip,
                              float *cjm, float ***d, float ***e, float ***dip,
                              float *hc,  float **vxx_1,float **vxx_2,float **vxx_3,float **vxx_4,float **vyy_1,float **vyy_2,float **vyy_3,float **vyy_4,float **vxy_1,float **vxy_2,float **vxy_3,float **vxy_4,float **vyx_1,float **vyx_2,float **vyx_3,float **vyx_4,float **svx_1,float **svx_2,float **svx_3,float **svx_4,float **svy_1,float **svy_2,float **svy_3,float **svy_4,float ***pr_2,float ***pr_3,float ***pr_4, float ***pp_2, float ***pp_3, float ***pp_4, float ***pq_2, float ***pq_3, float ***pq_4);

void update_s_visc_PML_4(int nx1, int nx2, int ny1, int ny2, int *gx, int *gy,int nt,
                         float   **vx, float    **vy, float    **sxx, float    **syy,
                         float    **sxy, float *hc, float ***r, float ***p, float ***q, float **fipjp, float **f, float **g,
                         float *bip, float *bjm, float *cip, float *cjm, float ***d, float ***e, float ***dip,
                         float *K_x, float *a_x, float *b_x, float *K_x_half, float *a_x_half, float *b_x_half,
                         float *K_y, float *a_y, float *b_y, float *K_y_half, float *a_y_half, float *b_y_half,
                         float **psi_vxx, float **psi_vyy, float **psi_vxy, float **psi_vyx,  float **vxx_1,float **vxx_2,float **vxx_3,float **vxx_4,float **vyy_1,float **vyy_2,float **vyy_3,float **vyy_4,float **vxy_1,float **vxy_2,float **vxy_3,float **vxy_4,float **vyx_1,float **vyx_2,float **vyx_3,float **vyx_4,float **svx_1,float **svx_2,float **svx_3,float **svx_4,float **svy_1,float **svy_2,float **svy_3,float **svy_4,float ***pr_2,float ***pr_3,float ***pr_4, float ***pp_2, float ***pp_3, float ***pp_4, float ***pq_2, float ***pq_3, float ***pq_4);

void PML_pro(float *d_x, float *K_x, float *alpha_prime_x, float *a_x, float *b_x,
             float *d_x_half, float *K_x_half, float *alpha_prime_x_half, float *a_x_half, float *b_x_half,
             float *d_y, float *K_y, float *alpha_prime_y, float *a_y, float *b_y,
             float *d_y_half, float *K_y_half, float *alpha_prime_y_half, float *a_y_half, float *b_y_half);


void update_v(int nx1, int nx2, int ny1, int ny2, int nt,
              float   **pvx, float **pvy, float **psxx, float **psyy,
              float **psxy, float **prho, float  **prip, float **prjp,
              float   **srcpos_loc, float **signals, int nsrc, float **absorb_coeff,
              float *hc);

void update_v_abs(int nx1, int nx2, int ny1, int ny2, int *gx, int *gy, int nt, float   **vx, float **vy,
                  float **sxx, float **syy, float **sxy,  float  **rip, float **rjp,
                  float **absorb_coeff,float *hc);

void update_v_abs_4(int nx1, int nx2, int ny1, int ny2, int *gx, int *gy, int nt,
                    float   **vx, float **vy, float **sxx, float **syy, float **sxy,
                    float  **rip, float **rjp, float **absorb_coeff,float *hc,float **svx_1,float **svx_2,float **svx_3,float **svx_4,float **svy_1,float **svy_2,float **svy_3,float **svy_4);

void update_v_interior(int nx1, int nx2, int ny1, int ny2, int *gx, int *gy, int nt,
                       float   **vx, float **vy, float **sxx, float **syy,
                       float **sxy, float **rho, float  **rip, float **rjp,
                       float   **srcpos_loc, float **signals, int nsrc,
                       float *hc);

void update_v_interior_4(int nx1, int nx2, int ny1, int ny2, int *gx, int *gy, int nt,
                         float   **vx, float **vy, float **sxx, float **syy,
                         float **sxy, float **rho, float  **rip, float **rjp,
                         float   **srcpos_loc, float **signals, int nsrc,float *hc,float **svx_1,float **svx_2,float **svx_3,float **svx_4,float **svy_1,float **svy_2,float **svy_3,float **svy_4);

void update_v_PML(int nx1, int nx2, int ny1, int ny2,int *gx, int *gy, int nt, float   **vx, float **vy,
                  float **sxx, float **syy, float **sxy, float  **rip, float **rjp,
                  float *hc, float *K_x, float *a_x, float *b_x, float *K_x_half,
                  float *a_x_half, float *b_x_half, float *K_y, float *a_y, float *b_y, float *K_y_half, float *a_y_half, float *b_y_half, float **psi_sxx_x, float **psi_syy_y,
                  float **psi_sxy_y, float **psi_syx_x);
void update_v_PML_4(int nx1, int nx2, int ny1, int ny2, int *gx, int *gy, int nt, float   **vx, float **vy,
                    float **sxx, float **syy, float **sxy,  float  **rip, float **rjp,
                    float *hc, float *K_x, float *a_x, float *b_x, float *K_x_half, float *a_x_half,
                    float *b_x_half, float *K_y, float *a_y, float *b_y, float *K_y_half, float *a_y_half, float *b_y_half,
                    float **psi_sxx_x, float **psi_syy_y, float **psi_sxy_y, float **psi_sxy_x,float **svx_1,float **svx_2,float **svx_3,float **svx_4,float **svy_1,float **svy_2,float **svy_3,float **svy_4);

void wavefield_update_s_el(int i, int j,float   vxx, float  vyx,float vxy,float  vyy, float **sxy,
                           float **sxx, float **syy, float **pi, float **u, float **uipjp);

void wavefield_update_s_visc(int i, int j,float   vxx, float  vyx,float vxy,float  vyy, float **sxy,
                             float **sxx, float **syy, float ***r, float ***p,
                             float ***q,float **fipjp, float **f, float **g, float *bip,
                             float *bjm,float *cip, float *cjm, float ***d, float ***e, float ***dip);

void wavefield_update_s_visc_4(int i, int j,float   vxx, float  vyx,float vxy,float  vyy, float **sxy,
                               float **sxx, float **syy, float ***r, float ***p,
                               float ***q,float **fipjp, float **f, float **g, float *bip,
                               float *bjm,float *cip, float *cjm, float ***d, float ***e, float ***dip,  float **vxx_1,float **vxx_2,float **vxx_3,float **vxx_4,float **vyy_1,float **vyy_2,float **vyy_3,float **vyy_4,float **vxy_1,float **vxy_2,float **vxy_3,float **vxy_4,float **vyx_1,float **vyx_2,float **vyx_3,float **vyx_4,float **svx_1,float **svx_2,float **svx_3,float **svx_4,float **svy_1,float **svy_2,float **svy_3,float **svy_4,float ***pr_2,float ***pr_3,float ***pr_4, float ***pp_2, float ***pp_3, float ***pp_4, float ***pq_2, float ***pq_3, float ***pq_4);

void wavefield_update_v(int i, int j,float   sxx_x, float  sxy_x,float sxy_y,float  syy_y, float **vx,
                        float **vy, float **rip, float **rjp);
void wavefield_update_v_4(int i, int j,float   sxx_x, float  sxy_x,float sxy_y,float  syy_y, float **vx,
                          float **vy, float **rip, float **rjp,float **svx_1,float **svx_2,float **svx_3,float **svx_4,float **svy_1,float **svy_2,float **svy_3,float **svy_4);

void wavefield_update_s_el_4(int i, int j,float   vxx, float  vyx,float vxy,float  vyy, float **sxy,
                             float **sxx, float **syy, float **pi, float **u, float **uipjp,float **vxx_1,float **vxx_2,float **vxx_3,float **vxx_4,float **vyy_1,float **vyy_2,float **vyy_3,float **vyy_4,float **vxy_1,float **vxy_2,float **vxy_3,float **vxy_4,float **vyx_1,float **vyx_2,float **vyx_3,float **vyx_4);

float **wavelet(float **srcpos_loc, int nsrc);

void writebufs(float **sxx, float **syy,
               float **sxy, float **bufferlef_to_rig, float **bufferrig_to_lef,
               float **buffertop_to_bot, float **bufferbot_to_top);

void writebufv(float **vx, float **vy,
               float **bufferlef_to_rig, float **bufferrig_to_lef,
               float **buffertop_to_bot, float **bufferbot_to_top);

void write_par(FILE *fp);

void writedsk(FILE *fp_out, float amp, int format);

void writemod(char modfile[STRING_SIZE], float **array, int format);

void zero_elastic(int nx1, int nx2, int ny1, int ny2, float **vx, float **vy, float **sxx, float **syy, float **sxy);

void zero_elastic_4(int nx1, int nx2, int ny1, int ny2, float **vxx_1,float **vxx_2,float **vxx_3,float **vxx_4,float **vyy_1,float **vyy_2,float **vyy_3,float **vyy_4,float **vxy_1,float **vxy_2,float **vxy_3,float **vxy_4,float **vyx_1,float **vyx_2,float **vyx_3,float **vyx_4,float **svx_1,float **svx_2,float **svx_3,float **svx_4,float **svy_1,float **svy_2,float **svy_3,float **svy_4);

void zero_visco_4(int nx1, int nx2, int ny1, int ny2, float **vxx_1,float **vxx_2,float **vxx_3,float **vxx_4,float **vyy_1,float **vyy_2,float **vyy_3,float **vyy_4,float **vxy_1,float **vxy_2,float **vxy_3,float **vxy_4,float **vyx_1,float **vyx_2,float **vyx_3,float **vyx_4,float **svx_1,float **svx_2,float **svx_3,float **svx_4,float **svy_1,float **svy_2,float **svy_3,float **svy_4,float ***pr_2,float ***pr_3,float ***pr_4, float ***pp_2, float ***pp_3, float ***pp_4, float ***pq_2, float ***pq_3, float ***pq_4);

void zero_visc(int nx1, int nx2, int ny1, int ny2, float **vx, float **vy, float **sxx, float **syy, float **sxy, float *** pr, float *** pp, float *** pq);

void zero_PML_elastic(int ny1, int ny2, int nx1, int nx2, float **vx, float **vy, float **sxx,
                      float **syy, float **sxy,
                      float **psi_sxx_x, float **psi_sxy_x, float **psi_vxx, float **psi_vyx, float **psi_syy_y, float **psi_sxy_y, float **psi_vyy, float **psi_vxy,float **psi_vxxs);

void zero_PML_visc(int ny1, int ny2, int nx1, int nx2, float **vx, float **vy, float **sxx,
                   float **syy, float **sxy,
                   float **psi_sxx_x, float **psi_sxy_x, float **psi_vxx, float **psi_vyx, float **psi_syy_y, float **psi_sxy_y, float **psi_vyy, float **psi_vxy, float **psi_vxxs,
                   float ***pr, float ***pp, float ***pq);

/* declaration of functions for parser*/

/* declaration of functions for json parser in json_parser.c*/
int read_objects_from_intputfile(FILE *fp, char input_file[STRING_SIZE],char **varname_list,char **value_list);

void print_objectlist_screen(FILE *fp, int number_readobject,char **varname_list,char **value_list);

int count_occure_charinstring(char stringline[STRING_SIZE], char teststring[]);

void copy_str2str_uptochar(char string_in[STRING_SIZE], char string_out[STRING_SIZE], char teststring[]);

int get_int_from_objectlist(char string_in[STRING_SIZE], int number_readobject, int *int_buffer,
                            char **varname_list,char **value_list);

int get_float_from_objectlist(char string_in[STRING_SIZE], int number_readobject, float *double_buffer,
                              char **varname_list,char **value_list);

int get_string_from_objectlist(char string_in[STRING_SIZE], int number_readobject, char string_buffer[STRING_SIZE],
                               char **varname_list,char **value_list);

int is_string_blankspace(char string_in[STRING_SIZE]);

void remove_blankspaces_around_string(char string_in[STRING_SIZE]);

void add_object_tolist(char string_name[STRING_SIZE],char string_value[STRING_SIZE], int *number_read_object,
                       char **varname_list,char **value_list);


/* utility functions */
void declare_error(char err_text[]);
void err2(char errformat[],char errfilename[]);
void warning(char warn_text[]);

double maximum(float **a, int nx, int ny);
float *vector(int nl, int nh);
int *ivector(int nl, int nh);
float **matrix(int nrl, int nrh, int ncl, int nch);
int **imatrix(int nrl, int nrh, int ncl, int nch);
float ** *f3tensor(int nrl, int nrh, int ncl, int nch,int ndl, int ndh);
void free_vector(float *v, int nl, int nh);
void free_ivector(int *v, int nl, int nh);
void free_matrix(float **m, int nrl, int nrh, int ncl, int nch);
void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch);
void free_f3tensor(float ***t, int nrl, int nrh, int ncl, int nch, int ndl,
                   int ndh);

