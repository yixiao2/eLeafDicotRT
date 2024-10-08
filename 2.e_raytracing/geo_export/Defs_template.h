/***********************************************
*    update 2021 June
*           - add refr_s
*           - use sac_o and sac_i instead of absorb
*           - use I_absorb_tri_o/i instead of I_absorb_tri
*           - use I_absorb_o/i instead of I_absorb
*
*    update 2021 May
*           - add flag_*
*  
*    version 11/11/11
**********************************************/

#ifndef DEFS_H
#define DEFS_H

typedef double Vec[3];
typedef double Light;	// first we will consider monolight

typedef struct t_tri
{
	int n[3];	// number of the three points
	Vec Nor;	// normal vector for the triangle
	double I_intersect; // not used?
	double I_absorb_tri_o;
	double I_absorb_tri_i;
	int count_hit_tri;
}Triangle;

typedef struct t_pt
{
	Vec coordinate;
	Vec normal;
}Point;

typedef struct t_sur
{
	int pnum;
	int tnum;
	Point *vertex;
	Triangle *tri;
	double z_buf;
}Surface;

typedef struct t_obj
{
	Surface *S;
	double refr_index_o;
	double refr_index_i;
	double refr_index_s; // refractive index of surface, currently active for cell wall only, inactive for chlo and vac surface
	double sac_o;	// index of absorption
	double sac_i;
	double scatter;	// index of scattering, NOT USED, to develop
	double I_absorb_o;
	double I_absorb_i;
	struct t_obj* subobj;
	struct t_obj* nextobj;
	struct t_obj* belongobj;
	int celltype;//0->leaf; 1->ms cellwall; 2->ms chlo; 4->ms vac; 3->nonms cellwall
	int cellname;// leaf surface=0
	int chlname;// cell surface=0
	int vacuolename;// chlo surface=0
	int count_hit;
}Object;

typedef struct t_isect
{
	Vec inter_pt;
	Vec inter_nor;
	Object *inter_obj;
	int inter_tri;
	double t;
}Isect;

typedef struct t_ray
{
	Vec P;	// origin point
	Vec D;	// direction
	Light I;	// light intensity or energy in each beam
	Object *P_obj;
	int P_tri;
	int interflag;
	int level;
	double R;
	double T;
	Object *trace_obj;
	double refr_index;
	int outin_P;// 1=P is on the inner side of obj; 0=P is on the outter surface of obj
}Ray;

#define Infinity 1000//1.7E+308

#define PI 3.141592653589793

#define discardI 0.0001
//Object *GLB_curobj;//,*curbelong,*curobj_refl,*curobj_trans;
//int TIR;
//Isect *GLB_intersect;

double I_discard;
double I_discard_Rf;
double I_discard_Tr;

//#define chl_concentration_MS (2.352941e4) //g m^-3
//#define chl_concentration_BS (2.352941e4) //g m^-3
//#define SAC_chlab 3.27 //m^2 g^-1
//#define SAC_water 0.41
//#define SAC_cellwall (1.47e6)*(4.38e-4)
//#define SAC_cyto 0.41
double chl_con_MS;
double chl_con_BS;
double SAC_chlab;
//double SAC_chlab_BS;
double SAC_water;
double SAC_cellwall;
double SAC_cyto;
double SAC_air;

#define air_refr_index 1
#define wall_refr_index 1.415
#define cyto_refr_index 1.353
#define chlo_refr_index 1.511
#define vacu_refr_index 1.333

#define DBL_EPSILON (2.2204460492503131e-16)
#define DIS_EPSILON (1.0e-11)//in eLeaf, model_ddis=0.1e-6, DIS_EPSILON=1e-4*model_ddis

int num_chl_hit;
double debugI;

#define flag_RTdebug_printf 0//=1 debugmode
#define flag_warningmsg 0//=1 debugmode
#define flag_errormsg 1 //keep 1 even under non-debugging mode
#define flag_RT_file4plot 0//debug mode is more useful
FILE *fout_file4plot;
#define flag_randseed 1 //change to 0 for debug
#define flag_pertube_dir 1 //pertubate incident direction within 1 degree
#define flag_opt_absorbevents 1 // output absorb events for RT training
FILE *fout_absorbevents;
#define flag_output_results_tri 0 //output light absorptance of each triangle in the geometry, notice it will cost a lot of space if run in parallel
#define flag_output_results_srf 1
#define flag_debug_precal 0 //useless now, once used in dev branch
#define flag_RTtrain_mod 1//=1 less output messages for ray tracing module training
int ray_i,ray_j;

//#define xmax 1.610000e-05
//#define xmin 0.000000e+00
//#define zmax 2.731000e-04
//#define zmin -1.710000e-05
//#define ymax 1.610000e-05
//#define ymin 0.000000e+00
//
//#define pal_num 4
//#define spo_num 21
//#define ms_num (pal_num+spo_num)
//int count_ms;
//#define ms_max_chl_num 52
//int ms_chl_num[ms_num];
//double ms_chl_con[ms_num];
//double ms_chl_con_v2[ms_num][ms_max_chl_num];
//int flag_format_count_chl;
//#define epl_l_num 3
//#define epl_u_num 3
//#define nonms_num 6
//int count_nonms;
//
//Object *p_cell_ms[ms_num];
//Object *p_chl_ms[ms_num][ms_max_chl_num];
//Object *p_vac_ms[ms_num];
//Object *p_leaf;
//Object *p_cell_ns[nonms_num];
//
//#endif

