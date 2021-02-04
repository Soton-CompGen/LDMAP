#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#define MAX_LOCI 1000000
#define MAX_LINE 50000 
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS2 1.2e-7
#define RNMX (1.0-EPS2)
#define SQR(a)  (a)*(a)
#define FMAX(a,b) ((a) > (b)? (a): (b))
/*>>>>>>>>DEFINE GLOBAL VARIABLES, STRUCTURES, & TYPES<<<<<<<<<*/

/*LOCI AND LDU DATA FOR MERGING  ------------------------------------------*/

struct ld{
char locus[20];
double xkb;
double xldu;
char kb[20];
char ldu[20];
int flag;
struct ld *nextPtr;
};

typedef struct ld lds;
typedef lds *ldsPtr;

/*FIELDS---------------------------------------------------*/

/*abcd table ----------------------------------------------------------*/

struct ab{
char locus[20];
char locus2[20];
double p;
double k;
double kb;
double kb1;
double kb2;
double ldu;
double ldu1;
double ldu2;
double n;
int i1;
int i2;
int flag;
struct ab *nextPtr;
};
typedef struct ab abcds;
typedef abcds *abcdPtr;

/****************************************************/
/*atoi table (of 9 counts for diallelic diplotypes) */

struct ai{
char locus1[20];
char locus2[20];
double kb1;
double kb2;
int aitab[15];
double n;
double freq1;
double freq2;
struct ai *nextPtr;
};
typedef struct ai ais;
typedef ais *aisPtr;
/***************************************************/
/****************************************************/
/*atoi table (of 9 counts for diallelic diplotypes) */

/*POINTER BY INTERVAL TO INT FILE */


struct inty{
abcdPtr p;
int    left;
int    right;
double k;
double ab2p;
double dd;
double ed;
struct inty *nextPtr;
};
typedef struct inty ints;
typedef ints *intsPtr;

/*POINTERS TO ABOVE */

struct hug{
intsPtr pp;
struct hug *nextPtr;
};
typedef struct hug hugs;
typedef hugs *hugsPtr;

/*Mapc data---------------------------------------------------*/

/*Mapc data---------------------------------------------------*/
struct map2{
char locus[20];
double kb;
int order;
double freq;
double chi;
struct map2 *nextPtr;
};
typedef struct map2 mapkb2;
typedef mapkb2 *mapPtr2;


/*ped data---------------------------------------------------*/
struct ped{
int GEN[MAX_LOCI][2];
int numg;
struct ped *nextPtr;
};
typedef struct ped peds;
typedef peds *pedPtr;

/*********************************************************************************/
/*FILES.....*/
FILE *output_f;
FILE *output_ft;
FILE *job_fp;
FILE *output_f2;
FILE *fped;

/*POINTERS.....*/
ldsPtr ldPtr,ldsstartPtr;
abcdPtr abcdstartPtr,g_abcd;
aisPtr aistartPtr,g_aiPtr,g_ai2Ptr;
mapPtr2 g_p12,map_startPtr2;
pedPtr g_ped1,ped_startPtr;
hugsPtr bigstart;

/********************/
double g_maf, g_hwp;
double g_flag, g_chi1,g_prob,global_p, g_lnl,g_location[MAX_LOCI][3];
double g_alllocation[MAX_LOCI][4];
double global_m,global_l,global_e;
double g_kel,g_kem,g_klm;
double g_uei,g_uli,g_umi,g_pred;
double g_startlnl;
double interv[MAX_LOCI][7],g_rho,g_rhoi;
double cinterv[MAX_LOCI][7];
double g_best,g_sumki2,g_sumki;
double g_max,g_sed;
double g_df,g_V,g_ekeep,g_mkeep;
/********************/
char g_temp[MAX_LINE];
char gg_temp[500];
char terfile[20],intefile[20],datfile[20];
char g_loci[MAX_LOCI][20];
char g_allloci[MAX_LOCI][20];

/********************/
int g_bad,g_n,g_nloci,g_nallloci,g_niter,g_finish;
int ldflag;
int itlast, iteast,itmast,g_nki;
int g_iter[5], global_fin,g_calls,g_ind;
int ite,itl,itm,g_sig;
int g_estimatel;
int g_opt,path;
int g_npar,g_int;
int g_nhole;
int g_initial;
float ran1(long * );
long idum;
/********************/
double **dmatrix(long nrl, long nrh, long ncl, long nch);
void jobin();
void getpa(char[500]);
int *ivector();
double bound(int,double[]);
void free_dmatrix();
void cal_SE2(double **,int);
void get_tab();
void getfields();
double func(double[]);
void dfunc(double[],double[]);
double *dvector(long nl, long nh);
void nrerror();
void runewt2();
void diplo(char *);
void ludcmp(double **,int,int *,double *);
void fill_gai(int, int, int, int, int, int, int, int, int, double, double, char[20], char[20],double, double);
void make_map4();
void readped();
void metricho(double,double,double,double);
void meankchi(char[20]);
void update();
void one_interval(int , intsPtr,double *, double *, double *);
void keep_it2();
void fast();
void quicklike(double,double,double);
void keep_it();
void get_predl();
void locus_list();
void directmap2();
void updat1();
void reorderall();
void runewt3();
void reorder();
void newentry();
void newentry2();
void file_inputnew();
void getfreqs();
void multallele();
void segments();
void lubksb(double **, int, int *, double []);
void intfile4();
void intfileseg();
void ldmapseg();
void writeterfinseg();
void intfile6();
void free_dvector(double *, long, long);
void dfpmin(double p[], int n, double gtol, int *iter, double *fret,
	    double(*func)(double []), void (*dfunc)(double [], double []));
double gammln(double);
void gcf(double *,double,double,double *,double *);
double gammp(double,double);
void gser(double *,double,double,double *,double *);
void chi(double,double);
