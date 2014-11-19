


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#define MAX_INDS 2000
#define MAX_LOCI 150000
#define MAX_LINE 500000 
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS2 1.2e-7
#define RNMX (1.0-EPS2)
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

struct field{
char locus[20];
int col1;
int col2;
struct field *nextPtr;
};

typedef struct field fields;
typedef fields *fieldsPtr;

/*Haplotypes : up to  loci  ------------------------------------------*/

struct hap{
char allele[MAX_LOCI][6];
double count;
int flag;
struct hap *nextPtr;
};

typedef struct hap haps;
typedef haps *hapsPtr;

/*For pooling haplotype counts - whole file  ------------------------------------------*/

struct all{
char locus1[15];
char locus2[15];
double kb1;
double kb2;
double a;
double b;
double c;
double d;
int flag;
struct all *nextPtr;
};

typedef struct all alls;
typedef alls *allsPtr;



/*For pooling haplotype counts  ------------------------------------------*/
/*
struct poo{
char locus1[15];
char locus2[15];
double kb1;
double kb2;
double a;
double b;
double c;
double d;
struct poo *nextPtr;
};

typedef struct poo poos;
typedef poos *poosPtr;
*/
/*Allele frequencies -------------------------------------------------*/

struct al{
char locus[20];
int index;
int idis;
char allele[6];
double count;
double count2;
double it;
double freq;
struct al *nextPtr;
};
typedef struct al als;
typedef als *alsPtr;

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
double chi;
double freq1;
double freq2;
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
char locus[20];
char locus1[20];
char locus2[20];
char ckb[25];
char ckb1[25];
char ckb2[25];
double kb1;
double kb2;
double aitab[50];
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
char locus[15];
double kb;
int order;
int col1;
int col2;
double freq;
struct map2 *nextPtr;
};
typedef struct map2 mapkb2;
typedef mapkb2 *mapPtr2;

/*Input fields---------------------------------------------------*/
struct fields{
char field[15];
int col1;
int col2;
struct fields *nextPtr;
};
typedef struct fields fields1;
typedef fields1 *fieldPtr;

/*ped data---------------------------------------------------*/
struct ped{
char FAM[15];
char ID[15];
int N;
int GEN[MAX_LOCI][2];
int numg;
struct ped *nextPtr;
};
typedef struct ped peds;
typedef peds *pedPtr;

/*********************************************************************************/
/*FILES.....*/
FILE *output_f;
FILE *output_fi;
FILE *output_f5;
FILE *output_f3;
FILE *output_f2;
FILE *fped;
FILE *output_ft;
FILE *output_pf;
FILE *job_fp;

/*POINTERS.....*/
ldsPtr ldPtr,ldsstartPtr;
fieldsPtr g_p1,fields_startPtr;
hapsPtr hap_startPtr,g_hp1;
alsPtr alstartPtr;
abcdPtr abcdstartPtr,g_abcd;
aisPtr aistartPtr,g_aiPtr,g_ai2Ptr;
allsPtr all_startPtr;
mapPtr2 g_p12,map_startPtr2;
pedPtr ped_startPtr;
fieldPtr g_f1,field_startPtr;

/********************/
double g_lnl,g_location[MAX_LOCI][3];
double g_alllocation[MAX_LOCI][4];
double global_m,global_l,global_e;
double g_kee,g_kll,g_ktt,g_kss,g_kel,g_ket,g_kem,g_kes,g_klt,g_klm,g_kls,g_kts,g_kms;
double g_uei,g_uli,g_umi,g_pred;
double g_meank,g_meanchi,g_startlnl;
double interv[MAX_LOCI][7],g_rho,g_rhoi;
double cinterv[MAX_LOCI][7];
double g_best,g_sumki2,g_sumki,g_cov,g_covi;
double g_seE,g_seL,g_seM,g_bestlik,g_bestE,g_bestL,g_bestM;
double g_bestseE,g_bestseL,g_bestseM;
/********************/
char ped_file[20];
char g_datfile[15];
char g_posneg[6];
char g_temp[MAX_LINE];
char gg_temp[500];
char terfile[20],intefile[20],datfile[20];
char g_loci[MAX_LOCI][20];
char g_allloci[MAX_LOCI][20];
char g_disease[15],g_jobfile[20];

/********************/
int g_bad,g_n,g_nloci,g_nallloci,g_ind,g_niter,g_finish;
int ldflag2,ldflag;
int itlast, iteast,itmast,g_nki;
int g_iter[5], global_fin,g_calls,g_ind;
int ite,itl,itm,itr,itd,g_sig;
int g_estimatel;
/********************/

void chi(double,double);
void final();
void jobin();
void file_input(char *);
void pool(char *);
void file_input2(char *);
void make_map();
void getfields();
void make_hap();
void get_tab();
void getpa(char[500]);
double gammln(double);
double *dvector(long nl, long nh);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
int **imatrix(long nrl, long nrh, long ncl, long nch);
void runewt();
void cal_SE2(double **,int);
void runewt2();
void intfile2();
void intfile3();
#define SQR(a)  (a)*(a)
#define FMAX(a,b) ((a) > (b)? (a): (b))
void nrerror();
double *dvector();
double **dmatrix();
int *ivector();
void free_dmatrix();
double func(double[]);
double bound(int,double[]);
void dfunc(double[],double[]);
void writeint();
void metricho(double,double,double,double);
void diplo(char *);
void matrix();
void fill_str2(char[20],char[20]);
void update();
void ludcmp(double **,int,int *,double *);
void meankchi(char[20]);
void get_predl();
void make_map4();
void readped();
void writ();
void keep_it();
void keep_it2();
void one_interval(int , intsPtr,double *, double *, double *);
void locus_list();
int g_nhole;
double g_max,g_pp,g_ppi,g_sed;
void fast();
hugsPtr bigstart;
int g_npar,g_int;
double g_df,g_V,g_ekeep,g_mkeep;
void directmap();
int g_code,g_opt,path;
void hapallele();
void haplo();
char g_SNP[20];
void quicklike(double,double,double);
int g_initial;
void writeter();
void writeterfin();
void printinputmap();
void reorder();
void reorderall();
void directmap2();
void updat1();
void runewt3();
void file_inputg();
void multallele1();
void diplo1(char *);
void getfreqsh();
void getfreqsh2();
void hapallele2();
void getfreqs();
void multallele();
void segments();
void ldmap2();
void mergeloc();
void lubksb(double **, int, int *, double []);
void intfile4();
void intfileseg();
void ldmapseg();
void writeterfinseg();
void intfile6();
void free_dvector(double *, long, long);
void printmap();
float ran1(long * );
char g_buffer[MAX_INDS][MAX_LINE];
int g_used[MAX_INDS];
long idum;
int g_nrec;
int g_number;

void dfpmin(double p[], int n, double gtol, int *iter, double *fret,
	    double(*func)(double []), void (*dfunc)(double [], double []));

void intfile7();
void file_inputg2();

void simboot();
