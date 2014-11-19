#include "allass.h"
int main(int argc, char *argv[])
{
char datfile2[50],temp[20],outputfile[20],jobfile[20],flagg[5];
/* time_t start; */
g_int=100;
g_max=500.;
g_finish=0;
g_calls=0;
ldflag2=0;
if(argc==2)strcpy(flagg,argv[1]);
if(argc==3)strcpy(flagg,argv[2]);
if(argc==4)strcpy(flagg,argv[3]);
if(argc==5)strcpy(flagg,argv[4]);
if(argc==6)strcpy(flagg,argv[5]);
if(argc==7)strcpy(flagg,argv[6]);

g_opt=atoi(flagg);

if(g_opt==1||g_opt==2) {strcpy(datfile,argv[1]);strcpy(intefile,argv[2]);}

if(g_opt!=1&&g_opt!=2&&g_opt!=3)
{
if(argc==3){strcpy(datfile,argv[1]);}
if(argc==4){strcpy(intefile,argv[1]);strcpy(outputfile,argv[2]);}
if(argc==5){strcpy(jobfile,argv[1]);strcpy(intefile,argv[2]);strcpy(outputfile,argv[3]);}
if(argc==6){strcpy(jobfile,argv[1]);strcpy(datfile,argv[2]);strcpy(intefile,argv[3]);strcpy(outputfile,argv[4]);}
}

/*////////////////////////////////////////////////////////////////////////////*/

/***OPTION 1 AND 2 ******/

/*Make int files ******/
if(g_opt==1||g_opt==2) 
  { 

printf("\n**NOTE**: The current settings of the program can accomodate up to %d loci and files up to %d characters wide ",MAX_LOCI,MAX_LINE);
strcpy(ped_file,datfile);
strcpy(datfile2,"temp.dat");
if((output_f3=fopen(datfile2,"w"))==NULL){printf("\nCannot open temp.dat file");exit(1);} 
fprintf(output_f3,"\n\n\nASSOCIATION IN SNP DIPLOTYPES (LDMAP)\n\n\n ");
fflush(output_f3);


if(g_opt==1)
{
printf("\nAre alleles coded as 0,1 or 1,2 ? (if 0,1 enter 0, if 1,2 enter 1)"); 
scanf("%s",temp);
if(temp[0]=='1')g_code=0;
if(temp[0]=='0')g_code=1;
}
file_inputg();

if(g_opt==1) { if(g_code==0){getfreqsh2();hapallele();} if(g_code==1){getfreqsh();hapallele2();} }

if(g_opt==2){getfreqs();multallele();}
fflush(output_f3);
fclose(output_f3);
file_input2(datfile2);
if(g_opt==1)haplo(intefile);
if(g_opt==2)diplo(intefile);
}
/*////////////////////////////////////////////////////////////////////////////*/
if(g_opt==3)
{
/*CONSTRUCT MAP IN SECTIONS */
strcpy(jobfile,argv[1]);
strcpy(g_jobfile,jobfile);
if ((job_fp=fopen(jobfile, "r"))==NULL){ printf("\nCannot open jobfile.\n\n"); exit(1);}

strcpy(datfile,argv[2]);
strcpy(intefile,argv[3]);
strcpy(outputfile,argv[4]);
  if ((output_f=fopen(outputfile, "w"))==NULL){ printf("\nCannot open outputfile.\n\n"); exit(1);}
segments();
}
/*////////////////////////////////////////////////////////////////////////////*/
/***OPTION 4 ******/
/*Construct LD map ******/
if(g_opt==4) 
  {
  strcpy(g_jobfile,jobfile);
  if ((job_fp=fopen(jobfile, "r"))==NULL){ printf("\nCannot open jobfile.\n\n"); exit(1);}
  if ((output_f=fopen(outputfile, "w"))==NULL){ printf("\nCannot open outputfile.\n\n"); exit(1);}
  printf("\n\n\nL D M A P   - VERSION 1.0 CREATED 12/11/04 \n\n\n"); 
  fprintf(output_f,"\n\n\nL D M A P   - VERSION 1.0 CREATED 12/11/04 \n\n\n"); 
  intfile3();
  locus_list(); 
  ldmap2(); 
  }

/*////////////////////////////////////////////////////////////////////////////*/
/***OPTION 5 ******/

if(g_opt==5)mergeloc();
/*////////////////////////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////////////////////*/
/***OPTION 6 ******/
/*Process pooled count file (hill.out) ******/
if(g_opt==6)
{
printf("\nEnter name of intermediate file to be constructed :");
scanf("%s",intefile);
if ((output_f=fopen(intefile, "w"))==NULL){ printf("\nCannot open intermed file.\n\n"); exit(1);}
  pool(datfile);
}

if(g_opt==8)
{
strcpy(datfile,outputfile);
intfile7();

}
/*////////////////////////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////////////////////*/
if(g_opt==9)
{
/*GENERATE A SUB-SAMPLE DROPPING ONE INDIVIDUAL */
strcpy(ped_file,argv[1]);strcpy(outputfile,argv[2]);
if ((output_f=fopen(outputfile, "w"))==NULL){ printf("\nCannot open outputfile.\n\n"); exit(1);}
file_inputg2();
fclose(output_fi);
}




return 0;
}
