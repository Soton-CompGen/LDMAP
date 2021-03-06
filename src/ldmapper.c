
#include "allass.h"
int main(int argc, char *argv[])
{
char outputfile[50],jobfile[20], temp1[20],temp2[20];
/*This is the maximum number of intervals used:*/
g_int=100;
/*This is the maximum distance between pairs of SNPs used:*/
g_max=500.;
/*This is the Hardy-Weinberg cut off */
g_hwp=0.001;
/*This is the MAF cut off */
g_maf=0.05;

g_finish=0;
g_calls=0;

if(argc != 8)
	{
	fprintf(stderr, "\nIncorrect number of options given \n\nExpected usage: ./ldmapper1 [input.tped] [intermediate.txt] [job] [out.ldmap] [out.log] [MAF] [HWE]\n\n");
	return 1;
	}

strcpy(datfile,argv[1]);strcpy(intefile,argv[2]); strcpy(jobfile,argv[3]);strcpy(terfile,argv[4]);strcpy(outputfile,argv[5]);
strcpy(temp1,argv[6]);strcpy(temp2,argv[7]);

g_maf=atof(temp1);
g_hwp=atof(temp2);

if ((output_f=fopen(outputfile, "w"))==NULL){ printf("\nCannot open outputfile.\n\n"); exit(1);}

/*////////////////////////////////////////////////////////////////////////////*/

/*Make int file ******/
printf("\n**NOTE**: The current settings of the program can accomodate up to %d loci and files up to %d characters wide ",MAX_LOCI,MAX_LINE);
file_inputnew();
getfreqs();
multallele();
diplo(intefile);
/*////////////////////////////////////////////////////////////////////////////*/
/*CONSTRUCT MAP IN SECTIONS */
if ((job_fp=fopen(jobfile, "r"))==NULL){ printf("\nCannot open jobfile.\n\n"); exit(1);}
segments();
/*////////////////////////////////////////////////////////////////////////////*/

return 0;
}
