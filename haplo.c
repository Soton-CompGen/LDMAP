#include "allass.h"
void haplo(char *outputfile)
{
/*Haplotype */
/*Needs 4 counts from a 4*4 genotype table */
char xflag=' ';
int numy=0,usemean,nrec=0,flag,kk;
double kb,kb2,kb1,D,temp,nh,n11,n12,n21,n22,Q,R;
double a,b,c,d;
output_f2=NULL;

usemean=0;
if((output_f2=fopen(outputfile,"w"))==NULL)
   {printf("\nCannot open intermediate file\n");exit(1);} 

/*hill.out file gives the haplotype counts A,B,C,D - can be used for pooling across
populations*/
g_aiPtr=aistartPtr;

while(g_aiPtr!=NULL)
{
kb1=atof(g_aiPtr->ckb1);
kb2=atof(g_aiPtr->ckb2);
kb=fabs(kb1-kb2);
nh=0;
for(kk=0;kk<4;kk++){ if(kk<4)nh=nh+g_aiPtr->aitab[kk]; }
flag=0;
xflag=' ';
if(nh==0)goto skip;
n11=g_aiPtr->aitab[0];
n12=g_aiPtr->aitab[1];
n21=g_aiPtr->aitab[2];
n22=g_aiPtr->aitab[3];
a=n11/nh;b=n12/nh;c=n21/nh;d=n22/nh;

numy++;

D=a*d-b*c;
Q=a+b;
R=a+c;
/**********************************************/
/*DO SWITCHING TO GET Q,R*/
/*3 CONSTRAINTS */
if(n11*n22<n12*n21){temp=a;a=c;c=temp;temp=b;b=d;d=temp;
                    Q=a+b; R=a+c;}

if(Q>R)            {temp=b;b=c;c=temp;
                    Q=a+b; R=a+c;}
                
if(Q>(1-R))        {temp=a;a=d;d=temp;
                    Q=a+b; R=a+c;}

if(fabs(D)<0.000001)
         {
           if((b>a)||(b>c)||(b>d))
            {
            temp=a;a=b;b=temp;
            temp=c;c=d;d=temp;
                    Q=a+b; R=a+c;
             }
         }
/*D=0 */
D=fabs(D);
/*SKIP RECORDS HERE !!!! */
if(Q<=0.000001)goto skip;
/*The metrics and their information */
/*nh is number of random haplotypes*/
metricho(Q,R,nh,D);
/***********************************/
if(nrec>0)
{
fprintf(output_f2,"\n%-15s %-15s %10.3f %10.3f %12.10f %10.3f %8.2f %5.0f %12.10f %12.10f %12.10f max_int %6d max_kb %8.2f %12.10f %12.10f", 
g_aiPtr->locus1,g_aiPtr->locus2,kb1,kb2,g_rho,g_rhoi,(g_rho*g_rho)*g_rhoi,nh,Q,R,D,g_int,g_max,g_aiPtr->freq1,g_aiPtr->freq2);
}
if(nrec==0)
{
fprintf(output_f2,"%-15s %-15s %10.3f %10.3f %12.10f %10.3f %8.2f %5.0f %12.10f %12.10f %12.10f max_int %6d max_Kb %8.2f %12.10f %12.10f", 
g_aiPtr->locus1,g_aiPtr->locus2,kb1,kb2,g_rho,g_rhoi,(g_rho*g_rho)*g_rhoi,nh,Q,R,D,g_int,g_max,g_aiPtr->freq1,g_aiPtr->freq2);
}

nrec++;
skip:g_aiPtr=g_aiPtr->nextPtr;
}
  fprintf(output_f2,"\n");
  fflush(output_f2); 
  fclose(output_f2);
return;
}
