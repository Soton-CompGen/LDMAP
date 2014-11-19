#include "allass.h"
void diplo(char *outputfile)
{
/*Diallelic diplotypes */
/*Needs 9 counts from a 3*3 genotype table */
char xflag=' ';
int numy=0,nrec=0,times,flag,kk;
double nd,kb,kb2,kb1,f11o,D,x11,x12,x21,x22,f11,f12,f21,f22,temp,sum,n11,n12,n13,n21,n22,n23,n31,n32,n33,Q,R;
double a,b,c,d;
output_f2=NULL;
output_f3=NULL;

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
sum=0;
for(kk=0;kk<9;kk++){ if(kk<9)sum=sum+g_aiPtr->aitab[kk]; }
flag=0;
xflag=' ';
if(sum==0)goto skip;
n11=g_aiPtr->aitab[0];
n12=g_aiPtr->aitab[1];
n13=g_aiPtr->aitab[2];
n21=g_aiPtr->aitab[3];
n22=g_aiPtr->aitab[4];
n23=g_aiPtr->aitab[5];
n31=g_aiPtr->aitab[6];
n32=g_aiPtr->aitab[7];
n33=g_aiPtr->aitab[8];
Q=(n11+n12+n13+0.5*(n21+n22+n23))/sum;
R=(n11+n21+n31+0.5*(n12+n22+n32))/sum;

x11=2*n11+n12+n21;
x12=2*n13+n12+n23;
x21=2*n31+n21+n32;
x22=2*n33+n23+n32;

/*STARTING VALUE */
f11=(1./(4.*sum))*(x11-x12-x21+x22)+0.5-(1-Q)*(1-R);
f12=Q-f11;
f21=R-f11;
f22=1-Q-R+f11;
f11o=f11;
/**********************************************/
/**********************************************/
/*ITERATION */
times=0;
top:if(  (f11*(1.-Q-R+f11)+(Q-f11)*(R-f11))>0.)
{
f11=(x11+n22*f11*(1.-Q-R+f11)/(f11*(1.-Q-R+f11)+(Q-f11)*(R-f11))) /(2.*sum);
f12=Q-f11;
f21=R-f11;
f22=1-Q-R+f11;
if(f11<0.)f11=0.;
if(f12<0.)f12=0.;
if(f21<0.)f21=0.;
if(f22<0.)f22=0.;
}

times++;
if(times>14999)xflag='*';
if(times<15000&&fabs(f11o-f11)>0.000000001) { f11o=f11; goto top; }

/**********************************************/

/*CHECK, GET f11 under E.M. algorithm*/
nd=sum;
numy++;
/**********************************************/
a=f11;
b=f12;
c=f21;
d=f22;
Q=a+b;
R=a+c;
if(Q<=0.000001)goto skip;
if(R<=0.000001)goto skip;
if((1.-Q)<=0.000001)goto skip;
if((1.-R)<=0.000001)goto skip;

/**********************************************/
/**********************************************/
/*DO SWITCHING TO GET Q,R*/

/* CONSTRAINTS */
/*We require Q to be the minimum allele frequency*/
if(Q>R)             {temp=b;b=c;c=temp; Q=a+b; R=a+c;}
if(Q>(1.-R))        {temp=a;a=d;d=temp; Q=a+b; R=a+c;}
if(Q>(1.-Q))        {temp=a;a=c;c=temp;temp=b;b=d;d=temp; Q=a+b; R=a+c;}

/*SWITCHES R and 1-R TO ENSURE ad>bc*/
if((a*d)<(b*c)){temp=a;a=b;b=temp;temp=c;c=d;d=temp; Q=a+b; R=a+c;}

f11=a;f12=b;f21=c;f22=d; 
D=f11*f22-f12*f21;

/*D=0 */
if((D<0.000001)&&(b>a)) {temp=a;a=b;b=temp; temp=c;c=d;d=temp; Q=a+b; R=a+c;}

f11=a; f12=b; f21=c; f22=d; 
D=f11*f22-f12*f21;
if(D<0.0)D=0.0;
/*The metrics and their information */
/*sum is number of random diplotypes*/
metricho(Q,R,sum,D);
if(g_rho>1.0)g_rho=1.0;
if(g_rho<0.0)g_rho=0.0;
/***********************************/

if(nrec>0)
{
fprintf(output_f2,"\n%-15s %-15s %10.3f %10.3f %12.10f %10.3f %8.2f %5.0f %12.10f %12.10f %12.10f max_int %6d max_kb %8.2f %12.10f %12.10f", 
g_aiPtr->locus1,g_aiPtr->locus2,kb1,kb2,g_rho,g_rhoi,(g_rho*g_rho)*g_rhoi,sum,Q,R,D,g_int,g_max,g_aiPtr->freq1,g_aiPtr->freq2);
}
if(nrec==0)
{
fprintf(output_f2,"%-15s %-15s %10.3f %10.3f %12.10f %10.3f %8.2f %5.0f %12.10f %12.10f %12.10f max_int %6d max_Kb %8.2f %12.10f %12.10f", 
g_aiPtr->locus1,g_aiPtr->locus2,kb1,kb2,g_rho,g_rhoi,(g_rho*g_rho)*g_rhoi,sum,Q,R,D,g_int,g_max,g_aiPtr->freq1,g_aiPtr->freq2);
}

nrec++;
skip:g_aiPtr=g_aiPtr->nextPtr;
}
  fprintf(output_f2,"\n");
  fflush(output_f2); 
  fclose(output_f2);

return;
}
