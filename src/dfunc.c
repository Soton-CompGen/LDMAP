
#include "allass.h" 
/********************************************************************/
	void dfunc(double x[],double fvec[])
{
/*FILL FVEC VECTOR */
/*Called by DFPMIN - RETURNS NEGATIVE FIRST DERIVS = GRADIENT */     

double smins,dldp,lnl,xx,uei,uli,umi,di,phat,E,L,M,pi,ki;
double kee,kel,kem,kes,kll,klm,kmm,kms,kss;
double **fjac,yy,ss2,expx2,negl2,M2,di2;
int j,cnt=1;
int n,row, col;
char blank[14];
abcdPtr abcd_p1;
strcpy(blank,"             ");
j=1;
n=g_n;

g_kel = 0;
g_kem = 0;
g_klm = 0;

g_sumki=0.;
g_sumki2=0.;
g_nki=0;
fjac=dmatrix(1,n,1,n);
n=0;
if(g_iter[1]==1){global_e=x[j]; j++;n++;}
if(g_iter[2]==1){global_l=x[j]; j++;n++;}
if(g_iter[3]==1){global_m=x[j]; j++;n++;}
g_npar=j-1;
E=global_e;
L=global_l;
M=global_m;

for(cnt=1; cnt<=n; cnt++){fvec[cnt]=0.0;}
for(row=1; row<=n; row++) { for(col=1; col<=n; col++) {fjac[row][col]=0.0;} }
lnl=0.;

abcd_p1=abcdstartPtr;
while(abcd_p1!=NULL)
{
if(abcd_p1->flag==1)
{

pi=abcd_p1->p; 
ki=abcd_p1->k;
g_sumki=g_sumki+ki;
g_nki++;
g_sumki2=g_sumki2+sqrt(ki);


smins=1.0;
di=abcd_p1->kb;

if(path==1||path==2||path==3) { di=abcd_p1->ldu; }

/*if LDU map does not yet exist, use kb */
if(g_initial==1) { di=abcd_p1->kb; }

xx=-E*di*smins;
phat=(1.-L)*M*exp(xx)+L;
dldp=(pi-phat)*ki;

xx=-E*di*smins;
M2=M*M;
di2=di*di;
negl2=(1.-L)*(1.-L);
ss2 = smins*smins;

yy=exp(xx);
expx2=yy*yy;
uei = -(pi-(1.-L)*M*exp(xx)-L)*ki*(1.-L)*M*di*smins*exp(xx);

umi = (pi-(1.-L)*M*exp(xx)-L)*ki*(1.-L)*exp(xx); 

uli = -(pi-(1.-L)*M*exp(xx)-L)*ki*(M*exp(xx)-1);



/*Elements of K matrix - for Stdrd Errors & printing only - not used by dfpmin */

kem = ki*negl2*M*di*smins*expx2-(pi-(1.-L)*M*exp(xx)-L)*ki*(1.-L)*di*smins*exp(xx);

kes = ki*negl2*M2*di2*smins*expx2*E+(pi-(1.-L)*M*yy-L)*ki*(1.-L)*M*di*yy-(pi-(1.-L)*M*yy-L)*ki*(1.-L)*M*di2*smins*E*yy; 

kll = -ki*((M*exp(xx)-1.)*(M*exp(xx)-1.));

kss = -ki*negl2*M2*(E*E)*di2*expx2+(pi-(1.-L)*M*exp(xx)-L)*ki*(1.-L)*M*(E*E)*di2*exp(xx); 

kel = -ki*(1.-L)*M*di*smins*yy*(M*yy-1)+(pi-(1.-L)*M*yy-L)*ki*M*di*smins*yy;

kee = -ki*negl2*M2*di2*ss2*expx2+(pi-(1.-L)*M*exp(xx)-L)*ki*(1.-L)*M*di2*ss2*exp(xx);

klm = ki*(M*yy-1.)*(1.-L)*yy-(pi-(1.-L)*M*yy-L)*ki*yy;

kmm =  -ki*negl2*expx2;

kms = -ki*negl2*expx2*M*E*di+(pi-(1.-L)*M*exp(xx)-L)*ki*(1.-L)*E*di*exp(xx);

/*Negative 1st derivs */
j=1;
if(g_iter[1]==1){fvec[j]-=uei;j++;}
if(g_iter[2]==1){fvec[j]-=uli;j++;}
if(g_iter[3]==1){fvec[j]-=umi;j++;}

/*All k matrix - for getting S.Errors & printing only - not used by function */

if(g_iter[1]==1&&g_iter[2]==1&&g_iter[3]==1&&g_iter[4]==0)
{ /*1,2,3*/
fjac[1][1]+=-kee;
fjac[1][2]+=-kel;
fjac[1][3]+=-kem;
fjac[2][1]+=-kel;
fjac[2][2]+=-kll;
fjac[2][3]+=-klm;
fjac[3][1]+=-kem;
fjac[3][2]+=-klm;
fjac[3][3]+=-kmm;

g_kel = fjac[1][2];
g_kem = fjac[1][3];
g_klm = fjac[2][3];
}

if(g_iter[1]==1&&g_iter[2]==1&&g_iter[3]==0&&g_iter[4]==0)
{ /* 1,2 */
fjac[1][1]+=-kee;
fjac[1][2]+=-kel;
fjac[2][1]+=-kel;
fjac[2][2]+=-kll;
g_kel = fjac[1][2];
}

if(g_iter[1]==1&&g_iter[2]==0&&g_iter[3]==1&&g_iter[4]==0)
{ /* 1,3 */
fjac[1][1]+=-kee;
fjac[1][2]+=-kem;
fjac[2][1]+=-kem;
fjac[2][2]+=-kmm;
g_kem = fjac[1][2];
}

if(g_iter[1]==0&&g_iter[2]==1&&g_iter[3]==1&&g_iter[4]==0)
{ /* 2,3 */
fjac[1][1]+=-kll;
fjac[1][2]+=-klm;
fjac[2][1]+=-klm;
fjac[2][2]+=-kmm;
g_klm = fjac[1][2];
}
if(g_iter[1]==1&&g_iter[2]==0&&g_iter[3]==0&&g_iter[4]==0) { fjac[1][1]+=-kee; }
if(g_iter[1]==0&&g_iter[2]==1&&g_iter[3]==0&&g_iter[4]==0) { fjac[1][1]+=-kll; }
if(g_iter[1]==0&&g_iter[2]==0&&g_iter[3]==1&&g_iter[4]==0) { fjac[1][1]+=-kmm; }

lnl+=-(((pi-phat)*(pi-phat))*ki)/2.;
}/*TEST ON FLAG - ONLY USE REDUCED DATA SET FOR TESTING SPECIFIC LOCI */
abcd_p1=abcd_p1->nextPtr;
}
/***********************************/
g_lnl=-2.*lnl;

if(ldflag!=2)
{
fprintf(output_f,"\n %10.2f  ",-2.*lnl);
fprintf(output_f,"%13.5f %13.5f %13.5f               ",E,L,M);
fflush(output_f);

if(g_iter[1]==1&&g_iter[2]==1&&g_iter[3]==1&&g_iter[4]==0)
{
  fprintf(output_f," %13.5f ",fjac[1][1]);
  fprintf(output_f," %13.5f ",fjac[2][2]);
  fprintf(output_f," %13.5f ",fjac[3][3]);
}


if(g_iter[1]==1&&g_iter[2]==1&&g_iter[3]==0&&g_iter[4]==0)
{  
  fprintf(output_f," %13.5f ",fjac[1][1]);
  fprintf(output_f," %13.5f ",fjac[2][2]);
}

if(g_iter[1]==1&&g_iter[2]==0&&g_iter[3]==1&&g_iter[4]==0)
{  
  fprintf(output_f," %13.5f ",fjac[1][1]);
  fprintf(output_f," %13s ",blank);
  fprintf(output_f," %13.5f ",fjac[2][2]);
}

if(g_iter[1]==0&&g_iter[2]==1&&g_iter[3]==1&&g_iter[4]==0)
{  
  fprintf(output_f," %13s ",blank);
  fprintf(output_f," %13.5f ",fjac[1][1]);
  fprintf(output_f," %13.5f ",fjac[2][2]);
}

if(g_iter[1]==1&&g_iter[2]==0&&g_iter[3]==0&&g_iter[4]==0)
 fprintf(output_f," %13.5f ",fjac[1][1]);
if(g_iter[1]==0&&g_iter[2]==1&&g_iter[3]==0&&g_iter[4]==0)
{ 
  fprintf(output_f," %13s ",blank);
fprintf(output_f," %13.5f ",fjac[1][1]);
}

if(g_iter[1]==0&&g_iter[2]==0&&g_iter[3]==1&&g_iter[4]==0)
{
  fprintf(output_f," %13s ",blank);
  fprintf(output_f," %13s ",blank);
 fprintf(output_f," %13.5f ",fjac[1][1]);
}

g_uei=g_uli=g_umi=0.;
j=1;
if(g_iter[1]==1){g_uei=fvec[j];j++;}
if(g_iter[2]==1){g_uli=fvec[j];j++;}
if(g_iter[3]==1){g_umi=fvec[j];j++;}

g_calls++;
if(g_finish==1&&g_calls>0){global_fin=1;}

if(global_fin==1)
{cal_SE2(fjac,n); 
}


fflush(output_f);
}
if(g_finish==1)exit(0);
}
