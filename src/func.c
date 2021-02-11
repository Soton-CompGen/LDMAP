#define NRANSI
#include "allass.h" 
/**************************************************************/
double func(double x[])
{
/*Function call by dfpmin - returns only -2lnL */

double smins,xxx,lnl,xx,di,phat,E,L,M,pi,ki;
int n,j;
abcdPtr abcd_p1;
n=0;
j=1;
if(g_iter[1]==1){global_e=x[j];j++;n++;}
if(g_iter[2]==1){global_l=x[j];j++;n++;}
if(g_iter[3]==1){global_m=x[j];j++;n++;}

E=global_e;
L=global_l;
M=global_m;

lnl=0.;

/* Go through each marker in turn */
abcd_p1=abcdstartPtr;
while(abcd_p1!=NULL)
{
if(abcd_p1->flag==1)
{
pi=abcd_p1->p;
ki=abcd_p1->k;
smins=1.0;
di = abcd_p1->kb;

if(path==1||path==2||path==3) { di=abcd_p1->ldu; }

if(g_initial==1) { di = abcd_p1->kb; }

xx=-E*di*smins;
phat=(1.-L)*M*exp(xx)+L;
lnl+=-(((pi-phat)*(pi-phat))*ki)/2.;

}/*TEST ON FLAG - ONLY USE REDUCED SET FOR TESTING SPECIFIC LOCI*/
abcd_p1=abcd_p1->nextPtr;
}
/***********************************/
g_lnl=-2.*lnl;
xxx=g_lnl;
return xxx;
}
