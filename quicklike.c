
#include "allass.h" 
/********************************************************************/
	void quicklike(double E,double L,double M)
{

double lnl,xx,di,phat,pi,ki;
abcdPtr abcd_p1;
lnl=0.;
di=0;
/* Go through each marker in turn */
abcd_p1=abcdstartPtr;
while(abcd_p1!=NULL)
{
pi=abcd_p1->p; 
ki=abcd_p1->k;
if(path==1)di=abcd_p1->kb;
if(path==3)di=abcd_p1->ldu;
xx=-E*di;
phat=(1.-L)*M*exp(xx)+L;
lnl+=-(((pi-phat)*(pi-phat))*ki)/2.;

abcd_p1=abcd_p1->nextPtr;
}
g_lnl=-2.*lnl;
}
