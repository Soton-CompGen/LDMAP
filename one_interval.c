
#include "allass.h"
void one_interval(int ii, intsPtr intop,double *u, double *bigk, double *lnl )
{
/*INPUT Malecot M and L, INPUT AND OUTPUT u,bigk,lnl*/
/*GLOBAL - NEEDS TO BE PASSED = interv ARRAY + global_l, global_m, g_nloci*/

double u1,bigk1,lnl1=0.,term0,term,l,m,ed,dd,pexp,k,p;

u1=0.;
bigk1=0.;


while(intop!=NULL)
    {
     
         dd=intop->dd;
         ed=intop->ed;
         l=global_l;
         m=global_m;
         k=intop->k;
         p=intop->ab2p;
         pexp=(1.-l)*m*exp(-ed)+l; 
         term0=k*(p-pexp); 
         term=-(1.-l)*m*dd*exp(-ed);
         /*u, k and lnl are returned */
         u1=u1+(term0*term); 
         bigk1=bigk1+(k*term*term); 
         lnl1=lnl1-(k*((p-pexp)*(p-pexp)))/2.;
       intop=intop->nextPtr;
       }
*u=u1;
*bigk=bigk1;
*lnl=*lnl+lnl1;
}
