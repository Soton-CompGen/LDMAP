#define NRANSI
#include "allass.h" 
/**************************************************************/
double bound(int n,double x[])
{
double xflag;
int j, i;
/*IMPLEMENT BOUNDS - THIS IS DONE AFTER 1 UNCONSTRAINED CYCLE */

xflag=0;
j=1;
if(g_iter[1]==1&&global_e<0.000001){global_e=0.000001;g_iter[1]=0;xflag=1;j++;}
if(g_iter[2]==1&&global_l>1.000000){global_l=1.000000;g_iter[2]=0;xflag=1;j++;}
if(g_iter[2]==1&&global_l<0.000001){global_l=0.000001;g_iter[2]=0;xflag=1;j++;}
if(g_iter[3]==1&&global_m>1.000000){global_m=1.000000;g_iter[3]=0;xflag=1;j++;}
if(g_iter[3]==1&&global_m<0.000001){global_m=0.000001;g_iter[3]=0;xflag=1;j++;}

for(i=0;i<5;i++)x[i]=0.;
j=1;
n=0;
if(g_iter[1]==1){x[j]=global_e;j++;n++;}
if(g_iter[2]==1){x[j]=global_l;j++;n++;}
if(g_iter[3]==1){x[j]=global_m;j++;n++;}
g_n=n;

return xflag;
}

