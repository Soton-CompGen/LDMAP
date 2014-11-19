
#include "allass.h"
/*************************/
void runewt3()
{
double gtol,x[5],fret;
double *fvec;
int iter,i,j,n;
gtol=0.0000001;
if(ldflag==0)
{
fprintf(output_f,"\n--------------------------------------------------------------------------------------------------------------------------------");
fprintf(output_f,"\nIT   -2lnlk          E             L             M             S            Kee            Kll            Kmm            Kss          ");  
fflush(output_f);
}
g_kel=0.; g_kem=0.; g_ket=0.; g_kes=0.; g_klm=0.; g_klt=0.; g_kls=0.; g_kms=0.; g_kts=0.;
n=0;
global_fin=0;

g_iter[1]=0;
if(iteast==1)g_iter[1]=1;n++;

g_iter[2]=0;
if(itlast==1)g_iter[2]=1;n++;

g_iter[3]=0;
if(itmast==1)g_iter[3]=1;n++;

g_iter[4]=0;
for(i=0;i<5;i++)x[i]=0.;
j=1;
n=0;
if(g_iter[1]==1){x[j]=global_e;j++;n++;}
if(g_iter[2]==1){x[j]=global_l;j++;n++;}
if(g_iter[3]==1){x[j]=global_m;j++;n++;}
g_bad=0;
g_n=n;
fvec=dvector(1,g_n); dfunc(x,fvec);
dfpmin(x,n,gtol,&iter,&fret,&func,&dfunc);
if(ldflag==0) fprintf(output_f,"\n--------------------------------------------------------------------------------------------------------------------------------");
return;
}
