
#include "allass.h"
/*************************/
void runewt2()
{
double gtol,x[5],fret;
double *fvec;
int keepflag,iter,i,j,n;
gtol=0.0000001;
if(ldflag!=1)
{
g_sig=0;
lab: fprintf(output_f,"\n--------------------------------------------------------------------------------------------------------------------------------");
jobin();

/*NEED TO RE-THINK THIS TO AVOID INTERACTION WITH END OF FILE CC  MARK ! */
if(g_sig==0)
{
fprintf(output_f,
"\nIT   -2lnlk          E             L             M                          Kee            Kll            Kmm                        ");  
fflush(output_f);
g_lnl = 999999.;
g_kel=0.; g_kem=0.; g_klm=0.; 
n=0;
global_fin=0;
g_iter[1]=0;
keepflag=0;
if(ite==1){g_iter[1]=1;n++;}

g_iter[2]=0;
if(itl==1){g_iter[2]=1;n++;}

g_iter[3]=0;
if(itm==1){g_iter[3]=1;n++;}
if(ite==1&&itm==1&&itl==0)keepflag=1;

for(i=0;i<5;i++)x[i]=0.;
j=1;
n=0;
if(g_iter[1]==1){x[j]=global_e;j++;n++;}
if(g_iter[2]==1){x[j]=global_l;j++;n++;}
if(g_iter[3]==1){x[j]=global_m;j++;n++;}

g_bad=0;
g_n=n;
fvec=dvector(1,g_n); dfunc(x,fvec);
/*THIS IS VARIABLE METRIC - 10.7 in Numerical Recipes - note that gradient 'g' is NEGATIVE first derivs. */

dfpmin(x,n,gtol,&iter,&fret,&func,&dfunc);

fflush(output_f);
}

if(g_sig==0)goto lab;

}/*ldflag!=1 */


/**********************************************************************/
if(ldflag==1)
{
ite=1;itl=1;itm=1;

fprintf(output_f,"\n--------------------------------------------------------------------------------------------------------------------------------");
fprintf(output_f,
"\nIT   -2lnlk          E             L             M                          Kee            Kll            Kmm                         ");  
n=0;
global_fin=0;
g_iter[1]=0;
if(ite==1){g_iter[1]=1;n++;}

g_iter[2]=0;
if(itl==1){g_iter[2]=1;n++;}

g_iter[3]=0;
if(itm==1){g_iter[3]=1;n++;}

for(i=0;i<5;i++)x[i]=0.;
j=1;
n=0;
if(g_iter[1]==1){x[j]=global_e;j++;n++;}
if(g_iter[2]==1){x[j]=global_l;j++;n++;}
if(g_iter[3]==1){x[j]=global_m;j++;n++;}

g_bad=0;
g_n=n;
fvec=dvector(1,g_n); dfunc(x,fvec);
fflush(output_f);
dfpmin(x,n,gtol,&iter,&fret,&func,&dfunc);
}
return;
}
