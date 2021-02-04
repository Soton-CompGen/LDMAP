#include "allass.h"

void chi(double df, double chi)
{
double thep,prob,chi1;
double xp,c0,c1,c2,d1,d2,d3,t,t2,t3;

t=gammp(df,chi);
prob=1.-global_p;
c0 = 2.515517;
c1 = 0.802853;
c2 = 0.010328;
d1 = 1.432788;
d2 = 0.189269;
d3 = 0.001308;
thep = prob/2;
if(g_flag==1)t=sqrt(log(1./(thep*thep)));
t2 = t*t;
t3 = t2*t;

xp = t - (c0+c1*t+c2*t2)/ (1.+d1*t+d2*t2+d3*t3) ;
chi1 = xp*xp; 

g_chi1=chi1;
g_prob=2.*thep;
}

/*-----------------------------------------------------------------------------------
Function name:		gammp
Date last modified:	7/2/95
Description:		Chi square function
Variables passed:	a and x
Global variables used:	-
Variables returned:	result
------------------------------------------------------------------------------------*/

double gammp(double a , double x)

{

void gcf(double *gammcf,double a, double x,double *gln,double *bigp);
void gser(double *gamser,double a, double x, double *gln, double *bigp);
double gamser,gammcf,gln,bigp;
a=a/2.;
x=x/2.;

if(x<0.0 || a<=0.0){printf("\nInvalid arguments in gammp");}
if(x<(a+1.0)){
gser(&gamser,a,x,&gln,&bigp);
global_p = bigp;
g_flag=1;
return gamser;

} else { 
   gcf(&gammcf,a,x,&gln,&bigp);
    global_p = 1.-bigp;  
    g_flag=0; 
    return gammcf;

}
}

#define ITMAX 50
#define EPS 3.0e-7
void gser(double *gamser,double a,double x,double *gln,double *bigp)
{
double gammln(double xx);
int n;
double sum,del,ap;
*gln=gammln(a);
if(x<=0.0){
if(x<0.0){printf("x less than 0 in routine gser");return;}

*gamser=0.0;
return;
}else{
ap=a;
del=sum=1.0/a;
for(n=1;n<=ITMAX;n++){
++ap;
del *= x/ap;
sum += del;
if(fabs(del)<fabs(sum)*EPS){
*bigp=sum*exp(-x+a*log(x)-*gln);
return;
}}
printf("a too large, ITMAX too small in routine gser");
return;
}
}

#define ITMAX 50
#define EPS 3.0e-7
#define FPMIN   1.0e-300
void gcf(double *gammcf,double a,double x,double *gln,double *bigp)
{
  double gammln(double xx);
  int i;
  double an,b,c,d,del,h;
  *gln=gammln(a);
  b=x+1.0-a;
  c=1.0/FPMIN;
  d=1.0/b;
h=d;
 for(i=1;i<=ITMAX;i++){
  an=-i*(i-a);
  b+=2.0;
  d=an*d+b;
 if(fabs(d)<FPMIN)d=FPMIN;
 c=b+an/c;
 if(fabs(c)<FPMIN)c=FPMIN;
 d=1.0/d;
 del=d*c;
 h*=del;
 if(fabs(del-1.0)<EPS) break;
}
if(i>ITMAX){printf("a too large, ITMAX too small in gcf");}
*bigp=exp(-x+a*log(x)-(*gln))*h;
*gammcf =  sqrt((-2.* (log(h)-x + a*log(x)- *gln  ))+log(4.));
}

double gammln(double xx)
{
 double x,y,tmp,ser;
static double cof[6]={76.18009172947146,-86.50532032941677,
24.01409824083091,-1.231739572450155,
0.1208650973866179e-2,-0.5395239384953e-5};
int j;
y=x=xx;
tmp=x+5.5;
tmp -=(x+0.5)*log(tmp);
ser=1.000000000190015;
for(j=0;j<=5;j++)ser+=cof[j]/++y;
return -tmp+log(2.5066282746310005*ser/x);
}
