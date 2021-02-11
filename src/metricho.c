
#include "allass.h"
void metricho(double Q,double R,double n,double DDD)
{
/*RANDOM SAMPLES !!! */
/*DEFINE DIFFERENT METRICS FOR ASSOCIATION UNDER H0 */
/*n is the number of random haplotypes or diplotypes*/

double D,C,kd;
D=DDD;

kd=n/(Q*(1.-Q)*R*(1.-R));
/*The metrics and their information */

/*Association rho */
C=(Q*(1.-R));
g_rho=D/C;
g_rhoi=(C*C)*kd;


}
