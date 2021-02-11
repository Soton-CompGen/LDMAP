
#include "allass.h" 
/********************************************************************/
	void get_predl()
{
abcdPtr abcd_p1;
double pi,ki;
g_sumki=0.;
g_sumki2=0.;
g_nki=0;
/* Go through each marker in turn */
abcd_p1=abcdstartPtr;
while(abcd_p1!=NULL)
{
pi=abcd_p1->p; 
ki=abcd_p1->k;
g_sumki=g_sumki+ki;
g_nki++;
g_sumki2=g_sumki2+sqrt(ki);
abcd_p1=abcd_p1->nextPtr;
}
g_pred=sqrt(2./3.14159265)*(g_sumki2/g_nki)/(g_sumki/g_nki); 
}
