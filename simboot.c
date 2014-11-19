#include "allass.h"

void simboot()
{
/*Generate a sub-sample dropping one individual*/
char temp[20];
/* time_t start;
   float xxx,largest;
   int icount=0,n0,n1,ic,j,ikeep,y,icc,ifin,inum,ilocus;
   double x,freq;
*/
 int i; 
printf("\nThe total number of individuals is %d ",g_nrec);

printf("\nEnter the index of the individual to drop from the sample (1 - %4d):",g_nrec);
scanf("\n%s",temp);
g_number=atoi(temp);
/*********************************************************/
for(i=0;i<g_nrec;i++)
{
if(i!=(g_number-1))fprintf(output_f,"\n%s",g_buffer[i]);
fflush(output_f);
}
fclose(output_f);
}
