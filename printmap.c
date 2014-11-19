#include "allass.h"
/*********************************************************************************************/
void printmap()
/*PRINTING A MAP !! */
{
double loc1;
char flank1[20];
int ii;

/*************************************************************/
fprintf(output_f,"\n\nMAP USED >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
fprintf(output_f,"\n      Locus     Location     \n    ");

for(ii=0;ii<g_nloci;ii++)
{
   strcpy(flank1,g_loci[ii]);
   loc1=g_location[ii][0];
  
fprintf(output_f,"\n%12s %12.6f ",flank1,loc1); 
}
/*************************************************************/
}
