#include "allass.h"
/*********************************************************************************************/
void printinputmap()
/*PRINT OUT A TRIAL MAP  */
{
double skb,sed,d,ed,e,kb1,kb2,ldu1,ldu2;
char flank1[20],flank2[20];
int jth,ii;
sed=0.;
skb=0.;
jth=0;
fprintf(output_f,"\n\n ***** THE INPUT MAP INTERVAL DATA ........... ");
fprintf(output_f,"\nInterval  locus1   location         locus2   location       dj            Ej     LDUs (Ejdj) ");
for(ii=0;ii<g_nloci;ii++)
{
   if(ii+1<g_nloci)
   { 
   strcpy(flank1,g_loci[ii]);
   strcpy(flank2,g_loci[ii+1]); 
   ldu1=g_location[ii][0];
   ldu2=g_location[ii+1][0]; 
   kb1=g_location[ii][1];
   kb2=g_location[ii+1][1]; 
   jth=jth+1; 
   d=interv[jth][1];
   e=interv[jth][2];
   ed=d*e;
   sed=sed+ed;
   skb=skb+d;
   fprintf(output_f,"\n%3d %12s %12.3f %12s %12.3f %12.6f %10.6f %10.6f ", jth,flank1,kb1,flank2,kb2,d,e,ed);  
    fflush(output_f);  
   }
}
   fprintf(output_f,"\nTOTAL:                                                  %12.6f            %10.6f ", skb,sed);  
}
