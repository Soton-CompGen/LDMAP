#include "allass.h"
/*************************************************************/
void writeterfinseg()
{
double ldu1,ldu2,kb1,kb2;
char flank1[20],copy1[20],flank2[20],copy2[20];
int jth,ii;

/* removed */
if((output_ft=fopen(terfile,"w"))==NULL){ printf("\nCannot open map (output) file");exit(1);}
fprintf(output_ft,"# ");
fprintf(output_ft,"\n# LDU MAP PARAMETERS..............................");
fprintf(output_ft,"\n# E= %10.6f  L= %10.6f  M=%10.6f -2*lnlk= %12.5f",global_e,global_l,global_m,g_lnl);
fprintf(output_ft,"\n# N(number of pairs)=%15d  m(number of SNPs)=%15d df=%14.1f V(error variance)=%14.5f", g_nki,g_nloci,g_df,g_V);

/*************************************************************/
/*************************************************************/
jth=0;
fprintf(output_ft,"\n#                                                                           ");
fprintf(output_ft,
"\n#           Locus     kb map       LDU map                                 ");
fprintf(output_ft,"\n# ");

for(ii=0;ii<g_nallloci;ii++)
{
   if(ii+1<g_nallloci)
   { 
   strcpy(flank1,g_allloci[ii]);
   strcpy(flank2,g_allloci[ii+1]); 
   kb1=g_alllocation[ii][1];
   kb2=g_alllocation[ii+1][1]; 
   ldu1=g_alllocation[ii][0];
   ldu2=g_alllocation[ii+1][0]; 

if(jth==0){
          strcpy(copy1,flank1); 
         fprintf(output_ft, "\n%5d %12s %12.5f %12.6f                                                     ",
jth+1,flank1,kb1,ldu1);
         }          
   jth=jth+1; 

if(jth>0){
         strcpy(copy2,flank2);
          fprintf(output_ft,"\n%5d %12s %12.5f %12.6f                                                     ",
jth+1,flank2,kb2,ldu2);
          }
}
}
fflush(output_ft);
fclose(output_ft);
}
