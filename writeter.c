#include "allass.h"
/*************************************************************/
void writeter()
{
double sumdf=0.,sumln=0.,freq1,freq2,dplus,d,sed,ed,e,kb1,kb2;
char flank1[20],copy1[20],flank2[20],copy2[20];
int jth,ii;

/*************************************************************/
/*************************************************************/
dplus=0;
sed=0.;
jth=0;
if(path==1||path==3)
{

fprintf(output_ft,"\n# ");
fprintf(output_ft,"\n#           Locus     kb map       LDU map      MAF                                           ");
fprintf(output_ft,"\n# ");
}

for(ii=0;ii<g_nloci;ii++)
{
   if(ii+1<g_nloci)
   { 
   strcpy(flank1,g_loci[ii]);
   strcpy(flank2,g_loci[ii+1]); 
   kb1=g_location[ii][1];
   kb2=g_location[ii+1][1]; 
   freq1=g_location[ii][2];
   freq2=g_location[ii+1][2];

if(jth==0){
          strcpy(copy1,flank1); 
          meankchi(copy1); 
          sumdf=sumdf+ g_meank;
          sumln=sumln+g_meanchi; 

          fprintf(output_ft,"\n%5d %12s %12.6f %12.6f  %10.8f                                          ",
jth+1,flank1,kb1,sed,freq1);
          
          }

   jth=jth+1; 
   d=interv[jth][1];
   dplus=dplus+d; 
   e=interv[jth][2];
   ed=d*e;
   sed=sed+ed;
if(jth>0){
         strcpy(copy2,flank2);
         meankchi(copy2); 
         sumdf=sumdf+ g_meank;
         sumln=sumln+g_meanchi; 
          
         fprintf(output_ft,"\n%5d %12s %12.6f %12.6f  %10.8f                                            ",
jth+1,flank2,kb2,sed,freq2);
         }
   }
fflush(output_ft);
}

fprintf(output_ft,"\n#                                                                                          ");

}
