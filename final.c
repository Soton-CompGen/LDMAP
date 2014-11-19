#include "allass.h"
/*********************************************************************************************/
void final()
/*PRINTING AND TESTS ON FINAL MAP !! */
{
double keepl,keepe,keepm,d,u,ed,e,kb1,kb2,k;
char flank1[20],flank2[20],star;
int jth,ii;
g_meank=0;
rewind(output_ft);
fprintf(output_ft,"# ");
fprintf(output_ft,"\n# Reading the intermediate data file: %s ",datfile);
fprintf(output_ft,"\n# Writing the LD map output file: %s ",intefile);
fprintf(output_ft,"\n# ");
fprintf(output_ft,"\n# FINAL LDU MAP PARAMETERS..............................");
fprintf(output_ft,"\n# E= %10.6f  L= %10.6f  M=%10.6f -2*lnlk= %12.5f",global_e,global_l,global_m,g_lnl);
fprintf(output_ft,"\n# N(number of pairs)=%15d  m(number of SNPs)=%15d df=%14.1f V(error variance)=%14.5f", g_nki,g_nloci,g_df,g_V);

/*************************************************************/
/*************************************************************/
jth=0;
fprintf(output_f,
"\n\nInterval  locus1   location         locus2   location       dj            Ej     LDUs (Ejdj)       UE          KE");
    fflush(output_f);  
for(ii=0;ii<g_nloci;ii++)
{
   if(ii+1<g_nloci)
   { 
   strcpy(flank1,g_loci[ii]);
   strcpy(flank2,g_loci[ii+1]); 
   kb1=g_location[ii][1];
   kb2=g_location[ii+1][1]; 
   jth=jth+1; 
   d=interv[jth][1];
   e=interv[jth][2];
   ed=d*e;
   u=interv[jth][3];
   k=interv[jth][4]; 
   star=' ';
   if(interv[jth][6]==0)star='*';
fprintf(output_f,"\n%3d %12s %12.3f %12s %12.3f %12.6f %10.6f %10.6f %12.3f %15.3f%c ", 
jth,flank1,kb1,flank2,kb2,d,e,ed,u,k,star);  
    fflush(output_f);  
   }
}
keepe=global_e;
keepm=global_m;
keepl=global_l;
writeterfin();
global_e=keepe;
global_m=keepm;
global_l=keepl;
fclose(output_ft);
/*************************************************************/
fflush(output_f);  
fprintf(output_f,"\n\n\n*****    FITTING MALECOT MODEL TO THE FINAL LINKAGE DISEQUILIBRIUM MAP (LDU).................\n ");
ldflag=0;

fprintf(output_f,"\n\nIT(E,M)  \n");
itmast=1;
iteast=1;
itlast=0;
runewt3();

fprintf(output_f,"\n\nIT(E,M,L)  \n");
itmast=1;
iteast=1;
itlast=1;
runewt3();

fflush(output_f);  
}
/*************************************************************/
