
#include "allass.h"
/*********************************************************************************************/
void keep_it()
/*KEEP AN INTERIM MAP !! */
{
double dplus,d,sed,u,ed,e,k;
char flank2[20],star;
int jth,ii;
rewind(output_ft);

fprintf(output_ft,"# ");
fprintf(output_ft,"\n# USING THE INTERMEDIATE DATA FILE %s ",datfile);
fprintf(output_ft,"\n# WRITING THE LD MAP FILE %s ",intefile);
fprintf(output_ft,"\n# ");
fprintf(output_ft,"\n# INTERIM LDU MAP PARAMETERS..............................");
fprintf(output_ft,"\n# iter=%5.0d E= %10.6f  L= %10.6f  M=%10.6f -2*lnlk= %12.5f",g_niter,global_e,global_l,global_m,g_lnl);
fprintf(output_ft,"\n# N(number of pairs)=%15d  m(number of SNPs)=%15d df=%14.1f V(error variance)=%14.5f", g_nki,g_nloci,g_df,g_V);
/*************************************************************/
g_abcd=abcdstartPtr;
while(g_abcd!=NULL)
{
g_abcd->ldu1=0.;
g_abcd->ldu2=0.;
g_abcd=g_abcd->nextPtr;
}
/*************************************************************/
dplus=0;
sed=0.;
jth=0;
for(ii=0;ii<g_nloci;ii++)
{
   if(ii+1<g_nloci)
   { 
   strcpy(flank2,g_loci[ii+1]); 
   jth=jth+1; 
   d=interv[jth][1];
   dplus=dplus+d; 
   e=interv[jth][2];
   ed=d*e;
   sed=sed+ed;
   u=interv[jth][3];
   k=interv[jth][4]; 
   star=' ';
   if(interv[jth][6]==0)star='*';
   g_abcd=abcdstartPtr;
   while(g_abcd!=NULL)
    {
     if(strcmp(flank2,g_abcd->locus2)==0) { g_abcd->ldu2=sed; }
     if(strcmp(flank2,g_abcd->locus)==0) { g_abcd->ldu1=sed; }
     g_abcd=g_abcd->nextPtr;
    }
   }
}
/*************************************************************/
/*NOTE THAT THE DIFFERENCE IS USED BY DFUNC*/
g_abcd=abcdstartPtr;
while(g_abcd!=NULL)
{
g_abcd->ldu=fabs(g_abcd->ldu1-g_abcd->ldu2);       
g_abcd=g_abcd->nextPtr;
}
/*************************************************************/

writeter();
g_sed=sed;
}
