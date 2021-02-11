
#include "allass.h"
/*********************************************************************************************/
void update()
{
int jth,ii;
double ed,d,sed,e;
char flank2[20];
/*This function updates M and E for the LD map that is being constructed*/
/*Replace locations with current LDU ones in g_abcd structure prior to re-fitting pairwise data*/
/*Note that interv has to hold the current LDU map data */


   g_abcd=abcdstartPtr;
   while(g_abcd!=NULL)
    {
      g_abcd->ldu1=0.;
      g_abcd->ldu2=0.;
      g_abcd=g_abcd->nextPtr;
    }

sed=0.;
jth=0;
g_nhole=0;
for(ii=0;ii<g_nloci;ii++)
{
   if(ii+1<g_nloci)
   { 
   strcpy(flank2,g_loci[ii+1]); 
   jth=jth+1; 
   d=interv[jth][1];
   e=interv[jth][2];
   ed=d*e;
   sed=sed+ed;
   if(interv[jth][6]==0)g_nhole++;
   /*THIS SETS THE MAP DISTANCES TO THE CURRENT ONES*/
   g_abcd=abcdstartPtr;
   while(g_abcd!=NULL)
    {
      if(strcmp(flank2,g_abcd->locus)==0) { g_abcd->ldu1=sed; }
      if(strcmp(flank2,g_abcd->locus2)==0) {g_abcd->ldu2=sed; }
      g_abcd=g_abcd->nextPtr;
    }
   }
}
/*NOTE THAT THE DIFFERENCE IS USED BY DFUNC*/
g_abcd=abcdstartPtr;
while(g_abcd!=NULL)
{
g_abcd->ldu=fabs(g_abcd->ldu1-g_abcd->ldu2);       
g_abcd=g_abcd->nextPtr;
}

rewind(job_fp);
/*FIT pairwise data to the current LDU map - estimating E, (L) and M */
ldflag=2;
iteast=1;
itmast=1;
itlast=1;
if(g_estimatel==0){itlast=0;global_l=g_pred;}

runewt3();
return;
}
