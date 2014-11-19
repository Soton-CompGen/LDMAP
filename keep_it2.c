
#include "allass.h"
/*********************************************************************************************/
void keep_it2()
/*KEEP AN INTERIM MAP !! */
{
double dplus,d,sed,u,ed,e,k;
char flank2[20],star;
int jth,ii;

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
g_sed=sed;
}
