#include "allass.h"
/*************************************************************/
void updat1()
{
double d,sed,ed,e;
int jth,ii;
            
sed=0.;
jth=0;

for(ii=0;ii<g_nloci;ii++)
{
   if(ii+1<g_nloci)
   { 
if(jth==0){ g_location[ii][0]=sed; }

   jth=jth+1; 
   d=interv[jth][1];
   e=interv[jth][2];
   ed=d*e;
   sed=sed+ed;
if(jth>0){ g_location[ii+1][0]=sed; }
   }
}
}
