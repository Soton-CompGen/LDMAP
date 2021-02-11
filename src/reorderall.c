
#include "allass.h"
/*********************************************************************************************/
/*********************************************************************************************/
void reorderall()
{
int largi=99999999;
double maxkb=0,minkb,larg=99999999.,large=99999999999., xmin,location[MAX_LOCI][3];
char loci[MAX_LOCI][30];
int ikeep, ic,i;

minkb=large;

/****************************************************************************************************/
/*This section to re-order the map and then define the interval flankers being studied */
/*Note that the kb map is stored in position 1. Since locations are unique, the map is best
ordered on these locations as below*/
ic=0;

top: xmin=larg;
ikeep=largi;

for(i=0;i<g_nallloci;i++)
  {
    if(g_alllocation[i][1]<=xmin&&g_alllocation[i][1]<larg)
     {
       ikeep=i;
       xmin=g_alllocation[i][1];
       if(xmin>maxkb&&xmin!=larg)maxkb=xmin;
       if(xmin<minkb&&xmin!=larg)minkb=xmin;
     }
   }

if(ikeep!=largi)
  {
    strcpy(loci[ic],g_allloci[ikeep]);
    location[ic][0]=g_alllocation[ikeep][0]; 
    location[ic][1]=g_alllocation[ikeep][1]; 
    location[ic][2]=g_alllocation[ikeep][2]; 
    g_alllocation[ikeep][1]=larg; 
    ic++; 
    goto top;
  }
for(i=0;i<g_nallloci;i++)
{
g_alllocation[i][0]=location[i][0];
g_alllocation[i][1]=location[i][1];
g_alllocation[i][2]=location[i][2];
strcpy(g_allloci[i],loci[i]);
}
}
