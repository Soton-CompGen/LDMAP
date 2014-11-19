
#include "allass.h"
/*********************************************************************************************/
/*********************************************************************************************/
void reorder()
{
int largi=99999999;
double maxkb=0,minkb,larg=99999999.,large=99999999999., xmin,location[MAX_LOCI][3];
char loci[MAX_LOCI][30];
int ikeep,ic,i;

minkb=large;

/****************************************************************************************************/
/*This section to re-order the map and then define the interval flankers being studied */
/*Note that the kb map is stored in position 1. Since locations are unique, the map is best
ordered on these locations as below*/
ic=0;

top: xmin=larg;
ikeep=largi;

for(i=0;i<g_nloci;i++)
  {
    if(g_location[i][1]<=xmin&&g_location[i][1]<larg)
     {
       ikeep=i;
       xmin=g_location[i][1];
       if(xmin>maxkb&&xmin!=larg)maxkb=xmin;
       if(xmin<minkb&&xmin!=larg)minkb=xmin;
     }
   }

if(ikeep!=largi)
  {
    strcpy(loci[ic],g_loci[ikeep]);
    location[ic][0]=g_location[ikeep][0]; 
    location[ic][1]=g_location[ikeep][1]; 
    g_location[ikeep][1]=larg; 
    ic++; 
    goto top;
  }
for(i=0;i<g_nloci;i++)
{
g_location[i][0]=location[i][0];
g_location[i][1]=location[i][1];
strcpy(g_loci[i],loci[i]);
}

}
