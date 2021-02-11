#include "allass.h"

void newentry()
{
int i;
/******************************************************************************************/
   if(g_ped1==NULL)
     {
      ped_startPtr=(pedPtr)malloc(sizeof(peds));
      ped_startPtr->nextPtr=NULL;
      g_ped1=ped_startPtr;
      for(i=0;i<MAX_LOCI;i++){g_ped1->GEN[i][0]=0;g_ped1->GEN[i][1]=0;} 
      } 
      else
      {
      g_ped1->nextPtr=(pedPtr)malloc(sizeof(peds));
      g_ped1=g_ped1->nextPtr;
      for(i=0;i<MAX_LOCI;i++){g_ped1->GEN[i][0]=0;g_ped1->GEN[i][1]=0;} 
     }
}
