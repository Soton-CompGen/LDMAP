
#include "allass.h"
/*********************************************************************************************/
/*********************************************************************************************/
void directmap2()
{
double kb,ldu;
char locus[20];
int i;
ldsstartPtr=NULL;
ldPtr=NULL;

for(i=0;i<g_nloci;i++)
{

   strcpy(locus,g_loci[i]);
   kb=g_location[i][1];
   ldu=g_location[i][0];

   /*SAVE TO STRUCTURE*/ 
         /****FILL STRUCTURE !!!! *****************/
         /****FILL STRUCTURE !!!! *****************/
                  if(ldsstartPtr==NULL)
                   {
                     ldsstartPtr=(ldsPtr)malloc(sizeof(lds));
                     ldsstartPtr->nextPtr=NULL;
                     ldPtr=ldsstartPtr;
                     strcpy(ldPtr->locus,locus); 
                     ldPtr->xkb=kb; 
                     ldPtr->xldu=ldu; 
                     ldPtr->nextPtr=NULL; 
                     }
                     else
                     {
                     ldPtr->nextPtr=(ldsPtr)malloc(sizeof(lds));
                     ldPtr=ldPtr->nextPtr;
                     strcpy(ldPtr->locus,locus); 
                     ldPtr->xkb=kb; 
                     ldPtr->xldu=ldu; 
                    ldPtr->nextPtr=NULL; 
                     }
           /****FILL STRUCTURE !!!! *****************/
           /****FILL STRUCTURE !!!! *****************/
}
}
