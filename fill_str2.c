
#include "allass.h"
         void fill_str2(char locus[20],char ldu[20])
         {

         /****FILL STRUCTURE !!!! *****************/
         /****FILL STRUCTURE !!!! *****************/
                  if(ldsstartPtr==NULL)
                   {
                     ldsstartPtr=(ldsPtr)malloc(sizeof(lds));
                     ldsstartPtr->nextPtr=NULL;
                     ldPtr=ldsstartPtr;
                     strcpy(ldPtr->locus,locus); 
                     strcpy(ldPtr->ldu,ldu); 
                     ldPtr->nextPtr=NULL; 
                     }
                     else
                     {
                     ldPtr->nextPtr=(ldsPtr)malloc(sizeof(lds));
                     ldPtr=ldPtr->nextPtr;
                     strcpy(ldPtr->locus,locus); 
                     strcpy(ldPtr->ldu,ldu); 
                     ldPtr->nextPtr=NULL; 
                     }
           /****FILL STRUCTURE !!!! *****************/
           /****FILL STRUCTURE !!!! *****************/

           }
