
#include "allass.h"
void fast()
{
/*ESTABLISH ALL OF THE POINTER RECORDS FOR FAST ACCESS TO PAIRWISE DATA*/
/*There are a set of pointers to relevant data for each interval (intsPtr)
and there are a set of pointers to the starting pointers for each of these
(hugsPtr)*/

intsPtr top=NULL, topstart=NULL;
abcdPtr abcd_p1;
hugsPtr big=NULL;
int ii,ikeep,left,right;

bigstart=NULL;
  for(ii=1;ii<g_nloci;ii++)
  { 
  topstart=NULL;
  abcd_p1=abcdstartPtr;
  while(abcd_p1!=NULL)
    {
     left=abcd_p1->i1; 
     right=abcd_p1->i2; 
     if(left>right){ikeep=left;left=right;right=ikeep;} 
     if( (left<ii) && (right>=(ii)) )
        {
         if(topstart==NULL)
            {
              topstart=(intsPtr)malloc(sizeof(ints));
              topstart->nextPtr=NULL;
              top=topstart;
              top->p=abcd_p1;
            }
           else
            {
              top->nextPtr=(intsPtr)malloc(sizeof(ints));
              top=top->nextPtr;
              top->p=abcd_p1;
              top->nextPtr=NULL; 
            }

        }
       abcd_p1=abcd_p1->nextPtr;
       }
         if(bigstart==NULL)
            {
              bigstart=(hugsPtr)malloc(sizeof(hugs));
              bigstart->nextPtr=NULL;
              big=bigstart;
              big->pp=topstart;
            }
           else
            {
              big->nextPtr=(hugsPtr)malloc(sizeof(hugs));
              big=big->nextPtr;
              big->pp=topstart;
              big->nextPtr=NULL; 
            }
     }
}
