#include "allass.h"

void hapallele()
{
/*SET UP 4 CELL TABLES FOR DIPLOTYPE DATA (ASSOCIATION) -FOR SNPS ONLY */
/*SNPS MUST BE CODED 1,2 */
double diff,kb,kb1,kb2,sum,maxkb;
char temp[20];
int theint,i1,i2,idiff,N,j,numg,i,n11,n12,n21,n22;
pedPtr ped1;
mapPtr2 map1,map2;


maxkb=g_max;
printf("\n\nThe maximum distance between a pair of SNPs to be used is: %10.4f kb, do you wish to change this? (y/n) ",maxkb); 
scanf("%s",temp);
if(temp[0]=='y')
{
printf("\nWhat is the maximum distance to be used ?");
scanf("%s",temp);
maxkb=atof(temp);
g_max=maxkb;
}

theint=g_int;
printf("\nThe maximum number of intervals between pairs of SNPs to be used is: %d, do you wish to change this? (y/n) ",theint); 
scanf("%s",temp);
if(temp[0]=='y')
{
printf("\nWhat is the maximum number of intervals to be used ?");
scanf("%s",temp);
theint=atoi(temp);
g_int=theint;
}
/*******************/
/*******************/
i=0;
map1=map_startPtr2;
while (map1!=NULL)
{
map1->order=i;
i++;
map1=map1->nextPtr;
}
/*******************/
/*******************/


ped1=ped_startPtr;
numg=ped1->numg;
      
/************************************************************************/
/************************************************************************/

fprintf(output_f3,"\n\n\nNumber of loci =%d ",numg+1);
fprintf(output_f3,"\n------------------------------------------------------------------------------------------------------------------------"); 

fprintf(output_f3,"\n(11,1-4) ");
fprintf(output_f3,"\n(12,5-9) ");
fprintf(output_f3,"\n(21,10-14) ");
fprintf(output_f3,"\n(22,15-19) ");
fprintf(output_f3,"\n(N,20-24) ");
fprintf(output_f3,"\n(locus1,26-37) ");
fprintf(output_f3,"\n(locus2,39-50) ");
fprintf(output_f3,"\n(kb1,52-63) ");
fprintf(output_f3,"\n(kb2,65-76) ");
fprintf(output_f3,"\n(freq1,78-90) ");
fprintf(output_f3,"\n(freq2,92-104) ");

fprintf(output_f3,"\n------------------------------------------------------------------------------------------------------------------------"); 
/************************************/
i=0;
map1=map_startPtr2;
while (map1!=NULL)
{
i1=map1->order;
j=i+1;
map2=map1->nextPtr;
    while(map2!=NULL) 
    {
kb1=map1->kb;
kb2=map2->kb;
i2=map2->order;

diff=fabs(kb1-kb2);
if(diff>maxkb)goto skip;

idiff=abs(i1-i2);
if(idiff>theint+1)goto skip;
/************************************************************************/
/*THIS SECTION REFERS TO A PAIR OF LOCI **/
/*FOR THIS PAIR HAVE TO LOOK AT EACH ALLELE COMBINATION */
/************************************************************************/
n11=0;
n12=0;
n21=0;
n22=0;
ped1=ped_startPtr;
while(ped1!=NULL)   
  {
   if(ped1->N==999)ped1->N=1;  
   if(ped1->N==0)ped1->N=1;  
   N=ped1->N; 

   if(ped1->GEN[i][0]!=0&&ped1->GEN[j][0]!=0)
    {
   if((ped1->GEN[i][0]==1))
     {
    if((ped1->GEN[j][0]==1)) n11=n11+N;
    if((ped1->GEN[j][0]==2)) n12=n12+N;
     } 

   if((ped1->GEN[i][0]==2))
     {
    if((ped1->GEN[j][0]==1)) n21=n21+N;
    if((ped1->GEN[j][0]==2)) n22=n22+N;
     } 

     }/*all missing */ 
  ped1=ped1->nextPtr;
  }
/************************************************************************/
/************************************************************************/
sum=n11+n12+n21+n22;
kb1=map1->kb;
kb2=map2->kb;
kb=fabs(kb1-kb2);
map1->locus[13]='\0';
map2->locus[13]='\0';
if(sum>0.)
{
fprintf(output_f3,"\n%4d %4d %4d %4d ",n11,n12,n21,n22);
fprintf(output_f3,"%4d %12s %12s %12.4f %12.4f  %12.10f  %12.10f", n11+n12+n21+n22,map1->locus,map2->locus,kb1,kb2,map1->freq,map2->freq);
} 
skip:map2=map2->nextPtr;
   j++; 
   }/*j index*/
i++;
map1=map1->nextPtr;
}
}
