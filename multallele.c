#include "allass.h"

void multallele()
{
/*SET UP 9 CELL TABLES FOR DIPLOTYPE DATA (ASSOCIATION) -FOR SNPS ONLY */
/*SNPS MUST BE CODED 1,2 */
double diff,kb1,kb2,sum,maxkb;
int theint,i1,i2,idiff,N,j,i,n11,n12,n13;
int n21,n22,n23;
int n31,n32,n33;
pedPtr ped1;
mapPtr2 map1,map2;


maxkb=g_max;
/*
printf("\n\nThe maximum distance between a pair of SNPs to be used is: %10.4f kb, do you wish to change this? (y/n) ",maxkb); 
scanf("%s",temp);
if(temp[0]=='y')
{
printf("\nWhat is the maximum distance to be used ?");
scanf("%s",temp);
maxkb=atof(temp);
g_max=maxkb;
}
*/
theint=g_int;
/*
printf("\nThe maximum number of intervals between pairs of SNPs to be used is: %d, do you wish to change this? (y/n) ",theint); 
scanf("%s",temp);
if(temp[0]=='y')
{
printf("\nWhat is the maximum number of intervals to be used ?");
scanf("%s",temp);
theint=atoi(temp);
g_int=theint;
}
*/
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
n11=0; n12=0; n13=0; n21=0; n22=0; n23=0; n31=0; n32=0; n33=0;
ped1=ped_startPtr;
while(ped1!=NULL)   
  {
   N=1; 
   if(ped1->GEN[i][0]!=0&&ped1->GEN[i][1]!=0&&ped1->GEN[j][0]!=0&&ped1->GEN[j][1]!=0)
    {
    if((ped1->GEN[i][0]==1)&&(ped1->GEN[i][1]==1))
     {
    if((ped1->GEN[j][0]==1)&&(ped1->GEN[j][1]==1)) n11=n11+N;
    if((ped1->GEN[j][0]==1)&&(ped1->GEN[j][1]==2)) n12=n12+N;
    if((ped1->GEN[j][0]==2)&&(ped1->GEN[j][1]==1)) n12=n12+N;
    if((ped1->GEN[j][0]==2)&&(ped1->GEN[j][1]==2)) n13=n13+N;
     } 

   if((ped1->GEN[i][0]==1)&&(ped1->GEN[i][1]==2))
     {
    if((ped1->GEN[j][0]==1)&&(ped1->GEN[j][1]==1)) n21=n21+N;
    if((ped1->GEN[j][0]==1)&&(ped1->GEN[j][1]==2)) n22=n22+N;
    if((ped1->GEN[j][0]==2)&&(ped1->GEN[j][1]==1)) n22=n22+N;
    if((ped1->GEN[j][0]==2)&&(ped1->GEN[j][1]==2)) n23=n23+N;
     } 

   if((ped1->GEN[i][0]==2)&&(ped1->GEN[i][1]==1))
     {
    if((ped1->GEN[j][0]==1)&&(ped1->GEN[j][1]==1)) n21=n21+N;
    if((ped1->GEN[j][0]==1)&&(ped1->GEN[j][1]==2)) n22=n22+N;
    if((ped1->GEN[j][0]==2)&&(ped1->GEN[j][1]==1)) n22=n22+N;
    if((ped1->GEN[j][0]==2)&&(ped1->GEN[j][1]==2)) n23=n23+N;
     } 

   if((ped1->GEN[i][0]==2)&&(ped1->GEN[i][1]==2))
     {
    if((ped1->GEN[j][0]==1)&&(ped1->GEN[j][1]==1)) n31=n31+N;
    if((ped1->GEN[j][0]==1)&&(ped1->GEN[j][1]==2)) n32=n32+N;
    if((ped1->GEN[j][0]==2)&&(ped1->GEN[j][1]==1)) n32=n32+N;
    if((ped1->GEN[j][0]==2)&&(ped1->GEN[j][1]==2)) n33=n33+N;
     } 
     }/*all missing */ 
  ped1=ped1->nextPtr;
  }
/************************************************************************/
/************************************************************************/
sum=n11+n12+n13+n21+n22+n23+n31+n32+n33;

kb1=map1->kb;
kb2=map2->kb;
if(sum>0.)
{
/*Select pairs of loci with sample MAF >0.05 for both*/

if(map1->chi>g_hwp&&map2->chi>g_hwp)
{ 
if(map1->freq>g_maf&&map2->freq>g_maf) 
{
fill_gai(n11,n12,n13,n21,n22,n23,n31,n32,n33,kb1,kb2,map1->locus, map2->locus, map1->freq, map2->freq);
}
}


} 
skip:map2=map2->nextPtr;
   j++; 
   }/*j index*/
i++;
map1=map1->nextPtr;
}

}
