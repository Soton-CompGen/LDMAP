#include "allass.h"

void multallele1()
{
/*SET UP 9 CELL TABLES FOR DIPLOTYPE DATA (ASSOCIATION) -FOR ONE SNP ONLY */
/*SNPS MUST BE CODED 1,2 */
double diff,kb,kb1,kb2,sum,maxkb;
char temp[20];
int idiff,theint,N,j,numg,i,n11,n12,n13;
int n21,n22,n23,i1,i2;
int n31,n32,n33;
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

fprintf(output_f3,"\n(1111,1-4) ");
fprintf(output_f3,"\n(1112,5-9) ");
fprintf(output_f3,"\n(1122,10-14) ");
fprintf(output_f3,"\n(1211,15-19) ");
fprintf(output_f3,"\n(1212,20-24) ");
fprintf(output_f3,"\n(1222,25-29) ");
fprintf(output_f3,"\n(2211,30-34) ");
fprintf(output_f3,"\n(2221,35-39) ");
fprintf(output_f3,"\n(2222,40-44) ");
fprintf(output_f3,"\n(N,45-49) ");
fprintf(output_f3,"\n(locus1,51-65) ");
fprintf(output_f3,"\n(locus2,67-82) ");
fprintf(output_f3,"\n(kb1,83-97) ");
fprintf(output_f3,"\n(kb2,98-110) ");
fprintf(output_f3,"\n------------------------------------------------------------------------------------------------------------------------"); 
/************************************/
i=0;
map1=map_startPtr2;
while (map1!=NULL)
{
i1=map1->order;
j=0;
if(strcmp(map1->locus,g_SNP)==0)
{
    map2=map_startPtr2;
    while(map2!=NULL) 
    {
    if(strcmp(map2->locus,g_SNP)==0)goto skip;
kb1=map1->kb;
kb2=map2->kb;
i2=map2->order;
diff=fabs(kb1-kb2);
if(diff>maxkb)goto skip;
idiff=abs(i1-i2);
if(idiff>theint+1)goto skip;

n11=0; n12=0; n13=0; n21=0; n22=0; n23=0; n31=0; n32=0; n33=0;


ped1=ped_startPtr;
while(ped1!=NULL)   
  {

   if(ped1->N==999)ped1->N=1;  
   if(ped1->N==0)ped1->N=1;  
   N=ped1->N; 
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
kb=fabs(kb1-kb2);
map1->locus[13]='\0';
map2->locus[13]='\0';
if(sum>0.)
{
fprintf(output_f3,"\n%4d %4d %4d %4d %4d %4d %4d %4d %4d ",n11,n12,n13,n21,n22,n23,n31,n32,n33);
fprintf(output_f3,"%4d  %12s     %12s     %12.4f %12.4f ", 
n11+n12+n13+n21+n22+n23+ n31+n32+n33,map1->locus,map2->locus,kb1,kb2);
} 
skip:map2=map2->nextPtr;
j++; 
}/*j index*/

} /*Test on selected SNP found */ 
i++;
map1=map1->nextPtr;
}
}
