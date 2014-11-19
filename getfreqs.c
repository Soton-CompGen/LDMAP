#include "allass.h"

void getfreqs()
{
int r,nr,N,
numg,nals[MAX_LOCI],i,altab[2][2][MAX_LOCI];
pedPtr ped1;
mapPtr2 p1;
double fmin,f;

/*numg is the number of loci in the file */
/*nals is the rth allele at locus i */
/*altab - r alleles, 0=allele label, 1 = allele count ith locus */

for(i=0;i<MAX_LOCI;i++)
{
 nals[i]=0;
 altab[0][0][i]=0;
 altab[0][1][i]=0;
 altab[1][0][i]=0;
 altab[1][1][i]=0;
}

ped1=ped_startPtr;
numg=ped1->numg;
/******************************************************************************/
for(i=0;i<=numg;i++)
{
ped1=ped_startPtr;
while(ped1!=NULL)   
  {
   if(ped1->N==999)ped1->N=1;  

      for(r=0;r<2;r++)
          {
          if(ped1->GEN[i][0]!=0)    
             {if(altab[r][0][i]==ped1->GEN[i][0]) {altab[r][1][i]=altab[r][1][i]+ped1->N;goto lab1;}}
          } 
           if(ped1->GEN[i][0]!=0)    
             {altab[nals[i]][0][i]=ped1->GEN[i][0];
              altab[nals[i]][1][i]=ped1->N; 
              nals[i]++;}
lab1:   for(r=0;r<2;r++)
          {
           if(ped1->GEN[i][1]!=0)    
             {if(altab[r][0][i]==ped1->GEN[i][1]) {altab[r][1][i]=altab[r][1][i]+ped1->N;goto lab2;}}
          } 
           if(ped1->GEN[i][1]!=0)    
             {altab[nals[i]][0][i]=ped1->GEN[i][1]; altab[nals[i]][1][i]=ped1->N; nals[i]++; }
lab2: ped1=ped1->nextPtr;
  }
}
/*****************************************************************************/
/*LOOP OVER ALL LOCI */
p1=map_startPtr2;
for(i=0;i<=numg;i++)
  {
    fmin=1.;    
    N=0;
    for(r=0;r<2;r++)N=N+altab[r][1][i];
for(r=0;r<2;r++)
     {
       nr=altab[r][1][i];
       f=(double)nr/(double)N;
       if (f<fmin)fmin=f; 

     }
p1->freq=fmin;
p1=p1->nextPtr;
}
  }
