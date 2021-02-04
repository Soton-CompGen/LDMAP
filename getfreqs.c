#include "allass.h"

void getfreqs()
{
int j,nr1,nr2,r,nr,N,
numg,nals[MAX_LOCI],i,altab[2][2][MAX_LOCI];
pedPtr ped1;
mapPtr2 p1;
double chis[MAX_LOCI],df=1.,a,b,c,xchi,etot,tot,f,n1=0,n2=0,f1,f2,e11,e12,e22,n11,n12,n22;

/*numg is the number of loci in the file */
/*nals is the rth allele at locus i */
/*altab - r alleles, 0=allele label, 1 = allele count ith locus */
/*N is the count of all alleles */
/*nr is the number of rare alleles */
for(i=0;i<MAX_LOCI;i++)
{
 nals[i]=0;
 chis[i]=0.;
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
xchi=tot=etot=0.;
n11=n12=n22=0;
ped1=ped_startPtr;
while(ped1!=NULL)   
  {
   /***Count genotypes******/
           if((ped1->GEN[i][0]!=0)&&(ped1->GEN[i][1]!=0))    
              {
           if((ped1->GEN[i][0]==1)&&(ped1->GEN[i][1]==1))n11=n11+1.;    
           if((ped1->GEN[i][0]==1)&&(ped1->GEN[i][1]==2))n12=n12+1.;    
           if((ped1->GEN[i][0]==2)&&(ped1->GEN[i][1]==1))n12=n12+1.;    
           if((ped1->GEN[i][0]==2)&&(ped1->GEN[i][1]==2))n22=n22+1.;    
              }

   /***Count genotypes******/

      for(r=0;r<2;r++)
          {
          if(ped1->GEN[i][0]!=0)    
             {if(altab[r][0][i]==ped1->GEN[i][0]) {altab[r][1][i]++;goto lab1;}}
          } 
           if(ped1->GEN[i][0]!=0)    
             {altab[nals[i]][0][i]=ped1->GEN[i][0];
              altab[nals[i]][1][i]=1; 
              nals[i]++;}

lab1:   for(r=0;r<2;r++)
          {
           if(ped1->GEN[i][1]!=0)    
             {if(altab[r][0][i]==ped1->GEN[i][1]) {altab[r][1][i]++;goto lab2;}}
          } 
           if(ped1->GEN[i][1]!=0)    
             {altab[nals[i]][0][i]=ped1->GEN[i][1]; 
              altab[nals[i]][1][i]=1; 
              nals[i]++; }

lab2: ped1=ped1->nextPtr;
  }
tot=n11+n12+n22;
if((altab[0][0][i])==1) {n1=altab[0][1][i]; n2=altab[1][1][i]; }
if((altab[1][0][i])==1) {n2=altab[0][1][i]; n1=altab[1][1][i]; }
f1=n1/(n1+n2);
f2=n2/(n1+n2);
e11=(f1*f1)*tot;
e12=2*(f1*f2)*tot;
e22=(f2*f2)*tot;
etot=e11+e12+e22;
if(e11==0)e11=0.5;
if(e12==0)e12=0.5;
if(e22==0)e22=0.5;
a=((n11-e11)*(n11-e11))/e11;
b=((n12-e12)*(n12-e12))/e12;
c=((n22-e22)*(n22-e22))/e22;
xchi=a+b+c;
g_prob=1.;
if(xchi>0.)chi(df,xchi);
/*fprintf(output_f,"\nExp: %10.3f %10.3f %10.3f %10.3f Obs: %10.3f %10.3f %10.3f %10.3f Chi: %10.3f P=%10.6f ", e11,e12,e22,etot,n11,n12,n22,tot,xchi,g_prob);
*/
chis[i]=g_prob;

} /*outer loop over all markers*/
/*****************************************************************************/
fprintf(output_f,"\n\nTotal number of loci in file=%d\n",numg+1);
fprintf(output_f,"\nCut-off used for MAF = %f",g_maf);
fprintf(output_f,"\nCut-off used for Hardy-Weinberg test (Pvalue) = %f ",g_hwp);
fprintf(output_f,"\nLoci excluded with low MAF and/or significant Hardy-Weinberg test ");
fprintf(output_f,"\n                Locus        Kb    N_min_als  N_als     MAF        HW_Pvalue ");
p1=map_startPtr2;
j=0;
for(i=0;i<=numg;i++)
  {
    N=0;
    for(r=0;r<2;r++)N=N+altab[r][1][i];
    nr1=altab[0][1][i];
    nr2=altab[1][1][i];
    nr=nr1;
    if(nr2<nr1)nr=nr2;
    f=(double)nr/(double)N;
    p1->freq=f;
    p1->chi=chis[i];
if((p1->freq<=g_maf)&&(p1->chi>g_hwp)){fprintf(output_f,"\n%20s %10.3f   %5d   %5d     %10.4f ",p1->locus, p1->kb,nr,N,p1->freq);}
if((p1->freq>g_maf)&&(p1->chi<=g_hwp)){fprintf(output_f,"\n%20s %10.3f   %5d   %5d                %10.4f ",p1->locus, p1->kb,nr,N,p1->chi); }
if((p1->freq<=g_maf)&&(p1->chi<=g_hwp)){fprintf(output_f,"\n%20s %10.3f   %5d   %5d     %10.4f %10.4f ",p1->locus, p1->kb,nr,N,p1->freq,p1->chi); }
p1=p1->nextPtr;
}
fflush(output_f);
  }
