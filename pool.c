
#include "allass.h"



void pool(char *datfile)
{
FILE *fp;
/*Local variables*/
char buffer[200];
char locus1[15],locus2[15],kb1[15],kb2[15],ca[15],cb[15],ccc[15],cd[15];
/* char dorh[10]; */
int ic,iccy=0,nrec,i,j;
allsPtr a2=NULL,a1=NULL;
double x,y,temp,a,b,c,d,D,Q,R,f11,f12,f21,f22,bign,kb11,kb22,aa,bb,cc,dd;
/*POOL FILES FROM EACH POPULATION (hill.out) BY COMBINING HAPLOTYPES COUNTS A,B,C,D
AND THEN RE-COMPUTE RHO AND K RHO TO GIVE A STANDARD (POOLED) INTERMEDIATE FILE */

nrec=0;

all_startPtr=NULL;

if ((fp=fopen(datfile, "r"))==NULL){
printf("\nCannot open datafile.\n\n"); exit(1);}

while(fgets(buffer, 199, fp)!=NULL)
{

j=0;
for(i=0;i<15;i++) { locus1[j]=buffer[i];j++; } locus1[14]='\0'; j=0;
for(i=16;i<31;i++) { locus2[j]=buffer[i];j++; } locus2[14]='\0'; j=0;
for(i=32;i<42;i++) { kb1[j]=buffer[i];j++; } kb1[10]='\0'; kb11=atof(kb1); j=0;
for(i=43;i<53;i++) { kb2[j]=buffer[i];j++; } kb2[10]='\0'; kb22=atof(kb2); j=0;
for(i=54;i<64;i++) { ca[j]=buffer[i];j++; } ca[10]='\0'; aa=atof(ca); j=0;
for(i=65;i<75;i++) { cb[j]=buffer[i];j++; } cb[10]='\0'; bb=atof(cb); j=0;
for(i=76;i<86;i++) { ccc[j]=buffer[i];j++; } ccc[10]='\0'; cc=atof(ccc); j=0;
for(i=87;i<97;i++) { cd[j]=buffer[i];j++; } cd[10]='\0'; dd=atof(cd);

                  if(all_startPtr==NULL)
                   {
                     all_startPtr=(allsPtr)malloc(sizeof(alls));
                     all_startPtr->nextPtr=NULL;
                     a1=all_startPtr;
                     strcpy(a1->locus1,locus1); 
                     strcpy(a1->locus2,locus2); 
                     a1->kb1=kb11; 
                     a1->kb2=kb22; 
                     a1->a=aa; 
                     a1->b=bb; 
                     a1->c=cc; 
                     a1->d=dd;
                     a1->flag=0; 
                     a1->nextPtr=NULL; 
                     }
                     else
                     {
                     a1->nextPtr=(allsPtr)malloc(sizeof(alls));
                     a1=a1->nextPtr;
                     strcpy(a1->locus1,locus1); 
                     strcpy(a1->locus2,locus2); 
                     a1->kb1=kb11; 
                     a1->kb2=kb22; 
                     a1->a=aa; 
                     a1->b=bb; 
                     a1->c=cc; 
                     a1->d=dd; 
                     a1->flag=0; 
                     a1->nextPtr=NULL; 
                     }

} /*end of big while*/



ic=0;
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
a1=all_startPtr;

while(a1!=NULL)
{
ic++;
iccy++;
/*GO THROUGH EACH RECORD AND THEN SEARCH REST OF STRUCTURE FOR ALL MATCHING RECORDS - 
ELIMINATE THOSE WHICH HAVE ALREADY BEEN USED (BY FLAG)*/
aa=bb=cc=dd=0;
if(a1->flag==0)
{
a1->flag=1;
strcpy(locus1,a1->locus1);
strcpy(locus2,a1->locus2);
kb11=a1->kb1;
kb22=a1->kb2;
aa=a1->a;
bb=a1->b;
cc=a1->c;
dd=a1->d;
/*NOW SEARCH THE REST OF THE DATA FOR MATCHES AND FORM SUMS OVER COUNTS */
/*______________________________________*/           
/*______________________________________*/           
   a2=a1->nextPtr;
   while(a2!=NULL)
    {
    if(a2->flag!=1)
      { 
      if ( ((strcmp(a2->locus1,locus1)==0)&& (strcmp(a2->locus2,locus2)==0)) ||
      ((strcmp(a2->locus2,locus1)==0)&& (strcmp(a2->locus1,locus2)==0)))
         {
          aa=aa+a2->a;
          bb=bb+a2->b;
          cc=cc+a2->c;
          dd=dd+a2->d;
          iccy++; 
          a2->flag=1; /*IF MATCH - SET FLAG TO 1 TO PREVENT RECORD BEING USED AGAIN*/ 
         }
      }
   a2=a2->nextPtr;
   }
/*______________________________________*/           
/*______________________________________*/           
/*WRITE INT FILE IN HERE*/
/***************************************************************************/
bign=(aa+bb+cc+dd);
f11=aa/bign;
f12=bb/bign;
f21=cc/bign;
f22=dd/bign;
D=f11*f22-f12*f21;
a=f11;
b=f12;
c=f21;
d=f22;
Q=a+b;
R=a+c;
/*DO SWITCHING TO GET Q,R*/
/*3 CONSTRAINTS */
if(Q>(1.-Q))        {temp=a;a=c;c=temp;temp=b;b=d;d=temp; Q=a+b; R=a+c;}
if(Q>R)            {temp=b;b=c;c=temp; Q=a+b; R=a+c;}
if(Q>(1.-R))        {temp=a;a=d;d=temp; Q=a+b; R=a+c;}
/*SWITCHES R and 1-R TO ENSURE ad>bc*/
if((a*d)<(b*c)){temp=a;a=b;b=temp;temp=c;c=d;d=temp; Q=a+b; R=a+c;}
f11=a;f12=b;f21=c;f22=d;
D=f11*f22-f12*f21;
/*D=0 */
if(D<0.000001) { if((b>a)||(b>c)||(b>d)) { temp=a;a=b;b=temp; temp=c;c=d;d=temp; Q=a+b; R=a+c; } }
f11=a;f12=b;f21=c;f22=d;
D=f11*f22-f12*f21;
if(Q<=0.000001)goto skip;
/*Need n for random diplotypes or random haplotypes*/
metricho(Q,R,bign,D);
if(nrec>0)
{
fprintf(output_f,"\n%-15s %-15s %10.3f %10.3f %12.10f %10.3f %8.2f %5.0f %12.10f %12.10f %12.10f ", 
locus1,locus2,kb11,kb22,g_rho,g_rhoi,(g_rho*g_rho)*g_rhoi,bign,Q,R,D);
fflush(output_f);
}
if(nrec==0)
{
fprintf(output_f,"%-15s %-15s %10.3f %10.3f %12.10f %10.3f %8.2f %5.0f %12.10f %12.10f %12.10f ", 
locus1,locus2,kb11,kb22,g_rho,g_rhoi,(g_rho*g_rho)*g_rhoi,bign,Q,R,D);
fflush(output_f);
}
nrec++;
/*WRITE INT FILE ABOVE HERE*/
skip:x=y;
} /*END a1->flag=0*/

printf("\n....Processing file record # %d...clusters so far: %d ",iccy,ic);
a1=a1->nextPtr;
} /*end of big while*/

/****************************************************************************/
/*MAKE hill.sum file */
}
