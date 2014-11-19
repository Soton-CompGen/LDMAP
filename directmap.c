
#include "allass.h"
/*********************************************************************************************/
/*********************************************************************************************/
void directmap()
{
/*Use an LDMAP output terse file to provide LDU locations directly for fitting the data */
double lduold,ldunew,x,bigx;
char buff[200],terse[20];
char toks[20][20],locus[20],ldu[20],kb[20];
FILE *fp2;
int nholes=0,ntoks,k,ic,iset,a,b,j,ii,jj;
ldsstartPtr=NULL;
ldPtr=NULL;
bigx=0.;
lduold=999999999.;

printf("\nTo run this option you will need to use a previously constructed LDMAP map output file");
printf("\nEnter name of map file: ");
scanf("%s",terse);
if((fp2=fopen(terse,"r"))==NULL){printf("\n Cannot open LDMAP map file ");exit(1);}
fprintf(output_f,"\n\n *** The LDU map used is in the file called: %s *** \n\n",terse);
/******************************************************/
/******************************************************/
while(fgets(buff,199,fp2)!=NULL)
{
/*!!!!!!!!!!!!!!!!!!*/
   if(buff[0]!='#') 
    {
    iset=0;
    j=0;
    ntoks=0;
    ic=0;
    for(ii=0;ii<10;ii++){for(jj=0;jj<20;jj++)toks[ii][jj]='\0'; } 
    for(k=0;k<200;k++)
      {
      if((buff[k]==' ')||(buff[k]=='\0')){ic=0;iset=0;goto skip;}
      if(buff[k]!=' ')
        {
         if(iset==1){toks[ntoks][ic]=buff[k];ic++;}
         if(iset==0){iset=1;ntoks++;toks[ntoks][ic]=buff[k];ic++;}
        }
      skip:a=b;
      }

   strcpy(locus,toks[2]);
   strcpy(kb,toks[3]);
   strcpy(ldu,toks[4]);
   x=atof(ldu);
   bigx=bigx+x;
   ldunew=atof(ldu);
   if(lduold<999999.)
   {
   if((ldunew-lduold)>2.5)nholes++;
   }
   lduold=ldunew;




   /*SAVE TO STRUCTURE*/ 
         /****FILL STRUCTURE !!!! *****************/
         /****FILL STRUCTURE !!!! *****************/
                  if(ldsstartPtr==NULL)
                   {
                     ldsstartPtr=(ldsPtr)malloc(sizeof(lds));
                     ldsstartPtr->nextPtr=NULL;
                     ldPtr=ldsstartPtr;
                     strcpy(ldPtr->locus,locus); 
                     strcpy(ldPtr->kb,kb); 
                     strcpy(ldPtr->ldu,ldu); 
                     ldPtr->nextPtr=NULL; 
                     }
                     else
                     {
                     ldPtr->nextPtr=(ldsPtr)malloc(sizeof(lds));
                     ldPtr=ldPtr->nextPtr;
                     strcpy(ldPtr->locus,locus); 
                     strcpy(ldPtr->kb,kb); 
                     strcpy(ldPtr->ldu,ldu); 
                    ldPtr->nextPtr=NULL; 
                     }
           /****FILL STRUCTURE !!!! *****************/
           /****FILL STRUCTURE !!!! *****************/
    }
/*!!!!!!!!!!!!!!!!!!*/
}
fclose(fp2);
 
if(bigx<0.01){printf("\nERROR INPUT LDU MAP NOT FOUND - EXITING ");exit(0);}

fprintf(output_f,"\n\n** The number of intervals with LDU > 2.5 in this LD map is:  %d \n\n",nholes);
fflush(output_f);
}
