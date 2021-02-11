#include "allass.h"

void file_inputnew()
{
/*READ IN GENOTYPE FILE*/
/*This is expecting a genotype file in PLINK TPED format which has chromosome number, SNP name, genetic distance in cM, and base pair position followed by geontypes for all individuals in one line. Therefore very long files of many SNPs can be handled. The chrom number and cM locations are not used by this program  */
char buffer[MAX_LINE];
int j,iset,tok, alkeep1, alkeep2, ii, ikeep=0,iii,nloci=0;
char temp1[50];
int i;

g_ped1=NULL;
g_p12=NULL;map_startPtr2=NULL;
nloci=0;
if((fped=fopen(datfile,"r"))==NULL){printf("\n Cannot open genotype file ");exit(1);}
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
while(fgets(buffer, MAX_LINE, fped)!=NULL)
  {
    for(iii=0;iii<MAX_LINE;iii++){if(buffer[iii]=='\0'){ikeep=iii;goto lab1;} }
      /*************************************************************************/
      /*************************************************************************/
lab1: buffer[ikeep+1]=' ';
        if(nloci==0)
   {
      for(ii=0;ii<50;ii++){temp1[ii]='\0';}
      j=0;
      iset=0;
      alkeep1=0;
      alkeep2=0;
      tok=0;
      /*****/
      for(i=0;i<=ikeep;i++)
        {
            if((buffer[i]!=' ')&&(buffer[i+1]!=' ')) { temp1[j]=buffer[i];j++; }
            if((buffer[i]!=' ')&&(buffer[i+1]==' '))
               {
                  temp1[j]=buffer[i];j++; 
                  tok++;
                  if(tok==2)
                    {
                    newentry2(); 
                    strcpy(g_p12->locus, temp1); 
                    }
                  if(tok==4)
                  {   
                  g_p12->kb=atof(temp1)/1000.;  
                  }
                  if(tok>4)
                    { if(iset==0){alkeep1=atoi(temp1);iset=1;goto skip;} 
                       if(iset==1){alkeep2=atoi(temp1);newentry(); 
                       g_ped1->GEN[nloci][0]=alkeep1; g_ped1->GEN[nloci][1]=alkeep2;iset=0; }
                     }
skip:;              for(ii=0;ii<50;ii++){temp1[ii]='\0';}
              j=0;
              }

         }/*for loop*/
/***IF FIRST LINE USE TO CREATE WHOLE PED STRUCTURE */
     } /*if nloci==0*/

      /*************************************************************************/
      /*************************************************************************/
     if(nloci>0)
     {
      g_ped1=ped_startPtr;
      for(ii=0;ii<50;ii++){temp1[ii]='\0';}
      j=0;
      iset=0;
      alkeep1=0;
      alkeep2=0;
      tok=0;
      /*****/
      for(i=0;i<=ikeep;i++)
        {
            if((buffer[i]!=' ')&&(buffer[i+1]!=' ')) { temp1[j]=buffer[i];j++; }
            if((buffer[i]!=' ')&&(buffer[i+1]==' '))
               {
                  temp1[j]=buffer[i];j++; 
                  tok++;
                  if(tok==2)
                    {
                    newentry2(); 
                    strcpy(g_p12->locus, temp1); 
                    }
                  if(tok==4)
                    {
                    g_p12->kb=atof(temp1)/1000.;  
                    }
                  if(tok>4)
                    {
                       if(iset==0){alkeep1=atoi(temp1);iset=1;goto skips;} 
                       if(iset==1){alkeep2=atoi(temp1);
                       g_ped1->GEN[nloci][0]=alkeep1; 
                       g_ped1->GEN[nloci][1]=alkeep2;g_ped1=g_ped1->nextPtr;iset=0; }
                     }
skips:;       for(ii=0;ii<50;ii++){temp1[ii]='\0';}
              j=0;
              }

         }/*for loop*/
     } /* if nloci >0*/
      /*************************************************************************/
      /*************************************************************************/
     nloci++;
    } /*fgets*/



g_ped1=ped_startPtr;
g_ped1->numg=nloci-1;
}

