
#include "allass.h"
/*********************************************************************************************/
/*********************************************************************************************/
void mergeloc()
{
/*Merge the LDU map locations with the kb ones on the raw data file */
char theldu[20],locs[20],buff[200],buffer[10000],terse[20];
char toks[10][20],temp[200];
FILE *fp,*fp2;
int ig,ncom,jcounter,ntoks,k,ic,iflag,iset,a,b,i,j,ii,jj;
ldsstartPtr=NULL;
ldPtr=NULL;
iflag=0;
printf("\nThis option requires a diplotype/haplotype data file with locations given in kb and a terse output file containing the LD map ");
printf("\nThe updated diplotype/haplotype file containing both kb and ldu maps will be written to a file called temp.gen");
printf("\nEnter name of raw genotype file: ");
scanf("%s",ped_file);
if((fped=fopen(ped_file,"r"))==NULL){printf("\n Cannot open genotype file ");exit(1);}

printf("\nEnter name of LDMAP terse output file: ");
scanf("%s",terse);
if((fp2=fopen(terse,"r"))==NULL){printf("\n Cannot open LDMAP terse file ");exit(1);}

if((fp=fopen("temp.gen","w"))==NULL){printf("\n Cannot open updated temp file ");exit(1);}
/******************************************************/
/******************************************************/
jcounter=0;
while(fgets(buff,199,fp2)!=NULL)
{
/*!!!!!!!!!!!!!!!!!!*/
   if(jcounter>10)/*first 10 lines of the file are skipped before the map is reached*/ 
    {
    iset=0;
    j=0;
    ntoks=0;
    ic=0;
    for(ii=0;ii<10;ii++){for(jj=0;jj<20;jj++)toks[ii][jj]='\0'; } 
    fflush(fp);
    for(k=0;k<200;k++)
      {
      if(buff[k]==' '||buff[k]=='\0'){ic=0;iset=0;goto skip;}
      if(buff[k]!=' ')
        {
         if(iset==1){toks[ntoks][ic]=buff[k];ic++;}
         if(iset==0){iset=1;ntoks++;toks[ntoks][ic]=buff[k];ic++;}
        }
      skip:a=b;
      }
   /*SAVE TO STRUCTURE*/ 
    fill_str2(toks[2],toks[4]);
    }
jcounter++;
/*!!!!!!!!!!!!!!!!!!*/
}
/******************************************************/
/******************************************************/
ldPtr=ldsstartPtr;
while(ldPtr!=NULL)
{

ldPtr=ldPtr->nextPtr;
}
/******************************************************/
/******************************************************/



while(fgets(buffer, 9999, fped)!=NULL)
  {
for(i=0;i<10000;i++){if(buffer[i]=='\n')buffer[i]='\0';}
 if(buffer[0]=='-'||buffer[0]=='_')iflag++;
         if(iflag<1){fprintf(fp,"\n%s",buffer);} 
/*****************************************/
     if(iflag==1)
       {
 if(buffer[0]=='-'||buffer[0]=='_'){fprintf(fp,"\n%s",buffer);} 
 if(buffer[0]!='-'&&buffer[0]!='_')
      { 
       /*NEED TO PARSE THESE BRACKET CONTROLS SO WE CAN ADD IN REVISED LDUs*/  
        for(k=0;k<20;k++)locs[k]='\0'; 
       ic=0; 
       for(k=0;k<10000;k++)
         {
         if(buffer[k]=='('){iset=1;goto skip3;}
         if(buffer[k]==',')goto skip4;
         locs[ic]=buffer[k];ic++; 


         skip3:a=b; 
         }

        skip4:a=b;
        for(k=0;k<20;k++)theldu[k]='\0'; 
        ldPtr=ldsstartPtr;
        while(ldPtr!=NULL)
        {
        if(strcmp(ldPtr->locus,locs)==0)strcpy(theldu,ldPtr->ldu);
        ldPtr=ldPtr->nextPtr;
        }
        
/*------------------------------------*/
        ncom=0;ig=0; 
        for(k=0;k<200;k++)temp[k]='\0';
        for(k=0;k<200;k++)
           {
            if(buffer[k]==')')buffer[k]=','; 
            if(buffer[k]==',')ncom++;
            temp[ig]=buffer[k];ig++; 
             if(ncom==3)
              {
              /*printf("\n********%s%s)",temp,theldu);*/
              fprintf(fp,"\n%s%s)",temp,theldu);
              ncom=4;
              }                  
           }
/*------------------------------------*/
       
    }   
    }   

/*****************************************/

     if(iflag>1)
       {
         fprintf(fp,"\n%s",buffer); 
       }   

    }
fprintf(fp,"\n");
fclose(fp);
 }

