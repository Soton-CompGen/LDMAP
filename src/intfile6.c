
#include "allass.h"
/*********************************************************************************************/
/*********************************************************************************************/
void intfile6()
{
/*Read from intermed file the values of k, p etc.... */
/*THIS IS FOR FITTING REGIONS, LOCI */ 
/*READ IN FULL INT FILE AND ALSO CONSTRUCT A LOCUS LIST - LATER RE-ORDERED TO A MAP*/

int ic=0,i,j;
char temp[20];
char buff[200];
char cbuff[200];
FILE *fp;
abcdPtr abcd_p1;
char locs[20],loc1[20],loc2[20];
double chi,n,kb1,kb2,p,k;
g_nloci=0;

if((fp=fopen(intefile,"r"))==NULL){ printf("\nCannot open intermediate file for reading");exit(1);}
fprintf(output_f,"\nReading the intermediate data file: %s ",intefile);
fprintf(output_f,"\nWriting the map file: %s ",terfile);

abcd_p1=NULL;

while(fgets(buff,200,fp)!=NULL)
{
  if(buff[0]=='\n')goto skip; 
  ic++;
strcpy(cbuff,buff);
/*******************/
for(i=0;i<20;i++) { loc1[i]='\0'; loc2[i]='\0'; }
kb1=kb2=p=k=chi=n=0.;

for(i=0;i<20;i++)temp[i]='\0';
j=0;
for(i=0;i<15;i++){if(buff[i]!=' '){temp[j]=buff[i];j++;}}
strcpy(loc1,temp);
if(loc1[0]==' '||loc1[0]=='\0'||loc1[0]=='\n'){goto skip;}

for(i=0;i<20;i++)temp[i]='\0';
j=0;
for(i=16;i<31;i++){if(buff[i]!=' '){temp[j]=buff[i];j++;}}
strcpy(loc2,temp);
if(loc2[0]==' '||loc2[0]=='\0'||loc2[0]=='\n'){goto skip;}

for(i=0;i<20;i++)temp[i]='\0';
j=0;
for(i=32;i<42;i++){temp[j]=buff[i];j++;}
kb1=atof(temp);

for(i=0;i<20;i++)temp[i]='\0';
j=0;
for(i=43;i<53;i++){temp[j]=buff[i];j++;}
kb2=atof(temp);

for(i=0;i<20;i++)temp[i]='\0';
j=0;
for(i=54;i<66;i++){temp[j]=buff[i];j++;}
p=atof(temp);

for(i=0;i<20;i++)temp[i]='\0';
j=0;
for(i=68;i<78;i++){temp[j]=buff[i];j++;}
k=atof(temp);
/*MAKE A LOCUS LIST*/
for(i=0;i<20;i++)locs[i]='\0';
for(i=0;i<20;i++) { if(loc1[i]=='.')goto gogo; locs[i]=loc1[i]; }
gogo:for(i=0;i<g_nloci;i++) { if(strcmp(locs,g_loci[i])==0)goto endit; }
strcpy(g_loci[g_nloci],locs);
g_nloci++;

endit: for(i=0;i<20;i++)locs[i]='\0';
for(i=0;i<20;i++) { if(loc2[i]=='.')goto gogo2; locs[i]=loc2[i];  }
gogo2:for(i=0;i<g_nloci;i++) { if(strcmp(locs,g_loci[i])==0)goto endit2; }
strcpy(g_loci[g_nloci],locs);
g_nloci++;


endit2:if(k<=0.0001)goto skip;

if(abcd_p1==NULL)
{
abcdstartPtr=(abcdPtr)malloc(sizeof(abcds));
abcdstartPtr->nextPtr=NULL;
abcd_p1=abcdstartPtr;
strcpy(abcd_p1->locus,loc1);
strcpy(abcd_p1->locus2,loc2);
abcd_p1->k=k;
abcd_p1->p=p;
abcd_p1->ldu=0.;
abcd_p1->ldu1=0.;
abcd_p1->ldu2=0.;
abcd_p1->kb1=kb1;
abcd_p1->kb2=kb2;
abcd_p1->kb=fabs(kb2-kb1);
abcd_p1->flag=1;
}
else
{
abcd_p1->nextPtr=(abcdPtr)malloc(sizeof(abcds));
abcd_p1=abcd_p1->nextPtr;
strcpy(abcd_p1->locus,loc1);
strcpy(abcd_p1->locus2,loc2);
abcd_p1->k=k;
abcd_p1->p=p;
abcd_p1->ldu=0.;
abcd_p1->ldu1=0.;
abcd_p1->ldu2=0.;
abcd_p1->kb1=kb1;
abcd_p1->kb2=kb2;
abcd_p1->kb=fabs(kb2-kb1);
abcd_p1->flag=1;
}

skip:i=1;
} /*end while*/

fclose(fp);
fflush(output_f);
}
