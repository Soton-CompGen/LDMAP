
#include "allass.h"
/*********************************************************************************************/
/*********************************************************************************************/
void intfile4()
{
/*Read from intermed file the map - used in segments logic only ... */
/*THIS IS FOR FITTING REGIONS, LOCI */ 
/*READ IN FULL INT FILE AND ALSO CONSTRUCT A LOCUS LIST - LATER RE-ORDERED TO A MAP*/
/*THIS ONE JUST GETS A FULL MAP FROM AN INT FILE*/
int ic=0,i,j;
char temp[20];
char buff[200];
char cbuff[200];
char mint[30]; 
FILE *fp;
abcdPtr abcd_p1;
char locs[20],loc1[20],loc2[20];
double f1,f2,kb1,kb2;
g_nallloci=0;


if((output_ft=fopen(intefile,"w"))==NULL){ printf("\nCannot open map (output) file");exit(1);}
if((fp=fopen(datfile,"r"))==NULL){ printf("\nCannot open intermediate file for reading");exit(1);}
fprintf(output_f,"\nReading the intermediate data file: %s ",datfile);
fprintf(output_f,"\nWriting the map file: %s ",intefile);
fprintf(output_ft,"\n#Reading the intermediate data file: %s ",datfile);
fprintf(output_ft,"\n#Writing the map file: %s ",intefile);

abcd_p1=NULL;

while(fgets(buff,200,fp)!=NULL)
{
  if(buff[0]=='\n')goto endit2; 
  ic++;
strcpy(cbuff,buff);
/*******************/
  if(ic==1){
     fprintf(output_f,"\n\nFirst line of the intermediate file:");
     fprintf(output_f,"\n%s",cbuff);
     for(i=0;i<30;i++)mint[i]='\0';

     for(i=0;i<200;i++)
       {
         if(cbuff[i]=='n'&&cbuff[i+1]=='t')
          {
            mint[0]=cbuff[i+2]; mint[1]=cbuff[i+3]; mint[2]=cbuff[i+4]; mint[3]=cbuff[i+5]; mint[4]=cbuff[i+6]; 
            mint[5]=cbuff[i+7]; mint[6]=cbuff[i+8]; mint[7]=cbuff[i+9]; mint[8]='\0';
          }
       }
       fprintf(output_f,"\n\nThe maximum number of intervals for computing each epsilon is:  %s ",mint);

       for(i=0;i<30;i++)mint[i]='\0';
       for(i=0;i<200;i++)
         {
            if(cbuff[i]=='K'&&cbuff[i+1]=='b')
              {
                mint[0]=cbuff[i+2]; mint[1]=cbuff[i+3]; mint[2]=cbuff[i+4]; mint[3]=cbuff[i+5]; mint[4]=cbuff[i+6]; 
                mint[5]=cbuff[i+7]; mint[6]=cbuff[i+8]; mint[7]=cbuff[i+9]; mint[8]=cbuff[i+10]; mint[9]=cbuff[i+11];mint[10]='\0';
              }
         }
         fprintf(output_f,"\n\nThe maximum distance (kb) between any pair is set to %s ",mint);
       }
/*******************/
for(i=0;i<20;i++) { loc1[i]='\0'; loc2[i]='\0'; }
kb1=kb2=0.;

for(i=0;i<20;i++)temp[i]='\0';
j=0;
for(i=0;i<15;i++){if(buff[i]!=' '){temp[j]=buff[i];j++;}}
strcpy(loc1,temp);
if(loc1[0]==' '||loc1[0]=='\0'||loc1[0]=='\n'){goto endit2;}

for(i=0;i<20;i++)temp[i]='\0';
j=0;
for(i=16;i<31;i++){if(buff[i]!=' '){temp[j]=buff[i];j++;}}
strcpy(loc2,temp);
if(loc2[0]==' '||loc2[0]=='\0'||loc2[0]=='\n'){goto endit2;}

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
for(i=163;i<175;i++){temp[j]=buff[i];j++;}
f1=atof(temp);
 
for(i=0;i<20;i++)temp[i]='\0';
j=0;
for(i=177;i<188;i++){temp[j]=buff[i];j++;}
f2=atof(temp);
     
/*MAKE A LOCUS LIST*/
for(i=0;i<20;i++)locs[i]='\0';
for(i=0;i<20;i++) { if(loc1[i]=='.')goto gogo; locs[i]=loc1[i]; }
gogo:for(i=0;i<g_nallloci;i++) { if(strcmp(locs,g_allloci[i])==0)goto endit; }
strcpy(g_allloci[g_nallloci],locs);
g_alllocation[g_nallloci][1]=kb1;
g_alllocation[g_nallloci][2]=f1;

g_nallloci++;

endit: for(i=0;i<20;i++)locs[i]='\0';
for(i=0;i<20;i++) { if(loc2[i]=='.')goto gogo2; locs[i]=loc2[i];  }
gogo2:for(i=0;i<g_nallloci;i++) { if(strcmp(locs,g_allloci[i])==0)goto endit2; }
strcpy(g_allloci[g_nallloci],locs);
g_alllocation[g_nallloci][1]=kb2;
g_alllocation[g_nallloci][2]=f2;
g_nallloci++;


endit2:i=1; 
} /*end while*/
fclose(fp);
reorderall();

fflush(output_f);
}
