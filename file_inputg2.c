#include "allass.h"

void file_inputg2()
{
char buffer[MAX_LINE];
int a,b,iii,ikeep=-1,iflag=0,j,i;

for(i=0;i<MAX_INDS;i++)g_used[i]=0;

if((fped=fopen(ped_file,"r"))==NULL){printf("\n Cannot open input file ");exit(1);}
fflush(output_f);
  
i=0;j=0; 

while(fgets(buffer, MAX_LINE-1, fped)!=NULL)
  {
    for(iii=0;iii<MAX_LINE;iii++){if(buffer[iii]=='\0'){ikeep=iii;goto lab1;} }
lab1:for(iii=ikeep-1;iii<MAX_LINE-1;iii++)buffer[iii]='\0';
buffer[MAX_LINE-1]='\0';
if(buffer[0]=='-'||buffer[0]=='_')iflag++;
if(iflag<2){fprintf(output_f,"\n%s",buffer);} 
if(iflag==2){fprintf(output_f,"\n%s",buffer);iflag=3;goto skip;} 
if(iflag>2)
{
strcpy(g_buffer[i],buffer);
i++;
}
skip:a=b;  }   
g_nrec=i;

simboot();
fclose(fped);
}
