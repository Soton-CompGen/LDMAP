#include <ctype.h>
#include "allass.h"



void jobin()
{
/*GET AND PROCESS A JOB FILE */
/*Local variables*/
char ibuffer[80];
char buffer[80];
int iset, c, z;
int i;

g_sig=1;

if(itm==1||ite==1||itl==1)
{
itlast=itl;
itmast=itm;
iteast=1;
}

ite=0;
itl=0;
itm=0;
/**TOP OF WHILE LOOP **/
ibuffer[0]='\0';
for(i=0;i<80;i++)ibuffer[i]=' ';

while(fgets(ibuffer, 79, job_fp)!=NULL)
{
for(i=0;i<80;i++)buffer[i]=' ';
   /*Removes blanks and makes all alphabetical characters uppercase*/
   c=0; z=0;
   while(ibuffer[c]!='\0')
      {if(ibuffer[c]==' '){c++;}
       else {buffer[z]=toupper(ibuffer[c]);c++;z++;}
      }
      buffer[z]='\0';
/***************************************/
if(buffer[0]=='P'&&buffer[1]=='A')
{
g_sig=0;
fprintf(output_f,"\n--------------------------------------------------------------------------------------------------------------------------------"); 
if(ldflag==0)fprintf(output_f,"\n>>> %s ",buffer);
getpa(buffer);
fflush(output_f);
}

iset=0;
/***************************************/
/* IT CONTROL */
if(buffer[0]=='I' && buffer[1]=='T')
    { 

g_sig=0;

fprintf(output_f,"\n\n>>> %s ",buffer);
fflush(output_f);
    for(i=3;i<79;i++)
     {
     if(buffer[i]=='E')ite=1;  
     if(buffer[i]=='L')itl=1;  
     if(buffer[i]=='M')itm=1;  
     }  /*end for*/
    

    goto lab;
    }
/***************************************/
if(buffer[0]=='C' && buffer[1]=='C'){g_sig=1;goto lab;}

} /*End of while */

lab:fflush(output_f);
return;
}
