#include "allass.h"

void file_input2(char *datfile)
{

/*Local variables*/
FILE *fp;
char buffer[500];
int ikeep=0,iii,iflag,iset,j;
int i;

g_ind=0;

if ((fp=fopen(datfile, "r"))==NULL){
printf("\nCannot open datafile.\n\n"); exit(1);}
iflag=0;

while(fgets(buffer, 499, fp)!=NULL)
  {
for(iii=0;iii<500;iii++){if(buffer[iii]=='\0'){ikeep=iii;goto labb;} }
 labb:for(iii=ikeep-1;iii<499;iii++)buffer[iii]=' ';
 buffer[499]='\0';
if(buffer[0]=='-'||buffer[0]=='_')iflag++;
         if(iflag<1){buffer[120]='\0';
               }
     if(iflag==1)
       {
        iset=0; 
         i=0;j=0; 
           
          while(buffer[i]!='\0')
           {
           if(buffer[i]=='('){iset++;goto lab;}
           if(buffer[i]==')'){iset++;goto lab;}
             
            if(iset==1)
             {
              gg_temp[j]=buffer[i];
              j++; 
             }
          lab: i++;
 
           } 
          gg_temp[j]='\0';

      
          getfields(); 
       }
  
       if(iflag==2)
        {
           if(buffer[0]!='-'&&buffer[0]!='_')
            {
            for(iii=0;iii<500;iii++)gg_temp[iii]=' '; 
            strcpy(gg_temp,buffer);
            get_tab();
            }
         }
  }
fclose(fp);

}
/****************************************************************************************/ 
/****************************************************************************************/ 
void getfields()
{
/*Process the format control - (locus name, field position) */
int jjj,jjjj,i,j,iset;
char locus[15],ccol1[15],ccol2[15];
int jj,col1, col2;



i=0,j=0;jj=0;jjj=0;jjjj=0;
iset=0;

i=0;
while(gg_temp[i]!='\0'){if(gg_temp[i]==',')goto lab;i++;}return;



lab:i=0;
while(gg_temp[i]!='\0')
  {
    
    if(gg_temp[i]==','||gg_temp[i]=='-'){i++;iset++;} 
     if(iset==0) 
     {
     locus[j]=gg_temp[i]; 
     j++; 
     }
     locus[j]='\0';
     if(iset==1) 
     {
     ccol1[jj]=gg_temp[i]; 
     jj++; 
     }
     ccol1[jj]='\0';
     if(iset==2) 
     {
     ccol2[jjj]=gg_temp[i]; 
     jjj++; 
     }
     ccol2[jjj]='\0';
    i++;
 }

col1=atoi(ccol1);
col2=atoi(ccol2);
/*We got the record parsed - add to structure */

   if(g_p1==NULL)
     {
      fields_startPtr=(fieldsPtr)malloc(sizeof(fields));
      fields_startPtr->nextPtr=NULL;
      g_p1=fields_startPtr;
      strcpy(g_p1->locus,locus);
      g_p1->col1=col1;
      g_p1->col2=col2;
      } 
      else
      {
      g_p1->nextPtr=(fieldsPtr)malloc(sizeof(fields));
      g_p1=g_p1->nextPtr;
      strcpy(g_p1->locus,locus);
      g_p1->col1=col1;
      g_p1->col2=col2;
      } 

}
/***********************************************************************************/
/***********************************************************************************/
void get_tab()
{
/*Process the file to get all the 9 counts (diplo) & fill the structure */
int j=0,k=0,kk=0,i=0,iii;
char temp[500];

for(iii=0;iii<500;iii++)temp[iii]=' ';


if(g_aiPtr==NULL)
  {
      aistartPtr=(aisPtr)malloc(sizeof(ais));
      aistartPtr->nextPtr=NULL;
      g_aiPtr=aistartPtr;
      
      g_p1=fields_startPtr; 
      while(g_p1!=NULL)
        {
                  
                  j=0; 
                  for(i=g_p1->col1-1;i<=g_p1->col2-1;i++) { if(gg_temp[i]!=' '){temp[j]=gg_temp[i]; j++;} }
                       temp[j]='\0';
                    if(g_p1->locus[0]=='k'&& g_p1->locus[1]=='b'&&g_p1->locus[2]=='1'){strcpy(g_aiPtr->ckb1,temp);g_aiPtr->kb1=atof(temp);goto l1;} 
                    if(g_p1->locus[0]=='k'&& g_p1->locus[1]=='b'&&g_p1->locus[2]=='2'){strcpy(g_aiPtr->ckb2,temp);g_aiPtr->kb2=atof(temp);goto l1;} 
                    if(g_p1->locus[0]=='N'){ g_aiPtr->n=atof(temp);goto l1;}
                    if(g_p1->locus[0]=='l'&& g_p1->locus[5]=='1'){strcpy(g_aiPtr->locus1,temp);goto l1;}
                    if(g_p1->locus[0]=='l'&& g_p1->locus[5]=='2'){strcpy(g_aiPtr->locus2,temp);goto l1;}
                    if(g_p1->locus[0]=='f'&& g_p1->locus[4]=='1')g_aiPtr->freq1=atof(temp);
                    if(g_p1->locus[0]=='f'&& g_p1->locus[4]=='2')g_aiPtr->freq2=atof(temp);
                    else {g_aiPtr->aitab[kk]=atof(temp);kk++;} 
l1:                 g_p1=g_p1->nextPtr;
        }

  }


else  /* If not the first record in the ais structure */
  {
      g_aiPtr->nextPtr=(aisPtr)malloc(sizeof(ais));
      g_aiPtr=g_aiPtr->nextPtr;
      g_p1=fields_startPtr;
      j=0; 
      while(g_p1!=NULL)
        {
             j=0; 
              for(i=g_p1->col1-1;i<=g_p1->col2-1;i++) { if(gg_temp[i]!=' '){temp[j]=gg_temp[i]; j++;} }
                    temp[j]='\0';
                    if(g_p1->locus[0]=='k'&& g_p1->locus[1]=='b'&&g_p1->locus[2]=='1'){strcpy(g_aiPtr->ckb1,temp);g_aiPtr->kb1=atof(temp);} 
                    if(g_p1->locus[0]=='k'&& g_p1->locus[1]=='b'&&g_p1->locus[2]=='2'){strcpy(g_aiPtr->ckb2,temp);g_aiPtr->kb2=atof(temp);} 
                    if(g_p1->locus[0]=='N') g_aiPtr->n=atof(temp);
                    if(g_p1->locus[0]=='l'&& g_p1->locus[5]=='1'){strcpy(g_aiPtr->locus1,temp);}
                    if(g_p1->locus[0]=='l'&& g_p1->locus[5]=='2'){strcpy(g_aiPtr->locus2,temp);}
                    if(g_p1->locus[0]=='f'&& g_p1->locus[4]=='1')g_aiPtr->freq1=atof(temp);
                    if(g_p1->locus[0]=='f'&& g_p1->locus[4]=='2')g_aiPtr->freq2=atof(temp);
                    else {g_aiPtr->aitab[kk]=atof(temp);kk++;} 
           g_p1=g_p1->nextPtr;
           k++; 
         }

  }
/*************************/
}
