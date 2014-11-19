#include "allass.h"

void file_inputg()
{
/*Local variables*/
char buffer[MAX_LINE],loco[50];
int col1o,col2o,ikeep=0,iii,iflag,iset,j;
int i;
double kbo;
mapPtr2 p1,pold;
g_p12=NULL;map_startPtr2=NULL;
g_f1=NULL;
g_ind=0;
if((fped=fopen(ped_file,"r"))==NULL){printf("\n Cannot open genotype file ");exit(1);}
fflush(output_f3);
iflag=0;
while(fgets(buffer, MAX_LINE-1, fped)!=NULL)
  {
    for(iii=0;iii<MAX_LINE;iii++){if(buffer[iii]=='\0'){ikeep=iii;goto lab1;} }
 lab1:for(iii=ikeep-1;iii<MAX_LINE-1;iii++)buffer[iii]='\0';
 buffer[MAX_LINE-1]='\0';

 if(buffer[0]=='-'||buffer[0]=='_')iflag++;
         if(iflag<1)fprintf(output_f3,"\n%s",buffer); 

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
              g_temp[j]=buffer[i];
              j++; 
             }
          lab: i++;
 
           } 
          g_temp[j]='\0';
          make_map4(); 
       }   

    }

/*******************************/
/*******************************/
/*CHECK MAP HERE */
/*SWITCH MIS_ORDERED LOCI */
lab2:p1=map_startPtr2;
pold=NULL;

while(p1!=NULL)
{

if(pold!=NULL)
{
if(pold->kb>p1->kb)
{
strcpy(loco,pold->locus);
kbo=pold->kb;
col1o=pold->col1;
col2o=pold->col2;

strcpy(pold->locus,p1->locus);
pold->col1=p1->col1;
pold->col2=p1->col2;
pold->kb=p1->kb;

strcpy(p1->locus,loco);
p1->col1=col1o;
p1->col2=col2o;
p1->kb=kbo;

goto lab2;
}

}
pold=p1;
p1=p1->nextPtr;
}


/*******************************/
g_f1=field_startPtr;
while(g_f1!=NULL)
 {
 if(g_f1->col2==0)g_f1->col2=g_f1->col1; 

if(strcmp(g_f1->field,"N")!=0)
{ 

if(g_f1->col1==0) { 
printf("\n*****PROBLEM READING COLUMN LOCATIONS FOR THIS FIELD : %s ",g_f1->field); 
fprintf(output_f3,"\n*****PROBLEM READING COLUMN LOCATIONS FOR THIS FIELD : %s ",g_f1->field); 
                  }
printf("\n*****WARNING UNRECOGNISED FIELD IN PEDIGREE FILE : %s ",g_f1->field);
fprintf(output_f3,"\n*****WARNING UNRECOGNISED FIELD IN PEDIGREE FILE : %s ",g_f1->field);
                  }
 g_f1=g_f1->nextPtr; 
 }

/********************************/
/********************************/

readped();

}
/****************************************************************************************/ 
void pastdash()
{ /*Skip past two lines of dashes to allow further processing of a file */
char buffer[MAX_LINE];
int iflag;
rewind(fped);
iflag=0;
while(fgets(buffer,MAX_LINE-1, fped)!=NULL)
  {
if(buffer[0]=='-'||buffer[0]=='_')iflag++;
         if(iflag==2)return;
  }
} 
/****************************************************************************************/ 
void make_map4()
{
/*Process the format control - (locus name, field position, distance in kb) */
int jjj,jjjj,i,j,iset;
char locus[50],ccol1[50],ccol2[50],kkb[50];
double kb;
int jj,col1, col2;

i=0,j=0;jj=0;jjj=0;jjjj=0;
iset=0;


i=0;
while(g_temp[i]!='\0'){if(g_temp[i]==',')goto lab;i++;}return;


lab:i=0;
while(g_temp[i]!='\0')
  {
    
    if(iset==0)if(g_temp[i]==','){i++;iset++;}/*Allow dash in locus name */ 
    if(iset>0)if(g_temp[i]==','||g_temp[i]=='-'){i++;iset++;} 
     if(iset==0) 
     {
     locus[j]=g_temp[i]; 
     j++; 
     }
     locus[j]='\0';
     if(iset==1) 
     {
     ccol1[jj]=g_temp[i]; 
     jj++; 
     }
     ccol1[jj]='\0';
     if(iset==2) 
     {
     ccol2[jjj]=g_temp[i]; 
     jjj++; 
     }
     ccol2[jjj]='\0';
     if(iset==3) 
     {
     kkb[jjjj]=g_temp[i]; 
     jjjj++; 
     }
     kkb[jjjj]='\0';
    i++;
 }

col1=atoi(ccol1);
col2=atoi(ccol2);
if(col2>MAX_LINE){printf("\nERROR DATA FILE RECORDS EXCEED MAX_LINE IN LENGTH ");exit(0);}

kb=999999999.;
if(jjjj>0)kb=atof(kkb);
/*We got the record parsed - add to structure */
if(jjjj>0)
  {
   if(g_p12==NULL)
     {
      map_startPtr2=(mapPtr2)malloc(sizeof(mapkb2));
      map_startPtr2->nextPtr=NULL;
      g_p12=map_startPtr2;
      strcpy(g_p12->locus,locus);
      g_p12->col1=col1;
      g_p12->col2=col2;
      g_p12->kb=kb;
      g_ind++; 
      } 
      else
      {
      g_p12->nextPtr=(mapPtr2)malloc(sizeof(mapkb2));
      g_p12=g_p12->nextPtr;
      strcpy(g_p12->locus,locus);
      g_p12->col1=col1;
      g_p12->col2=col2;
      g_p12->kb=kb;
      g_ind++; 
      /*g_loci++;*/
      } 
   }
if(jjjj==0)
  {
   if(g_f1==NULL)
     {
      field_startPtr=(fieldPtr)malloc(sizeof(fields1));
      field_startPtr->nextPtr=NULL;
      g_f1=field_startPtr;
      strcpy(g_f1->field,locus);
      g_f1->col1=col1;
      g_f1->col2=col2;
      g_f1->nextPtr=NULL; 
      } 
      else
      {
      g_f1->nextPtr=(fieldPtr)malloc(sizeof(fields1));
      g_f1=g_f1->nextPtr;
      strcpy(g_f1->field,locus);
      g_f1->col1=col1;
      g_f1->col2=col2;
      g_f1->nextPtr=NULL; 
      } 
   }

}

/*************************************************************/
void readped()
{
/*READ IN COMPLETE PED FILE */
fieldPtr f1;
char temp[50];
char blan[50],dot[50];
char temp1[50];
char temp2[50];
int iii,ikeep=0,ic,i,j,k,iset;
mapPtr2 p1;
pedPtr ped1=NULL;
char buffer[MAX_LINE]; 
char bufferc[MAX_LINE]; 

pastdash();
while(fgets(buffer, (MAX_LINE-1), fped)!=NULL)
  {
    for(iii=0;iii<MAX_LINE;iii++){if(buffer[iii]=='\0'){ikeep=iii;goto lab1;} }
 lab1:for(iii=ikeep-1;iii<(MAX_LINE-1);iii++)buffer[iii]='\0';
 buffer[MAX_LINE-1]='\0';
/****CREATE a new record in genotypes structure */
   strcpy(bufferc,buffer);

   if(ped1==NULL)
     {
      ped_startPtr=(pedPtr)malloc(sizeof(peds));
      ped_startPtr->nextPtr=NULL;
      ped1=ped_startPtr;
      ped1->N=999; 
      } 
      else
      {
      ped1->nextPtr=(pedPtr)malloc(sizeof(peds));
      ped1=ped1->nextPtr;
      ped1->N=999; 
     }
/****END CREATE **********************************/

/*****Get the genotypes ! *********/
ic=0;
  p1=map_startPtr2;
  while(p1!=NULL)
    {
       
       strcpy(buffer,bufferc);

       iset=1; 
       j=0; 
       k=0; 
       for(i=0;i<50;i++){blan[i]='\0';temp1[i]='\0';temp2[i]='\0';dot[i]='\0';}
       strcpy(dot,".\0"); 
       strcpy(blan," \0"); 
       if(g_opt==1&&g_code==1){strcpy(temp1,"9");}
     
       for(i=(p1->col1-1);i<p1->col2;i++)
         {
           if(buffer[i]==' ') 
             { 
               if(j>0)iset=2; 
             } 
           else
             {
               if(iset!=2)
                 {
                   temp1[j]=buffer[i];blan[j]=' ';j++;
                 } 
               if(iset==2)
                {
                   temp2[k]=buffer[i];k++;
                } 
             } 
         } 
      if(g_opt==1&&g_code==0){ if(strcmp(temp1,dot)==0){strcpy(temp1,"0");}}
      if(g_opt==1&&g_code==1){ if(strcmp(temp1,blan)==0){strcpy(temp1,"9");}}

      if(strcmp(temp2,dot)==0){strcpy(temp2,"0");}
      
      if(g_opt==1&&g_code==0)
      {
      ped1->GEN[ic][0]=0; 
      ped1->GEN[ic][1]=0; 
      ped1->GEN[ic][0]=atoi(temp1); 
      }
 
      if(g_opt==1&&g_code==1)
      {
      ped1->GEN[ic][0]=9; 
      ped1->GEN[ic][1]=0; 
      ped1->GEN[ic][0]=atoi(temp1); 
      } 
      
      if(g_opt!=1)
      {

      ped1->GEN[ic][0]=atoi(temp1); 
      ped1->GEN[ic][1]=atoi(temp2); 
      if(ped1->GEN[ic][0]==0)ped1->GEN[ic][1]=0; 
      if(ped1->GEN[ic][1]==0)ped1->GEN[ic][0]=0; 
      }  

  ped1->numg=ic; 
   ic++;
   if(ic==(MAX_LOCI-5)){printf("\n ERROR INCREASE MAX_LOCI, CURRENTLY SET AT : %d ",MAX_LOCI);exit(0);}

   p1=p1->nextPtr;
    }
/*****End get the genotypes ! *********/
/****GET THE REST OF PED DATA *********/
f1=field_startPtr;
while(f1!=NULL)
  { 
for(i=0;i<50;i++){temp[i]='\0';}
if(f1->col2==0)f1->col2=f1->col1; 
j=0;for(i=f1->col1-1;i<f1->col2;i++){if(buffer[i]!=' '){temp[j]=buffer[i];j++;}}
 if(strcmp(f1->field,"FAM")==0)strcpy(ped1->FAM,temp);
 if(strcmp(f1->field,"ID")==0)strcpy(ped1->ID,temp); 
 if(strcmp(f1->field,"N")==0){ped1->N=atoi(temp);} 
  f1=f1->nextPtr;
  }
/****END GET THE REST OF PED DATA *****/
  } /*End of fgets */
}
