
#include "allass.h"
/*********************************************************************************************/
/*********************************************************************************************/
void segments()
{
/*READ IN A PORTION OF AN INT FILE - FOR SEGMENTAL MAP CONSTRUCTION */
char temp[10];
int j,i;
int nseg=-1,bign;
int overlap=-1,pos,segsize,a,b,ifin=0;
double big1,big2,offset,large=999999999.;
/*intfile4 just gets a full locus list and stores it in g_allloci and g_alllocations
- the element [i][2] or the latter has a flag which labels (by 1) loci to be used in 
the segmental construction of an ld map*/

printf("\n....Compiling locus list....");
intfile4();
printf("\nThe input map has %d loci in total ",g_nallloci);
fprintf(output_f,"\nThe input map has %d loci in total ",g_nallloci);

/*!!!!!!!!!!!!!!!!!!!!*/
/*MODIFY FOR IRIDIS - JULY 09 */
printf("\nFor files with >500 loci this program will assemble the map in overlapping segments ");

nseg=g_nallloci/500;
if(nseg<1)nseg=1;
fprintf(output_f,"\nCompiling the map in %d segments ",nseg);
printf("\nCompiling the map in %d segments ",nseg);

/*
printf("\nFor assembling the LD map how many segments are required ?"); 
scanf("%s", temp);
nseg=atoi(temp);
*/



for(i=0;i<g_nallloci;i++) { g_alllocation[i][0]=large; }

bign=g_nallloci;
segsize=bign/nseg;
segsize=segsize+1; /*ADD 1 TO MANAGE ROUNDING */

printf("\nThis has approximately %d loci per segment",segsize-1);
fprintf(output_f,"\nThis has approximately %d loci per segment",segsize-1);


/*
printf("\nWhat overlap is required ? "); 
scanf("%s", temp);
overlap=atoi(temp);
*/

overlap=25;

fprintf(output_f,"\nThe overlap is %d ",overlap);

if(overlap<=0)overlap=1;
pos=0;
a=pos-overlap;
b=pos+segsize+overlap;
/*****************************************/
top:for(i=0;i<g_nallloci;i++) { g_alllocation[i][3]=0.; }

if(a<0)a=0;
if(b>g_nallloci){b=g_nallloci;ifin=1;}
/*SET A FLAG WHICH SELECTS A SUB SET OF LOCI */
for(i=a;i<b;i++)
{
 g_alllocation[i][3]=1.;
}
fprintf(output_f,"\n------------------------------------------------------------------------------------------------------------------------");  
fprintf(output_f,"\n------------------------------------------------------------------------------------------------------------------------\n\n");  
fprintf(output_f,"\nCONSTRUCTING A MAP OF THIS SEGMENT WITH THE FOLLOWING MARKERS AS THE LIMITS : ");
fprintf(output_f,"\nBEGIN: %16s  %12.6f ",g_allloci[a],g_alllocation[a][1]);
fprintf(output_f,"\n  END: %16s  %12.6f ",g_allloci[b-1],g_alllocation[b-1][1]);

intfileseg();
path=1;
ldmapseg();
fflush(output_f);
pos=pos+segsize;
a=pos-overlap;
b=pos+segsize+overlap;

/******************************************************************************************/
/*APPLY FIRST OFFSET - THIS IS THE LOCATION IN THE CURRENT ALL MAP
WHICH CORRESPONDS TO THAT FOR THE FIRST LOCUS IN THE NEW SEGMENT*/

offset=0.;
for(i=0;i<g_nallloci;i++)
  {
    if(g_alllocation[i][0]<large)
      {

        /*First locus in the new segment is g_loci[0] - which has a location at zero*/
       	if (strcmp(g_allloci[i],g_loci[0])==0)offset=g_alllocation[i][0];
      }
  }

/*NOW ADD THIS TO THE NEW MAP*/
for(j=0;j<g_nloci;j++)
{
g_location[j][0]=g_location[j][0]+offset;

}
/******************************************************************************************/
/*NOW DEAL WITH THE OVERLAP SECTION*/
/*THIS VERSION AVERAGES THE OVERLAP*/
big1=0.;
big2=0.;
for(i=0;i<g_nallloci;i++)
  {
    if(g_alllocation[i][0]<large)
      {
        for(j=0;j<g_nloci;j++)
          {
            if(strcmp(g_allloci[i],g_loci[j])==0) { 
              if(g_location[j][0]>big2)big2=g_location[j][0]; 

            g_alllocation[i][0]=(g_location[j][0]+g_alllocation[i][0])/2.;
            if(g_alllocation[i][0]>big1)big1=g_alllocation[i][0]; 

           }
         }
      }
   }
/****************************************************************/
/******************************************************************************************/
/*NOW DEAL WITH THE NON OVERLAP SECTION WHICH FOLLOWS - NEED A SECOND SMALL OFFSET TO ACCOUNT
FOR SMALL DISCREPANCIES IN SCALES IN THE OVERLAP SECTION - MAY BE NEGATIVE*/
offset=big1-big2;
/*fprintf(output_f,"\nSECOND OFFSET = %f  big1 = %f big2 = %f ",offset,big1,big2);
*/

for(i=0;i<g_nallloci;i++)
  {
    if(g_alllocation[i][0]==large)
      {
        for(j=0;j<g_nloci;j++)
          {
            if(strcmp(g_allloci[i],g_loci[j])==0) { 
            g_alllocation[i][0]=g_location[j][0]+offset; 
                                                  }
         }
    }  
  }
/******************************************************************************************/
writeterfinseg();
if(ifin==0)goto top;

/*NOW NEED TO GET ALL THE DATA TO COMPUTE -2Lnl AGAINST BOTH Kb and LDU MAPS*/


fprintf(output_f,"\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"); 
fprintf(output_f,"\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"); 
fprintf(output_f,"\n\n\nTESTING THE FIT OF THE ASSEMBLED MAP...........");
fprintf(output_f,"\nFIRSTLY TESTING  THE FIT OF THE Kb MAP....\n\n");
for(i=0;i<g_nallloci;i++) { 
    g_location[i][0]=g_alllocation[i][0]; 
    g_location[i][1]=g_alllocation[i][1]; 
    g_location[i][2]=g_alllocation[i][2]; 
    strcpy(g_loci[i],g_allloci[i]);
                          }

g_nloci=g_nallloci;

intfile6();
path=0;
rewind(job_fp);
ldmapseg();

fprintf(output_f,"\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"); 
fprintf(output_f,"\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"); 
fprintf(output_f,"\nNOW TESTING  THE FIT OF THE LDU MAP....\n\n");
path=2;
rewind(job_fp);
ldmapseg();
writeterfinseg();
}
