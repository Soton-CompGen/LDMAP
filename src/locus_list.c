
#include "allass.h"
/*********************************************************************************************/
/************************************************************************/
void locus_list()
{
int i,ii;
double kb1,kb2,ldu1,ldu2;
abcdPtr abcd_p1;
char locs[20];
/*************GET THE LOCATIONS FOR ALL IN THE LOCUS LIST ***/
/*NOTE THAT THIS LIST IS NOT SORTED BY LOCATION **/
/*When called from ldmap2.c the g_location[ ][0] array holds the LD map and the
g_location[ ][1] array holds the kb map */
/*USES INT FILE DATA TO COMPILE LOCUS LIST - EACH LOCUS FORMS PART OF A PAIR SO DEAL WITH
THE FIRST ONE AND THEN MOVE ON TO SECOND, NOTE THAT INT FILE DOES NOT CONTAIN ANY LDU DATA
THE FIRST TIME THIS IS CALLED - SO ZEROS ARE ENTERED IN INTFILE3.C*/
/*THE TWO ARRAYS OUTPUT ARE g_loci AND g_location */

kb1=0;
kb2=0;
ldu1=0;
ldu2=0;

abcd_p1=abcdstartPtr;
while(abcd_p1!=NULL)
{
for(i=0;i<20;i++)locs[i]='\0';
for(i=0;i<20;i++) { 
if(abcd_p1->locus[i]=='.')goto gogo4; strcpy(locs,abcd_p1->locus);kb1=abcd_p1->kb1;ldu1=abcd_p1->ldu1;}

gogo4: for(ii=0;ii<g_nloci;ii++)
{
if(strcmp(g_loci[ii],locs)==0){g_location[ii][0]=ldu1;g_location[ii][1]=kb1;}
}
abcd_p1=abcd_p1->nextPtr;
}

abcd_p1=abcdstartPtr;
while(abcd_p1!=NULL)
{
for(i=0;i<20;i++)locs[i]='\0';
for(i=0;i<20;i++) { if(abcd_p1->locus2[i]=='.')goto gogo3; strcpy(locs,abcd_p1->locus2);kb2=abcd_p1->kb2;ldu2=abcd_p1->ldu2;}
gogo3: for(ii=0;ii<g_nloci;ii++)
{
if(strcmp(g_loci[ii],locs)==0){g_location[ii][0]=ldu2;g_location[ii][1]=kb2;}
}
abcd_p1=abcd_p1->nextPtr;
}
/*************GET THE LOCATIONS FOR ALL IN THE LOCUS LIST ***/
reorder();

}

