

#include "allass.h"
/*********************************************************************************************/
/*********************************************************************************************/
void ldmap3()
{
/*Fit Malecot to kb or LD map */

char temp[20];

ldflag2=1;
itlast=0;
iteast=0;
itmast=0;
ldflag=0;

g_abcd=abcdstartPtr;
get_predl();
printf("\nARE YOU USING A kb or LDU MAP (k/l) ? ");
scanf("%s",temp);
if(temp[0]=='k')fprintf(output_f,"\n\n\n*****   FITTING MALECOT MODEL TO THE PHYSICAL MAP (kb).................\n ");
if(temp[0]=='l')fprintf(output_f,"\n\n\n*****   FITTING MALECOT MODEL TO THE LDU MAP .................\n ");
runewt2();
fprintf(output_f,"\n\n>>>>>> STANDARD ERRORS ON THE FULL MODEL >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n\n");
itlast=1;
iteast=1;
itmast=1;
g_finish=1;
runewt3();
}
