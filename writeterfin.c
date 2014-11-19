#include "allass.h"
/*************************************************************/
void writeterfin()
{
double diff,keepe,keepm,keepl,sumdf=0.,sumln=0.,freq1,freq2,dplus,d,sed,ed,e,kb1,kb2;
char flank1[20],copy1[20],flank2[20],copy2[20];
int jth,ii;
abcdPtr p1;
            
keepe=global_e;
keepm=global_m;
keepl=global_l;

/*************************************************************/
/*************************************************************/
dplus=0;
sed=0.;
jth=0;
if(path==1||path==3)
{

fprintf(output_ft,"\n# ");
fprintf(output_ft,"\n#                                                                                                                                  ");
fprintf(output_ft,
"\n#           Locus     kb map       LDU map      MAF                                                                                ");
fprintf(output_ft,"\n# ");
}





for(ii=0;ii<g_nloci;ii++)
{
   if(ii+1<g_nloci)
   { 
   strcpy(flank1,g_loci[ii]);
   strcpy(flank2,g_loci[ii+1]); 
   kb1=g_location[ii][1];
   kb2=g_location[ii+1][1]; 
   freq1=g_location[ii][2];
   freq2=g_location[ii+1][2];

if(jth==0){
          global_e=keepe;
          global_m=keepm;
          global_l=keepl;
          strcpy(copy1,flank1); 
          meankchi(copy1); 
          sumdf=sumdf+ g_meank;
          sumln=sumln+g_meanchi; 

          p1=abcdstartPtr;
          while(p1!=NULL)
           {
           p1->flag=0;
           if((strcmp(flank1,p1->locus)==0)||(strcmp(flank1,p1->locus2)==0))p1->flag=1;
           p1=p1->nextPtr;
           }

           g_location[ii][0]=sed;

           ldflag=2;



/**********FIX E FIRST*/
           global_e=keepe;
           global_m=keepm;
           global_l=keepl; 
           itmast=1;
           iteast=0;
           itlast=1;
           if(g_estimatel==0){itlast=0;global_l=g_pred;}
           runewt3();

           itmast=1;
           iteast=1;
           itlast=1;
           if(g_estimatel==0){itlast=0;global_l=g_pred;}
           runewt3();
           itmast=1;
           iteast=1;
           itlast=1;
           if(g_estimatel==0){itlast=0;global_l=g_pred;}
           runewt3();
          diff=(g_meanchi/g_meank)-(g_lnl/g_meank);
          fprintf(output_ft,"\n%5d %12s %12.5f %12.6f  %10.8f                                                       ",
jth+1,flank1,kb1,sed,freq1);
          
          
          p1=abcdstartPtr;
          while(p1!=NULL)
           {
           p1->flag=1;
           p1=p1->nextPtr;
           }
           global_e=keepe;
           global_m=keepm;
           global_l=keepl; 
          }

   jth=jth+1; 
   d=interv[jth][1];
   dplus=dplus+d; 
   e=interv[jth][2];
   ed=d*e;
   sed=sed+ed;
if(jth>0){
          global_e=keepe;
          global_m=keepm;
          global_l=keepl;
         strcpy(copy2,flank2);
         meankchi(copy2); 
         sumdf=sumdf+ g_meank;
         sumln=sumln+g_meanchi; 
          
          p1=abcdstartPtr;
          while(p1!=NULL)
           {
           p1->flag=0;
           if((strcmp(flank2,p1->locus)==0)||(strcmp(flank2,p1->locus2)==0))p1->flag=1;
           p1=p1->nextPtr;
           }

           g_location[ii+1][0]=sed;

           ldflag=2;

/**********FIX E FIRST*/
           global_e=keepe;
           global_m=keepm;
           global_l=keepl; 
           itmast=1;
           iteast=0;
           itlast=1;
           if(g_estimatel==0){itlast=0;global_l=g_pred;}
           runewt3();

           itmast=1;
           iteast=1;
           itlast=1;
           if(g_estimatel==0){itlast=0;global_l=g_pred;}
           runewt3();
           itmast=1;
           iteast=1;
           itlast=1;
           if(g_estimatel==0){itlast=0;global_l=g_pred;}
           runewt3();
          diff=(g_meanchi/g_meank)-(g_lnl/g_meank);
         fprintf(output_ft,"\n%5d %12s %12.5f %12.6f  %10.8f                                                        ",
jth+1,flank2,kb2,sed,freq2);
          p1=abcdstartPtr;
          while(p1!=NULL)
           {
           p1->flag=1;
           p1=p1->nextPtr;
           }
           global_e=keepe;
           global_m=keepm;
           global_l=keepl; 
         }
   }
fflush(output_ft);
}

fprintf(output_ft,"\n#                                                                                          ");

}
