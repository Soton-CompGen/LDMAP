#include "allass.h"


void getpa(char buffer[500])
{
  char ieval[20],ilval[20],imval[20]; 
  int i,n,lp;

get_predl();
global_l=g_pred;
 
  /* PA CONTROL ****/
      for(n=0; n<500; n++)
        {
         if(buffer[n]=='E' && buffer[n+1]=='=')
           { i=0; lp=n+2;
             while(buffer[lp]!=',' && buffer[lp]!=')') {ieval[i]=buffer[lp];i++;lp++;}
             ieval[i]='\0';
             if(ldflag==0)global_e=atof(ieval);
            }

         if(buffer[n]=='L' && buffer[n+1]=='=')
           {i=0; lp=n+2;
            while(buffer[lp]!=',' && buffer[lp]!=')') {ilval[i]=buffer[lp];i++;lp++;}
            ilval[i]='\0';
            if(ldflag==0)global_l=atof(ilval);
           }
         if(buffer[n]=='M' && buffer[n+1]=='=')
           { i=0; lp=n+2;
             while(buffer[lp]!=',' && buffer[lp]!=')') {imval[i]=buffer[lp];i++;lp++;}
             imval[i]='\0';
             if(ldflag==0)global_m=atof(imval);
            }
        }
}
