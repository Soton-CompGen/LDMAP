

#include "allass.h"
         void meankchi(char theloc1[20])
         {
/*GET LOCUS SPECIFIC CHI SQUARE ETC */

         double lnl=0., xx,d,k,n=0,po,pe,E,L,M;
         g_abcd=abcdstartPtr;
         while(g_abcd!=NULL)
         {
          if ( (strcmp(theloc1,g_abcd->locus)==0)||(strcmp(theloc1,g_abcd->locus2)==0))
            {
            L=global_l;
            M=global_m;
            E=global_e;
            d=g_abcd->ldu;
            k=g_abcd->k;
            xx=-E*d;
            pe=(1.-L)*M*exp(xx)+L;
            po=g_abcd->p;
           
            lnl+=-(((po-pe)*(po-pe))*k)/2.;
            n=n+1.;
             }
         g_abcd=g_abcd->nextPtr;
         }
           }
