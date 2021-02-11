#include "allass.h"

void cal_SE2(double **fjac, int n)
{

/* Inverts K matrix prior to deriving standard errors  */
int i=1, j=1, *indx;
double big=9999999.99999,**temp, **y, d, *col;
char blank[14];

if(ldflag==2)return;
strcpy(blank,"            ");
indx=ivector(1,n);
col=dvector(1,n);
y=dmatrix(1,n,1,n);
temp=dmatrix(1,n,1,n);
for(i=1;i<=n;i++){ for(j=1;j<=n;j++){ temp[i][j]=fjac[i][j];} }

fprintf(output_f,"\n\n                     Ue            Ul            Um                          ");
fprintf(output_f,"\n             %13.5f %13.5f %13.5f               ",g_uei,g_uli,g_umi); 

ludcmp(temp,n,indx, &d);
if(g_bad==1){fprintf(output_f,"\n\nCannot get standard errors ! ");return;}
for(j=1; j<=n; j++)
{
  
  for(i=1; i<=n; i++) col[i]=0.0;
      col[j]=1.0;
      lubksb(temp,n,indx,col);
  for(i=1; i<=n; i++) y[i][j]=col[i];
}   

j=1;
fprintf(output_f,"\n\nStdrd Errors:");
if(g_iter[1]==1&&(y[j][j]>=0.)&&(y[j][j]<big)){fprintf(output_f,"%13.5f",sqrt(y[j][j]));j++;
                }   
                else
                {
                fprintf(output_f," %s",blank);
                }

if(g_iter[2]==1&&(y[j][j]>=0.)&&(y[j][j]<big)){fprintf(output_f," %13.5f",sqrt(y[j][j]));j++;}   
                else
                {
                fprintf(output_f,"  %s",blank);
                }
if(g_iter[3]==1&&(y[j][j]>=0.)&&(y[j][j]<big)){fprintf(output_f," %13.5f",sqrt(y[j][j]));j++;}   
                else
                {
                fprintf(output_f,"  %s",blank);
                }

fprintf(output_f,"\n\nN (Number of pairs)= %7d     Predicted L (Lp)= %13.5f",
g_nki,sqrt(2./3.14159265)*(g_sumki2/g_nki)/(g_sumki/g_nki)); 
}
