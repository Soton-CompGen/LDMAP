

#include "allass.h"
/*********************************************************************************************/
/*********************************************************************************************/
void ldmap2()
{
/*Construct Ld map according to algorithm for locus oriented mapping found in map+ */
/*This version always uses predicted L and stops iterating epsilon once at 0 */
time_t start1,end1;
int stope,istop=0,intervals,if1,if2,useld;
double xxx,BESTLNL,x,y,allldu,maxkb=0,minkb,large=99999999999.,swept,diff,LDU,orig_e,conv,
old,new,bigk,u,lnl,ldu1,ldu2,kb1,kb2;
char flank1[30],flank2[30];
char temp[20];
abcdPtr abcd_p1;
hugsPtr hug1;
intsPtr intop;
int ntimes=0,maxit,ix,i,jth,j,ii;

/* set the vv variables..*/ 
/* double vvthee,vved,vvdd,vved_old,vved_new; */
double vvthee,vved,vvdd;
int vvjj,vvleft,vvright;
intsPtr vvintop;
hugsPtr vvhug1;

/** aman added   0926 ***/

g_estimatel=0;
g_initial=0;
conv=0.1;
path=0;
minkb=large;
g_best=large;
maxit=20000;
ldflag2=1;
itlast=0; iteast=0; itmast=0;
useld=0;
intervals=g_int;
stope=1;
/************************************************************************************************************/
/*strcpy(temp,"                   ");
printf("\nDo you want to change the maximum number of iterations, currently: =%d  (y/n) ?",maxit);
fprintf(output_f,"\nDo you want to change the maximum number of iterations, currently: =%d  (y/n) ?",maxit);
scanf("\n%s",temp);
fprintf(output_f, " %s ",temp);
if(temp[0]=='y')
{
strcpy(temp,"                   ");
printf("\nEnter the maximum number of iterations (default 20000): ");
scanf("\n%s",temp);
maxit=atoi(temp);
fprintf(output_f, "\nThe maximum number of iterations is now %d ",maxit);
}*/
/************************************************************************************************************/
/*strcpy(temp,"                   ");
printf("\nDo you want to continue iterating in an interval at epsilon=0 (y/n) ?");
fprintf(output_f,"\n\nDo you want to continue iterating in an interval at epsilon=0 (y/n) ?");
scanf("\n%s",temp);
fprintf(output_f, " %s ",temp);
if(temp[0]=='y')stope=0;
*/
if(stope==1)fprintf(output_f,"\nSo....epsilon is NOT iterated once at lower boundary");
if(stope==0)fprintf(output_f,"\nSo....epsilon IS iterated, even at lower boundary");
/************************************************************************************************************/
printf("\nSelect an option from 0, 1, 2 or 3 ");
strcpy(temp,"                   ");
printf("\n *0* Fit the pairwise data to the kb map (Enter 0) ");
printf("\n *1* Create an LD map from a kb map (Enter 1) ");
printf("\n *2* Fit the pairwise data to an LD map (Enter 2) ");
printf("\n *3* Create an LD map from an existing LD map (Enter 3)\n ");
scanf("\n%s",temp);
path=atof(temp);
if(path!=0&&path!=1&&path!=2&&path!=3){printf("\nERROR - YOU MUST SELECT AN OPTION ");
                             fprintf(output_f,"\nERROR - YOU MUST SELECT AN OPTION ");exit(0);}

/************************************************************************************************************/
if(path==1||path==3)
{

/*strcpy(temp,"                   ");
printf("\nDo you want to estimate L during LD map construction (y/n) ? "); 
fprintf(output_f,"\n\nDo you want to estimate L during LD map construction (y/n) ? "); 
scanf("\n%s",temp);
fprintf(output_f, " %s ",temp);
if(temp[0]=='y'){
fprintf(output_f,"\nEstimating L during LD map construction ");
 g_estimatel=1; }
*/
fprintf(output_f,"\nUsing predicted L during LD map construction ");
}

/******************************************************************************************************/
if(path==0)fprintf(output_f,"\n\n\n******** Fitting pairwise data to kb map ******** ");
if(path==1)fprintf(output_f,"\n\n\n******** Creating an LD map from a Kb map ******** ");
if(path==2)fprintf(output_f,"\n\n\n******** Fitting pairwise data to LD map ******** ");
if(path==3)fprintf(output_f,"\n\n\n******** Creating an LD map from an existing LD map ******** ");
/******************************************************************************************************/
start1=time(NULL);
if(path==2)istop=1;
if((path==2)||(path==3)) { useld=1; directmap();}

/*IF NEEDED INSERT EXISTING LDU MAP AT THIS POINT*/
if(useld==1)
{
  abcd_p1=abcdstartPtr;
  while(abcd_p1!=NULL)
    {
    if1=0;if2=0;
    ldPtr=ldsstartPtr;
    while(ldPtr!=NULL)
        {  
           if(strcmp(abcd_p1->locus,ldPtr->locus)==0){abcd_p1->kb1=atof(ldPtr->kb);abcd_p1->ldu1=atof(ldPtr->ldu);if1=1;}
           if(strcmp(abcd_p1->locus2,ldPtr->locus)==0){abcd_p1->kb2=atof(ldPtr->kb);abcd_p1->ldu2=atof(ldPtr->ldu);if2=1;}
           ldPtr=ldPtr->nextPtr;
        }
if((if1==0)||(if2==0)){printf("\nWARNING: TRYING TO USE AN EXISTING LD MAP BUT  MISSING LOCUS FROM THE INTERMED FILE...EXITING"); exit(1);}

       x=abcd_p1->kb2-abcd_p1->kb1;
       if(x<0.)x=-x;
       abcd_p1->kb=x; 
       abcd_p1->ldu=fabs(abcd_p1->ldu2-abcd_p1->ldu1);  
      abcd_p1=abcd_p1->nextPtr;
    }
locus_list(); /*this list is normally compiled in intfile3 - need to update again here, now the map has changed*/
}/*ABOVE SECTION ONLY FOR INSERTING EXISTING LDU MAP */

/******************************************************************************************************/
maxkb=0;
minkb=large; 
ldflag=0;
g_abcd=abcdstartPtr;
abcd_p1=abcdstartPtr;

while(abcd_p1!=NULL)
{

if(abcd_p1->kb1>maxkb)maxkb=abcd_p1->kb1;     
if(abcd_p1->kb2>maxkb)maxkb=abcd_p1->kb2;      
if(abcd_p1->kb1<minkb)minkb=abcd_p1->kb1;     
if(abcd_p1->kb2<minkb)minkb=abcd_p1->kb2;      

abcd_p1=abcd_p1->nextPtr;
       }
get_predl();

fflush(output_f);
if(path==0||path==1)g_initial=1;
    

runewt2();
g_initial=0;

/*******************************************/
orig_e=global_e;
swept=1./orig_e;

/****************************************************************************************************/
/*
reorder();
*/
/****************************************************************************************************/

maxkb=maxkb-minkb;
fprintf(output_f,"\n\n***THE SWEPT RADIUS (MEAN EXTENT OF 'USEFUL' LD IN THIS SAMPLE) IS %12.6f    ***",1./global_e);

if((global_l>=(g_pred-0.01))&&(global_l<=(g_pred+0.01))){g_df=g_nki-2;}
else{g_df=g_nki-3;}
g_V=g_lnl/g_df;

fprintf(output_f,"\n\nN(number of pairs)=%15d  m(number of SNPs)=%15d df=%14.1f V(error variance)=%14.5f\n\n", g_nki,g_nloci,g_df,g_V);


fflush(output_f);
if(path==0)return;

if(istop==1) { printmap(); goto fin; }

/****************************************************************************************************/
jth=0;
for(ii=0;ii<g_nloci;ii++)
{
   if(ii+1<g_nloci)
   { 
   strcpy(flank1,g_loci[ii]);
   strcpy(flank2,g_loci[ii+1]); 

/**g_location[][0] holds the LD map*/
   ldu1=g_location[ii][0];
   ldu2=g_location[ii+1][0]; 

/**g_location[][1] holds the Kb map*/
   kb1=g_location[ii][1];
   kb2=g_location[ii+1][1]; 

   jth=jth+1; 
   interv[jth][1] = kb2-kb1; 
   interv[jth][2]=global_e; 
   cinterv[jth][1] = kb2-kb1; 
   cinterv[jth][2]=global_e; 

   /*If using an existing LD map replace starting epsilons with LDU/Kb in each interval*/
   if(useld==1){x=fabs(ldu1-ldu2);y=fabs(kb1-kb2); interv[jth][2]=0.001;if(y>0.)interv[jth][2]=x/y; }  
   if(useld==1){x=fabs(ldu1-ldu2);y=fabs(kb1-kb2); cinterv[jth][2]=0.001;if(y>0.)cinterv[jth][2]=x/y; }  
   interv[jth][5]=1.0; 
   interv[jth][6]=1.0; 
   cinterv[jth][5]=1.0; 
   cinterv[jth][6]=1.0; 
   }
}

if(path==3)printinputmap();

fprintf(output_f,"\n\n\n*****   CONSTRUCTING LDU MAP.................\n ");
if(path==1)
{
if(g_estimatel==0)global_l=g_pred;
quicklike(global_e,global_l,global_m);
fprintf(output_f,"\n\n\nIter     E          L          M       n_holes   LDU_length      -2lnlk"); 
fprintf(output_f,"\n     %10.6f %10.6f %10.6f      %16.8f %12.5f ",
 global_e,global_l,global_m,maxkb*global_e,g_lnl);
}

if(path==3)
{
if(g_estimatel==0)global_l=g_pred;
quicklike(global_e,global_l,global_m);

fprintf(output_f,"\n\n\nIter     E          L          M       n_holes   LDU_length      -2lnlk"); 
fprintf(output_f,
"\n     %10.6f %10.6f %10.6f                       %12.5f ", 
global_e,global_l,global_m,g_lnl);
}

g_startlnl=g_lnl;
fflush(output_f);

/***************************************************************************************************************/
      abcd_p1=abcdstartPtr;
      while(abcd_p1!=NULL)
       {
       for(ii=0;ii<g_nloci;ii++)
         {  
          if(strcmp(abcd_p1->locus,g_loci[ii])==0)abcd_p1->i1=ii;
         }
       for(ii=0;ii<g_nloci;ii++)
         {  
          if(strcmp(abcd_p1->locus2,g_loci[ii])==0)abcd_p1->i2=ii;
         }
       abcd_p1=abcd_p1->nextPtr;
        }

g_niter=0;
keep_it();
BESTLNL=g_lnl;

/*Build the pointer relationships to the int file data - 
allows more rapid access to records informative for each interval*/
fast();

/***************************************************************************************************************/
/***************************************************************************************************************/
/*M A I N  S E C T I O N *********************/
/***************************************************************************************************************/
/***************************************************************************************************************/
g_niter=0;
old=large;
new=large;
ix=0;

topp:lnl=0.;
hug1=bigstart;
ii=0;
allldu=0;


/********** aman added ...computing the first value of ed  0925  *************/
/*THIS IS NOT LAYERED ANY MORE */  

vvhug1=bigstart;
while(vvhug1!=NULL)
{

  vvintop=vvhug1->pp;
  while(vvintop!=NULL)
    {
     vvleft=vvintop->p->i1; 
     vvright=vvintop->p->i2;
     vved=0.; 
     vvdd=0.; 
       /*Sum up adjacent interval distances and epsilon*d's for this pair*/ 
     vvjj=vvleft+1;
     while(vvjj<=vvright && vvjj<=g_nloci)
       {
        vvdd=vvdd+interv[vvjj][1];
        vvthee=interv[vvjj][2];
        vved=vved+interv[vvjj][1]*vvthee;
        vvjj=vvjj+1;
       }
      vvintop->left=vvintop->p->i1; 
      vvintop->right=vvintop->p->i2;
      vvintop->dd=vvdd;
      vvintop->ed=vved;
      vvintop->k=vvintop->p->k;
      vvintop->ab2p=vvintop->p->p;

      vvintop=vvintop->nextPtr;

      }        
      vvhug1=vvhug1->nextPtr;

    }
/**********SEARCH ALL INTERVALS!!!*****************/
while(hug1!=NULL)
{
ii++;
allldu=allldu+interv[ii][2]*interv[ii][1];
/*Point to the subset of pairs informative for this interval*/
if(interv[ii][5]>0.99)
{
if(interv[ii][6]>0.99)
{
intop=hug1->pp;
u=0.;
bigk=0.;
one_interval(ii,intop,&u,&bigk,&lnl);
if(g_niter>0)
       {
       if(bigk<0.00001)bigk=0.00001;
       interv[ii][2]=interv[ii][2]+(u/bigk); 
       if(interv[ii][2]<0.)interv[ii][2]=0.;
       if(g_niter>50){if((interv[ii][2]<0.00000001)&&(stope==1))interv[ii][5]=0.;} /*stop iterating those at bound */ 
       
       LDU=interv[ii][2]*interv[ii][1];
       if(g_niter<=50){  if(LDU>2.)interv[ii][2]=2./interv[ii][1]; } 


       if(g_niter>50)
       { 
       if(LDU>3.){interv[ii][2]=3./interv[ii][1];interv[ii][5]=0.;interv[ii][6]=0.;} /*stop iterating at upper bound */
       } 

       interv[ii][3]=u; 
       interv[ii][4]=bigk; 
       }
}/*interv[ii][5] test*/
}/*interv[ii][6] test*/
   hug1=hug1->nextPtr;
} /*END OF INTERVAL LOOP*/
/**********END OF SEARCH ALL INTERVALS!!!*****************/

old=new;      
new=-2.*lnl;
g_niter++;
diff=old-new;
/*------------*/
/* print out every 100th iteration */
/*ippy=ippy+1;
if(ippy==100) { fprintf(output_f,"\nIter=%4d diff=%16.6f LDU=%16.8f",
g_niter,fabs(old-new),allldu); fflush(output_f); }
if(ippy==100)ippy=0;
*/
/*------------*/

/*KEEP THE BEST RESULT */
if(new<g_best) { ntimes=0; g_best=new; }

fflush(output_ft);

if(g_niter==25||g_niter==50||g_niter==100||g_niter==200||g_niter==400||g_niter==800||g_niter==1600||g_niter==3200||g_niter==6400)
{
ldflag=0;
update();
if(g_lnl<BESTLNL)keep_it();

fflush(output_f);
fprintf(output_f,"\n%4d %10.6f %10.6f %10.6f %4d %16.8f %12.5f", g_niter,global_e,global_l,global_m,g_nhole,g_sed,g_lnl);
fflush(output_f);
fflush(output_ft);

if(g_lnl<BESTLNL)
{
/*MAKE A COPY OF THE BEST MAP TO DATE*/
BESTLNL=g_lnl;
for(i=0;i<MAX_LOCI;i++)
{
for(j=1;j<7;j++)cinterv[i][j]=interv[i][j];
}
}

if(g_lnl>g_startlnl){printf("\nWARNING: POORER FIT THAN STARTING -2lnlk"); 
           fprintf(output_f," WARNING: POORER FIT THAN STARTING -2lnlk");} 
fflush(output_f);
}
if(g_niter<100)goto topp;
if(g_niter<maxit){if(fabs(old-new)>conv)goto topp;}

/*FINAL UPDATE OF M */
ldflag=0;
update();
keep_it();
fprintf(output_f,"\n%4d %10.6f %10.6f %10.6f %4d %16.8f %12.5f\n\n", g_niter,global_e,global_l,global_m,g_nhole,g_sed,g_lnl);
if(g_lnl<BESTLNL)
{
/*MAKE A COPY OF THE BEST MAP TO DATE*/
BESTLNL=g_lnl;
for(i=0;i<MAX_LOCI;i++)
{
for(j=1;j<7;j++)cinterv[i][j]=interv[i][j];
}
}

/******RESTORE BEST MAP AT THIS POINT */
for(i=0;i<MAX_LOCI;i++)
{
for(j=1;j<7;j++)interv[i][j]=cinterv[i][j];
}
update();

if((global_l>=(g_pred-0.01))&&(global_l<=(g_pred+0.01))){g_df=g_nki-g_nloci;}
else{g_df=g_nki-g_nloci-1;}
g_V=g_lnl/g_df;
fprintf(output_f,"\n\nN(number of pairs)=%15d  m(number of SNPs)=%15d df=%14.1f V(error variance)=%14.5f\n\n", g_nki,g_nloci,g_df,g_V);
final();
end1=time(NULL);
xxx=end1-start1;
xxx=xxx/60.;
fprintf(output_f,"\n\n\n******ELAPSED TIME IN MINUTES FOR THIS RUN:  %f ******\n\n\n",xxx);

fin:fclose(output_f);
}
