/*MODIFIED SEPT 2001 - NOW TESTS BOUNDS AFTER 1 CYCLE. THEN IT TAKES BOUNDED PARAMETERS OUT
AND CONTINUES ITERATION */

#include "allass.h" 
#define NRANSI
#define ITMAX 50 
#define EPS 3.0e-8
#define TOLX (4*EPS)
#define STPMX 10.000

#define FREEALL free_dvector(xi,1,n);free_dvector(pnew,1,n); \
free_dmatrix(hessin,1,n,1,n);free_dvector(hdg,1,n);free_dvector(g,1,n); \
free_dvector(dg,1,n);

void dfpmin(double p[], int n, double gtol, int *iter, double *fret,
	double(*func)(double []), void (*dfunc)(double [], double []))
{
	void lnsrch(int n, double xold[], double fold, double g[], double p[], double x[],
		 double *f, double stpmax, int *check, double (*func)(double []));
	int check,i,itss,j;
	double xall,den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test;
	double **copy,*dg,*g,*hdg,**hessin,*pnew,*xi;
top:	dg=dvector(1,n);
	g=dvector(1,n);
	hdg=dvector(1,n);
	copy=dmatrix(1,n,1,n);
	hessin=dmatrix(1,n,1,n);
	pnew=dvector(1,n);
	xi=dvector(1,n);
	fp=(*func)(p);
	(*dfunc)(p,g);
        for (i=1;i<=n;i++) {
		for (j=1;j<=n;j++){ hessin[i][j]=0.0;copy[i][j]=0.0;}
		hessin[i][i]=1.0;
                copy[i][i]=1.0;
	 /*xi here are the starting U scores !! */	
                xi[i] = -g[i];
        /*p[i] are the parameter estimates!!! */        
                sum += p[i]*p[i];
	}

         stpmax=STPMX*FMAX(sqrt(sum),(double)n);
         for (itss=1;itss<=ITMAX;itss++) {
		*iter=itss;
                 lnsrch(n,p,fp,g,xi,pnew,fret,stpmax,&check,func);
/*fp is the value of -2lnL **************************************************/
                fp = *fret;
                for (i=1;i<=n;i++) {
/*xi[i] is now the difference between the old and new parameter estimates!! */			
                     xi[i]=pnew[i]-p[i];
			p[i]=pnew[i];

                }
		/*test=0.0;*/
		test=0.0;
		for (i=1;i<=n;i++) {
			temp=fabs(xi[i])/FMAX(fabs(p[i]),1.0);
                        if (temp > test) test=temp;
		                  }
	  	if (test < TOLX) {
                /**********************/ 
                xall=bound(n,p); 
                /*FREEALL*/
                n=g_n; 
                if(xall==1){
                            FREEALL
                            goto top; 
                            } 
                /**********************/ 
	        global_fin=1;	
                (*dfunc)(p,g);
                FREEALL
	         return;
		}
		for (i=1;i<=n;i++){ dg[i]=g[i];} /* dg is the old gradient (u score) - saved here !!! */
		(*dfunc)(p,g);
                test=0.0;
		den=FMAX(*fret,1.0);
		for (i=1;i<=n;i++) {
			temp=fabs(g[i])*FMAX(fabs(p[i]),1.0)/den;
			if (temp > test) test=temp;
		}
		if (test < gtol) {
                /**********************/ 
                xall=bound(n,p); 
                /*FREEALL*/
                n=g_n; 
                if(xall==1){
                            FREEALL
                            goto top; 
                            } 
	        /**********************/ 
                global_fin=1;	
                (*dfunc)(p,g);
                    FREEALL
                    return;
		}
		for (i=1;i<=n;i++) dg[i]=g[i]-dg[i];
		for (i=1;i<=n;i++) {
			hdg[i]=0.0;
			for (j=1;j<=n;j++){if(i==j)hessin[i][j]=fabs(hessin[i][j]);copy[i][j]=hessin[i][j]; hdg[i] += hessin[i][j]*dg[j];}
		}
		fac=fae=sumdg=sumxi=0.0;
		for (i=1;i<=n;i++) {
			fac += dg[i]*xi[i];
			fae += dg[i]*hdg[i];
			sumdg += SQR(dg[i]);
			sumxi += SQR(xi[i]);
		}
		if (fac*fac > EPS*sumdg*sumxi) {
			fac=1.0/fac;
			fad=1.0/fae;
			for (i=1;i<=n;i++) dg[i]=fac*xi[i]-fad*hdg[i];
                        for (i=1;i<=n;i++) {
				for (j=1;j<=n;j++) {
					hessin[i][j] += fac*xi[i]*xi[j]
					-fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
			                if(i==j)hessin[i][j]=fabs(hessin[i][j]);	
                                        copy[i][j]=hessin[i][j];          
                                                        }
			}
		
                }
		for (i=1;i<=n;i++) {
			xi[i]=0.0;
			for (j=1;j<=n;j++){if(i==j)hessin[i][j]=fabs(hessin[i][j]); copy[i][j]=hessin[i][j];xi[i] -= hessin[i][j]*g[j];}
		}
	}
	nrerror("too many iterations in dfpmin");
                /**********************/ 
                xall=bound(n,p); 
                /*FREEALL*/
                n=g_n; 
                if(xall==1){
                            FREEALL
                            goto top; 
                            } 
                /**********************/ 
	        global_fin=1;	
                (*dfunc)(p,g);
        FREEALL
}
#undef ITMAX
#undef EPS
#undef TOLX
#undef STPMX
#undef FREEALL
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software ;1[V3. */
