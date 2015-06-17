/* rSalvador: An R tool for the Luria-Delbruck fluctuation assay */
/* Qi Zheng, Department of epidemiology and Biostatistics */
/* Texas A&M School of Public Health */
/* Version 1.0: April 20, 2014 */
/* Version 1.2: June 2, 2015 */


/* July 9, 2013: to rewrite Salvador in R */

#include <math.h>
#include <malloc.h>

#include <stdio.h>  /* added Dec 26, 2013, EOS needs this! */

#include <stdlib.h>  /* added for the exit() function, December 22, 2013 */

#ifdef _WIN32   /* April 20, 2014, Eric */
    #define RAND() rand()
#else
    #define RAND() random()
#endif

void pdfluria(double m, double phi, int k, double* problist){
double cumuphi;
int n,j;

for(n=1;n<=k;n++) problist[n]=0;
problist[0]=exp(-m); 
for(n=1;n<=k;n++) {
  
   cumuphi=1.0;
   for(j=1;j<=n;j++){
      
      problist[n]=problist[n]+cumuphi*(1-j*phi/(j+1))*problist[n-j];
      cumuphi=cumuphi*phi;}
   problist[n]=problist[n]*m/n;};
}


/* March 6, 2001: to convolute 2 lists of PDFs; July 9, 2013 */
  
void pdfconv(double* x, int xlen, double* y, int ylen, double* z){ 
int i, k; 

for (i=0;i<=xlen-1;i++) z[i]=0;
for(k=1;k<=xlen;k++){
for(i=1;i<=k;i++){

z[k-1]+=x[i-1]*y[k-i];};
};

} 

/* ------------ plating efficiency, July 13, 2013 -------------- */


void pmfLDPlat(double m, double e, double* eta, int etaLen, double* prob){
double odds;
int n,j;

odds=e/(1-e);
for(n=1;n<=etaLen-1;n++) prob[n]=0;
prob[0]=exp(m*odds*log(e)); 
for(n=1;n<=etaLen-1;n++) {

   for(j=1;j<=n;j++){

      prob[n]=prob[n]+j*eta[j]*prob[n-j]; }
   prob[n]=prob[n]*m*odds/n;};


}


/* Taken from SALVADOR 2.3, December 22, 2013 */
/* binomial distribution */
double binProb(int n, double p, int x){
double pp,f;
int xx,j0, j1, j2;
if ( (x>n) | (x<0) ) return 0.0; /* added Jan 17, 2008 */
if (2*x>n) {pp=1-p; xx=n-x; }
 else {pp=p;xx=x;}
j0=j1=j2=0;
f=1.0;
while ((j0<xx) | (j1<xx) | (j2<n-xx))
{ if ((j0<xx) && (f<1) )
 { j0++;
   f*=(double)(n-xx+j0)/(double)j0;
}
 else
 { if (j1<xx) { j1++; f*=pp; }
        else { j2++; f*= 1-pp;}
}
}
return(f);
}


/* the Akj sequence as described in Mathematical Biosciences */
double* getAkj(int g, int k, double mu, int N0)
{
int j, U, Ng, minJ;
double* Akj;
Ng=N0*pow(2,g);
if (k>Ng) minJ=k-Ng; else minJ=0;
U=floor(k/2);
Akj=(double*) malloc(sizeof(double)*(U+1));
if (Akj==NULL) {fprintf(stderr,"insufficient memory!"); exit(1);}
for (j=0;j<=U;j++) Akj[j]=0;

if (k%2 == 0) Akj[U]=pow(1-mu, Ng-k/2);
else Akj[U]=(Ng-(k-1)/2)*mu*pow(1-mu,Ng-(k+1)/2);

for(j=U-1;j>=minJ;j--){
/* this re-arrangement is necessary to ensure numerical stability, 1-14-08 */
 Akj[j]=(mu*(Ng-j))*(mu*(Ng-k+j+1))/((k-2*j)*(k-2*j-1)*(1-mu))*Akj[j+1];} 
return Akj;
}


/* Adapted from SALVADOR 2.3, December 22, 2013 */

void pmfHald(int gen, double mu, int n, int N0, double* prob)
{
int i,j,k,g,minJ, nn[gen+1];
double p[n+1],  *A;

for (i=0;i<=n;i++) prob[i]=0;

for (i=0;i<=n;i++) p[i]=0;

nn[gen]=n;
for (i=gen-1;i>=1;i--)  nn[i]=floor(nn[i+1]/2);

for (i=0; i<=nn[1]; i++) prob[i]=binProb(N0,mu,i);

for (g=2; g<=gen; g++) {
   for (i=0;i<=nn[g];i++) {p[i]=prob[i]; prob[i]=0; }

   for (k=0;k<=nn[g];k++) {
      
      A=getAkj(g-1,k,mu,N0);
      if(k>N0*pow(2,g-1)) minJ=k-N0*pow(2,g-1); else minJ=0;
      for(j=minJ;j<=floor(k/2);j++) prob[k]+=A[j]*p[j];
      free(A); 
    }
   }

}

/* newton for Haldane, adapted from SALVADOR 2.3, December 23, 2013 */

/* -------------------- derivatives for Haldane -------------------- */

double a(int g, int k, int j, double mu, int N0){
return (k-2*j-(N0*pow(2,g)-j)*mu)/(mu*(1-mu));
}

double ap(int g, int k, int j, double mu, int N0){
return (2*j-k+2*(k-2*j)*mu-(N0*pow(2,g)-j)*mu*mu)/(mu*mu*(1-mu)*(1-mu));
}


void derivHald(int gen, double mu, int n, int N0, double* prob, double* prob1, double* prob2)
{
int i,j,k,g,minJ, *nn, Ng;
double  *prob0, *prob10, *prob20, *A;

for (i=0;i<=n;i++) prob[i]=0;

prob0=(double*) malloc(sizeof(double)*(n+1)); 
if (prob0==NULL) {fprintf(stderr,"insufficient memory!"); exit(1);}
for (i=0;i<=n;i++) prob0[i]=0;

for (i=0;i<=n;i++) prob1[i]=0;

prob10=(double*) malloc(sizeof(double)*(n+1)); 
if (prob10==NULL) {fprintf(stderr,"insufficient memory!"); exit(1);}
for (i=0;i<=n;i++) prob10[i]=0;

for (i=0;i<=n;i++) prob2[i]=0;

prob20=(double*) malloc(sizeof(double)*(n+1)); 
if (prob20==NULL) {fprintf(stderr,"insufficient memory!"); exit(1);}
for (i=0;i<=n;i++) prob20[i]=0;

nn=(int*) malloc(sizeof(int)*(gen+1)); /* nn[i] is l[i] in the paper */
if (nn==NULL) {fprintf(stderr,"insufficient memory!"); exit(1);}
/* computing the nn[i]=l[i], nn[0] is not used */
nn[gen]=n;
for (i=gen-1;i>=1;i--)  nn[i]=floor(nn[i+1]/2);

/* initialization */
for (j=0; j<=nn[1]; j++) prob[j]=binProb(N0,mu,j);

for (j=0; j<=nn[1]; j++) prob1[j]=(j-N0*mu)/(mu*(1-mu))*prob[j];

for (j=0; j<=nn[1]; j++) 
  prob2[j]=(j*(j-1)-2*j*(N0-1)*mu+N0*mu*(N0-1)*mu)*prob[j]/(mu*mu*(1-mu)*(1-mu));

/* iteration from g=2 */

for (g=2; g<=gen; g++) {
   for (i=0;i<=nn[g];i++) {prob0[i]=prob[i]; prob[i]=0; 
                           prob10[i]=prob1[i]; prob1[i]=0; 
                           prob20[i]=prob2[i]; prob2[i]=0; }

   for (k=0;k<=nn[g];k++) {
      A=getAkj(g-1,k,mu,N0);
      if(k>N0*pow(2,g-1)) minJ=k-N0*pow(2,g-1); else minJ=0;

      for(j=minJ;j<=floor(k/2);j++){
        
        prob[k]+=A[j]*prob0[j];
        prob1[k]+=( A[j]*a(g-1,k,j,mu,N0)*prob0[j] + A[j]*prob10[j] );
        prob2[k]+=( A[j]*(a(g-1,k,j,mu,N0)*a(g-1,k,j,mu,N0)+ap(g-1,k,j,mu,N0))*\
                  prob0[j]+2*A[j]*a(g-1,k,j,mu,N0)*prob10[j]+A[j]*prob20[j] );
        }
       free(A);
    }
   }

free(nn);

free(prob0);

free(prob10);

free(prob20);

}


/*  derivHald without the second derivative, for CI, to save time, Dec 24, 2013 */

void pAndP1Hald(int gen, double mu, int n, int N0, double* prob, double* prob1)
{
int i,j,k,g,minJ, *nn, Ng;
double  *prob0, *prob10, *A;

for (i=0;i<=n;i++) prob[i]=0;

prob0=(double*) malloc(sizeof(double)*(n+1)); 
if (prob0==NULL) {fprintf(stderr,"insufficient memory!"); exit(1);}
for (i=0;i<=n;i++) prob0[i]=0;

for (i=0;i<=n;i++) prob1[i]=0;

prob10=(double*) malloc(sizeof(double)*(n+1)); 
if (prob10==NULL) {fprintf(stderr,"insufficient memory!"); exit(1);}
for (i=0;i<=n;i++) prob10[i]=0;

nn=(int*) malloc(sizeof(int)*(gen+1)); /* nn[i] is l[i] in the paper */
if (nn==NULL) {fprintf(stderr,"insufficient memory!"); exit(1);}
/* computing the nn[i]=l[i], nn[0] is not used */
nn[gen]=n;
for (i=gen-1;i>=1;i--)  nn[i]=floor(nn[i+1]/2);

/* initialization */
for (j=0; j<=nn[1]; j++) prob[j]=binProb(N0,mu,j);

for (j=0; j<=nn[1]; j++) prob1[j]=(j-N0*mu)/(mu*(1-mu))*prob[j];

/* iteration from g=2 */

for (g=2; g<=gen; g++) {
   for (i=0;i<=nn[g];i++) {prob0[i]=prob[i]; prob[i]=0; 
                           prob10[i]=prob1[i]; prob1[i]=0;}

   for (k=0;k<=nn[g];k++) {
      A=getAkj(g-1,k,mu,N0);
      if(k>N0*pow(2,g-1)) minJ=k-N0*pow(2,g-1); else minJ=0;

      for(j=minJ;j<=floor(k/2);j++){
        
        prob[k]+=A[j]*prob0[j];
        prob1[k]+=( A[j]*a(g-1,k,j,mu,N0)*prob0[j] + A[j]*prob10[j] );
        }
       free(A);
    }
   }

free(nn);

free(prob0);

free(prob10);

}



/* --------------- Haldane simulation December 23, 2013 -------------- */

void simuHald(int gen, double mu, int* mut) {
int wild=1, newWild, newMut, i, j;

*mut=0;

for (i=1; i<=gen; i++) {
  newMut=newWild=0;
  for (j=1;j<=wild; j++) {
    if ( RAND()/( (double)RAND_MAX+1) < mu) newMut++;  /* from rand() to random(),Dec 23, 2013 */
    else newWild++;
  }
    mut[0]=2* (*mut)+newMut;
    wild=wild+newWild;
}
}

/* ---------------- Kimmel model as described by Sinden, new to rSalvador, Dec 24, 2013 */

void simuKimmel(int nGen, double mutRate, int* mutants) {

int  cell=0, total, gen, j;

total=1;

/* mutants[0]=-1; */

for(gen=0; gen < nGen; gen++){
cell=total;
total=0;

for(j=0; j<(2*cell); j++){
double mu=RAND()/( (double)RAND_MAX+1);

if(mu > mutRate){

total++;
}
} 
}

*mutants=pow(2,nGen)-total;
}

/* ----------------- Dec 30, 2013, Kimmel model with the last two generations exported ---- */


void simuKimm2(int nGen, double mutRate, int* mutants0, int* mutants) {

int  cell=0, total, gen, j;

total=1;

/* mutants[0]=-1; */

for(gen=0; gen < nGen; gen++){
cell=total;
total=0;

for(j=0; j<(2*cell); j++){
double mu=RAND()/( (double)RAND_MAX+1);

if(mu > mutRate){

total++;
}
} 
}

*mutants0=pow(2,nGen-1)-cell;
*mutants=pow(2,nGen)-total;
}



/* --------------- Haldane simulation for the last two generations --------- */
/* to report the last two generations, mut0 is for the next to last generation, Dec 31, 2013 */

void simuHald2(int gen, double mu, int* mut0, int* mut) {
int wild=1, newWild, newMut, i, j;

*mut0=0;
*mut=0;

for (i=1; i<=gen; i++) {
  newMut=newWild=0;
  for (j=1;j<=wild; j++) {
    if ( RAND()/( (double)RAND_MAX+1) < mu) newMut++;  /* from rand() to random(),Dec 23, 2013 */
    else newWild++;
  }
    mut0[0]=*mut; /* to keep the previous generation, Dec 31, 2013 */
    mut[0]=2* (*mut)+newMut;
    wild=wild+newWild;
}
}


/* ------------ fitness, May 15, 2015 -------------- */

void pmfMK(double m, double w, double* B, int betaLen, double* prob){

int n,j;

for(n=1;n<=betaLen-1;n++) prob[n]=0;
prob[0]=exp(-m); 
for(n=1;n<=betaLen-1;n++) {

   for(j=1;j<=n;j++){

      prob[n]=prob[n]+j*B[j-1]*prob[n-j]; }

   prob[n]=prob[n]*m/n/w;};

}




