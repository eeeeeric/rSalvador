/* rSalvador: An R tool for the Luria-Delbruck fluctuation assay */
/* Qi Zheng, Department of epidemiology and Biostatistics */
/* Texas A&M School of Public Health */
/* Version 1.0: April 20, 2014 */
/* Version 1.2: June 2, 2015 */


/* double* pdfluria(double, double, int, double*); */

void pdfluria(double, double, int, double*); 

void pdfconv(double*,int, double*, int, double*);

void pmfLDPlat(double, double, double*, int, double*);

/* ------- Haldane -------- */

void pmfHald(int, double, int, int, double*);

void derivHald(int, double, int, int, double*, double*, double*);

void pAndP1Hald(int, double, int, int, double*, double*);

void simuHald(int, double, int*);

void simuKimmel(int, double, int*);

void simuKimm2(int, double, int*, int*);

void simuHald2(int, double, int*, int*);




/* ------ Wrappers ------------------ */

void luria_R_wrapper(double *m, double *phi, int *k, double *result) {

 pdfluria(*m, *phi, *k, result);

};



void pdfconv_R_wrapper(double* x, int* xlen, double* y, int* ylen, double* result) {

pdfconv(x, *xlen, y, *ylen, result);

};


void pmfLDPlat_R_wrapper(double* m, double* e, double* eta, int* etaLen, double* prob){

pmfLDPlat(*m, *e, eta, *etaLen, prob);

};

/* December 22, 2013, Haldane formulation */

void pmfHald_R_wrapper(int* gen, double* mu, int* n, int* N0, double* prob){

pmfHald(*gen, *mu, *n, *N0, prob);

};

void derivHald_R_wrapper(int* gen, double* mu, int* n, int* N0,
     double* prob, double* prob1, double* prob2){

derivHald(*gen, *mu, *n, *N0, prob, prob1, prob2);

};

void pAndP1Hald_R_wrapper(int* gen, double* mu, int* n, int* N0, double* prob, double* prob1){

pAndP1Hald(*gen, *mu, *n, *N0, prob, prob1);

};

void simuHald_R_wrapper(int* gen, double* mu, int* mut) {

simuHald(*gen, *mu, mut); 

};

void simuKimmel_R_wrapper(int* nGgen, double* mutRate, int* mutants) {

simuKimmel(*nGgen, *mutRate, mutants);

};


void simuKimm2_R_wrapper(int* nGgen, double* mutRate, int* mutants0, int* mutants) {

simuKimm2(*nGgen, *mutRate, mutants0, mutants);

};

/* Jan 1, 2014 */
void simuHald2_R_wrapper(int* gen, double* mu, int* mut0, int* mut) {

simuHald2(*gen, *mu, mut0, mut); 

};



/* ------------ taking care of fitness, May 15, 2015 --------------- */

void pmfMK(double, double, double*, int, double*);

void pmfMK_R_wrapper(double* m, double* w, double* beta, int* betaLen, double* prob){

pmfMK(*m, *w, beta, *betaLen, prob);

};

