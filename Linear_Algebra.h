#ifndef _LINEAR_ALGEBRAHEADER
#define _LINEAR_ALGEBRAHEADER

METHOD_STATE   State;
static int check = 0;

//parameter initial value
static double xsi[MAX_ORDER] = { 20000, 18700, 8900, 17900 };
static double ysi[MAX_ORDER] = { 19400, 1800, 17900, 21000 };
static double zsi[MAX_ORDER] = { 19740, 18500, 20000, 18500 };
static double ri[MAX_ORDER] = { 19992.498593, 18992.103621, 19463.491788, 19383.574838 };

//Newton Raphon's Method
static double initNR[MAX_ORDER][TRUE] = { 16273.0, 10924.0, 2012.0, 2.4263 * 0.00000001 };
static double preNR[MAX_ORDER][TRUE] = { 0.0 };
static double curNR[MAX_ORDER][TRUE] = { 0.0, };

static double Jacobian[MAX_ORDER][MAX_ORDER] = { 0, };
static double inv_Jacobian[MAX_ORDER][MAX_ORDER] = { 0, };
static double output[MAX_ORDER][TRUE] = { 0, };


void Find_GPS                    (IN int _method);

double* matrix_Jacobian          (IN double mat1[MAX_ORDER][TRUE]);

double Det_3                     (IN double Det[MAX_ORDER][MAX_ORDER]);

double Det_4                     (IN double Det[MAX_ORDER][MAX_ORDER]);

double* matrix_inverse_4         (IN double inv[MAX_ORDER][MAX_ORDER]);

double* matrix_mul               (IN double mat1[][MAX_ORDER], IN double mat2[][TRUE]);

double* matrix_substract         (IN double mat1[][TRUE], IN double mat2[][TRUE]);

double* ans_function             (IN double mat1[MAX_ORDER][TRUE]);

double RootFinding_NEWTONRAPHON  (IN double initNR[MAX_ORDER][TRUE]);
#endif

