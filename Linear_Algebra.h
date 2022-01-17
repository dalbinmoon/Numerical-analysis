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

//HW3-1 parameter
static double H[11][2] = { {0, 1}, {10, 1}, {20, 1}, {30, 1}, {40, 1}, {50, 1}, {60, 1}, {70, 1}, {80, 1}, {90, 1}, {100, 1} };
static float P[11] = { 0.94, 0.96, 1.00, 1.05, 1.07, 1.09, 1.14, 1.17, 1.21, 1.24, 1.28 };

//HW3-2 parameter
static double R = 5000000;
static double H2[15][2] = { 0, };
static float V[15] = { 9.7, 8.1, 6.6, 5.1, 4.4, 3.7, 2.8, 2.4, 2.0, 1.6, 1.4, 1.1, 0.85, 0.69, 0.6 };

void Find_GPS                    (IN int _method);

double* matrix_Jacobian          (IN double mat1[MAX_ORDER][TRUE]);

float Det_3                     (IN float Det[READ_COLUMN][READ_COLUMN]);

double Det_4                     (IN double Det[MAX_ORDER][MAX_ORDER]);

double* matrix_inverse_4         (IN double inv[MAX_ORDER][MAX_ORDER]);

void matrix_inverse_3            (IN float inv[READ_COLUMN][READ_COLUMN], OUT float out[READ_COLUMN][READ_COLUMN]);

void matrix_mul                  (IN float mat1[READ_COLUMN][READ_COLUMN], IN float mat2[READ_COLUMN][READ_COLUMN], OUT float out[READ_COLUMN][READ_COLUMN]);

void matrix_add                   (IN float mat1[READ_COLUMN][READ_COLUMN], IN float mat2[READ_COLUMN][READ_COLUMN], OUT float out[READ_COLUMN][READ_COLUMN]);

double* matrix_substract         (IN double mat1[][TRUE], IN double mat2[][TRUE]);

double* ans_function             (IN double mat1[MAX_ORDER][TRUE]);

double RootFinding_NEWTONRAPHON  (IN double initNR[MAX_ORDER][TRUE]);

void transpose_matrix(IN float input[READ_COLUMN][READ_COLUMN], OUT float result[READ_COLUMN][READ_COLUMN]);


//HW3_1
void HW3_1();
void transpose_matrix_3_1    (IN double input[11][2], double result[2][11]);
void matrix_mul_3_1          (IN double input[2][11],  IN double input2[11][2], double output[2][2]);
void matrix_inverse_2        (IN double inv[2][2], OUT double output[2][2]);

//HW3_2
void HW3_2();
void transpose_matrix_3_2    (IN double input[15][2], double result[2][15]);
void matrix_mul_3_2          (IN double input[2][15], IN double input2[15][2], double output[2][2]);


#endif

