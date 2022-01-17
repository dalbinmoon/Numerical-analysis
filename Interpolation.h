#ifndef _INTERPOLATION
#define _INTERPOLATION

INTERPOLATION INTERPOL;

//Hermite Interpolation
static double init_x[MAX_INTER] = { -2163.22447, -2165.29848, -2167.26210, -2169.41539, -2171.45828 };
//static double init_x[MAX_INTER] = { 1, 2, 4, 4, 5 };
static double init_y[MAX_INTER] = { 4415.78095, 4426.45691, 4437.11169, 4447.74523, 4458.35748 };
static double init_z[MAX_INTER] = { 4879.67927, 4869.10118, 4858.49973, 4847.87496, 4837.22693 };
static double init_vx[MAX_INTER] = { -1.03959, -1.03441, -1.02922, -1.02404, -1.01885 };
//static double init_vx[MAX_INTER] = { 2, 2,2,2,2 };
//static double time[MAX_INTER] = { 2, 2,2,2,2 };
static double init_vy[MAX_INTER] = { 5.3427, 5.33269, 5.32208, 5.31145, 5.30080 };
static double init_vz[MAX_INTER] = { -5.28319, -5.29489, -5.30656, -5.31820, -5.32983 };
static double time[MAX_INTER] = { 2123, 2125, 2127, 2129, 2131 };

static double init_pose[MAX_INTER] = { -1577.12682, 4207.54609, 4511.39410 };


double Lagrange_Coef(double coef[MAX_INTER], int in);

double Langrange(double intput[MAX_INTER], double input2[MAX_INTER]);

double PI(double input[MAX_INTER], int size, double in);

double G(double input[MAX_INTER], int size, double in);

double Hermite(double input[MAX_INTER], double prime[MAX_INTER], double in);

double Hermite_Prime(double input[MAX_INTER], double prime[MAX_INTER], double in);

void HW2_1();

void InterPolation(int INTERPOL);

#endif
