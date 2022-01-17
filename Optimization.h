#ifndef _OPTIMIZATION
#define _OPTIMIZATION

#define no_hw4     (int) 2

//Initial value
static double init_hw4[no_hw4] = { -2, 1 };
static double A[no_hw4][no_hw4] = { { 3, 1 }, { 1, 2 } };
static double ans_x[no_hw4] = { 0, };
static double pre_x[no_hw4] = { 0, };
static double grad[no_hw4] = { 0, };
static double tk = 0.1;
static double output_grad[100][no_hw4] = { 0, };
static double output_steep[100][no_hw4] = { 0, };

void Gradient_Descent();

double NORM2(double x1, double x2);             

void Gradient(double input[no_hw4]);

void Next_GD(double input[no_hw4], double tk);

double FindStepsize(double input[no_hw4]);

void Steepest_Gradient_Descent();


#endif