#ifndef _OPTIMIZATION_
#define _OPTIMIZATION_

#include "basicHeader.h"
#include "UserDefine.h"
#include "Optimization.h"

double FindStepsize(double input[no_hw4]) {
	double temp[no_hw4] = { 0, };
	double den = 0, num = 0;

	num = grad[0] * grad[0] + grad[1] * grad[1];

	for (int i = 0; i < no_hw4; i++) {
		for (int j = 0; j < no_hw4; j++) {
			temp[i] = temp[i] + grad[j] * A[i][j];
		}
		den = den + temp[i] * grad[i];
	}

	
	double output = num / den;

	return output;
}


void Gradient_Descent() {
	
	int num = 0; 
	printf("\n\nGradient Method.......\n");

	Gradient(init_hw4);
	
	pre_x[0] = init_hw4[0];
	pre_x[1] = init_hw4[1];

	while (NORM2(grad[0], grad[1]) >= TOL_ERROR) {
		Next_GD(pre_x, tk);
		Gradient(pre_x);
		output_grad[num][0] = ans_x[0];
		output_grad[num][1] = ans_x[1];
		printf("%2d %3.2f %3.2f\n", num+1, output_grad[num][0], output_grad[num][1]);
		num++;
	}
	
}


void Steepest_Gradient_Descent() {

	int num = 0;
	printf("\n\nSteepest Method.......\n");

	Gradient(init_hw4);
	tk = FindStepsize(grad);

	pre_x[0] = init_hw4[0];
	pre_x[1] = init_hw4[1];

	while (NORM2(grad[0], grad[1]) >= TOL_ERROR) {
		Next_GD(pre_x, tk);
		Gradient(pre_x);
		tk = FindStepsize(grad);
		output_steep[num][0] = ans_x[0];
		output_steep[num][1] = ans_x[1];
		printf("%2d %3.2f %3.2f\n", num + 1, output_steep[num][0], output_steep[num][1]);
		num++;
	}

}


void Gradient(double input[no_hw4]){
	double x = input[0];
	double y = input[1];

	grad[0] = 3 * x + y - 4;
	grad[1] = x + 2 * y + 2;
}


double NORM2(double x1, double x2) {
	double temp = pow(x1, 2) + pow(x2, 2);
	return sqrt(temp);
}


void Next_GD(double input[no_hw4], double tk) {
	ans_x[0] = input[0] - tk * grad[0];
	ans_x[1] = input[1] - tk * grad[1];

	pre_x[0] = ans_x[0];
	pre_x[1] = ans_x[1];
}

#endif
