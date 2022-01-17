#ifndef _LINEAR_ALGEBRA
#define _LINEAR_ALGEBRA

#include "basicHeader.h"
#include "UserDefine.h"
#include "Linear_Algebra.h"
#include "Sol_NonLinearEqns.h"

//extern void Error_Message(IN unsigned int _Errortype);
void Find_GPS(IN int _method) {
	printf("\nFinding GPS.....\n");
	RootFinding(NEWTONRAPHON);
	printf("Solution = { %10lf %10lf %10lf %10lf }\n", curNR[0][0], curNR[1][0], curNR[2][0], curNR[3][0]);
}

double* matrix_Jacobian(IN double mat1[MAX_ORDER][TRUE]) {
	
	for (int i = 0 ; i < MAX_ORDER ; i++){
		for (int j = 0; j < MAX_ORDER; j++) {
			double jacob_den = sqrt(pow(mat1[0][0] - xsi[j], 2.0) + pow(mat1[0][1] - ysi[j], 2.0) + pow(mat1[0][2] - zsi[j], 2.0));
			if (j == 0) Jacobian[i][j] = (mat1[0][0] - xsi[j]) / jacob_den;
			if (j == 1) Jacobian[i][j] = (mat1[0][1] - ysi[j]) / jacob_den;
			if (j == 2) Jacobian[i][j] = (mat1[0][2] - zsi[j]) / jacob_den;
			if (j == 3) Jacobian[i][j] = -SPEED_LIGHT;
		}
	}

	return Jacobian;
}

double Det_3(IN double Det[MAX_ORDER][MAX_ORDER]) {
	double det =
		Det[0][0] * Det[1][1] * Det[2][2] + Det[1][0] * Det[2][1] * Det[0][2] + Det[2][0] * Det[0][1] * Det[1][2] - \
		Det[0][0] * Det[2][1] * Det[1][2] - Det[2][0] * Det[1][1] * Det[0][2] - Det[1][0] * Det[0][1] * Det[2][2];

	if (det == 0.0) { State = NO_DET; Error_Message(State); return 0.0; }
	return det;
}

double Det_4(IN double Det[MAX_ORDER][MAX_ORDER]) {
	double det = Det[0][0] * Det[1][1] * Det[2][2] * Det[3][3] + Det[0][0] * Det[1][2] * Det[2][3] * Det[3][1] + Det[0][0] * Det[1][3] * Det[2][1] * Det[3][2] + Det[0][1] * Det[1][0] * Det[2][3] * Det[3][2] + Det[0][1] * Det[1][2] * Det[2][0] * Det[3][3] + Det[0][1] * Det[1][3] * Det[2][2] * Det[3][0] + Det[0][2] * Det[1][0] * Det[2][1] * Det[3][3] + Det[0][2] * Det[1][1] * Det[2][3] * Det[3][0] + Det[0][3] * Det[1][3] * Det[2][0] * Det[3][1] + Det[0][3] * Det[1][0] * Det[2][2] * Det[3][1] + Det[0][3] * Det[1][1] * Det[2][0] * Det[3][2] + Det[0][3] * Det[1][2] * Det[2][1] * Det[3][0] - Det[0][0] * Det[1][1] * Det[2][3] * Det[3][2] - Det[0][0] * Det[1][2] * Det[2][1] * Det[3][3] - Det[0][0] * Det[1][3] * Det[2][2] * Det[3][1] - Det[0][1] * Det[1][0] * Det[2][2] * Det[3][3] - Det[0][1] * Det[1][2] * Det[2][3] * Det[3][0] - Det[0][1] * Det[1][3] * Det[2][0] * Det[3][2] - Det[0][2] * Det[1][0] * Det[2][3] * Det[3][1] - Det[0][2] * Det[1][1] * Det[2][0] * Det[3][3] - Det[0][2] * Det[1][3] * Det[2][1] * Det[3][0] - Det[0][3] * Det[1][0] * Det[2][1] * Det[3][2] - Det[0][3] * Det[1][1] * Det[2][2] * Det[3][0] - Det[0][3] * Det[1][2] * Det[2][0] * Det[3][1];

	if (det == 0.0) { State = NO_DET; Error_Message(State); return 0.0; }
	return det;
}

double* matrix_inverse_4(IN double inv[MAX_ORDER][MAX_ORDER]) {
	double parameter = 1 / Det_4(Jacobian);

	inv_Jacobian[0][0] = parameter * (Jacobian[1][1] * Jacobian[2][2] * Jacobian[3][3] + Jacobian[1][2] * Jacobian[2][3] * Jacobian[3][1] + Jacobian[1][3] * Jacobian[2][1] * Jacobian[3][2] - Jacobian[1][1] * Jacobian[2][3] * Jacobian[3][2] - Jacobian[1][2] * Jacobian[2][1] * Jacobian[3][3] - Jacobian[1][3] * Jacobian[2][2] * Jacobian[3][1]);
	inv_Jacobian[0][1] = parameter * (Jacobian[0][1] * Jacobian[2][3] * Jacobian[3][2] + Jacobian[0][2] * Jacobian[2][1] * Jacobian[3][3] + Jacobian[0][3] * Jacobian[2][2] * Jacobian[3][1] - Jacobian[0][1] * Jacobian[2][2] * Jacobian[3][3] - Jacobian[0][2] * Jacobian[2][3] * Jacobian[3][1] - Jacobian[0][3] * Jacobian[2][1] * Jacobian[3][2]);
	inv_Jacobian[0][2] = parameter * (Jacobian[0][1] * Jacobian[1][2] * Jacobian[3][3] + Jacobian[0][2] * Jacobian[1][3] * Jacobian[3][1] + Jacobian[0][3] * Jacobian[1][1] * Jacobian[3][2] - Jacobian[0][1] * Jacobian[1][3] * Jacobian[3][2] - Jacobian[0][2] * Jacobian[1][1] * Jacobian[3][3] - Jacobian[0][3] * Jacobian[1][2] * Jacobian[3][1]);
	inv_Jacobian[0][3] = parameter * (Jacobian[0][1] * Jacobian[1][3] * Jacobian[2][2] + Jacobian[0][2] * Jacobian[1][1] * Jacobian[2][3] + Jacobian[0][3] * Jacobian[1][2] * Jacobian[2][1] - Jacobian[0][1] * Jacobian[1][2] * Jacobian[2][3] - Jacobian[0][2] * Jacobian[1][3] * Jacobian[2][1] - Jacobian[0][3] * Jacobian[1][1] * Jacobian[2][2]);
	
	inv_Jacobian[1][0] = parameter * (Jacobian[1][0] * Jacobian[2][3] * Jacobian[3][2] + Jacobian[1][2] * Jacobian[2][0] * Jacobian[3][3] + Jacobian[1][3] * Jacobian[2][2] * Jacobian[3][0] - Jacobian[1][0] * Jacobian[2][2] * Jacobian[3][3] - Jacobian[1][2] * Jacobian[2][3] * Jacobian[3][0] - Jacobian[1][3] * Jacobian[2][0] * Jacobian[3][2]);
	inv_Jacobian[1][1] = parameter * (Jacobian[0][0] * Jacobian[2][2] * Jacobian[3][3] + Jacobian[0][2] * Jacobian[2][3] * Jacobian[3][1] + Jacobian[0][3] * Jacobian[2][0] * Jacobian[3][2] - Jacobian[0][0] * Jacobian[2][3] * Jacobian[3][2] - Jacobian[0][2] * Jacobian[2][0] * Jacobian[3][3] - Jacobian[0][3] * Jacobian[2][2] * Jacobian[3][0]);
	inv_Jacobian[1][2] = parameter * (Jacobian[0][0] * Jacobian[1][3] * Jacobian[3][2] + Jacobian[0][2] * Jacobian[1][0] * Jacobian[3][3] + Jacobian[0][3] * Jacobian[1][2] * Jacobian[3][0] - Jacobian[0][0] * Jacobian[1][2] * Jacobian[3][3] - Jacobian[0][2] * Jacobian[1][3] * Jacobian[3][0] - Jacobian[0][3] * Jacobian[1][0] * Jacobian[3][2]);
	inv_Jacobian[1][3] = parameter * (Jacobian[0][0] * Jacobian[1][2] * Jacobian[2][3] + Jacobian[0][2] * Jacobian[1][3] * Jacobian[2][1] + Jacobian[0][3] * Jacobian[1][0] * Jacobian[2][2] - Jacobian[0][0] * Jacobian[1][3] * Jacobian[2][2] - Jacobian[0][2] * Jacobian[1][0] * Jacobian[2][3] - Jacobian[0][3] * Jacobian[1][2] * Jacobian[2][0]);
	
	inv_Jacobian[2][0] = parameter * (Jacobian[1][0] * Jacobian[2][1] * Jacobian[3][3] + Jacobian[1][1] * Jacobian[2][3] * Jacobian[3][0] + Jacobian[1][3] * Jacobian[2][0] * Jacobian[3][1] - Jacobian[1][0] * Jacobian[2][3] * Jacobian[3][1] - Jacobian[1][1] * Jacobian[2][0] * Jacobian[3][3] - Jacobian[1][3] * Jacobian[2][2] * Jacobian[3][0]);
	inv_Jacobian[2][1] = parameter * (Jacobian[0][0] * Jacobian[2][3] * Jacobian[3][1] + Jacobian[0][1] * Jacobian[2][0] * Jacobian[3][3] + Jacobian[0][3] * Jacobian[2][1] * Jacobian[3][0] - Jacobian[0][0] * Jacobian[2][1] * Jacobian[3][3] - Jacobian[0][1] * Jacobian[2][3] * Jacobian[3][0] - Jacobian[0][3] * Jacobian[2][0] * Jacobian[3][1]);
	inv_Jacobian[2][2] = parameter * (Jacobian[0][0] * Jacobian[1][1] * Jacobian[3][3] + Jacobian[0][1] * Jacobian[1][3] * Jacobian[3][0] + Jacobian[0][3] * Jacobian[1][0] * Jacobian[3][1] - Jacobian[0][0] * Jacobian[1][3] * Jacobian[3][1] - Jacobian[0][1] * Jacobian[1][0] * Jacobian[3][3] - Jacobian[0][3] * Jacobian[1][2] * Jacobian[3][0]);
	inv_Jacobian[2][3] = parameter * (Jacobian[0][0] * Jacobian[1][3] * Jacobian[2][1] + Jacobian[0][1] * Jacobian[1][0] * Jacobian[2][3] + Jacobian[0][3] * Jacobian[1][1] * Jacobian[2][0] - Jacobian[0][0] * Jacobian[1][1] * Jacobian[2][3] - Jacobian[0][1] * Jacobian[1][3] * Jacobian[2][0] - Jacobian[0][3] * Jacobian[1][0] * Jacobian[2][1]);
	
	inv_Jacobian[3][0] = parameter * (Jacobian[1][0] * Jacobian[2][2] * Jacobian[3][1] + Jacobian[1][1] * Jacobian[2][0] * Jacobian[3][2] + Jacobian[1][2] * Jacobian[2][1] * Jacobian[3][0] - Jacobian[1][0] * Jacobian[2][1] * Jacobian[3][2] - Jacobian[1][1] * Jacobian[2][2] * Jacobian[3][0] - Jacobian[1][2] * Jacobian[2][0] * Jacobian[3][1]);
	inv_Jacobian[3][1] = parameter * (Jacobian[0][0] * Jacobian[2][1] * Jacobian[3][2] + Jacobian[0][1] * Jacobian[2][2] * Jacobian[3][0] + Jacobian[0][2] * Jacobian[2][0] * Jacobian[3][1] - Jacobian[0][0] * Jacobian[2][2] * Jacobian[3][1] - Jacobian[0][1] * Jacobian[2][0] * Jacobian[3][2] - Jacobian[0][2] * Jacobian[2][1] * Jacobian[3][0]);
	inv_Jacobian[3][2] = parameter * (Jacobian[0][0] * Jacobian[1][2] * Jacobian[3][1] + Jacobian[0][1] * Jacobian[1][0] * Jacobian[3][2] + Jacobian[0][2] * Jacobian[1][1] * Jacobian[3][0] - Jacobian[0][0] * Jacobian[1][1] * Jacobian[3][2] - Jacobian[0][1] * Jacobian[1][2] * Jacobian[3][0] - Jacobian[0][2] * Jacobian[1][0] * Jacobian[3][1]);
	inv_Jacobian[3][3] = parameter * (Jacobian[0][0] * Jacobian[1][1] * Jacobian[2][2] + Jacobian[0][1] * Jacobian[2][2] * Jacobian[2][0] + Jacobian[0][2] * Jacobian[1][0] * Jacobian[2][1] - Jacobian[0][0] * Jacobian[1][2] * Jacobian[2][1] - Jacobian[0][1] * Jacobian[1][0] * Jacobian[2][2] - Jacobian[0][2] * Jacobian[1][1] * Jacobian[2][0]);

	return inv_Jacobian;
}

double* matrix_mul(IN IN double mat1[][MAX_ORDER], IN double mat2[][TRUE]) {
	int output_col = sizeof(mat2[0]) / sizeof(double);  
	int output_row = sizeof(mat1) / sizeof(mat1[0]);
	//double** output = malloc(sizeof(double*)*output_row);
	//for (int i = 0; i < output_row; i++) {
	//	output[i] = malloc(sizeof(double) * output_col);
	//}

	for (int i = 0; i < output_row; i++) {
		for (int j = 0; j < output_col; j++) {
			output[i][j] = output[i][j] + mat1[i][j] * mat2[j][i];
		}
	}

	return output;
}

double* matrix_substract(IN double mat1[][TRUE], IN double mat2[][TRUE]) {
	if (check == 0) {
		for (int i = 0; i < MAX_ORDER; i++) {
			preNR[i][0] = mat1[i][0] + mat2[i][0];
		}
		check = 1;
	}

	if (check == 1) {
		for (int i = 0; i < MAX_ORDER; i++) {
			curNR[i][0] = mat1[i][0] + mat2[i][0];
		}
	}



}

double* ans_function(IN double mat1[MAX_ORDER][TRUE]) {
	double output_ans[MAX_ORDER][TRUE] = { 0, };

	for (int i = 0; i < MAX_ORDER; i++) {
		double temp = pow(output[0][0] - ri[i], 2.0) + pow(output[1][0] - ri[i], 2.0) + pow(output[2][0] - ri[i], 2.0);
		output_ans[i][0] = sqrt(temp) - output[3][0] * SPEED_LIGHT;
	}

	return output;
}

double RootFinding_NEWTONRAPHON(IN double initNR[MAX_ORDER][TRUE]) {

	matrix_substract(initNR, matrix_mul(matrix_inverse_4(matrix_Jacobian(initNR)), ans_function(initNR)));

	for (; seq < MAX_ITER; seq++) {

		if (seq > 1) {
			matrix_substract(output, matrix_mul(matrix_inverse_4(matrix_Jacobian(output)), ans_function(output)));
		}


		if (fabs(curNR[0][0] - preNR[0][0]) < TOL_RESIDUAL) {
			State = CONVERGE; Error_Message(State);	break;
		}

		for (int i = 0; i < MAX_ORDER; i++) { preNR[i][0] = curNR[i][0]; }
	}

	if (seq >= MAX_ITER) { State = ERR_ITER; Error_Message(State); }
}
#endif