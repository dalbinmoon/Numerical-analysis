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

float Det_3(IN float Det[READ_COLUMN][READ_COLUMN]){
	float det = Det[0][0] * Det[1][1] * Det[2][2] - Det[0][0] * Det[1][2] * Det[2][1] + Det[0][1] * Det[1][2] * Det[2][0] - Det[0][1] * Det[1][0] * Det[2][2] + Det[0][2] * Det[1][0] * Det[2][1] - Det[0][2] * Det[1][1] * Det[2][0];

	if (det == 0.0) { State = NO_DET; Error_Message(State); return 0.0; }
	return det;
}

double Det_4(IN double Det[MAX_ORDER][MAX_ORDER]){
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

void matrix_inverse_3(IN float inv[READ_COLUMN][READ_COLUMN], OUT float out[READ_COLUMN][READ_COLUMN]){
	float det = 1/Det_3(inv);

	out[0][0] = det * (inv[1][1] * inv[2][2] - inv[1][2] * inv[2][1]);
	out[0][1] = det * (inv[0][2] * inv[2][1] - inv[0][1] * inv[2][2]);
	out[0][2] = det * (inv[0][1] * inv[1][2] - inv[0][2] * inv[1][1]);

	out[1][0] = det * (inv[1][2] * inv[2][0] - inv[1][0] * inv[2][2]);
	out[1][1] = det * (inv[0][0] * inv[2][2] - inv[0][2] * inv[2][0]);
	out[1][2] = det * (inv[0][2] * inv[1][0] - inv[0][0] * inv[1][2]);

	out[2][0] = det * (inv[1][0] * inv[2][1] - inv[1][1] * inv[2][0]);
	out[2][1] = det * (inv[0][1] * inv[2][0] - inv[0][0] * inv[2][1]);
	out[2][2] = det * (inv[0][0] * inv[1][1] - inv[0][1] * inv[1][0]);

	return 0; 
}

void matrix_mul(IN float mat1[READ_COLUMN][READ_COLUMN], IN float mat2[READ_COLUMN][READ_COLUMN], OUT float out[READ_COLUMN][READ_COLUMN]) {

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				out[i][j] = out[i][j] + mat1[i][k] * mat2[k][j];
			}
		}
	}

	return 0;
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

void matrix_add(IN float mat1[READ_COLUMN][READ_COLUMN], IN float mat2[READ_COLUMN][READ_COLUMN], OUT float out[READ_COLUMN][READ_COLUMN]) {

	for (int i = 0; i < READ_COLUMN; i++) {
		for (int j = 0; j < READ_COLUMN; j++) {
			out[i][j] = mat1[i][j] + mat2[i][j];
			
		}
	}
	
	return 0;
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

	//matrix_substract(initNR, matrix_mul(matrix_inverse_4(matrix_Jacobian(initNR)), ans_function(initNR)));

	for (; seq < MAX_ITER; seq++) {

		if (seq > 1) {
			//matrix_substract(output, matrix_mul(matrix_inverse_4(matrix_Jacobian(output)), ans_function(output)));
		}


		if (fabs(curNR[0][0] - preNR[0][0]) < TOL_RESIDUAL) {
			State = CONVERGE; Error_Message(State);	break;
		}

		for (int i = 0; i < MAX_ORDER; i++) { preNR[i][0] = curNR[i][0]; }
	}

	if (seq >= MAX_ITER) { State = ERR_ITER; Error_Message(State); }
}

void transpose_matrix(IN float input[READ_COLUMN][READ_COLUMN], OUT float result[READ_COLUMN][READ_COLUMN]) {

	for (int i = 0; i < READ_COLUMN; i++) {
		for (int j = 0; j < READ_COLUMN; j++) {
			result[i][j] = input[j][i];
		}
	}

	return 0;
}

void HW3_1() {
	printf("\nProblem 1\nFinding coefficient a1, ao...\n");
	printf("The equation is Hx = y\n");
	printf("H = ");
	
	for (int i = 0; i < 11; i++) {
		if(i > 0) printf("     ");
		for (int j = 0; j < 2; j++) {
			printf("%4.f ", H[i][j]);
		}
		printf("\n");
	}

	printf("Y = ");
	for (int i = 0; i < 11; i++) printf("%3.2f ", P[i]);
	printf("\n");

	double result[2][11] = { 0 };
	transpose_matrix_3_1(H, result);

	double result1[11][11] = { 0 };
	matrix_mul_3_1(result, H, result1);

	double result2[2][2] = { 0 };
	matrix_inverse_2(result1, result2);

	double result3[2][1] = { 0 };
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 1; j++) {
			for (int k = 0; k < 11; k++) {
				result3[i][j] = result3[i][j] + result[i][k] * P[k];
			}
		}
	}

	double x[2] = { 0 };
	for (int i = 0; i < 2; i++) {
		for (int k = 0; k < 2; k++) {
			x[i] = x[i] + result2[i][k] * result3[k][0];
		}
	}

	printf("\na1 : %lf a0 : %lf \n", x[0], x[1]);
	
}

void transpose_matrix_3_1(IN double input[11][2], double result[2][11]) {

	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 11; j++) {
			result[i][j] = input[j][i];
		}
	}

}

void matrix_mul_3_1(IN double input[2][11], IN double input2[11][2], double output[2][2]) {

	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			for (int k = 0; k < 11; k++) {
				output[i][j] = output[i][j] + input[i][k] * input2[k][j];
			}
		}
	}

	return 0; 
}

void matrix_inverse_2(IN double inv[2][2], OUT double output[2][2]) {
	double det = 1 / (inv[0][0] * inv[1][1] - inv[0][1] * inv[1][0]);

	output[0][0] = det * inv[1][1];
	output[0][1] = -det * inv[0][1];
	output[1][0] = -det * inv[1][0];
	output[1][1] = det * inv[0][0];

	return 0;
}

void HW3_2() {
	//initialize H2 & V
	int t = 2;
	for (int i = 0; i < 15; i++) {
		H2[i][0] = -t/R; H2[i][1] = 1;
		t = t + 2;
		V[i] = log(V[i]);
	}

	printf("\n\nProblem 2\nFinding coefficient c, v...\n");
	printf("The equation is Hx = y\n");
	printf("H = ");

	for (int i = 0; i < 15; i++) {
		if (i > 0) printf("    ");
		for (int j = 0; j < 2; j++) {
			printf("%lf ", H2[i][j]);
		}
		printf("\n");
	}

	printf("Y = ");
	for (int i = 0; i < 11; i++) printf("%3.2f ", V[i]);
	printf("\n");

	double result[2][15] = { 0 };
	transpose_matrix_3_2(H2, result);

	double result1[15][15] = { 0 };
	matrix_mul_3_2(result, H2, result1);

	double result2[2][2] = { 0 };
	matrix_inverse_2(result1, result2);

	double result3[2][1] = { 0 };
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 1; j++) {
			for (int k = 0; k < 15; k++) {
				result3[i][j] = result3[i][j] + result[i][k] * V[k];
			}
		}
	}

	double x[2] = { 0 };
	for (int i = 0; i < 2; i++) {
		for (int k = 0; k < 2; k++) {
			x[i] = x[i] + result2[i][k] * result3[k][0];
		}
	}

	printf("\n1/C : %lf V : %lf \n", (x[0]), exp(x[1]));


}

void transpose_matrix_3_2(IN double input[15][2], double result[2][15]) {

	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 15; j++) {
			result[i][j] = input[j][i];
		}
	}

}

void matrix_mul_3_2(IN double input[2][15], IN double input2[15][2], double output[2][2]) {

	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			for (int k = 0; k < 15; k++) {
				output[i][j] = output[i][j] + input[i][k] * input2[k][j];
			}
		}
	}

	return 0;
}
#endif