#ifndef _INTERPOLATION_
#define _INTERPOLATION_

#include "basicHeader.h"
#include "UserDefine.h"
#include "Interpolation.h"

double Lagrange_Coef(double coef[MAX_INTER], int in) {

}

double Langrange(double intput[MAX_INTER], double input2[MAX_INTER]) {

}

void InterPolation(int INTERPOL) {
	switch (INTERPOL) {

	case LAGRANGE:
		//RootFinding_BISECT(lowerB, upperB);
		break;

	case NEWTON_POL:
		//RootFinding_NEWTON(initN);
		break;

	default:
		break;
	}
}

double PI(double input[MAX_INTER], int size, double in) {
	double output = 1;

	for (int i = 0; i < MAX_INTER ; i++) {
		double temp = in - input[i];
		if (i == size) temp = 1;
		output = output * temp;
	}
	//printf("\n");

	return output;
}

double G(double input[MAX_INTER], int size, double in) {
	double output = 0, temp;
	//printf("size : %d\n", size);

	for (int i = 0; i < MAX_INTER; i++) {
		temp = in - input[i];
		temp = (double)1 / temp;
		if (i == size) temp = 0;
		output = output + temp;
		//printf("output = %.2f input = %.2f in = %.2f temp = %.2f\n", output, input[i], in, temp);
		
	}

	return output;
}

double Hermite(double input[MAX_INTER], double prime[MAX_INTER], double in) {
	double output = 0.0;

	//for (int i = 0; i < MAX_INTER; i++) {
		int i = 0;
		double a = 2 * G(input, i, time[i]) * sqrt(PI(input, i, in));
		double b = 2 * G(input, i, in) * sqrt(PI(input, i, time[i]));

		double c = sqrt(PI(input, i, in)) * (in - time[i]);
		double d = sqrt(PI(input, i, time[i]));

		output = output + a / b * input[i] + c / d * prime[i];
		//printf("N = %d output : %lf\n", i, output);
	//}

	return output;
}

double Hermite_Prime(double input[MAX_INTER], double prime[MAX_INTER], double in) {
	double output = 0;

	//for (int i = 0; i < MAX_INTER; i++) {
     	int i = 0;
		double a = sqrt(PI(input, i, in)) / sqrt(PI(input, i, time[i]));

		double b = 2 * (G(input, i, in) - G(input, i, time[i]));
		double c = 4 * G(input, i, in) * G(input, i, time[i]) * (time[i] - in);

		double d = 2 * G(input, i, in) * (in - time[i]) + 1;
		output = output + a*(b+c)*input[i] + a*d*prime[i];
		//printf("N = %d output : %lf\n", i, output);
	//}

	return output;
}


void HW2_1() {
	double hermitx = 0, hermity = 0, hermitz = 0, hermitvx = 0, hermitvy = 0, hermitvz = 0;
	printf("Problem 1\nHermite Interpolation...\n");
	printf("Input   x_position    y_position    z_position   x_velocity   y_velocity  z_velocity\n");

	for (int i = 0; i < MAX_INTER; i++) {
		hermitx = Hermite(init_x, init_vx, time[i]);
		hermity = Hermite(init_y, init_vy, time[i]);
		hermitz = Hermite(init_z, init_vz, time[i]);

		hermitvx = Hermite_Prime(init_x, init_vx, time[i]);
		hermitvy = Hermite_Prime(init_y, init_vy, time[i]);
		hermitvz = Hermite_Prime(init_z, init_vz, time[i]);
		
		printf("%5d  %lf  %lf   %lf   %lf    %lf    %lf\n", (int)time[i], hermitx, hermity, hermitz, hermitvx, hermitvy, hermitvz);

	}
	
}
#endif