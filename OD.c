#ifndef _ODE_
#define _ODE_

#include "OD.h"
#include "UserDefine.h"
#include "basicHeader.h"

//initial condition & defualt condition
double kp = 27.0;
double kd = 0.003;
double ts = 0.0005;
double tf = 2.0;
double zeta = 0.7;

double omega_a = 2 * M_PI * 20.5;
double theta_c = M_PI / 4;
double init_theta = M_PI / 4 + 1.0 * (M_PI / 180);
double omega_s = 2 * M_PI * 100.0;

double Lp = 3.2;
double Lpi = 1200;
double Ltheta = 16000;

//RK4
double dstate1[RKnum][sysnum] = { 0, };
double dstate2[RKnum][sysnum] = { 0, };
double state1_old[RKnum][sysnum] = { 0, };
double state1[RKnum][sysnum] = { 0, };
double state2[RKnum][sysnum] = { 0, };
double k[RKnum][sysnum] = { 0, };

void choose_init() {
	int wa, thetac = 0;
	printf("/*******************************chooose initial value********************************/\n");
	printf("/**actuator BW : 1 = 2*pi*20.5[rad/s]   2 = 2*pi*16.4[rad/s]  3 = 2*pi*24.6[rad/s]**/\n");
	printf("enter the number...");
	scanf("%d", &wa); getchar();
	printf("/*roll attitude cmd & initial roll angle : 1 = pi/4[rad]  2 = pi/8[rad]  3 = 0[rad]*/\n");
	printf("enter the number...");
	scanf("%d", &thetac); getchar();
	printf("/***********************************************************************************/\n");
	printf("start simulating....\n");

	switch (wa) {
	case 1:
		omega_a = 2 * M_PI * 20.5;
		break;
	case 2:
		omega_a = 2 * M_PI * 16.4;
		break;
	case 3:
		omega_a = 2 * M_PI * 24.6;
		break;
	default:
		printf("ERROR...\n");
		break;
	}

	switch (thetac) {
	case 1:
		theta_c = M_PI / 4;
		break;
	case 2:
		theta_c = M_PI / 8;
		break;
	case 3:
		theta_c = 0;
		break;
	default:
		printf("ERROR...\n");
		break;
	}

	init_theta = theta_c + 1.0 * (M_PI / 180);

	for (int i = 0; i < sysnum; i++) state1[2][i] = init_theta;

}

void sys() {
	double t = 0;
	copy_state(state1, state1_old);

	do {
		
		for (int i = 0; i < RKnum; i++) {
			dstate2[i][2] = sys3(state1[i][1], state1_old[i][2], dstate1[i][2]);
			dstate2[i][3] = sys4(state1[i][2], dstate1[i][3], state1[i][3]);
			dstate2[i][4] = sys4(dstate1[i][3], dstate1[i][4], state1[i][4]);
			dstate2[i][1] = sys2(state1[i][0], dstate1[i][1], state1[i][1]);
			dstate2[i][0] = sys1(state1[i][0], state1[i][4], state1[i][3]);
		}

		copy_state(state1, state1_old);

		RK4(ts, dstate2, dstate1);
		RK4(ts, dstate1, state1);
		

		t = t + ts;
		

	} while (t < tf);
	
	printf("\nThe output is : %lf\n\n", state1[0][2]);

}

double sys1(double input_x, double input_y, double input_y2) {
	double output = ((theta_c - input_y) * kp - input_y2)*kd;
	return output;
}


double sys2(double input_x, double input_y, double input_y2) {
	double output = pow(omega_a, 2) * input_x - 2 * zeta * input_y - pow(omega_a, 2) * input_y2;
	return output;
}

double sys3(double input_x, double input_y, double input_y2) {
	double output = Ltheta * input_x + Lpi * sin(4 * input_y) - Lp * input_y2;
	return output;
}

double sys4(double input_x, double input_y, double input_y2) {
	double output = pow(omega_s, 2) * input_x - 2 * zeta * input_y - pow(omega_s, 2) * input_y2;
	return output;
}

void RK4(double ts, double input[RKnum][sysnum], double output[RKnum][sysnum]) {

	for (int i = 0; i < RKnum; i ++ ) {
		for (int j = 0; j < sysnum; j++) {
			k[i][j] = ts * input[i][j];
		}
	}

	for (int h = 0; h < sysnum; h++) {
		output[0][h] = output[0][h] + (k[0][h] + 2 * k[1][h] + 2 * k[2][h] + k[3][h]) / 6;
		output[1][h] = output[0][h] + 0.5 * k[0][h];
		output[2][h] = output[0][h] + 0.5 * k[1][h];
		output[3][h] = output[0][h] + k[2][h];
	}

	
}

double copy_state(double input[RKnum][sysnum], double output[RKnum][sysnum]) {
	for (int i = 0; i < RKnum; i++) {
		for (int j = 0; j < sysnum; j++) {
			output[i][j] = input[i][j];
		}
	}
}
#endif