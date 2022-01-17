#ifndef _SOL_TRANSFERFUNCTION
#define _SOL_NONLINEAREQUATION

#include "basicHeader.h"
#include "UserDefine.h"
#include "Sol_TranserFunc.h"
#include "Sol_NonLinearEqns.h"

void print_GoOpt(void) {

	achieve_Go();

	printf("\nThe inserted open loop system model is as follows:\n");
	print_Go(denuminator_size+1, numerator_size+1);

	char check[TRUE];
	printf("\n\nIs it correct?[y/n] : ");
	scanf("%*[\n]%c", check);
	if (check[0] == 'y') {
		printf("\nCalculating Gain Crossover and Phase Crossover Frequencies...\n");
		
	}
	if (check[0] == 'n') {
		printf("\nAchieve parameter Again\n");
		for (int i = 0; i < MAX_ORDER + 1; i++) {
			denuminator[i] = 0.0;
			numerator[i] = 0.0;
		}
		print_GoOpt();
	}

}

void achieve_Go(void) {
	char check[TRUE];

	printf("\nInput the order of denminator             : ");
	scanf("%d", &denuminator_size);
	if (denuminator_size > MAX_ORDER) {
		State = OVER_FLOW; Error_Message(State); exit(0);
	}

	printf("What are the coefficients of denominator? : ");
	scanf("%*[\n]%c", check);
	for (int i = denuminator_size; i > -1; i--) {
		scanf("%lf", &denuminator[i]);
	}
	scanf("%c", check);

	printf("Input the order of numerator              : ");
	scanf("%d", &numerator_size);
	if (numerator_size > MAX_ORDER) {
		State = OVER_FLOW; Error_Message(State); exit(0);
	}

	printf("What are the coefficients of numerator?   : ");
	scanf("%*[\n]%c", check);
	for (int i = numerator_size; i > -1 ; i--) {
		scanf("%lf", &numerator[i]);
	}
	scanf("%c", check);
}

void print_Go(IN int num1, IN int num2) {
	printf("         ");
	for (int i = MAX_ORDER-1; i > -1; i--) {
		if (numerator[i] != 0) {
			printf("%5.2f", numerator[i]);

			if (i == 1)  printf("s ", i);
			else if (i > 1) printf("s^%d ", i);

			if (i > 0 && numerator[i - 1] != 0) printf(" + ");
		}
		else printf("            ");
	}

	printf("\nGo(s) = --------------------------------------------\n");
	printf("         ");

	for (int i = MAX_ORDER-1; i > -1; i--) {
		if (denuminator[i] != 0) {
			printf("%5.2f", denuminator[i]);

			if (i == 1)  printf("s ", i);
			else if (i > 1
				) printf("s^%d ", i);

			if (i > 0 && denuminator[i - 1] != 0) printf(" + ");
		}
		else printf("            ");
	}
}

void calculate_feq(void) {	
	printf("\n[NEWTON METHOD]");
	printf("\nFinding Gain over Frequency....\n");
	switch_freq(WGC);
	printf("\nFinding Phase over Frequency....\n");
	switch_freq(WPC);

	if (output_wgc == -1) printf("\n\nGain Crossover Frequency Wgc  =    -     [rad/s]", output_wgc);
	else printf("\n\nGain Crossover Frequency Wgc  = %6.2f[rad/s]", output_wgc);

	if(output_wpc == -1) printf("\nPhase Crossover Frequency Wpc =    -     [rad/s]\n", output_wpc);
	else printf("\nPhase Crossover Frequency Wpc = %6.2f[rad/s]\n", output_wpc);
}

void switch_freq(freq) {
	switch (freq) {

	case WGC:
		WFinding_NEWTON(9.5, WGC);
		output_wgc = curwN;
		break;

	case WPC:
		WFinding_NEWTON(12.8, WPC);
		output_wpc = curwN;
		break;

	default:
		break;
	}
}

float sDiff0_Func(IN float x, IN int _freq) {
	float y = 0.0;
	if (_freq == WGC) {
		y = (pow(numerator[3], 2) - 1) * pow(x, 6.0) + (pow(numerator[2], 2.0) - pow(denuminator[2], 2.0) - 2 * numerator[1] * numerator[3] + 2.0 * denuminator[1]) * pow(x, 4.0) + (pow(numerator[1], 2.0) - pow(denuminator[1], 2.0) - 2 * numerator[0] * numerator[2] + 2.0 * denuminator[0] * denuminator[2]) * pow(x, 2.0) + pow(numerator[0], 2.0) - pow(denuminator[0], 2.0);
	}

	else if (_freq == WPC) {
		float a = -numerator[2] * pow(x, 2.0) + numerator[0];
		float b = -numerator[3] * pow(x, 3.0) + numerator[1] * x;
		float c = -denuminator[2] * pow(x, 2.0) + denuminator[0];
		float d = -denuminator[3] * pow(x, 3.0) + denuminator[1] * x;
		float input = (b * c - a * d) / (a * c + b * d);

		y = atan(input) - M_PI;
	}

	return y;

}

float sDiff1_Func(IN float x, IN int _freq) {
	float y = 0.0;
	if (_freq == WGC) {
		y = 6.0 * (pow(numerator[3], 2) - 1) * pow(x, 5.0) + 4.0 * (pow(numerator[2], 2.0) - pow(denuminator[2], 2.0) - 2 * numerator[1] * numerator[3] + 2.0 * denuminator[1]) * pow(x, 3.0) + 2.0 * (pow(numerator[1], 2.0) - pow(denuminator[1], 2.0) - 2 * numerator[0] * numerator[2] + 2.0 * denuminator[0] * denuminator[2]) * pow(x, 1.0) ;
		
	}
	else if (_freq == WPC) {
		float left = 5 * (denuminator[2] * numerator[3] - denuminator[3] * numerator[2]) * pow(x, 4.0) + 3 * (-denuminator[0] * numerator[3] - denuminator[2] * numerator[1] + numerator[1] * numerator[2] + numerator[0] * numerator[3]) * pow(x, 2.0) + denuminator[0] * numerator[1] - numerator[0] * numerator[1];
		float right = tan(M_PI)* (6.0*denuminator[3] * numerator[3] * pow(x, 5.0) - 4.0*(denuminator[3] * numerator[1] + numerator[1] * numerator[3] - denuminator[2] * numerator[2]) * pow(x, 3.0) - 2*(denuminator[0] * numerator[2] + denuminator[2] * numerator[0] - numerator[1] * numerator[1]) * pow(x, 1.0));

		y = left - right;
	}
	return y;
}

float WFinding_NEWTON(IN float initwN, IN int freq) {

	prewN = initwN - sDiff0_Func(initwN, freq) / sDiff1_Func(initwN, freq);


	for (; sequence < MAX_ITER; sequence++) {

		
		curwN = prewN - sDiff0_Func(prewN, freq) / sDiff1_Func(prewN, freq);

		Disp(WFINDINGNEWTON, sequence, curwN, 0.0, 0.0, 0.0, 0.0);

		if (fabs(curwN - prewN) < TOL_RESIDUAL) {

			State = CONVERGE; Error_Message(State); 
			
			break;
		}

		prewN = curwN;
	}

	if (sequence >= MAX_ITER) {
		State = ERR_ITER; Error_Message(State);
		if (freq == WGC) { output_wgc = FALSE; }
		else if (freq == WPC) { output_wpc = FALSE; }
	}
}
#endif