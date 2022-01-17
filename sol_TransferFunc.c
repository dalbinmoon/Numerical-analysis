#ifndef _SOL_TRANSFERFUNCTION
#define _SOL_NONLINEAREQUATION

#include "basicHeader.h"
#include "UserDefine.h"
#include "C:\Users\MoonSunBeen\Desktop\HGU 20-1\수치해석\Exp\myFunc\sol_TranserFunc.h"

extern void Error_Message(IN unsigned int _Errortype);
extern void Disp(IN int _Method, IN int seq, IN float disp1, IN float disp2, IN float disp3, IN float disp4, IN float disp5);

void print_GoOpt(void) {

	achieve_Go();

	printf("\nThe inserted open loop system model is as follows:\n");
	print_Go(denuminator_size, numerator_size);

	char check[TRUE];
	printf("\nIs it correct?[y/n] : ");
	scanf("%*[\n]%c", check);
	if (check[0] == 'y') {
		printf("\nCalculating Gain Crossover and Phase Crossover Frequencies...\n");
		
	}
	if (check[0] == 'n') {
		printf("Achieve parameter Again\n");
		print_GoOpt();
	}

}

void achieve_Go(void) {

	printf("\nInput the order of denminator             : ");
	scanf("%d", &denuminator_size);
	if (denuminator_size > 10) {
		State = OVER_FLOW; Error_Message(State); exit(0);
	}
	denuminator_size++;

	char check[TRUE];
	printf("What are the coefficients of denominator? : ");
	scanf("%*[\n]%c", check);
	for (int i = 0; i < denuminator_size; i++) {
		scanf("%lf", &denuminator[i]);
	}
	scanf("%c", check);

	printf("Input the order of numerator              : ");
	scanf("%d", &numerator_size);
	if (numerator_size > 10) {
		State = OVER_FLOW; Error_Message(State); exit(0);
	}
	numerator_size++;

	printf("What are the coefficients of numerator?   : ");
	scanf("%*[\n]%c", check);
	for (int i = 0; i < numerator_size; i++) {
		scanf("%lf", &numerator[i]);
	}
	scanf("%c", check);
}

void print_Go(IN int num1, IN int num2) {
	printf("          ");
	int num_num = num2 - 1;
	for (int i = 0; i < num2; i++) {
		printf("%.2f", numerator[i]);
		if (num_num > 1) printf("s^%d ", num_num);
		if (num_num == 1) printf("s ", num_num);
		if (num_num > 0) printf(" + ");
		num_num--;
	}

	printf("\nGo(s) = --------------------------\n");
	printf("          ");

	int num_den = num1 - 1;
	for (int i = 0; i < num1; i++) {
		printf("%.2f", denuminator[i]);
		if (num_den > 1) printf("s^%d ", num_den);
		if (num_den == 1) printf("s ", num_den);
		if (num_den > 0) printf(" + ");
		num_den--;
	}
	printf("\n");
}

void calculate_feq(void) {

	for (int i = 0; i < denuminator_size; i++) {
		int order = denuminator_size - i - 1;
		if (i == (denuminator_size - 1)) {
			den_w_ideal[i][0] = denuminator[i];
			den_w_ideal[i][1] = order;
		}
		else if (order % 2 == 0) {
			if (order % 4 == 3 || order % 4 == 2) den_w_ideal[i][0] = (-denuminator[i]);
			else den_w_ideal[i][0] = (denuminator[i]);
			den_w_ideal[i][1] = (double)order;
		}
		else if (order % 2 == 1) {
			if (order % 4 == 3 || order % 4 == 2) den_w_img[i][0] = (-denuminator[i]);
			else den_w_img[i][0] = (denuminator[i]);
			den_w_img[i][1] = (double)order;
		}
	}

	for (int i = 0; i < numerator_size; i++) {
		int order = numerator_size - i - 1;
		if (i == (numerator_size - 1)) {
			num_w_ideal[i][0] = numerator[i];
			num_w_ideal[i][1] = order;
		}
		else if (order % 2 == 0) {
			if (order % 4 == 3 || order % 4 == 2) num_w_ideal[i][0] = (-numerator[i]);
			else num_w_ideal[i][0] = (numerator[i]);
			num_w_ideal[i][1] = (double)order;
		}
		else if (order % 2 == 1) {
			if (order % 4 == 3 || order % 4 == 2) num_w_img[i][0] = (-numerator[i]);
			else num_w_img[i][0] = (numerator[i]);
			num_w_img[i][1] = (double)order;
		}

		printf("\n num  : %lf", num_w_img[i][0]);
	}

	//calculate_wpc(); 
	calculate_wgc(); 


	freq = WGC;
	WFinding_NEWTON(4, freq);
	//freq = WPC;
	//WFinding_NEWTON(12, freq);

}

void calculate_wpc(void) {
	for (int i = 0; i < numerator_size; i++) {

		for (int j = 0; j < denuminator_size; j++) {
			int order = num_w_ideal[i][1] + den_w_img[j][1];
			wpc_img_order[order] = order;
			wpc_img_val[order] = wpc_img_val[order] + num_w_ideal[i][0] * den_w_img[j][0];
		}
		for (int j = 0; j < denuminator_size; j++) {
			int order = num_w_img[i][1] + den_w_ideal[j][1];
			wpc_img_order[order] = order;
			wpc_img_val[order] = wpc_img_val[order] + num_w_img[i][0] * den_w_ideal[j][0];
		}
	}
}

void calculate_wgc(void) {

	for (int i = 0; i < numerator_size; i++) {
		
		for (int j = 0; j < numerator_size; j++) {
			int ordera = num_w_ideal[i][1] + num_w_ideal[j][1];
			num_wgcordera[ordera] = ordera;
			num_wgcvala[ordera] = num_wgcvala[ordera] + num_w_ideal[i][0] * num_w_ideal[j][0];

			int orderb = num_w_img[i][1] + num_w_img[j][1];
			num_wgcorderb[orderb] = orderb;
			num_wgcvalb[orderb] = num_wgcvala[orderb] + num_w_img[i][0] * num_w_img[j][0];
			printf("orderb : %d\n", orderb);
			printf("orvalb : %lf\n", num_wgcvala[orderb]);
		}

		

	}

	for (int i = 0; i < denuminator_size; i++) {
		
		for (int j = 0; j < denuminator_size; j++) {
			int orderc = den_w_ideal[i][1] + den_w_ideal[j][1];
			den_wgcorderc[orderc] = orderc;
			den_wgcvalc[orderc] = den_wgcvalc[orderc] + den_w_ideal[i][0] * den_w_ideal[j][0];
		}

		for (int j = 0; j < denuminator_size; j++) {
			int orderd = den_w_img[i][1] + den_w_img[j][1];
			den_wgcorderd[orderd] = orderd;
			den_wgcvald[orderd] = den_wgcvald[orderd] + den_w_img[i][0] * den_w_img[j][0];
		}

	}

	printf("1 : %lf\n", num_wgcvala[1]);
	printf("2 : %lf\n", num_wgcvala[2]);


	for (int i = 0; i < MAX_ORDER; i++) {
		wgc_val_img[i] = num_wgcvala[i] - den_wgcvalc[i];
		wgc_order_img[i] = den_wgcorderc[i];
		printf(" j : %lf %d\n", wgc_val_img[i], wgc_order_img[i]);
		wgc_val_ideal[i] = num_wgcvalb[i] - den_wgcvald[i];
		wgc_order_ideal[i] = den_wgcorderd[i];
	}

}

float sDiff0_Func(IN float x, IN int w) {
	float y = 0;
	if (w == wgc) {
		for (int i = 0; i < sizeof(wgc_order_img) / sizeof(int); i++) {
			y = y + wgc_val_img[i] * pow(x, wgc_order_img[i]);
		}
		return y;
	}
	if (w == wpc) {
		for (int i = 0; i < sizeof(wgc_order_img) / sizeof(int); i++) {
			y = y + wgc_val_img[i] * pow(x, wgc_val_img[i]);
		}
		return y;
	}
}

float sDiff1_Func(IN float x, IN int w) {
	float y = 0;
	if (w == wgc) {
		for (int i = 0; i < sizeof(wgc_order_img) / sizeof(int)-1; i++) {
			y = y + wgc_order_img[i] * wgc_order_img[i] * pow(x, wgc_order_img[i]-1);
		}
		return y;
	}
	if (w == wpc) {
		for (int i = 1; i < sizeof(wpc_img_order) / sizeof(int); i++) {
			y = y + wpc_img_val[i] * wpc_img_order[i] * pow(x, wpc_img_order[i]-1);
		}
		return y;
	}
}

float WFinding_NEWTON(IN float initN) {

	prewN = initwN - sDiff0_Func(initN, freq) / sDiff1_Func(initN, freq);

	for (; sequence < MAX_ITER; sequence++) {

		if (sequence > 1) curwN = prewN - sDiff0_Func(prewN, freq) / sDiff1_Func(prewN, freq);

		Disp(NEWTON, sequence, curwN, sDiff1_Func(curwN, freq), sDiff0_Func(curwN, freq), 0.0, 0.0);

		if (fabs(curwN - prewN) < TOL_RESIDUAL) {
			State = CONVERGE; Error_Message(State); break;
		}

		prewN = curwN;
	}

	if (sequence >= MAX_ITER) { State = ERR_ITER; Error_Message(State); }
}
#endif