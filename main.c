/*
2020_05_08 21700242 Sun-Been Moon
Numerical Analysis - HW4
*/

#include "basicHeader.h"
#include "UserDefine.h" 
#include "Sol_NonLinearEqns.h"
#include "Sol_TranserFunc.h"
#include "Linear_Algebra.h"
#include "Interpolation.h"
#include "Optimization.h"

void main() {

	/*******************************
	HW4 - problem2
	*******************************/

	printf("/******* HW4 *******/");
	Gradient_Descent();
	Steepest_Gradient_Descent();

	system("pause");

}