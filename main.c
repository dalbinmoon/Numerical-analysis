/*
2020_04_14 21700242 Sun-Been Moon
Numerical Analysis - HW3
*/

#include "basicHeader.h"
#include "UserDefine.h" 
#include "Sol_NonLinearEqns.h"
#include "Sol_TranserFunc.h"
#include "Linear_Algebra.h"
#include "Interpolation.h"

void main() {

	/*******************************
	HW3 - problem1
	*******************************/
	
	HW3_1();
	
	/*****************************
	HW3 - problem 2
	*****************************/

	HW3_2();
	
	/*****************************
	HW3 - problem 3
	*****************************/

	RootFinding(RECURSIVE_LMS);

	system("pause");

}