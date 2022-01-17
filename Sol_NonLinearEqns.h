#ifndef _SOL_NONLINEAREQNS
#define _SOL_NONLINEAREQNS

METHOD_EQN     Method;
METHOD_STATE   State;
static int     printflag = 0;

//Common
static int seq = 1;
static float residual;
static float bound;

//Bisection Method
static float lowerB = 0.0;
static float upperB = 1.0;
static float halfB;

//Newton's Method
static float initN = 1.0;
static float preN  = 0.0;
static float curN  = 0.0;

//Secant Method
static float initS0 = 0.0;
static float initS1 = 1.0;
static float nextS  = 0.0;
static float curS   = 0.0;
static float preS   = 0.0;

void   print_MethodOpt            (OUT unsigned int *_Mehtod);

float  Diff0_Func			      (IN float _x);

float  Diff1_Func                 (IN float _x);

void   RootFinding                (IN int _Method);

float  RootFinding_BISECT         (IN float _lowerB, IN float _upperB);

float  RootFinding_NEWTON         (IN float _initN);

float  RootFinding_SECANT         (IN float _initS0, IN float _initS1);

void   Disp                       (IN int _Method, IN int seq, IN float disp1, IN float disp2, IN float disp3, IN float disp4, IN float disp5);

void   Error_Message              (IN unsigned int _Errortype);

#endif
