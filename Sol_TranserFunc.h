#ifndef _SOL_TRANSFERFUNC
#define _SOL_TRANSFERFUNC

METHOD_EQN     Method;
METHOD_STATE   State;
FREQUENCY      freq;

//wNewton Method 
static int sequence = 0;
static float prewN = 0.0;
static float curwN = 0.0;

//Go parameter
static int denuminator_size;
static int numerator_size;
static double denuminator[MAX_ORDER] = { 0, };
static double numerator[MAX_ORDER] = { 0, };

static double output_wgc = 0.0;
static double output_wpc = 0.0;


void   print_GoOpt        (void);

void   achieve_Go         (void);

void   print_Go           (IN int num1, IN int num2);

void   calculate_feq      (void);

void   switch_freq        (IN int _freq);

float  sDiff0_Func        (IN float x, IN int _freq);

float  sDiff1_Func        (IN float x, IN int _freq);

float  WFinding_NEWTON    (IN float initwN, IN int freq);

void   Error_Message      (IN unsigned int _Errortype);

#endif
