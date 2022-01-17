#ifndef _SOL_TRANSFERFUNC
#define _SOL_TRANSFERFUNC

METHOD_EQN     Method;
METHOD_STATE   State;
FREQUENCY freq;

//wNewton Method 
static int sequence = 0;
static float initwN = 1.0;
static float prewN = 0.0;
static float curwN = 0.0;

#define wgc (int)( 1)
#define wpc (int)( 2)

//Go parameter
static int denuminator_size;
static int numerator_size;
static double denuminator[MAX_ORDER] = { 0, };
static double numerator[MAX_ORDER] = { 0, };

//w parameter
static double den_w_ideal[MAX_ORDER][2];
static double den_w_img[MAX_ORDER][2];

static double num_w_ideal[MAX_ORDER][2];
static double num_w_img[MAX_ORDER][2];

//w parameter
static double wpc_img_val[MAX_ORDER];
static int wpc_img_order[MAX_ORDER];

static int num_wgcordera[MAX_ORDER];
static double num_wgcvala[MAX_ORDER];
static int num_wgcorderb[MAX_ORDER];
static double num_wgcvalb[MAX_ORDER];
static int den_wgcorderc[MAX_ORDER];
static double den_wgcvalc[MAX_ORDER];
static int den_wgcorderd[MAX_ORDER];
static double den_wgcvald[MAX_ORDER];

static double wgc_val_img[MAX_ORDER];
static int wgc_order_img[MAX_ORDER];

static double wgc_val_ideal[MAX_ORDER];
static int wgc_order_ideal[MAX_ORDER];

//Newton's Method
extern float initN;
extern float preN;
extern float curN;


void   print_MethodOpt(OUT unsigned int* _Mehtod);

void   print_GoOpt(void);

void   achieve_Go(void);

void   print_Go(IN int num1, IN int num2);

void   calculate_feq(void);

void   calculate_wpc(void);

void   calculate_wgc(void);

float sDiff0_Func(IN float x, IN int w);

float sDiff1_Func(IN float x, IN int w);

float WFinding_NEWTON(IN float initN);

void   Error_Message(IN unsigned int _Errortype);

#endif
