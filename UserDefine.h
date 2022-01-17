#ifndef _USERDEFINE
#define _USERDEFINE

#define ZERO_CROSSING(fx1, fx2)   (fx1*fx2 < 0) ? 1 : 0
#define SIGN(x)                   (x > 0) ? 1 : -1

#define MAX_ITER                  (unsigned int) ( 100  )
#define MAX_ORDER                 (unsigned int) (  10  )

#define TOL_RESIDUAL              (double)       ( 0.00001 )

typedef enum _METHOD_EQN          { BISECT = 1, NEWTON, SECANT }METHOD_EQN;
typedef enum _MESSAGE             { CONVERGE = 1, ERR_ITER, ERR_DIVERGE, OVER_FLOW }METHOD_STATE;
typedef enum _FREQUENCY {WGC = 1, WPC }FREQUENCY;

#endif