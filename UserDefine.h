#ifndef _USERDEFINE
#define _USERDEFINE

#define ZERO_CROSSING(fx1, fx2)   (fx1*fx2 < 0) ? 1 : 0
#define SIGN(x)                   (x > 0) ? 1 : -1
#define ABSOLUTE(x)               (x >= 0 ) ? 1 : 0      

#define MAX_ITER                  (unsigned int) ( 100  )
#define MAX_ORDER                 (unsigned int) (   4  )
#define MAX_INTER                 (unsigned int) (   5  )

#define SPEED_LIGHT               (int)          (299792458) //m/s
#define TOL_RESIDUAL              (double)       ( 0.00001 )
#define TOL_ERROR                 (double)       ( 0.00001 )

#define READ_ROW                  (int)          (100)
#define READ_COLUMN               (int)          (3)

typedef enum _METHOD_EQN          { BISECT = 1, NEWTON, SECANT, WFINDINGNEWTON, NEWTONRAPHON, RECURSIVE_LMS}METHOD_EQN;
typedef enum _MESSAGE             { CONVERGE = 1, ERR_ITER, ERR_DIVERGE, OVER_FLOW, NO_DET }METHOD_STATE;
typedef enum _FREQUENCY           { WGC = 1, WPC }FREQUENCY;
typedef enum _INTERPOLATION       { LAGRANGE = 1, NEWTON_POL }INTERPOLATION;

#endif