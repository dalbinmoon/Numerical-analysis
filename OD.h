#ifndef _ORDINARY_
#define _ORDINARY_

#define sysnum (int) 5
#define RKnum  (int) 4

//initial condition & defualt condition
extern double kp;
extern double kd;
extern double ts;
extern double tf;
extern double zeta;

extern double omega_a;
extern double init_theta;
extern double theta_c;

extern double omega_s;
extern double Lp;
extern double Lpi;
extern double Ltheta;

//RK4
extern double dstate1[RKnum][sysnum];
extern double dstate2[RKnum][sysnum];
extern double state1_old[RKnum][sysnum];
extern double state1[RKnum][sysnum];
extern double state2[RKnum][sysnum];
extern double k[RKnum][sysnum];



void choose_init();

double sys1(double input_x, double input_y, double input_y2);

double sys2(double input_x, double input_y, double input_y2);

double sys3(double input_x, double input_y, double input_y2);

double sys4(double input_x, double input_y, double input_y2);

double copy_state(double input[RKnum][sysnum], double output[RKnum][sysnum]);

void sys();

void RK4(double ts, double input[RKnum][sysnum], double output[RKnum][sysnum]);


#endif
#pragma once
