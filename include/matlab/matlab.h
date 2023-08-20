//
// Created by Catigeart on 2023/6/21.
//

#ifndef CUBEGPA6SAB_MATLAB_H
#define CUBEGPA6SAB_MATLAB_H

void fix(double* X, int* Y, int len);

void deye(int n, double* a);

void dkron(int x_m, int x_n, const double* x, int y_m, int y_n, const double* y, double* z);

void cnorm2(int n, const double complex* x, int incx, double complex* norm_);

#endif //CUBEGPA6SAB_MATLAB_H
