#ifndef __RANDOM_H_
#define __RANDOM_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#define PI 3.1415926535897932384626

long setseed(void);
double ran2(double min, double max);
int ran2(int min, int max);
int ran2(int max);
double ran2(void);
double gammln(double xx);
double normrnd(void);
double gamrnd(double ia, double ib);
int poissrnd(double xm);
double exprnd(double mu);
int bnlrnd(double pp, int n);

#endif
