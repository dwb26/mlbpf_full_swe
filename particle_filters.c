#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <assert.h>
#include "solver.h"
#include "particle_filters.h"

double sigmoid(double x, double a, double b) {
	return a / (1.0 + exp(-0.01 * M_PI * x)) + b;
}


double sigmoid_inv(double x, double a, double b) {
	return log((x - b) / (a + b - x)) / (0.01 * M_PI);
}