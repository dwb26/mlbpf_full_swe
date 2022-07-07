#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include "solver.h"
#include <assert.h>

const double GR = 9.81;
const double ZERO_THRESH = 1e-16;
const double GLOB_EPS = 1e-10;
const double C_CUTOFF = 1.35;

void gen_Z_drain(int nx, double * xs, double * Z, double height, double centre) {
	for (int j = 0; j < nx + 2; j++) {
		Z[j] = height - 0.05 * pow(xs[j] - centre, 2);
		Z[j] = Z[j] > 0 ? Z[j] : 0.0;
	}
}


void apply_topography(double * Z_arr, double * xs, double k, double gamma_of_k, double theta, int nx) {
	double gamma_term = -gamma_of_k * pow(theta, k) / 25.0;
	for (int j = 0; j < nx + 2; j++)
		Z_arr[j] = pow(xs[j], k - 1) * exp(-xs[j] / theta) / gamma_term;
}


void prescribe_left_bcs(double ** W, double dx, double t) {
	W[0][0] = W[1][0], W[0][1] = 0.0;
}


void prescribe_right_bcs(double ** W, int nx) {

	double uN, hN = W[nx][0];
	if (hN <= ZERO_THRESH)
		uN = 0.0;
	else
		uN = W[nx][1] / hN;
	double hR = 1.0 / (9.0 * GR) * pow(uN + 2 * sqrt(GR * hN), 2);
	hR = hR < hN ? hR : hN;
	double qR = hR / 3.0 * (uN + 2 * sqrt(GR * hN));
	W[nx + 1][0] = hR, W[nx + 1][1] = qR;

}


void output_data(double ** W, double * Z, int nx, double t, double sig_theta, FILE * CURVE_DATA, FILE * TOP_DATA, FILE * TIMES) {

	for (int j = 1; j < nx + 1; j++)
		fprintf(CURVE_DATA, "%.16e ", W[j][0] + Z[j]);
	fprintf(CURVE_DATA, "\n");
	for (int j = 1; j < nx + 1; j++)
		fprintf(TOP_DATA, "%.16e ", Z[j]);
	fprintf(TOP_DATA, "\n");
	fprintf(TIMES, "%.5lf ", t);
}


void zero_cutoff(double * W_L, double * W_R) {

	double h_L = W_L[0], h_R = W_R[0];
	if ( (h_L < ZERO_THRESH) && (h_R < ZERO_THRESH) )
		W_L[0] = 0.0, W_L[1] = 0.0, W_R[0] = 0.0, W_R[1] = 0.0;
	else if (h_L < ZERO_THRESH)
		W_L[0] = 0.0, W_L[1] = 0.0;
	else if (h_R < ZERO_THRESH)
		W_R[0] = 0.0, W_R[1] = 0.0;
}


void compute_wave_speeds(double * W_L, double * W_R, double * lmbda_neg, double * lmbda_pos, int j) {
	/**
	 * Compute the negative travelling (< 0) and positive travelling (> 0) wave speeds
	*/
	double h_L = W_L[0], h_R = W_R[0], c_L, u_L, c_R, u_R, lmbda_L, lmbda_R, C;
	assert( (h_L >= 0.0) && (h_R >= 0.0) );

	/* Extract the velocities and the speed of sounds */
	if (fabs(h_L) < ZERO_THRESH)
		c_L = 0.0, u_L = 0.0;
	else
		c_L = sqrt(GR * h_L), u_L = W_L[1] / h_L;
	if (fabs(h_R) < ZERO_THRESH)
		c_R = 0.0, u_R = 0.0;
	else
		c_R = sqrt(GR * h_R), u_R = W_R[1] / h_R;

	/* Compute the left and right wave speeds */
	lmbda_L = -fabs(u_L) - c_L;
	C = -fabs(u_R) - c_R;
	lmbda_L = C > lmbda_L ? lmbda_L : C;
	lmbda_L = -GLOB_EPS > lmbda_L ? lmbda_L : -GLOB_EPS;

	lmbda_R = fabs(u_L) + c_L;
	C = fabs(u_R) + c_R;
	lmbda_R = lmbda_R > C ? lmbda_R : C;
	lmbda_R = lmbda_R > GLOB_EPS ? lmbda_R : GLOB_EPS;

	assert( (lmbda_L < 0) && (lmbda_R > 0) );
	lmbda_neg[j] = lmbda_L;
	lmbda_pos[j] = lmbda_R;

}


void flux(double * W, double * W_flux) {
	/**
	 * FLux function for the system. This is used in the computation of the HLL terms.	
	*/
	double h = W[0], q = W[1];
	if (h >= ZERO_THRESH) {
		W_flux[0] = q;
		W_flux[1] = q * q / h + 0.5 * GR * h * h;
	}
}


void compute_W_HLL(double * W_HLL, double lmbda_L, double lmbda_R, double * W_L, double * W_R, double * W_flux_L, double * W_flux_R) {
	flux(W_L, W_flux_L);
	flux(W_R, W_flux_R);
	// W_HLL[0] = (lmbda_R * W_R - lmbda_L * W_L + flux(W_L, W_flux) - flux(W_R, W_flux)) / (lmbda_R - lmbda_L);
	W_HLL[0] = (lmbda_R * W_R[0] - lmbda_L * W_L[0] + W_flux_L[0] - W_flux_R[0]) / (lmbda_R - lmbda_L);
	W_HLL[1] = (lmbda_R * W_R[1] - lmbda_L * W_L[1] + W_flux_L[1] - W_flux_R[1]) / (lmbda_R - lmbda_L);
	assert(W_HLL[0] >= 0.0);
}


int sgn(double x) {
	if (x > 0)
		return 1;
	else if (x < 0)
		return -1;
	else
		return 0;
}


double h_cutoff(double h_L, double h_R, double dx) {
	/** 
	 * Cutoff function to ensure the scheme is well-balanced when according to smooth steady states. However, it means that the source term approximation does not vanish when the topography is flat.
	 * */
	if (fabs(h_R - h_L) <= C_CUTOFF * dx)
		return h_R - h_L;
	else
		return sgn(h_R - h_L) * C_CUTOFF * dx;
}


double compute_S_dx(double h_L, double h_R, double Z_L, double Z_R, double dx) {
	double Z_diff  = Z_R - Z_L;
	if ( (h_L < ZERO_THRESH) && (h_R < ZERO_THRESH) )
		return 0.0;
	else if ( (h_L) < ZERO_THRESH || (h_R < ZERO_THRESH) )
		return -0.5 * GR * Z_diff * (h_L + h_R);
	else {
		double h_C = h_cutoff(h_L, h_R, dx);
		double a = -GR * Z_diff * 2 * h_L * h_R / (h_L + h_R);
		double b = 0.5 * GR * pow(h_C, 3) / (h_L + h_R);
		return a + b;
	}

}


void WB_solver(double ** W, double ** W_forward, double ** W_L_stars, double ** W_R_stars, double * lmbda_neg, double * lmbda_pos, double * Z, int nx, double dx, double T_stop, double * xs, double k, double gamma_of_k, double sig_theta, FILE * CURVE_DATA, FILE * TOP_DATA, FILE * TIMES) {

	/**
	 * This is the well-balanced solver within which we apply the conservative formula. To do this we need to compute 	the time increment, the wave speeds and the intermediate states at each time iterate. We apply the conservative formula until T_stop is exceeded.
	 * 	nx :: number of cells
		The space discretisation consists in cells (x_{i - 1/2}, x_{i + 1/2}) of volume dx and centered at x_i = x_{i - 1/2} + dx / 2
	*/

	int n = 0;
	double t = 0, dt = 0.1;
	double Z_L, Z_R, h_L, h_R, q_L, q_R, h_HLL, q_HLL, S_dx, q_star, S_dx_by_alph, h_L_star, h_R_star, lmbda_L, lmbda_R, lmbda_jump;
	double * W_L = (double *) malloc(2 * sizeof(double));
	double * W_R = (double *) malloc(2 * sizeof(double));
	double * W_HLL = (double *) malloc(2 * sizeof(double));
	double * W_flux_L = (double *) calloc(2, sizeof(double));
	double * W_flux_R = (double *) calloc(2, sizeof(double));
	double height = 0.2, centre = 10.0; // Temporary. These will be the random walk
	gen_Z_drain(nx, xs, Z, height, centre);
	for (int j = 0; j < nx + 2; j++)
		memcpy(W_forward[j], W[j], 2 * sizeof(double));
	
	while (t <= T_stop) {

		for (int j = 0; j < nx + 2; j++)
			memcpy(W[j], W_forward[j], 2 * sizeof(double));

		/* Assign the values of the ghost cells wrt the Neumann boundary conditions */
		printf("t = %lf\n", t);
		prescribe_left_bcs(W, dx, t);
		prescribe_right_bcs(W, nx);
		output_data(W, Z, nx, t, sig_theta, CURVE_DATA, TOP_DATA, TIMES);

		/* Iterate along each boundary (Riemann problem) */
		for (int j = 0; j < nx + 1; j++) {

			/* Define the left and right states of the Riemann problem */
			zero_cutoff(W[j], W[j + 1]);
			memcpy(W_L, W[j], 2 * sizeof(double));
			memcpy(W_R, W[j + 1], 2 * sizeof(double));
			Z_L = Z[j], Z_R = Z[j + 1], h_L = W_L[0], h_R = W_R[0], q_L = W_L[1], q_R = W_R[1];

			/* Compute the wave speeds and the HLL terms */
			compute_wave_speeds(W_L, W_R, lmbda_neg, lmbda_pos, j);
			lmbda_L = lmbda_neg[j], lmbda_R = lmbda_pos[j], lmbda_jump = lmbda_R - lmbda_L;
			compute_W_HLL(W_HLL, lmbda_neg[j], lmbda_pos[j], W_L, W_R, W_flux_L, W_flux_R);
			h_HLL = W_HLL[0], q_HLL = W_HLL[1];

			/* Compute the intermediate states */
			S_dx = compute_S_dx(h_L, h_R, Z_L, Z_R, dx);


		}


		t += dt, n++;

	}
	fprintf(CURVE_DATA, "%d\n", n);

}




























