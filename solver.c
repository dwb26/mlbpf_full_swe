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
const double m_glob = 0.5;
const double M_glob = 1e-10;


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


double compute_q_star(double q_HLL, double S_dx, double lmbda_L, double lmbda_R) {
	/**
	 * This is the discharge component of the intermediate states. Note that there is only one since we set q_L* = q_R*
	 */
	return q_HLL * S_dx / (lmbda_R - lmbda_L);
}


double compute_S_dx_by_alph(double h_L, double h_R, double Z_L, double Z_R, double S_dx, double q_star) {

	if ( (h_L < ZERO_THRESH) && (h_R < ZERO_THRESH) )
		return 0.0;
	else if ( (h_L < ZERO_THRESH) || (h_R < ZERO_THRESH) )
		return -(Z_R - Z_L);
	else {
		double alph = -(q_star * q_star) / (h_L * h_R) + 0.5 * GR * (h_L + h_R);
		return S_dx / alph;
	}

}


void compute_h_stars(double * h_stars, double h_HLL, double lmbda_L, double lmbda_R, double S_dx_by_alph, double h_L, double h_R) {
	/**
	 * Computes the lefy and right intermediate height states
	 */
	double eps = h_L < h_R ? h_L : h_R;
	eps = eps < h_HLL ? eps: h_HLL;
	double lmbda_jump = lmbda_R - lmbda_L;
	double lmbda_R_rat = lmbda_R / lmbda_L;
	double lmbda_L_rat = lmbda_L / lmbda_R;
	assert(eps >= 0);

	/* Compute the left intermediate height states */
	double term = h_HLL - lmbda_R * S_dx_by_alph / lmbda_jump;
	double a = term > eps ? term : eps;
	double b = (1 - lmbda_R_rat) * h_HLL + lmbda_R_rat * eps;
	h_stars[0] = a < b ? a : b;

	/* Compute the right intermediate height states */
	term = h_HLL - lmbda_L * S_dx_by_alph / lmbda_jump;
	a = term < eps ? term : eps;
	b = (1 - lmbda_L_rat) * h_HLL + lmbda_L_rat * eps;
	h_stars[1] = a < b ? a : b;

	double diff = lmbda_R * h_stars[1] - lmbda_L * h_stars[0] - lmbda_jump * h_HLL;
	assert(diff <= 2.0 * ZERO_THRESH);

}


double compute_timestep(double lmbda_L_max, double lmbda_R_max, double dx) {
	/**
	 * Compute the time increment using the maximum magnitude wavespeed and a CFL type condition
	 */
	double lmbda_max = fabs(lmbda_L_max) > lmbda_R_max ? fabs(lmbda_L_max) : lmbda_R_max;
	assert (lmbda_max > 0);
	double dt = 0.5 * dx / lmbda_max;
	assert (dt / dx <= 1);
	return dt;

}


double compute_theta(double ** W, double * Z, double dx, int j) {
	/**
	 * This is a parameter that describes the convexity of the reconstruction. Depends only on h and q.
	 */

	double h_L, h_M, h_R, q_L, q_M, q_R, Z_L, Z_M, Z_R, S_dx_LM, S_dx_MR, psi_LM, psi_MR, a, b, varphi, theta;
	double * W_L = (double *) malloc(2 * sizeof(double));
	double * W_M = (double *) malloc(2 * sizeof(double));
	double * W_R = (double *) malloc(2 * sizeof(double));

	memcpy(W_L, W[j - 1], 2 * sizeof(double));
	memcpy(W_M, W[j], 2 * sizeof(double));
	memcpy(W_R, W[j + 1], 2 * sizeof(double));
	Z_L = Z[j - 1], Z_M = Z[j], Z_R = Z[j + 1];
	h_L = W_L[0], q_L = W_L[1];
	h_M = W_M[0], q_M = W_M[1];
	h_R = W_R[0], q_R = W_R[1];
	S_dx_LM = compute_S_dx(h_L, h_M, Z_L, Z_M, dx);
	S_dx_MR = compute_S_dx(h_M, h_R, Z_M, Z_R, dx);

	/* Compute the boundary psi terms */
	if (h_L < ZERO_THRESH) {
		if (h_M < ZERO_THRESH)
			psi_LM = -S_dx_LM;
		else
			psi_LM = q_M * q_M / h_M + 0.5 * GR * h_M * h_M - S_dx_LM;
	}
	else if (h_M < ZERO_THRESH) {
		if (h_L < ZERO_THRESH)
			psi_LM = -S_dx_LM;
		else
			psi_LM = -q_L * q_L / h_L - 0.5 * GR * h_L * h_L - S_dx_LM;
	}
	else
		psi_LM = q_M * q_M / h_M - q_L * q_L / h_L + 0.5 * GR * (h_M * h_M - h_L * h_L) - S_dx_LM;

	if (h_M < ZERO_THRESH) {
		if (h_R < ZERO_THRESH)
			psi_MR = -S_dx_MR;
		else
			psi_MR = q_R * q_R / h_R + 0.5 * GR * h_R * h_R - S_dx_MR;
	}
	else if (h_R < ZERO_THRESH) {
		if (h_M < ZERO_THRESH)
			psi_MR = -S_dx_MR;
		else
			psi_MR = -q_M * q_M / h_M - 0.5 * GR * h_M * h_M - S_dx_MR;
	}
	else
		psi_MR = q_R * q_R / h_R - q_M * q_M / h_M + 0.5 * GR * (h_R * h_R - h_M * h_M) - S_dx_MR;


	/* Compute the centred-cell varphi term */
	a = sqrt((q_M - q_L) * (q_M - q_L) + psi_LM * psi_LM);
	b = sqrt((q_R - q_M) * (q_R - q_M) + psi_MR * psi_MR);
	varphi = a + b;

	/* Compute the centred-cell theta term */
	if (varphi < m_glob * dx)
		theta = 0.0;
	else if (varphi <= M_glob * dx)
		theta = (varphi - m_glob * dx) / (M_glob * dx - m_glob * dx);
	else
		theta = 1.0;

	assert( (theta >= 0) && (theta <= 1) );

	free(W_L);
	free(W_M);
	free(W_R);

	return theta;

}


double minmod(double a, double b) {

	if ( (fabs(a) < fabs(b)) && (a * b > 0) )
		return a;
	else if ( (fabs(b) < fabs(a)) && (a * b > 0) )
		return b;
	else
		return 0.0;

}


double compute_slope(double * v, double dx, int j) {
	/**
	 * Computes the slope of the linear reconstruction corresponding to the same cell w_M is considered on.
	 */
	double v_L, v_M, v_R;
	v_L = v[j - 1], v_M = v[j], v_R = v[j + 1];
	return minmod((v_R - v_M) / dx, (v_M - v_L) / dx);

}


void MUSCL_forward_solution(double ** W_forward, double ** W, double dx, double dt, int nx, double * Z, double * lmbda_neg, double * lmbda_pos, double ** W_L_stars, double ** W_R_stars) {


	/* -------------------------------------------------------------------------------------------------------------- */
	/*  																											  */
	/* Array definition 																							  */
	/*  																											  */
	/* -------------------------------------------------------------------------------------------------------------- */
	double h_LM, q_LM, h_ML, q_ML, h_MR, q_MR, h_RM, q_RM;
	double * h = (double *) malloc((nx + 2) * sizeof(double));
	double * q = (double *) malloc((nx + 2) * sizeof(double));
	double * h_plus_Z = (double *) malloc((nx + 2) * sizeof(double));
	for (int j = 0; j < nx + 2; j++) {
		h[j] = W[j][0];
		q[j] = W[j][1];
		h_plus_Z[j] = h[j] + Z[j];
	}
	double * thetas = (double *) calloc((nx + 2), sizeof(double));
	double ** slopes = (double **) malloc(3 * sizeof(double *));
	for (int m = 0; m < 3; m++)
		slopes[m] = (double *) calloc((nx + 2), sizeof(double));
	double ** W_tilde = (double **) malloc((nx + 2) * sizeof(double *));
	double ** gs = (double **) malloc((nx + 2) * sizeof(double *));
	for (int m = 0; m < nx + 2; m++) {
		W_tilde[m] = (double *) malloc(2 * sizeof(double));
		memcpy(W_tilde[m], W[m], 2 * sizeof(double));
		gs[m] = (double *) malloc(2 * sizeof(double));
	}
	double * h_MUSCL = (double *) malloc(4 * sizeof(double));
	double * q_MUSCL = (double *) malloc(4 * sizeof(double));
	double * Z_MUSCL = (double *) malloc(4 * sizeof(double));
	double * W_LM = (double *) malloc(2 * sizeof(double));
	double * W_ML = (double *) malloc(2 * sizeof(double));
	double * W_MR = (double *) malloc(2 * sizeof(double));
	double * W_RM = (double *) malloc(2 * sizeof(double));
	double * h_tilde = (double *) malloc((nx + 2) * sizeof(double));
	double * q_tilde = (double *) malloc((nx + 2) * sizeof(double));
	double * h_plus_Z_tilde = (double *) malloc((nx + 2) * sizeof(double));



	/* -------------------------------------------------------------------------------------------------------------- */
	/*  																											  */
	/* Basic MUSCL reconstruction																					  */
	/*  																											  */
	/* -------------------------------------------------------------------------------------------------------------- */
	for (int j = 1; j < nx + 1; j++) {

		/* Compute the theta parameters for the convex MUSCL reconstruction */
		thetas[j] = compute_theta(W, Z, dx, j);

		/* Compute the respective slopes */
		slopes[0][j] = compute_slope(h, dx, j);
		slopes[0][j] = compute_slope(q, dx, j);
		slopes[0][j] = compute_slope(h_plus_Z, dx, j);

	}

	/* Apply the forward formula */






	/* -------------------------------------------------------------------------------------------------------------- */
	/*  																											  */
	/* Stable forward time MUSCL reconstruction																		  */
	/*  																											  */
	/* -------------------------------------------------------------------------------------------------------------- */





	free(h);
	free(q);
	free(h_plus_Z);
	free(thetas);
	free(slopes);
	free(W_tilde);
	free(gs);
	free(h_MUSCL);
	free(q_MUSCL);
	free(Z_MUSCL);
	free(W_LM);
	free(W_ML);
	free(W_MR);
	free(W_RM);
	free(h_tilde);
	free(q_tilde);
	free(h_plus_Z_tilde);

}


void WB_solver(double ** W, double ** W_forward, double ** W_L_stars, double ** W_R_stars, double * lmbda_neg, double * lmbda_pos, double * Z, int nx, double dx, double T_stop, double * xs, double k, double gamma_of_k, double sig_theta, FILE * CURVE_DATA, FILE * TOP_DATA, FILE * TIMES) {

	/**
	 * This is the well-balanced solver within which we apply the conservative formula. To do this we need to compute 	the time increment, the wave speeds and the intermediate states at each time iterate. We apply the conservative formula until T_stop is exceeded.
	 * 	nx :: number of cells
		The space discretisation consists in cells (x_{i - 1/2}, x_{i + 1/2}) of volume dx and centered at x_i = x_{i - 1/2} + dx / 2
	*/

	int n = 0;
	double t = 0, dt = 0.1;
	double Z_L, Z_R, h_L, h_R, q_L, q_R, h_HLL, q_HLL, S_dx, q_star, S_dx_by_alph, lmbda_L, lmbda_R;
	double lmbda_L_max, lmbda_R_max;
	double * W_L = (double *) malloc(2 * sizeof(double));
	double * W_R = (double *) malloc(2 * sizeof(double));
	double * W_HLL = (double *) malloc(2 * sizeof(double));
	double * W_flux_L = (double *) calloc(2, sizeof(double));
	double * W_flux_R = (double *) calloc(2, sizeof(double));
	double * h_stars = (double *) calloc(2, sizeof(double));
	double height = 0.2, centre = 10.0; // Temporary. These will be the random walk
	gen_Z_drain(nx, xs, Z, height, centre);
	for (int j = 0; j < nx + 2; j++)
		memcpy(W_forward[j], W[j], 2 * sizeof(double));
	
	while (t <= T_stop) {

		lmbda_L_max = 0.0, lmbda_R_max = 0.0;
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
			lmbda_L = lmbda_neg[j], lmbda_R = lmbda_pos[j];
			lmbda_L_max = lmbda_L_max < lmbda_L ? lmbda_L_max : lmbda_L;
			lmbda_R_max = lmbda_R_max > lmbda_R ? lmbda_R_max : lmbda_R;
			compute_W_HLL(W_HLL, lmbda_neg[j], lmbda_pos[j], W_L, W_R, W_flux_L, W_flux_R);
			h_HLL = W_HLL[0], q_HLL = W_HLL[1];

			/* Compute the intermediate states */
			S_dx = compute_S_dx(h_L, h_R, Z_L, Z_R, dx);
			q_star = compute_q_star(q_HLL, S_dx, lmbda_L, lmbda_R);
			S_dx_by_alph = compute_S_dx_by_alph(h_L, h_R, Z_L, Z_R, S_dx, q_star);
			compute_h_stars(h_stars, h_HLL, lmbda_L, lmbda_R, S_dx_by_alph, h_L, h_R);
			W_L_stars[j][0] = h_stars[0], W_L_stars[j][1] = q_star;
			W_R_stars[j][0] = h_stars[1], W_R_stars[j][1] = q_star;

		}
		dt = compute_timestep(lmbda_L_max, lmbda_R_max, dx);
		MUSCL_forward_solution(W_forward, W, dx, dt, nx, Z, lmbda_neg, lmbda_pos, W_L_stars, W_R_stars);
		t += dt, n++;

	}
	fprintf(CURVE_DATA, "%d\n", n);

}




























