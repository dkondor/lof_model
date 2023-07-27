/*
 * pop_model3.cpp -- -- numerically solve analytic population models with
 * 	three components
 * 
 * Copyright 2023 Daniel Kondor <kondor@csh.ac.at>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */


#ifndef HAVE_INLINE
#define HAVE_INLINE
#endif

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <random>
#include <time.h>
#include <stdexcept>



struct params {
	enum class type {
		NO_TYPE,
		LOF /* landscape of fear model */
	};
	type t = type::NO_TYPE;
	double r = 0.01; // base population growth rate
	double K = 1000.0; // total environment carrying capacity
	double a = 0.1; // warfare growth rate
	double b = 0.01; // warfare decline rate
	double c = 0.01; // migration rate from the general environment to the reguge
	double d = 0.01; // migration rate from the refuge to the general environment
	double rho = 0.1; // ratio for the refuge capacity
	
	static bool try_parse_type(const char* x, type& t) {
		if(!x) return false;
		
		if(x[0] == 'l' || x[0] == 'L') {
			t = type::LOF;
			return true;
		}
		return false;
	}
	
	
	int dfunc2(const double y1[3], double dydt[3]) const {
		double x = y1[0];
		double y = y1[1];
		double z = y1[2];
		
		switch(t) {
			case type::LOF:
			{
				/* x: warfare level, y: population, z: refuge population */
				double R = rho * K;
				dydt[0] = x * (a * y / K - b);
				dydt[1] = r * y * (1.0 - y / K) - c * x * y + d * z;
				dydt[2] = r * z * (1.0 - z / R) + c * x * y - d * z;
				break;
			}
			default:
				return GSL_EBADFUNC;
		}
		
		return GSL_SUCCESS;
	}
};


int dfunc(double t, const double y1[], double dydt[], void* params1) {
	const params& pars = *(const params*)params1;
	return pars.dfunc2(y1, dydt);
}


/* run one instance of the simulation from the given initial conditions using the given parameters */
int do_one_run(const params& pars, const double y0[3], double tmax, double dt, bool fixed_step, double noise_d[3],
		bool discrete_noise, FILE* fout, std::mt19937_64& rng) {
	
	std::normal_distribution<double> nd[3];
	bool have_noise = false;
	for(int i = 0; i < 3; i++) if(noise_d[i] > 0.0) {
		nd[i] = std::normal_distribution<double>(0.0, noise_d[i]);
		have_noise = true;
	}
	
	gsl_odeiv2_system sys = {dfunc, nullptr, 3, const_cast<params*>(&pars)};
	gsl_odeiv2_driver* d = nullptr;
	if(fixed_step <= 0.0) d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);
	
	double t = 0.0;
	double y1[] = {y0[0], y0[1], y0[2]};
	double dydt[] = {0.0, 0.0, 0.0};
	
	if(fout) fprintf(fout, "%f\t%f\t%f\t%f\n", t, y1[0], y1[1], y1[2]);
	
	// keep track of the derivatives in dydt
	pars.dfunc2(y1, dydt);
	
	while(t < tmax) {
		// add noise sampled from a Gaussian to either or both variables
		for(int i = 0; i < 3; i++) if(noise_d[i] > 0.0) y1[i] = std::max(0.0, y1[i] + nd[i](rng) * y1[i]);
		
		if(fixed_step) {
			for(int i = 0; i < 3; i++) if(discrete_noise && dydt[i] != 0.0) {
				bool xn = false;
				if(dydt[i] < 0.0) {
					xn = true;
					dydt[i] *= -1.0;
				}
				std::exponential_distribution<double> ed1(1.0 / dydt[i]);
				dydt[i] = ed1(rng);
				if(xn) dydt[i] *= -1.0;
			}
			for(int i = 0; i < 3; i++) {
				y1[i] += dydt[i];
				if(y1[i] < 0.0) y1[i] = 0.0; // all variables are constrained to >= 0
			}
			t += dt;
		}
		else {
			if(have_noise) gsl_odeiv2_driver_reset(d);
			double t1 = t + dt;
			int r = gsl_odeiv2_driver_apply(d, &t, t1, y1);
			if(r != GSL_SUCCESS) {
				fprintf(stderr, "Error advancing the solution (%d)!\n", r);
				gsl_odeiv2_driver_free(d);
				return r;
			}
		}
		
		if(fout) fprintf(fout, "%f\t%f\t%f\t%f\n", t, y1[0], y1[1], y1[2]);
		
		if(fixed_step) pars.dfunc2(y1, dydt);
	}
	
	if(d) gsl_odeiv2_driver_free(d);
	
	return 0;
}



int main(int argc, char **argv) {
	params pars;
	double x = 0.1;
	double y = 0.1 * pars.K;
	double z = 0.1 * y;
	double tmax = 500.0; // maximum time
	double dt = 1.0; // time step for output
	bool fixed_step = false; // if true, instead of running a solver,
	// use discrete-time steps (of e.g. 1 year), essentially evolving
	// the system as a difference equation
	double noise_d[] = {0.0, 0.0, 0.0}; // add noise from a Gaussian with this relative scale to the variables each year
	bool discrete_noise = false; // add noise in the discrete model directly to the yearly differences
	uint64_t seed = time(0);
	
	for(int i = 1; i < argc; i++) if(argv[i][0] == '-') switch(argv[i][1]) {
		case 'a':
			pars.a = atof(argv[i+1]);
			i++;
			break;
		case 'b':
			pars.b = atof(argv[i+1]);
			i++;
			break;
		case 'c':
			pars.c = atof(argv[i+1]);
			i++;
			break;
		case 'd':
			pars.d = atof(argv[i+1]);
			i++;
			break;
		case 'r':
			pars.r = atof(argv[i+1]);
			i++;
			break;
		case 'R':
			pars.rho = atof(argv[i+1]);
			i++;
			break;
		case 'K':
			pars.K = atof(argv[i+1]);
			i++;
			break;
		case 't':
			if(!params::try_parse_type(argv[i+1], pars.t)) {
				fprintf(stderr, "Unknown model type: %s!\n", argv[i+1]);
				return 1;
			}
			i++;
			break;
		case 'x':
			x = atof(argv[i+1]);
			i++;
			break;
		case 'y':
			y = atof(argv[i+1]);
			i++;
			break;
		case 'z':
			z = atof(argv[i+1]);
			i++;
			break;
		case 'T':
			tmax = atof(argv[i+1]);
			i++;
			break;
		case 'D':
			dt = atof(argv[i+1]);
			i++;
			break;
		case 'F':
			fixed_step = true;
			break;
		case 'N':
			switch(argv[i][2]) {
				case 'x':
					noise_d[0] = atof(argv[i+1]);
					i++;
					break;
				case 'y':
					noise_d[1] = atof(argv[i+1]);
					i++;
					break;
				case 'z':
					noise_d[2] = atof(argv[i+1]);
					i++;
					break;
				case 0:
					noise_d[0] = atof(argv[i+1]);
					noise_d[1] = noise_d[0];
					noise_d[2] = noise_d[0];
					i++;
					break;
				case 'd':
					discrete_noise = true;
					break;
				default:
					fprintf(stderr, "Unknown parameter: %s!\n", argv[i]);
					return 1;
			}
			break;
		case 's':
			seed = strtoul(argv[i+1], nullptr, 10);
			i++;
			break;
		default:
			fprintf(stderr, "Unknown parameter: %s!\n", argv[i]);
			return 1;
	}
	else {
		fprintf(stderr, "Unknown parameter: %s!\n", argv[i]);
		return 1;
	}
	
	if(pars.t == params::type::NO_TYPE) {
		fprintf(stderr, "Type of system not given!\n");
		return 1;
	}
	
	std::mt19937_64 rng(seed);
	FILE* fout = stdout;
	
	double y0[] = {x, y, z};
	
	return do_one_run(pars, y0, tmax, dt, fixed_step, noise_d, discrete_noise, fout, rng);
}


