/*
 * cell.hpp -- a cell represents a local area where demographic simulations take place
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

#ifndef CELL_HPP
#define CELL_HPP

#include <stdint.h>
#include <cmath>
#include <algorithm>
#include <random>
#include <stdexcept>

/**
 * A cell is the basic unit of the simulation. Demographics in cells are
 * simulated using a logistic equation (individual people are not considered).
 * This class contains the variables defining the current state of a cell
 * and basic functionality to update the local population.
 * 
 * Cells are organized in a landscape in the container classes:
 * cellnet and cellgrid (in their respective header)
 * Typically they are not used directly.
 */
struct cell {
	/* possible states of a cell */
	enum class state {
		EMPTY = 0, /* cell is empty, no agricultural village present */
		FARMING, /* farming village is present */
		RAIDERS /* group of warriors who live by raiding nearby farming villages */
		/* additional states could represent foragers or warrior groups */
	};
	
	uint32_t N = 0; // current population
	double K = 0.0; // carrying capacity (number of people) -- this can change in time generally, or multiplied by factors dependant on local conditions
	state st = state::EMPTY; // current state (note: st == EMPTY implies N == 0 and vice versa)
	
	double conflict_level = 0.0; // current perception of the level of conflict around this cell
	double defensibility = 0.0; // (perceived) defensibility of this cell
	// note: if conflict_level > defensibility, this cell is abandoned
	
	// some simulation parameters (note: these are not stored for each cell, but defined here since they belong here logically)
	struct cell_pars {
		double r = 0.0135; // birth rate (per year)
		double delta = -1.0; // population collapse if larger than carrying capacity (< 0: no collapse, use logistic model)
	};
	
	// get the expected population after one step
	static double get_expected_pop(double t, double K, double N, const cell_pars& pars) {
		double target = K; /* note: actual carrying capacity could be a function of K */
		double r = pars.r;
		double delta = pars.delta;
		
		/* population collapses to delta * target */
		if(N > target && delta > 0.0) return delta * target;
		else {
			/* logistic growth / decline */
			double tmp = N;
			tmp = (target - tmp) / tmp;
			tmp *= exp(-1.0 * r * t);
			tmp = target / (1 + tmp);
			/* tmp is the expected population based on the solution of the 
			 * logistic differential equation */
			return tmp;
		}
	}
	
	// advance the demographic model by t time (years); above parameters are given as a variable
	template<class RNG>
	void step(double t, RNG& rng, const cell_pars& pars = default_pars());
	
	// default cell parameters
	static constexpr cell_pars default_pars() {
		return cell_pars();
	}
	
	/* cells can be created by specifying their carrying capacity */
	explicit cell(double K_ = 0.0) : K(K_) { }
};


template<class RNG>
void cell::step(double t, RNG& rng, const cell_pars& pars) {
	if(N == 0) throw std::runtime_error("cell::step(): empty cell found!\n");
	
	double tmp = get_expected_pop(t, K, N, pars);
	if(tmp < N && pars.delta > 0.0) {
		/* population collapses to delta * target */
		N = std::round(tmp);
		if(N > K) N = std::floor(K);
	}
	else {
		/* logistic growth / decline;
		 * tmp is the expected population based on the solution of the 
		 * logistic differential equation
		 * we calculate actual growth based on a Poisson distribution */
		double gr = tmp - N;
		bool growth = (gr > 0.0) ? true : false;
		gr = fabs(gr);
		if(gr > 0.0) {
			std::poisson_distribution<unsigned int> pd(gr);
			unsigned int N1 = pd(rng);
			if(growth) {
				N += N1;
				if(N > K) N = std::floor(K);
			}
			else N -= std::min(N, N1);
		}
	}
}

#endif

