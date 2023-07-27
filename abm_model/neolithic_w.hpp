/*
 * neolithic_w.hpp -- simulation of neolithic farming communities with
 * 	simple model of warfare and conflicts
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



#ifndef NEOLITHIC_W_HPP
#define NEOLITHIC_W_HPP

#include "cell.hpp"
#include "neolithic_base.hpp"
#include <vector>
#include <unordered_set>
#include <utility>
#include <memory>
#include <time.h>

template<class cell_collection>
class neolithic_w : public nbase2<cell_collection> {
protected:
	/* main parameters -- group migration */
	bool stationary_raiders_split = true; /* if true, even for stationary raiders, the origin cell is split for the first conflict */
	bool defenders_create_raiders = false; /* if true, raiders are created if a cell is successfully defended from an attack as well */
	bool raiders_do_revert = true; /* if true, raiders revert to farming if they cannot find a suitable target, or if their attack fails */
	bool raider_dist_real_pop = false; /* if true, raiders consider the real population of target cells for attack */
	double raider_dist = 5.0; /* characteristic distance for raiders */
	double raider_max_dist = -1.0; /* separate variable controlling the maximum distance raiders can look for targets */
	double attack_success_prob = 0.5; /* probability of an attacker being successfull */
	double survivor_ratio = 0.8; /* rate of survivors after (1) a successfull attack by non-mobile raiders; (2) an unsuccessful attack of raiders */
	double praiders = 0.5; /* chance of attackers converting to raiders */
	double raiders_revert_rate = 0.0; /* rate of raiders reverting (only for the new stationary warrior model) */
	double pattack = 1.0; /* rate at which raiders attack in the new stationary model */
	unsigned int raider_revert_threshold = 0; /* raiders revert back to farming if their population falls below this threshold (after an attack by other raiders) */
	
	/* landscape of fear model parameters */
	bool discrete_conflict = false; /* if true, conflict level in a cell jumps discretely from 0 to 1 in case of a conflict */
	double conflict_effect_base = 1.0; /* base effect of conflict on a cell (conflict_level increases by this amount or 
		chance of jumping to 1 if discrete_conflict = true) */
	bool conflict_effect_only_on_success = true; /* if true, conflicts only have an effect if attackers win */
	double refuge_bonus = 0.0; /* bonus for refuge cells for defense */
	double refugee_threshold = 0.0; /* a cell is abandoned if its population falls below this share of K (after an attack) */
	bool raiders_tolerate_conflict = false; /* if true, raiders do not abandon a cell with conflict */
	bool raiders_tolerate_self = false; /* if true, raiders are not affected by conflicts initiated by themselves */
	bool raiders_attack_refuges = true; /* if false, raiders do not attack cells with defensibility == 1.0 */
	double conflict_dist = -1.0; /* characteristic distance for conflict effects (< 0 means use raider_dist) */
	/* note: raiders_revert_rate is used for decreaseing the conflict levels as well */
	double refugee_dist = -1.0; /* characteristic distance for refugees (< 0 means use raider_dist) */
	bool refugee_target_def = false; /* refugees only flee to defensible locations */
	bool only_instant_conflict = false; /* conflicts do not change cells' "conflict_level", only have immediate effects */
	bool limit_split_by_raiders = false; /* if a (farmer) cells' neighbors include raiders, it does not initiate split-offs at all */
	double scale_refuge_targets = -1.0; /* relative weight for attacking refuges */
	double limit_split_by_raiders_dist = -1.0; /* distance of LoF effect for limiting split-offs */
	
	using nbase<cell_collection>::rng;
	using nbase<cell_collection>::pempty;
	using nbase<cell_collection>::use_global_neighbors;
	using nbase<cell_collection>::group_mig_dist;
	using nbase<cell_collection>::group_mig_pow;
	using nbase<cell_collection>::cph;
	using nbase<cell_collection>::group_mig_prob;
	using nbase<cell_collection>::use_prob_helper;
	using group_mig_order = typename nbase<cell_collection>::group_mig_order;
	using cell_prob_helper_type = typename nbase<cell_collection>::cell_prob_helper_type;
	
public:
	using nbase<cell_collection>::g;
	// using nbase<cell_collection>::orig_K; -- not sure if this is required
	
	neolithic_w() : nbase2<cell_collection>(false) {
		use_prob_helper = cell_prob_helper_type::BTCB;
	}

	unsigned int attack_success = 0; /* count the number of successful and failed attacks in each time step */
	unsigned int attack_fail = 0;
	unsigned int raiders_created = 0; /* count the number of raider cells created and reverted */
	unsigned int raiders_reverted = 0;
	unsigned int total_killed = 0; /* total population directly killed in conflicts */
	unsigned int total_dec1 = 0; /* total population decrease due to natural causes (when N < K) */
	unsigned int total_dec2 = 0; /* total population decrease due to overcrowding (when N > K) */
	unsigned int total_refugees = 0; /* total number of people moving to different cells as refugees */
	unsigned int cancelled_refugee_cells = 0; /* cells where refugees could not find any target */
	unsigned int cancelled_refugees = 0; /* number of refugees who could not find a target */
	
	/* add potential parse command line options for parameters in this class to parser */
	void add_options(option_parser& op) override {
		nbase<cell_collection>::add_options(op);
		
		op.add_numeric_option("a", &attack_success_prob);
		op.add_bool_option("Rc", &defenders_create_raiders);
		op.add_bool_option("Rr0", &raiders_do_revert);
		op.add_numeric_option("RR", &raiders_revert_rate);
		op.add_numeric_option("Rp", &praiders);
		op.add_bool_option("RP", &raider_dist_real_pop);
		op.add_numeric_option("RA", &pattack);
		op.add_numeric_option("Rs", &survivor_ratio);
		op.add_numeric_option("R", &raider_dist);
		op.add_numeric_option("Rm", &raider_max_dist);
		op.add_numeric_option("RT", &raider_revert_threshold);
		op.add_bool_option("RS", &stationary_raiders_split, false);
		op.add_option("b", [this] (int, int, char**) {
			use_prob_helper = cell_prob_helper_type::CDF;
			return 0;
		});
		op.add_bool_option("LD", &discrete_conflict);
		op.add_numeric_option("Lb", &conflict_effect_base);
		op.add_bool_option("LR", &raiders_tolerate_conflict);
		op.add_bool_option("LRS", &raiders_tolerate_self);
		op.add_bool_option("La", &raiders_attack_refuges, false);
		op.add_bool_option("Ls", &conflict_effect_only_on_success, false);
		op.add_numeric_option("LB", &refuge_bonus);
		op.add_numeric_option("Lc", &conflict_dist);
		op.add_numeric_option("Lt", &refugee_threshold);
		op.add_bool_option("Lrd", &refugee_target_def);
		op.add_numeric_option("Lrm", &refugee_dist);
		op.add_bool_option("LI", &only_instant_conflict);
		op.add_bool_option("LS", &limit_split_by_raiders);
		op.add_numeric_option("LSd", &limit_split_by_raiders_dist);
		op.add_numeric_option("LA", &scale_refuge_targets);
		op.add_post_callback([this]() {
			if(this->def_base < 0.0) this->def_base = std::min(conflict_effect_base, 1.0);
			if(conflict_dist <= 0.0) conflict_dist = raider_dist;
			if(refugee_dist <= 0.0) refugee_dist = raider_dist;
			if(scale_refuge_targets < 0.0) scale_refuge_targets = raiders_attack_refuges ? 1.0 : 0.0;
			else if(scale_refuge_targets > 0.0) raiders_attack_refuges = true;
			else raiders_attack_refuges = false;
			if(limit_split_by_raiders_dist > 0.0) limit_split_by_raiders = true;
			return true;
		});
	}
	
	std::vector<group_mig_order> orders;
	unsigned int orders_rejected = 0; // number of cells where a split-off was rejected due to not finding a target
	std::vector<typename cell_collection::point> raiders_revert; // aggressor cells that revert to being peaceful
	
	struct refugee_order {
		typename cell_collection::point src; // source cell
		double r; // ratio of refugees
		unsigned int seq;
		
		static void sort(std::vector<refugee_order>& orders) {
			std::sort(orders.begin(), orders.end(),
				[](const refugee_order& o1, const refugee_order& o2) {
					return o1.seq < o2.seq;
				});
		}
	};
	std::vector<refugee_order> refugees; // cells where (part) of the population runs away as refugees
	
	/* Create a set of "orders", i.e. migration choices for split-off groups,
	 * to be carried out later. */
	void create_group_orders() {
		orders.clear();
		raiders_revert.clear();
		
		this->ensure_omp_threads();
		
		/* init orders and reverted raiders to the maximum possible size */
		size_t tmpsize = g.size();
		orders.resize(tmpsize);
		raiders_revert.resize(tmpsize);
		size_t norders = 0;
		size_t nrr = 0;
		orders_rejected = 0;

#ifdef _OPENMP		
#pragma omp parallel for
#endif
		for(size_t j = 0; j < tmpsize; j++) {
			auto& rng2 = this->prng();
			const auto& pt1 = g.ptid(j);
			auto& c = g.at_ix(j);
			if(c.st == cell::state::EMPTY) continue;
			unsigned int N3 = 0; // faction size (in case of farming cell)
			
			/* estimate split-off groups and their targets, apply these changes */
			std::uniform_int_distribution<unsigned int> osd;
			
			std::vector<typename cell_collection::point> targets;
			std::vector<double> w;
			
			if(c.st == cell::state::FARMING) {
				double ps;
				N3 = this->cell_group_mig(j, rng2, ps);
				if(!N3) continue;
				
				if(limit_split_by_raiders) {
					/* check that any neighbor (up to raider_max_dist) is a raider */
					double max_dist;
					if(limit_split_by_raiders_dist > 0.0) max_dist = limit_split_by_raiders_dist;
					else max_dist = raider_max_dist > 0.0 ? raider_max_dist : raider_dist * 10.0;
					bool found = false;
					for(const auto& tmp1 : g.neighbors(pt1, max_dist)) {
						const auto& pt2 = tmp1.first;
						const auto& c2 = g.at(pt2);
						if(c2.st == cell::state::RAIDERS) {
							found = true;
							break;
						}
					}
					if(found) continue;
				}
				
				/* note: only farmers choose from the long-range distribution,
				 * raiders make separate choices among only farming cells (see below) */
				group_mig_order o;
				unsigned int max_retries = this->target_max_retries;
				if(this->scale_max_retries) {
					max_retries = (unsigned int)std::round(max_retries * ps);
					if(!max_retries) max_retries = 1;
				}
				if(this->choose_dst_helper(pt1, o, max_retries, rng2, N3)) {
					o.group_size = N3;
					o.seq = osd(rng2);
					size_t tmpo;
#ifdef _OPENMP		
#pragma omp atomic capture
#endif
					tmpo = norders++;
					orders[tmpo] = o;
				}
				else {
#ifdef _OPENMP
#pragma omp atomic update
#endif
					orders_rejected++;
				}
			}
			else {
				/* possibility for a raider cell reverting */
				bool revert = false;
				if(raiders_revert_rate > 0.0 && c.conflict_level <= c.defensibility) {
					std::uniform_real_distribution<double> rrd(0.0, 1.0);
					double x = rrd(rng2);
					if(x < raiders_revert_rate) revert = true;
				}
				
				if(!revert && pattack < 1.0) {
					/* in this case, raiders not necessarily attack in each year */
					std::uniform_real_distribution<double> da(0.0, 1.0);
					double x = da(rng2);
					if(x >= pattack) continue;
				}
				
				/* here, this is a raider cell, so we have to choose
				 * from a "local" distribution */
				double max_dist = raider_max_dist > 0.0 ? raider_max_dist : raider_dist * 10.0;
				if(!revert) for(const auto& tmp1 : g.neighbors(pt1, max_dist)) {
					const auto& pt2 = tmp1.first;
					const auto& c2 = g.at(pt2);
					double p1 = 1.0;
					if(c2.defensibility == 1.0) {
						if(!raiders_attack_refuges) continue; /* raiders do not attack defensible targets */
						p1 = scale_refuge_targets;
					}
					p1 *= exp(-1.0 * tmp1.second / raider_dist);
					if(raider_dist_real_pop) p1 *= c2.N;
					else p1 *= c2.K;
					
					if(p1 > 0.0) {
						targets.push_back(pt2);
						w.push_back(p1);
					}
				}
				
				if(w.size()) {
					/* raiding this way is only possible if there are possible target cells */
					std::discrete_distribution<unsigned int> td(w.cbegin(), w.cend());
					group_mig_order o;
					o.src = pt1;
					unsigned int tix = td(rng2);
					o.dst = targets[tix];
					o.group_size = 0;
					o.seq = osd(rng2);
					size_t tmpo;
#ifdef _OPENMP		
#pragma omp atomic capture
#endif
					tmpo = norders++;
					orders[tmpo] = o;
					targets.clear();
					w.clear();
				}
				else if(revert || raiders_do_revert) {
					size_t tmpr;
#ifdef _OPENMP		
#pragma omp atomic capture
#endif
					tmpr = nrr++;
					raiders_revert[tmpr] = pt1;
				}
			}
		}
		
		orders.resize(norders);
		raiders_revert.resize(nrr);
	}
	
	template<class CB>
	void apply_group_mig(CB&& cb) {
		std::unordered_set<typename cell_collection::point, typename cell_collection::pointhash> already_targeted;
		std::unordered_set<typename cell_collection::point, typename cell_collection::pointhash> refugees_pts;
		
		attack_success = 0;
		attack_fail = 0;
		raiders_created = 0;
		raiders_reverted = 0;
		total_killed = 0;
		total_refugees = 0;
		cancelled_refugee_cells = 0;
		cancelled_refugees = 0;
		
		/* process reverted raiders */
		for(const auto& pt1 : raiders_revert) {
			auto& c = g.at(pt1);
			c.st = cell::state::FARMING;
			raiders_reverted++;
		}
		
		/* sort the orders randomly */
		group_mig_order::sort_by_seq(orders);
		
		/* process migration events */
		std::uniform_real_distribution<double> pf1;
		std::uniform_int_distribution<unsigned int> osd;
		for(size_t i = 0; i < orders.size(); i++) {
			const auto& src = orders[i].src;
			const auto& dst = orders[i].dst;
			if(already_targeted.count(src)) {
				/* source cell was attacked already, simply skip this order */
				continue;
			}
			auto& c1 = g.at(src);
			auto& c2 = g.at(dst);
			if(already_targeted.count(dst) || 
				(c2.st == cell::state::EMPTY && c1.st == cell::state::RAIDERS) )  {
				/* target was already attacked, or it was originally a farming cell, but
				 * everyone migrated away -- cancel this order */
				if(c1.st == cell::state::RAIDERS && raiders_do_revert) {
					raiders_reverted++;
					c1.st = cell::state::FARMING;
				}
				continue;
			}
			already_targeted.insert(dst);
			
			uint32_t N3 = orders[i].group_size;
			cb(orders[i].src, orders[i].dst, N3);
			
			cell::state dsts = c2.st;
			
			bool conflict = false;
			bool migration = false;
			
			if(c2.st != cell::state::EMPTY) conflict = true;
			else migration = true;
			
			if(conflict) {
				bool success = false;
				double prob1 = attack_success_prob;
				if(c2.defensibility == 1.0) prob1 = std::max(0.0, prob1 - refuge_bonus);
				double x1 = (prob1 < 1.0) ? pf1(rng) : 0.0;
				if(x1 < prob1) {
					/* successful attack */
					attack_success++;
					success = true;
					if(c1.st == cell::state::FARMING) migration = stationary_raiders_split;
					if(!migration) {
						/* target is depopulated */
						uint32_t N2 = (unsigned int)std::round(survivor_ratio * c2.N);
						if(N2 < c2.N) {
							total_killed += (c2.N - N2);
							c2.N = N2;
						}
						if(!c2.N) c2.st = cell::state::EMPTY;
						else if(c2.st == cell::state::RAIDERS && c2.N < raider_revert_threshold) {
							c2.st = cell::state::FARMING;
						}
						else if(c2.N < c2.K * refugee_threshold) {
							/* target cell's inhabitant run away */
							if(refugees_pts.insert(dst).second) {
								refugee_order tmp2;
								tmp2.src = dst;
								tmp2.r = 1.0;
								tmp2.seq = osd(rng);
								refugees.push_back(tmp2);
							}
						}
						if(c1.st == cell::state::FARMING) {
							/* chance to convert to raiders in this case */
							double x1 = pf1(rng);
							if(x1 < praiders) {
								c1.st = cell::state::RAIDERS;
								raiders_created++;
							}
						}
					}
				}
				else {
					conflict = false;
					attack_fail++;
					if(defenders_create_raiders) {
						x1 = pf1(rng);
						if(x1 < praiders) {
							c2.st = cell::state::RAIDERS;
							raiders_created++;
						}
					}
				}
				
				if(conflict_effect_base > 0.0 && (success || !conflict_effect_only_on_success)) {
					/* effect of conflict on the target and nearby cells */
					double max_dist = conflict_dist * 10.0;
					
					for(const auto& tmp1 : g.neighbors(dst, max_dist, true)) {
						const auto& pt2 = tmp1.first;
						if(raiders_tolerate_self && src == pt2) continue;
						if(refugees_pts.count(pt2)) continue;
						
						auto& c3 = g.at(pt2);
						bool rtmp = false;
						double p1 = conflict_effect_base * exp(-1.0 * tmp1.second / conflict_dist);
						if(only_instant_conflict) {
							if(c3.st != cell::state::EMPTY && c3.defensibility < 1.0 &&
									(c3.st == cell::state::FARMING || !raiders_tolerate_conflict)) {
								if(pf1(rng) < p1) rtmp = true;
							}
						}
						else {
							if(!discrete_conflict) c3.conflict_level = std::min(1.0, c3.conflict_level + p1);
							else if(pf1(rng) < p1) c3.conflict_level = 1.0;
							if(c3.conflict_level > c3.defensibility && c3.st != cell::state::EMPTY &&
								(c3.st == cell::state::FARMING || !raiders_tolerate_conflict)) rtmp = true;
						}
						
						if(rtmp) {
							refugee_order tmp2;
							tmp2.src = pt2;
							tmp2.r = 1.0;
							tmp2.seq = osd(rng);
							refugees.push_back(tmp2);
							refugees_pts.insert(pt2);
						}
					}
				}
			}
			
			if(migration) {
				/* in this case, c1.st == cell::state::FARMING always
				 * (we only have stationary raiders, migration happens
				 * for farming source cells only) */
				if(c1.st != cell::state::FARMING) throw std::runtime_error("Inconsistent migration!\n");
				/* this is the case of mobile split-off from a farming cell */
				c1.N -= N3;
				c2.N += N3;
				c2.st = cell::state::FARMING;
				if(conflict) {
					double x1 = pf1(rng);
					if(x1 < praiders) {
						c2.st = cell::state::RAIDERS;
						raiders_created++;
					}
				}
			}
			
			/* in (the very unlikely) case that everyone migrated away, mark the source cell as empty */
			if(c1.N == 0) c1.st = cell::state::EMPTY;
			
			if(c2.st != dsts && dsts == cell::state::EMPTY && !this->external_helper) {
				/* cell became occupied might need to create the helper distribution for it */
				if(cph.count(dst) == 0) this->create_cell_prob_helper_one(dst);
			}
			
		}
		
		/* process refugees (in random order) */
		refugee_order::sort(refugees);
		for(const auto& o : refugees) {
			//!! TODO: should we use the "regular" migration distances here?
			double max_dist = refugee_dist * 10.0;
			auto& c = g.at(o.src);
			std::vector<double> w;
			std::vector<typename cell_collection::point> targets;
			double sum1 = 0.0;
			for(const auto& tmp1 : g.neighbors(o.src, max_dist, true)) {
				const auto& pt2 = tmp1.first;
				const auto& c2 = g.at(pt2);
				
				if(c2.st == cell::state::EMPTY) continue;
				if(refugee_target_def) { if(c2.defensibility < 1.0) continue; }
				else if(c2.conflict_level > c2.defensibility) continue;
				
				double p1 = exp(-1.0 * tmp1.second / refugee_dist) * c2.K / (double) c2.N;
				w.push_back(p1);
				targets.push_back(pt2);
				sum1 += p1;
			}
			
			if(targets.size() == 0) {
				cancelled_refugee_cells++;
				cancelled_refugees += c.N;
				continue; //!! TODO: should this be an error? there should always be a refuge close by
			}
			
			
			double Ntot = c.N;
			double Norig = Ntot * o.r;
			for(size_t i = 0; i < targets.size(); i++) {
				double N1 = w[i] * Norig / sum1;
				Ntot -= N1;
				unsigned int N2 = (unsigned int)std::round(Ntot);
				if(N2 < c.N) {
					unsigned int tmp = c.N - N2;
					g.at(targets[i]).N += tmp;
					c.N = N2;
					total_refugees += tmp;
				}
			}
			
			// we should have c.N == 0 at this point
			if(o.r == 1.0 && c.N != 0) throw std::runtime_error("Error distributiong refugees!\n");
			if(c.N == 0) c.st = cell::state::EMPTY;
		}
		
		refugees.clear();
	}
	
	template<class point>
	struct noop_cb {
		void operator ()(const point& pt1, const point& pt2, unsigned int N) const { }
	};
	
	void apply_group_mig() { apply_group_mig(noop_cb<typename cell_collection::point>()); }
	
	void do_group_mig() {
		create_group_orders();
		apply_group_mig(noop_cb<typename cell_collection::point>());
	}
	
	template<class CB>
	void do_group_mig(CB&& cb) {
		create_group_orders();
		apply_group_mig(cb);
	}
	
	/* step the demographic simulation in each cell */
	void cells_step(double t) {
		total_dec1 = 0;
		total_dec2 = 0;
		this->ensure_omp_threads();
		size_t tmpsize = g.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for(size_t i = 0; i < tmpsize; i++) {
			auto& c = g.at_ix(i); /* cell to update */
			/* decrease conflict levels */
			if(c.conflict_level > 0.0) c.conflict_level = std::max(0.0, c.conflict_level - raiders_revert_rate);
			if(c.st == cell::state::EMPTY) continue;
			unsigned int tmp = c.N;
			c.step(t, this->prng(), this->pars);
			if(tmp > c.N) {
				if(tmp > c.K) {
#ifdef _OPENMP		
#pragma omp atomic
#endif
					total_dec2 += (tmp - c.N);
				}
				else {
#ifdef _OPENMP		
#pragma omp atomic
#endif
					total_dec1 += (tmp - c.N);
				}
			}
			/* in the case this cell became empty, we need to update the state */
			if(c.N == 0) c.st = cell::state::EMPTY;
		}
	}
};

#endif

