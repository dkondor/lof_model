/*
 * neolithic_common.hpp -- common helper to run the simulation with or
 * 	without raiders
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
 */


#ifndef NEOLITHIC_COMMON_HPP
#define NEOLITHIC_COMMON_HPP

#include "neolithic_base.hpp"
#include "neolithic_w.hpp"
#include "cellnet.hpp"
#include "cellnet_options.hpp"
#include <type_traits>
#include <stdio.h>


/* output to files -- need to be specialized to be used! */
struct output_files {
	FILE* summary_out = stdout;
	FILE* detail_out = nullptr;
	unsigned int out_period = 10;
	bool track_avg_distance = false;
	bool write_progress = true;
	unsigned int steps = 0; /* total number of steps (only used for writing the progress) */
	
	using nclass = neolithic_w<cellnet>;
	
	void operator()(const nclass& nn, unsigned int i, const step_res<nclass>& sr) {
		FILE* out1 = summary_out;
		FILE* outf = detail_out;
		
		if(out1) {
			fprintf(out1, "%u\t%u\t%" PRIu64 "\t%u\t%" PRIu64 "\t%u\t%" PRIu64 "\t%u\t%" PRIu64 "\t%u\t%u", i,
				sr.t.farming_cells, sr.t.farmers, sr.t.raider_cells, sr.t.raiders,
				sr.t.farming_cells_r, sr.t.farmers_r, sr.t.raider_cells_r, sr.t.raiders_r,
				 nn.raiders_created, nn.raiders_reverted);
			if(track_avg_distance) {
				fprintf(out1, "\t%f\t%u", sr.sum_dist / (double)sr.cnt_dist, sr.total_mig);
				fprintf(out1, "\t%f", sr.sum_dist_new / (double)sr.cnt_dist_new);
			}
			fprintf(out1, "\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n", sr.t.conflict_cells, sr.t.excess_pop,
				nn.total_killed, nn.total_dec1, nn.total_dec2, nn.total_refugees,
				nn.cancelled_refugee_cells, nn.cancelled_refugees);
		}
		
		if(outf && (i % out_period == 0)) nn.g.write_pop(outf, i);
		
		if(write_progress) {
			fprintf(stderr, "%u / %u steps complete\r", i, steps);
			fflush(stderr);
		}
	}
};


class neolithic_sim {
	protected:
		neolithic_w<cellnet> nn;
		cellnet_options opts;
	
	public:
		int run(int argc, char** argv) {
			option_parser op;
			nn.add_options(op);
			opts.add_options(op);
			
			if(!op.parse_options(argc, argv)) return 1;
			
			if(opts.weather_allow_missing && opts.weather_cut) {
				fprintf(stderr, "Incompatible options!\n");
				return 1;
			}
			
			/* read the simulation space */
			if(opts.helper_matrix_fn) nn.read_combined_matrix(opts.helper_matrix_fn, opts.helper_only_net);
			else {
				opts.read_net_data(nn);
				if(opts.use_distance_matrix) nn.g.create_distance_matrix();
			}
			opts.do_adjust_K(nn, (opts.helper_matrix_fn != nullptr));
			opts.ensure_start_cell(nn);
			
			/* add the starting population */
			nn.set_start_pop(opts.start_cell, opts.starting_pop, false, opts.start_pop_fill_factor, opts.start_pop_pchoice);

			/* set defensibilities */
			nn.set_cell_defensibility();
			
			/* open the detailed output file (if needed) */
			FILE* outf = nullptr;
			if(opts.out_base) {
				/* optionally allow writing detailed results to the standard output */
				if(opts.out_base[0] == '-' && opts.out_base[1] == 0) outf = stdout;
				else {
					if(opts.out_base_zip) {
						char* tmp1 = new char[strlen(opts.out_base) + 30];
						sprintf(tmp1, "/bin/gzip -c -1 > %s", opts.out_base);
						outf = popen(tmp1, "w");
						delete[]tmp1;
					}
					else outf = fopen(opts.out_base, "w");
					if(!outf) {
						fprintf(stderr, "Error opening output file %s!\n", opts.out_base);
						return 1;
					}
				}
			}
			
			
			output_files of;
			of.detail_out = outf;
			/* note: if we are writing the detailed results to stdout,
			 * we need to supress the summary output that would normally
			 * go there */
			if(outf == stdout) of.summary_out = nullptr;
			of.steps = opts.steps;
			of.write_progress = true;
			of.out_period = opts.out_period;
			of.track_avg_distance = opts.track_avg_distance;
			
			if(opts.weather_file) {
				weather_update_file wf;
				opts.read_weather_mapping(wf, nn);
				wf.open_file(opts.weather_file);
				do_one_run(nn, opts.steps, opts.start_cell, wf, of, true);
			}
			else do_one_run(nn, opts.steps, opts.start_cell, weather_update_noop(), of, true);
		
			if(outf && outf != stdout) {
				if(opts.out_base_zip) pclose(outf);
				else fclose(outf);
			}
			return 0;
		}
};

#endif

