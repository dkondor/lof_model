/*
 * cellnet_options.hpp -- options only relevant when running the
 * 	simulation with a network of cells
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

#ifndef CELLNET_OPTIONS_HPP
#define CELLNET_OPTIONS_HPP
#include "neolithic_base.hpp"
#include "cellnet.hpp"

struct cellnet_options : public common_options<cellnet> {	
	char* net_file = nullptr; /* read a network from this file */
	bool net_file_edges = false; /* set to true if the file with the network is an edgelist */
	bool net_file_bin = false; /* set to true if the cellnet object should be read from a binary file (without the probability helpers) */
	bool net_file_write = false; /* set to true if the cellnet object should be written to a binary file */
	bool net_recalculate_distances = false; /* set to true to recalculate distances along edges after reading the network */
	
	std::vector<const char*> extra_edges_files; /* files containing additional edges */
	
	char* coords_file = nullptr; /* read cell (center) coordinates from this file */
	bool coords_file_csv = false;
	bool coords_file_has_types = false; /* set to true if the above file includes information about coastal cells (as cell type) */
	bool coords_file_has_factors = false; /* set to true if the above file includes land travel scaling factors */
	double coastal_travel_factor = 1.0; /* factor used for coastal travel */
	
	bool have_K = false; /* set to true if we should load carrying capacities from the above input file */
	std::vector<const char*> scale_K_fn; /* scale productivity values by these given as inputs */
	bool scale_K_geo = false; /* scale productivity values by a factor corresponding to the latitude */
	
	/* use a distance matrix */
	bool use_distance_matrix = false;
	
	/* add specialized options to parser; note: this is NOT a virtual function,
	 * so this needs to be called on an instance of this object and not
	 * on the base class */
	void add_options(option_parser& op) {
		common_options<cellnet>::add_options(op);
		op.add_bool_option("Kg", &scale_K_geo);
		op.add_option("Ks", [this] (int i, int argc, char** argv) {
			if(i + 1 >= argc) return -1;
			scale_K_fn.push_back(argv[i+1]);
			return 1;
		});
		op.add_bool_option("k", &have_K);
		op.add_bool_option("ct", &coords_file_has_types);
		op.add_bool_option("cl", &coords_file_has_factors);
		op.add_string_option("c", &coords_file);
		op.add_string_option("cc", &coords_file, &coords_file_csv);
		op.add_bool_option("nr", &net_recalculate_distances);
		op.add_numeric_option("nC", &coastal_travel_factor);
		op.add_option("nE", [this] (int i, int argc, char** argv) {
			if(i + 1 >= argc) return -1;
			extra_edges_files.push_back(argv[i+1]);
			return 1;
		});
		op.add_bool_option("nw", &net_file_write);
		op.add_string_option("n", &net_file);
		op.add_string_option("ne", &net_file, &net_file_edges);
		op.add_string_option("nb", &net_file, &net_file_bin);
		op.add_bool_option("z", &use_distance_matrix);
		
		op.add_post_callback([&]() { if(net_file_bin) net_file_edges = true; return true; });
	}
	
	/* load a network data (only if cell_collection supports it) */
	void read_net_data(nbase<cellnet>& nn) {
		if(helper_matrix_fn) throw std::runtime_error("common_options::read_net_data() should only be called if no external helper was given!\n");
		if(!(net_file && coords_file)) throw std::runtime_error("common_options::read_net_data(): missing input filenames!\n");
		
		uint32_t flags = 0;
		if(have_K) flags |= cellnet::load_net_flags::read_K;
		if(net_file_edges) flags |= cellnet::load_net_flags::read_edges;
		if(coords_file_has_types) flags |= cellnet::load_net_flags::read_types;
		if(coords_file_has_factors) flags |= cellnet::load_net_flags::read_travel_factors;
		
		nn.g.load_net(net_file, coords_file, flags, coords_file_csv);
		for(auto fn : extra_edges_files) nn.g.read_edges(fn);
		if(net_recalculate_distances) nn.g.calculate_edge_distances(coastal_travel_factor);
		
		if(have_K) {
			if(scale_K_geo) nn.g.scale_K_lat();
			for(const char* fn1 : scale_K_fn) nn.g.scale_K(fn1, true);
		}
		else for(const auto& pt : nn.g) nn.g.at(pt).K = K1;
	}
};

#endif

