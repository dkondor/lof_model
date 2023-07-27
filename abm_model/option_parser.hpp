/*
 * option_parser.hpp -- helper class to organize a set of command line
 * 	options that might belong to different subcomponents
 * ensures that there are no conflicting option flags (only at runtime)
 * and that each command line option was processed
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

#include <unordered_map>
#include <string>
#include <memory>
#include <utility>
#include <functional>
#include <stdexcept>
#include <stdio.h>
	

class option_parser {
	public:
		/* actions to be taken have the following signature:
			int action(int i, int argc, char** argv);
		   where i is the current position in argv (of the matched option)
		   and argc and argv correspond to the command line arguments
		   the function should return the additional arguments consumed
		   (not counting the matched parameter, i.e. flags should return 0)
		   or -1 on error
		 */
		using action_t = std::function<int(int, int, char**)>;
		/* action for unknown parameters -- can be set externally
		 * if not set, the default action is to write an error
		 * message to stderr */
		action_t default_action;
		
	protected:
		std::unordered_map<std::string, action_t> options;
		std::vector<std::function<bool()> > post_cbs;
		
	public:
		/* add the given action to the given option (without the leading '-') */
		template<class F>
		void add_option(std::string opt, F&& action) {
			auto res = options.emplace(opt, std::move(action));
			if(!res.second) throw std::runtime_error("option_parser::add_option(): duplicate option!\n");
		}
		
		template<class num>
		static void parse_num(num* res, const char* str);
		
		/* convenience wrapper for adding a simple numeric option
		 * (with a bool indicator in flag set if it was found) */
		template<class num>
		void add_numeric_option(std::string opt, num* target, bool* flag = nullptr) {
			add_option(opt, [target, flag](int i, int argc, char** argv) {
				if(i + 1 < argc) {
					parse_num(target, argv[i+1]);
					if(flag) *flag = true;
					return 1;
				}
				else return -1;
			});
		}
		
		/* convenience wrapper for adding a simple numeric option
		 * (with a bool indicator in flag set if it was found) */
		void add_string_option(std::string opt, char** target, bool* flag = nullptr) {
			add_option(opt, [target, flag](int i, int argc, char** argv) {
				if(i + 1 < argc) {
					*target = argv[i+1];
					if(flag) *flag = true;
					return 1;
				}
				else return -1;
			});
		}
		
		/* convenience wrapper for adding a simple bool option set to
		 * the given value if the target is found */
		void add_bool_option(std::string opt, bool* target, bool value = true) {
			add_option(opt, [target, value](int, int, char**) { *target = value; return 0; });
		}
		
		/* Add a callback that will be executed after all command line
		 * options have been succesfully processed; these can be used
		 * e.g. for consistence checks. Callbacks should return a bool
		 * which will be propagated to the caller. If any callback
		 * returns false, the rest is not called.
		 * Note: if there is an error processing any option, callbacks
		 * are not called. */
		template<class F>
		void add_post_callback(F&& cb) { post_cbs.emplace_back(std::move(cb)); }
		
		/* parse all options on the command line, starting from argv[1]
		 * returns true on success, false on any error encountered */
		bool parse_options(int argc, char** argv) {
			for(int i = 1; i < argc; i++) {
				action_t* action = nullptr;
				
				if(argv[i][0] == '-') {
					std::string tmp(argv[i] + 1);
					auto it = options.find(tmp);
					if(it != options.end()) action = &(it->second);
				}
				
				action_t& action2 = action ? *action : default_action;
				if(action2) {
					int tmp = action2(i, argc, argv);
					if(tmp >= 0) i += tmp;
					else return false;
				}
				else {
					fprintf(stderr, "Unknown parameter: %s!\n", argv[i]);
					return false;
				}
			}
			
			for(auto& x : post_cbs) if(!x()) return false;
			return true;
		}
};


template<> void option_parser::parse_num<int>(int* res, const char* str) { *res = atoi(str); }
template<> void option_parser::parse_num<unsigned int>(unsigned int* res, const char* str) { *res = atoi(str); }
template<> void option_parser::parse_num<uint64_t>(uint64_t* res, const char* str) { *res = strtoul(str, 0, 10); }
template<> void option_parser::parse_num<double>(double* res, const char* str) { *res = atof(str); }


