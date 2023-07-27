/*
 * runmultiple.cpp -- Execute a set of shell commands in parallel
 * 
 * Copyright 2023 Daniel Kondor <kondor@csh.ac.at>
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 * * Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above
 *   copyright notice, this list of conditions and the following disclaimer
 *   in the documentation and/or other materials provided with the
 *   distribution.
 * * Neither the name of the  nor the names of its
 *   contributors may be used to endorse or promote products derived from
 *   this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */


#include <iostream>
#include <fstream>
#include <istream>
#include <thread>
#include <atomic>
#include <string>
#include <vector>
#include <stdlib.h>

int main(int argc, char **argv)
{
	char* fn = 0;
	int nthreads = -1;
	for(int i = 1; i < argc; i++) {
		if(argv[i][0] == '-') switch(argv[i][1]) {
			case 'f':
				fn = argv[i+1];
				i++;
				break;
			case 't':
				nthreads = atoi(argv[i+1]);
				if(!nthreads) {
					fprintf(stderr, "Invalid value for the number of threads!\n");
					return 1;
				}
				i++;
				break;
			default:
				fprintf(stderr, "Unknown parameter: %s!\n", argv[i]);
				break;
		}
		else fprintf(stderr, "Unknown parameter: %s!\n", argv[i]);
	}
	
	if(nthreads < 0) {
		nthreads = std::thread::hardware_concurrency();
		if(nthreads <= 0) nthreads = 1;
	}
	
	std::ifstream fs;
	if(fn) fs.open(fn);
	
	std::istream& is = fn ? fs : std::cin;
	
	std::vector<std::string> cmds;
	for(std::string l; std::getline(is, l);) cmds.push_back(std::move(l));
	
	if(nthreads == 1) {
		for(const auto& cmd : cmds)
			if(system(cmd.c_str()) < 0) throw std::runtime_error("Error running the child process!\n");
	}
	else {
		std::atomic<size_t> j{0};
		std::vector<std::thread> t;
		t.reserve(nthreads);
		auto fn = [&cmds, &j]() {
			while(true) {
				size_t x = j++;
				if(x >= cmds.size()) break;
				if(system(cmds[x].c_str()) < 0) throw std::runtime_error("Error running the child process!\n");
			}
		};
		for(int i = 0; i < nthreads; i++) t.emplace_back(fn);
		for(int i = 0; i < nthreads; i++) t[i].join();
	}
	
	return 0;
}

