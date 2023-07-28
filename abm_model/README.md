# LoF agent-based model

Agent-based simulation of density-dependent conflict with landscape of fear effects.

## Examples for the simulation

Examples of a few cases are given below that were used to generate Fig. 2 in our paper. More simulation runs are in the script [lof_runs_rep.fsh](lof_runs_rep.fsh) that performs repeated realizations of the simulation and analysis of temporal patterns in the results, generating the main results in our paper. These results are then analyzed with the script [lof_acf_rep_plot.r](lof_acf_rep_plot.r), creating the figures in the paper.

The below commands are provided to be run in a standard UNIX (Linux, OSX) shell (such as `bash`). They require a `C++17` compiler such as a recent version of GCC or CLang (to use CLang, replace the `g++` commands with `clang++`). The code for creating images and videos (`spout2png.cpp`) uses the [Cairo](https://www.cairographics.org/) graphics libraries. [FFMPEG](https://ffmpeg.org/) is used for creating videos. Running the below code will require at least 8GB free memory.

The simulation relies on the same base data and perprocessing steps as our previous model available [here](https://github.com/dkondor/neolithic_simulation/) (but without the climate data).

All code is assumed to be run from the current directory, and assumes a shell variable `$base_dir` set to the base directory where all data files are stored after the preprocessing steps. Files needed for the simulation and results are also written here, under the subdirectory `simulation`; if this does not exist, it should be created before running any of the below commands:
```
mkdir -p $base_dir/simulation
```


### 1. Compile the simulation code

Manual compilation from the command line with GCC:
```
g++ -o cm create_helper_matrix.cpp read_table_cpp.cpp -std=gnu++17 -O3 -march=native -lm -fopenmp -pthread
g++ -o n2w neolithic2w.cpp read_table_cpp.cpp -std=gnu++17 -O3 -march=native -lm -fopenmp
g++ -o spas spaggr_simple.cpp read_table_cpp.cpp -std=gnu++17 -O3 -march=native
g++ -o rm runmultiple.cpp -std=gnu++17 -pthread
g++ -o sp2png spout2png.cpp read_table_cpp.cpp cellnet_img.cpp -O3 -march=native -std=gnu++17 -I/usr/include/cairo -I/usr/include/glib-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include -I/usr/include/pixman-1 -I/usr/include/uuid -I/usr/include/freetype2 -I/usr/include/libpng16 -lcairo
```

Notes:
 - The `-fopenmp` switch is optional and will enable parallelization in some parts of the code; it requires support for OpenMP.
 - The `rm` program is only needed for running multiple simulations in parallel by the analysis [script](lof_runs_rep.fsh) and can be skipped if that script is not used.
 - Compiling the `sp2png` program (the last command) above requires options locating the Cairo graphics libraries. The development headers for this needs to be installed, e.g. on Ubuntu, the `libcairo2-dev` package. The command line above includes default locations on typical Linux distributions; this likely requires adjustments on other systems. On many systems, the `pkg-config --libs --cflags cairo` command can be used in determining the options needed. Note that this program is only needed for visualization of the results.


Alternatively, the [Meson](https://mesonbuild.com/) build system can be used to compile the codes:
```
meson -Dbuildtype=release build
ninja -C build
```
This will place the executables under the newly created `build` subdirectory. Note that the following examples assumes the output of manual compilation, i.e. that all executables are under the current directory. Copy them or create symbolic links as follows:
```
ln -s build/cm
ln -s build/n2w
ln -s build/spas
ln -s build/sp2png
```


### 2. Create the base probability distributions for the migrations

These are used by the simulation later. The main parameter is given with the `-G` switch that determines the characteristic distance in kilometers (the parameter of the exponential function that gives the decrease of migration probability with distance). In the current exapmles and all analysis in our paper G = 80 km is used (note: the `-t` parameter controls how many CPU cores to use in parallel, here it is set to 4).


```
./cm -o $base_dir/simulation/matrix_dg_G80zC10.bin -G 80 -k -K 150 -Kt 50 -cc $base_dir/gaez_crop/eu_dggrid_coords_K_land4a.csv -ct -n $base_dir/dggrid/isea3h12eun.nbr -nr -nC 0.1 -nE $base_dir/eu_dggrid_extra_edges.dat -Ks $base_dir/gaez_crop/eu_dggrid_land_share.csv -z -t 4
```

These result in output files with an approximate size of 6.2 GB (assuming a study area with 36400 cells)


### 3. Run the simulation

In the following, we use the `$outdir` variable as the target of simulation results. E.g.
```
outdir=$base_dir/simulation/lof_ex
mkdir -p $outdir
```

Two examples are given here, for a set of default parameters and an 1% or 10% share of refuge cells. These create the results shown in Fig. 2 of our paper.

10% share of refuges (R = 0.1):
```
$program_dir/n2w -Hf $base_dir/simulation/matrix_dg_G80zC10.bin -LI -LS -Ld 0 -LD -Lr 0.1 -Lb 1 -Lt 0 -Rp 1 -RP -RA 0.1 -R 8 -RR 0.01 -Rc -Rp 1 -Rs 1 -RS -d 1 -s 899705504 -S 10000 -E 5 -a 0.5 -o $outdir/res_G80_C10_E5_A10_S1_L1_R0.1_2sp.out.gz -oz -op 1 -i 1596250 -I 5000 -If 0.5 -Ip 0.5 -Ce 0 -Cb 0 -De 1 -Dm 0.5 -Db 0.25 -Dl 200 > $outdir/res_G80_C10_E5_A10_S1_L1_R0.1_2.out
```

1% share of refuges (R = 0.01):
```
$program_dir/n2w -Hf $base_dir/simulation/matrix_dg_G80zC10.bin -LI -LS -Ld 0 -LD -Lr 0.01 -Lb 1 -Lt 0 -Rp 1 -RP -RA 0.1 -R 8 -RR 0.01 -Rc -Rp 1 -Rs 1 -RS -d 1 -s 899705504 -S 10000 -E 5 -a 0.5 -o $outdir/res_G80_C10_E5_A10_S1_L1_R0.01_2sp.out.gz -oz -op 1 -i 1596250 -I 5000 -If 0.5 -Ip 0.5 -Ce 0 -Cb 0 -De 1 -Dm 0.5 -Db 0.25 -Dl 200 > $outdir/res_G80_C10_E5_A10_S1_L1_R0.01_2.out
```


### 4. Spatial aggregation

This step aggregates population numbers in the tiling that was created for spatial aggregation in the preprocessing step (examples are given for the above two cases of the simulation):
```
zcat $outdir/res_G80_C10_E5_A10_S1_L1_R0.1_2sp.out.gz | ./spas -m $base_dir/new_grid_aggr_flat_r500.csv -f $base_dir/new_dgid_aggr_filter.dat | gzip -c > $outdir/res_G80_C10_E5_A10_S1_L1_R0.1_2spaggr.out.gz
zcat $outdir/res_G80_C10_E5_A10_S1_L1_R0.01_2sp.out.gz | ./spas -m $base_dir/new_grid_aggr_flat_r500.csv -f $base_dir/new_dgid_aggr_filter.dat | gzip -c > $outdir/res_G80_C10_E5_A10_S1_L1_R0.01_2spaggr.out.gz
```


### 5. Compile videos

This uses the `sp2png` program that has two modes of operation. The first one outputs a series of PNG images that can be later used by any video editing software to create video files:
```
mkdir $outdir/video_example
zcat $outdir/res_G80_C10_E5_A10_S1_L1_R0.1_2sp.out.gz | ./sp2png -oy -7000 -op 10 -om 200 -ow 1608 -oh 900 -o $outdir/video_example/res_G80_C10_E5_A10_S1_L1_R0.1 -n ../../../DGGRID/test1/isea3h12eun.nbr -cp $base_dir/dggrid/dggrid_poly.csv
```
Note that the input is always read from the standard input; this should be redirected from the previously created output file. This will create a series of PNG images under `$outdir/video_example` that can be used e.g. as an input to `ffmpeg` or any similar software.

An alternative (and preferred) way is to run ffmpeg directly while generating the images. This can be achieved with the `-of` option that will try to run ffmpeg with the correct options and the given file as the final output, supplying a set of images to encode as a video directly:

```
zcat $outdir/res_G80_C10_E5_A10_S1_L1_R0.1_2sp.out.gz | ./sp2png -oy 0 -op 10 -om 200 -ow 1608 -oh 900 -of $outdir/res_G80_C10_E5_A10_S1_L1_R0.1_2.mp4 -oc "h264 -pix_fmt yuv420p" -n ../../../DGGRID/test1/isea3h12eun.nbr -cp $base_dir/dggrid/dggrid_poly.csv
zcat $outdir/res_G80_C10_E5_A10_S1_L1_R0.01_2sp.out.gz | ./sp2png -oy 0 -op 10 -om 200 -ow 1608 -oh 900 -of $outdir/res_G80_C10_E5_A10_S1_L1_R0.01_2.mp4 -oc "h264 -pix_fmt yuv420p" -n ../../../DGGRID/test1/isea3h12eun.nbr -cp $base_dir/dggrid/dggrid_poly.csv
```
Note that the output container format will depend on the filename given (multiple formats are supported by ffmpeg; the above example uses the MP4 container format). The actual video codec and additional options can be selected with the `-oc` parameter; in the above example, H264 is used which works reasonably fast, and also the YUV420P image format is selected which makes the result compatible with more video players.



## Program parameters

A more detailed description of the command line arguments of the individual programs.

## cm / create_helper_matrix.cpp

This program creates a binary file that contains the cumulative distribution function of migration probabilities to be used during the simulations. Essentially, the distributions form an NxN matrix (for N cells); given that any one probability can be calculated on the fly, it is sufficient to store NxN/2 values. Optionally, a matrix of distances can be stored as well; this is useful if non-trivial adjustments such as scaling factors or long-range connections are used.

Command line parameters used:

 - `-o fn` output file name
 - `-K num` scale carrying capacities so that the average capacity of a cell is the value given here (recommended to use, since the input files contain arbitrary values)
 - `-Kg` additionally scale the capacity in each cell by the cosine of the latitude (to account for cells whose size depends on the latitude -- this is not used with DGGRID, where size of cells is constant)
 - `-Ks fn` read additional scaling factors for each cell from the given file (this is used for taking into account the share of actual usable soil in each cell); this option can be used multiple times, resulting in the product of all scaling factors used
 - `-c fn` or `-cc fn` name of the file with a list of cell IDs, coordinates and (optionally) base carrying capacities; by default, this file is expected in a TSV format without header; use the `-cc` form if the input is a CSV with header
 - `-k` has to be given if the above file contains carrying capacities (e.g. GAEZ base estimates); otherwise a constant value is assumed (note: this is used in all cases)
 - `-n fn` or `-ne fn` name of the file containing the neighbor relations among cells; use `-n` if this file is in the DGGRID output format (all neighbors of a cell are listed on a line), and the `-ne` variant if this file is an edgelist (each line contains two node IDs); the input file should be symmetric, i.e. it should contain each edge both ways
 - `-nC num` scaling factor for travel distance between coastal cells; calculated distance is multiplied by this, so giving a number < 1 allows faster / further travel along coastal cells
 - `-nE fn` name of a file containing additional edges (as an edgelist) that should be added to the cell neighbor relations; this is currently used to add links between non-adjacent cells
 - `-nr` if given, distances are calculated along neighbor edges taking into account any scaling factor 
 - `-z` calculate distances along shortest paths in the cell network, taking into account any scaling, and store a matrix of such distances in the output
 - `-C` instead of creating a binary file, read one and check that the probability distributions are correct in it
 - `-G dist` characteristic distance used in the migration probabilities (exponential decay)
 - `-t num` number of threads to use for distance and probability calculations

## Main simulation: n2w (neolithic2w.cpp)


### Initial conditions and basic settings:
 - `-i ID` ID of the cell where to place a starting population (a location in Anatolia is used in the examples: cells 1596250 and 1607919 for the cases with and without conflict respectively)
 - `-I pop` size of the starting population to use (default: 100); if this is larger than the carrying capacity of the starting cell, the additional population is distributed in neighboring cells
 - `-If factor` when distributing the initial population, cells are only filled up to this factor (e.g. 0.5); this option is useful to avoid an initial "explosion" of migrations that would happen if the initial conditions include cells close to their carrying capacity
 - `-Ip factor` when distributing the initial population, any neighbor cell is chosen with this probability (e.g. 0.5); this can be used to start from a less tightly packed configuration of occupied cells
 - `-S n` number of steps (years) to run the simulation for
 - `-K num` scale carrying capacities so that the average over all cells is this value (can be adjusted regardless of the value given when creating the migration distributions, since this means scaling the whole distribution with the same number)
 - `-T` if given, the average distance of newly settled cells (from the starting cell) is tracked and output
 - `-o fn` detailed output filename; the population in each cell is written here every 10 years
 - `-oz` if present, the above output file is compressed with `/bin/gzip`
 - `-s num` number to use to seed the random number generator (current time is used by default)
 - `-t num` number of parallel OpenMP threads to use for evaluating population growth and split-offs (default: 1, i.e. no parallelism as it matters very little in this case)
 - `-r r` base population growth rate (per year)
 - `-d r` population collapse factor: if population is above carrying capacity, it collapses to this ratio of it (values below one can model a larger collapse); a negative number means use the logistic equation even if population is above carrying capacity

### Group fission dynamics:
 - `-Hf fn` name of the binary file containing the pre-computed probability distributions (created with `cm`)
 - `-Dl num` Dunbar number to use (default: 150)
 - `-Db p` probability of a group splitting off if the population of a cell is equal to the Dunbar number
 - `-Dm r` minimum ratio of population to the Dunbar number to consider a possibility of split-off
 - `-De e` exponent used to calculate the probability of split-off due to approaching the Dunbar limit
 - `-Cb p` probability of a group splitting off if the population of a cell is equal to the currnt local carrying capacity (note: this is set to zero in all simulations)
 - `-Cm r` minimum ratio of population to the carrying capacity to consider a possibility of split-off
 - `-Ce e` exponent used to calculate the probability of split-off due to approaching the carrying capacity limit

Note: each year, the possiblity of a group splitting off (and migrating elsewhere) is considered for each cell. This is done by first estimating a probability for this to happen based on "scalar stress" considering both the absolute group size and relative group size compared to the current carrying capacity. In each case, this is done in the following way:
```
P = p * ((x - r) / (1 - r))**e
```
where `x` is the current population relative to the limit considered (either the "Dunbar number", given by the `-Dl` parameter, in case of considering absolute group size; or the local carrying capacity), and `p`, `r` and `e` are the parameters given above. By default, both mechanisms are considered, and a migration occurs if either of them yields it by random choice (given the probability computed as above). Setting the `-Db` or `-Cb` value to zero disables that mechanism.

### Climate data (for variation in agriculture):

Note: this functionality was not used in our analysis, but is carried over from the [original model implementation](https://github.com/dkondor/neolithic_simulation).

 - `-w fn1 fn2` specifies the names of two files with the data to use; `fn1` should be the mapping between weather cells and simulation cells (i.e. `cell_ids_dggrid.csv` in this case), `fn2` should be the file with the actual variations (`crop_data_cmb.dat` in this case)
 - `-wm` if given, missing cells are allowed; without this, every cell has to appear in all years in the input file
 - `-wc` if given, cut the simulation area to cells that have climatic variation data
 - `-wf` if given, the mapping between weather cells and simulation cells contains weighting factors; this means that multiple weather cells overlap with one simulation cell and a weighted average is used to calculate any actual effect of climate variation; this option is necessary when using the hexagon grid (maybe this should be made the default)
 - `-wC` if given the mapping file is in CSV format (this is necessary for the above files)
 - `-wF year` if given, assume that the start year of the simulation is the one given here (as an absolute number, in calendar years, i.e. BCE / CE); data before this date in the input file is skipped; if not given, the first year is the input file is assumed to be the beginning of the simulation
 - `-ws factor` scale any variation in yield by this factor; this can be used to exaggerate the effect of climatic variability
 

### Command line parameters specific to agressors:

 - `-R dist` characteristic distance in the probability of choosing a target for aggressors (this can be varied independently from the migration probability distribution parameter)
 - `-Rc` if a cell is succesfully defended from an attack, it can alos turn into aggressors (with a probability given by the `-Rp` parameter)
 - `-Rp prob` probability that a winner of a conflict converts to aggressors
 - `-RP` if given, aggressors take into account the current farmer population when choosing a target (they prefer cells with larger populations)
 - `-RR rate` yearly rate at which aggressors revert back to farmers
 - `-RA rate` yearly rate for aggressors to engage in conflicts
 - `-Rs prob` share of survivors after a successful attack (i.e. 1 means that everyone survives, 0 means everyone is killed)
 - `-Rm dist` maximum distance that aggressors are able to attack (default is 10x the distance given by the `-R` option)
 - `-m` aggressors are considered mobile: after each successfull attack, they take over the attacked cell
 - `-a prob` probability of an attack (by aggressors) being successful (e.g. setting this to 1 means that they always win; this is used in the default version of the simulation)
 - `-Ra dist` if given, the probability of successful attack (by aggressors) depends on the distance exponentially, using the characteristic distance given here (i.e. as `a * exp(- d / dist)`, where `a` is the base success rate given by the `-a` parameter, `d` is the distance of raiders to the target cell, and `dist` is the paremeter value given here; not used in the main analysis)
 - `-E num` preference to choose an empty cell for group migrations (typical value: 10, default value is 1); i.e. a larger value will result in less aggressors to be created until the landscape is more saturated; note: this is the inverse of the p_E parameter in the paper

### Command line parameters specific to the landscape of fear effect:

 - `-Lb num` base effect of conflicts; i.e. this gives the probability of inhabitants of the attacked cells becoming refugees; probabilities for neighbor cells decrease with distance
 - `-Lc dist` characteristic distance of conflict effects (i.e. decline of effect by an exponential); if not given, aggressors' characteristic distance is used instead (`-R`)
 - `-Lr prob` probability of any cell to be a refuge
 - `-LD` "discrete" version of model: a cell is considered either safe / unsafe
 - `-LI` conflicts only have an immediate effect (refugee probabilities are evaluated only immediately after the conflict, but cells are not marked unsafe)
 - `-LS` presence of aggressors inhibits cell fission, i.e. cells that are close to an agressor cell will not try to fission and settle empty land
 - `-LSd dist` maximum distance used for the above effect (if not given 10x the aggressor attack characteristic distance is used)
 - `-Lt ratio` if a cells' population falls below this after an attack (relative to its carrying capacity), its inhabitants become refugees (set to zero in all of our analysis)
 - `-La` if set, aggressors cannot attack refuge cells (not used)
 - `-Ls` if set, unsuccesful attacks also contribute to the landscape of fear effect (by default, only successful attacks have an effect)
 - `-LB prob` increase the defensive capacity of refuge cells by the given probability
 - `-LR` aggressors are unaffected by conflict in their neighborhood


## sp2png / spout2png.cpp

This program converts the detailed output of the simulation run among arbitrary regions to a series of PNG images that can be used to create a video. It also supports directly running ffmpeg to create a video output file.

Parameters:

 - `-ow num` width of the output images (note that images are scaled in a way so that the whole area fits in the requested size)
 - `-oh num` height of the output images
 - `-op num` output images in this intervals (e.g. 10 means every 10 years)
 - `-c fn` or `-cc fn` name of the file with a list of cell IDs, coordinates and (optionally) base carrying capacities; by default, this file is expected in a TSV format without header; use the `-cc` form if the input is a CSV with header
 - `-n fn` or `-ne fn` name of the file containing the neighbor relations among cells; use `-n` if this file is in the DGGRID output format (all neighbors of a cell are listed on a line), and the `-ne` variant if this file is an edgelist (each line contains two node IDs); the input file should be symmetric, i.e. it should contain an edge both ways
 - `-cp fn` name of the file containing the actual cell polygons in a simple CSV format (note: needs to be ordered in a way that corresponds to how they should be displayed)
 - `-om num` maximum population value expected in the input (needs to be given since the data is not read in advance)
 - `-og num` gamma correction to apply to the output color scale (essentially scales output values when selecting the color of cells)
 - `-o fn` base output filename when creating a set of PNG images; image numbers are appended based on time step
 - `-of fn` output video filename; if this option is used, this program will try to run `ffmpeg` directly to create a video
 - `-oc codec` change the video codec used (only matters when using `ffmpeg`, i.e. the `-of` option); this should be something that `ffmpeg` supports, the default is `vp9` (setting this to `h264` can often result in faster encoding); additional options can also be given as part of this
 - `oq num` change the video encoding quality (only matters when using `ffmpeg`, i.e. the `-of` option), `num` should be a small integer (default: 15); lower numbers result in better quality and larger file sizes; setting a too low quality (a high number) will result in noticable issues, e.g. blurry or choppy videos
 

Currently, the color scale used is hardcoded (goes from dark to light blue, bandits are displayed in red), but this should be made configurable in the future. If the `-of` option is used, no PNG images are used. The `ffmpeg` executable should be in PATH to run it successfully. It is given the following options:
```
ffmpeg -y -f rawvideo -pix_fmt bgra -s %ux%u -r 24 -i - -codec:v vp9 -crf 15 %s
```
(where the size and the output filename are filled in based on the options)


