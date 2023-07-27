#!/usr/bin/fish
# re-run simulations with 100 realizations,
# using the original processing methodology

# !! TODO: set the correct working directory !!
cd /mnt/smallext/kondor
# !! TODO: set this to the base data directory !!
set data_dir .
# !! TODO: set this to the directory with the simulation code !!
set code_dir /home/kondor/holosim_cpp/neolithic_cpp
set out_base lof
# !! TODO: set this to the maximum number of CPU cores to use !!
set ncores 50

# simulation time
set tmax 5000
# simulation realizations
set N 100


# 1. explore different L and S values
set outdir $out_base/runs_rep1
mkdir -p $outdir

echo G,A,E,L,S,i,seed > $outdir/seeds_run1.csv
set d1 (date)

for G in 80
for A in 10
set pa2 (math 1/$A)
for E in 5
cat $data_dir/matrix_dg_G"$G"zC10.bin > /dev/null
for L in 0.2 1 5
for S in 0.7 1
set fnbase $outdir/res_G"$G"_C10_E"$E"_A"$A"_S"$S"_L$L
for i in (seq $N)
set seed (head -c 4 /dev/urandom | hexdump -e "\"%u\"")
echo $G,$A,$E,$L,$S,$i,$seed >> $outdir/seeds_run1.csv
echo $code_dir/n2ws -Hf $data_dir/matrix_dg_G"$G"zC10.bin -LI -LS -Ld 0 -LD -Lr 0.01 -Lb $L -Lt 0 -Rp 1 -RP -RA $pa2 -R 8 -RR 0.01 -Rc -Rp 1 -Rs $S -RS -d 1 -s $seed -S $tmax -E $E -a 0.5 -o - -op 1 -i 1596250 -I 5000 -If 0.5 -Ip 0.5 -Ce 0 -Cb 0 -De 1 -Dm 0.5 -Db 0.25 -Dl 200 \| $code_dir/spas -m new_grid_aggr_flat_r500.csv -f new_dgid_aggr_filter.dat \| gzip -c \> $fnbase"_"$i.out.gz
end
end
end | $code_dir/rm -t $ncores

date
echo \n
for L in 0.2 1 5
for S in 0.7 1
set fnbase $outdir/res_G"$G"_C10_E"$E"_A"$A"_S"$S"_L$L
echo Rscript --verbose --no-save acf_new_aggr_rep_orig.r $fnbase $N
end
end | /home/kondor/holosim_cpp/neolithic/rm -t 2

rm $outdir/*out.gz

echo $d1
date

end
end
end




# 2. explore different R values
set outdir $out_base/runs_rep2
mkdir -p $outdir

echo G,A,E,L,S,R,i,seed > $outdir/seeds_run1.csv

set d2 (date)

for G in 80
for A in 10
set pa2 (math 1/$A)
for E in 5
for L in 1
cat $data_dir/matrix_dg_G"$G"zC10.bin > /dev/null
for S in 1
for R in 0.01 0.02 0.04 0.06 0.08 0.1
set fnbase $outdir/res_G"$G"_C10_E"$E"_A"$A"_S"$S"_L"$L"_R$R
for i in (seq $N)
set seed (head -c 4 /dev/urandom | hexdump -e "\"%u\"")
echo $G,$A,$E,$L,$S,$R,$i,$seed >> $outdir/seeds_run1.csv
echo $code_dir/n2ws -Hf $data_dir/matrix_dg_G"$G"zC10.bin -LI -LS -Ld 0 -LD -Lr $R -Lb $L -Lt 0 -Rp 1 -RP -RA $pa2 -R 8 -RR 0.01 -Rc -Rp 1 -Rs $S -RS -d 1 -s $seed -S $tmax -E $E -a 0.5 -o - -op 1 -i 1596250 -I 5000 -If 0.5 -Ip 0.5 -Ce 0 -Cb 0 -De 1 -Dm 0.5 -Db 0.25 -Dl 200 \| $code_dir/spas -m new_grid_aggr_flat_r500.csv -f new_dgid_aggr_filter.dat \| gzip -c \> $fnbase"_"$i.out.gz
end
end
end | $code_dir/rm -t $ncores

date
echo \n
for R in 0.02 0.04 0.06 0.08 0.1
for S in 1
set fnbase $outdir/res_G"$G"_C10_E"$E"_A"$A"_S"$S"_L"$L"_R$R
echo Rscript --verbose --no-save acf_new_aggr_rep_orig.r $fnbase $N
end
end | /home/kondor/holosim_cpp/neolithic/rm -t 2

rm $outdir/*out.gz

date

end
end
end
end


# warm-up time period
set tmin 5000
# simulation maximum time
set tmax 10000



# 3. explore different R values, no detrending, use a warm-up period
set outdir $out_base/runs_rep3
mkdir -p $outdir

echo G,A,E,L,S,R,i,seed > $outdir/seeds_run1.csv

set d2 (date)

for G in 80
for A in 10
set pa2 (math 1/$A)
for E in 5
for L in 1
cat $data_dir/matrix_dg_G"$G"zC10.bin > /dev/null
for S in 1
for R in 0.01 0.02 0.04 0.06 0.08 0.1
set fnbase $outdir/res_G"$G"_C10_E"$E"_A"$A"_S"$S"_L"$L"_R$R
for i in (seq $N)
set seed (head -c 4 /dev/urandom | hexdump -e "\"%u\"")
echo $G,$A,$E,$L,$S,$R,$i,$seed >> $outdir/seeds_run1.csv
echo $code_dir/n2ws -Hf $data_dir/matrix_dg_G"$G"zC10.bin -LI -LS -Ld 0 -LD -Lr $R -Lb $L -Lt 0 -Rp 1 -RP -RA $pa2 -R 8 -RR 0.01 -Rc -Rp 1 -Rs $S -RS -d 1 -s $seed -S $tmax -E $E -a 0.5 -o - -op 1 -i 1596250 -I 5000 -If 0.5 -Ip 0.5 -Ce 0 -Cb 0 -De 1 -Dm 0.5 -Db 0.25 -Dl 200 \| mawk \'\{if\(\$1 \>= $tmin\) print \$0\;\}\' \| $code_dir/spas -m new_grid_aggr_flat_r500.csv -f new_dgid_aggr_filter.dat \| gzip -c \> $fnbase"_"$i.out.gz
end
end
end | $code_dir/rm -t $ncores

echo \n
for R in 0.02 0.04 0.06 0.08 0.1
for S in 1
set fnbase $outdir/res_G"$G"_C10_E"$E"_A"$A"_S"$S"_L"$L"_R$R
echo Rscript --verbose --no-save acf_new_aggr_rep_lof.r $fnbase $N
end
end | /home/kondor/holosim_cpp/neolithic/rm -t 2

rm $outdir/*out.gz

end
end
end
end


echo $d1
echo $d2
date



