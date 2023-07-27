# Basic setup
se xr [0:3000]
se xl 'Time [years]'
se yl 'Population' off 1

# run the simulation (R = 0.1, 5% noise)
!./pm3 -t l -a 0.2 -b 0.03 -d 0.01 -c 1 -r 0.0135 -K 10000 -R 0.1 -x 0.1 -y 100 -z 10 -T 3000 -s 123 -Nx 0.05 -Ny 0.05 -Nz 0.05 > /tmp/pm1.out

se title 'R = 0.1, 5% noise'

# create the base plot
p '/tmp/pm1.out' u 1:(10000*$2) w l lw 3 lc 7 not
rep '/tmp/pm1.out' u 1:3 w l lw 3 lc -1 not
rep '/tmp/pm1.out' u 1:4 w l lw 3 lc 2 not   

se te post eps color solid size 3.2, 2.2
se out '/tmp/pm11.eps'
rep
se out
se te wxt


# re-run the simulation (R = 0.1, no noise, overwrite the previous results)
!./pm3 -t l -a 0.2 -b 0.03 -d 0.01 -c 1 -r 0.0135 -K 10000 -R 0.1 -x 0.1 -y 100 -z 10 -T 3000 > /tmp/pm1.out

# redo the plot
se title 'R = 0.1, without noise'

se te post eps color solid size 3.2, 2.2
se out '/tmp/pm10.eps'
rep
se out
se te wxt


# re-run the simulation (R = 0.01, 5% noise, overwrite the previous results)
!./pm3 -t l -a 0.2 -b 0.03 -d 0.01 -c 1 -r 0.0135 -K 10000 -R 0.01 -x 0.1 -y 100 -z 10 -T 3000 -s 123 -Nx 0.05 -Ny 0.05 -Nz 0.05 > /tmp/pm1.out

se title 'R = 0.01, 5% noise'

rep

# redo the plot
se te post eps color solid size 3.2, 2.2
se out '/tmp/pm21.eps'
rep
se out
se te wxt


# re-run the simulation (R = 0.01, no noise, overwrite the previous results)
!./pm3 -t l -a 0.2 -b 0.03 -d 0.01 -c 1 -r 0.0135 -K 10000 -R 0.01 -x 0.1 -y 100 -z 10 -T 3000 > /tmp/pm1.out

se title 'R = 0.01, without noise'

# redo the plot
se te post eps color solid size 3.2, 2.2
se out '/tmp/pm20.eps'
rep
se out
se te wxt



