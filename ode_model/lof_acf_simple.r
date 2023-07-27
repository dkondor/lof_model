# some plots of model behavior with noise

library(ggplot2)


# set the working directory to where the code is
setwd('/home/dkondor/CSH/HoloSim/lof_model/ode_models')

# run the simulation for the two main cases (run the below commands in the command line):
# ./pm3 -t l -a 0.2 -b 0.03 -d 0.01 -c 1 -r 0.0135 -K 10000 -R 0.1 -x 0.1 -y 100 -z 10 -T 3000 -s 123 -Nx 0.05 -Ny 0.05 -Nz 0.05 > lof_R0.1.out
# ./pm3 -t l -a 0.2 -b 0.03 -d 0.01 -c 1 -r 0.0135 -K 10000 -R 0.01 -x 0.1 -y 100 -z 10 -T 3000 -s 123 -Nx 0.05 -Ny 0.05 -Nz 0.05 > lof_R0.01.out


####################################################################
# ACFs with noise

tmp1 = read.table("lof_R0.01.out", header=FALSE)
names(tmp1) = c("t", "W", "N", "P")
tmp1$sum = tmp1$N + tmp1$P

p1 = ggplot(tmp1) + geom_line(aes(x=t, y=sum), color="black")
p1 = p1 + geom_line(aes(x=t, y=N), color="blue")
p1 = p1 + geom_line(aes(x=t, y=P), color="red")


tmp2 = read.table("lof_R0.1.out", header=FALSE)
names(tmp2) = c("t", "W", "N", "P")
tmp2$sum = tmp2$N + tmp2$P

p2 = ggplot(tmp2) + geom_line(aes(x=t, y=sum), color="black")
p2 = p2 + geom_line(aes(x=t, y=N), color="blue")
p2 = p2 + geom_line(aes(x=t, y=P), color="red")

# do ACF, leave out the first 500 years
tmp1$s2 = tmp1$sum - mean(tmp1$sum[tmp1$t > 500])
acf1 = acf(tmp1$s2[tmp1$t > 500], lag.max = 2000, plot = FALSE)
acf1 = data.frame(lag = acf1$lag, acf = acf1$acf)
tmp2$s2 = tmp2$sum - mean(tmp2$sum[tmp2$t > 500])
acf2 = acf(tmp2$s2[tmp2$t > 500], lag.max = 2000, plot = FALSE)
acf2 = data.frame(lag = acf2$lag, acf = acf2$acf)

p1 = ggplot(acf1) + geom_line(aes(x=lag, y=acf))
p2 = ggplot(acf2) + geom_line(aes(x=lag, y=acf))

p1 = p1 + theme_bw(6) + xlab("Lag [years]") + ylab("ACF")
p2 = p2 + theme_bw(6) + xlab("Lag [years]") + ylab("ACF")

ggsave("lof_acf_R0.01.pdf", p1, width=3.2, height=2)
ggsave("lof_acf_R0.01.png", p1, width=3.2, height=2, dpi=300)
ggsave("lof_acf_R0.1.pdf", p2, width=3.2, height=2)
ggsave("lof_acf_R0.1.png", p2, width=3.2, height=2, dpi=300)


