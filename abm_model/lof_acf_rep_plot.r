# lof_acf_rep_plot.r -- aggregate repeated simulation results
#   create plots

library(ggplot2)

#!! TODO: set this to the base data directory (with results for 14C dates)
data_dir = '~/CSH/HoloSim/data'
# !! TODO: set this to the simulation results
setwd('~/CSH/HoloSim/data/simulation')
base_dir0 = "lof/" # base directory with results

# main simulation parameters -- these do not change 
G = 80
C = 10
E = 5
A = 10


##################################################################
# 1. helper functions for loading and processing data
binwidth = 100 # binwidth for ACF peak distributions
cv_bw = 0.04 # binwidth for CV distributions
get_bin = function(x, bw) {
  return(bw * floor(x / bw) + bw / 2)
}
create_hist = function(peaks1, bw, var) {
  tmp1 = peaks1
  tmp1["bin"] = get_bin(tmp1[var], bw)
  tmp1$count = 1
  tmp1 = aggregate(tmp1["count"], by=tmp1["bin"], FUN=sum)
  scnt = sum(tmp1$count)
  tmp1$p.mean = tmp1$count / scnt
  return(tmp1)
}


# read the results of running multiple realizations of the simulation
# (distribution of ACF minima or CV)
read_histograms = function(fn, r = 500, varname1 = "all_res_min_dt",
                           varname2 = "lag", bw = binwidth, all_res1 = NULL) {
  if(is.na(fn)) {
    if(is.null(all_res1)) stop("Need filename or result object!")
    all_res = all_res1
  }
  else load(fn) # loads variable all_res
  N = length(all_res)
  
  histograms = data.frame()
  for(i in 1:N) {
    tmp1 = all_res[[i]]$res
    peaks = data.frame()
    for(j in 1:length(tmp1)) {
      if(tmp1[[j]]$r %in% r) {
        tmp2 = tmp1[[j]][[varname1]]
        if(length(r) > 1) tmp2$r = tmp1[[j]]$r
        peaks = rbind(peaks, tmp2)
      }
    }
    if(is.null(peaks)) stop('Cannot find result!\n')
    if(nrow(peaks) == 0) next
    peaks$bin = get_bin(peaks[[varname2]], bw)
    peaks$count = 1
    if(length(r) > 1) {
      peaks = aggregate(peaks["count"], by=peaks[c("bin", "r")], FUN=sum)
      sum1 = aggregate(peaks["count"], by=peaks["r"], FUN=sum)
      names(sum1)[names(sum1) == "count"] = "scnt"
      peaks = merge(peaks, sum1, by="r")
      peaks$p = peaks$count / peaks$scnt
    } else {
      peaks = aggregate(peaks["count"], by=peaks["bin"], FUN=sum)
      sum1 = sum(peaks$count)
      peaks$p = peaks$count / sum1
    }
    peaks$i = i
    histograms = rbind(histograms, peaks)
  }
  rm(all_res)
  return(histograms)
}

# aggregate histograms from multiple simulation realizations
aggr_hist = function(histograms) {
  all_i = unique(histograms$i)
  all_bins = unique(histograms$bin)
  tmp1 = data.frame(i = rep(all_i, each=length(all_bins)),
                    bin = rep(all_bins, length(all_i)))
  tmp1 = merge(tmp1, histograms, by=c("i", "bin"), all.x = TRUE)
  tmp1[is.na(tmp1)] = 0
  tmp1 = aggregate(tmp1["p"], by=tmp1["bin"], FUN=function(x) {
    return(c(mean = mean(x), sd = sd(x)))
  })
  tmp1 = do.call(data.frame, tmp1)
  return(tmp1)
}


###################################################################
# 2. load the "baseline" data (results from processing 14C dates,
#   see here:  https://github.com/dkondor/neolithic_simulation/tree/main/analysis)
# ACF minimum distribution
c14_peaks = read.csv(paste0(data_dir, '/population/new_analysis/run2r2/res_min_orig.csv'))
# only use results without a moving average
c14_peaks = c14_peaks[c14_peaks$rm == 0,]
c14_peaks = create_hist(c14_peaks, binwidth, "m")

# CV distribution
c14_cv = read.csv(paste0(data_dir, '/population/new_analysis/run2r2/res_cv.csv'))
c14_cv = c14_cv[c14_cv$rm == 0,]
c14_cv = create_hist(c14_cv, cv_bw, "cv")


# base figures (used later as the background)
p1 = ggplot(mapping = aes(x=bin, y=p.mean)) + theme_bw(6)
p1 = p1 + xlab("Lag [years]") + ylab("Relative frequency")
p1 = p1 + geom_col(data=c14_peaks, color="grey", fill="grey")
p1 = p1 + scale_x_continuous(limits=c(0,2000))

p2 = ggplot(mapping = aes(x=bin, y=p.mean)) + theme_bw(6)
p2 = p2 + xlab("CV") + ylab("Relative frequency")
p2 = p2 + geom_col(data=c14_cv, color="grey", fill="grey")
p2 = p2 + scale_x_continuous(limits=c(0,1.5))



#############################################################
# 3. base results (previous simulation, from our previous paper)
fn1 = 'revisions2/runs_rep2/res_G80_C10_E5_A10_s0_all_res.RData'
tmp1 = read_histograms(fn1)
acf_base = aggr_hist(tmp1)
tmp1 = read_histograms(fn1, 500, "cv", "cv", cv_bw)
cv_base = aggr_hist(tmp1)


# create figures of these
p12 = p1 + geom_col(data=acf_base, color="red", fill="red", alpha=0.4)
p22 = p2 + geom_col(data=cv_base, color="red", fill="red", alpha=0.4)

# save these -- Fig. 3, top two panels (note: label is added manually in Inkscape)
fn1 = paste0(base_dir0, "res_orig")
ggsave(paste0(fn1, "_acf.pdf"), p12, width=1.82, height=1.14)
ggsave(paste0(fn1, "_acf.png"), p12, width=1.82, height=1.14, dpi=300)
ggsave(paste0(fn1, "_cv.pdf"), p22, width=1.82, height=1.14)
ggsave(paste0(fn1, "_cv.png"), p22, width=1.82, height=1.14, dpi=300)



#################################################################x
# 4. Main results: varying S and L (Fig. 3, main panels)
acf2 = data.frame()
cv2 = data.frame()

base_dir = 'lof/runs_rep1/'

for(S in c(0.7,1)) for(L in c(0.2,1,5)) {
  fn0 = paste0(base_dir, "res_G", G, "_C", C, "_E", E, "_A", A,
               "_S", S, "_L", L)
  fn1 = paste0(fn0, "_all_res.RData")
  tmp1 = read_histograms(fn1)
  acf1 = aggr_hist(tmp1)
  tmp1 = read_histograms(fn1, 500, "cv", "cv", cv_bw)
  cv1 = aggr_hist(tmp1)
  
  acf1$L = L
  acf1$S = S
  acf2 = rbind(acf2, acf1)
  cv1$L = L
  cv1$S = S
  cv2 = rbind(cv2, cv1)
}

# create the figures
p12 = p1 + geom_col(data=acf2, color="red", fill="red", alpha=0.4)
p22 = p2 + geom_col(data=cv2, color="red", fill="red", alpha=0.4)
p12 = p12 + facet_grid(L~S, labeller = label_both)
p22 = p22 + facet_grid(L~S, labeller = label_both)

# save these
fn1 = paste0(base_dir, "res_cmb")
ggsave(paste0(fn1, "_acf.pdf"), p12, width=3.2, height=3)
ggsave(paste0(fn1, "_acf.png"), p12, width=3.2, height=3, dpi=300)
ggsave(paste0(fn1, "_cv.pdf"), p22, width=3.2, height=3)
ggsave(paste0(fn1, "_cv.png"), p22, width=3.2, height=3, dpi=300)


###################################################################
# 5. Variation in R (Fig. 4)
S = 1
L = 1
base_dir = "lof/runs_rep2/"

LL = c(1)
acf3 = acf2[acf2$L %in% LL & acf2$S == 1,]
acf3$R = 0.01
cv3 = cv2[cv2$L %in% LL & cv2$S == 1,]
cv3$R = 0.01

acf3 = data.frame()
cv3 = data.frame()

for(R in c(0.01, 0.02, 0.04, 0.06, 0.08, 0.1)) {
  fn0 = paste0(base_dir, "res_G", G, "_C", C, "_E", E, "_A", A,
               "_S", S, "_L", L, "_R", R)
  fn1 = paste0(fn0, "_all_res.RData")
  tmp1 = read_histograms(fn1)
  acf1 = aggr_hist(tmp1)
  tmp1 = read_histograms(fn1, 500, "cv", "cv", cv_bw)
  cv1 = aggr_hist(tmp1)
  
  acf1$L = L
  acf1$S = S
  acf1$R = R
  acf3 = rbind(acf3, acf1)
  cv1$L = L
  cv1$S = S
  cv1$R = R
  cv3 = rbind(cv3, cv1)
}

# figure with all panels
p12 = p1 + geom_col(data=acf3, color="red", fill="red", alpha=0.4)
p22 = p2 + geom_col(data=cv3, color="red", fill="red", alpha=0.4)
p12 = p12 + facet_grid(~R, labeller = label_both)
p22 = p22 + facet_grid(~R, labeller = label_both)

fn1 = paste0(base_dir, "res_cmbR")
w = 6.4; h = 1.8
ggsave(paste0(fn1, "_acf.pdf"), p12, width=w, height=h)
ggsave(paste0(fn1, "_acf.png"), p12, width=w, height=h, dpi=300)
ggsave(paste0(fn1, "_cv.pdf"), p22, width=w, height=h)
ggsave(paste0(fn1, "_cv.png"), p22, width=w, height=h, dpi=300)

# figure with only four panels -> Fig. 4
acf3 = acf3[acf3$R %in% c(0.01, 0.02, 0.06, 0.1),]
cv3 = cv3[cv3$R %in% c(0.01, 0.02, 0.06, 0.1),]

p12 = p1 + geom_col(data=acf3, color="red", fill="red", alpha=0.4)
p22 = p2 + geom_col(data=cv3, color="red", fill="red", alpha=0.4)
p12 = p12 + facet_grid(~R, labeller = label_both)
p22 = p22 + facet_grid(~R, labeller = label_both)

fn1 = paste0(base_dir, "res_cmbR2")
w = 5.6; h = 1.8
ggsave(paste0(fn1, "_acf.pdf"), p12, width=w, height=h)
ggsave(paste0(fn1, "_acf.png"), p12, width=w, height=h, dpi=300)
ggsave(paste0(fn1, "_cv.pdf"), p22, width=w, height=h)
ggsave(paste0(fn1, "_cv.png"), p22, width=w, height=h, dpi=300)



#############################################################
# 6. Fig S2 -- repeat the analysis for comparing different
#  R values, leave out a warm-up period
base_dir = 'lof/runs_rep3/'

# load the data
acf3 = data.frame()
cv3 = data.frame()

for(R in c(0.01, 0.02, 0.04, 0.06, 0.08, 0.1)) {
  fn0 = paste0(base_dir, "res_G", G, "_C", C, "_E", E, "_A", A,
               "_S", S, "_L", L, "_R", R)
  fn1 = paste0(fn0, "_all_res.RData")
  tmp1 = read_histograms(fn1)
  acf1 = aggr_hist(tmp1)
  tmp1 = read_histograms(fn1, 500, "cv", "cv", cv_bw)
  cv1 = aggr_hist(tmp1)
  
  acf1$L = L
  acf1$S = S
  acf1$R = R
  acf3 = rbind(acf3, acf1)
  cv1$L = L
  cv1$S = S
  cv1$R = R
  cv3 = rbind(cv3, cv1)
}

# plot without 14C data under it -- since the simulation time
# is much longer than the historical period with the 14C data,
# it would not make sense to do a direct compatison
p12 = ggplot(mapping = aes(x=bin, y=p.mean)) + theme_bw(6)
p12 = p12 + xlab("Lag [years]") + ylab("Relative frequency")
p12 = p12 + scale_x_continuous(limits=c(0,2000))

p22 = ggplot(mapping = aes(x=bin, y=p.mean)) + theme_bw(6)
p22 = p22 + xlab("CV") + ylab("Relative frequency")
p22 = p22 + scale_x_continuous(limits=c(0,1.5))

p12 = p12 + geom_col(data=acf3, color="red", fill="red", alpha=0.4)
p22 = p22 + geom_col(data=cv3, color="red", fill="red", alpha=0.4)

p12 = p12 + facet_grid(~R, labeller = label_both)
p22 = p22 + facet_grid(~R, labeller = label_both)

# full figure
fn1 = paste0(base_dir, "res_cmbR")
w = 6.4; h = 1.8
ggsave(paste0(fn1, "_acf.pdf"), p12, width=w, height=h)
ggsave(paste0(fn1, "_acf.png"), p12, width=w, height=h, dpi=300)
ggsave(paste0(fn1, "_cv.pdf"), p22, width=w, height=h)
ggsave(paste0(fn1, "_cv.png"), p22, width=w, height=h, dpi=300)

# repeat, filter only four panels
acf3 = acf3[acf3$R %in% c(0.01, 0.02, 0.06, 0.1),]
cv3 = cv3[cv3$R %in% c(0.01, 0.02, 0.06, 0.1),]
fn1 = paste0(base_dir, "res_cmbR2")
w = 5.6; h = 1.8

p12 = ggplot(mapping = aes(x=bin, y=p.mean)) + theme_bw(6)
p12 = p12 + xlab("Lag [years]") + ylab("Relative frequency")
p12 = p12 + scale_x_continuous(limits=c(0,2000))

p22 = ggplot(mapping = aes(x=bin, y=p.mean)) + theme_bw(6)
p22 = p22 + xlab("CV") + ylab("Relative frequency")
p22 = p22 + scale_x_continuous(limits=c(0,1.5))

p12 = p12 + geom_col(data=acf3, color="red", fill="red", alpha=0.4)
p22 = p22 + geom_col(data=cv3, color="red", fill="red", alpha=0.4)

p12 = p12 + facet_grid(~R, labeller = label_both)
p22 = p22 + facet_grid(~R, labeller = label_both)

ggsave(paste0(fn1, "_acf.pdf"), p12, width=w, height=h)
ggsave(paste0(fn1, "_acf.png"), p12, width=w, height=h, dpi=300)
ggsave(paste0(fn1, "_cv.pdf"), p22, width=w, height=h)
ggsave(paste0(fn1, "_cv.png"), p22, width=w, height=h, dpi=300)

