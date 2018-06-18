
library(outbreaker)
library(adegenet)

set.seed(2310)

dna_outb = fasta2DNAbin("tree/sclade_sequence_alignment.fa")

# <as per phybreak etc
time_years = c(0.015, 0.001, 0.273, 0.733, 2.018, 2.36, 2.64, 3.41, 4.878, 4.968, 5.037, 5.834, 6.579)
time_days = 365*time_years
first_date = 2005.956
names(time_years) = c("Case27_2005-12-20", "Case28_2005-12-15", "Case29_2006-03-24", "Case30_2006-09-08", "Case37_2007-12-21", 
                      "Case38_2008-04-24", "Case40_2008-08-04", "Case47_2009-05-12", "Case52_2010-10-30", "Case54_2010-12-02", 
                      "Case55_2010-12-27", "Case68_2011-10-14", "Case81_2012-07-12")
short_names = c("Case27", "Case28", "Case29", "Case30", "Case37", "Case38", "Case40", "Case47", "Case52", "Case54", "Case55", "Case68", "Case81")
very_short_names= c("27", "28", "29", "30", "37", "38", "40", "47", "52", "54", "55", "68", "81")
# as per phybreak etc>

#!WARNING! Initial state of the chain has a likelihood of zero. The chain may never converge.

w = seq(from=0,to=15,by=1) 
density_outb = dgamma(w,shape=1.3,scale=1/0.3)

vegard_outbreak = outbreaker(dna=dna_outb, dates=time_days, 
                     w.dens=density_outb, f.dens=density_outb, 
                     #init.tree="random",
                     #init.mu1 = 1e-4, init.mu2=init.mu1, 
                     n.iter = 10000, sample.every=50, tune.every=50, burnin=2000,# one tenth of defaults
                     #import.method="full", find.import.n=50,
                     pi.prior1 = 10, pi.prior2 = 1, spa1.prior = 1, move.mut = TRUE,
                     move.ances = TRUE, move.kappa = TRUE, move.Tinf = TRUE,
                     move.pi = TRUE, move.spa = TRUE, outlier.threshold = 5,
                     max.kappa = 10, quiet = TRUE, res.file.name = "chains.txt",
                     tune.file.name = "tuning.txt", seed = NULL)
