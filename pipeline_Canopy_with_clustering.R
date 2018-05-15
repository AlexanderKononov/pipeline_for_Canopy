#######################################################
#######################################################
#######                                         #######
#######             CNA and SNA input           #######
#######                                         #######
#######################################################
#######################################################
library(Canopy)
arguments <- commandArgs(TRUE)
setwd(arguments[1])
projectname = arguments[2] ## name of project
R <- read.table('Rout.txt', header = T) ## mutant allele read depth (for SNAs)
rownames(R) <- R[,1]
R$mut <- NULL
R <- as.matrix(R)
R <- na.omit(R)

X <- read.table('Xout.txt', header = T) ## total depth (for SNAs)
rownames(X) <- X[,1]
X$mut <- NULL
X <- as.matrix(X)
X <- na.omit(X)

WM <- read.table('WMout.txt', header = T)  ## observed major copy number (for CNA regions)
rownames(WM) <- WM[,1]
WM$CNA_mutation <- NULL
WM <- as.matrix(WM)

Wm <- read.table('Wmout.txt', header = T) ## observed minor copy number (for CNA regions)
rownames(Wm) <- Wm[,1]
Wm$CNA_mutation <- NULL
Wm <- as.matrix(Wm)

epsilonM = 0.01 ## standard deviation of WM, pre-fixed here
epsilonm = 0.01 ## standard deviation of Wm, pre-fixed here
## whether CNA regions harbor specific CNAs (only needed for overlapping CNAs)
C <- read.table('Cout.txt', header = T, check.names = F)
rownames(C) <- C[,1]
C$CNAs <- NULL
C <- as.matrix(C)

Y <- read.table('Yout.txt', header = T, check.names = F) ## whether SNAs are affected by CNAs
rownames(Y) <- Y[,1]
Y$mut <- NULL
Y <- as.matrix(Y)

save.image(file = paste0(projectname,'.Rdata'))
#######################################################
#######################################################
#######                                         #######
#######              SNV Clustering             #######
#######                                         #######
#######################################################
#######################################################
num_cluster=3:6 # Range of number of clusters to run
num_run=15 # How many EM runs per clustering step for each mutation cluster wave
canopy.cluster=canopy.cluster(R = R, X = X, num_cluster = num_cluster, num_run = num_run)

#bic_output=canopy.cluster$bic_output
sna_cluster=canopy.cluster$sna_cluster

save.image(file = paste0(projectname,'.Rdata'))
#######################################################
#######################################################
#######                                         #######
#######      MCMC sampling by clusters          #######
#######                                         #######
#######################################################
#######################################################
#######################################################
K = 2:4 # number of subclones
numchain = 4 # number of chains with random initiations
sampchain = canopy.sample.cluster(R = R, X = X, sna_cluster = sna_cluster, 
                                  WM = WM, Wm = Wm, epsilonM = epsilonM, 
                                  epsilonm = epsilonm, C = C, Y = Y, K = K, 
                                  numchain = numchain, max.simrun = 100000,
                                  min.simrun = 20000, writeskip = 200,
                                  projectname = projectname, cell.line = TRUE,
                                  plot.likelihood = TRUE)

save.image(file = paste0(projectname,'.Rdata'))
#######################################################
#######################################################
#######                                         #######
#######   BIC to determine number of subclones  #######
#######                                         #######
#######################################################
#######################################################

burnin = 50
thin = 5
# If pdf = TRUE, a pdf will be generated.
bic = canopy.BIC(sampchain = sampchain, projectname = projectname, K = K,
                 numchain = numchain, burnin = burnin, thin = thin, pdf = TRUE)
optK = K[which.max(bic)]

save.image(file = paste0(projectname,'.Rdata'))
#######################################################
#######################################################
#######                                         #######
#######         posterior tree evaluation       #######
#######                                         #######
#######################################################
#######################################################
post = canopy.post(sampchain = sampchain, projectname = projectname, K = K,
                   numchain = numchain, burnin = burnin, thin = thin, 
                   optK = optK, C = C, post.config.cutoff = 0.05)
samptreethin = post[[1]]   # list of all post-burnin and thinning trees
samptreethin.lik = post[[2]]   # likelihoods of trees in samptree
config = post[[3]]
config.summary = post[[4]]
print(config.summary)
# first column: tree configuration
# second column: posterior configuration probability in the entire tree space
# third column: posterior configuration likelihood in the subtree space
# note: if modes of posterior probabilities aren't obvious, run sampling longer.

save.image(file = paste0(projectname,'.Rdata'))
#######################################################
#######################################################
#######                                         #######
#######          Tree output and plot           #######
#######                                         #######
#######################################################
#######################################################
# choose the configuration with the highest posterior likelihood
config.i = config.summary[which.max(config.summary[,3]),1]
cat('Configuration', config.i, 'has the highest posterior likelihood.\n')
output.tree = canopy.output(post, config.i, C)
pdf.name = paste(projectname, '_config_highest_likelihood.pdf', sep='')
canopy.plottree(output.tree, pdf = TRUE, pdf.name = pdf.name)

save.image(file = paste0(projectname,'.Rdata'))