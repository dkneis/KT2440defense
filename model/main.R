rm(list=ls())

library(rodeo)    # builds R/Fortran source code from tabular model definition 
library(deSolve)  # provides numerical ODE solvers

source("functions.R")  # imports workhorse R functions

# generate the model source code and run the Fortran compiler
singleStepModel <- buildSingleStepModel("./rodeo")

# table of experiments with particular initial or boundary conditions
experiments <- as.matrix(read.table("./experiments.tsv", sep="\t",
  header=TRUE, row.names=1))

# get numeric parameter values as a vector
parameters <- read.table("./rodeo/parameters.tsv", sep="\t", header=TRUE)
parameters <- setNames(parameters[,"default"], parameters[,"name"])

# initialize output container
out <- NULL

# simulate all experiments
for (xps in colnames(experiments)) {
  vNames <- singleStepModel$namesVars()
  sim <- multiStep(
    n=experiments["n", xps],  # number of transfers
    dt=experiments["dt", xps],  # interval between transfers
    pars=parameters,  # biological parameters
    initial=experiments[vNames, xps],  # initial concentrations in system
    inputs=experiments[paste0(vNames,"_in"), xps],  # external input concentrations
    fracRepl=experiments["fracRepl", xps],  # fraction of liquid replaced (with external inputs) in a transfer
    singleStepModel
  )
  # add to output container
  out <- rbind(out, cbind(experiment=xps, as.data.frame(sim)))
}

# very basic plotting example
clr <- c(F0_50="darkorange", F3_50="skyblue", F5_50="steelblue4")
omar <- par("mar")
par(mfrow=c(4,2), mar=c(4.5, 4.5, 0.2, 0.2))
for (v in singleStepModel$namesVars()) {
  plot(range(out[,"day"]), range(out[,v]), type="n", xlab="Day", ylab=v)
  f <- function(x) {
    lines(x[,"day"], x[,v], col=clr[unique(x[,"experiment"])])
  }
  by(out, out[,"experiment"], f)
}
legend("topright", bty="n", lty=1, col=clr, legend=names(clr), title="Experiment")
par(mfrow=c(1,1), mar=omar)

