# Creates the basic ODE model which is applicable to a single time step
# given constant external forcings

buildSingleStepModel <- function(dir) {
  rd <- function(f, ...) { read.table(file=f, header=TRUE, sep="\t", ...) }
  vars <- rd(paste0(dir,"/variables.tsv"))
  prost <- rd(paste0(dir,"/processesAndStoichiometry.tsv"))
  stoi <- as.matrix(prost[,vars[,"name"]])
  rownames(stoi) <- prost[,"name"]
  pros <- prost[,c("name","unit","expression","description")]
  x <- rodeo$new(
    vars= vars,
    pars= rd(paste0(dir,"/parameters.tsv")),
    funs= rd(paste0(dir,"/functions.tsv")),
    pros= pros,
    stoi= stoi,
    asMatrix=TRUE, dim=1
  )
  x$compile(paste0(dir,"/functions.f95"))
  return(x)
}


# Implements the model for multiple time steps; specifically, a series
# of distinct transfers

multiStep <- function(
  n,         # number of transfers
  dt,        # interval between transfers
  pars,      # biological parameters
  initial,   # initial concentrations in system
  inputs,    # external input concentrations
  fracRepl,  # fraction of liquid replaced (with external inputs) in a transfer
  singleStepModel
) {
  if (!identical(paste0(names(initial),"_in"), names(inputs))) {
    print(paste("inputs:", paste(names(inputs), collapse=", ")))
    print(paste("initial:", paste(names(initial), collapse=", ")))
    stop("inputs not consistent with state variables")
  }
  singleStepModel$setPars(pars)
  out <- NULL
  state <- initial
  singleStepModel$setVars(state)
  mobile <- as.logical(singleStepModel$getVarsTable()[,"mobile"])
  for (i in 1:n) {
    if (i > 1) {
      state[mobile] <- state[mobile] * (1-fracRepl) + inputs[mobile] * fracRepl
      singleStepModel$setVars(state)
    }
    dyn <- singleStepModel$dynamics(times=seq(from=(i-1)*dt, to=i*dt, by=3),
      fortran=TRUE, rtol=1e-10, atol=1e-10)
    out <- rbind(out, dyn[2:nrow(dyn),])
    state <- dyn[nrow(dyn), names(state)]
  }
  out[,"time"] <- out[,"time"] / 24
  colnames(out)[colnames(out)=="time"] <- "day"
  out
}
