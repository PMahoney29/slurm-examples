#! /usr/bin/env Rscript

# get environment variables
MYSCRATCH <- Sys.getenv('MYSCRATCH')
DATADIR <- Sys.getenv('DATADIR')
RESULTDIR <- Sys.getenv('RESULTDIR')
STEPSIZE <- as.numeric(Sys.getenv('STEPSIZE'))
TASKID <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# set defaults if nothing comes from environment variables
MYSCRATCH[is.na(MYSCRATCH)] <- '.'
DATADIR[is.na(DATADIR)] <- '.'
RESULTDIR[is.na(RESULTDIR)] <- '.'
STEPSIZE[is.na(STEPSIZE)] <- 1
TASKID[is.na(TASKID)] <- 0

# Required packages
# Change to fit the modeling requirements.  In this case,
# fitting a step-selection function using clogit.
library(survival)
library(MuMIn)

# Any additional required functions...code for estimating QIC used in SSF
QIC.coxph <- function(mod, details = FALSE)
{
  if (!exists("naive.var", mod))
    stop("QIC can be computed only if robust variances are estimated.")
  trace <- sum(diag(solve(mod$naive.var) %*% mod$var))
  quasi <- mod$loglik[2]
  if (details)
    return(data.frame(QICR = -2 * quasi + 2 * trace, QICI = -2 *
                        quasi + 2 * length(mod$coefficients), QuasiLL = quasi,
                      n = mod$n, nevent = mod$nevent, K = length(mod$coefficients),
                      Trace = trace))
  else return(-2 * quasi + 2 * trace)
}

# get command lines arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1){
  stop("Not enough arguments. Please use args 'listsize', 'prepare', 'run <itemsize>' or 'merge'")
}

## Determine the number of iterations in a loop (mylistsize), typically the length of a matrix,
## array, dataframe or vector.
## In this case, we are importing a set of models to run as a model matrix.  Thus,
## DATADIR should contain two files, the model matrix and the datafile to be analysed
## The model matrix was created in advance using the following code.
#ParamMatrix <- function(params, maxp=length(params)) {
  ## All possible combinations of parameters
  #paramsL <- length(params)
  #Par <- matrix("",ncol=paramsL,nrow=(2^paramsL))
  #j=2
  #for (i in 1:paramsL) {
  #  x <- combinations(paramsL,i,v=params)
  #  Par[j:(j+nrow(x)-1),1:ncol(x)] <- x
  #  j <- j+nrow(x)
  #}
  
  #if (maxp < paramsL) Par <- Par[Par[,(maxp+1)]=="", 1:maxp]  
  #return(Par)
#}
#params <- c('aspectF', 'slope', 'tpi', 'dist_rds')
#MODELMATRIX <- ParamMatrix(params)
#MODELMATRIX[1,1] <- 1  ## If the model has an intercept
#MODELMATRIX[-1,]       ## If the model doesn't have an intercept  (clogit...or rather Coxph does not)

MODELMATRIX <- as.matrix(read.csv(paste0(DATADIR, '/modMatrix.csv'), as.is=T))
mylistsize <- nrow(MODELMATRIX)

# Storage for coefficients
nms <- unique(c(MODELMATRIX))
nms <- nms[-which(nms == "")]
if (!any(nms == "1")) nms <- c("1", nms)  # for intercept (if needed)
oParams <- as.data.frame(matrix(NA, ncol=length(nms)))
names(oParams) <- nms

# get the list size #########
if (args[1] == 'listsize') {
  cat(mylistsize)
}

# execute prepare job ##################
if (args[1] == 'prepare') {
  # Bring in data from a CSV and prep for sending out to the cluster
  inputdata <- read.csv(paste0(DATADIR, '/inDAT_SSF.csv'))
  save(inputdata, file=paste0(MYSCRATCH,'/input.dat'))
  print(paste0('initial data saved to: ', MYSCRATCH, '/input.dat'))
}

# execute parallel job #################################################
if (args[1] == 'run') {
  if (length(args) < 2) {
    stop("Not enough arguments. 'run' needs a second argument 'id'")
  }
  id<-as.numeric(args[2])
  inputdata <- get(load(paste0(MYSCRATCH,'/input.dat')))
  print(paste(Sys.time(), "arrid:" , id, "TASKID:",
              TASKID, "STEPSIZE:", STEPSIZE))
  for (i in (id+TASKID):(id+TASKID+STEPSIZE-1)) {
    if (i > mylistsize) {break}
    
    ####################
    ##### Do stuff
    ####################
      # Build formula from MODELMATRIX row
    model.d <- MODELMATRIX[i,]
    test <- sum(model.d != "")
    model <- as.formula(paste("Surv(Time, Used) ~",paste(model.d[1:test],collapse="+"),
                              '+ strata(Stratum) + cluster(ID)'))
      
      # Run a clogit using coxph()
    m <- coxph(model, data = inputdata) 
    
    # Build data frame
    Output <- data.frame(ModelN = i, Model = as.character(model)[[3]], 
                         AICc = AICc(m), QIC = QIC.coxph(m))
    Output <- cbind(Output, oParams)
    
    # Add coefficients (NA absent, 1's present in the model)
      # could adjust to store actual coefficient values and SEs
    Output[1, which(names(Output) %in% model.d)] <- 1
    
    
    # Print progress
    print(paste(Sys.time(), "i:" , i))
    # save to a temp file and then rename it as last action !
    save(Output, file=paste0(MYSCRATCH,'/run/',i,"-run.dat.tmp"))
    file.rename(paste0(MYSCRATCH,'/run/',i,"-run.dat.tmp"),
                paste0(MYSCRATCH,'/run/',i,"-run.dat"))
  }
}

# merge job ###########################
if (args[1] == 'merge') {
  outputdata <- c()  # Store very basic model output
  incomplete <- c()  # Store models that didn't run
  for (i in 1:mylistsize) {
    print(paste(Sys.time(), "i:" , i))
    
    if (file.exists(paste0(MYSCRATCH,'/run/',i,"-run.dat"))) {
      outputi <- get(load(paste0(MYSCRATCH,'/run/',i,"-run.dat")))
      outputdata <- rbind(outputdata, outputi)}
    else {
      incomplete <- rbind(incomplete, cbind(i, t(MODELMATRIX[i,])))
    }
  }
  mysum <- nrow(outputdata)
  expectsum <- nrow(MODELMATRIX)
  print(paste('result:', mysum, 'expected:', expectsum))
  
  # Export as a CSV, could also be saved as a .dat file as above.
  write.csv(outputdata, paste0(RESULTDIR, '/result.csv'), row.names=F)
  write.csv(incomplete, paste0(RESULTDIR, '/incomplete.csv'), row.names=F)
  print(paste0('saved result to: ', RESULTDIR, '/result.csv'))
}
