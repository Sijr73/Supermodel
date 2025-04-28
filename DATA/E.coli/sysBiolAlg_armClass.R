#  sysBiolAlg_armClass.R
#
#  Copyright (C) 2014 Claus Jonathan Fritzemeier
#  Institute for Informatics, Heinrich-Heine-University, Duesseldorf, Germany.
#  All right reserved.
#  Email: clausjonathan.fritzemeier@hhu.de
#
#  Sybil and SybilSWITCH are free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Sybil and SybilSWITCH are distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with sybilSWITCH.  If not, see <http://www.gnu.org/licenses/>.
#
#------------------------------------------------------------------------------#
# A ctive
# R eaction
# M inimization
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#                   definition of the class sysBiolAlg_arm                     #
#------------------------------------------------------------------------------#

setClass(Class = "sysBiolAlg_arm",
		 representation(
             additionalReact  = "numeric"
         ),
         contains = "sysBiolAlg"
)

#------------------------------------------------------------------------------#
# default constructor
#------------------------------------------------------------------------------#

# contructor for class sysBiolAlg_arm
setMethod(f = "initialize",
          signature = "sysBiolAlg_arm",
          definition = function(.Object,
                                model,
                                additionalReact,
                                biomassThreshold=SYBIL_SETTINGS("TOLERANCE")*10,
                                absMAX=SYBIL_SETTINGS("MAXIMUM"),
                                lpdir = "min",
                                scaling = NULL, ...) {

			if ( ! missing(model) ) {
				stopifnot(is(model, "modelorg"),
						is(additionalReact, "numeric")
				)
				
				MAX <- absMAX
				
#				remark to the binary variables:
#				1 -> reaction is ON
#				0 -> reaction is OFF
#                                                                               
#				 _______________________
#				|           |           |                                       
#				|     S     |     0     | = 0                                   
#				|___________|___________|
#  				|1          |MAX        |                                       
#				|  \        |  \        |                                      
#				|    \      |    \      | >= 0                                 
#				|      \    |      \    |                                       
#				|        \  |        \  |                                       
#				|          1|        MAX|                                     
#				|___________|___________|                                     
#				|1          |-MAX       |                                       
#				|  \        |  \        |                                      
#				|    \      |    \      | <= 0                                 
#				|      \    |      \    |                                       
#				|        \  |        \  |                                       
#				|          1|       -MAX|                                       
#				|___________|___________|                                       
#				|           |           |
#				|    1      |           | > biomassThreshold
#				|___________|___________|
#				     ^FBAobj_coef==1
#				                                                                
#		ctype= 	       C           B                                                
#		obj=           0     1---------1                                 

				
				mat <- cbind(S(model), Matrix(0, nrow=met_num(model), ncol=react_num(model)))
				
				
				#mat <- rBind(mat, cBind(diag(react_num(model)), diag(react_num(model))*MAX))
				#mat <- rBind(mat, cBind(diag(react_num(model)), diag(react_num(model))*-MAX))
				
				mat <- rbind(mat, c(obj_coef(model), rep(0, react_num(model))))
				
				fi <- 1:react_num(model)
				
				nRows <- dim(mat)[1]
				nCols <- dim(mat)[2]
				
				obj <- c(rep(0, react_num(model)), rep(1, react_num(model)))
				
				rub <- c(	rep(0, 		met_num(model)),
#							rep(0, 		react_num(model)),
#							rep(0, 		react_num(model)),
							MAX
						)
				rlb <- c(	rep(0,		met_num(model)),
#							rep(0, 	react_num(model)), 
#							rep(0, 	react_num(model)),
							biomassThreshold
						)
				rtype <- c(	rep("E", 	met_num(model)),
				#			rep("L", 	react_num(model)),
				#			rep("U", 	react_num(model)),
								"L"
						)
				
				ub <- c(
						uppbnd(model),
						rep(1, react_num(model))
					   ) 
				lb <- c(lowbnd(model),
						rep(1, react_num(model))
					   ) # all reactions are always on
				lb[react_num(model)+ additionalReact] <- 0 # allow additionalReact to be switched off
				
				ctype <- c(	rep("C", react_num(model)),
							rep("B", react_num(model))
						  )
				
				# use the default constructor for class sysBiolAlg
				.Object <- callNextMethod(.Object,
								        sbalg      = "arm",
								        pType      = "mip",
								        scaling    = scaling,
								        fi         = fi,
								        nCols      = nCols,
								        nRows      = nRows,
								        mat        = mat,
								        ub         = ub,
								        lb         = lb,
								        obj        = obj,
								        rlb        = rlb,
								        rtype      = rtype,
								        lpdir      = lpdir,
								        rub        = rub,
								        ctype      = ctype,
								        algPar     = list(),
								        ...)
				.Object@additionalReact <- additionalReact
				validObject(.Object)
				
				
				#binary = 0 ==> continuus = 0
				
				for(i in 1:react_num(model)){
					e <- addIndConstrCPLEX(.Object@problem@oobj@env, .Object@problem@oobj@lp,
									complemented = TRUE,
									sense = "E",
									rhs = 0,
									indvar = (i+react_num(model))-1,
									nzcnt=1, linind=i-1, linval=1,
									indname = paste0("indConst", i-1))
					if(e !=0){
						stop(paste0("addIndConstrCPLEX had non-zero exit code", e))
					}
				}
				
              }
              return(.Object)
          }
)

suggestedArmSolverSettings <- function(threads=1, timelimit=3600, workMemLimit = 2, treeMemLimit=2){
  
  workMemLimit <- floor(workMemLimit*1024)
  treeMemLimit <- floor(treeMemLimit*1024)
  threads <- as.integer(threads)
  
  sp <- list(	
    CPX_PARAM_EPINT=1e-6,	 		# tolerance to integer values
    CPX_PARAM_EPRHS=1e-6, 			# tolerance to bounds
    CPX_PARAM_TILIM=timelimit, 		# timelimit in seconds
    CPX_PARAM_WORKMEM=workMemLimit, # workspace memory limit
    CPX_PARAM_TRELIM=treeMemLimit, 	# tree memory limit
    CPX_PARAM_THREADS=threads, 		# only use 1 thread
    CPX_PARAM_PARALLELMODE=CPX_PARALLEL_OPPORTUNISTIC # no deterministic search
  )
  
  return(sp)
}

optimizeProbArmLP <- function(model=NULL, sysBiolAlg=NULL, additionalReact=NULL, fixedReact=NULL,
                              react=NULL, ub=NULL, lb=NULL,
                              biomassThreshold=1, absMAX=SYBIL_SETTINGS("MAXIMUM"),
                              randomAfterOptimizationRounds = 4, maxOptimizationRounds = 8,
                              printDebug=T, returnSBA=F, threads=1, restarts=2, algorithm="armLP", solverParm=NULL, ...){
  if(printDebug){
    cat("checking Parameters\n")
  }
  
  if(!is.null(react)){
    if(!is.null(ub)){
      if(length(ub)!=length(react)){
        stop("ub and react have to have the same length!")
      }
    }
    if(!is.null(lb)){
      if(length(lb)!=length(react)){
        stop("lb and react have to have the same length!")
      }
    }
  }
  
  nc <- react_num(model)
  
  lower  <- c(lowbnd(model), rep(0, nc))
  upper  <- c(uppbnd(model), rep(10*absMAX, nc))
  
  if(!is.null(ub)){
    upper[react] <- ub
  }
  if(!is.null(lb)){
    lower[react] <- lb
  }
  
  stopifnot(!is.null(model))
  stopifnot(is(model, "modelorg"))
  
  if(is.null(additionalReact)){
    if(is.null(fixedReact)){
      stop("you have to give either 'additionalReact' or 'fixedReact'")
    }else{
      additionalReact <- setdiff(1:react_num(model), fixedReact)
    }
  }
  
  solverTolerance <- SYBIL_SETTINGS("TOLERANCE")/10
  tolerance <- SYBIL_SETTINGS("TOLERANCE")
  
  
  SYBIL_SETTINGS("SOLVER", "cplexAPI")
  SYBIL_SETTINGS("METHOD", "primopt")
  
  sp <- list(	CPX_PARAM_EPRHS=solverTolerance, 		# tolerance to bounds
              CPX_PARAM_NUMERICALEMPHASIS = CPX_ON,	# more accurate results
              CPX_PARAM_TILIM=60*60, 					# timelimit in seconds
              CPX_PARAM_THREADS=as.integer(threads) 					# only use 1 thread
  )
  if(!is.null(solverParm)){
    sp[names(solverParm)] <- solverParm
  }
  
  if(printDebug){
    cat("calculating Costs\n")
  }
  
  initialcost <- rep(0, react_num(model))
  initialcost[additionalReact] <- 1
  
  
  sba <- NULL
  if(!is.null(sysBiolAlg)){
    if(is(sysBiolAlg, "sysBiolAlg_armLP")){
      sba <- sysBiolAlg
    }
  }else{
    if(printDebug){
      cat(">>> building LP\n")
    }
    sba <- sysBiolAlg(model=model, algorithm=algorithm, biomassThreshold=biomassThreshold, costs=initialcost, absMAX=absMAX, solverParm=sp, ...)
    if(printDebug){
      cat(">>> finished building LP\n")
    }
  }
  
  # initilize run
  cost <- initialcost
  epsilon <- 1
  #lastRes <- NULL
  
  opt <- NULL
  solStat <- NULL
  bestSolution <- list()
  
  if(printDebug){
    cat(">>> starting optimization Rounds\n")
  }
  i <- 1
  while(i <= maxOptimizationRounds){
    if(!printDebug){
      cat(".")
    }
    opt <- optimizeProb(sba, react=1:(nc*2), ub=upper, lb=lower, obj_coef=c(rep(0, nc), cost), resetChanges = FALSE)
    
    solStat <- c(solStat, opt$stat)
    
    if(length(checkSolStat(opt$stat))>0){
      if(printDebug){
        cat(paste0(getMeanStatus(opt$stat), "\n"))
      }
      if(opt$stat!=5){
        rsba <- NULL
        if(returnSBA){
          rsba <- sba
        }
        return(list(solution=NULL, fluxes=NULL, solStat=solStat, sba=rsba))
      }
      else{
        if(i == maxOptimizationRounds && is.null(bestSolution)){
          rsba <- NULL
          if(returnSBA){
            rsba <- sba
          }
          return(list(solution=NULL, fluxes=NULL, solStat=solStat, sba=rsba))
        }
        if(printDebug){
          cat("restarting round\n")
        }
        if(restarts>0){
          i <- (i %/% randomAfterOptimizationRounds * 4)
          restarts <- restarts - 1
        }
        else{
          if(printDebug){
            cat("no more restarts.\n")
          }
        }
      }
    }else{
      # solution was ok!
      if(abs(opt$obj) < tolerance){
        rsba <- NULL
        cat("No active reaction!\n")
        if(returnSBA){
          rsba <- sba
        }
        return(list(solution=NULL, fluxes=NULL, solStat=solStat, sba=rsba))
      }
      
      a <- opt$fluxes[additionalReact]
      res <- sum(abs(a) > tolerance)
      
      if(printDebug){
        cat(res)
        cat("\t")
        cat(epsilon)
      }
      
      if(is.null(bestSolution$res)){
        bestSolution$opt <- opt
        bestSolution$res <- res
      }
      else{
        if(res <= bestSolution$res){
          bestSolution$opt <- opt
          bestSolution$res <- res
        }
        else{
          if(printDebug){
            cat("\tprevious solution was better.")
          }
        }
      }
      #			lastRes <- append(lastRes, res)
    }
    
    
    # recalculate costs
    if(i < maxOptimizationRounds){
      if(i %% randomAfterOptimizationRounds == 0 || opt$stat==5){
        if(printDebug){
          cat("\tset weights randomly and start over.\n")
        }
        cost <- rep(0, nc)
        cost[additionalReact] <- abs(rnorm(length(additionalReact)))+1
        epsilon <- 1
      }else{
        epsilon <- epsilon / 10
        cost <- rep(0, nc)
        f <- opt$fluxes[additionalReact]
        cost[additionalReact] <- 1 / pmax(abs(f), epsilon)
        if(printDebug){
          cat("\n")
        }
      }
    }
    i <- i+1
  }
  if(printDebug){
    cat("\n")
  }
  rsba <- NULL
  if(returnSBA){
    rsba <- sba
  }
  
  if(printDebug){
    cat(">>> checking Parameters\n")
  }else{
    cat("\n")
  }
  
  # fallback if no solution found.
  if(is.null(bestSolution$res)){
    return(list(solution=NULL, fluxes=NULL, solStat=solStat, sba=rsba))
  }
  return(list(solution=which(abs(bestSolution$opt$fluxes[1:react_num(model)]) > tolerance), fluxes=bestSolution$opt$fluxes[1:react_num(model)], solStat=solStat, sba=rsba))
}






















#  sysBiolAlg_armLPClass.R
#
#  Copyright (C) 2014 Claus Jonathan Fritzemeier
#  Institute for Informatics, Heinrich-Heine-University, Duesseldorf, Germany.
#  All right reserved.
#  Email: clausjonathan.fritzemeier@hhu.de
#
#  Sybil and SybilSWITCH are free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Sybil and SybilSWITCH are distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with sybilSWITCH.  If not, see <http://www.gnu.org/licenses/>.
#
#------------------------------------------------------------------------------#
# A ctive
# R eaction
# M inimization
#   with
# L inear
# P rogramming
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#                   definition of the class sysBiolAlg_armLP                   #
#------------------------------------------------------------------------------#

setClass(Class = "sysBiolAlg_armLP",
         representation(
           maxobj = "numeric"
         ),
         contains = "sysBiolAlg"
)


#------------------------------------------------------------------------------#
#                            default constructor                               #
#------------------------------------------------------------------------------#

# contructor for class sysBiolAlg_armLP
setMethod(f = "initialize",
          signature = "sysBiolAlg_armLP",
          definition = function(.Object,
                                model,
                                biomassThreshold = SYBIL_SETTINGS("TOLERANCE")*10,
                                react = NULL, lb = NULL, ub = NULL,
                                costs = NULL,
                                additionalReact  = "numeric",
                                absMAX = SYBIL_SETTINGS("MAXIMUM"),
                                useNames = SYBIL_SETTINGS("USE_NAMES"),
                                cnames = NULL,
                                rnames = NULL,
                                pname = NULL,
                                scaling = NULL,
                                writeProbToFileName = NULL, ...) {
            
            if ( ! missing(model) ) {
              
              if (is.null(biomassThreshold)) {
                stop("no biomassThreshold given")
              }
              
              stopifnot(is(model, "modelorg"), is(biomassThreshold, "numeric"))
              
              if (length(biomassThreshold) != 1) {
                stop("biomassThreshold has to be length 1")
              }
              else {
                maxobj <- NULL
                currmo <- biomassThreshold[1]
              }
              
              
              #  the problem: minimize:
              #
              #            |      |
              #         S  |  0   |  = b
              #            |      |
              #       -------------
              #            |      |
              #         1  |  1   | >= 0
              #            |      |
              #       -------------
              #            |      |
              #         -1 |  1   | >= 0
              #            |      |
              #       -------------
              #       c_wt |  0   |>= c^T * v_wt
              #            |      
              #  lb   wt_lb|  0   |
              #  ub   wt_ub|10000 |
              #            |      |
              #  obj    0  | cost |
              
              
              # ---------------------------------------------
              # problem dimensions
              # ---------------------------------------------
              
              nc     <- react_num(model)
              nr     <- met_num(model)
              
              nCols  <- 2*nc
              nRows  <- nr + 2*nc + 1
              
              
              
              # ---------------------------------------------
              # constraint matrix
              # ---------------------------------------------
              
              # the initial matrix dimensions
              LHS <- Matrix::Matrix(0, 
                                    nrow = nRows,
                                    ncol = nCols,
                                    sparse = TRUE)
              
              # rows for the mutant strain
              LHS[1:nr,1:nc] <- S(model)
              
              # location of the mutant strain
              fi <- c(1:(nc*2))
              
              # rows for the delta match matrix
              diag(LHS[(nr+1)   :(nr+nc)  ,1       :nc    ]) <-  1
              diag(LHS[(nr+1)   :(nr+nc)  ,(nc+1)  :(2*nc)]) <-  1
              diag(LHS[(nr+nc+1):(nr+2*nc),1       :nc    ]) <- -1
              diag(LHS[(nr+nc+1):(nr+2*nc),(nc+1)  :(2*nc)]) <-  1
              
              # fix the value of the objective function
              LHS[(nr+2*nc+1),1:nc] <- obj_coef(model)
              
              
              # ---------------------------------------------
              # columns
              # ---------------------------------------------
              
              lower  <- c(lowbnd(model), rep(0, nc))
              upper  <- c(uppbnd(model), rep(absMAX, nc))
              
              
              # ---------------------------------------------
              # rows
              # ---------------------------------------------
              
              rlower <- c(rep(0, nr), rep(0, 2*nc), currmo)
              rupper <- c(rep(0, nr), rep(absMAX, 2*nc + 1))
              rtype  <- c(rep("E", nr), rep("L", 2*nc + 1))
              
              # ---------------------------------------------
              # objective function
              # ---------------------------------------------
              stopifnot(!is.null(costs))
              cobj <- c(rep(0, nc), costs)
              
              
              # ---------------------------------------------
              # row and column names for the problem object
              # ---------------------------------------------
              
              if (isTRUE(useNames)) {
                if (is.null(cnames)) {
                  cn <- c(react_id(model),
                          paste("cost", react_id(model), sep = "_")
                  )
                  colNames <- sybil:::.makeLPcompatible(cn, prefix = "x")
                }
                else {
                  stopifnot(is(cnames, "character"),
                            length(cnames) == nCols)
                  colNames <- cnames
                }
                
                if (is.null(rnames)) {
                  rn <- c(met_id(model),
                          paste("bw", 1:nc, sep = "_"),
                          paste("fw", 1:nc, sep = "_"),
                          "obj_wt"
                  )
                  rowNames <- sybil:::.makeLPcompatible(rn, prefix = "r")
                }
                else {
                  stopifnot(is(rnames, "character"),
                            length(rnames) == nRows)
                  rowNames <- rnames
                }
                
                if (is.null(pname)) {
                  probName <- sybil:::.makeLPcompatible(
                    paste("armLP", mod_id(model), sep = "_"),
                    prefix = "")
                }
                else {
                  stopifnot(is(pname, "character"),
                            length(pname) == 1)
                  probName <- pname
                }
              }
              else {
                colNames <- NULL
                rowNames <- NULL
                probName <- NULL
              }
              
              
              # ---------------------------------------------
              # build problem object
              # ---------------------------------------------
              
              .Object <- callNextMethod(.Object,
                                        sbalg      = "armLP",
                                        pType      = "lp",
                                        scaling    = scaling,
                                        fi         = fi,
                                        nCols      = nCols,
                                        nRows      = nRows,
                                        mat        = LHS,
                                        ub         = upper,
                                        lb         = lower,
                                        obj        = cobj,
                                        rlb        = rlower,
                                        rub        = rupper,
                                        rtype      = rtype,
                                        lpdir      = "min",
                                        ctype      = NULL,
                                        cnames     = colNames,
                                        rnames     = rowNames,
                                        pname      = probName,
                                        algPar     = list("biomassThreshold" = biomassThreshold,
                                                          "costs" = costs),
                                        ...)
              
              .Object@maxobj <- as.numeric(maxobj)
              
              if (!is.null(writeProbToFileName)) {
                writeProb(problem(.Object),
                          fname = as.character(writeProbToFileName))
              }
            }
            return(.Object)
          }
)


#  sysBiolAlg_armLPClass.R
#
#  Copyright (C) 2014 Claus Jonathan Fritzemeier
#  Institute for Informatics, Heinrich-Heine-University, Duesseldorf, Germany.
#  All right reserved.
#  Email: clausjonathan.fritzemeier@hhu.de
#
#  Sybil and SybilSWITCH are free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Sybil and SybilSWITCH are distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with sybilSWITCH.  If not, see <http://www.gnu.org/licenses/>.
#
#------------------------------------------------------------------------------#
# A ctive
# R eaction
# M inimization
#   with
# L inear
# P rogramming
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#           definition of the class sysBiolAlg_armLPEasyConstraint             #
#------------------------------------------------------------------------------#

setClass(Class = "sysBiolAlg_armLPEasyConstraint",
         representation(
           maxobj = "numeric",
           easyConstraint = "list"
         ),
         contains = "sysBiolAlg"
)


#------------------------------------------------------------------------------#
#                            default constructor                               #
#------------------------------------------------------------------------------#

# contructor for class sysBiolAlg_armLP
setMethod(f = "initialize",
          signature = "sysBiolAlg_armLPEasyConstraint",
          definition = function(.Object,
                                model,
                                biomassThreshold = SYBIL_SETTINGS("TOLERANCE")*10,
                                react = NULL, lb = NULL, ub = NULL,
                                costs = NULL,
                                additionalReact  = "numeric",
                                absMAX = SYBIL_SETTINGS("MAXIMUM"),
                                useNames = SYBIL_SETTINGS("USE_NAMES"),
                                cnames = NULL,
                                rnames = NULL,
                                pname = NULL,
                                scaling = NULL,
                                easyConstraint = NULL,
                                writeProbToFileName = NULL, ...) {
            
            if ( ! missing(model) ) {
              
              if (is.null(biomassThreshold)) {
                stop("no biomassThreshold given")
              }
              
              stopifnot(is(model, "modelorg"), is(biomassThreshold, "numeric"))
              
              if (length(biomassThreshold) != 1) {
                stop("biomassThreshold has to be length 1")
              }
              else {
                maxobj <- NULL
                currmo <- biomassThreshold[1]
              }
              
              
              #  the problem: minimize:
              #
              #            |      |
              #         S  |  0   |  = b
              #            |      |
              #       -------------
              #            |      |
              #         1  |  1   | >= 0
              #            |      |
              #       -------------
              #            |      |
              #         -1 |  1   | >= 0
              #            |      |
              #       -------------
              #       c_wt |  0   |>= c^T * v_wt
              #            |      
              #  lb   wt_lb|  0   |
              #  ub   wt_ub|10000 |
              #            |      |
              #  obj    0  | cost |
              
              
              # ---------------------------------------------
              # problem dimensions
              # ---------------------------------------------
              
              nc     <- react_num(model)
              nr     <- met_num(model)
              
              nCols  <- 2*nc
              nRows  <- nr + 2*nc + 1
              
              
              
              # ---------------------------------------------
              # constraint matrix
              # ---------------------------------------------
              
              # the initial matrix dimensions
              LHS <- Matrix::Matrix(0, 
                                    nrow = nRows,
                                    ncol = nCols,
                                    sparse = TRUE)
              
              # rows for the mutant strain
              LHS[1:nr,1:nc] <- S(model)
              
              # location of the mutant strain
              fi <- c(1:(nc*2))
              
              # rows for the delta match matrix
              diag(LHS[(nr+1)   :(nr+nc)  ,1       :nc    ]) <-  1
              diag(LHS[(nr+1)   :(nr+nc)  ,(nc+1)  :(2*nc)]) <-  1
              diag(LHS[(nr+nc+1):(nr+2*nc),1       :nc    ]) <- -1
              diag(LHS[(nr+nc+1):(nr+2*nc),(nc+1)  :(2*nc)]) <-  1
              
              # fix the value of the objective function
              LHS[(nr+2*nc+1),1:nc] <- obj_coef(model)
              
              
              # ---------------------------------------------
              # columns
              # ---------------------------------------------
              
              lower  <- c(lowbnd(model), rep(0, nc))
              upper  <- c(uppbnd(model), rep(absMAX, nc))
              
              
              # ---------------------------------------------
              # rows
              # ---------------------------------------------
              
              rlower <- c(rep(0, nr), rep(0, 2*nc), currmo)
              rupper <- c(rep(0, nr), rep(absMAX, 2*nc + 1))
              rtype  <- c(rep("E", nr), rep("L", 2*nc + 1))
              
              # ---------------------------------------------
              # objective function
              # ---------------------------------------------
              stopifnot(!is.null(costs))
              cobj <- c(rep(0, nc), costs)
              
              
              # ---------------------------------------------
              # row and column names for the problem object
              # ---------------------------------------------
              
              if (isTRUE(useNames)) {
                if (is.null(cnames)) {
                  cn <- c(react_id(model),
                          paste("cost", react_id(model), sep = "_")
                  )
                  colNames <- sybil:::.makeLPcompatible(cn, prefix = "x")
                }
                else {
                  stopifnot(is(cnames, "character"),
                            length(cnames) == nCols)
                  colNames <- cnames
                }
                
                if (is.null(rnames)) {
                  rn <- c(met_id(model),
                          paste("bw", 1:nc, sep = "_"),
                          paste("fw", 1:nc, sep = "_"),
                          "obj_wt"
                  )
                  rowNames <- sybil:::.makeLPcompatible(rn, prefix = "r")
                }
                else {
                  stopifnot(is(rnames, "character"),
                            length(rnames) == nRows)
                  rowNames <- rnames
                }
                
                if (is.null(pname)) {
                  probName <- sybil:::.makeLPcompatible(
                    paste("armLP", mod_id(model), sep = "_"),
                    prefix = "")
                }
                else {
                  stopifnot(is(pname, "character"),
                            length(pname) == 1)
                  probName <- pname
                }
              }
              else {
                colNames <- NULL
                rowNames <- NULL
                probName <- NULL
              }
              #add easyConstraints:
              if(!is.null(easyConstraint)){
                if(		length(easyConstraint$react) != length(easyConstraint$x)
                     | 	length(easyConstraint$react) != length(easyConstraint$rtype)
                ){
                  stop("easyConstraint elements have to have equal lengths")
                }
                stopifnot(is.list(easyConstraint$react))
                stopifnot(is.list(easyConstraint$x))
                stopifnot(all(easyConstraint$rtype %in% c("F", "L", "U", "D", "E")))
                
                # setting and checking rlb
                if(is.null(easyConstraint$lb)){
                  rlower <- c(rlower, rep(0, length(easyConstraint$react)))
                }else{
                  if(length(easyConstraint$react) != length(easyConstraint$lb)){
                    stop("easyConstraint$lb length has to match length of react argument")
                  }else{
                    stopifnot(is.numeric(easyConstraint$lb))
                    rlower <- c(rlower, easyConstraint$lb)
                  }
                }
                
                # setting and checking rub
                if(is.null(easyConstraint$ub)){
                  rupper <- c(rupper, rep(0, length(easyConstraint$react)))
                }else{
                  if(length(easyConstraint$react) != length(easyConstraint$ub)){
                    stop("easyConstraint$ub length has to match length of react argument")
                  }else{
                    stopifnot(is.numeric(easyConstraint$ub))
                    rupper <- c(rupper, easyConstraint$ub)
                  }
                }
                
                m <- Matrix(0, ncol=nCols, nrow=length(easyConstraint$react))
                
                for(i in 1:length(easyConstraint$react)){
                  m[i, easyConstraint$react[[i]]] <- easyConstraint$x[[i]]
                }
                
                
                LHS <- rbind2(LHS, m)
                rtype <- c(rtype, easyConstraint$rtype)
                nRows <- nRows + length(easyConstraint$react)
                if(!is.null(rowNames)){
                  rowNames <- c(rowNames, paste0("easyConstraint", 1:length(easyConstraint$react)))
                }
                
              }
              
              # ---------------------------------------------
              # build problem object
              # ---------------------------------------------
              
              .Object <- callNextMethod(.Object,
                                        sbalg      = "armLP",
                                        pType      = "lp",
                                        scaling    = scaling,
                                        fi         = fi,
                                        nCols      = nCols,
                                        nRows      = nRows,
                                        mat        = LHS,
                                        ub         = upper,
                                        lb         = lower,
                                        obj        = cobj,
                                        rlb        = rlower,
                                        rub        = rupper,
                                        rtype      = rtype,
                                        lpdir      = "min",
                                        ctype      = NULL,
                                        cnames     = colNames,
                                        rnames     = rowNames,
                                        pname      = probName,
                                        algPar     = list("biomassThreshold" = biomassThreshold,
                                                          "costs" = costs),
                                        ...)
              
              .Object@maxobj <- as.numeric(maxobj)
              .Object@easyConstraint <- easyConstraint
              
              if (!is.null(writeProbToFileName)) {
                writeProb(problem(.Object),
                          fname = as.character(writeProbToFileName))
              }
            }
            return(.Object)
          }
)

#  sysBiolAlg_armLPClass.R
#
#  Copyright (C) 2014 Claus Jonathan Fritzemeier
#  Institute for Informatics, Heinrich-Heine-University, Duesseldorf, Germany.
#  All right reserved.
#  Email: clausjonathan.fritzemeier@hhu.de
#
#  Sybil and SybilSWITCH are free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Sybil and SybilSWITCH are distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with sybilSWITCH.  If not, see <http://www.gnu.org/licenses/>.
#
#------------------------------------------------------------------------------#
# A ctive
# R eaction
# M inimization
#   with
# L inear
# P rogramming
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#           definition of the class sysBiolAlg_armLPEasyConstraint             #
#------------------------------------------------------------------------------#

setClass(Class = "sysBiolAlg_armLPEasyConstraint",
         representation(
           maxobj = "numeric",
           easyConstraint = "list"
         ),
         contains = "sysBiolAlg"
)


#------------------------------------------------------------------------------#
#                            default constructor                               #
#------------------------------------------------------------------------------#

# contructor for class sysBiolAlg_armLP
setMethod(f = "initialize",
          signature = "sysBiolAlg_armLPEasyConstraint",
          definition = function(.Object,
                                model,
                                biomassThreshold = SYBIL_SETTINGS("TOLERANCE")*10,
                                react = NULL, lb = NULL, ub = NULL,
                                costs = NULL,
                                additionalReact  = "numeric",
                                absMAX = SYBIL_SETTINGS("MAXIMUM"),
                                useNames = SYBIL_SETTINGS("USE_NAMES"),
                                cnames = NULL,
                                rnames = NULL,
                                pname = NULL,
                                scaling = NULL,
                                easyConstraint = NULL,
                                writeProbToFileName = NULL, ...) {
            
            if ( ! missing(model) ) {
              
              if (is.null(biomassThreshold)) {
                stop("no biomassThreshold given")
              }
              
              stopifnot(is(model, "modelorg"), is(biomassThreshold, "numeric"))
              
              if (length(biomassThreshold) != 1) {
                stop("biomassThreshold has to be length 1")
              }
              else {
                maxobj <- NULL
                currmo <- biomassThreshold[1]
              }
              
              
              #  the problem: minimize:
              #
              #            |      |
              #         S  |  0   |  = b
              #            |      |
              #       -------------
              #            |      |
              #         1  |  1   | >= 0
              #            |      |
              #       -------------
              #            |      |
              #         -1 |  1   | >= 0
              #            |      |
              #       -------------
              #       c_wt |  0   |>= c^T * v_wt
              #            |      
              #  lb   wt_lb|  0   |
              #  ub   wt_ub|10000 |
              #            |      |
              #  obj    0  | cost |
              
              
              # ---------------------------------------------
              # problem dimensions
              # ---------------------------------------------
              
              nc     <- react_num(model)
              nr     <- met_num(model)
              
              nCols  <- 2*nc
              nRows  <- nr + 2*nc + 1
              
              
              
              # ---------------------------------------------
              # constraint matrix
              # ---------------------------------------------
              
              # the initial matrix dimensions
              LHS <- Matrix::Matrix(0, 
                                    nrow = nRows,
                                    ncol = nCols,
                                    sparse = TRUE)
              
              # rows for the mutant strain
              LHS[1:nr,1:nc] <- S(model)
              
              # location of the mutant strain
              fi <- c(1:(nc*2))
              
              # rows for the delta match matrix
              diag(LHS[(nr+1)   :(nr+nc)  ,1       :nc    ]) <-  1
              diag(LHS[(nr+1)   :(nr+nc)  ,(nc+1)  :(2*nc)]) <-  1
              diag(LHS[(nr+nc+1):(nr+2*nc),1       :nc    ]) <- -1
              diag(LHS[(nr+nc+1):(nr+2*nc),(nc+1)  :(2*nc)]) <-  1
              
              # fix the value of the objective function
              LHS[(nr+2*nc+1),1:nc] <- obj_coef(model)
              
              
              # ---------------------------------------------
              # columns
              # ---------------------------------------------
              
              lower  <- c(lowbnd(model), rep(0, nc))
              upper  <- c(uppbnd(model), rep(absMAX, nc))
              
              
              # ---------------------------------------------
              # rows
              # ---------------------------------------------
              
              rlower <- c(rep(0, nr), rep(0, 2*nc), currmo)
              rupper <- c(rep(0, nr), rep(absMAX, 2*nc + 1))
              rtype  <- c(rep("E", nr), rep("L", 2*nc + 1))
              
              # ---------------------------------------------
              # objective function
              # ---------------------------------------------
              stopifnot(!is.null(costs))
              cobj <- c(rep(0, nc), costs)
              
              
              # ---------------------------------------------
              # row and column names for the problem object
              # ---------------------------------------------
              
              if (isTRUE(useNames)) {
                if (is.null(cnames)) {
                  cn <- c(react_id(model),
                          paste("cost", react_id(model), sep = "_")
                  )
                  colNames <- sybil:::.makeLPcompatible(cn, prefix = "x")
                }
                else {
                  stopifnot(is(cnames, "character"),
                            length(cnames) == nCols)
                  colNames <- cnames
                }
                
                if (is.null(rnames)) {
                  rn <- c(met_id(model),
                          paste("bw", 1:nc, sep = "_"),
                          paste("fw", 1:nc, sep = "_"),
                          "obj_wt"
                  )
                  rowNames <- sybil:::.makeLPcompatible(rn, prefix = "r")
                }
                else {
                  stopifnot(is(rnames, "character"),
                            length(rnames) == nRows)
                  rowNames <- rnames
                }
                
                if (is.null(pname)) {
                  probName <- sybil:::.makeLPcompatible(
                    paste("armLP", mod_id(model), sep = "_"),
                    prefix = "")
                }
                else {
                  stopifnot(is(pname, "character"),
                            length(pname) == 1)
                  probName <- pname
                }
              }
              else {
                colNames <- NULL
                rowNames <- NULL
                probName <- NULL
              }
              #add easyConstraints:
              if(!is.null(easyConstraint)){
                if(		length(easyConstraint$react) != length(easyConstraint$x)
                     | 	length(easyConstraint$react) != length(easyConstraint$rtype)
                ){
                  stop("easyConstraint elements have to have equal lengths")
                }
                stopifnot(is.list(easyConstraint$react))
                stopifnot(is.list(easyConstraint$x))
                stopifnot(all(easyConstraint$rtype %in% c("F", "L", "U", "D", "E")))
                
                # setting and checking rlb
                if(is.null(easyConstraint$lb)){
                  rlower <- c(rlower, rep(0, length(easyConstraint$react)))
                }else{
                  if(length(easyConstraint$react) != length(easyConstraint$lb)){
                    stop("easyConstraint$lb length has to match length of react argument")
                  }else{
                    stopifnot(is.numeric(easyConstraint$lb))
                    rlower <- c(rlower, easyConstraint$lb)
                  }
                }
                
                # setting and checking rub
                if(is.null(easyConstraint$ub)){
                  rupper <- c(rupper, rep(0, length(easyConstraint$react)))
                }else{
                  if(length(easyConstraint$react) != length(easyConstraint$ub)){
                    stop("easyConstraint$ub length has to match length of react argument")
                  }else{
                    stopifnot(is.numeric(easyConstraint$ub))
                    rupper <- c(rupper, easyConstraint$ub)
                  }
                }
                
                m <- Matrix(0, ncol=nCols, nrow=length(easyConstraint$react))
                
                for(i in 1:length(easyConstraint$react)){
                  m[i, easyConstraint$react[[i]]] <- easyConstraint$x[[i]]
                }
                
                
                LHS <- rbind2(LHS, m)
                rtype <- c(rtype, easyConstraint$rtype)
                nRows <- nRows + length(easyConstraint$react)
                if(!is.null(rowNames)){
                  rowNames <- c(rowNames, paste0("easyConstraint", 1:length(easyConstraint$react)))
                }
                
              }
              
              # ---------------------------------------------
              # build problem object
              # ---------------------------------------------
              
              .Object <- callNextMethod(.Object,
                                        sbalg      = "armLP",
                                        pType      = "lp",
                                        scaling    = scaling,
                                        fi         = fi,
                                        nCols      = nCols,
                                        nRows      = nRows,
                                        mat        = LHS,
                                        ub         = upper,
                                        lb         = lower,
                                        obj        = cobj,
                                        rlb        = rlower,
                                        rub        = rupper,
                                        rtype      = rtype,
                                        lpdir      = "min",
                                        ctype      = NULL,
                                        cnames     = colNames,
                                        rnames     = rowNames,
                                        pname      = probName,
                                        algPar     = list("biomassThreshold" = biomassThreshold,
                                                          "costs" = costs),
                                        ...)
              
              .Object@maxobj <- as.numeric(maxobj)
              .Object@easyConstraint <- easyConstraint
              
              if (!is.null(writeProbToFileName)) {
                writeProb(problem(.Object),
                          fname = as.character(writeProbToFileName))
              }
            }
            return(.Object)
          }
)

getArmAddedReactions <- function(model=NULL, res=NULL, additionalReact=NULL){
  if(is(res, "optsol_optimizeProb")){
    f <- fluxes(res)
  }
  else{
    f <- res$fluxes
  }
  intersect(which(getArmReactionState(model, f)), additionalReact)
}
getArmReactionFluxes <- function(model, res){
  # returns flux values of reactions
  if(is(res, "optsol_optimizeProb")){
    f <- fluxes(res)
  }
  else{
    f <- res$fluxes
  }
  stopifnot(length(f)==react_num(model)*2)
  f[1:react_num(model)]
}

getArmReactionState <- function(model, f){
  # returns TRUE for reactions switched on
  round(f[react_num(model)+(1:react_num(model))], digits=-log10(SYBIL_SETTINGS("TOLERANCE")))==1
}







