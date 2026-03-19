#' Run the MCMC for the Bayesian model
#'
#' The main function to run the MCMC for the Bayesian model, to obtain individual dynamics, model parameters such as infection probability, boosting, waning, and measurement error.
#' @param inputdata The data for running MCMC, in dataframe format. It should be in the same format as the data in the package. It includes: 1) age_group (0: children, 1: adults, 2: older adults), 2) start_time: start of follow-up, 3) end_time: end of follow-up, 4) time1: date for first serum collection, 5) time2: date for second serum collection, 6) time3: date for third serum collection, 7) HAI_titer_1: HAI titer for first serum collection, 8) HAI_titer_2: HAI titer for second serum collection, 9) HAI_titer_3: HAI titer for third serum collection.
#' @param inputILI The data for influenza activity used in the inference. The row number should match with the date in the inputdata.
#' @param n_iteration The number of iterations of the MCMC.
#' @param burnin The iteration for burn-in for MCMC.
#' @param thinning The number of thinning in MCMC.
#' @return A list object stores: 1) posterior samples for the model parameter, 2) posterior samples for the parameter for baseline HAI titer, 3) posterior samples for the infection status for each individual, 4) posterior samples for the infection time for each individual (0 for uninfected individuals), 5) posterior samples for the individual waning, 6) posterior samples for individual boosting, 7) posterior samples for the baseline HAI titer, 8) input data, 9) input influenza activity data.
#' @examples
#' \donttest{
#' a1 <- sero_reconstruct(inputdata, flu_activity,
#'                         n_iteration = 2000, burnin = 1000, thinning = 1)
#' }
#' @export
sero_reconstruct <- function(inputdata,inputILI,n_iteration = 2000,burnin = 1000,thinning = 1){
aaaaa1 <- Sys.time()

keep_iteration <- burnin + 1:((n_iteration - burnin)/thinning)*thinning

## first ensure the order
inputdata <- inputdata[,c("age_group", "start_time", "end_time", "time1", "time2", "time3", "HAI_titer_1", "HAI_titer_2", "HAI_titer3")]
## add the column for hhID and member
inputdata <- cbind(0,0,inputdata)
## each of the time need to -14 to represent the boosting delay
inputdata[,4:8] <- inputdata[,4:8]-14
inputdata[is.na(inputdata)] <- -1
inputdata <- cbind(inputdata,2,0,2)
inputILI[inputILI<0] <- 1e-11

inputdata <- as.matrix(inputdata)
inputILI <- as.matrix(inputILI)

#################################################

## try simple format, i.e. gamma with second parameter = 1, beta with second parameter = 1
## model parameter
# 1-18  1.random,  2.1-fold, 3-2fold error*6
# 19-42     1-4 children, adult, older adult boost and waning
# 43-63     infection para      
# 64-75    HAI protection
## first create a function for simulation

int_para <-  c(0.005,rep(0.6,17),rep(c(3.5,0.5),12),rep(c(0.4,0.2,0.2),7),rep(-0.1,12))

# children, adults
int_para2 <- c(0.82,0.08,0.065,0.005,0.005,0.005,0.005,0.005,0.005,0.005,               
               0.87,0.04,0.04,0.015,0.01,0.005,0.005,0.005,0.005,0.005,
               0.82,0.08,0.065,0.005,0.005,0.005,0.005,0.005,0.005,0.005,               
               0.87,0.04,0.04,0.015,0.01,0.005,0.005,0.005,0.005,0.005,
               0.82,0.08,0.065,0.005,0.005,0.005,0.005,0.005,0.005,0.005,               
               0.87,0.04,0.04,0.015,0.01,0.005,0.005,0.005,0.005,0.005,
               0.82,0.08,0.065,0.005,0.005,0.005,0.005,0.005,0.005,0.005,               
               0.87,0.04,0.04,0.015,0.01,0.005,0.005,0.005,0.005,0.005,
               0.82,0.08,0.065,0.005,0.005,0.005,0.005,0.005,0.005,0.005,               
               0.87,0.04,0.04,0.015,0.01,0.005,0.005,0.005,0.005,0.005,
               0.82,0.08,0.065,0.005,0.005,0.005,0.005,0.005,0.005,0.005,               
               0.87,0.04,0.04,0.015,0.01,0.005,0.005,0.005,0.005,0.005) 

### run this to create a template to generate the input for MCMC
t <- sim_data(inputdata,inputILI,int_para,int_para2)

input1 <- t[[1]]
input2 <- t[[2]]
input3 <- t[[3]]

## seting prior
sigma <- abs(int_para)/10
move <- rep(1,length(int_para))
## seting initial value

int_para3 <- rep(1,12)
sigma3 <- int_para3/10
## the original
move[c(4:18,35:42,64:65,67,69,71,73,75)] <- 0
## further for 1 season
move[c(3,19:26,31:51,55:67,69:75)] <- 0
int_para[64:65] <- 0

## season paralist to increase the speed
paraseason <- c(rep(0,42),rep(1,6),rep(2:6,each=3),rep(1:6,each=2))


tt <- mcmc(input1,inputdata,input3,inputILI,n_iteration,int_para,int_para2,int_para3,paraseason,move,sigma,sigma3,burnin,thinning)
aaaaa2 <- Sys.time()

output <- list( 'posterior_model_parameter' = tt[[1]][keep_iteration,which(move==1)],
                'posterior_baseline_HAI_titer' = tt[[2]][keep_iteration,41:60],
                'posterior_inf_status' = tt[[13]],
                'posterior_inf_time' = tt[[14]],
                'posterior_waning' = tt[[15]],
                'posterior_boosting' = tt[[16]],
                'posterior_baseline_titer' = tt[[17]],
                'data' = inputdata,
                'ILI_data' = inputILI)

print('A list containing MCMC results is created; please use the relevant function to obtain the needed information.')
print(paste0('The running time is ',round(difftime(aaaaa2,aaaaa1,units = "secs")), ' seconds'))

return(output)
}

#################################################
## function to compute the mcmc output
## internal use, not output
para_summary <- function(mcmc,a,b,print){
  y <- matrix(NA,ncol(mcmc),4)
  for (i in 1:ncol(mcmc)){
    y[i,1:3] <- quantile(mcmc[,i],c(0.5,0.025,0.975),na.rm=T)
    y[i,4] <- sum(diff(mcmc[,i])!=0)/nrow(mcmc)
    y[i,1] <- mean(mcmc[,i],na.rm=T)
  }
  layout(matrix(1:(a*b),nrow=a,byrow=T))
  par(mar=c(2,4,1,1))
  if (print==1){
    for (i in 1:ncol(mcmc)){
      plot(mcmc[,i],type="l")
    }
  }
  return(y)
}

#################################################
#' Extract the model estimates from the fitted MCMC
#'
#' The function to obtain the estimates of the infection probabilities, boosting, waning and measurement error
#' @param fitted_MCMC The fitted MCMC
#' @param period A vector indicate the start and the end of a season to compute the infection probabilities. If empty, the start and the end of a season was assumed to be the minimum of the start of follow-up and the maximum of the end of follow-up among individuals.
#' @return A data frame that store the estimates from a fitted MCMC.
#' @examples
#' \donttest{
#' a1 <- sero_reconstruct(inputdata, flu_activity,
#'                         n_iteration = 2000, burnin = 1000, thinning = 1)
#' fitted_result <- output_model_estimate(a1)
#' }
#' @export
output_model_estimate <- function(fitted_MCMC,period){
  a1 <- fitted_MCMC
  z1 <- para_summary(a1$posterior_model_parameter,4,3,0)
  z1[1,] <- z1[1,]*10*100
  z1[2,] <- 0.25/exp(z1[2,])*100
  z1[3:6,] <- 2^(z1[3:6,])
  
  
  ILI <- a1$ILI_data
  mcmc1 <- a1$posterior_model_parameter
  mcmc2 <- cbind(a1$posterior_baseline_HAI_titer,a1$posterior_baseline_HAI_titer[,11:20])
  
  xvec <- 0:100
  
  d1list <- list(NA)
  d1 <- matrix(NA,nrow(mcmc1),101)
  
  ## create the period to compute
  if (missing(period)){
  indrow <- min(a1$data[,"start_time"]):max(a1$data[,"end_time"])
  }
  else{
    indrow <- period[1]:period[2]
  }
  
  # compute the risk of the period
    for (u in 1:3){
    for (i in 1:nrow(mcmc1)){
      d1[i,] <- 1-exp(-(mcmc1[i,6+u]*sum(ILI[indrow,1]))*exp(xvec/10*mcmc1[i,10]))
    }
      d1list[[u]] <- d1
    }
  
  ## here need to generate the computation of attack rate
  # 1:2 overall, 3:4 < 1:10
  # 1,2 risk, 2,4 ratio
  plotmatrix <- matrix(NA,4,9)


  mean_ci_function <- function(comp_vec){
    c(mean(comp_vec),quantile(comp_vec,c(0.025,0.975)))
  }
  
  vec <- list(NA)
  for (uu in 1:3){
    vec[[uu]] <- rowSums(d1list[[uu]][,c(6+0:9*10)]*mcmc2[,1:10+10*(u-1)])
  }
    for (v in 0:2){
    plotmatrix[1,1:3+3*v] <- mean_ci_function(vec[[v+1]])
    plotmatrix[2,1:3+3*v] <- mean_ci_function(vec[[v+1]]/vec[[2]])
    plotmatrix[3,1:3+3*v] <- mean_ci_function(d1list[[v+1]][,1])
    plotmatrix[4,1:3+3*v] <- mean_ci_function(d1list[[v+1]][,1]/d1list[[2]][,1] )
    }

  
  output <- data.frame( matrix(NA,16,4) )
  output[1:6,] <- z1[1:6,c(4,1:3)]
  output[7:16,2:4] <- rbind(plotmatrix[,1:3],plotmatrix[c(1,3),4:6],plotmatrix[,7:9])[c(1,5,7,3,6,9,2,4,8,10),]
  
  names(output) <- c("Variable","Point estimate","Lower bound","Upper bound")
  output[,1] <- c("Random error (%)",
                  "Two-fold error (%)",
                  "Fold-increase after infection for children (Boosting)",
                  "Fold-decrease after 1 year for children (Waning)",
                  "Fold-increase after 1 year for adults (Boosting)",
                  "Fold-decrease after 1 year for adults (Waning)",
                  "Infection probability for children",
                  "Infection probability for adults",
                  "Infection probability for older adults",
                  "Infection probability for children with pre-epidemic HAI titer < 10",
                  "Infection probability for adults with pre-epidemic HAI titer < 10",
                  "Infection probability for older adults with pre-epidemic HAI titer < 10",
                  "Relative risk for children (Ref: Adults)",
                  "Relative risk for older adults (Ref: Adults)",
                  "Relative risk for children with pre-epidemic HAI titer < 10 (Ref: Adults with pre-epidemic HAI titer < 10)",
                  "Relative risk for older adults with pre-epidemic HAI titer < 10 (Ref: Adults with pre-epidemic HAI titer < 10)")
  
  output_print <- output
  output_print[,-1] <- round(  output_print[,-1],2 )
  print(output_print)
  
  return(output)
}

#################################################
#' Simulation of the dataset of the Bayesian model
#'
#' The function to simulate the dataset, for validation or other purpose. 
#' @param inputdata The data with the same format that for running MCMC, in dataframe format. 
#' @param inputILI The data for influenza activity used in the inference. The row number should match with the date in the inputdata.
#' @param para1 The parameter vector for the model parameters, in the following format: 1) the parameter value of the random measurement error, 2) the parameter value of the 2-fold error, 3) the boosting in HAI titer after infection for children (in log2 unit), 4) the waning in HAI titer for children (in log2 unit), 5) the boosting in HAI titer after infection for adults (in log2 unit), 6) the waning in HAI titer for adults (in log2 unit), 7) the scale parameter for children, 8) the scale parameter for adults, 9) the scale parameter for older adults, 10) the log of risk ratio of 2-fold increase in baseline HAI titer 
#' @param para2 The parameter vector for the parameter fot the HAI titer distribution. The first 10 elements are the probability that the HAI titer is equal to 0-9 for children, and the elements of 11-20 is the probability that the HAI titer is equal to 0-9 for adults.
#' @return A simulated data based on the input parameter vectors, with the format equal to the input data.
#' @examples
#' simulated <- simulate_data(inputdata, flu_activity, para1, para2)
#' @export
simulate_data <- function(inputdata,inputILI,para1,para2){

## first ensure the order
inputdata <- inputdata[,c("age_group", "start_time", "end_time", "time1", "time2", "time3", "HAI_titer_1", "HAI_titer_2", "HAI_titer3")]
## add the column for hhID and member
inputdata <- cbind(0,0,inputdata)
## each of the time need to -14 to represent the boosting delay
inputdata[,4:8] <- inputdata[,4:8]-14
inputdata[is.na(inputdata)] <- -1
inputdata <- cbind(inputdata,2,0,2)
inputILI[inputILI<0] <- 1e-11

inputdata <- as.matrix(inputdata)
inputILI <- as.matrix(inputILI)

#################################################

## try simple format, i.e. gamma with second parameter = 1, beta with second parameter = 1
## model parameter
# 1-18  1.random,  2.1-fold, 3-2fold error*6
# 19-42     1-4 children, adult, older adult boost and waning
# 43-63     infection para      
# 64-75    HAI protection
## first create a function for simulation

int_para <-  c(0.005,rep(0.6,17),rep(c(3.5,0.5),12),rep(c(0.4,0.2,0.2),7),rep(-0.1,12))
int_para[64:65] <- 0
move <- rep(1,length(int_para))
## the original
move[c(4:18,35:42,64:65,67,69,71,73,75)] <- 0
## further for 1 season
move[c(3,19:26,31:51,55:67,69:75)] <- 0
int_para[which(move==1)] <- para1

# children, adults
int_para2 <- para2

### run this to create a template to generate the input for MCMC
t <- sim_data(inputdata,inputILI,int_para,int_para2)

return(data.frame(t[[2]][,2+1:9]))
}


