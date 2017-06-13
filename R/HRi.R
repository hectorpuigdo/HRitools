devtools::use_package("nls2")
devtools::use_package("ggplot2")
devtools::use_package("rtracklayer")
devtools::use_package("curl")

#' @title HRi estimators
#'
#' @description 
#' It allows you to calculate HRi curvilinear functions over an adapatation-recombination dataset ponderating with 0-fold sites.
#' 
#' @details 
#' This is a first release. It might not work with really complicated datasets, as initial parameters still need to be improved.
#' 
#' @param x Dataframe consisting of recombination, adaptation and 0-fold sites count data.
#' @keywords Hill-Robertson
#' @return A summary of the analysis of the dataframe consisting of: 
#' the curvilinear model, 
#' LHRi of the model, 
#' the cutoff value where LHRi is below 0.05 from the initial LHRi value,
#' the Mean Ideal Adaptation or MIA,
#' the Maximum Detrimental Value or MDV,
#' the Curvature of the function or Recovery Strength,
#' the formula of the model as a Formula object,
#' the remaining LHRi values when recombination increases,
#' a graph with a representation of these values,
#' a general graph with the model
#' and a NLS object with the results of appliying the NLS computation to the data.
#' @export
HRi <- function(x) {
  
  errmsg <- "The input must be a dataframe with three columns: recombination, adaptation and counts on 0-fold sites."
  if (is.data.frame(x)==FALSE) stop(errmsg)
  if (dim(x)[2]!=3) stop(errmsg)
    
  dat <- data.frame(x)
  dat <- dat[ order(dat[,1]), ]
  
  X <- dat[,1]
  Y <- dat[,2]
  Z <- dat[,3]
  
  fo <- as.formula("Y~a-b*exp(-c*X)")
  
  st <- expand.grid(a=seq(0,1,length.out=30), b=seq(0,20,length.out=30), c=seq(0,20,length.out=30))
  first_res<-nls2::nls2(fo, start=st, algorithm="brute-force", control=nls.control(maxiter=length(X)))
  res <- nls2::nls2(formula=fo, start=first_res,control=nls.control(maxiter=200))
  
  A <- summary(res)$coefficients[1]
  B <- summary(res)$coefficients[2]
  C <- summary(res)$coefficients[3]
  A_err <- summary(res)$coefficients[4]
  sigma <- summary(res)$sigma
  
  min_x <- min(X)
  max_x <- max(X)
  
  y <- function(x) {A-B*exp(-C*x)}
  
  min_y <- y(min_x)
  Acorr <- A - min_y
  
  z <- function(x) {(Acorr)-B*exp(-C*x)}
  integral <- function(x) {(Acorr)*x+((B*exp(-C*x))/C)}
  L <- length(X)
  
  Load <- function(X,Y,Z) {
    
    areas <- c()
    size <- max(X) - min(X)
    L <- length(X)
    total_Z <- sum(Z)
    
    for (iter in 100:150) {
      part <- numeric(iter)
      sep <- numeric(iter)
      local_area <- numeric(iter)
      
      for (i in 1:iter) {
        
        sub <- which( X > ((i-1)/iter*size) & X<= (i/iter*size))
        part[i] <- sum(Z[sub])/total_Z
        sep[i] <- i/iter*size
        
        start <- (i-1)/iter*size
        end <- i/iter*size
        
        under <- integral(end)-integral(start)
        total <- Acorr*end - Acorr*start
        
        load <- (total-under)/total
        local_area[i] <- part[i] * load
        
      }
      
      areas <- c(areas,sum(local_area))
    }
    
    return(mean(areas))
  }
  
  loads <- numeric(L)
  for (i in 1:L) loads[i] <- Load(X[i:L],Y[i:L],Z[i:L])
  
  model_v <- c(A,B,C,sigma)
  names(model_v) <- c("a","b","c","sigma")
  
  maxload <- loads[1]
  percentloads <- loads/maxload
  
  ropt005  <- X[which(percentloads<0.05)[1]]
  ropt001  <- X[which(percentloads<0.01)[1]]
  ropt0005 <- X[which(percentloads<0.005)[1]]

  loads_graph <- ggplot2::qplot(X,loads,geom="line") + ggplot2::ylab(expression(L[HRi])) + ggplot2::xlab("Recombination")
  general_graph <- ggplot2::ggplot(dat[,1:2],ggplot2::aes(X,Y)) +
    ggplot2::ylab(expression("Adaptation")) + ggplot2::xlab("Recombination") +
    ggplot2::geom_point(ggplot2::aes(colour=loads), size=I(0.6)) +
    ggplot2::scale_color_gradientn(colours = rainbow(7)) +
    ggplot2::geom_vline(xintercept = ropt005, col=2) +
    ggplot2::geom_vline(xintercept = ropt001,col=3) +
    ggplot2::geom_vline(xintercept = ropt0005,col=4) +
    ggplot2::stat_function(fun=y)
  
  results <- list(model_v,
                  loads[[1]],
                  ropt005,
                  A,B,C,y,
                  loads,loads_graph,
                  general_graph,res)
  
  names(results) <- c("Model",
                      "LHRi","Ropt0.05",
                      "MIA","MDV","Curvature",
                      "Formula","Loads","LoadsGraph",
                      "Graph","nls2Results")
  
  return(results)
  
}

#' @title Linear and curvilinear model comparison
#'
#' @description 
#' It allows you to calculate HRi curvilinear and linear functions and compares both models.
#' 
#' @details 
#' This is a first release. It might not work with really complicated datasets, as initial parameters for curvilinear modeling still need to be improved.
#' 
#' @param x Dataframe consisting of recombination and adaptation, and maybe 0-fold sites count data.
#' @keywords Hill-Robertson
#' @return A list consisting of a data frame with both models' AICs and two objects with both models.
#' @export
LVNLtest <- function(x) {
  
  errmsg <- "The input must be a dataframe with, at least, two columns: recombination and adaptation."
  if (is.data.frame(x)==FALSE) stop(errmsg)
  if (dim(x)[2]<2) stop(errmsg)
  
  dat <- data.frame(x)
  dat <- dat[ order(dat[,1]), ]
  X <- dat[,1]
  Y <- dat[,2]
  
  fo <- as.formula("Y~a-b*exp(-c*X)")
  
  st <- expand.grid(a=seq(0,1,length.out=30), b=seq(0,20,length.out=30), c=seq(0,20,length.out=30))
  first_res<-nls2::nls2(fo, start=st, algorithm="brute-force", control=nls.control(maxiter=length(X)))
  resNL <- nls2::nls2(formula=fo, start=first_res,control=nls.control(maxiter=200))
  
  resL <- lm(Y~X)
  
  aics <- AIC(resL,resNL)
  rownames(aics) <- c("Linear","Curvilinear")
  
  results <- list(aics,resL,resNL)
  names(results) <- c("AIC","Linear model","Curvilinear model")
  return(results)
  
}

#' @title PopFly population Rho per kb values
#' 
#' @param population String: the desired population as encoded in PopFly. The default value is RAL.
#' @param window String: the desired window size. Allowed window sizes: 10kb, 50kb, 100kb, 
#' @keywords PopFly
#' @keywords Recombination
#' @return A data frame with recombination values in rho/kb, as well as the chromosome, start and end positions of each value.
#' @export
rhokbPopFly <- function(population="RAL",window="100kb") {
  
  feature="rhokb"
  template <- "http://popfly.uab.cat/files/"
  subfolder <- "wig/"
  extension <- ".bw"
  
  tmp <- tempfile()
  filename <- paste(c(template,subfolder,population,"_",feature,"_",window,extension),collapse = "")
  curl::curl_download(filename,tmp)
  big <- rtracklayer::import(tmp, format="bw")
  
  result <- as.data.frame(big)
  return(result)
  
}
