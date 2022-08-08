#' Function to normalize raw Ct values using a global mean method.
#' @param D KxJ matrix of raw Ct values for miRNAs across samples, where row k corresponds to k-th miRNA, and column j corresponds to j-th sample. Missing values can be represented with NA.
#' @return
#'  A list containing
#'  A J-length vector of normalizers for each sample \code{u_j}
#'  A KxJ matrix of normalized Ct values \code{D_normalized}
#' @export

get_global_mean_norm<-function(D){
  KK = nrow(D)
  E_Dj <- apply(D,2, mean,na.rm=T)
  E_D <- mean(E_Dj,na.rm = T)
  u_j <- E_Dj-E_D
  D_normalized <- D-rep(1,KK)%*%t(u_j)
  return(list(u_j=u_j,D_normalized=D_normalized))
}

#' Function to normalize raw Ct values using a global fit method.
#' @param D KxJ matrix of raw Ct values for miRNAs across samples, where row k corresponds to k-th miRNA, and column j corresponds to j-th sample. Missing values can be represented with NA.
#' @return
#'  A list containing
#'  A Jx3 matrix of parameter estimate: slope,bias,dropout \code{estimate}
#'  A J-length vector of normalizers for each sample \code{u_j}
#'  A KxJ matrix of normalized Ct values \code{D_normalized}
#' @export

get_global_fit_norm<-function(D){
  KK<-nrow(D)
  JJ<-ncol(D)
  drop_out<-colSums(is.na(D))
  D_min <- min(D,na.rm=T)
  D_max <- max(D,na.rm=T)
  ct_lb<-max(5,D_min-5)
  ct_rb<-max(45,D_max+5)
  ###starting value####
  slope_start<-0.5
  drop_start<-mean(drop_out,na.rm=T)
  bias_start<-0.5*(ct_lb+ct_rb)
  ##
  y_median<- KK - 0.5 * (KK - mean(drop_out))
  X_j_median<-rep(NA,JJ)
  param_mt<-matrix(rep(NA,JJ*3),nrow = JJ)
  for (j in 1:JJ) {
    X_j<-seq(ct_lb,ct_rb)
    D_j<-D[,j]
    y_j<- KK- sapply(X_j, function(x){sum(D_j<x,na.rm = TRUE)})
    ##Regression
    fitdata<-coef(nlsLM(y_j ~ (KK-c)/(1+ exp(a * (X_j-b))) + c, start=list(a=slope_start,b=bias_start, c=drop_start)))
    a<-fitdata[1]
    b<-fitdata[2]
    c<-fitdata[3]
    param_mt[j,]<-c(a,b,c)
    X_j_median[j]<- b + 1/a * log((KK-c)/(y_median - c) - 1)
  }
  X_median<-mean(X_j_median,na.rm=T)
  u_j<-X_j_median - X_median
  D_normalized <- D-rep(1,KK)%*%t(u_j)
  return(list(estimate=param_mt,u_j=u_j,D_normalized=D_normalized))
}


#' Function to conduct correction (secondary normalization) by linear regression normalized Ct matrix.
#' @param D KxJ matrix of normalized Ct values for miRNAs across samples, where row k corresponds to k-th miRNA, and column j corresponds to j-th sample. Missing values can be represented with NA.
#' @param u_j a J length vector of normalizers across J samples obtained from global mean/fit normalization
#' @param non_missing_thres the minimum number of samples required for conducting linear regression correction
#' @return
#'  A list containing
#'  A KxJ matrix of corrected/secondary normalized Ct values \code{D_normalized}
#'  A A Kx2 matrix of parameter estimate for fitted simple linear regression models for K miRNAs
#' @export

lr_correct_norm<-function(D,u_j,non_missing_thres=5){
  KK<-nrow(D)
  JJ<-ncol(D)
  param_mt<-matrix(rep(0,2*KK),ncol = 2)
  for (k in 1:KK) {
    non.na.counts<-sum(!is.na(D[k,]))
    if (non.na.counts<non_missing_thres) {
      next
    }
    param_mt[k,]<-coefficients(lm(as.numeric(D[k,]) ~ u_j))
    D[k,]<-D[k,]-u_j*param_mt[k,2]
  }
  return(list(D_normalized=D,param_mt=param_mt))
}



# 
# get_global_fit_norm_new<-function(D){
#   KK<-nrow(D)
#   JJ<-ncol(D)
#   drop_out<-colSums(is.na(D))
#   D_min <- min(D,na.rm=T)
#   D_max <- max(D,na.rm=T)
#   ct_lb<-max(5,D_min-5)
#   ct_rb<-max(45,D_max+5)
#   ###starting value####
#   slope_start<-0.5
#   drop_start<-mean(drop_out,na.rm=T)
#   bias_start<-0.5*(ct_lb+ct_rb)
#   ##
#   y_median<- KK - 0.5 * (KK - drop_out)
#   #print(y_median)
#   X_j_median<-rep(NA,JJ)
#   param_mt<-matrix(rep(NA,JJ*3),nrow = JJ)
#   for (j in 1:JJ) {
#     X_j<-seq(ct_lb,ct_rb)
#     D_j<-D[,j]
#     y_j<- KK- sapply(X_j, function(x){sum(D_j<x,na.rm = TRUE)})
#     ##Regression
#     fitdata<-coef(nlsLM(y_j ~ (KK-c)/(1+ exp(a * (X_j-b))) + c, start=list(a=slope_start,b=bias_start, c=drop_start)))
#     a<-fitdata[1]
#     b<-fitdata[2]
#     c<-fitdata[3]
#     param_mt[j,]<-c(a,b,c)
#     X_j_median[j]<- b + 1/a * log((KK-c)/(y_median[j] - c) - 1)
#   }
#   X_median<-mean(X_j_median,na.rm=T)
#   u_j<-X_j_median - X_median
#   D_normalized <- D-rep(1,KK)%*%t(u_j)
#   return(list(estimate=param_mt,u_j=u_j,D_normalized=D_normalized))
# }
# 


