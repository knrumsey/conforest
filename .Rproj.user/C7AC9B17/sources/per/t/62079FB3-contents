#' Conformal Random Forest with Out-of-bag KNN
#'
#' A conformal RF based on Johansson et al. (2014). Rather than using a traditional calibration set, conformal intervals are constructed using out-of-bag samples for each tree. K-nearest neighbors is used to calibrate the prediction intervals for each prediction location.
#'
#' @param X a data frame or a matrix of predictors
#' @param y a response vector
#' @param k Number of nearest neighbors to use for out-of-bag error estimaters
#' @param beta Sensitivity parameter
#' @param eps Tolerance for nearest neighbor approach. When \code{eps = 0}, exact nearest neighbors are used.
#' @param ... additional parameters passed to randomForest
#' @return An object with class "rfok"
#' @references Johansson, U., Boström, H., Löfström, T., & Linusson, H. (2014). Regression conformal prediction with random forests. Machine learning, 97, 155-176.
#' @examples
#' # Friedman function
#' f <-  function(x){
#'   10 * sin(pi * x[1] * x[2]) + 20 * (x[3] - 0.5)^2 + 10 * x[4] + 5 * x[5]
#' }
#' X <- matrix(runif(250), nrow=50, ncol=5)
#' y <- apply(X, 1, f)
#' fit <- rfok(X, y)
#'
#' Xnew <- matrix(runif(250), nrow=50, ncol=5)
#' predict(fit, Xnew)
#' @importFrom graphics par abline segments
#' @importFrom stats predict quantile rbinom runif sd
#' @export
rfok <- function(X, y, k=5, beta=sd(y)/30, eps=0, ...){
  n <- nrow(X)

  # Exact nn search?
  if(eps == 0){
    searchtype = "standard"
  }else{
    searchtype = "priority"
  }

  # Fit RF model
  fit <- randomForest::randomForest(X, y, keep.inbag=TRUE, ...)

  # Get predictions
  preds <- predict(fit, newdata=X, predict.all=TRUE)
  yhat <- preds$aggregate
  preds <- preds$individual

  # Get out of bag predictions
  out_of_bag_sets <- apply(fit$inbag, 1, function(xx) which(xx==0))
  yhat_oob <- rep(NA, n)
  for(i in 1:n){
    yhat_oob[i] <- mean(preds[i, out_of_bag_sets[[i]]])
  }

  # Get k nearest neighbors for each point
  neighbors <- RANN::nn2(X, X, k=min(k, nrow(X)), searchtype=searchtype, eps=eps)$nn.idx

  # Estimate the mus
  mu <- rep(NA, n)
  for(i in 1:n){
    ind <- neighbors[i,]
    mu[i] <- mean(abs(y[ind] - yhat_oob[ind]))
  }

  # Calculate non-conformity scores
  alpha <- abs(y - yhat_oob)/(mu + beta)
  #alpha_s <- quantile(alpha, conf)

  # Return object
  object <- list(fit=fit, k=k, beta=beta, alpha=alpha, oob_error=yhat_oob-y, X=X, y=y)
  class(object) <- "rfok"
  return(object)
}


#' Predict method for rfok
#'
#' A conformal RF based on Johansson et al. (2014). Rather than using a traditional calibration set, conformal intervals are constructed using out-of-bag samples for each tree. K-nearest neighbors is used to calibrate the prediction intervals for each prediction location.
#'
#' @param object Object returned by \code{rfok()}
#' @param newdata a data frame or matrix of new data
#' @param samples Number of samples from the predictive distribution.
#' @param conf a vector of desired confidence intervals. Ignored unless samples is NULL.
#' @param smoothing logical; should conformity scores be sampled or smoothed (using quantiles)
#' @param eps Tolerance for nearest neighbor approach. When \code{eps = 0}, exact nearest neighbors are used.
#' @param ... Additional, ignored, parameters
#' @return An object with class "rfok"
#' @examples
#' # Friedman function
#' f <-  function(x){
#'   10 * sin(pi * x[1] * x[2]) + 20 * (x[3] - 0.5)^2 + 10 * x[4] + 5 * x[5]
#' }
#' X <- matrix(runif(250), nrow=50, ncol=5)
#' y <- apply(X, 1, f)
#' fit <- rfok(X, y)
#'
#' Xnew <- matrix(runif(250), nrow=50, ncol=5)
#' predict(fit, Xnew)
#' @importFrom graphics par abline segments
#' @importFrom stats predict quantile rbinom runif sd
#' @export
predict.rfok <- function(object, newdata=NULL, samples=1000, conf=NULL, smoothing=TRUE, eps=0, ...){
  if(is.null(newdata)){
    newdata <- object$X
  }

  # Exact nn search?
  if(eps == 0){
    searchtype = "standard"
  }else{
    searchtype = "priority"
  }


  pred <- predict(object$fit, newdata)
  n <- length(pred)

  if(is.null(conf) & is.null(samples)){
    return(pred)
  }

  # Get nearest neighbors
  neighbors <- RANN::nn2(object$X, newdata, k=min(object$k, nrow(object$X)),
                         searchtype=searchtype, eps=eps)$nn.idx

  # Estimate mus
  mu <- rep(NA, n)
  for(i in 1:n){
    ind <- neighbors[i,]
    mu[i] <- mean(abs(object$oob_error[ind]))
  }

  if(is.null(samples)){
    # Generate predictions and confidence intervals
    moe <- matrix(NA, nrow=length(pred), ncol=length(conf))
    colnames(moe) <- gsub("\\.", "_", paste(conf))
    for(i in seq_along(conf)){
      moe[,i] <- quantile(object$alpha, conf[i])*(mu + object$beta)
    }
    out <- list(pred=pred, moe=moe)
    return(out)
  }else{
    # Generate predictive samples
    preds <- matrix(NA, nrow=samples, ncol=n)
    for(i in 1:n){
      if(smoothing){
        alpha_delta <- quantile(c(object$alpha, -object$alpha), runif(samples))
      }else{
        alpha_delta <- sample(object$alpha, samples, replace=TRUE) * (2*rbinom(samples,1,0.5)- 1)
      }
      error_term <- alpha_delta*(mu[i] + object$beta)
      preds[,i] <- pred[i] + error_term
    }
  }
  return(preds)
}

#' Plot method for rfok
#'
#' A conformal RF based on Johansson et al. (2014). Rather than using a traditional calibration set, conformal intervals are constructed using out-of-bag samples for each tree. K-nearest neighbors is used to calibrate the prediction intervals for each prediction location.
#'
#' @param x Object returned by \code{rfok()}
#' @param conf Confidence level for plotting
#' @param ... Optional parameters to pass to \code{plot()}
#' @return An object with class "rfok"
#' @examples
#' # Friedman function
#' f <-  function(x){
#'   10 * sin(pi * x[1] * x[2]) + 20 * (x[3] - 0.5)^2 + 10 * x[4] + 5 * x[5]
#' }
#' X <- matrix(runif(250), nrow=50, ncol=5)
#' y <- apply(X, 1, f)
#' fit <- rfok(X, y)
#'
#' plot(fit)
#' @importFrom graphics par abline segments
#' @importFrom stats predict quantile rbinom runif sd
#' @export
plot.rfok <- function(x, conf=0.95, ...){
  object <- x
  a <- 1-conf
  xx <- object$X
  pp <- predict(object)
  yhat <- colMeans(pp)
  bounds <- apply(pp, 2, quantile, probs = c(a/2, 1-a/2))

  hold_par <- par(no.readonly = TRUE)
  par(mfrow=c(1,2))
  plot(object$y, yhat, ylim=range(pp), main="Predictions", ...)
  abline(0, 1, lwd=2, col='orange')
  segments(object$y, bounds[1,], object$y, bounds[2,])

  plot(object$fit, main="MSE")
  par(hold_par)
}




