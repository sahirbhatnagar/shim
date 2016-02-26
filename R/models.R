#' Fit Strong Heredity Interaction Model
#'
#' @description function to fit the Strong Heredity Interaction Model for a
#'   sequence of tuning parameters. This is a penalized regression method that
#'   ensures the interaction term is non-zero only if its corresponding
#'   main-effects are non-zero.
#' @param x Design matrix of dimension \code{n x q}, where \code{n} is the
#'   number of subjects and q is the total number of variables; each row is an
#'   observation vector. This must include all main effects and interactions as
#'   well, with column names corresponding to the names of the main effects
#'   (e.g. \code{x1, x2, E}) and their interactions (e.g. \code{x1:E, x2:E}).
#'   All columns should be scaled to have mean 0 and variance 1; this is done
#'   internally by the \code{\link{shim}} function.
#' @param y response variable (matrix form) of dimension \code{n x 1}
#' @param main.effect.names character vector of main effects names. MUST be
#'   ordered in the same way as the column names of \code{x}. e.g. if the column
#'   names of \code{x} are \code{"x1","x2"} then \code{main.effect.names =
#'   c("x1","x2")}
#' @param interaction.names character vector of interaction names. MUST be
#'   separated by a colon (e.g. x1:x2), AND MUST be ordered in the same way as
#'   the column names of \code{x}
#' @param lambda.beta sequence of tuning parameters for the main effects. If
#'   \code{NULL} (default), this function will automatically calculate a
#'   sequence using the \code{\link{shim_once}} function which will be over a
#'   grid of tuning parameters for gamma as well. If the user specifies a
#'   sequence then this function will not automatically perform the serach over
#'   a grid. You will need to create the grid yourself e.g. repeat the
#'   lambda.gamma for each value of lambda.beta
#' @param lambda.gamma sequence of tuning parameters for the interaction
#'   effects. Default is \code{NULL} which means this function will
#'   automatically calculate a sequence of tuning paramters. See
#'   \code{\link{shim_once}} for details on how this sequence is calculated.
#' @param nlambda.gamma number of tuning parameters for gamma. This needs to be
#'   specified even for user defined inputs
#' @param nlambda.beta number of tuning parameters for beta. This needs to be
#'   specified even for user defined inputs
#' @param nlambda total number of tuning parameters. If \code{lambda.beta =
#'   NULL} and \code{lambda.gamma = NULL} then \code{nlambda} should be equal to
#'   \code{nlambda.beta x nlambda.gamma}. This is important to specify
#'   especially when a user defined sequence of tuning parameters is set.
#' @param threshold Convergence threshold for coordinate descent. Each
#'   coordinate-descent loop continues until the change in the objective
#'   function after all coefficient updates is less than threshold. Default
#'   value is \code{1e-4}.
#' @param max.iter Maximum number of passes over the data for all tuning
#'   parameter values; default is 100.
#' @param cores The number of cores to use for certain calculations in the
#'   \code{\link{shim}} function, i.e. at most how many child processes will be
#'   run simultaneously using the \code{parallel} package. Must be at least one,
#'   and parallelization requires at least two cores. Default is 2.
#' @param initialization.type The procedure used to estimate the regression
#'   coefficients and used in the \code{\link{uni_fun}} function. If
#'   \code{"univariate"} then a series of univariate regressions is performed
#'   with the response variable \code{y}. If \code{"ridge"} then ridge
#'   regression is performed using the \code{\link[glmnet]{cv.glmnet}} function
#'   and the tuning parameter is chosen using 10 fold cross validation. The
#'   default is \code{"ridge"}.
#' @param intercept Should \code{x} and \code{y} be centered. Default is
#'   \code{TRUE}
#' @param normalize Should \code{x} be scaled to have unit variance. Default is
#'   \code{TRUE}
#' @param verbose Should iteration number and vector of length \code{nlambda} be
#'   printed to console? Default is \code{TRUE}. 0 represents the algorithm has
#'   not converged for the pair of tuning parameters lambda.beta and
#'   lambda.gamma and 1 means it has converged
#' @details the index of the tuning parameters is as follows. If for example
#'   there are 10 lambda_gammas, and 20 lambda_betas, then the first
#'   lambda_gamma gets repeated 20 times. So the first twenty entries of tuning
#'   parameters correspond to 1 lambda_gamma and the 20 lambda_betas
#' @return An object with S3 class "shim" \describe{ \item{b0}{Intercept
#'   sequence of length \code{nlambda}} \item{beta}{A nvars x \code{nlambda}
#'   matrix of main effects (\eqn{\beta}) coefficients, stored in sparse column
#'   format \code{("CsparseMatrix")}} \item{alpha}{A nvars x \code{nlambda}
#'   matrix of interaction effects (\eqn{\alpha}) coefficients, stored in sparse
#'   column format \code{("CsparseMatrix")}} \item{gamma}{A nvars x
#'   \code{nlambda} matrix of (\eqn{\gamma}) coefficients, stored in sparse
#'   column format \code{("CsparseMatrix")}} \item{lambda.beta}{The sequence of
#'   tuning parameters used for the main effects} \item{lambda.gamma}{The
#'   sequence of tuning parameters used for the interaction effects}
#'   \item{tuning.parameters}{2 x nlambda matrix of tuning parameters. The first
#'   row corresponds to \code{lambda.beta} and the second row corresponds to
#'   \code{lambda.gamma}} \item{dfbeta}{list of length \code{nlambda} where each
#'   element gives the index of the nonzero \eqn{\beta} coefficients}
#'   \item{dfalpha}{list of length \code{nlambda} where each element gives the
#'   index of the nonzero \eqn{\alpha} coefficients} \item{x}{x matrix }
#'   \item{y}{response data} \item{bx}{column means of x matrix} \item{by}{mean
#'   of response} \item{sx}{column standard deviations of x matrix}
#'   \item{intercept}{logical of if the data was centered prior to fitting}
#'   \item{normalize}{logical of if the data was normalized prior to fitting}
#'   \item{call}{the call to the function} \item{nlambda.gamma}{nlambda.gamma}
#'   \item{nlambda.beta}{nlambda.beta} \item{nlambda}{nlambda}
#'   \item{interaction.names}{interaction names} \item{main.effect.names}{main
#'   effect names} }
#' @note if the user specifies lambda.beta and lambda.gamma then they this will
#'   not take all possible combinations of lambda.beta and lambda.gamma. It will
#'   be the first element of each as a pair, and so on. This is done on purpose
#'   for use with the cv.shim function which uses the same lambda sequences for
#'   each fold.
#' @seealso \code{\link{shim_once}}
#' @examples
#' # number of observations
#' n <- 100
#'
#' # number of predictors
#' p <- 5
#'
#' # environment variable
#' e <- sample(c(0,1), n, replace = T)
#'
#' # main effects
#' x <- cbind(matrix(rnorm(n*p), ncol = p), e)
#'
#' # need to label columns
#' dimnames(x)[[2]] <- c("x1","x2","x3","x4","x5","e")
#'
#' # design matrix without intercept (can be user defined interactions)
#' X <- model.matrix(~(x1+x2+x3)*e+x1*x4+x3*x5-1, data = as.data.frame(x))
#'
#' # names must appear in the same order as X matrix
#' interaction_names <- grep(":", colnames(X), value = T)
#' main_effect_names <- setdiff(colnames(X), interaction_names)
#'
#' # response
#' Y <- X %*% rbinom(ncol(X), 1, 0.6) + 3*rnorm(n)
#'
#' # standardize data
#' data_std <- standardize(X,Y)
#'
#' result <- shim(x = data_std$x, y = data_std$y,
#'             main.effect.names = main_effect_names,
#'             interaction.names = interaction_names)
#' @author Sahir Bhatnagar
#'
#'   Maintainer: Sahir Bhatnagar \email{sahir.bhatnagar@mail.mcgill.ca}
#' @export

shim <- function(x, y, main.effect.names, interaction.names,
                 lambda.beta = NULL, lambda.gamma = NULL,
                 nlambda.gamma = 10,
                 nlambda.beta = 10,
                 nlambda = 100,
                 threshold = 1e-4, max.iter = 100,
                 initialization.type = c("ridge","univariate"),
                 intercept = TRUE, normalize = TRUE, verbose = TRUE,
                 cores = 2) {

  initialization.type <- match.arg(initialization.type)
  this.call <- match.call()
  obj <- standardize(x = x, y = y,intercept = intercept, normalize = normalize)
  x <- obj$x
  y <- obj$y
  bx <- obj$bx
  by <- obj$by
  sx <- obj$sx

  if (is.null(lambda.gamma) & is.null(lambda.beta)) {

    tuning_params <- shim_once(x = x, y = y,
                               main.effect.names = main.effect.names,
                               interaction.names = interaction.names,
                               initialization.type = initialization.type,
                               nlambda.gamma = nlambda.gamma,
                               nlambda.beta = nlambda.beta)

    # convert to a list. each element corresponds to a value of lambda_gamma
    lambda_gamma_list <- rep(lapply(seq_len(length(tuning_params$lambda_gamma)),
                                    function(i) tuning_params$lambda_gamma[i]),
                             each = nlambda.beta)

    lambda_beta_list <- lapply(seq_len(length(unlist(tuning_params$lambda_beta))),
                               function(i) unlist(tuning_params$lambda_beta)[i])
  } else {

    # convert to a list. each element corresponds to a value of lambda_gamma
    # these are already of the proper length i.e., if the user specifies
    # lambda.beta and lambda.gamma then they this will not take all possible
    # combinations of lambda.beta and lambda.gamma. It will be the first element
    # of each as a pair, and so on. This is done on purpose for use with
    # the cv.shim function which uses the same lambda sequences for each fold...
    lambda_gamma_list <- lapply(seq_len(length(lambda.gamma)),
                                function(i) lambda.gamma[i])

    lambda_beta_list <- lapply(seq_len(length(unlist(lambda.beta))),
                               function(i) unlist(lambda.beta)[i])
  }

  adaptive.weights <- ridge_weights(x = x, y = y,
                                    main.effect.names = main.effect.names,
                                    interaction.names = interaction.names)
  adaptive.weights.mat <- replicate(nlambda,adaptive.weights, simplify = "matrix")
  rownames(adaptive.weights.mat) <- rownames(adaptive.weights)
  adaptive_weights_list <- lapply(seq_len(ncol(adaptive.weights.mat)),
                                  function(i) adaptive.weights.mat[,i, drop = F])

  # initialization
  betas_and_alphas <- uni_fun(variables = colnames(x), x = x, y = y,
                              include.intercept = F,
                              type = initialization.type)

  # this converts the alphas to gammas
  uni_start <- convert(betas_and_alphas, main.effect.names = main.effect.names,
                       interaction.names = interaction.names)

  # need to create a matrix here instead of a 1 column vector
  # dim1: # of variables,
  # dim2: # of lambdas

  beta_hat_previous <- replicate(nlambda, uni_start[main.effect.names, , drop = F],
                                 simplify = "matrix")
  rownames(beta_hat_previous) <- main.effect.names

  gamma_hat_previous <- replicate(nlambda, uni_start[interaction.names, , drop = F],
                                  simplify = "matrix")
  rownames(gamma_hat_previous) <- interaction.names

  # convert gamma and beta previous to lists each element corresponds to the
  # coefficients for each combination of lambda_gamma and lambda_beta
  beta_hat_previous_list <- lapply(seq_len(ncol(beta_hat_previous)),
                                   function(i) beta_hat_previous[,i, drop = F])
  gamma_hat_previous_list <- lapply(seq_len(ncol(gamma_hat_previous)),
                                    function(i) gamma_hat_previous[,i, drop = F])

  # store likelihood values at each iteration in a matrix Q
  # piping using magrittr::set_colnames is slower here
  # rows are the iterations, columns are the index of the sequence of
  # lambda_gammas and lambda_betas
  # rows: iteration number
  # columns: tuning parameter
  Q <- matrix(nrow = max.iter + 1, ncol = nlambda)

  # store the value of the likelihood at the 0th iteration
  Q[1,] <- parallel::mcmapply(Q_theta,
                              beta = beta_hat_previous_list,
                              gamma = gamma_hat_previous_list,
                              lambda.beta = lambda_beta_list,
                              lambda.gamma = lambda_gamma_list,
                              weights = adaptive_weights_list,
                              MoreArgs = list(x = x, y = y,
                                              main.effect.names = main.effect.names,
                                              interaction.names = interaction.names),
                              mc.cores = cores)

  m <- 1 # iteration counter
  delta <- 1 # threshold initialization
  # to see which lambdas have converged: 0=not converged, 1=converged
  converged <- rep(0, nlambda)
  y_tilde_list <- vector("list", nlambda)
  x_tilde_list <- vector("list", nlambda)
  gamma_hat_next_list <- vector("list", nlambda)
  y_tilde_2_list_temp <- vector("list", nlambda)
  term_2_temp_list <- vector("list", nlambda)
  term_1_list <- vector("list", nlambda)
  term_2_list <- vector("list", nlambda)


  # for all the x_tilde in zero_x_tilde, return the following matrix with 0 for each coefficient
  # this is like a place holder.
  coef_zero_gamma_matrix <- matrix(data = 0,
                                   nrow = length(interaction.names),
                                   ncol = 1,
                                   dimnames = list(interaction.names))

  # index data.frame to figure out which j < j'
  index <- data.frame(main.effect.names, seq_along(main.effect.names),
                      stringsAsFactors = F)
  colnames(index) <- c("main.effect.names","index")

  while (any(converged == 0) && m < max.iter){

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # update gamma (interaction parameter)
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # this is a nsubjects x lambda matrix for each tuning parameter stored in a list
    # each element of the list corresponds to a tuning parameter
    # need to keep y_tilde_list and x_tilde_list of length nlambda

    not_converged <- which(converged == 0)
    for (j in not_converged) {
      y_tilde_list[[j]] <- y - x[,main.effect.names,drop = F] %*% beta_hat_previous_list[[j]]
    }

    for (j in not_converged) {

      x_tilde_list[[j]] <- xtilde(interaction.names = interaction.names,
                                  data.main.effects = x[,main.effect.names, drop = F],
                                  beta.main.effects = beta_hat_previous_list[[j]])
    }


    # indices of the x_tilde matrices that have all 0 columns
    zero_x_tilde <- which(sapply(x_tilde_list,
                                 function(i) is.null(colnames(check_col_0(i)))))

    # this will store the results but will be shorter than nlambda
    gamma_hat_next_list_not_converged <-
      parallel::mclapply(seq_len(nlambda)[not_converged],
                         function(i) {
                           if (i %in% zero_x_tilde) coef_zero_gamma_matrix else
                             as.matrix(coef(glmnet::glmnet(
                               x = x_tilde_list[[i]],
                               y = y_tilde_list[[i]],
                               penalty.factor = adaptive_weights_list[[i]][interaction.names,,drop=F],
                               lambda = lambda_gamma_list[[i]],
                               standardize = F, intercept = F))[-1,,drop = F])
                         },
                         mc.cores = cores)

    gamma_hat_next_list <- replace(gamma_hat_next_list, not_converged,
                                   gamma_hat_next_list_not_converged)

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # update beta (main effect parameter) step 4 of algortihm in Choi et al
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    beta_hat_next_list <- beta_hat_previous_list

    for (j in main.effect.names) {

      # determine the main effects not in j
      j_prime_not_in_j <- setdiff(main.effect.names,j)

      for (notconverged in not_converged) {
        y_tilde_2_list_temp[[notconverged]] <- y -
          x[,j_prime_not_in_j, drop = F] %*%
          beta_hat_next_list[[notconverged]][j_prime_not_in_j, , drop = F]
      }

      # mclapply is faster than lapply even with just two cores
      term_2_temp_list_not_converged <-
        parallel::mclapply(seq_len(nlambda)[not_converged],
                           function(i)
                             as.matrix(
                               rowSums(
                                 xtilde_mod(beta.main.effects = beta_hat_next_list[[i]][j_prime_not_in_j, , drop = F],
                                            gamma.interaction.effects = gamma_hat_next_list[[i]],
                                            interaction.names = interaction.names[-grep(j, interaction.names)],
                                            data.main.effects = x[,j_prime_not_in_j, drop = F])
                               ),
                               ncol = 1),
                           mc.cores = cores)

      term_2_temp_list <- replace(term_2_temp_list, not_converged, term_2_temp_list_not_converged)

      # this is of length nlambda.beta*nlambda.gamma i.e. one set of y's for each tuning parameter
      y_tilde_2_list <- mapply("-", y_tilde_2_list_temp, term_2_temp_list, SIMPLIFY = F)

      # j' less than j
      j.prime.less <- index[which(index[,"index"] < index[which(index$main.effect.names == j),2]),
                            "main.effect.names"]

      # need to make sure paste(j.prime.less,j,sep=":") are variables in x matrix
      # this is to get around situations where there is only interactions with the E variable
      j.prime.less.interaction <- intersect(paste(j.prime.less,j, sep = ":"), colnames(x))

      # need to get the main effects that are in j.prime.greater.interaction
      j.prime.less <- gsub("\\:(.*)", "", j.prime.less.interaction)


      # the if conditions in term1 and term2 are to check if there are
      # any variables greater or less than j
      # lapply is faster than mclapply here
      term_1_list_not_converged <- if (length(j.prime.less.interaction) != 0 ) {
          lapply(seq_len(nlambda)[not_converged], function(i)
            x[,j.prime.less.interaction] %*%
              (gamma_hat_next_list[[i]][j.prime.less.interaction,, drop = F] *
                 beta_hat_next_list[[i]][j.prime.less, , drop = F]))} else
                   matrix(rep(0,length(beta_hat_next_list[not_converged])), ncol = 1)

      term_1_list <- replace(term_1_list, not_converged, term_1_list_not_converged)

      # j' greater than j
      j.prime.greater <- index[which(index[,"index"] >
                                       index[which(index$main.effect.names == j),2]),
                               "main.effect.names"]

      # need to make sure j.prime.greater is a variable in x matrix
      # this is to get around situations where there is only interactions with the E variable
      j.prime.greater.interaction <- intersect(paste(j,j.prime.greater,sep = ":"), colnames(x))

      # need to get the main effects that are in j.prime.greater.interaction
      j.prime.greater <- if (all(gsub("\\:(.*)", "", j.prime.greater.interaction) == j))
        gsub("(.*)\\:", "", j.prime.greater.interaction) else gsub("\\:(.*)", "", j.prime.greater.interaction)

      term_2_list_not_converged <- if (length(j.prime.greater) != 0) {
        lapply(seq_len(nlambda)[not_converged], function(i)
          x[,j.prime.greater.interaction] %*%
            (gamma_hat_next_list[[i]][j.prime.greater.interaction,, drop = F] *
               beta_hat_next_list[[i]][j.prime.greater,,drop = F])) } else
                 matrix(rep(0,length(beta_hat_next_list[not_converged])), ncol = 1)


      term_2_list <- replace(term_2_list, not_converged, term_2_list_not_converged)


      # lapply is faster than mclapply
      x_tilde_2_list <- lapply(seq_len(length(term_1_list)),
                               function(i) x[,j, drop = F] +
                                 term_1_list[[i]] + term_2_list[[i]])

      # glmnet is giving weired results for this... and is slower than using my
      # soft function. use this. non-parallel version is faster
      # the result of this should give 1 beta for each tuningn parameter
      # This calculates for all tuning parameters
      beta_hat_next_list_j <- lapply(seq_len(length(x_tilde_2_list)), function(i)
        soft(x = x_tilde_2_list[[i]],
             y = y_tilde_2_list[[i]],
             weight = adaptive_weights_list[[i]][j,,drop=F],
             lambda = lambda_beta_list[[i]]))

      # update beta_j for each tuning parameter but only those that
      # have not converged
      for (i in seq_len(nlambda)[not_converged]) {
        beta_hat_next_list[[i]][j,] <- beta_hat_next_list_j[[i]]
      }
    }

    Q[m + 1, not_converged] <- parallel::mcmapply(Q_theta,
                                               beta = beta_hat_next_list[not_converged],
                                               gamma = gamma_hat_next_list[not_converged],
                                               lambda.beta = lambda_beta_list[not_converged],
                                               lambda.gamma = lambda_gamma_list[not_converged],
                                               weights = adaptive_weights_list[not_converged],
                                               MoreArgs = list(x = x, y = y,
                                                               main.effect.names = main.effect.names,
                                                               interaction.names = interaction.names),
                                               mc.cores = cores)


    delta <- abs(Q[m,] - Q[m + 1, ])/abs(Q[m,])
    # if delta is NA, this means Q wasnt calculated for the previous iteration
    # because the algorithm converged, therefore replace with threshold
    # so that it stays as converged
    delta[is.na(delta)] <- threshold
    converged <- as.numeric(delta <= threshold)
    if (verbose) print(paste("Iteration:", m))
    if (verbose) print(converged)

    m <- m + 1

    beta_hat_previous_list <- beta_hat_next_list

    # adaptive weight for each tuning parameter. currently this is the
    # same for iterations, but I am coding it here
    # for flexibility in case we want to change the weights at each iteration

    adaptive_weights_list_not_converged <- lapply(seq_len(nlambda)[not_converged],
                                                  function(i)
                                                    update_weights(betas = beta_hat_previous_list[[i]],
                                                                   gammas = gamma_hat_previous_list[[i]],
                                                                   main.effect.names = main.effect.names,
                                                                   interaction.names = interaction.names))

    adaptive_weights_list <- replace(adaptive_weights_list, not_converged,
                                     adaptive_weights_list_not_converged)


  }

  # convert to original scale
  betas_original_scale_list <- lapply(beta_hat_next_list, function(i) i / sx[main.effect.names])
  gammas_original_scale_list <- lapply(gamma_hat_next_list, function(i) i / sx[interaction.names])


  # convert gammas to alphas
  betas_alphas_original_scale <- mapply(convert2,
                                        beta = betas_original_scale_list,
                                        gamma = gammas_original_scale_list,
                                        MoreArgs = list(main.effect.names = main.effect.names,
                                                        interaction.names = interaction.names))

  dimnames(betas_alphas_original_scale) <- list(c(main.effect.names, interaction.names),
                                                paste0("s",1:nlambda))

  betas_original_scale <- betas_alphas_original_scale[main.effect.names,]
  alphas_original_scale <- betas_alphas_original_scale[interaction.names,]

  b0 <- vector(length = nlambda)
  for (lam in seq_len(nlambda)) {
    b0[lam] <- by - sum(betas_original_scale[,lam,drop = F] * bx[main.effect.names]) -
      sum(alphas_original_scale[,lam,drop=F] * bx[interaction.names])
  }
  names(b0) <- paste0("s",1:nlambda)

  gamma_final <- as(matrix(unlist(gammas_original_scale_list, use.names = F),
                           ncol = nlambda,
                           byrow = T,
                           dimnames = list(interaction.names, paste0("s",1:nlambda))),
                    "dgCMatrix")
  beta_final <- as(betas_original_scale,"dgCMatrix")
  alpha_final <- as(alphas_original_scale,"dgCMatrix")

  lambda.beta <- unlist(lambda_beta_list)
  names(lambda.beta) <- paste0("s",1:nlambda,".beta")
  lambda.gamma <- unlist(lambda_gamma_list)
  names(lambda.gamma) <- paste0("s",1:nlambda, ".gamma")

  tuning.parameters <- matrix(nrow = 2, ncol = nlambda,
                              dimnames = list(c("lambda.beta", "lambda.gamma"),
                                              paste0("s",1:nlambda)))

  tuning.parameters["lambda.beta",] <- lambda.beta
  tuning.parameters["lambda.gamma",] <- lambda.gamma

  out <- list(b0 = b0,
             beta = beta_final,
             alpha = alpha_final,
             gamma = gamma_final,
             lambda.beta = lambda.beta,
             lambda.gamma = lambda.gamma,
             tuning.parameters = tuning.parameters,
             dfbeta = apply(beta_final, 2, eclust::nonzero),
             dfalpha = apply(alpha_final, 2, eclust::nonzero),
             converged = converged, x = x, y = y, bx = bx, by = by, sx = sx,
             intercept = intercept, normalize = normalize, call = this.call,
             nlambda.gamma = nlambda.gamma,
             nlambda.beta = nlambda.beta,
             nlambda = nlambda,
             interaction.names = interaction.names,
             main.effect.names = main.effect.names)
  class(out) <- "shim"
  return(out)
}


#' Cross-validation for shim
#'
#' @description Does k-fold cross-validation for shim and determines the optimal
#'   pair of tuning parameters (\eqn{\lambda_\beta} and \eqn{\lambda_\gamma})
#'
#' @inheritParams shim
#' @param parallel If \code{TRUE}, use parallel \code{foreach} to fit each fold.
#'   Must register parallel before hand using the
#'   \code{\link[doMC]{registerDoMC}} function from the \code{doMC} package. See
#'   the example below for details.
#' @param type.measure loss to use for cross-validation. Currently only 1
#'   option. The default is \code{type.measure="mse"}, which uses squared-error
#'   for gaussian models
#' @param nfolds number of folds - default is 10. Although nfolds can be as
#'   large as the sample size (leave-one-out CV), it is not recommended for
#'   large datasets. Smallest value allowable is \code{nfolds=3}
#' @details The function runs shim nfolds+1 times; the first to get the tuning
#'   parameter sequences, and then the remainder to compute the fit with each of
#'   the folds omitted. The error is accumulated, and the average error and
#'   standard deviation over the folds is computed. Note also that the results
#'   of cv.shim are random, since the folds are selected at random using the
#'   \code{\link{createfolds}} function. Users can reduce this randomness by
#'   running cv.shim many times, and averaging the error curves.
#' @author Sahir Bhatnagar
#'
#'   Maintainer: Sahir Bhatnagar \email{sahir.bhatnagar@mail.mcgill.ca}
#' @export

cv.shim <- function(x, y, main.effect.names, interaction.names,
                    lambda.beta = NULL, lambda.gamma = NULL,
                    nlambda.gamma = 10,
                    nlambda.beta = 10,
                    nlambda = 100,
                    threshold = 1e-4, max.iter = 100,
                    initialization.type = c("ridge","univariate"),
                    intercept = TRUE, normalize = TRUE, verbose = TRUE,
                    parallel = TRUE, type.measure = c("mse"),
                    nfolds = 10) {


  # x = X; y = Y; main.effect.names = main_effect_names;
  # interaction.names = interaction_names;
  # lambda.beta = NULL ; lambda.gamma = NULL
  # threshold = 1e-4 ; max.iter = 500 ; initialization.type = "ridge";
  # nlambda.gamma = 5; nlambda.beta = 20; cores = 1;
  # intercept=TRUE; normalize=TRUE
  #
  # nfolds = 5
  # grouped = TRUE; keep = FALSE; parallel = TRUE

  initialization.type <- match.arg(initialization.type)
  type.measure <- match.arg(type.measure)

  if (!is.null(lambda.beta) && length(lambda.beta) < 2)
    stop("Need more than one value of lambda.beta for cv.shim")
  if (!is.null(lambda.gamma) && length(lambda.gamma) < 2)
    stop("Need more than one value of lambda.gamma for cv.shim")

  N <- nrow(x)
  y <- drop(y)

  shim.object <- shim(x = x, y = y,
                      main.effect.names = main.effect.names,
                      interaction.names = interaction.names,
                      lambda.beta = lambda.beta, lambda.gamma = lambda.gamma,
                      nlambda = nlambda,
                      threshold = threshold, max.iter = max.iter,
                      initialization.type = initialization.type,
                      nlambda.gamma = nlambda.gamma,
                      nlambda.beta = nlambda.beta, cores = 1,
                      verbose = verbose)

  nz.main = sapply(predict(shim.object, type = "nonzero")[["main"]], length)
  nz.interaction = sapply(predict(shim.object, type = "nonzero")[["interaction"]], length)

  foldid <- createfolds(y, k = nfolds)
  if (nfolds < 3)
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  outlist = as.list(seq(nfolds))
  if (parallel) {
    outlist = foreach(i = seq(nfolds), .packages = c("glmnet")) %dopar%
    {
      which = foldid == i
      if (is.matrix(y)) y_sub = y[!which, ] else y_sub = y[!which]
      print(paste("Foldid = ",i))
      shim(x = x[!which, , drop = FALSE],
           y = y_sub,
           main.effect.names = main.effect.names,
           interaction.names = interaction.names,
           lambda.beta = shim.object$lambda.beta,
           lambda.gamma = shim.object$lambda.gamma,
           nlambda = shim.object$nlambda,
           threshold = threshold,
           max.iter = max.iter ,
           initialization.type = initialization.type,
           nlambda.gamma = nlambda.gamma,
           nlambda.beta = nlambda.beta, cores = 1, verbose = verbose)
    }
  } else {
    for (i in seq(nfolds)) {
      which = foldid == i

      if (is.matrix(y))
        y_sub = y[!which, ] else y_sub = y[!which]

      outlist[[i]] = shim(x[!which, , drop = FALSE],
                          y_sub,
                          main.effect.names = main.effect.names,
                          interaction.names = interaction.names,
                          lambda.beta = shim.object$lambda.beta,
                          lambda.gamma = shim.object$lambda.gamma,
                          nlambda = shim.object$nlambda,
                          threshold = threshold,
                          max.iter = max.iter ,
                          initialization.type = initialization.type,
                          nlambda.gamma = nlambda.gamma,
                          nlambda.beta = nlambda.beta, cores = 1,
                          verbose = verbose)
    }
  }

  lambda.beta = shim.object$lambda.beta
  lambda.gamma = shim.object$lambda.gamma

  cvstuff <- do.call(cv_lspath, list(outlist = outlist,
                                    x = x, y = y, foldid = foldid,
                                    nlambda = shim.object$nlambda,
                                    nlambda.beta = shim.object$nlambda.beta,
                                    nlambda.gamma = shim.object$nlambda.gamma))

  cvm <- cvstuff$cvm
  cvsd <- cvstuff$cvsd
  # the following checks is any of the tunining parameter pairs
  # have a cvsd==NA or did not converge. cvstuff$converged should be
  # equal to the number of folds because it is the sum of booleans
  # that converged over all folds.
  nas <- (is.na(cvsd) + (cvstuff$converged != nfolds)) != 0
  if (any(nas)) {
    lambda.beta = lambda.beta[!nas]
    lambda.gamma = lambda.gamma[!nas]
    cvm = cvm[!nas]
    cvsd = cvsd[!nas]
    # this is the total number of non-zero parameters (both betas and alphas)
    nz.main = nz.main[!nas]
    nz.interaction = nz.interaction[!nas]
  }
  cvname <- cvstuff$name

  df <- as.data.frame(cbind(lambda.beta = lambda.beta,
                            lambda.gamma = lambda.gamma,
                            mse = cvm,
                            upper = cvm + cvsd,
                            lower = cvm - cvsd,
                            nz.main = nz.main,
                            nz.interaction = nz.interaction,
                            log.gamma = round(log(lambda.gamma),2)))

  rownames(df) <- gsub("\\.(.*)", "",rownames(df))

  out <- list(lambda.beta = lambda.beta, lambda.gamma = lambda.gamma,
             cvm = cvm, cvsd = cvsd, cvup = cvm + cvsd,
             cvlo = cvm - cvsd, nz.main = nz.main, name = cvname,
             nz.interaction = nz.interaction,
             shim.fit = shim.object, converged = cvstuff$converged,
             cvm.mat.all = cvstuff$cvm.mat.all,
             df = df)

  lamin.beta = if (cvname == "AUC")
    getmin_type(lambda.beta, -cvm, cvsd, type = "beta") else getmin_type(lambda.beta, cvm, cvsd, type = "beta")

  obj <- c(out, as.list(lamin.beta))
  class(obj) <- "cv.shim"
  obj
}



