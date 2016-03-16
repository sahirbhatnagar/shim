#' Gaussian Response fitting function
#'
lspath <- function(x, y, main.effect.names, interaction.names,
               lambda.beta, lambda.gamma,
               weights,
               lambda.factor,
               nlambda.gamma,
               nlambda.beta,
               nlambda,
               threshold, max.iter,
               initialization.type,
               center, normalize, verbose,
               cores) {

  obj <- standardize(x = x, y = y, center = center, normalize = normalize)
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
                               nlambda.beta = nlambda.beta,
                               lambda.factor = lambda.factor)

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

  betas_original_scale <- betas_alphas_original_scale[main.effect.names, , drop = F]
  alphas_original_scale <- betas_alphas_original_scale[interaction.names, , drop = F]

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
              center = center, normalize = normalize,
              nlambda.gamma = nlambda.gamma,
              nlambda.beta = nlambda.beta,
              nlambda = nlambda,
              interaction.names = interaction.names,
              main.effect.names = main.effect.names)
  class(out) <- "lspath"
  return(out)

}



