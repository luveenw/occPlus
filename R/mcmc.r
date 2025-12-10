library(tidyverse)
library(Rcpp)
library(here)
library(MASS)

sourceCpp(here("src","functions.cpp"))
source(here("R","mcmcfun.R"))

# starting info
{
  threshold <- F

  data_info <- as.data.frame(data$info)
  OTU <- data$OTU


  occCovariates = c()
  detCovariates = c()
  ordCovariates = c()

  # data checks
  {
    if(!all(c(occCovariates, ordCovariates, detCovariates) %in% colnames(data$info))){
      stop("Covariate names provided not in data$info")
    }

    if(any(is.na(data_info$Site)) | any(is.na(data_info$Sample)) | any(is.na(data_info$Primer))){
      stop("NA in Site, Sample or Primer columns")
    }
  }

  # sort the data
  {
    data_info <- data_info %>%
      dplyr::arrange(Site, Sample, Primer)
  }

  # samples per site
  {
    M_df <- data_info %>%
      dplyr::group_by(Site, Sample) %>%
      slice(1) %>%
      dplyr::group_by(Site) %>%
      dplyr::summarise(M = n())

    M <- M_df$M
    names(M) <- M_df$Site

    n <- length(M)

    sumM <- c(0, cumsum(M)[-n])

    siteNames <- unique(data_info$Site)

  }

  # marker per samples
  {
    # number of markers

    L_df <- data_info %>%
      dplyr::group_by(Site, Sample, Primer) %>%
      dplyr::slice(1) %>%
      dplyr::group_by(Site, Sample) %>%
      dplyr::summarise(L = n())

    L <- L_df$L
    names(L) <- L_df$Sample

    # number of observations
    {
      L_all <- data_info %>%
        dplyr::group_by(Site, Sample, Primer) %>%
        dplyr::slice(1) %>%
        dplyr::group_by(Site, Sample) %>%
        dplyr::summarise(L_m = n())
      L_all <- L_all$L_m

      sumL <- c(0, cumsum(L_all)[-length(L_all)])

    }

    primerNames <- unique(data_info$Primer)
    # sumL <- c(0, cumsum(L)[-length(L)])

  }

  # pcr per marker
  {
    data_K <- data_info %>%
      dplyr::group_by(Site, Sample, Primer) %>%
      dplyr::summarise(K = n(),
                       dplyr::across(contains("Species"),function(x){sum(x > 0)})
      ) %>%
      dplyr::ungroup()#%>%
    # ungroup() %>%
    # mutate(across(contains("Species"), function(x){sum(x > 0)}))

    K <- data_K$K

    numL <- length(K)

    sumK <- c(0, cumsum(K)[-numL])
  }

  maxL = max(L)

  idx_z <- rep(1:n, M)
  idx_k <- data_info$Sample

  primerIdx <- data_info$Primer

  y <- OTU

  if(threshold != F){

    y[OTU >= threshold] <- 1
    y[OTU < threshold] <- 0

  }

  # priors
  {
    prior_beta_psi <- 0
    prior_beta_psi_sd <- 1
    prior_beta_theta <- 0
    prior_beta_theta_sd <- 1
    a_theta0 <- 1
    b_theta0 <- 20
    a_p <- 5
    b_p <- 1
    a_q <- 1
    b_q <- 20
    a_sigma0 <- 10
    b_sigma0 <- 10
    a_sigma1 <- 1
    b_sigma1 <- 1
    # a_pi0 <- 1
    # b_pi0 <- 1

    # library(MCMCpack)
    #
    # log(dinvgamma(2, 10, 10))

    mu_mu1 <- 3
    sd_mu1 <- 3
    mu_mu0 <- 1
    sd_mu0 <- .5
  }

  # prior
  {
    b_betatheta <- rep(1, ncov_theta)
    B_betatheta <- diag(1, nrow = ncov_theta)

    b_betatheta[1] <- prior_beta_theta
    B_betatheta[1,1] <- prior_beta_theta_sd
  }

}

trueStartingPoint <- F

# starting points
{

  # true starting points
  if (trueStartingPoint){

    z <- z_true
    w <- delta
    theta <- theta_true
    theta0 <- theta0_true
    beta_psi <- beta_psi_true
    beta_ord <- beta_ord_true
    beta_theta <- beta_theta_true
    U <- U_true
    E <- E_true
    LL <- L_true
    p <- p_true
    q <- q_true
    c_imk <- cimk_true
    mu1 <- mu1_true
    sigma1 <- sd1_true
    mu0 <- mu0_true
    sigma0 <- sd0_true

  } else {

    c_imk <- y > 0

    w <- matrix(NA, N, S)

    for (s in 1:S) {
      for (i in 1:n) {
        for (m in 1:M[i]) {
          idx_im1 <- sumK[sumL[sumM[i] + m] + 1] + 1
          idx_im2 <- sumK[sumL[sumM[i] + m] + L[sumM[i] + m]] +
            K[sumL[sumM[i] + m] + L[sumM[i] + m]]
          w[sumM[i] + m,s] <- as.numeric(any(y[idx_im1:idx_im2,s] > 0))
        }
      }
    }

    z <- matrix(NA, n, S)

    for (s in 1:S) {
      for (i in 1:n) {
        idx_i1 <- sumM[i] + 1
        idx_i2 <- sumM[i] + M[i]

        z[i,s] <- as.numeric(any(w[idx_i1:idx_i2, s] > 0))
      }
    }

    beta_psi <- matrix(0, ncov_psi, S)
    beta_theta <- matrix(0, ncov_theta, S)
    beta_ord <- matrix(0, ncov_ord, d)
    E <- matrix(0, n, d)
    LL <- matrix(1, d, S)

    U <- computeU(X_ord, beta_ord, E)
    psi <- computePsi(X_psi, beta_psi, U, LL)
    theta <- computeTheta(X_theta, beta_theta)

    p <- matrix(.9, maxL, S)
    q <- matrix(.05, maxL, S)
    theta0 <- rep(.05, S)

    mu1 <- 7
    sigma1 <- 3
    mu0 <- 1
    sigma0 <- 1
  }

  U <- computeU(X_ord, beta_ord, E)
  psi <- computePsi(X_psi, beta_psi, U, LL)
}

nburn <- 2000
niter <- 2000

# output
{
  beta_ord_output <- array(NA, dim = c(ncov_ord, d, niter))
  beta_psi_output <- array(NA, dim = c(ncov_psi, S, niter))
  beta_theta_output <- array(NA, dim = c(ncov_theta, S, niter))
  LL_output <- array(NA, dim = c(d, S, niter))
  E_output <- array(NA, dim = c(n, d, niter))
  U_output <- array(NA, dim = c(n, d, niter))
  UL_output <- array(NA, dim = c(n, S, niter))
  z_output <- array(NA, dim = c(n, S, niter))
  p_output <- array(NA, dim = c(maxL, S, niter))
  q_output <- array(NA, dim = c(maxL, S, niter))
  theta0_output <- array(NA, dim = c(S, niter))
}

for (iter in 1:(niter + nburn)) {

  print(iter)

  # sample z
  z <- sample_z_cpp(w, psi, theta, theta0, M, sumM)

  # sample psi

  {
    # k <- z - .5

    # sample Omega
    {
      # U <- (X_ord %*% beta_ord + E)
      # Xbeta_coef <- X_psi %*% beta_psi + U %*% LL
      #
      # Omega <- samplePGvariables(Xbeta_coef)
    }

    # sample beta psi and LL
    {
      # list_betapsiLL <- sample_betapsiLL(beta_psi, LL, Omega, X_psi, k, U,
      #                                    prior_beta_psi, prior_beta_psi_sd)
      # beta_psi <- list_betapsiLL$beta_psi
      # LL <- list_betapsiLL$LL
    }

    # beta_psi <- sample_betapsi(beta_psi, LL, Omega, X_psi, k, U,
    #                            prior_beta_psi, prior_beta_psi_sd)
    # LL <- sample_LL(beta_psi, LL, Omega, X_psi, k, U,
    #                 prior_beta_psi, prior_beta_psi_sd)
    #
    # E <- sample_E(k, Omega, X_psi, beta_psi, X_ord, beta_ord, LL)

    if(ncov_ord > 0){
      # list_betaordE <- sample_betaordE(k, Omega, X_psi, beta_psi, X_ord, LL)
      # beta_ord <- list_betaordE$beta_ord
      # E <- list_betaordE$E
      # beta_ord <- sample_betaord(k, LL, X_ord, X_psi, beta_psi, E, Omega)
    }
  }

  {
    list_betapsiLL <- sample_psivars(z, X_psi, beta_psi, X_ord, beta_ord, E, LL,
                                     prior_beta_psi, prior_beta_psi_sd)
    beta_psi <- list_betapsiLL$beta_psi
    beta_ord <- list_betapsiLL$beta_ord
    LL <- list_betapsiLL$LL
    E <- list_betapsiLL$E
    U <- computeU(X_ord, beta_ord, E)
    psi <- computePsi(X_psi, beta_psi, U, LL)
  }

  # sample w
  w <- sample_w_cpp(logy1, mu0, sigma0, mu1, sigma1, theta, theta0, p, q,
                M, K, sumL, sumM, sumK, maxL, z)

  # sample cimk

  {
    w_all <- w[idx_k,,drop=F]

    c_imk <- ifelse(w_all == 1 & logy1 > 0, 1,
                    ifelse(w_all == 0 & logy1 > 0, 2, 0))

    # c_imk <- sample_cimk(logy1, mu1, sigma1, pi0, sigma0,
    # p, q, idx_k, primerIdx)

    # c_imk <- sample_cimk_cpp(y, logy1, w, mu1, sigma1, pi0, sigma0,
    #                          p, q, idx_k, primerIdx, maxL)
  }

  # sample theta
  beta_theta <- sample_betatheta(w, z, beta_theta, idx_z, X_theta,
                                 b_betatheta, B_betatheta)
  theta <- computeTheta(X_theta, beta_theta)

  # sample pq
  list_pq <- sample_pq(y, w, primerIdx, idx_k, maxL, a_p, b_p, a_q, b_q)
  list_pq <- sample_pq_cpp(c_imk, w, primerIdx, idx_k, maxL, a_p, b_p, a_q, b_q)
  (p <- list_pq$p)
  (q <- list_pq$q)

  # sample theta0
  theta0 <- sample_theta0(z, w, idx_z, a_theta0, b_theta0)

  # sample mu1sigma1
  list_mu1sigma1 <- sample_musigma(sigma1,
                               logy1, c_imk,
                               mu_mu1, sd_mu1,
                               a_sigma1, b_sigma1, tpfp = T)
  mu1 <- list_mu1sigma1$mu1
  sigma1 <- list_mu1sigma1$sigma1

  # sample mu0sigma0
  list_mu1sigma1 <- sample_musigma(sigma0,
                               logy1, c_imk,
                               mu_mu0, sd_mu0,
                               a_sigma0, b_sigma0, tpfp = F)
  mu0 <- list_mu1sigma1$mu1
  sigma0 <- list_mu1sigma1$sigma1

  {

    if(iter > nburn){
      currentIter <- iter - nburn
      beta_psi_output[,,currentIter] <- beta_psi
      beta_ord_output[,,currentIter] <- beta_ord
      beta_theta_output[,,currentIter] <- beta_theta
      LL_output[,,currentIter] <- LL
      U_output[,,currentIter] <- U
      E_output[,,currentIter] <- E
      # UL_output[,,currentIter] <- U %*% LL
      p_output[,,currentIter] <- p
      q_output[,,currentIter] <- q
      theta0_output[,currentIter] <- theta0
      z_output[,,currentIter] <- z
    }

  }

}

# reparam
{

  L_output_reparam <- LL_output
  U_output_reparam <- U_output
  E_output_reparam <- E_output
  beta_ord_output_reparam <- beta_ord_output

  d <- dim(LL_output)[1]

  for (iter in 1:niter) {
    # print(iter)

    if(d == 1){

      L1 <- LL_output[1,1, iter]
      L_output_reparam[1,,iter] <- LL_output[1,,iter] / L1
      U_output_reparam[,1,iter] <- U_output[,1,iter] * L1
      E_output_reparam[,1,iter] <- E_output[,1,iter] * L1
      beta_ord_output_reparam[,1,iter] <- beta_ord_output[,1,iter] * L1

    } else {

      L_current <- LL_output[,,iter]
      # E_current <- E_output[iter,,]
      U_current <- U_output[,,iter]
      beta_ord_current <- beta_ord_output[,,iter]

      qr_decomp <- qr(L_current)
      Q_current <- qr.Q(qr_decomp)
      R_current <- qr.R(qr_decomp)

      Q2 <- Q_current %*% diag(diag(R_current), nrow = d)
      invQ2 <- diag(1 / diag(R_current), nrow = d) %*% t(Q_current)

      betapsiord_new <- beta_ord_current %*% Q2
      # E_new <- E_current %*% Q2
      L_new <- invQ2 %*% L_current
      U_new <- U_current %*% Q2

      L_output_reparam[,,iter] <- L_new
      # E_output_reparam[iter,] <- E_new
      U_output_reparam[,,iter] <- U_new
      beta_ord_output_reparam[,,iter] <- betapsiord_new
    }

  }
}

E_qtl <- apply(E_output[,1,], 1, function(x){
  quantile(x, probs = c(0.025, 0.975))
})

ggplot(data = NULL, aes(x = 1:n,
                        ymin = E_qtl[1,],
                        ymax = E_qtl[2,],
                        y = E_true[,1])) + geom_point() + geom_errorbar()

# beta psi
{
  # beta_ord_output_vec <- apply(beta_ord_output_reparam, 3, as.vector)
  beta_psi_output_vec <- apply(beta_psi_output, 3, as.vector)
  betaord_qtl <- apply(beta_psi_output_vec, 1, function(x){
    quantile(x, probs = c(0.025, 0.975))
  })

  ggplot(data = NULL, aes(x = 1:(ncov_psi * S),
                          ymin = betaord_qtl[1,],
                          ymax = betaord_qtl[2,],
                          y = as.vector(beta_psi_true))) + geom_point() +
    geom_errorbar() +
    geom_hline(aes(yintercept = 0), color = "red")

}

# beta ord
{
  # beta_ord_output_vec <- apply(beta_ord_output_reparam, 3, as.vector)
  beta_ord_output_vec <- apply(beta_ord_output, 3, as.vector)
  betaord_qtl <- apply(beta_ord_output_vec, 1, function(x){
    quantile(x, probs = c(0.025, 0.975))
  })

  ggplot(data = NULL, aes(x = 1:(ncov_ord * d),
                          ymin = betaord_qtl[1,],
                          ymax = betaord_qtl[2,],
                          y = as.vector(beta_ord_true))) +
    geom_point() + geom_errorbar() +
    geom_hline(aes(yintercept = 0))

}

# beta theta
{
  # beta_ord_output_vec <- apply(beta_ord_output_reparam, 3, as.vector)
  beta_theta_output_vec <- apply(beta_theta_output, 3, as.vector)
  betatheta_qtl <- apply(beta_theta_output_vec, 1, function(x){
    quantile(x, probs = c(0.025, 0.975))
  })

  ggplot(data = NULL, aes(x = 1:(ncov_theta * S),
                          ymin = betatheta_qtl[1,],
                          ymax = betatheta_qtl[2,],
                          y = as.vector(beta_theta_true))) +
    geom_point() + geom_errorbar() +
    geom_hline(aes(yintercept = 0))

}

# L
{
  L_output_vec <- apply(LL_output, 3, as.vector)
  L_qtl <- apply(L_output_vec, 1, function(x){
    quantile(x, probs = c(0.025, 0.975))
  })

  L_true

  ggplot(data = NULL, aes(x = 1:(S * d),
                          ymin = L_qtl[1,],
                          ymax = L_qtl[2,],
                          y = as.vector(L_true))) + geom_point() + geom_errorbar()

}

# E
{
  E_output_vec <- apply(E_output, 3, as.vector)
  E_qtl <- apply(E_output_vec, 1, function(x){
    quantile(x, probs = c(0.025, 0.975))
  })

  ggplot(data = NULL, aes(x = 1:(n * d),
                          ymin = E_qtl[1,],
                          ymax = E_qtl[2,],
                          y = as.vector(E_true))) + geom_point() + geom_errorbar()

}

# U
{
  U_output_vec <- apply(U_output, 3, as.vector)
  U_qtl <- apply(U_output_vec, 1, function(x){
    quantile(x, probs = c(0.025, 0.975))
  })

  ggplot(data = NULL, aes(x = 1:(n * d),
                          ymin = U_qtl[1,],
                          ymax = U_qtl[2,],
                          y = as.vector(U_true))) + geom_point() + geom_errorbar()

}

# UL
{
  i <- 1:50
  s <- 1:2
  UL_output_vec <- apply(UL_output[i,s,], 3, as.vector)
  UL_qtl <- apply(UL_output_vec, 1, function(x){
    quantile(x, probs = c(0.025, 0.975))
  })

  UL_true <- U_true %*% L_true
  UL_true <- UL_true[i,s]

  ggplot(data = NULL, aes(x = 1:(length(i) * length(s)),
                          ymin = UL_qtl[1,],
                          ymax = UL_qtl[2,],
                          y = as.vector(UL_true))) +
    geom_point() + geom_errorbar()

}

# p q theta0
{
  # beta_ord_output_vec <- apply(beta_ord_output_reparam, 3, as.vector)
  # p_output_vec <- apply(p_output, c(1,2), as.vector)
  p_qtl <- apply(p_output, c(1,2), function(x){
    quantile(x, probs = c(0.025, 0.975))
  })

  q_qtl <- apply(q_output, c(1,2), function(x){
    quantile(x, probs = c(0.025, 0.975))
  })

  theta0_qtl <- apply(theta0_output, 1, function(x){
    quantile(x, probs = c(0.025, 0.975))
  })

}

# diagnostics
{
  View(
    cbind(logy1[,1], c_imk[,1], w[idx_k,1])
  )
}
