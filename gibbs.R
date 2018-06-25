gibbs <- function(y, G, dirich.prior = NULL, niter = 7500, 
                  n.burn = 2500, n.thin = 10, relabel = TRUE, verbatim = TRUE)
  {
  if (! all(y == 0 | y == 1))
    stop("y must be code 0 or 1") # stop if y is not coded 0,1
  if (niter <= n.burn)   # niter has to be greater than n.burn, error if not
    stop(paste("niter=",niter,",must be greater than n.burn=", n.burn))
  ###
  # loading packages needed to run Gibbs sampling and basic settings
  ###
  require(gtools)  # rdirichlet()
  K <- ncol(y)   # number of items
  N <- nrow(y)  # number of respondents
  G <- G      # number of latent groups
  done.burn <- FALSE # burn.in is not yet done
  
  ###
  # MCMC basic setup, number of iterations and storages for chains
  ###
  Pi <- matrix(NA, nrow = niter, ncol = G) # storage for class membership
  Pjk <- array(NA, dim = c(niter, G, K))   # storage for item resp prob
  dimnames(Pjk) <- list(NULL, paste("G", 1:G, sep = ""), colnames(y))
  Cij <- array(NA, dim = c(niter, N, G)) # storage for discrete classes
  Cij.pr <- array(NA, dim = c(niter, N, G)) # storage for latent class prob
  labelStor <- matrix(NA, nrow = niter, ncol = G) # storage for relabeling
  #
  ## Storage for simulated parameters pjk, C, at iteration t+1
  #
  pjk.t1 <- matrix(NA, nrow = G, ncol = K) # latest p_jk^(t+1) stored here
  dimnames(pjk.t1) <- list(paste("G", 1:G, sep = ""), colnames(y))
  #N*G(people by group) matrix of each person's class membership prob
  Clp.t1 <- matrix(NA, nrow = N, ncol = G)
  dimnames(Clp.t1) <- list(paste("N", 1:N, sep = ""), paste("G", 1:G, sep=""))
  ###
  # Priors
  ###
  if (is.null(dirich.prior))
    dirich.prior <- rep(1,G) # flat Dirichlet by default
  # beta prior, alpha=1 and beta=1 for a flat prior
  alpha <- 1
  beta <-1
  ###
  # starting values of pi and pjk, drawn randomly from dirichlet, beta priors
  ###
  start.pi <- rdirichlet(n=1, alpha = dirich.prior)
  start.item.p <- matrix(NA, nrow = G, ncol = K)
  for (g in 1:G) {
    start.item.p[g,] <- rbeta(K, shape1 = alpha, shape2 = beta)
  }
  ###
  pi.t <- start.pi # membership distr [pi1=0.78, pi2=0.11, pi3=0.11]
  pjk.t <- start.item.p #item response probability pjk on iteration t
  # used later to address the label switch problem
  perm <- gtools::permutations(n=G, r=G) # 24 total permutations when G=4
  trace.num <- numeric(nrow(perm)) # trace of match between t0 and t+1
  ##########
  # Main MCMC simulation
  #########
  iter <- 1 # begins with iteration number 1
  while(iter <= niter)
    { # loop until niter is reached
    #each person's class membership prob, using Eq7
    # [c|y,pi,p], first calculate p and 1-p
    pr.p <- t(pjk.t) #transpose to K by G matrix for apply()
    pr.q <- 1-pr.p
    #step 2: binomial item response probability per Eq 2
    A <- apply(y, MAR =1, FUN = function(yv){pr.p^yv*pr.q^(1-yv)})
    A <- array(A, dim = c(K, G, N))
    A <- aperm(A, c(3,2,1)) #reshape into N*G*K
    eq2 <- apply(A, MARGIN = c(1,2), prod) # multiply across K, keeping N*G
    # step 3: each binomial item resp prob weighted by class distr prob pi[j]
    eq2 <- sweep(eq2, MARGIN = 2, STATS = pi.t, FUN = "*")
    # calculate total probability for each person, per Eq 5
    p.total <-apply(eq2, MARGIN = 1, sum)
    # finally, 'divided-by-total' yields latent class membership prob
    Clp.t1 <- eq2/p.total
    #
    #Clp.t1 gives us the probability of each person's latent class membership,
    # e.g., person 1 has (0.30,0.20,0.15, 0.35) of being in class 1, 2, 3, and 4.
    # so latest class membership can be c=[1,0,0,0] with 15% chance, and c=[0,0,0,1]
    # with 35% chance. Next we use these probes to draw a single specific smaple
    # of c from any of the 4 possibilities above. Each person has one and only one class
    # out of G latent classes.
    Cl.t1 <- apply(Clp.t1,1, function(prob){rmultinom(n = 1, size = 1, prob = prob)})
    Cl.t1 <- t(Cl.t1)
    ##
    # next, update pi (per eq 10) and pjk (per eq 11) using the newly calculated N*G matrix of 
    # discrete latent class membership
    ##
    # sample $\pi_j^{(t+1)}$, percentages of latent classes in the population 
    # Eq 10 shows posterior = data by colSums(C.t)+prior sample sizes
    pi.t1 <- rdirichlet(n = 1, alpha = colSums(Cl.t1) + dirich.prior)
    # sample item response probability, one latent class at a time, sum over N
    for (g in 1:G)
      { #eahc column not guaranteed to add up to 1
      # Eq 11 shows posterior beta(y*c + alpha, (1-y)*c + beta)
      pjk.t1[g,] <- rbeta(K, shape1 = alpha + colSums(Cl.t1[,g]*y), 
                          shape2 = beta + colSums(Cl.t1[,g]*(1-y)))
    }
    # simulated values in current iteraction are added into chain storages
    Pi[iter,] <- pi.t1
    Pjk[iter, , ] <- pjk.t1
    Cij[iter, , ] <- Cl.t1
    Cij.pr[iter, , ] <- Clp.t1
    # 'label switching' problem to match latent classes at end of burn-in 
    if (relabel && done.burn)
    {
      match.tab <- t(Cl.t1) %*% Cl.0 # match between t+1 and t0 latent classes
      for (1 in 1:nrow(perm))
        # across G! permutations, where matches are?
        trace.num[1] <- sum(diag(match.tab[, perm[1,]]))
        
        relabel.ord <- perm[which.max(trace.num), ] # relabel by best match
        labelStor[iter, ] <- relabel.ord
      
    }
    # current simulated values will be used to draw the next iteration
    pi.t <- pi.t1
    pjk.t <- pjk.t1
    # print a message if b.burn iterations done
    if (iter == n.burn){
      done.burn <- TRUE
      cat("\nburn-in completed\n")
      Cl.0 <- Cl.t1 # latent classes immediately after burn-in
    }
    # verbatim can be set by the user to print iteration count every 500th iter
    if (verbatim)
      if(iter == 1) cat("iteration(s) completed:", iter, "")
      if((iter %% 500) < 10^(-7)) {cat(iter, "")}
    iter <- iter + 1  # last thing before repeating is to increment iter by 1
  }
  # end while (iter <= niter)
  cat("\n")
  
  # discard burn-in iterations, thin by n.thin
  ti <- seq(from = n.burn +1, to = niter, by = n.thin)
  Pi.chain <- Pi[ti, ] #pi chain after burn-in and thinning 
  Pjk.chain <- Pjk[ti, ,]  # pjk after burn-in and thinning
  labelStor <- labelStor[ti,]
  Pi.mean <- apply(Pi.chain, 2, mean) # average pi
  Pjk.mean <- apply(Pjk.chain, c(2,3), mean) # average p[jk]
  # put the results together in a list and return() the results
  ans <- list(Pi.mean = Pi.mean, Pjk.mean = Pjk.mean, Pi.chain = Pi.chain,
              Pjk.chain = Pjk.chain, Pjk.mean = Pjk.mean, Pi.chain = Pi.chain,
              Pjk.chain = Pjk.chain, Cpr.chain = Cij.pr, relabelOrd = labelStor)
  return(ans)
}