####Function to calculate the FVI####
#Dr. Emigdio Marin-Enriquez (emarin@cibnor.mx)
FVI_MC <- function(
    K, KLow, KUp,
    Tmax, TmaxLow, TmaxUp,
    L50, L50Low, L50Up,
    Linf, LinfLow, LinfUp,
    cor_K_Linf = NULL,
    n = 10000,
    seed = NULL
) {
  
  #-----------------------------#
  # 1. Input validation
  #-----------------------------#
  
  if (n <= 0 || n != as.integer(n)) {
    stop("n must be a positive integer.")
  }
  
  if (KLow > K || K > KUp) {
    stop("K must be between KLow and KUp.")
  }
  
  if (TmaxLow > Tmax || Tmax > TmaxUp) {
    stop("Tmax must be between TmaxLow and TmaxUp.")
  }
  
  if (L50Low > L50 || L50 > L50Up) {
    stop("L50 must be between L50Low and L50Up.")
  }
  
  if (LinfLow > Linf || Linf > LinfUp) {
    stop("Linf must be between LinfLow and LinfUp.")
  }
  
  if (any(c(
    K, KLow, KUp,
    Tmax, TmaxLow, TmaxUp,
    L50, L50Low, L50Up,
    Linf, LinfLow, LinfUp
  ) <= 0)) {
    stop("All parameters and uncertainty bounds must be positive.")
  }
  
  # Check correlation if provided
  if (!is.null(cor_K_Linf)) {
    
    if (length(cor_K_Linf) != 1 ||
        !is.finite(cor_K_Linf) ||
        cor_K_Linf < -1 ||
        cor_K_Linf > 1) {
      stop("cor_K_Linf must be a single value between -1 and 1.")
    }
    
    if (cor_K_Linf > 0) {
      warning(
        "Positive correlation between K and Linf is unusual for ",
        "von Bertalanffy growth parameters. Please verify the ",
        "original source and parameterization."
      )
    }
  }
  
  # Set random seed
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  
  #-----------------------------#
  # 2. Triangular random sampler
  #-----------------------------#
  
  rtriangle <- function(n, min, mode, max) {
    
    if (max == min) {
      return(rep(min, n))
    }
    
    u <- runif(n)
    Fc <- (mode - min) / (max - min)
    
    ifelse(
      u < Fc,
      min + sqrt(
        u * (max - min) * (mode - min)
      ),
      max - sqrt(
        (1 - u) * (max - min) * (max - mode)
      )
    )
  }
  
  
  #-----------------------------#
  # 3. Simulate K and Linf
  #-----------------------------#
  
  if (is.null(cor_K_Linf)) {
    
    # Independent sampling
    K_sim <- rtriangle(
      n, KLow, K, KUp
    )
    
    Linf_sim <- rtriangle(
      n, LinfLow, Linf, LinfUp
    )
    
  } else {
    
    # Correlation matrix
    R <- matrix(
      c(
        1, cor_K_Linf,
        cor_K_Linf, 1
      ),
      nrow = 2,
      byrow = TRUE
    )
    
    # Cholesky decomposition
    L <- chol(R)
    
    # Independent standard normal variables
    Z <- matrix(
      rnorm(2 * n),
      nrow = n,
      ncol = 2
    )
    
    # Correlated normal variables
    Z_cor <- Z %*% L
    
    # Convert to uniform variables
    U <- pnorm(Z_cor)
    
    # Triangular quantile function
    triangular_quantile <- function(
    u, min, mode, max
    ) {
      
      if (max == min) {
        return(rep(min, length(u)))
      }
      
      Fc <- (mode - min) / (max - min)
      
      ifelse(
        u < Fc,
        min + sqrt(
          u * (max - min) * (mode - min)
        ),
        max - sqrt(
          (1 - u) * (max - min) * (max - mode)
        )
      )
    }
    
    # Transform correlated uniforms
    # into triangular distributions
    K_sim <- triangular_quantile(
      U[, 1],
      KLow,
      K,
      KUp
    )
    
    Linf_sim <- triangular_quantile(
      U[, 2],
      LinfLow,
      Linf,
      LinfUp
    )
  }
  
  
  #-----------------------------#
  # 4. Simulate Tmax and L50
  #-----------------------------#
  
  Tmax_sim <- rtriangle(
    n,
    TmaxLow,
    Tmax,
    TmaxUp
  )
  
  L50_sim <- rtriangle(
    n,
    L50Low,
    L50,
    L50Up
  )
  
  
  #-----------------------------#
  # 5. Calculate FVI
  #-----------------------------#
  
  FVI_sim <- (
    Tmax_sim / K_sim
  ) * (
    L50_sim / Linf_sim
  )
  
  
  #-----------------------------#
  # 6. Summary statistics
  #-----------------------------#
  
  summary <- c(
    median = median(FVI_sim),
    lower_95 = unname(
      quantile(FVI_sim, 0.025)
    ),
    upper_95 = unname(
      quantile(FVI_sim, 0.975)
    )
  )
  
  
  #-----------------------------#
  # 7. Return results
  #-----------------------------#
  
  list(
    FVI = summary,
    simulations = FVI_sim,
    K_sim = K_sim,
    Linf_sim = Linf_sim,
    Tmax_sim = Tmax_sim,
    L50_sim = L50_sim,
    correlation = cor(K_sim, Linf_sim)
  )
}