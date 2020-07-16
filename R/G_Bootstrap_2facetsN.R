calcVarComp <- function(x) {
  ## Calculate 2-facet variance components via ANOVA from 'long' data Get sample sizes
  np <- nlevels(x$p)
  ni <- nlevels(x$i)
  no <- nlevels(x$o)
  # Compute Sums of Squared Mean Scores (T), Brennan p.69
  xbar <- mean(x$Score)
  Tp <- no * ni * sum(aggregate(x[, 4], list(x$p), mean)[, 2]^2)
  To <- np * ni * sum(aggregate(x[, 4], list(x$o), mean)[, 2]^2)
  Tpo <- ni * sum(aggregate(x[, 4], list(x$p, x$o), mean)[, 3]^2)
  Tio <- np * sum(aggregate(x[, 4], list(x$o, x$i), mean)[, 3]^2)
  Tpio <- sum(x$Score^2)
  Tmu <- np * ni * no * xbar^2
  # Compute Sum of Squares (SS)
  SSp <- Tp - Tmu
  SSo <- To - Tmu
  SSpo <- Tpo - Tp - To + Tmu
  SSio <- Tio - To
  SSpio <- Tpio - Tpo - Tio + To
  # Comupute Mean squares (MS)
  MSp <- SSp/(np - 1)
  MSo <- SSo/(no - 1)
  MSpo <- SSpo/((np - 1) * (no - 1))
  MSio <- SSio/(no * (ni - 1))
  MSpio <- SSpio/(no * (np - 1) * (ni - 1))
  # Compute variance components
  var_pio <- MSpio
  var_po <- (MSpo - MSpio)/ni
  var_io <- (MSio - MSpio)/np
  var_p <- (MSp - ni * var_po - MSpio)/(ni * no)
  var_o <- (MSo - ni * var_po - np * var_io - MSpio)/(np * ni)
  return(c(var_p, var_o, var_po, var_io, var_pio))
}

calcAdjustedVar <- function(.data, .var_p,  .var_o, .var_po, .var_io, .var_pio,
                            type) {
  # Calculate adjusted variance components, given data and original components
  np <- nlevels(.data$p)
  ni <- nlevels(.data$i)
  no <- nlevels(.data$o)
  if (type == "p") {
    adj_var_p <- (np/(np - 1)) * .var_p
    adj_var_o <- .var_o - (1/(np - 1)) * .var_po
    adj_var_po <- (np/(np - 1)) * .var_po
    adj_var_io <- .var_io - (1/(np - 1)) * .var_pio
    adj_var_pio <- (np/(np - 1)) * .var_pio
  }
  if (type == "o") {
    adj_var_p <- .var_p - (1/(no - 1)) * .var_po
    adj_var_o <- (no/(no - 1)) * .var_o
    adj_var_po <- (no/(no - 1)) * .var_po
    adj_var_io <- .var_io
    adj_var_pio <- .var_pio
  }
  if (type == "po") {
    adj_var_p <- (np/(np - 1)) * .var_p - (np/((np - 1) * (no - 1))) * .var_po
    adj_var_o <- (no/(no - 1)) * .var_o - (no/((np - 1) * (no - 1))) * .var_po
    adj_var_po <- (np * no/((np - 1) * (no - 1))) * .var_po
    adj_var_io <- .var_io - (1/(np - 1)) * .var_pio
    adj_var_pio <- np/(np - 1) * .var_pio
  }
  if (type == "io") {
    adj_var_p <- .var_p - (1/(no * (ni - 1))) * .var_pio
    adj_var_o <- (no/(no - 1)) * .var_o - 1/(ni - 1) * .var_io
    adj_var_po <- (no/(no - 1)) * .var_po - 1/(ni - 1) * .var_pio
    adj_var_io <- (ni/(ni - 1)) * .var_io
    adj_var_pio <- (ni/(ni - 1)) * .var_pio
  }
  if (type == "pio") {
    adj_var_p <- (np/(np - 1)) * .var_p - (np/((np - 1) * (no - 1))) * .var_po
    adj_var_o <- (no/(no - 1)) * .var_o - (no/((np - 1) * (no - 1))) * .var_po - 1/(ni - 1) * .var_io + 1/((np - 1) * (ni - 1)) * .var_pio
    adj_var_po <- (np * no/((np - 1) * (no - 1))) * .var_po - np/((np - 1) * (ni - 1)) * .var_pio
    adj_var_io <- (ni/(ni - 1)) * .var_io - (ni/((np - 1) * (ni - 1))) * .var_po
    adj_var_pio <- np * ni/((np - 1) * (ni - 1)) * .var_pio
  }
  return(c(adj_var_p, adj_var_o, adj_var_po, adj_var_io, adj_var_pio))
}

CalcGTheoryCI <- function(Data, B = 1000, type) {
  colnames(Data) <- c("p", "i", "o", "Score")
  Data$p <- factor(Data$p)
  Data$i <- factor(Data$i)
  Data$o <- factor(Data$o)
  np <- nlevels(Data$p)
  ni <- nlevels(Data$i)
  no <- nlevels(Data$o)

  results <- NULL
  registerDoParallel(cores = 4)
  r <- foreach(icount(B), .combine = rbind) %dopar% {
    if (type == "p") {
      ## Get indices for p boot-p
      boot_indices <- sample(1:np, np, replace = T)
      boot_sample <- unlist(sapply(boot_indices, function(x) which(Data$p == x)))
      boot_data <- Data[boot_sample, ]
      # re-index person
      boot_data$p <- Data[order(Data$p), ]$p
      # Compute estimated variance components
      boot_var <- calcVarComp(boot_data)
      # Compute adjusted variance component estimates
      boot_AdjVar <- calcAdjustedVar(boot_data, boot_var[1], boot_var[2], boot_var[3], boot_var[4],
                                     boot_var[5], type = "p")
    }
    if (type == "o") {
      ## Get indices for o boot-o
      boot_indices <- sample(1:no, no, replace = T)
      boot_sample <- unlist(sapply(boot_indices, function(x) which(Data$o == x)))
      boot_data <- Data[boot_sample, ]
      # re-index occasion
      boot_data$o <- Data[order(Data$o), ]$o
      # Compute estimated variance components
      boot_var <- calcVarComp(boot_data)
      # Compute adjusted variance component estimates
      boot_AdjVar <- calcAdjustedVar(boot_data, boot_var[1], boot_var[2], boot_var[3], boot_var[4],
                                     boot_var[5], type = "o")
    }
    if (type == "po") {
      ## Get indices for p and o
      boot_indices_p <- sample(1:np, np, replace = T)
      boot_indices_o <- sample(1:no, no, replace = T)
      # draw person
      boot_sample_p <- unlist(sapply(boot_indices_p, function(x) which(Data$p == x)))
      boot_data_p <- Data[boot_sample_p, ]
      # re-index person
      boot_data_p$p <- Data[order(Data$p), ]$p
      # draw occasion
      boot_sample <- unlist(sapply(boot_indices_o, function(x) which(boot_data_p$o == x)))
      boot_data <- boot_data_p[boot_sample, ]
      # re-index item
      boot_data$o <- Data[order(Data$o), ]$o
      # Compute estimated variance components
      boot_var <- calcVarComp(boot_data)
      # Compute adjusted variance component estimates
      boot_AdjVar <- calcAdjustedVar(boot_data, boot_var[1], boot_var[2], boot_var[3], boot_var[4],
                                     boot_var[5], type = "po")
    }
    if (type == "io") {
      ## Get indices for o
      boot_indices_o <- sample(1:no, no, replace = T)
      # draw occasion
      boot_sample_o <- unlist(sapply(boot_indices_o, function(x) which(Data$o == x)))
      boot_data_o <- Data[boot_sample_o, ]
      # re-index occasion
      boot_data_o$o <- Data[order(Data$o), ]$o
      boot_data <- NULL
      for (t in 1:no){
        boot_indices_i <- sample(1:ni, ni, replace = T)
        # draw item within each occasion
        boot_sample <- unlist(sapply(boot_indices_i, function(x) which(boot_data_o[which(boot_data_o[,'o'] == t), ]$i == x)))
        boot_data <- rbind(boot_data, boot_data_o[which(boot_data_o[,'o'] == t), ][boot_sample, ])
      }
      # re-index item
      boot_data$i <- Data[order(Data$i), ][order(Data$o), ]$i
      # Compute estimated variance components
      boot_var <- calcVarComp(boot_data)
      # Compute adjusted variance component estimates
      boot_AdjVar <- calcAdjustedVar(boot_data, boot_var[1], boot_var[2], boot_var[3], boot_var[4],
                                     boot_var[5], type = "io")
    }
    if (type == "pio") {
      ## Get indices for p
      boot_indices_p <- sample(1:np, np, replace = T)
      # draw person
      boot_sample_p <- unlist(sapply(boot_indices_p, function(x) which(Data$p == x)))
      boot_data_p <- Data[boot_sample_p, ]
      # re-index person
      boot_data_p$p <- Data[order(Data$p), ]$p
      ## Get indices for o
      boot_indices_o <- sample(1:no, no, replace = T)
      # draw occasion
      boot_sample_o <- unlist(sapply(boot_indices_o, function(x) which(boot_data_p$o == x)))
      boot_data_o <- boot_data_p[boot_sample_o, ]
      # re-index occasion
      boot_data_o$o <- Data[order(Data$o), ]$o
      boot_data <- NULL
      for (t in 1:no){
        boot_indices_i <- sample(1:ni, ni, replace = T)
        # draw item within each occasion
        boot_sample <- unlist(sapply(boot_indices_i, function(x) which(boot_data_o[which(boot_data_o[,'o'] == t), ]$i == x)))
        boot_data <- rbind(boot_data, boot_data_o[which(boot_data_o[,'o'] == t), ][boot_sample, ])
      }
      # re-index item
      boot_data$i <- Data[order(Data$i), ][order(Data$o), ]$i
      # Compute estimated variance components
      boot_var <- calcVarComp(boot_data)
      # Compute adjusted variance component estimates
      boot_AdjVar <- calcAdjustedVar(boot_data, boot_var[1], boot_var[2], boot_var[3], boot_var[4],
                                     boot_var[5], type = "pio")
    }

    ## Compute D-study coefficients, B p. 16
    boot_AbsErrVar <- boot_AdjVar[4]/(ni * no) + boot_AdjVar[5]/(ni * no)
    boot_GenCoef <- (boot_AdjVar[1] + boot_AdjVar[3]/no)/(boot_AdjVar[1] + boot_AdjVar[3]/no + boot_AdjVar[5]/(ni * no))
    boot_DepCoef <- (boot_AdjVar[1] + boot_AdjVar[3]/no)/(boot_AdjVar[1] + boot_AdjVar[3]/no + boot_AbsErrVar)

    # Accumulate results such that cols 1-8 are varp,varo,varpo,vario,varpio,
    # AbsErrVar,GenCoef,DepCoef
    result <- c(boot_AdjVar, boot_AbsErrVar, boot_GenCoef, boot_DepCoef)
    result
  }
  results <- rbind(results, r)
  return(results)
}

summaryCI <- function(.result, ConfLevel, rounding) {
  # GstudyEstimates (means and standard error)
  means <- colMeans(.result)
  sds <- apply(.result, 2, sd)
  gstudy_est <- round(rbind(c(means[1:5], sds[1:5])), rounding)
  colnames(gstudy_est) <- c("p Var", "o Var", "po Var",
                            "i:o Var", "ResidVar", "p Var_SE", "o Var_SE",
                            "po Var_SE", "i:o Var_SE", "ResidVar_SE")
  # CIs
  lb_g <- round(apply(.result, 2, quantile, probs = (1 - ConfLevel)/2, names = F), rounding)
  ub_g <- round(apply(.result, 2, quantile, probs = (1 - (1 - ConfLevel)/2), names = F), rounding)
  gstudy_cis_vec <- noquote(paste0("(", lb_g, ", ", ub_g, ")"))
  gstudy_cis <- rbind(gstudy_cis_vec[1:5])
  colnames(gstudy_cis) <- c("p Var", "o Var", "po Var", "i:o Var", "ResidVar")
  # DstudyEstimates (mean and standard error)
  dstudy_est <- round(rbind(c(means[6:8], sds[6:8])), 4)
  colnames(dstudy_est) <- c("AbsErrVar", "GenCoef", "DepCoef", "AbsErrVar_SE", "GenCoef_SE",
                            "DepCoef_SE")
  # DstudyCIs
  dstudy_cis <- rbind(gstudy_cis_vec[6:8])
  colnames(dstudy_cis) <- c("AbsErrVar_SE", "GenCoef_SE", "DepCoef_SE")
  # Return four matrices in a named list
  return(list(Gstudy_Estimates = gstudy_est, Gstudy_Intervals = gstudy_cis, Dstudy_Estimates = dstudy_est,
              Dstudy_Intervals = dstudy_cis))
}
