calcVarComp <- function(x) {
    ## Calculate 2-facet variance components via ANOVA from 'long' data Get sample sizes
    np <- nlevels(x$p)
    ni <- nlevels(x$i)
    no <- nlevels(x$o)
    # Compute Sums of Squared Mean Scores (T), Brennan p.69
    xbar <- mean(x$Score)
    Tp <- no * ni * sum(aggregate(x[, 4], list(x$p), mean)[, 2]^2)
    Ti <- np * no * sum(aggregate(x[, 4], list(x$i), mean)[, 2]^2)
    To <- np * ni * sum(aggregate(x[, 4], list(x$o), mean)[, 2]^2)
    Tpi <- no * sum(aggregate(x[, 4], list(x$p, x$i), mean)[, 3]^2)
    Tpo <- ni * sum(aggregate(x[, 4], list(x$p, x$o), mean)[, 3]^2)
    Tio <- np * sum(aggregate(x[, 4], list(x$o, x$i), mean)[, 3]^2)
    Tpio <- sum(x$Score^2)
    Tmu <- np * ni * no * xbar^2
    # Compute Sum of Squares (SS)
    SSp <- Tp - Tmu
    SSi <- Ti - Tmu
    SSo <- To - Tmu
    SSpi <- Tpi - Tp - Ti + Tmu
    SSpo <- Tpo - Tp - To + Tmu
    SSio <- Tio - Ti - To + Tmu
    SSpio <- Tpio - Tpi - Tpo - Tio + Tp + Ti + To - Tmu
    # Comupute Mean squares (MS)
    MSp <- SSp/(np - 1)
    MSi <- SSi/(ni - 1)
    MSo <- SSo/(no - 1)
    MSpi <- SSpi/((np - 1) * (ni - 1))
    MSpo <- SSpo/((np - 1) * (no - 1))
    MSio <- SSio/((ni - 1) * (no - 1))
    MSpio <- SSpio/((np - 1) * (ni - 1) * (no - 1))
    # Compute variance components
    var_pio <- MSpio
    var_p <- (MSp - MSpi - MSpo + MSpio)/(ni * no)
    var_i <- (MSi - MSpi - MSio + MSpio)/(np * no)
    var_o <- (MSo - MSpo - MSio + MSpio)/(np * ni)
    var_pi <- (MSpi - MSpio)/no
    var_po <- (MSpo - MSpio)/ni
    var_io <- (MSio - MSpio)/np
    return(c(var_p, var_i, var_o, var_pi, var_po, var_io, var_pio))
}

calcAdjustedVar <- function(.data, .var_p, .var_i, .var_o, .var_pi, .var_po, .var_io, .var_pio,
    type) {
    # Calculate adjusted variance components, given data and original components
    np <- nlevels(.data$p)
    ni <- nlevels(.data$i)
    no <- nlevels(.data$o)
    if (type == "p") {
        adj_var_p <- (np/(np - 1)) * .var_p
        adj_var_i <- .var_i - (1/(np - 1)) * .var_pi
        adj_var_o <- .var_o - (1/(np - 1)) * .var_po
        adj_var_pi <- (np/(np - 1)) * .var_pi
        adj_var_po <- (np/(np - 1)) * .var_po
        adj_var_io <- .var_io - (1/(np - 1)) * .var_pio
        adj_var_pio <- (np/(np - 1)) * .var_pio
    }
    if (type == "i") {
        adj_var_p <- .var_p - (1/(ni - 1)) * .var_pi
        adj_var_i <- (ni/(ni - 1)) * .var_i
        adj_var_o <- .var_o - (1/(ni - 1)) * .var_io
        adj_var_pi <- (ni/(ni - 1)) * .var_pi
        adj_var_po <- .var_po - (1/(ni - 1)) * .var_pio
        adj_var_io <- (ni/(ni - 1)) * .var_io
        adj_var_pio <- (ni/(ni - 1)) * .var_pio
    }
    if (type == "o") {
        adj_var_p <- .var_p - (1/(no - 1)) * .var_po
        adj_var_i <- .var_i - (1/(no - 1)) * .var_io
        adj_var_o <- (no/(no - 1)) * .var_o
        adj_var_pi <- .var_pi - (1/(no - 1)) * .var_pio
        adj_var_po <- (no/(no - 1)) * .var_po
        adj_var_io <- (no/(no - 1)) * .var_io
        adj_var_pio <- (no/(no - 1)) * .var_pio
    }
    if (type == "pi") {
        adj_var_p <- (np/(np - 1)) * .var_p - (np/((np - 1) * (ni - 1))) * .var_pi
        adj_var_i <- (ni/(ni - 1)) * .var_i - (ni/((np - 1) * (ni - 1))) * .var_pi
        adj_var_o <- .var_o - (1/(np - 1)) * .var_po - (1/(ni - 1)) * .var_io + (1/((np - 1) *
            (ni - 1))) * .var_pio
        adj_var_pi <- (np * ni/((np - 1) * (ni - 1))) * .var_pi
        adj_var_po <- (np/(np - 1)) * .var_po - (np/((np - 1) * (ni - 1))) * .var_pio
        adj_var_io <- (ni/(ni - 1)) * .var_io - (ni/((np - 1) * (ni - 1))) * .var_pio
        adj_var_pio <- (np * ni/((np - 1) * (ni - 1))) * .var_pio
    }
    if (type == "po") {
        adj_var_p <- (np/(np - 1)) * .var_p - (np/((np - 1) * (no - 1))) * .var_po
        adj_var_i <- .var_i - (1/(np - 1)) * .var_pi - (1/(no - 1)) * .var_io + (1/((np - 1) *
            (no - 1))) * .var_pio
        adj_var_o <- (no/(no - 1)) * .var_o - (no/((np - 1) * (no - 1))) * .var_po
        adj_var_pi <- (np/(np - 1)) * .var_pi - (np/((np - 1) * (no - 1))) * .var_pio
        adj_var_po <- (np * no/((np - 1) * (no - 1))) * .var_po
        adj_var_io <- (no/(no - 1)) * .var_io - (no/((np - 1) * (no - 1))) * .var_pio
        adj_var_pio <- (np * no/((np - 1) * (no - 1))) * .var_pio
    }
    if (type == "io") {
        adj_var_p <- .var_p - (1/(ni - 1)) * .var_pi - (1/(no - 1)) * .var_po + (1/((ni - 1) *
            (no - 1))) * .var_pio
        adj_var_i <- (ni/(ni - 1)) * .var_i - (ni/((ni - 1) * (no - 1))) * .var_io
        adj_var_o <- (no/(no - 1)) * .var_o - (no/((ni - 1) * (no - 1))) * .var_io
        adj_var_pi <- (ni/(ni - 1)) * .var_pi - (ni/((ni - 1) * (no - 1))) * .var_pio
        adj_var_po <- (no/(no - 1)) * .var_po - (no/((ni - 1) * (no - 1))) * .var_pio
        adj_var_io <- (ni * no/((ni - 1) * (no - 1))) * .var_io
        adj_var_pio <- (ni * no/((ni - 1) * (no - 1))) * .var_pio
    }
    if (type == "pio") {
        adj_var_p <- (np/(np - 1)) * .var_p - (np/((np - 1) * (ni - 1))) * .var_pi - (np/((np -
            1) * (no - 1))) * .var_po + (np/((np - 1) * (ni - 1) * (no - 1))) * .var_pio
        adj_var_i <- (ni/(ni - 1)) * .var_i - (ni/((np - 1) * (ni - 1))) * .var_pi - (ni/((ni -
            1) * (no - 1))) * .var_io + (ni/((np - 1) * (ni - 1) * (no - 1))) * .var_pio
        adj_var_o <- (no/(no - 1)) * .var_o - (no/((np - 1) * (no - 1))) * .var_po - (no/((ni -
            1) * (no - 1))) * .var_io + (no/((np - 1) * (ni - 1) * (no - 1))) * .var_pio
        adj_var_pi <- (np * ni/((np - 1) * (ni - 1))) * .var_pi - (np * ni/((np - 1) * (ni - 1) *
            (no - 1))) * .var_pio
        adj_var_po <- (np * no/((np - 1) * (no - 1))) * .var_po - (np * no/((np - 1) * (ni - 1) *
            (no - 1))) * .var_pio
        adj_var_io <- (ni * no/((ni - 1) * (no - 1))) * .var_io - (ni * no/((np - 1) * (ni - 1) *
            (no - 1))) * .var_pio
        adj_var_pio <- (np * ni * no/((np - 1) * (ni - 1) * (no - 1))) * .var_pio
    }
    return(c(adj_var_p, adj_var_i, adj_var_o, adj_var_pi, adj_var_po, adj_var_io, adj_var_pio))
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
                boot_var[5], boot_var[6], boot_var[7], type = "p")
        }
        if (type == "i") {
            ## Get indices for i boot-i
            boot_indices <- sample(1:ni, ni, replace = T)
            boot_sample <- unlist(sapply(boot_indices, function(x) which(Data$i == x)))
            boot_data <- Data[boot_sample, ]
            # re-index item
            boot_data$i <- Data[order(Data$i), ]$i
            # Compute estimated variance components
            boot_var <- calcVarComp(boot_data)
            # Compute adjusted variance component estimates
            boot_AdjVar <- calcAdjustedVar(boot_data, boot_var[1], boot_var[2], boot_var[3], boot_var[4],
                boot_var[5], boot_var[6], boot_var[7], type = "i")
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
                boot_var[5], boot_var[6], boot_var[7], type = "o")
        }
        if (type == "pi") {
            ## Get indices for p and i
            boot_indices_p <- sample(1:np, np, replace = T)
            boot_indices_i <- sample(1:ni, ni, replace = T)
            # draw person
            boot_sample_p <- unlist(sapply(boot_indices_p, function(x) which(Data$p == x)))
            boot_data_p <- Data[boot_sample_p, ]
            # re-index person
            boot_data_p$p <- Data[order(Data$p), ]$p
            # draw item
            boot_sample <- unlist(sapply(boot_indices_i, function(x) which(boot_data_p$i == x)))
            boot_data <- boot_data_p[boot_sample, ]
            # re-index item
            boot_data$i <- Data[order(Data$i), ]$i
            # Compute estimated variance components
            boot_var <- calcVarComp(boot_data)
            # Compute adjusted variance component estimates
            boot_AdjVar <- calcAdjustedVar(boot_data, boot_var[1], boot_var[2], boot_var[3], boot_var[4],
                boot_var[5], boot_var[6], boot_var[7], type = "pi")
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
                boot_var[5], boot_var[6], boot_var[7], type = "po")
        }
        if (type == "io") {
            ## Get indices for i and o
            boot_indices_i <- sample(1:ni, ni, replace = T)
            boot_indices_o <- sample(1:no, no, replace = T)
            # draw item
            boot_sample_i <- unlist(sapply(boot_indices_i, function(x) which(Data$i == x)))
            boot_data_i <- Data[boot_sample_i, ]
            # re-index item
            boot_data_i$i <- Data[order(Data$i), ]$i
            # draw occasion
            boot_sample <- unlist(sapply(boot_indices_o, function(x) which(boot_data_i$o == x)))
            boot_data <- boot_data_i[boot_sample, ]
            # re-index occasion
            boot_data$o <- Data[order(Data$o), ]$o
            # Compute estimated variance components
            boot_var <- calcVarComp(boot_data)
            # Compute adjusted variance component estimates
            boot_AdjVar <- calcAdjustedVar(boot_data, boot_var[1], boot_var[2], boot_var[3], boot_var[4],
                boot_var[5], boot_var[6], boot_var[7], type = "io")
        }
        if (type == "pio") {
            ## Get indices for pio
            boot_indices_p <- sample(1:np, np, replace = T)
            boot_indices_i <- sample(1:ni, ni, replace = T)
            boot_indices_o <- sample(1:no, no, replace = T)
            # draw person
            boot_sample_p <- unlist(sapply(boot_indices_p, function(x) which(Data$p == x)))
            boot_data_p <- Data[boot_sample_p, ]
            # re-index person
            boot_data_p$p <- Data[order(Data$p), ]$p
            # draw item
            boot_sample_pi <- unlist(sapply(boot_indices_i, function(x) which(boot_data_p$i ==
                x)))
            boot_data_pi <- boot_data_p[boot_sample_pi, ]
            # re-index item
            boot_data_pi$i <- Data[order(Data$i), ]$i
            # draw occasion
            boot_sample <- unlist(sapply(boot_indices_o, function(x) which(boot_data_pi$o == x)))
            boot_data <- boot_data_pi[boot_sample, ]
            # re-index occasion
            boot_data$o <- Data[order(Data$o), ]$o
            # Compute estimated variance components
            boot_var <- calcVarComp(boot_data)
            # Compute adjusted variance component estimates
            boot_AdjVar <- calcAdjustedVar(boot_data, boot_var[1], boot_var[2], boot_var[3], boot_var[4],
                boot_var[5], boot_var[6], boot_var[7], type = "pio")
        }

        ## Compute D-study coefficients
        boot_AbsErrVar <- boot_AdjVar[2]/ni + boot_AdjVar[3]/no + boot_AdjVar[4]/ni + boot_AdjVar[5]/no +
            boot_AdjVar[6]/(ni * no) + boot_AdjVar[7]/(ni * no)
        boot_GenCoef <- boot_AdjVar[1]/(boot_AdjVar[1] + boot_AdjVar[4]/ni + boot_AdjVar[5]/no +
            boot_AdjVar[7]/(ni * no))
        boot_DepCoef <- boot_AdjVar[1]/(boot_AdjVar[1] + boot_AbsErrVar)

        # Accumulate results such that cols 1-10 are varp,vari,varo,varpi,varpo,vario,varpio,
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
    gstudy_est <- round(rbind(c(means[1:7], sds[1:7])), rounding)
    colnames(gstudy_est) <- c("p Var", "i Var", "o Var", "pi Var", "po Var",
        "io Var", "ResidVar", "p Var_SE", "i Var_SE", "o Var_SE", "pi Var_SE",
        "po Var_SE", "io Var_SE", "ResidVar_SE")
    # CIs
    lb_g <- round(apply(.result, 2, quantile, probs = (1 - ConfLevel)/2, names = F), rounding)
    ub_g <- round(apply(.result, 2, quantile, probs = (1 - (1 - ConfLevel)/2), names = F), rounding)
    gstudy_cis_vec <- noquote(paste0("(", lb_g, ", ", ub_g, ")"))
    gstudy_cis <- rbind(gstudy_cis_vec[1:7])
    colnames(gstudy_cis) <- c("p Var", "i Var", "o Var", "pi Var", "po Var",
                              "io Var", "ResidVar")
    # DstudyEstimates (mean and standard error)
    dstudy_est <- round(rbind(c(means[8:10], sds[8:10])), 4)
    colnames(dstudy_est) <- c("AbsErrVar", "GenCoef", "DepCoef", "AbsErrVar_SE", "GenCoef_SE",
        "DepCoef_SE")
    # DstudyCIs
    dstudy_cis <- rbind(gstudy_cis_vec[8:10])
    colnames(dstudy_cis) <- c("AbsErrVar_SE", "GenCoef_SE", "DepCoef_SE")
    # Return four matrices in a named list
    return(list(Gstudy_Estimates = gstudy_est, Gstudy_Intervals = gstudy_cis, Dstudy_Estimates = dstudy_est,
        Dstudy_Intervals = dstudy_cis))
}
