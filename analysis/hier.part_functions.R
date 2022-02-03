# Randy Peterson made modifications to these hier.part functions to "hot wire" them in to accepting output from MASS::glm.nb
# 12/01/2020

all.regs <-
function (y, xcan, family = c("gaussian", "binomial", 
    "Gamma", "inverse.gaussian", "poisson", 
    "quasi", "quasibinomial", "quasipoisson", 
    "beta", "ordinal", "glm.nb"), link = c("logit", "probit", 
    "cloglog", "cauchit", "loglog"), gof = c("Rsqu", 
    "RMSPE", "logLik"), print.vars = FALSE, ...) 
{
    if (length(family) > 1) 
        family <- family[1]
    if (length(link) > 1) 
        link <- link[1]
    if (length(gof) > 1) 
        gof <- gof[1]
    if (sum(is.na(xcan)) > 0) {
        missing <- is.na(apply(xcan, 1, FUN = sum))
        xcan <- xcan[!missing, ]
        y <- y[!missing]
        warning(paste(sum(missing), "observations deleted due to missingness in xcan\n"), 
            call. = FALSE)
    }
    if (sum(is.na(y)) > 0) {
        missing <- is.na(y)
        xcan <- xcan[!missing, ]
        y <- y[!missing]
        warning(paste(sum(missing), "observations deleted due to missingness in y\n"), 
            call. = FALSE)
    }
    if (!family %in% c("gaussian", "binomial", "Gamma", 
        "inverse.gaussian", "poisson", "quasi", 
        "quasibinomial", "quasipoisson", "beta", 
        "ordinal","glm.nb")) {
        stop("The 'family' argument must equal one of 'gaussian', 'binomial',\n              'Gamma', 'inverse.gaussian','poisson', 'quasi', 'quasibinomial',\n              'quasipoisson', 'beta', or 'ordinal'", 
            call. = FALSE)
    }
    if (family != "gaussian" & gof == "Rsqu") {
        stop("The 'gof' argument can only equal R-squared\n             if family = 'gaussian'", 
            call. = FALSE)
    }
    if (!is.vector(y) && dim(y)[2] != 1) {
        cat("\ny must be a vector or a single column data frame")
    }
    pcan <- dim(xcan)[2]
    n <- (2^pcan) - 1
    combs <- combos1(pcan)$ragged
    if (gof != "RMSPE" && gof != "logLik" && gof != 
        "Rsqu") {
        stop(paste("\n gof (goodness of fit measure) must equal", 
            "\n \"RMSPE\" (Root-mean-square \"prediction\" error", 
            "\n \"logLik\" (Log-Likelihood) or", "\n \"Rsqu\" (R-squared)\n\n"), 
            call. = FALSE)
    }
    else {
        if (gof == "RMSPE") {
            #neg binom added
            if (family == "glm.nb") {
                gfs <- sqrt(sum((MASS::glm.nb(y ~ 1, ...)$fitted.values - y)^2))
            }
            if (family == "beta") {
                gfs <- sqrt(sum((betareg::betareg(y ~ 1, family = family, 
                  link = link, ...)$fitted.values - y)^2))
            }
            if (family == "ordinal") {
                gfs <- sqrt(sum((stats::glm(y ~ 1, family = family, 
                  ...)$fitted.values - y)^2))
            }
            if (family == "ordinal") {
                gfs <- sqrt(sum((MASS::polr(y ~ 1, family = family, 
                  method = ifelse(is.null(link), "logistic", 
                    gsub("logit", "logistic", link)), 
                  ...)$fitted.values - y)^2))
            }
            if (!family %in% c("beta", "ordinal", "glm.nb")) {
                gfs <- sqrt(sum((stats::glm(y ~ 1, family = family, 
                  ...)$fitted.values - y)^2))
            }
        }
        if (gof == "logLik") {
            #neg binom added
            if (family == "glm.nb") {
                gfs <- as.vector(stats::logLik(MASS::glm.nb(y ~ 1, ...)))
            }
            if (family == "beta") {
                gfs <- as.vector(stats::logLik(betareg::betareg(y ~ 
                  1, family = family, link = link, ...)))
            }
            if (family == "ordinal") {
                gfs <- as.vector(stats::logLik(MASS::polr(y ~ 
                  1, family = family, method = ifelse(is.null(link), 
                  "logistic", gsub("logit", "logistic", 
                    link)), ...)))
            }
            if (!family %in% c("beta", "ordinal", "glm.nb")) {
                gfs <- as.vector(stats::logLik(stats::glm(y ~ 
                  1, family = family, ...)))
            }
        }
        if (gof == "Rsqu") 
            gfs <- 0
    }
    for (i in 1:n) {
        if (i%%500 == 0) 
            cat(i, "regressions calculated:", n - i, "to go...\n")
        current.comb <- as.vector(combs[i, ][combs[i, ] > 0])
        combn <- paste(names(data.frame(xcan)[current.comb]), 
            "", collapse = "")
        if (gof == "RMSPE") 
            new.line <- current.model(y, current.comb, xcan, 
                family = family, gof = "RMSPE")
        if (gof == "logLik") 
            new.line <- current.model(y, current.comb, xcan, 
                family = family, gof = "logLik")
        if (gof == "Rsqu") 
            new.line <- current.model(y, current.comb, xcan, 
                gof = "Rsqu")
        gfs <- c(gfs, new.line)
    }
    if (print.vars) {
        cat("regressions done: formatting results\n")
        var.names <- "Theta"
        for (i in 1:n) {
            current.comb <- as.vector(combs[i, ][combs[i, ] > 
                0])
            combn <- paste(names(data.frame(xcan)[current.comb]), 
                "", collapse = "")
            new.line <- combn
            var.names <- c(var.names, new.line)
        }
        gfs <- data.frame(`variable combination` = var.names, 
            gof = gfs)
    }
    gfs
}

hier.part <-
function (y, xcan, family = c("gaussian", "binomial", 
    "Gamma", "inverse.gaussian", "poisson", 
    "quasi", "quasibinomial", "quasipoisson", 
    "beta", "ordinal", "glm.nb"), link = c("logit", "probit", 
    "cloglog", "cauchit", "loglog"), gof = c("Rsqu", 
    "RMSPE", "logLik"), barplot = TRUE, ...) 
{
    if (length(family) > 1) 
        family <- family[1]
    if (length(link) > 1) 
        link <- link[1]
    if (length(gof) > 1) 
        gof <- gof[1]
    pcan <- dim(xcan)[2]
    if (pcan > 12) 
        stop("Number of variables must be < 13 for current implementation", 
            call. = FALSE)
    else {
        gfs <- all.regs(y, xcan, family = family, gof = gof, 
            link = link, ...)
        hp <- partition(gfs, pcan, var.names = names(data.frame(xcan)))
        if (barplot) {
            ymin <- min(c(0, floor(min(hp$I.perc) * 0.1) * 10))
            ymax <- ceiling(max(hp$I.perc) * 0.1) * 10
            barplot(t(hp$I.perc), col = c(1), ylim = c(ymin, 
                ymax), ylab = "% Independent effects (%I)")
        }
        params <- list(full.model = paste("y ~", paste(names(xcan), 
            collapse = " + ")), family = family, link = ifelse(family %in% 
            c("beta", "ordinal"), link, "default"), 
            gof = gof)
        list(gfs = gfs, IJ = hp$IJ, I.perc = hp$I.perc, params = params)
    }
}

current.model <-
function (y, current.comb, xcan, family = c("gaussian", 
    "binomial", "Gamma", "inverse.gaussian", 
    "poisson", "quasi", "quasibinomial", "quasipoisson", 
    "beta", "ordinal", "glm.nb"), link = c("logit", "probit", 
    "cloglog", "cauchit", "loglog"), gof = c("Rsqu", 
    "RMSPE", "logLik"), ...) 
{
    if (length(family) > 1) 
        family <- family[1]
    if (length(link) > 1) 
        link <- link[1]
    if (length(gof) > 1) 
        gof <- gof[1]
    if (sum(is.na(xcan)) > 0) {
        missing <- is.na(apply(xcan, 1, FUN = sum))
        xcan <- xcan[!missing, ]
        y <- y[!missing]
    }
    if (sum(is.na(y)) > 0) {
        missing <- is.na(y)
        xcan <- xcan[!missing, ]
        y <- y[!missing]
    }
    if (family != "gaussian" & gof == "Rsqu") {
        stop("R-squared is only appropriate if family = 'gaussian'", 
            call. = FALSE)
    }
    comb.data <- data.frame(xcan[, current.comb])
    colnames(comb.data) <- colnames(xcan)[current.comb]
    data <- data.frame(y, comb.data)
    depv <- names(data)[1]
    n.comb <- dim(comb.data)[2]
    xs <- vector("character", n.comb)
    for (i in 1:(n.comb - 1)) xs[i] <- paste(names(comb.data)[i], 
        "+", sep = "")
    xs[n.comb] <- names(comb.data)[n.comb]
    xss <- paste(xs, collapse = " ", sep = "")
    formu <- stats::formula(paste(depv, "~", xss, sep = ""))
    if (gof == "RMSPE") {
        if (family == "glm.nb") {
            gf <- sqrt(sum((MASS::glm.nb(formu, ...)$fitted.values - y)^2))
        }
        if (family == "beta") {
            gf <- sqrt(sum((betareg::betareg(formu, family = family, 
                link = link, ...)$fitted.values - y)^2))
        }
        if (family == "ordinal") {
            gf <- sqrt(sum((MASS::polr(formu, family = family, 
                method = ifelse(is.null(link), "logistic", 
                  gsub("logit", "logistic", link)), 
                ...)$fitted.values - y)^2))
        }
        if (!family %in% c("beta", "ordinal", "glm.nb")) {
            gf <- sqrt(sum((stats::glm(formu, data = data, family = family, 
                ...)$fitted.values - y)^2))
        }
    }
    if (gof == "logLik") {
        if (family == "glm.nb") {
            gf <- as.vector(stats::logLik(MASS::glm.nb(formu, data = data, ...)))
        }
        if (family == "beta") {
            gf <- as.vector(stats::logLik(betareg::betareg(formu, 
                data = data, family = family, link = link, ...)))
        }
        if (family == "ordinal") {
            gf <- as.vector(stats::logLik(MASS::polr(formu, data = data, 
                family = family, method = ifelse(is.null(link), 
                  "logistic", gsub("logit", "logistic", 
                    link)), ...)))
        }
        if (!family %in% c("beta", "ordinal", "glm.nb")) {
            gf <- as.vector(stats::logLik(stats::glm(formu, data = data, 
                family = family, ...)))
        }
    }
    if (gof == "Rsqu") 
        gf <- summary(stats::lm(formu, data = data))$r.squared
    gf
}