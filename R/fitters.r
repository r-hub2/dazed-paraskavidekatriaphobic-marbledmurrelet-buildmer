fit.GLMMadaptive <- function (p,formula) {
	fixed <- reformulas::nobars(formula)
	bars <- reformulas::findbars(formula)
	if (is.null(bars)) {
		return(fit.buildmer(p,formula))
	}
	if (length(bars) > 1) {
		# could be a ZCP term
		groups <- sapply(bars,function (x) as.character(list(x[[3]])))
		if ((lug <- length(unique(groups))) > 1) {
			stop(paste0('mixed_model can only handle a single random-effect grouping factor, yet you seem to have specified ',lug))
		}
		terms <- sapply(bars,function (x) as.character(list(x[[2]])))
		# findbars always adds in '0 +'
		if ('1' %in% terms) {
			terms <- sapply(terms,function (x) {
				tab <- tabulate.formula(mkForm(x))
				if (!nrow(tab)) {
					return(x)
				}
				tab$term
			})
		}
		random <- mkForm(paste0(paste0(terms,collapse='+'),'||',groups[[1]]))
	} else {
		random <- mkForm(as.character(bars))
	}
	progress(p,'Fitting via mixed_model: ',fixed,', random=',random)
	patch.GLMMadaptive(p,GLMMadaptive::mixed_model,c(list(fixed=fixed,random=random,data=p$data,family=p$family),p$args))
}

fit.bam <- function (p,formula) {
	re <- re2mgcv(formula,p$data)
	formula <- re$formula
	p$data <- re$data
	if (length(attr(stats::terms(formula),'term.labels')) == 0) {
		# bam is unable to fit intercept-only models
		formula <- add.terms(formula,'intercept')
		nr <- NROW(p$data)
		p$data$intercept <- cbind(rep(1,nr),rep(0,nr))
	}
	if (p$reml) {
		method <- 'fREML'
	} else {
		method <- 'ML'
		p$args$discrete <- FALSE
	}
	progress(p,'Fitting via bam, with ',method,': ',formula)
	patch.lm(p,mgcv::bam,c(list(formula=formula,family=p$family,data=p$data,method=method),p$args))
}

fit.buildmer <- function (p,formula) {
	reml <- p$reml && p$is.gaussian
	if (is.null(reformulas::findbars(formula))) {
		p$args$control <- NULL
		if (reml) {
			# gls has issues with weights
			p$args <- p$args[names(p$args) %in% c('weights','subset','na.action','offset')]
			return(fit.gam(p,formula))
		}
		if (p$is.gaussian) {
			p$args <- p$args[names(p$args) %in% names(formals(stats::lm))]
			progress(p,'Fitting via lm: ',formula)
			patch.lm(p,stats::lm,c(list(formula=formula,data=p$data),p$args))
		} else {
			p$args <- p$args[names(p$args) %in% names(formals(stats::glm))]
			progress(p,'Fitting via glm: ',formula)
			patch.lm(p,stats::glm,c(list(formula=formula,family=p$family,data=p$data),p$args))
		}
	} else {
		if (p$is.gaussian) {
			progress(p,'Fitting via lmer, with ',ifelse(reml,'REML','ML'),': ',formula)
			patch.lmer(p,lme4::lmer,c(list(formula=formula,data=p$data,REML=reml),p$args))
		} else {
			progress(p,'Fitting via glmer, with ',ifelse(reml,'REML','ML'),': ',formula)
			patch.lmer(p,lme4::glmer,c(list(formula=formula,data=p$data,family=p$family),p$args))
		}
	}
}

fit.clmm <- function (p,formula) {
	clm.control <- p$args$clm.control
	clmm.control <- p$args$clmm.control
	p$args$clm.control <- p$args$clmm.control <- NULL
	if (is.null(reformulas::findbars(formula))) {
		if (length(p$args)) {
			p$args <- p$args[names(p$args) %in% names(formals(ordinal::clm))]
			p$args$control <- clm.control
		}
		p$control.name <- p$control.names$clm
		patch.lm(p,ordinal::clm,c(list(formula=formula,data=p$data),p$args))
	} else {
		if (length(p$args)) {
			p$args <- p$args[names(p$args) %in% names(formals(ordinal::clmm))]
			p$args$control <- clmm.control
		}
		p$control.name <- p$control.names$clmm
		patch.lm(p,ordinal::clmm,c(list(formula=formula,data=p$data),p$args))
	}
}

fit.gam <- function (p,formula) {
	re <- re2mgcv(formula,p$data)
	formula <- re$formula
	p$data <- re$data
	if (length(attr(stats::terms(formula),'term.labels')) == 0) {
		# gam is sometimes unable to fit intercept-only models
		formula <- add.terms(formula,'intercept')
		p$data$intercept <- 1
	}
	if (p$quickstart > 0) {
		data <- p$data
		method <- if (p$reml || p$quickstart > 1) 'fREML' else 'ML'
		args <- p$args[names(p$args) %in% names(formals(mgcv::bam))]
		if (method == 'fREML' && p$quickstart > 2 && !'discrete' %in% names(args)) args$discrete <- TRUE
		if (p$quickstart > 3) {
			samfrac <- p$quickstart - floor(p$quickstart)
			if (samfrac == 0) samfrac <- .1
			n <- NROW(data)
			data <- data[sample.int(n,n*samfrac),]
		}
		if (p$quickstart > 4) {
			args$control <- c(p$args$control,list(epsilon=.02))
		}
		progress(p,'Quickstart fit with bam/',method,': ',formula)
		qs <- patch.lm(p,mgcv::bam,c(list(formula=formula,family=p$family,data=data,method=method),args))
		if (!inherits(qs,'try-error')) {
			p$args$in.out <- list(sp=unname(qs$sp),scale=qs$sig2)
			if (startsWith(qs$family$family,'Scaled t')) {
				if (utils::packageVersion('mgcv') < '1.8.32') {
					progress(p,paste0('Starting values: ',paste0(p$args$in.out$sp,collapse=' '),', excluding scaled-t theta values as mgcv version < 1.8.32'))
				} else {
					# set up starting values for theta
					th.notrans <- qs$family$getTheta(FALSE)
					th.trans <- qs$family$getTheta(TRUE)
					min.df <- th.trans[1] - exp(th.notrans[1])
					progress(p,paste0('Starting values: ',paste0(p$args$in.out$sp,collapse=' '),' with theta values ',paste0(th.trans,collapse=' '),' and min.df ',min.df))
					p$family <- mgcv::scat(theta=-th.trans,link=qs$family$link,min.df=min.df)
				}
			} else {
				progress(p,paste0('Starting values: ',paste0(p$args$in.out$sp,collapse=' '),' with scale parameter ',p$args$in.out$scale))
			}
		}
	}
	method <- if (p$reml) 'REML' else 'ML'
	p$args <- p$args[names(p$args) %in% names(formals(mgcv::gam))]
	progress(p,'Fitting via gam, with ',method,': ',formula)
	patch.lm(p,mgcv::gam,c(list(formula=formula,family=p$family,data=p$data,method=method),p$args))
}

fit.gamm <- function (p,formula) {
	fixed <- reformulas::nobars(formula)
	bars <- reformulas::findbars(formula)
	if (is.null(bars)) {
		if (!has.smooth.terms(formula)) {
			p$args <- p$args[names(p$args) %in% names(formals(mgcv::gam))]
			p$quickstart <- 0
			return(fit.gam(p,formula))
		}
		random <- NULL
	} else {
		random <- lapply(bars,function (x) mkForm(as.character(x[2])))
		names(random) <- sapply(bars,function (x) as.character(x[[3]]))
	}
	method <- if (p$reml) 'REML' else 'ML'
	progress(p,'Fitting via gamm, with ',method,': ',fixed,', random=',random)
	m <- patch.lm(p,mgcv::gamm,c(list(formula=fixed,random=random,family=p$family,data=p$data,method=method),p$args))
	if (inherits(m,'try-error') || p$finalize) m else m$lme
}

fit.gamm4 <- function (p,formula) {
	reml <- p$reml && p$is.gaussian
	fixed <- reformulas::nobars(formula)
	bars <- reformulas::findbars(formula)
	random <- if (length(bars)) mkForm(paste('(',sapply(bars,function (x) as.character(list(x))),')',collapse=' + ')) else NULL
	if (is.null(random) && !has.smooth.terms(formula)) return(fit.buildmer(p,formula))
	progress(p,'Fitting via gamm4, with ',ifelse(reml,'REML','ML'),': ',fixed,', random=',random)
	model <- patch.gamm4(p,gamm4::gamm4,c(list(formula=fixed,random=random,family=p$family,data=p$data,REML=reml),p$args))
	if (inherits(model,'try-error') || p$finalize) model else model$mer
}

fit.glmmTMB <- function (p,formula) {
	if (p$reml && is.null(reformulas::findbars(formula))) {
		# work around bug in glmmTMB: REML only works if at least one non-f.e. parameter is specified
		family <- p$family
		if (is.character(family)) family <- get(family)
		if (is.function (family)) family <- family()
		if (family$family %in% c('poisson','binomial')) {
			p$args$control <- NULL
			p$quickstart <- 0
			return(fit.gam(p,formula))
		}
	}
	if ('offset' %in% names(p$args)) {
		# glmmTMB issue #612
		buildmer_offset <- p$args$offset
		p$args$offset <- NULL
		fun <- function (...) glmmTMB::glmmTMB(...,offset=buildmer_offset)
	} else {
		fun <- glmmTMB::glmmTMB
	}
	progress(p,'Fitting via glmmTMB, with ',ifelse(p$reml,'REML','ML'),': ',formula)
	patch.lm(p,fun,c(list(formula=formula,data=p$data,family=p$family,REML=p$reml),p$args))
}

fit.gls <- function (p,formula) {
	method <- if (p$reml) 'REML' else 'ML'
	progress(p,'Fitting via gls, with ',method,': ',formula)
	# gls cannot handle rank-deficient fixed effects --- work around the problem
	dep <- as.character(formula[[2]])
	y <- p$data[[dep]]
	y <- y[!is.na(y)]
	X <- model.matrix(formula,p$data)
	newform <- y ~ 0+X
	newdata <- list(y=y,X=X)
	na <- is.na(coef(stats::lm(newform,newdata)))
	if (ndrop <- sum(na)) {
		progress(p,'gls model is rank-deficient, so dropping ',ndrop,if (ndrop > 1) ' columns/coefficients' else ' column/coefficient','. If this is the final model, the resulting summary may look a bit strange.')
		newdata$X <- newdata$X[,!na]
		return(patch.lm(p,nlme::gls,c(list(newform,data=newdata,method=method),p$args)))
	}
	patch.lm(p,nlme::gls,c(list(formula,data=p$data,method=method),p$args))
}

fit.lme <- function (p,formula) {
	fixed <- reformulas::nobars(formula)
	bars <- reformulas::findbars(formula)
	if ((length(bars) + !is.null(p$args$random)) > 1) stop(paste0('lme can only handle a single random-effect grouping factor, yet you seem to have specified ',length(bars)))
	if (!is.null(bars)) {
		random <- mkForm(as.character(bars))
		# and continue with lme
	} else {
		if (!is.null(p$args$random)) {
			random <- p$args$random
			p$args$random <- NULL
			# and continue with lme
		} else {
			p$args <- p$args[names(p$args) %in% names(c(formals(stats::lm),formals(nlme::gls)))]
			p$args$control <- NULL
			return((if (!is.null(p$args$correlation)) fit.gls else fit.buildmer)(p,formula))
		}
	}
	method <- if (p$reml) 'REML' else 'ML'
	progress(p,'Fitting via lme, with ',method,': ',fixed,', random=',random)
	patch.lm(p,nlme::lme,c(list(fixed,data=p$data,random=random,method=method),p$args))
}

fit.mertree <- function (p,formula) {
	fixed <- reformulas::nobars(formula)
	bars <- reformulas::findbars(formula)
	if (is.null(bars)) {
		ftext <- paste0(as.character(list(fixed)),' | ',p$partitioning,sep='',collapse=' + ')
		f <- stats::as.formula(ftext)
		if (p$is.gaussian) {
			progress(p,'Fitting via lmtree: ',f)
			p$args <- p$args[names(p$args) %in% names(formals(partykit::lmtree))]
			patch.lm(p,partykit::lmtree,c(list(formula=f,data=p$data),p$args))
		} else {
			progress(p,'Fitting via glmtree: ',f)
			p$args <- p$args[names(p$args) %in% names(formals(partykit::glmtree))]
			patch.lm(p,partykit::glmtree,c(list(formula=f,data=p$data,family=p$family),p$args))
		}
	} else {
		random <- paste0('(',sapply(bars,function (x) as.character(list(x))),')',collapse=' + ')
		ftext <- paste0(as.character(list(fixed)),' | ',random,' | ',p$partitioning,collapse=' + ')
		f <- stats::as.formula(ftext)
		if (p$is.gaussian) {
			progress(p,'Fitting via lmertree: ',f)
			patch.mertree(p,glmertree::lmertree,c(list(formula=f,data=p$data),p$args))
		} else {
			progress(p,'Fitting via glmertree: ',f)
			patch.mertree(p,glmertree::glmertree,c(list(formula=f,data=p$data,family=p$family),p$args))
		}
	}
}

fit.multinom <- function (p,formula) {
	progress(p,'Fitting via multinom: ',formula)
	patch.lm(p,nnet::multinom,c(list(formula=formula,data=p$data),p$args))
}

fit.nb <- function (p,formula) {
	# We can't rely on matching formals between glm.nb and glmer.nb, because glmer.nb accepts dots arguments
	in.glm   <- names(p$args %in% names(formals(MASS::glm.nb)))
	in.glmer <- names(p$args %in% names(formals(lme4::glmer.nb)))
	if (is.null(reformulas::findbars(formula))) {
		if (length(p$args)) {
			p$args <- p$args[!(in.glmer & !in.glm)]
		}
		progress(p,'Fitting via glm.nb: ',formula)
		patch.lm(p,MASS::glm.nb,c(list(formula=formula,data=p$data),p$args))
	} else {
		if (length(p$args)) {
			p$args <- p$args[!(in.glm & !in.glmer)]
		}
		progress(p,'Fitting via glmer.nb: ',formula)
		patch.lmer(p,lme4::glmer.nb,c(list(formula=formula,data=p$data),p$args))
	}
}
