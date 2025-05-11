backward <- function (p) {
	fit.references.parallel <- function (p) {
		progress(p,'Fitting ML and REML reference models')
		plist <- list(within(p, reml <- TRUE),within(p, reml <- FALSE))
		repeat {
			res <- p$parply(plist,p$fit,p$formula)
			conv <- lapply(res,converged,p$singular.ok,p$grad.tol,p$hess.tol)
			if (all(unlist(conv))) {
				p$cur.reml <- res[[1]]
				p$cur.ml <- res[[2]]
				return(p)
			}
			p <- reduce.model(p,conv)
			if (!any(!is.na(p$tab$block))) return(p)
		}
	}

	if (is.null(p$tab)) {
		p$tab <- tabulate.formula(p$formula)
	}
	if (!p$crit.name %in% colnames(p$tab)) {
		if (nrow(p$tab)) {
			p$tab$Iteration <- 0
			p$tab[,p$crit.name] <- rep(NA,nrow(p$tab))
		}
	}

	iter <- 0
	check.random <- function (fun,tab) {
		fun(sapply(unique(tab$block),function (b) {
			i <- which(tab$block == b)
			fun(!is.na(tab[i,'grouping']))
		}))
	}
	repeat {
		need.ml   <- !p$force.reml || !check.random(all,p$tab)
		need.reml <-  p$force.reml ||  check.random(any,p$tab)
		if (need.ml && need.reml && is.null(p$cur.ml) && is.null(p$cur.reml)) {
			p <- fit.references.parallel(p)
		} else {
			if (need.ml && is.null(p$cur.ml)) {
				progress(p,'Fitting ML reference model')
				p$reml <- FALSE
				p$cur.ml <- p$fit(p,p$formula)
				conv <- converged(p$cur.ml,p$singular.ok,p$grad.tol,p$hess.tol)
				if (!conv) {
					p <- reduce.model(p,conv)
					if (!any(!is.na(p$tab$block))) {
						return(p)
					}
					p <- fit.references.parallel(p)
				}
			}
			if (need.reml && is.null(p$cur.reml)) {
				progress(p,'Fitting REML reference model')
				p$reml <- TRUE
				p$cur.reml <- p$fit(p,p$formula)
				conv <- converged(p$cur.reml,p$singular.ok,p$grad.tol,p$hess.tol)
				if (!conv) {
					p <- reduce.model(p,conv)
					if (!any(!is.na(p$tab$block))) {
						return(p)
					}
					p <- fit.references.parallel(p)
				}
			}
		}

		if (!nrow(p$tab)) {
			progress(p,"There's nothing left!")
			return(p)
		}
		progress(p,'Testing terms')
		results <- p$parply(unique(p$tab$block),function (b) {
			if (is.na(b)) {
				# cannot remove term because of 'include'
				return(list(val=rep(NA,sum(is.na(p$tab$block)))))
			}
			i <- which(p$tab$block == b)
			if (!can.remove(p$tab,i) || any(paste(p$tab[i,'term'],p$tab[i,]$grouping) %in% paste(p$include$term,p$include$grouping))) {
				# cannot remove term due to marginality
				return(list(val=rep(NA,length(i))))
			}
			if (p$force.reml) {
				p$reml <- TRUE
				m.cur <- p$cur.reml
			} else {
				p$reml <- p$can.use.reml && all(!is.na(p$tab[i,]$grouping))
				m.cur <- if (need.reml && p$reml) p$cur.reml else p$cur.ml
			}
			f.alt <- build.formula(p$dep,p$tab[-i,],p$env)
			m.alt <- p$fit(p,f.alt)
			val <- if (converged(m.alt,p$singular.ok,p$grad.tol,p$hess.tol)) p$crit(p,m.alt,m.cur) else NaN
			val <- rep(val,length(i))
			list(val=val,model=m.alt)
		})
		results <- unlist(sapply(results,`[[`,1))
		p$tab[,p$crit.name] <- results
		p$tab$Iteration <- iter <- iter+1
		p$results <- rbind(p$results,p$tab)
		progrep <- p$tab
		progrep$index <- progrep$code <- progrep$ok <- NULL
		if (p$crit.name %in% c('LRT','LRT2')) progrep[,p$crit.name] <- exp(results)
		if (p$crit.name %in% c('deviance','devexp')) progrep[,p$crit.name] <- -progrep[,p$crit.name]
		if (!p$quiet) {
			print(progrep)
		}
		remove <- p$elim(results)
		remove <- which(!is.na(remove) & !is.nan(remove) & remove)
		if (length(remove) == 0) {
			progress(p,'All terms are significant')
			p$model <- if (need.reml) p$cur.reml else p$cur.ml
			return(p)
		}

		# Remove the worst offender(s) and continue
		remove <- remove[p$tab[remove,p$crit.name] == max(p$tab[remove,p$crit.name])]
		p$tab <- p$tab[-remove,]
		p$formula <- build.formula(p$dep,p$tab,p$env)
		p$cur.ml <- p$cur.reml <- NULL
		if (length(results) == 1) {
			# Recycle the current model as the next reference model
			p[[if (p$reml) 'cur.reml' else 'cur.ml']] <- results[[remove]]$model
		}
		progress(p,'Updating formula: ',p$formula)
	}
}

can.remove <- function (tab,i) {
	unravel2 <- function (x) unravel(stats::as.formula(paste0('~',x))[[2]])
	t <- tab[i,'term']
	g <- tab[i,'grouping']
	fx <- which(is.na(tab$g))
	tfx <- tab[intersect(i,fx),'term']

	if ('1' %in% t) {
		# If fixed intercept: never remove it
		if (any(is.na(g))) return(FALSE)
		# If random intercept: do not remove if there are subsequent terms
		for (x in g) if (x %in% tab[-c(fx,i),'grouping']) return(FALSE)
	}

	# Do not remove fixed effects that have corresponding random effects
	if (any(tfx %in% tab$term[-fx])) return(FALSE)

	for (x in g) {
		# Do not remove effects participating in interactions
		scope <- if (is.na(x)) fx else which(tab$grouping == x)
		scope <- scope[!scope %in% i]
		for (t in tab[i,'term']) {
			t <- unravel2(t)
			if (any(sapply(tab[scope,'term'],function (x) all(t %in% unravel2(x))))) return(FALSE)
		}
	}

	TRUE
}

forward <- function (p) {
	if (p$ordered != p$crit.name) {
		p <- order(p)
	} else if (p$ordered == 'custom') {
		warning("Assuming, but not checking, that direction='order' had used the same elimination criterion as requested for forward stepwise. If this is not the case, add an explicit 'order' step before the 'forward' step using the desired criterion.")
	}
	progrep <- p$tab
	progrep$index <- progrep$code <- progrep$ok <- NULL
	if (p$crit.name %in% c('LRT','LRT2')) progrep$score <- exp(progrep$score)
	if (p$crit.name %in% c('deviance','devexp')) progrep[,p$crit.name] <- -progrep[,p$crit.name]
	if (!p$quiet) {
		print(progrep)
	}
	remove <- p$elim(p$tab$score)
	# Retain all terms up to the last significant one, even if they were not significant themselves
	# This happens if they hade a smallest crit in the order step, but would still be subject to elimination by the elimination function
	keep <- which(!remove)
	remove[1:length(keep)] <- FALSE
	remove.ok <- sapply(1:nrow(p$tab),function (i) {
		if (is.na(p$tab[i,]$block)) return(FALSE)
		if (!can.remove(p$tab,i)) return(FALSE)
		TRUE
	})
	p$tab[,p$crit.name] <- p$tab$score
	p$results <- p$tab
	p$tab <- p$tab[!(remove & remove.ok),]
	p$formula <- build.formula(p$dep,p$tab,p$env)
	p$reml <- p$can.use.reml
	p$model <- p$fit(p,p$formula)
	p
}

order <- function (p) {
	reorder <- function (p,tab) {
		# Test for marginality
		can.eval <- function (tab) {
			my.ddply <- function (data,split,fun) {
				vec <- data[[split]]
				res <- lapply(unique(vec),function (x) {
					take <- if (is.na(x)) is.na(vec) else !is.na(vec) & vec == x
					fun(data[take,])
				})
				Reduce(rbind,res)
			}

			# 0. Initialize
			tab$ok <- TRUE
			# 1. If there are random effects, evaluate them as a group
			mine <- is.na(tab$grouping)
			my <- tab[mine,]
			tab[!mine,] <- my.ddply(tab[!mine,],'grouping',function (my) {
				g <- my$grouping
				my$grouping <- NA
				my <- can.eval(my)
				my$grouping <- g
				my
			})

			if (nrow(my)) {
				# 2. The intercept should always come first
				if (any(my$term == '1')) {
					my$ok <- my$term == '1'
					return(my)
				}

				# 3. Evaluate marginality. We cannot take the terms already in the formula into account, because that will break things like nesting.
				# Thus, we have to define marginality as ok if there is no lower-order term whose components are a proper subset of the current term.
				if (length(my[my$ok,'term']) > 1) {
					all.components <- lapply(my[my$ok,'term'],function (x) {
						x <- stats::as.formula(paste0('~',x))[[2]]
						if (is.smooth.term(x)) unpack.smooth.terms(x) else unravel(x)
					})
					check <- function (i) {
						test <- all.components[[i]]
						for (x in all.components[-i]) { #walk all other terms' components
							if (any(x == '1')) return(FALSE) #intercept should always come first
							if (all(x %in% test)) return(FALSE)
						}
						TRUE
					}
					my[my$ok,'ok'] <- sapply(1:length(all.components),check)
				}
				tab[mine,] <- my
			}

			# 4. If any term belonging to a single block could not be selected, disqualify the whole block
			tab <- my.ddply(tab,'block',function (x) within(x,{ if (!all(ok)) ok <- FALSE }))

			tab
		}

		p$ordered <- p$crit.name
		repeat {
			have <- p$tab
			cur <- p$fit(p,build.formula(p$dep,have,p$env))
			conv <- converged(cur,p$singular.ok,p$grad.tol,p$hess.tol)
			if (conv) break
			p <- reduce.model(p,conv)
			if (!any(!is.na(p$tab$block))) return(p)
		}
		repeat {
			check <- tab[!tab$code %in% have$code,]
			if (!nrow(check)) {
				p$tab <- have
				p$model <- cur
				return(p)
			}
			check <- can.eval(check)
			check <- check[check$ok,]
			if (!nrow(check)) {
				progress(p,'Could not proceed ordering terms')
				p$tab <- have
				p$model <- cur
				return(p)
			}
			progress(p,paste0('Currently evaluating ',p$crit.name,' for: ',paste0(ifelse(is.na(check$grouping),check$term,paste(check$term,'|',check$grouping)),collapse=', ')))
			mods <- p$parply(unique(check$block),function (b,check,have,p) {
				check <- check[check$block == b,]
				tab <- rbind(have[,c('index','grouping','term')],check[,c('index','grouping','term')])
				form <- build.formula(p$dep,tab,p$env)
				mod <- list(p$fit(p,form))
				rep(mod,nrow(check))
			},check,have,p)
			mods <- unlist(mods,recursive=FALSE)
			check$score <- sapply(mods,function (mod) if (converged(mod,p$singular.ok,p$grad.tol,p$hess.tol)) p$crit(p,cur,mod) else NaN)
			ok <- Filter(function (x) !is.na(x) & !is.nan(x),check$score)
			if (!length(ok)) {
				statuses <- sapply(mods,function (mod) attr(converged(mod,p$singular.ok,p$grad.tol,p$hess.tol),'reason'))
				statuses <- statuses[!is.na(statuses)]
				progress(p,'Ending the ordering procedure due to having reached the maximal feasible model - all higher models failed to converge. The types of convergence failure are:\n',paste(unique(statuses),collapse='\n    '))
				p$tab <- have
				p$model <- cur
				return(p)
			}
			winner <- which(check$score == min(ok))
			check <- check[winner,]
			have <- rbind(have,check)
			if (length(unique(check[winner,'block'] == 1))) cur <- mods[[winner[1]]] else {
				# In principle, there should be only one winner. If there are multiple candidates which happen to add _exactly_ the same amount of information to the model, this is
				# suspicious. Probably the reason is that this is an overfitted model and none of the candidate terms add any new information. The solution is to add both terms, but this
				# needs an extra fit to obtain the new 'current' model.
				form <- build.formula(p$dep,have,p$env)
				cur <- p$fit(p,form)
				conv <- conv(cur,p$singular.ok,p$grad.tol,p$hess.tol)
				if (!conv) {
					progress(p,'The reference model for the next step failed to converge - giving up ordering attempt.\nThe failure was: ',attr(conv,'reason'))
					return(p)
				}
			}
			progress(p,'Updating formula: ',build.formula(p$dep,have,p$env))
		}
	}

	progress(p,'Determining predictor order')
	p$tab$ok <- p$tab$score <- NULL
	tab <- p$tab
	if (!is.null(p$include)) {
		p$tab <- transform(p$include,block=NA,ok=TRUE,score=NA)
		# Remove possible duplicate terms (predictors specified both in 'formula' and in 'include') from consideration in reordering
		fxd.tab <- is.na(tab$grouping)
		fxd.inc <- is.na(p$include$grouping)
		overlap <- NULL
		if (any(fxd.tab) && any(fxd.inc)) {
			overlap <- which(fxd.tab & tab$term %in% p$include$term[fxd.inc])
		}
		if (any(!fxd.tab) && any(!fxd.inc)) {
			for (g in unique(tab$grouping[!fxd.tab])) {
				overlap <- c(overlap,which(tab$grouping == g & tab$term %in% p$include$term[p$include$grouping == g]))
			}
		}
		if (length(overlap)) {
			tab <- tab[-overlap,]
		}
	} else {
		p$tab <- cbind(tab[0,],ok=logical(),score=numeric())
	}
	fxd <- is.na(tab$grouping)
	if ('1' %in% tab[fxd,'term']) { #always keep the intercept
		where <- tab$block == tab[fxd & tab$term == '1',]$block
		p$tab <- rbind(p$tab,transform(tab[where,],ok=TRUE,score=NA))
		tab <- tab[!where,]
		fxd <- is.na(tab$grouping)
	}

	p$reml <- p$force.reml
	if (any(fxd)) {
		p <- reorder(p,tab[fxd,])
	}
	if (any(!fxd)) {
		p$reml <- p$can.use.reml
		p <- reorder(p,tab[!fxd,])
	}
	p$formula <- build.formula(p$dep,p$tab,p$env)
	p
}

reduce.model <- function (p,conv) {
	if (length(conv) == 1) {
		progress(p,'Convergence failure. Reducing terms and retrying...\nThe failure was: ',attr(conv,'reason'))
	} else {
		statuses <- sapply(conv[!unlist(conv)],function (x) attr(x,'reason'))
		progress(p,'Convergence failure. Reducing terms and retrying...\nThe failures were: ',paste(unique(statuses),collapse='\n    '))
	}
	cands <- p$tab$block[!is.na(p$tab$block)]
	if (length(unique(cands)) < 2) {
		stop('No terms left for reduction, giving up')
	}
	p$tab <- p$tab[is.na(p$tab$block) | p$tab$block != cands[length(cands)],]
	p$formula <- build.formula(p$dep,p$tab,p$env)
	p
}
