
library(rstan)
library(coda)
library(tidyverse)
rm(list=ls())


assign_labels <- function(E){
	
	if(!is_tibble(E) | !is.data.frame(E)){
		stop("E must be a tibble or data frame")
	}
	
	# get the species names and num sp
	sp_names <- sort(names(E))
	nspp <- length(sp_names)
	
	# sort the cols alphabetically
	E <- E %>% select(sp_names)
	
 	# create label names	
	labels <- expand.grid(replicate(nspp, 0:1, simplify = FALSE)) %>% filter(rowSums(.)>0) %>% as_tibble %>% setNames(sp_names) %>% 
		mutate(comm_id = as.numeric((.>0)%*%(2^(0:(ncol(.)-1))))) %>% gather(species,abundance,-comm_id) %>% filter(abundance>0) %>% 
		select(-abundance) %>% group_by(comm_id) %>% summarize(label=paste(species,collapse="-"))
	
	# assigne labels and label names
 	E <- E %>% filter(rowSums(.)>0) %>% mutate(tot_abund=rowSums(.), comm_id = as.numeric((.>0)%*%(2^(0:(ncol(.)-1))))) %>% 
			arrange(comm_id, tot_abund) %>% select(-tot_abund) %>% left_join(labels,"comm_id")
 
 	# rep numbers
 	rep_ids <- E %>% group_by(comm_id) %>% tally() %>% ungroup %>% group_by(comm_id) %>% expand(count = seq(1:n)) %>% ungroup %>% arrange(comm_id, count) %>% select(-comm_id) %>% rename(rep_num = count)

 	E <- E %>% bind_cols(rep_ids) %>% select(-comm_id) %>% rename(comm_id=label)
	return(E)
}


fit_stan <- function(dt, stan_file, B_upper, B_lower, omit_comms = "none", 
					 chains, cores, iter, thin, warmup, control=NULL,  seed=122457){
	omit_comms <- c("as","as-fa")
	
	
	options(mc.cores = parallel::detectCores())	
	Sys.setenv(USE_CXX14 = 1)	
	rstan_options(auto_write = TRUE)
		
	# make sure the comm labels are in the right order
	omit_comms <- as.character(sapply(omit_comms, function(x) paste(sort(strsplit(x,"-")[[1]]),collapse="-")))
	
	
	E <- dt %>% filter(!comm_id%in%omit_comms)
	
	Emean <- E %>% gather(species, abundance, -comm_id, -rep_num) %>% filter(abundance>0) %>% group_by(species, comm_id) %>% 
		summarize(abundance = median(abundance)) %>% spread(species, abundance, fill=0) %>% arrange(desc(comm_id)) %>% select(-comm_id) %>% as.matrix()
	
	E <- E %>% select(-comm_id,-rep_num) %>% as.matrix()
	
	nspp <- ncol(E)
	sp_names <- colnames(E)
	
	# get the naive regression estimate of B
	B_lm <- matrix(NA,nspp,nspp)
	for(i in 1:nspp){
		myE <- Emean[Emean[,i]>0,]
		B_lm[i,] <- as.numeric(lm(rep(-1,nrow(myE))~-1+myE)$coefficients)
	}
	
	if(any(is.na(B_lm))){
		stop("Not enough endpoints to fit the model")
	}
	
	
	if(missing(B_upper)){
		B_upper <- matrix(1,nspp,nspp)
	}
	if(missing(B_lower)){
		B_lower <- matrix(-1,nspp,nspp)
		diag(B_upper) <- 0
	}

	
	B_init <- -diag(1/apply(E %>% as.matrix(), 2, function(x) quantile(x[x>0],0.9)))
	
	maxB <- max(abs(B_init))		
	
	# build stan data list
	stan_data <- list(  N = ncol(E),
						M = nrow(E),
						E = E,
						B_upper = B_upper,
						B_lower = B_lower,
						y = rep(-1,nspp),
						maxB = maxB)
	
	init_list <- replicate(chains,list(sigmax = rep(0.25,nspp), B = B_init), simplify = FALSE)	
	
	stan_fit <- stan(stan_file,
								  data=stan_data, cores = cores, iter=iter, thin = thin, warmup=warmup, 
								  chains = chains, seed=122457, init=init_list,
								  control = control)

	return(list(E_orig = dt, E = E, stan_fit=stan_fit, omit_comms = omit_comms, B_upper = B_upper, B_lower = B_lower, B_init = B_init, seed = seed))
}
			
plot_stan_results <- function(stan_results){
	
	
	with(stan_results, {
		sp_names <- colnames(E)
		nspp <- length(sp_names)
			
		for(i in 1:length(stan_fit@model_pars)){
			if(stan_fit@model_pars[i]!="Blim"){
				post_plot<-As.mcmc.list(stan_fit,pars=stan_fit@model_pars[i]) # this samples from the posterior
				if(ncol(post_plot[[1]])>1){
						show(plot(post_plot)) # this shows the convergence
				}
				else{
					show(plot(post_plot,main=stan_fit@model_pars[i]))
				}
			}
		}
	
		# wrap the chains into a matrix
		stan_fit <- as.matrix(stan_fit)
		
		
		# grab the individual components and rename
		B_stan <- stan_fit[,paste0("B[",t(outer(1:nspp,1:nspp,paste,sep=",")),"]")] %>% data.frame() %>% setNames(as.character(t(outer(sp_names,sp_names,paste,sep="_")))) %>% as.matrix()
		sigma_stan <- stan_fit[,grepl("sigma",colnames(stan_fit)) | grepl("lp",colnames(stan_fit))]
		# get the predicted and observed eq abundances, using the median B
		pred_eq <- NULL
		my_B <- t(matrix(apply(B_stan,2,median),nspp,nspp))
		for(i in 1:nrow(E)){
			subs <- (1:nspp)[E[i,]>0]
			pred_x <-  -rowSums(solve(my_B[subs,subs]))
			if(any(pred_x<0)){
				pred_x <- rep(NA,length(pred_x))
			}
			pred_eq <- rbind(pred_eq, data_frame(Predicted=pred_x,Observed=E[i,subs],species=sp_names[subs],comm_id=as.numeric(matrix((E[i,]>0),nrow=1)%*%(2^(0:(nspp-1)))), rep=i))
		}
		# plot
		show(ggplot(pred_eq,aes(x=Observed, y=Predicted,color=species))+geom_point()+scale_y_log10()+scale_x_log10()+geom_abline(intercept=0, slope=1)+ylab("Predicted equilibrium abundances") + xlab("Observed abundances"))
		# plot the histograms 
		show(ggplot(B_stan %>% data.frame() %>% gather("B_ij", "value"), aes(x = value)) + geom_histogram() + geom_vline(xintercept=0, col="darkred", linetype=2)+facet_wrap(~B_ij, scales = "free") + 
					 	theme(axis.text.x = element_text(size=5))+theme_bw())
	})
	
}


dt <- read_csv("~/Git/endpoints/data/Kuebbing_plants/cleaned_data/invasives.csv")

stan_file <- "~/Git/endpoints/github/Stan_file.stan"


dt <- assign_labels(dt)

stan_results <- fit_stan(dt, stan_file = stan_file, omit_comms = "none", chains = 2, cores = 2, iter = 1000, warmup = 500, thin = 20, seed=10)
	
plot_stan_results(stan_results)


bootstrap_results <- function(stan_results){
	
	
	
	
	