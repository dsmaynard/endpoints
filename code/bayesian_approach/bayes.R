
library(rstan)
library(coda)
library(tidyverse)


# assign a label to each community based on the presence/absence of each species
assign_labels <- function(E){
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
	# clean up
 	E <- E %>% bind_cols(rep_ids) %>% select(-comm_id) %>% rename(comm_id=label)
	return(E)
}

# estimate B using a Bayesian appraoch via Stan
fit_stan <- function(E_full, stan_file = "../code/Stan_file.stan", B_upper, B_lower, exclude = NULL, 
					 chains, cores, iter, thin, warmup, seed=122457){
	# make sure it has a names() field present
	if(!is_tibble(E) | !is.data.frame(E)){
		stop("E must be a tibble or data frame")
	}	
	# assign comm_ids, for outputting
	E_full <- assign_labels(E)
	E <- E_full %>% select(-comm_id, -rep_id)
	# stan options
	options(mc.cores = parallel::detectCores())	
	Sys.setenv(USE_CXX14 = 1)	
	rstan_options(auto_write = TRUE)
	# make sure the comm labels are in the right order
	exclude <- as.character(sapply(exclude, function(x) paste(sort(strsplit(x,"-")[[1]]),collapse="-")))
	if(any(!exclude%in%unique(E_full$comm_id))){
		stop("`exclude' contains unobserved communities")
	}
	nspp <- ncol(E)
	sp_names <- colnames(E)
	# create indicators for upper/lower bounds
	if(missing(B_upper)){
		B_upper <- matrix(1,nspp,nspp)		
		diag(B_upper) <- 0
		if("fa"%in%sp_names){
			# constrain the plant community to be competitive
			B_upper <- matrix(0,nspp,nspp) 
		}
	}
	if(missing(B_lower)){
		B_lower <- matrix(-1,nspp,nspp)
	}
	# get the intial B matrix, assuming a diagonal matrix
	B_init <- -diag(1/apply(E %>% as.matrix(), 2, function(x) quantile(x[x>0],0.9)))
	# set the s.d. for the prior to be the largest diagonal element
	maxB <- max(abs(B_init))		
	# build stan data list
	stan_data <- list(  N = ncol(E),
						M = nrow(E),
						E = E,
						B_upper = B_upper,
						B_lower = B_lower,
						y = rep(-1,nspp),
						maxB = maxB)
	# initialization vector
	init_list <- replicate(chains,list(sigmax = rep(0.25,nspp), B = B_init), simplify = FALSE)	
	# fit the stan model!
	stan_fit <- stan(stan_file, data=stan_data, cores = cores, iter=iter, thin = thin, warmup=warmup, 
					  chains = chains, seed=seed, init=init_list)
	return(list(E_full = E_full, E = E, stan_fit=stan_fit, exclude = exclude, B_upper = B_upper, B_lower = B_lower, B_init = B_init, seed = seed))
}
			
# estimate B via the naive linear regression approach
fit_lm <- function(E){
	# make sure it has a names() field present
	if(!is_tibble(E) | !is.data.frame(E)){
		stop("E must be a tibble or data frame")
	}	
	# assign comm_ids, for outputting
	E_full <- assign_labels(E)	
	# get the median value across 
	Emean <- E_full %>% gather(species, abundance, -comm_id, -rep_num) %>% filter(abundance>0) %>% group_by(species, comm_id) %>% 
		summarize(abundance = median(abundance)) %>% spread(species, abundance, fill=0) %>% arrange(desc(comm_id)) %>% select(-comm_id) %>% as.matrix()
	# get the naive regression estimate of B
	B_lm <- matrix(NA,nspp,nspp)
	for(i in 1:nspp){
		myE <- Emean[Emean[,i]>0,]
		B_lm[i,] <- as.numeric(lm(rep(-1,nrow(myE))~-1+myE)$coefficients)
	}
	return(B_lm)
}
	
# plot the mcmc results, along with histograms and obs vs. pred
plot_diagnostics <- function(stan_results){
	with(stan_results, {
		sp_names <- colnames(E)
		nspp <- length(sp_names)
		# plot the mcmc runs
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
		B_stan <- stan_fit[,1:nspp^2]
		sigma_stan <- stan_fit[,grepl("sigma",colnames(stan_fit)) | grepl("lp",colnames(stan_fit))]
		# get the predicted and observed eq abundances, using the median B
		pred_eq <- NULL
		my_B <- matrix(apply(B_stan,2,median),nspp,nspp)
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

# boostrap the mcmc results and calculate predicted abundances and probability of coexistence
bootstrap_results <- function(stan_results, nboot){
	with(stan_results, {
		nspp <- ncol(E)
		sp_names <- colnames(E)
		# create a matrix of communities presence/absence and label name
		labels <- expand.grid(replicate(nspp, 0:1, simplify = FALSE)) %>% filter(rowSums(.)>0) %>% as_tibble %>% setNames(sp_names) %>% 
			mutate(comm_id = as.numeric((.>0)%*%(2^(0:(ncol(.)-1))))) %>% gather(species,abundance,-comm_id) %>% filter(abundance>0) %>% 
			select(-abundance) %>% group_by(comm_id) %>% summarize(label=paste(species,collapse="-")) %>% left_join(expand.grid(replicate(nspp, 0:1, simplify = FALSE)) %>% filter(rowSums(.)>0) %>% as_tibble %>% setNames(sp_names) %>% 
			mutate(comm_id = as.numeric((.>0)%*%(2^(0:(ncol(.)-1))))), "comm_id") %>% select(-comm_id) %>% rename(comm_id = label)
		# if doing out of fit, only predict the unobserved
		if(any(exclude%in%labels$comm_id)){
			labels <- labels %>% filter(comm_id%in%exclude)
		}
		# conver the stan results to a matrix
		stan_mat <- as.matrix(stan_fit)
		# extract the B entries
		B_stan <- stan_mat[,1:nspp^2]
		# extract the sigma entries
		sig_stan <- stan_mat[,paste0("sigmax[",1:nspp,"]")]
		abund_mat <- coexist_mat <- tibble()
		# cycle through all the communities
		for(i in 1:nrow(labels)){
			print(paste("Bootstrapping community",i,"of", nrow(labels)))
			# get the current community id and species locations
			my_id <- labels$comm_id[i]
			my_sp <- (1:nspp)[labels %>% select(-comm_id) %>% slice(i) %>% unlist %>% as.logical()]
			n_coex <- 0
			for(j in 1:nboot){
				# sample an index
				ind <- sample(1:nrow(B_stan),1)
				# get the B matrix
				my_B <- matrix(B_stan[ind,],nspp,nspp)[my_sp, my_sp]
				# get the sigmas
				my_sig <- as.numeric(sig_stan[ind,my_sp])
				# calculate the solution
				my_x <- -rowSums(solve(my_B))
				my_x_noise <- rep(NA, length(my_x))
				if(all(my_x>0)){
					# if all species are present, sample from lognormal with sigma					
					my_x_noise <- as.numeric(apply(cbind(my_x, my_sig), 1, function(x) rlnorm(1, log(x[1]), x[2])))
					# add one to the coexistence counter
					n_coex <- n_coex+1
				}
				# add to abund_mat
				abund_mat <- bind_rows(abund_mat,tibble(comm_id = my_id, species = sp_names[my_sp], abundance = my_x, abundance_noise = my_x_noise, coexist = as.numeric(all(my_x>0)), out_fit = as.numeric(my_id%in%exclude)))			
			}
			# get the proportion that go extinct
			coexist_mat <- bind_rows(coexist_mat, tibble(comm_id = my_id, num_coexist = n_coex, num_comms = nboot) %>% mutate(prop_extinct = (num_comms - n_coex)/num_comms))
		}
		return(append(stan_results, list(boot_X = abund_mat, boot_C = coexist_mat)))
	})
}


# plot the boostrap results, both violin and obs vs. pred
plot_boot_results <- function(boot_results){
	with(boot_results, {
		# relevel the comm_ids by size
		boot_X <- boot_X %>% mutate(comm_id = factor(comm_id, level=boot_X %>%  select(comm_id) %>% distinct() %>% mutate(len = nchar(comm_id)) %>% arrange(len, comm_id) %>% select(-len) %>% unlist %>% as.character()))
		# gather the observed datapoints, filter to the current fitted values, and relevel comm_id
		obs_X_filt <- E_full %>% filter(comm_id%in%unique(boot_X$comm_id)) %>% select(-rep_num)  %>% gather(species, abundance, -comm_id) 
		obs_X <- obs_X_filt %>% mutate(comm_id = factor(comm_id, level  = E_full %>%  select(comm_id) %>% distinct() %>% mutate(len = nchar(comm_id)) %>%
				arrange(len, comm_id) %>% select(-len) %>% unlist %>% as.character())) %>% filter(abundance>0)
		# plot violins
		show(g1 <- ggplot(boot_X %>% filter(coexist==1), aes(y=abundance_noise, x= species, fill=species))+geom_violin(alpha=0.8)+
			geom_point(data=obs_X, aes(x=species, y = abundance), color="black",size=1,shape=5, alpha=0.6) +
			facet_grid(.~comm_id, scales="free_x", space = "free_x")+scale_y_log10()+theme_bw()+ylab("Abundance"))+
				theme(axis.title.x = element_blank(), legend.position = "none")
		# get median obs vs. pred values
		obs_pred <- obs_X_filt %>% group_by(comm_id, species) %>% filter(abundance>0) %>% summarize(Observed = median(abundance)) %>% 
				left_join(boot_X %>% mutate(comm_id = as.character(comm_id)) %>% filter(coexist==1, abundance>0) %>% group_by(comm_id, species) %>% summarize(Predicted = median(abundance)), c("comm_id","species"))
		# plot
		show(g2 <- ggplot(obs_pred, aes(x=Observed, y=Predicted, fill=species)) + geom_point(shape=23, alpha=0.6, size=3) + geom_abline(intercept=0, slope=1)+theme_bw()+scale_y_log10()+scale_x_log10())
			 
		return(list(g1 = g1, g2 = g2))
	})
}

