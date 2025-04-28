compareReactionsFromDifferentModels <- function(models, id){
	
	if(length(id)==1){
		id <- rep(id, length(models))
	}
	names(id) <- models
	stopifnot(length(models)==length(id))
	
	
	s <- lapply(models,
		function(m){
			sm <- shrinkMatrix(allModels[[m]], j=id[m])
			sm[order(rownames(sm)), ,drop=F]
		}
	)
	
	# are rownames equal!?
	rn <- lapply(s, rownames)
	if(!all(table(unlist(rn))==length(models))){
		print(models)
		print(id)
		print(table(unlist(rn)))
		stop("reaction definition differ in checked models")
	}
	
	
	s <- do.call(cBind, s)
	
#	print(s)
	
	check <- apply(s, 1, function(x){
		length(unique(x))==1
	})
	
	return(all(check))
	
}
