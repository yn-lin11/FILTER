
library(MASS)
library(dummies)
library(rpart)
library(glmnet)
library(AUC)
library(randomForest)
library(scoringRules)
library(doSNOW)
library(parallel)


## this version is the negative Brier score in 2007 JASA paper
Brier = function(y, p){
	return(1 + p^2 + (1-p)^2 - 2*dbinom(y, 1, p))
}

compute_scores = function(y, p)
{
	if(!is.numeric(y)) y = as.numeric(y)-1

	scores_logs = -logs_binom(y=y, size=1, prob=p)
	scores_crps = -crps_binom(y=y, size=1, prob=p)
	scores_brier = -Brier(y=y, p=p)

	results = cbind(scores_logs, scores_crps, scores_brier)
	colnames(results) = c("Logs", "CRPS", "Brier")
	Inf_index = which(scores_logs==-Inf | scores_crps==-Inf | scores_brier==-Inf)
	if(length(Inf_index)>0)	results = results[-Inf_index,]
	results = colMeans(results)	

	return(results)
}

##evaluate function for prediction performance
eva = function(y_test, y_pred, p_pred)
{
	library(AUC)

	table = table(y_test, y_pred)
	if(ncol(table)==1){
		if(colnames(table)=="0"){
			table = cbind(table, c(0, 0))
		}else if(colnames(table)=="1"){
			table = cbind(c(0, 0), table)
		}
	}
	acc = sum(diag(table)) / length(y_test)

	##Recall, sensitivity
	Sensitivity = table[2,2] / (table[2,2] + table[2,1])
		
	##Precision
	if((table[2,2] + table[1,2])==0){
		Precision = 0
	}else{
		Precision = table[2,2] / (table[2,2] + table[1,2])
	}
		
	##FDR, equals to (1 - precision)
	FDR = 1 - Precision
	
	##True negative rate, specificity
	Specificity = table[1,1] / (table[1,1] + table[1,2])
		
	##False Positive Rate, equals to (1 - specificity)
	FPR = table[1,2] / (table[1,2] + table[1,1])
	
	##AUC by ROC
	AUC = auc(roc(p_pred, labels=y_test))

	# proper scoring rules
	scores = compute_scores(y_test, p_pred)
	
	index = c(ACC=acc, AUC=AUC, scores, Sensitivity=Sensitivity, Precision=Precision, FDR=FDR, Specificity=Specificity, FPR=FPR)
	return(index)
}


rpart_splits <- function(fit, digits = getOption("digits")) {
  splits <- fit$splits
  if (!is.null(splits)) {
    ff <- fit$frame
    is.leaf <- ff$var == "<leaf>"
    n <- nrow(splits)
    nn <- ff$ncompete + ff$nsurrogate + !is.leaf
    ix <- cumsum(c(1L, nn))
    ix_prim <- unlist(mapply(ix, ix + c(ff$ncompete, 0), FUN = seq, SIMPLIFY = F))
    type <- rep.int("surrogate", n)
    type[ix_prim[ix_prim <= n]] <- "primary"
    type[ix[ix <= n]] <- "main"
    left <- character(nrow(splits))
    side <- splits[, 2L]
    for (i in seq_along(left)) {
      left[i] <- if (side[i] == -1L)
                   paste("<", format(signif(splits[i, 4L], digits)))
                 else if (side[i] == 1L)
                   paste(">=", format(signif(splits[i, 4L], digits)))
                 else {
                   catside <- fit$csplit[splits[i, 4L], 1:side[i]]
                   paste(c("L", "-", "R")[catside], collapse = "", sep = "")
                 }
    }
    cbind(data.frame(var = rownames(splits),
                     type = type,
                     node = rep(as.integer(row.names(ff)), times = nn),
                     ix = rep(seq_len(nrow(ff)), nn),
                     left = left),
          as.data.frame(splits, row.names = F))
  }
}





###############################################################################
##set the parameters
###############################################################################

##set the main path to restore the results
main_path = "/home/linyn/FILTER/Simulations/predictions/"
if(!dir.exists(main_path)) dir.create(main_path, recursive=TRUE)
setwd(main_path)
tmp_path = paste0(main_path, "tmp/")
if(!dir.exists(tmp_path)) dir.create(tmp_path, recursive=TRUE)


args = commandArgs(T)
#args = c("100", "0.5", "5", "threshold", "aligned")

# n: sample size
# rho: parameter for the covariance matrix of the noise vector
# ncore: number of cores used in parallel computations
# type: model family, 'threshold' for Table 1 of Model (I);'piecewise' for Table 2 of Model (II)
# beta.type: coefficient type, 'aligned' for Aligned Case; 'cross' for Alternating Case



##set the model parameters
n = as.numeric(args[1])
p0 = 6				# number of non-zero elements in the true coefficient vector
p = 500				# dimension of the true coefficient vector
rho = as.numeric(args[2])	# {0, 0.5}
M = 500				# number of experiments' replications
coef_level = 1			# signal strength, default is 1
CV = 10
one_prob = 0.5
ncore = as.numeric(args[3])
type = args[4]			# {threshold, piecewise}
beta.type = args[5]		# {aligned, cross, random}
X.type = "Gaussian"		# {Gaussian, Tdistr}: type of covariates

tmp_list = c("X", "y", "tab.index", "centers", "n", "p", "p0", "rho", "type", "beta.type", "X.type", 
		"p_fusion", "y_fusion", "index_fusion",
		"p_cart", "thresh_cart", "y_cart", "index_cart",
		"y_rf", "p_rf", "index_rf", 
		"y_la", "p_la", "index_la", 
		"y_glm", "p_glm", "index_glm", 
		"y_glm_only", "p_glm_only", "index_glm_only")


##main effects, generated from a multivariates normal distribution,
##which has mu0=0 mean and (Sigma0)_{i,j}=\rho^{|i-j|} covariance matrix
mu0 = rep(0, p)

Sigma0 = diag(rep(0.5,p))
for(i in 1:(p-1))
{
	for(j in (i+1):p)
	{
		Sigma0[i,j] = rho^(abs(i-j))
	}
}
Sigma0 = Sigma0 + t(Sigma0)



start = proc.time()
fun.list = as.vector(lsf.str())
ncore = min(ncore,detectCores())
cl = snow::makeCluster(ncore, type="SOCK")  
snow::clusterExport(cl, list=fun.list)
registerDoSNOW(cl)

ndiv = M/5
progress = function(x) if( x %% ndiv == 0 ) cat(sprintf("task %d is completed.\n", x))
opts = list(progress=progress)
res = foreach(ik=1:M, 
		.packages=c("MASS", "scoringRules", "caret", "cluster", "rpart", "glmnet", "AUC", "randomForest", "dummies"),
		.options.snow=opts, .combine="rbind") %dopar%{

	setwd(tmp_path)
	tmp_file = paste0(type, "_", beta.type, "_", X.type, "_n", n, "_p", p, "_p0", p0, "_rho", rho*10, "_seed", ik, ".RData")
	if( file.exists(tmp_file) ){
		load(tmp_file)
		tab.index
	}else{
		print(paste0("Loop: ", ik))
		tab.index = NULL

		#ik=1
		set.seed(ik*1)
		if( X.type=="Gaussian" ){
			X = mvrnorm(n=n, mu0, Sigma0)
		}else if( X.type=="Tdistr" ){
			X = mvtnorm::rmvt(n, sigma=1/3*Sigma0, df=3, delta=mu0)
		}

		##generate the cut-points
		class_number = 4
		cut_quantiles = seq(1/class_number, by=1/class_number, length.out=class_number-1)
		cut_points = qnorm(cut_quantiles)
		X_factor = data.frame(lapply(data.frame(X), function(x){cut(x, breaks=c(-Inf, cut_points, Inf), labels=1:class_number)}))
		cnames = sapply(1:(p), function(x){paste("X", x, sep='')})
		colnames(X_factor) = cnames
		colnames(X) = cnames

		##encode the X_factor to dummy variables
		XX_true = dummy.data.frame(data=X_factor, names=cnames)
		XX_int = XX_true

		##genereate the true coefficient
		betas = coef_level * c(0, 1, 2, 1)
		beta_true = NULL
		if( beta.type=="cross" ){
			for(i in 1:p0) beta_true = c(beta_true, (-1)^i*2*i*betas)
		}else if( beta.type=="aligned" ){
			#beta_true = rep(betas, p0)
			for(i in 1:p0) beta_true = c(beta_true, 2*i*betas)
		}else if( beta.type=="random" ){
			set.seed(500+ik*11)
			for(i in 1:p0) beta_true = c(beta_true, sample(c(-1,1), size=1)*betas)
		}
		beta_true = c(beta_true, rep(0, (p-p0)*class_number))

		##generate the ys by the logistic model
		if( type=="threshold" ){
			linear_part = apply(XX_int, 1, function(x){sum(x*beta_true)})
		}else if( type=="piecewise" ){
			linear_part = NULL
			for(i in 1:nrow(X))
			{
				x = X[i,]
				z = XX_int[i,]
				x.vec = as.vector(sapply(x, function(t)c(1, sin(pi*t), (t-0.5)^2, t)))
				linear_part = c(linear_part, sum(z*beta_true*x.vec))
			}
		}
		intercept = -mean(linear_part)

		p_y = intercept + linear_part
		p_y = exp(p_y) / (1+exp(p_y))

		set.seed(1000+ik*111)
		y = as.factor(sapply(p_y, function(x){rbinom(1, size=1, prob=x)}))


		####################################
		##Start training
		####################################
		set.seed(1500+ik*1111)
		train_index = sample(x=n, size=n*0.8, replace=F)

		X_train = data.frame(X[train_index,])
		y_train = y[train_index]
		X_test = data.frame(X[-train_index,])
		y_test = y[-train_index]

		##############################################################
		##Using CART+Bagging to get the cut points estimation
		##############################################################
		factor_vars = NULL
		if(is.null(factor_vars)){
			vars = c(1:ncol(X_train))
		}else{
			vars = c(1:ncol(X_train))[-factor_vars]
		}
	
		# parameters for CART
		cv.cart = 10
		complex.cart = 1e-3	#0.01
		crit.cart = "gini"

		##with some variables always in the model for estimation's robust
		y_b = y_train
		X_b = X_train

		names = colnames(X_b)
		var_split = NULL
		for(i in 1:length(vars)) var_split = c(var_split, list(NULL))
		names(var_split) = names[vars]

		removed_index = NULL
		for(i in 1:length(vars))
		{
			dt = data.frame(X_b[,i], y_b)
			colnames(dt) = c(names(var_split)[i], "y")
			model = rpart(y~., method="class", data=dt, parms=list(split=crit.cart), control=rpart.control(xval=cv.cart, cp=complex.cart))
			model.split = rpart_splits(model)
			
			if( is.null(model.split) ){
				removed_index = c(removed_index, i)
			}else{
				rows.idx = which(model.split$var==names(var_split)[i])
				var_split[[i]] = model.split[rows.idx,]$index
			}
		}
	
		if( !is.null(removed_index) ){
			var_split[removed_index] = NULL
			vars = setdiff(vars, removed_index)
		}
		#table(sapply(var_split, length))

		centers = var_split
		##encode the X into dummy variables via estimated cut point for each variable
		XX = data.frame(X[,vars])
		for(i in 1:length(vars))
		{
			x = XX[,i]
			class = rep(0, length(x))
			splits = centers[[i]]
			for(j in 1:length(splits)) class[which(x<=splits[j])] = class[which(x<=splits[j])] + 1
			XX[,i] = as.factor(length(splits) + 1 - class)
		}
		XXX = dummy.data.frame(data=data.frame(XX), names=colnames(XX))
		variable_levels = NULL
		for(i in 1:ncol(XX)) variable_levels = c(variable_levels, length(levels(XX[,i])))

		################################################################################
		##set the baseline and the levels of each variable, the default using the first
		##level as the baseline
		################################################################################
		base = rep(1, ncol(XX))
		base = c(base[1], cumsum((variable_levels))[-length(variable_levels)]+base[-1])
	
		D = matrix(0, sum(variable_levels), sum(variable_levels))
		last_length = 0
	
		for(i in 1:length(variable_levels))
		{
			sub_matrix = matrix(0, nrow=variable_levels[i], ncol=variable_levels[i])
			
			if(i==1){
				base_index = base[i]
			}else{
				base_index = base[i] - sum(variable_levels[1:(i-1)]+1)
			}   
	
			if(base_index==1){
				for(j in 1:variable_levels[i])
				{
		
					row_number = j
					if(j==1){
						sub_matrix[row_number, j] = 1
						next
					}
					
					sub_matrix[row_number, j-1] = -1
					sub_matrix[row_number, j] = 1
				}
			}else if(base_index==(variable_levels[i]+1)){
				for(j in 1:variable_levels[i])
				{
			
					row_number = j
					if(j==variable_levels[i]){
						sub_matrix[row_number, j] = 1
						next
					}
					
					sub_matrix[row_number, j] = -1
					sub_matrix[row_number, j+1] = 1
				}       
			}else{
				for(j in 1:variable_levels[i])
				{
			
					row_number = j
					if(j==base_index | j==(base_index-1)){
						sub_matrix[row_number, j] = 1
						next
					}
					
					sub_matrix[row_number, j-1] = -1
					sub_matrix[row_number, j] = 1
				}       
			}
			
			inv_submatrix = solve(sub_matrix)
	
			D[(last_length+1):(last_length+variable_levels[i]),(last_length+1):(last_length+variable_levels[i])] = inv_submatrix
			last_length = last_length + variable_levels[i]
		}

		X_tilde = as.matrix(XXX) %*% D
		X_tilde = apply(X_tilde, 2, function(x){x-mean(x)})

		X_model = X_tilde[train_index,]
		y_model = y[train_index]
		X_test2 = X_tilde[-train_index, ]
		y_test2 = y[-train_index]

		##############################################################
		##Using penalized likelihood to estimate the coefficients
		##with fusion penalty
		##############################################################

		##using accuracy as the criteria
		cv = cv.glmnet(x=X_model, y=y_model, family="binomial", type.measure="class", nfolds=5)
		lambda = cv$lambda.min

		##using glmnet with lasso penalty to fit the model for transformed data
		model = glmnet(x=X_model, y=y_model, family="binomial", alpha=1)
		gamma_hat = coef(model, s=lambda)

		##beta_hat is the estimated coefficients
		beta_hat = as.vector(D%*%gamma_hat[-1])
		names(beta_hat) = colnames(XXX)
		beta0_hat = gamma_hat[1]

		##compute the index for prediction performance
		p_fusion = as.vector(predict(model, newx=X_test2, type="response", s=lambda))
		y_fusion = as.factor(ifelse(p_fusion>one_prob, 1, 0))
		index_fusion = eva(y_test2, y_fusion, p_fusion)
		tab.index = rbind(tab.index, index_fusion)


		##############################################################
		##Other methods' performance
		##############################################################
		X_train = X[train_index,]
		y_train = y[train_index]
		X_test = X[-train_index,]
		y_test = y[-train_index]

		##CART
		dt = data.frame(X_train, y_train)
		model_cart = rpart(y_train~., method="class", data=dt, parms=list(split=crit.cart), control=rpart.control(xval=cv.cart, cp=complex.cart))

		p_cart = predict(model_cart, newdata=data.frame(X_test), type="prob")[,2]
		roc_cart = roc(p_cart, y_test)
		thresh_cart = roc_cart$cutoff[which.max(abs(roc_cart$fpr - roc_cart$tpr))]
		y_cart = as.factor(ifelse(p_cart>thresh_cart, 1, 0))
		index_cart = eva(y_test, y_cart, p_cart)
		tab.index = rbind(tab.index, index_cart)

		##RF
		dt = data.frame(X_train, y_train)
		prior = c(1-one_prob, one_prob)
		names(prior) = c("0", "1")
		set.seed(3000+ik*1)
		model_rf = randomForest(y_train~., data=dt, ntree=500, cutoff=prior, xtest=data.frame(X_test), ytest=y_test)

		y_rf = model_rf$test$predicted
		p_rf = model_rf$test$votes[,2]
		index_rf = eva(y_test, y_rf, p_rf)
		roc_rf = roc(p_rf, y_test)
		tab.index = rbind(tab.index, index_rf)

		##Logistic(l1)
		cv = cv.glmnet(x=X_train, y=y_train, family="binomial", type.measure="class", nfolds=5)
		model = glmnet(x=X_train, y=y_train, family="binomial", alpha=1, lambda=cv$lambda.min)
		beta = as.numeric(model$beta)

		p_la = predict(model, newx=X_test, type="response")
		y_la = as.factor(ifelse(p_la>one_prob, 1, 0))
		index_la = eva(y_test, y_la, p_la)
		roc_la = roc(p_la, y_test)
		tab.index = rbind(tab.index, index_la)

		##Logistic(refit)
		glm_data = data.frame(X=X_train[,which(beta!=0)], y=y_train)
		glm_model = glm(y~., data=glm_data, family="binomial")

		p_glm = predict(glm_model, newdata=data.frame(X=X_test[,which(beta!=0)]), type="response")
		y_glm = as.factor(ifelse(p_glm>one_prob, 1, 0))
		index_glm = eva(y_test, y_glm, p_glm)
		roc_glm = roc(p_glm, y_test)
		tab.index = rbind(tab.index, index_glm)

		###Logistic(only)
		glm_data = data.frame(X_train, y=y_train)
		glm_model = glm(y~., data=glm_data, family="binomial")

		p_glm_only= predict(glm_model, newdata=data.frame(X_test), type="response")
		roc_glm_only = roc(p_glm_only, y_test)
		y_glm_only = as.factor(ifelse(p_glm_only>one_prob, 1, 0))
		index_glm_only = eva(y_test, y_glm_only, p_glm_only)
		tab.index = rbind(tab.index, index_glm_only)

		filter.rows = 1
		tab.index = rbind(tab.index[-filter.rows,], tab.index[filter.rows,])
		rownames(tab.index) = c("cart", "rf", "lasso", "glm", "glm_only", "FILTER")

		save(list=tmp_list, file=tmp_file)
		tab.index
	}
}

snow::stopCluster(cl)
end = proc.time()
print(end-start)


##############################################################
##Compute the average performance for different methods
##############################################################
setwd(main_path)
dig = 3

n.M = (nrow(res)/M)
method.names = rownames(res)[1:n.M]
index.names = colnames(res)
tab.results = NULL
for(i in 1:n.M)
{
	#i = 1
	res.cur = res[((1:M)-1)*n.M+i,]
	means.cur = colMeans(res.cur)
	sds.cur = apply(res.cur, 2, sd)
	tab.cur = sapply(1:length(means.cur), function(j)sprintf(paste0("%.",dig,"f(%.",dig,"f)"), means.cur[j], sds.cur[j]))
	tab.results = rbind(tab.results, tab.cur)
}
rownames(tab.results) = method.names
colnames(tab.results) = index.names


#Results
file_name = paste0("prediction_", type, "_", beta.type, "_", X.type, "_n", n, "_p", p, "_p0", p0, "_rho", rho*10, "_M", M, ".csv")
write.csv(tab.results, file_name, row.names=TRUE)

file_name = paste0("prediction_", type, "_", beta.type, "_", X.type, "_n", n, "_p", p, "_p0", p0, "_rho", rho*10, "_M", M, ".RData")
save_list = c("res", "tab.results", "n", "p", "p0", "rho", "type", "beta.type", "X.type", "M", "coef_level")
save(list=save_list, file=file_name)

print(tab.results)

