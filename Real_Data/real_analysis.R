
#args.input = commandArgs(T)
args.input = c("tertLower", "FALSE", "TRUE", "RSEM", "brca_tcga", "zscore")

response = args.input[1]
favorable = as.logical(args.input[2])
BRCA1.include = as.logical(args.input[3])
mrna.method = args.input[4]
id = args.input[5]
score.type = args.input[6]

# evaluate function for prediction performance
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
	Precision = table[2,2] / (table[2,2] + table[1,2])
		
	##FDR, equals to (1 - precision)
	FDR = table[1,2] / (table[1,2] + table[2,2])	
	
	##True negative rate, specificity
	Specificity = table[1,1] / (table[1,1] + table[1,2])
		
	##False Positive Rate, equals to (1 - specificity)
	FPR = table[1,2] / (table[1,2] + table[1,1])
	
	##AUC by ROC
	AUC = auc(roc(p_pred, labels=y_test))

	if( is.factor(y_test) ) y_test = as.numeric(y_test) - 1
	scores_logs = mean(-logs_binom(y=y_test, size=1, prob=p_pred))
	scores_crps = mean(-crps_binom(y=y_test, size=1, prob=p_pred))
	scores_brier = mean(-Brier(y=y_test, p=p_pred))
#	if( scores_logs==-Inf ) scores_logs = NULL

	index = c(Logs=scores_logs, CRPS=scores_crps, Brier=scores_brier, ACC=acc, AUC=AUC, Sensitivity=Sensitivity, Precision=Precision, FDR=FDR, Specificity=Specificity, FPR=FPR)
	return(index)
}


## this version is the negative Brier score in 2007 JASA paper
Brier = function(y, p){
	return(1 + p^2 + (1-p)^2 - 2*dbinom(y, 1, p))
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


library(scoringRules)
library(caret)
library(cluster)
library(dummies)
library(rpart)
library(glmnet)
library(AUC)
library(randomForest)
library(doSNOW)
library(parallel)


ncore = 5

paths = "C:/Users/linyn/Desktop/FILTER_codes/Real_Data/"
main.path = paths[dir.exists(paths)]
setwd(main.path)
out.path = paste0(main.path, "outputs/")
if( !dir.exists(out.path) ) dir.create(out.path)
data.set = paste0(id, "_prognostic-genes", ifelse(favorable, "-favorable", ""), "_", response, ifelse(BRCA1.include, "", "_woBRCA1"), "_", mrna.method, "_", score.type)

# M: number of experiments' replications
# cv_train: logical, indicate whether do the cross-validation analysis for the data or not
M = 10
cv_train = TRUE
fold = M	# for real data, fold=M
file_name = paste0(data.set, "_Folds", fold, ".RData")
print(file_name)

##read data
file = paste0(data.set, ".csv")
data = read.csv(file)
pp = ncol(data)
cols.idx = (pp-7):(pp-1)
data = data[,-cols.idx]
data = data[complete.cases(data),]

setwd(out.path)
n = nrow(data)
p = ncol(data)-1

##the index of factor variables, which are no need to be discreted
factor_vars = (p-1):p
one_probs = NULL

ps_la_out = NULL
ps_la_in = NULL
index_las = NULL
index_las_in = NULL
threshes_la = NULL
rocs_la = list()

ps_glm_out = NULL
ps_glm_in = NULL
index_glms = NULL
index_glms_in = NULL
threshes_glm = NULL
rocs_glm = list()

ps_glm_only_out = NULL
ps_glm_only_in = NULL
index_glms_only = NULL
index_glms_only_in = NULL
threshes_glm_only = NULL
rocs_glm_only = list()

ps_cart_out = NULL
ps_cart_in = NULL
index_carts = NULL
index_carts_in = NULL
threshes_cart = NULL
rocs_cart = list()

ps_rf_out = NULL
ps_rf_in = NULL
index_rfs = NULL
index_rfs_in = NULL
threshes_rf = NULL
rocs_rf = list()

ys_out = NULL
ys_in = NULL

loop = (1):(M)
print(loop)

savelist = c("ik", "ps_fusion_out", "ps_fusion_in", "index_fusions", "rocs_fusion",
		"ps_la_in", "ps_la_out", "index_las", "rocs_la", "ps_glm_out", "ps_glm_in", "index_glms", "rocs_glm",
		"ps_glm_only_out", "ps_glm_only_in", "index_glms_only", "rocs_glm_only", "ys_in", "ys_out",
		"ps_cart_out", "ps_cart_in", "index_carts", "rocs_cart", "ps_rf_out", "ps_rf_in", "index_rfs", "rocs_rf",
		"ps_CB_out", "ps_CB_in", "index_CBs", "rocs_CB")

if(FALSE){
## create fold id
i.fold = 0
while(i.fold<1e3)
{
	if( i.fold %% 100 == 0 ) print(i.fold)

	#set.seed(fold*(1+i.fold))
	set.seed(fold*(i.fold))	#111+i.fold
	flds = createFolds(data$y, k=fold, list=TRUE, returnTrain=FALSE)

	y_test_sums <- valids <- NULL
	for( ik in loop )
	{
		cv_loop = (1:fold)[-ik]
		train_index = sort(unlist(sapply(cv_loop, function(i){flds[[i]]})))
		num.y.test = sum(data$y[-train_index])
		y_test_sums = c(y_test_sums, num.y.test)
		valid = (num.y.test>0 & num.y.test<length(data$y[-train_index]))
		valids = c(valids, valid)
	}
	#print(y_test_sums)

	if( all(valids) ){
		print("Valid folds founded.")
		print(y_test_sums)
		break
	}
	i.fold = i.fold + 1
}

y_test_sums = NULL
for( ik in loop )
{
	cv_loop = (1:fold)[-ik]
	train_index = sort(unlist(sapply(cv_loop, function(i){flds[[i]]})))
	y_test_sums = c(y_test_sums, sum(data$y[-train_index]))
}
rbind(print(y_test_sums),sapply(flds, length))
}

set.seed(fold+111)
flds = createFolds(data$y, k=fold, list=TRUE, returnTrain=FALSE)

start = proc.time()
fun.list = as.vector(lsf.str())
ncore = min(ncore,detectCores())
cl = snow::makeCluster(ncore, type="SOCK")  
snow::clusterExport(cl, list=fun.list)
registerDoSNOW(cl)


ndiv = M/5
progress = function(x) if( x %% ndiv == 0 ) cat(sprintf("task %d is completed.\n", x))
opts = list(progress=progress)
res.cur = foreach(ik=loop, 
		.packages=c("scoringRules", "caret", "cluster", "rpart", "glmnet", "AUC", "randomForest", "dummies"),
		.options.snow=opts, .combine="rbind") %dopar% {
			
	#ik=1
	print(paste0("Loop: ", ik))
	tab.index = NULL

	if(cv_train){
		# use cross-validation to evaluate methods
		print("CV_train")	

		cv_loop = (1:fold)[-ik]
		train_index = sort(unlist(sapply(cv_loop, function(i){flds[[i]]})))

		X = data.frame(data[,-1])
		y = as.factor(data$y)
		if( !is.null(factor_vars) ) for(i in factor_vars) X[,i] = as.factor(X[,i])

		set.seed(111+ik)
		X_train_all = data.frame(X[train_index,])
		y_train_all = y[train_index]
		
		X_test = data.frame(X[-train_index,])
		y_test = y[-train_index]
	}else{
		##divide the data into training and testing data
		set.seed(1000+(ik-1)*2)
		train_index = sample(x=n, size=n*0.8, replace=F)
		
		X = data.frame(data[,-1])
		y = as.factor(data$y)
		for(i in factor_vars)
		{
			X[,i] = as.factor(X[,i])
		}
		
		X_train_all = data.frame(X[train_index,])
		y_train_all = y[train_index]
		
		X_test = data.frame(X[-train_index,])
		y_test = y[-train_index]
	}

	var_index = 1:ncol(X)
	one_prob = 0.5


	##data used for other methods
	X_train = X_train_all[,var_index]
	y_train = y_train_all
	X_test = X_test[, var_index]


	##############################################################
	##Using CART+Bagging to get the cut points estimation
	##############################################################
	if(is.null(factor_vars)){
		vars = c(1:ncol(X_train))
	}else{
		vars = c(1:ncol(X_train))[-factor_vars]
	}
	
	# parameters for CART
	cv.cart = 10
	complex.cart = 1e-3	#1e-3
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
	XXX = dummy.data.frame(data=data.frame(XX), names=colnames(XX), sep="_")
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

	print("Start fusion CV")
	##using accuracy as the criteria
	cv = cv.glmnet(x=X_model, y=y_model, family="binomial", type.measure="class", nfolds=5)
	lambda = cv$lambda.min

	##using glmnet with lasso penalty to fit the model for transformed data
	model = glmnet(x=X_model, y=y_model, family="binomial", alpha=1)

	p_fusion = as.vector(predict(model, newx=X_test2, type="response", s=lambda))
	fit_fusion = as.vector(predict(model, newx=X_model, type="response", s=lambda))
	roc_fusion = roc(p_fusion, y_test2)
	thresh_fusion = one_prob
	y_fusion = as.factor(ifelse(p_fusion>thresh_fusion, 1, 0))
	index_fusion = eva(y_test2, y_fusion, p_fusion)
	y_fusion_in = as.factor(ifelse(fit_fusion>thresh_fusion, 1, 0))
	index_fusion_in = eva(y_model, y_fusion_in, fit_fusion)
	tab.index = rbind(tab.index, index_fusion)

	##############################################################
	##Other methods' performance
	##############################################################

	X_train = X[train_index,]
	y_train = y[train_index]
	X_test = X[-train_index,]
	y_test = y[-train_index]

	print("Start Others")
	##CART
	dt = data.frame(X_train, y_train)
	model_cart = rpart(y_train~., method="class", data=dt, parms=list(split=crit.cart), control=rpart.control(xval=cv.cart, cp=complex.cart))

	p_cart = predict(model_cart, newdata=data.frame(X_test), type="prob")[,2]
	fit_cart = predict(model_cart, newdata=data.frame(X_train), type="prob")[,2]
	roc_cart = roc(p_cart, y_test)
	thresh_cart = one_prob
	y_cart = as.factor(ifelse(p_cart>thresh_cart, 1, 0))
	index_cart = eva(y_test, y_cart, p_cart)
	y_cart_in = as.factor(ifelse(fit_cart>thresh_cart, 1, 0))
	index_cart_in = eva(y_train, y_cart_in, fit_cart)
	tab.index = rbind(tab.index, index_cart)

	##RF
	dt = data.frame(X_train, y_train)
	prior = c(1-one_prob, one_prob)
	names(prior) = c("0", "1")

	set.seed(ik*333)
	model_rf = randomForest(y_train~., data=dt, ntree=500, cutoff=prior)
	p_rf = predict(object=model_rf, newdata=data.frame(X_test), type="prob")[,2]
	fit_rf = predict(object=model_rf, newdata=data.frame(X_train), type="prob")[,2]
	roc_rf = roc(p_rf, y_test)
	y_rf = predict(object=model_rf, newdata=data.frame(X_test), type="response")
	index_rf = eva(y_test, y_rf, p_rf)
	y_rf_in = predict(object=model_rf, newdata=data.frame(X_train), type="response")
	index_rf_in = eva(y_train, y_rf_in, fit_rf)
	tab.index = rbind(tab.index, index_rf)


	##Logistic(l1)
	fold.la = 5
	set.seed(ik*555)
	cv = cv.glmnet(x=data.matrix(X_train), y=y_train, family="binomial", type.measure="class", nfolds=fold.la)
	lambda = cv$lambda.min


	model = glmnet(x=data.matrix(X_train), y=y_train, family="binomial", alpha=1)
	beta = as.numeric(coef(model, s=lambda)[-1])

	p_la = predict(model, newx=data.matrix(X_test), type="response", s=lambda)
	fit_la = predict(model, newx=data.matrix(X_train), type="response", s=lambda)
	roc_la = roc(p_la, y_test)
	thresh_la = one_prob
	y_la = as.factor(ifelse(p_la>thresh_la, 1, 0))
	index_la = eva(y_test, y_la, p_la)
	y_la_in = as.factor(ifelse(fit_la>thresh_la, 1, 0))
	index_la_in = eva(y_train, y_la_in, fit_la)
	tab.index = rbind(tab.index, index_la)

	##Logistic(refit)
	X_temp = data.frame(X_train[,which(beta!=0)])
	glm_data = data.frame(X_temp, y=y_train)
	glm_model = glm(y~., data=glm_data, family="binomial")

	X_test_temp = data.frame(X_test[,which(beta!=0)])
	colnames(X_test_temp) = colnames(glm_data)[-ncol(glm_data)]
	p_glm= predict(glm_model, newdata=X_test_temp, type="response")
	fit_glm= predict(glm_model, newdata=X_temp, type="response")
	roc_glm = roc(p_glm, y_test)
	thresh_glm = one_prob
	y_glm = as.factor(ifelse(p_glm>thresh_glm, 1, 0))
	index_glm = eva(y_test, y_glm, p_glm)
	y_glm_in = as.factor(ifelse(fit_glm>thresh_glm, 1, 0))
	index_glm_in = eva(y_train, y_glm_in, fit_glm)
	tab.index = rbind(tab.index, index_glm)

	###Logistic(only)
	glm_data = data.frame(X_train, y=y_train)
	glm_model = glm(y~., data=glm_data, family="binomial")

	p_glm_only= predict(glm_model, newdata=data.frame(X_test), type="response")
	fit_glm_only= predict(glm_model, newdata=data.frame(X_train), type="response")
	roc_glm_only = roc(p_glm_only, y_test)
	thresh_glm_only = one_prob
	y_glm_only = as.factor(ifelse(p_glm_only>thresh_glm_only, 1, 0))
	index_glm_only = eva(y_test, y_glm_only, p_glm_only)
	y_glm_only_in = as.factor(ifelse(fit_glm_only>thresh_glm_only, 1, 0))
	index_glm_only_in = eva(y_train, y_glm_only_in, fit_glm_only)
	tab.index = rbind(tab.index, index_glm_only)


	filter.rows = 1
	tab.index = rbind(tab.index[-filter.rows,], tab.index[filter.rows,])
	rownames(tab.index) = c("cart", "rf", "lasso", "glm", "glm_only", "FILTER")
	tab.index
		}
snow::stopCluster(cl)
end = proc.time()
print(end-start)




M = nrow(res.cur) / length(loop)
tab.means = NULL
for(i in loop)
{
	rows = ((i-1)*M+1):(i*M)
	tab.cur = res.cur[rows,1:5]
	tab.cur
	if( is.null(tab.means) ){
		tab.means = tab.cur
		div.mat = matrix(1, nrow=nrow(tab.means), ncol=ncol(tab.means))
		if(any(tab.cur==-Inf)){
			div.mat[tab.means==-Inf] = 0
			tab.means[tab.means==-Inf] = 0
		}
	}else{
		div.mat = div.mat + 1
		if(any(tab.cur==-Inf)){
			div.mat[tab.cur==-Inf] = div.mat[tab.cur==-Inf] - 1
			tab.cur[tab.cur==-Inf] = 0
		}
		tab.means = tab.means + tab.cur
	}
}
tab.means = tab.means / div.mat


dig = 3
results = round(tab.means, dig)[,c(4,5,1)]
results

file.name = paste0("Real_result_", data.set, "_Folds", 10, ".csv")
write.csv(results, file=file.name)


