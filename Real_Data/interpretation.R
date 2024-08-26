

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



args.input = c("tertLower", "FALSE", "TRUE", "RSEM", "brca_tcga", "zscore")

response = args.input[1]
favorable = as.logical(args.input[2])
BRCA1.include = as.logical(args.input[3])
mrna.method = args.input[4]
id = args.input[5]
score.type = args.input[6]


library(caret)
library(cluster)
library(dummies)
library(rpart)
library(glmnet)


###############################################################################
##set the parameters
###############################################################################

paths = "C:/Users/linyn/Desktop/FILTER_codes/Real_Data/"
main.path = paths[dir.exists(paths)]
setwd(main.path)
out.path = paste0(main.path, "outputs/")

fold = 10

##read the data
data.set = paste0(id, "_prognostic-genes", ifelse(favorable, "-favorable", ""), "_", response, ifelse(BRCA1.include, "", "_woBRCA1"), "_", mrna.method, "_", score.type)
save.list = c("combined_results", "selected", "beta_hat", "beta0_hat", "selected_index",
		"model", "lambda", "cv", "tune_type", "X_model", "y_model", "tune", "non.genes",
		"variable_levels", "XX", "centers", "var_split", "one_prob","non.mono.idx",
		"X", "y", "factor_vars", "data.set", "response", "mrna.method", "K0", "B", "data.set", "pred.p")


##read data
file = paste0(data.set, ".csv")
data = read.csv(file)
setwd(out.path)


##the index of factor variables, which are no need to be discreted
pp = ncol(data)
cols.idx = (pp-7):(pp-1)
X = data[,-c(1, cols.idx)]
data.het = data[,cols.idx]
valid.idx = which(complete.cases(X))
data.het = data.het[valid.idx,]
X = X[valid.idx,]
y = as.factor(data$y[valid.idx])
n = nrow(X)
p = ncol(X)
K = 100	#floor(sqrt(n))



Rfile = paste0(data.set, ".RData")
if( file.exists(Rfile) ) file.remove(Rfile)


if( file.exists(Rfile) ){
	load(Rfile)
}else{

factor_vars = (p-1):p
one_probs = NULL
if( !is.null(factor_vars) ) for(i in factor_vars) X[,i] = as.factor(X[,i])

var_index = 1:ncol(X)
one_prob = 0.5

##using all the data to get the score system
X_train = X
y_train = y



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
		tmp = model.split[rows.idx,]$improve
		tmp.idx = sort(tmp, decreasing=T, index.return=T)$ix[1:min(K, length(tmp))]
		var_split[[i]] = sort(model.split[tmp.idx,]$index)
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

X_model = X_tilde
y_model = y



##############################################################
##Using penalized likelihood to estimate the coefficients
##with fusion penalty
##############################################################


##cross validation to choose best lambda
nlambda = 100
fold = 10
tune = "acc"


## using ACC as the criteria
set.seed(124)	# 124, 111, 124, 8
if(tune=="auc"){
	##using AUC as the criteria
	cv = cv.glmnet(x=X_model, y=y_model, family="binomial", type.measure="auc", nfolds=fold, nlambda=nlambda)
}else if(tune=="likelihood"){
	##using likelihood as the criteria
	cv = cv.glmnet(x=X_model, y=y_model, family="binomial", type.measure="deviance", nfolds=fold, nlambda=nlambda)
}else if(tune=="acc"){
	##using accuracy as the criteria
	cv = cv.glmnet(x=X_model, y=y_model, family="binomial", type.measure="class", nfolds=fold, nlambda=nlambda)
}


tune_type = "optimal"
# choose lambda with given criteria
least = which(cv$lambda==cv$lambda.min)
if(tune_type=="onese"){
	lambda = cv$lambda.1se
}else if(tune_type=="twose"){
	rule = cv$cvm[least] - 2*cv$cvsd[least] / sqrt(fold)
	nonzeros = cv$nzero[which(cv$cvm>rule)]
	jump = which(names(cv$nzero)==names(nonzeros)[1])
	lambda = cv$lambda[jump]
}else if(tune_type=="threese"){
	rule = cv$cvm[least] - 3*cv$cvsd[least] / sqrt(fold)
	nonzeros = cv$nzero[which(cv$cvm>rule)]
	jump = which(names(cv$nzero)==names(nonzeros)[1])
	lambda = cv$lambda[jump]
}else if(tune_type=="threesigma"){
	rule = cv$cvm[least] - 3*cv$cvsd[least]
	nonzeros = cv$nzero[which(cv$cvm>rule)]
	jump = which(names(cv$nzero)==names(nonzeros)[1])
	lambda = cv$lambda[jump]
}else if(tune_type=="optimal"){
	lambda = cv$lambda.min
}

# fit the model
model = glmnet(x=X_model, y=y_model, family="binomial", alpha=1, nlambda=nlambda)
gamma_hat = coef(model, s=lambda)

# get the estimated coefficients
beta_hat = as.vector(D%*%gamma_hat[-1])
names(beta_hat) = colnames(XXX)
beta0_hat = gamma_hat[1]

pred.p = as.vector(predict(model, newx=X_model, type="response", s=lambda))

##################################################
##transform to a score system
##################################################


##combine data with the same adjoint coefficient levels
combined_results = list()
start_tag = 1
#with_centers = sort(c(factor_vars, as.numeric(sapply(names(centers), function(x){gsub(pattern='x', replacement='', x)}))))
for(i in 1:length(variable_levels))
{
	levels = variable_levels[i]
	end_tag = start_tag + levels - 1
	betas = beta_hat[start_tag:end_tag]
	var.name.cur = unlist(strsplit(names(betas[1]), split="_"))[1]
	thr = abs(diff(range(betas)))/6

	# there are some variable without any split point by CART procedure
	if(i %in% factor_vars){
		center = NULL
	}else{
		center = centers[[var.name.cur]]
	}

	combined_beta = betas[1]
	combined_center = NULL

	for(b in 2:length(betas))
	{
		#if(betas[b]!=betas[b-1]){
		if(abs(betas[b]-betas[b-1])>thr){
			combined_beta = c(combined_beta, betas[b])
			if(!is.null(center)){
				combined_center = c(combined_center, center[b-1])
				#combined_center = c(combined_center, center[b])
			}
		}
	}

	if( !is.null(combined_center) & FALSE ){
		thr.centers = abs(diff(range(combined_center)))/3
		dropped.idx = which(abs(diff(combined_center))<thr.centers)
		if( length(dropped.idx)>0 ){
			combined_center = combined_center[-dropped.idx]	
			combined_beta = combined_beta[-(dropped.idx+1)]
		}
	}

	sign = if( length(combined_beta)>1 ){
		ifelse(all(diff(combined_beta)<0), -1, ifelse(all(diff(combined_beta)>0), 1, 2))
	}else{
		0
	}

	#combined_beta = combined_beta + abs(min(combined_beta))

	result = list(beta=combined_beta, center=combined_center, sign=sign)
	combined_results = c(combined_results, list(result))

	start_tag = start_tag + levels
}
names(combined_results) = colnames(XX)

combined_results$MATK
combined_results$RACGAP1


# find the selected variables
selected = list()
for(i in 1:length(combined_results))
{
	if(length(combined_results[[i]]$beta)!=1){
		selected = c(selected, list(combined_results[[i]]))
		names(selected)[length(selected)] = names(combined_results)[i]
	}
}
length(selected) 	# number of all selected variables


# find the variables with more than 2 levels
thresh = 2
selected_index = NULL
for(i in 1:length(combined_results)) if(length(combined_results[[i]]$beta)>thresh) selected_index = c(selected_index, i)
length(selected_index)  # number of selected variables with more than 2 levels

non.genes = NULL
non.mono.idx = which(sapply(selected, function(x)x$sign)==2)
for(i in non.mono.idx) non.genes = c(non.genes, names(selected)[i])

save.list = ls()
save(list=save.list, file=Rfile)

}# end of file.exists(Rfile)



if(TRUE){
cat(sprintf("Score Type: %s,\n Tuning Criterion: %s, Lambda: %s,\n p: %d, n: %d,\n selected: %d, selected_2: %d, non-monotone: %d,\n non-genes: %s\n", 
	score.type, tune, tune_type, p, n, length(selected), length(selected_index), length(non.mono.idx),
	paste0(non.genes, collapse=", ")))
}

# heterogeneity analysis
library(survival)
library(survminer)
library("scales")
library(ggsci)
library(ggplot2)

ranks <- vars <- last.tmp <- NULL
r = 1
for(i in 1:ncol(model$beta))
{
	beta.cur = model$beta[,i]
	tmp = names(beta_hat)[which(beta.cur!=0)]
	tmp = unique(as.vector(sapply(tmp, function(s)unlist(strsplit(s, split="_"))[1])))
	tmp.diff = setdiff(tmp, last.tmp)

	
	if( any(!(tmp.diff %in% vars)) ){
		vars.cur = tmp.diff[!(tmp.diff %in% vars)]
		vars = c(vars, vars.cur)
		ranks = c(ranks, rep(r, length(vars.cur)))
		r = r + 1
	}
	last.tmp = tmp
}


tmp = cbind(vars, ranks)
selected.order = which(vars %in% names(selected))
tmp = tmp[selected.order,]
head(tmp, 10)
write.csv(head(tmp, 10), file=paste0(id, "_selected_genes.csv"))


# KM plots
genes.first = c("RACGAP1", "MATK")
count = 1
for(gene in genes.first)
{
#gene = genes.first[count]
	gene.idx = which(names(combined_results)==gene)

	# extract data
	result.cur = combined_results[[gene.idx]]
	result.cur
	center.cur = result.cur$center

	levels = X[,gene.idx]
	groups = cut(levels, breaks=c(-Inf, center.cur, Inf), labels=1:(length(center.cur)+1))
	table(groups)

	# all subjects
	df = data.frame(y=data.het$os, time=data.het$os.month, group=factor(groups))
	df = df[complete.cases(df),]
	df$y = as.numeric(sapply(df$y, function(s)unlist(strsplit(s, split=":"))[1]))
	summary(df$time)
	df = df[df$time<(12*10),]

	fp = 2
	pdf.fp = 1.5
	plot.name = paste0(id, "_10yKM-", gene, ".pdf")
	print(plot.name)
	pdf(plot.name, onefile=FALSE, width=7*pdf.fp, height=7*pdf.fp)
	#par(mar=c(5.1+5,4.1,2.1,2.1))
	km.fit = survfit(Surv(time, y)~group, data=df)
	pair.p = pairwise_survdiff(Surv(time, y)~group, data=df, p.adjust.method="none")
	tab.pair.p = pair.p$p.value
	vec.pair.p <- vec.pair.p.names <- NULL
	for(i in 1:ncol(tab.pair.p))
	{
		for(j in i:nrow(tab.pair.p))
		{
			vec.pair.p = c(vec.pair.p, tab.pair.p[j,i])
			vec.pair.p.names = c(vec.pair.p.names, paste0(colnames(tab.pair.p)[i], rownames(tab.pair.p)[j]))
		}
	}

	cols.pool = pal_jco("default")(10)
	if( count==1 ){
		cols.set = cols.pool[1:3]
	}else if( count==2 ){
		cols.set = cols.pool[c(4,6,7)]
	}

	p = ggsurvplot(km.fit,
	           conf.int=FALSE, # add confidence intervals
	           pval=FALSE, # show the p-value for the log-rank test
	           risk.table=FALSE, # show a risk table below the plot
	           legend.labs=levels(df$group), # change group labels
	           legend.title="Group",  # add legend title
	           palette=cols.set, 
	           lwd=3*fp, pval.size=8*fp, 
		     font.legend=c(size=20*fp, color="black", face="bold"),
		     font.x=c(size=15*fp, color="black", face="bold"),
		     font.y=c(size=15*fp, color="black", face="bold"),
		     font.tickslab=c(size=15*fp, color="black", face="bold"),
	           title="", xlab="Month")
	vec.pair.p = sprintf("%.3f", round(vec.pair.p, 3))

	if( count==1 ){
		pval.1 = bquote(P[.(vec.pair.p.names[1])]^R==.(vec.pair.p[1]))
		pval.2 = bquote(P[.(vec.pair.p.names[2])]^R==.(vec.pair.p[2]))
		pval.3 = bquote(P[.(vec.pair.p.names[3])]^R==.(vec.pair.p[3]))
	}else if( count==2 ){
		pval.1 = bquote(P[.(vec.pair.p.names[1])]^M==.(vec.pair.p[1]))
		pval.2 = bquote(P[.(vec.pair.p.names[2])]^M==.(vec.pair.p[2]))
		pval.3 = bquote(P[.(vec.pair.p.names[3])]^M==.(vec.pair.p[3]))
	}
	pval.gap = 0.075*fp
	pval.size = 8*fp
	pval.x = 10*fp+2
	y.base = 0.015*fp
	print(
	p$plot+theme(axis.title.y=element_text(vjust=1*fp, size=20*fp),
			axis.text.y=element_text(margin=margin(r=5*fp)),
			axis.title.x=element_text(vjust=-0.05*fp, size=20*fp),
			axis.text.x=element_text(margin=margin(t=5*fp)),
			axis.ticks=element_line(size=1.5*fp),
			axis.ticks.length=unit(0.15*fp, "cm"))+
		theme(plot.margin = margin(2.1,2.1+10,15-10,15))+
		annotate("text", x=pval.x, y=2*pval.gap+y.base, label=pval.1, size=pval.size)+
		annotate("text", x=pval.x, y=1*pval.gap+y.base, label=pval.2, size=pval.size)+
		annotate("text", x=pval.x, y=0*pval.gap+y.base, label=pval.3, size=pval.size)
	)
	dev.off()

	count = count + 1
}



# plot coefs
L = 5*1.5
cex = 3
lwd = 15
font = 2
count = 1
for(i in selected_index)
{
	#i = selected_index[1]
	cur = combined_results[[i]]
	name.cur = names(combined_results)[i]
	range = range(X[,i])
	ys = cur$beta
	x = sort(cur$center)
	if(name.cur=="ATP5F1B"){
		xs = c(x[1]-10*diff(x)[1], x, tail(x,1)+tail(diff(x),1))
	}else{
		xs = c(x[1]-2*diff(x)[1], x, tail(x,1)+2*tail(diff(x),1))
	}
	ylim = range(ys)
	ylim[1] = ylim[1]-diff(ylim)/5
	ylim[2] = ylim[2]+diff(ylim)/5
	xlim = range(xs)

	color = "black"
	lty = 1
	plot.name = paste0(id, "_coefs_", name.cur, ".pdf")
	pdf(plot.name, height=L, width=L)
	#par(mar=c(5.1+3,4.1+2.2,2.1,2.1), mgp=c(6,1.5,0))
	par(mar=c(5.1-1,4.1,2.1,2.1), mgp=c(3, 1, 0))
	plot(0, xlim=xlim, ylim=ylim, col="white", xlab="", ylab="",
		cex.axis=cex, cex.lab=cex, font.lab=font, font.axis=font, xaxt='n')
	x.axis.ticks = axis(side = 1, at=NULL, labels=FALSE)
	axis(1, at=x.axis.ticks, labels=x.axis.ticks, tick=FALSE, line=1.2, cex.axis=cex, font.axis=font)
	segments(xs[-length(xs)], ys, xs[-1], ys, lwd=lwd, col=color, lty=lty)
	dev.off()

	
	if( name.cur %in% c("RACGAP1", "MATK")  ){
		if( name.cur=="RACGAP1" ){
			count = 1
			cols.pool = pal_jco("default")(10)
			color = cols.pool[1:3]
			lty = 1
		}else if( name.cur=="MATK" ){
			count = 2
			cols.pool = pal_jco("default")(10)
			color = cols.pool[c(4,6,7)]
			lty = 1
		}	

		plot.name = paste0(id, "_coefs_", name.cur, "_colored.pdf")
		pdf(plot.name, height=L, width=L)
		par(mar=c(5.1+3,4.1+2.2,2.1,2.1), mgp=c(6,1.5,0))
		plot(0, xlim=xlim, ylim=ylim, col="white", xlab=name.cur, ylab="",
			cex.axis=cex, cex.lab=cex, font.lab=font, font.axis=font, xaxt='n')
		x.axis.ticks = axis(side = 1, at=NULL, labels=FALSE)
		axis(1, at=x.axis.ticks, labels=x.axis.ticks, tick=FALSE, line=1.2, cex.axis=cex, font.axis=font)
		segments(xs[-length(xs)], ys, xs[-1], ys, lwd=lwd, col=color, lty=lty)
		dev.off()
	}
}

selected.morethan2 = sapply(selected_index, function(i)names(combined_results)[i])
write.csv(selected.morethan2, file="selected_3.csv")





