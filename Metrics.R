
### Get Percent Agreement ###
pcagree = function(xmat, within=0){  
  output = matrix(NA, ncol=3, nrow=1)
  colnames(output) = c("n.agree","n.total","percent.agree")
  
  output[1,"n.total"] = nrow(xmat)
  output[1,"n.agree"] = sum(apply(xmat, 1, max)-apply(xmat, 1, min) <= within)
  output[1,"percent.agree"] = 100*sum(apply(xmat, 1, max)-apply(xmat, 1, min) <= within) / nrow(xmat)
  
  return(output)
}


### Get Correlation ###
cor.me = function(xmat) {
  output = matrix(NA, ncol=4, nrow=1)
  colnames(output)=c("items", "min", "average", "max")
  
  n = ncol(xmat)
  output[,"items"] = n 
  cormat = cor(xmat, use="pairwise")
  output[,"min"] = min(cormat)
  diag(cormat) = 0
  output[,"max"] = max(cormat)
  output[,"average"] = sum(cormat) / (n^2-n)
  
  return(output)
}


### Get Raw & Std Alpha ###
alpha.me = function(xdat, digits=4) { # called by 
  # compute raw and standardized alpha
  n = ncol(xdat); k=(n^2-n)/2 # number of correlation/covariance to be estimated
  cov = cov(xdat, use="pairwise.complete.obs")
  totalvar = sum(cov)
  commonvar = totalvar - sum(diag(cov)) # sum of inter-item covariance
  
  rawalpha = (n/(n-1))*(commonvar/totalvar) # raw Cronbach's alpha
  avgcov = commonvar/(2*k) # average iter-item (pair-wise) covariance
  r=(sum(cor(xdat, use="pairwise.complete.obs"))-n)/(2*k) # average inter-item correlation
  stdalpha = n*r/(1+(n-1)*r) # standardized Cronbach's alpha
  
  output = round(c(n, rawalpha, stdalpha, r, avgcov), digits)
  names(output) = c("Num.Items","Raw.Alpha", "Std.Alpha", "Avg.Cor", "Avg.Cov")
  return(output)
}

### Convertion between Different Test/Scale Lengths ###
newalpha = function(n1, n2, alphaold) {
  k = n2 / n1
  alphanew = k*alphaold / (1+(k-1)*alphaold)
  return(alphanew)
}

moreitems = function(alphaold, alphanew, n1=1) {
  k = alphanew*(1-alphaold) / (alphaold*(1-alphanew))
  n2 = k*n1
  more = n2 - n1
  
  tab = vector()
  tab = c(n1, k, n2, more)
  names(tab) = c("num.before", "kfold", "num.after", "num.more")
  
  return(tab)
}

### Convertion between Stdalpha, Avg.Iter-Item.Cor, Num.Item, NewStdAlpha ###
avgr = function(n, stdalpha){ stdalpha / (n - (n-1)*stdalpha) }

stdalpha = function(n, avgr) { n*avgr/(1+(n-1)*avgr) }

nitems = function(avgr, stdalpha) { stdalpha*(1-avgr) / (avgr*(1-stdalpha)) }

newstdalpha = function(n1, n2, stdalpha) {
  avgr = avgr(n1, stdalpha)
  output = stdalpha(n2, avgr)
  return(output)
}

### New Alpha With One Item Omitted ###
# Full Version from Stata #
alpha.item = function(xmat, deci=4) {
  
  tab=as.data.frame(matrix(NA, ncol=10, nrow=ncol(xmat)))
  colnames(tab)=c("Missing","Reversed","Raw.Alpha","Std.Alpha",
                  "Avg.Cor","Avg.Cov","r.Item.Rest","r.Item.Test",
                  "Mean","Var")
  rownames(tab)=names(xmat)
  
  tab[,"Reversed"]="..."
  for (i in 1:ncol(xmat)){
    tab[i,"Missing"]=sum(is.na(xmat[,i]))
    
    tab[i,c("Raw.Alpha","Std.Alpha","Avg.Cor","Avg.Cov")]=alpha.me(xmat[,-i])[c("Raw.Alpha","Std.Alpha","Avg.Cor","Avg.Cov")] # Raw & Std Alpha with one item omitted
    
    scale.rest = apply(xmat[,-i], 1, mean, na.rm=T) # scale formed by all other items omitting one
    tab[i,"r.Item.Rest"]=round(cor(xmat[,i], scale.rest, use="pairwise.complete.obs"), deci) # item-rest correlation btw item and scale
    
    scale.test = apply(xmat, 1, mean, na.rm=T) # scale formed by all other items omitting one
    tab[i,"r.Item.Test"]=round(cor(xmat[,i], scale.test, use="pairwise.complete.obs"), deci) # item-test correlation btw item and scale
    
    tab[i,"Mean"] = round(mean(xmat[,i],na.rm=T), deci)
    tab[i,"Var"] = round(var(xmat[,i],na.rm=T), deci)
  }
  
  return(tab)
}


# A Simplified Version #
alpha.omit = function(xmat, deci=4) {
  
  tab=as.data.frame(matrix(NA, ncol=8, nrow=ncol(xmat)))
  colnames(tab)=c("Missing", "Reversed","Raw.Alpha","Std.Alpha",
                  "Avg.Cor","r(item-rest)","Mean","Var")
  rownames(tab)=names(xmat)
  
  tab[,"Reversed"]=paste("...")
  for (i in 1:ncol(xmat)){
    tab[i,"Missing"]=sum(is.na(xmat[,i]))
    
    tab[i,c("Raw.Alpha","Std.Alpha","Avg.Cor")]=alpha.me(xmat[,-i])[c("Raw.Alpha","Std.Alpha","Avg.Cor")] # Raw & Std Alpha with one item omitted
    
    scale.rest=apply(xmat[,-i], 1, mean, na.rm=T) # scale formed by all other items omitting one
    tab[i,"r(item-rest)"]=round(cor(xmat[,i], scale.rest, use="pairwise.complete.obs"), deci) # item-rest correlation btw item and scale
    
    tab[i,"Mean"] = round(mean(xmat[,i],na.rm=T), deci)
    tab[i,"Var"] = round(var(xmat[,i],na.rm=T), deci)
  }
  
  return(tab)
}


### Get Alpha and New Alpha ###

alpha.go = function(xmat, sort="Yes", deci=4) {
  
  alpha=alpha.me(xmat)
  item = alpha.omit(xmat)
  bad=item[item[,"Raw.Alpha"]>alpha["Raw.Alpha"],]
  if (sort == "Yes") {
    item = item[order(item[,"Raw.Alpha"]),]
    bad =bad[order(bad[,"Raw.Alpha"], decreasing=T),]
  }
  
  baditems = rownames(bad)
  gooditems = setdiff(rownames(item), baditems)
  
  alpha.rm = alpha.me(xmat[,gooditems])
  
  alpha.raw = newalpha(alpha.rm["Num.Items"], alpha["Num.Items"], alpha["Raw.Alpha"])
  alpha.std = stdalpha(alpha["Num.Items"], alpha.rm["Avg.Cor"])
  alpha.inflate = c(alpha["Num.Items"], alpha.raw, alpha.std, alpha.rm["Avg.Cor"], alpha.rm["Avg.Cov"])
  names(alpha.inflate) = names(alpha)
  
  list(item.omit=item, baditems=bad, alpha=alpha, 
       alpha.rm = alpha.rm, alpha.inflate = round(alpha.inflate, deci), 
       gooditems=gooditems, baditems=baditems)
}




### Standardize Xmat ###
std.xmat = function(xmat) {
  colmean = colMeans(xmat, na.rm=T)
  colsd = sapply(xmat, sd, na.rm=T) # equivalent to apply(mat, 2, sd)
  
  xmat.mn = t(t(xmat)-colmean) # mn for mean-normalization
  xmat.sd = t((t(xmat)-colmean)/colsd) # convert all cases into z-scores
  
  return(xmat.sd)
}




### Cohen's Kappa ###

kappa.opt = function(xdat, weight=c("none","equal","squared"), alpha=0.05){
  
  require(psych)
  
  m = nrow(xdat); n = ncol(xdat); weight = match.arg(weight)
  
  # Prepare Variables #
  if(n==2 && m!=n) { # if it's not a table
    lvls = levels(as.factor(as.matrix(xdat))) # ensure all levels will be displayed from table()
    x1 = factor(xdat[,1], levels=lvls)
    x2 = factor(xdat[,2], levels=lvls)
    tab = table(x1, x2)
  }

  # Get Weighted Kappa #
  tab = tab/sum(tab) # use joint probabilities instead of counts
  rs = rowSums(tab) # rs for row sum; conditional probabilities
  cs = colSums(tab) # cs for column sum; conditional probabilities
  prob = rs%o%cs

  wt = matrix(0, ncol=length(lvls), nrow=length(lvls))
  
  if(weight=="none")
  {wt = diag(1, length(lvls))}
  
  if(weight=="squared")
  {wt = 1 - abs((col(wt) - row(wt)))^2 / (length(lvls)-1)^2}
  
  if(weight=="equal")
  {wt = 1 - abs((col(wt) - row(wt))) / (length(lvls)-1)}
   
  wpo = sum(wt*tab) # weighted probability of observation
  wpe = sum(wt*prob) # weighted probability of error # rs%*%t(cs)
  
  wkappa = (wpo - wpe)/(1 - wpe)
  
  # Get Variance #
  wrs = colSums(wt * rs)
  wcs = colSums(wt * cs)
  
  var = prob * (wt - wcs %+% t(wrs))^2
  varkappa = (sum(var) - wpe^2) / (m * (1-wpe)^2)

  z = wkappa / sqrt(varkappa)
  p = 2*(1-pnorm(abs(z)))
  
  ### Produce Output ###
  output = matrix(NA, ncol=6, nrow=1)
  colnames(output) = c("Kappa", "se", "z stat", "p-value", "lower CI", "upper CI")
  rownames(output) = paste("Weight:", weight)
  output[1,] = c(wkappa, sqrt(varkappa), z, p, wkappa+c(1,-1)*qnorm(alpha/2)*sqrt(varkappa))
  
  return(output)
}


kappa.me = function(xdat, alpha=0.05){
  
  require(psych) #getAnywhere("%+%")
  
  m = nrow(xdat); n = ncol(xdat)
  
  # Setup Variables #
  if(n==2 && m!=n) {
    lvls = levels(as.factor(as.matrix(xdat)))
    x1 = factor(xdat[,1], levels=lvls)
    x2 = factor(xdat[,2], levels=lvls)
    tab = table(x1, x2)
  }
  
  tot = sum(tab)
  tab = tab/tot
  rs = rowSums(tab)
  cs = colSums(tab)
  prob = rs%o%cs
  
  wt = matrix(0, ncol=length(lvls), nrow=length(lvls))
  
  # Get Raw Kappa #
  wt.n = diag(1, length(lvls))
  
  po.n = sum(diag(tab)) # probability of observation
  pe.n = sum(diag(prob)) # probability of error # rs%*%t(cs)
  kappa = (po.n - pe.n)/(1 - pe.n)
  
  var = prob * (wt.n - cs %+% t(rs))^2
  varkappa = (sum(var) - pe.n^2) / (m * (1-pe.n)^2)
  
  z = kappa / sqrt(varkappa)
  p = 2*(1-pnorm(abs(z)))
  
  
  # Get Equally Weighted Kappa #
  wt.e = 1 - abs((col(wt) - row(wt))) / (length(lvls)-1)
  po.e = sum(wt.e*tab) # weighted probability of observation
  pe.e = sum(wt.e*prob) # weighted probability of error # rs%*%t(cs)
  
  kappa.e = (po.e - pe.e)/(1 - pe.e)
  
  rs.e = colSums(wt.e * rs)
  cs.e = colSums(wt.e * cs)
  
  var.e = prob * (wt.e - cs.e %+% t(rs.e))^2
  varkappa.e = (sum(var.e) - pe.e^2) / (m * (1-pe.e)^2)
  
  z.e = kappa.e / sqrt(varkappa.e)
  p.e = 2*(1-pnorm(abs(z.e)))
  
  
  # Get Squaredly Weighted Kappa #
  wt.s = 1 - abs((col(wt) - row(wt)))^2 / (length(lvls)-1)^2
  po.s = sum(wt.s*tab) # weighted probability of observation
  pe.s = sum(wt.s*prob) # weighted probability of error # rs%*%t(cs)
  
  kappa.s = (po.s - pe.s)/(1 - pe.s)
  
  rs.s = colSums(wt.s * rs)
  cs.s = colSums(wt.s * cs)
  
  var.s = prob * (wt.s - cs.s %+% t(rs.s))^2
  varkappa.s = (sum(var.s) - pe.s^2) / (m * (1-pe.s)^2)
  
  z.s = kappa.s / sqrt(varkappa.s)
  p.s = 2*(1-pnorm(abs(z.s)))
  
  ### Output ###
  output = matrix(NA, ncol=6, nrow=3)
  colnames(output) = c("Kappa", "se", "z stat", "p-value", "lower CI", "upper CI")
  rownames(output) = c("No Weight", "Equal Weight", "Squared Weight")
  output["No Weight",] = c(kappa, sqrt(varkappa), z, p, kappa+c(1,-1)*qnorm(alpha/2)*sqrt(varkappa))
  output["Equal Weight",] = c(kappa.e, sqrt(varkappa.e), z.e, p.e, kappa.e+c(1,-1)*qnorm(alpha/2)*sqrt(varkappa.e))
  output["Squared Weight",] = c(kappa.s, sqrt(varkappa.s), z.s, p.s, kappa.s+c(1,-1)*qnorm(alpha/2)*sqrt(varkappa.s))
  
  return(output)
}


##################################
### Non-parametric Correlation ###
kendall = function(master, student) {
  
  yvec = rank(master, ties.method="first") # get the rank of master
  xvec = rank(student, ties.method="first") # get the rank of student
  xvec = xvec[order(yvec)] # sort xvec based on the rank of yvec from min to max 
  
  concor=NULL; discor=NULL; n = length(xvec)
  
  for (i in 1:(n-1)) {
    concor[i] = sum(xvec[i] < xvec[(i+1):length(xvec)])
    discor[i] = sum(xvec[i] > xvec[(i+1):length(xvec)]) 
  }
  
  c = sum(concor); d = sum(discor)
  S = c-d # S is known as Kendall's Score
  D = c+d # D is the number of all possible pairs = choose(n,2)
  tau.a = S/D
  
  z = (3*tau.a*sqrt(n*(n-1))) / (sqrt(4*n+10))
  
  p = 2*pnorm(z, lower.tail=F)
  
  names=c("tau.a", "p.a", "S")
  output = vector(length=length(names))
  names(output) = names
  output[1]=tau.a; output[2]=p; output[3]=S
  
  return(output)
  
}
############################################



spearman = function(master, student) {
  yvec = rank(master, ties.method="first") # get the rank of master
  xvec = rank(student, ties.method="first") # get the rank of student
  
  n = length(xvec) # total number of case pairs
  d2 = (yvec-xvec)^2
  
  rho = 1 - (6*sum(d2)/(n^3-n)) # rho for spearman's rho  
  # rho ignores minor rank differences between master and student
  # however, it punishes heavily on egregious differences
  # rho cares more about magnitude than the number of differences
  
  # tau is the opposite (number > magnitude)
  # tau cares a lot about the number of differences, but less about magnitude 
  
  return(rho)
}


###########
### ICC ###
icc.me = function(xdat) {
  # Setup Data #
  m = nrow(xdat)
  n = ncol(xdat)
  dat = stack(xdat)
  dat = data.frame(dat, id = factor(rep(1:m, 2)))
  
  # Get ANOVA Results #
  mod = summary(aov(values ~ id + ind, data=dat))
  stats = matrix(unlist(mod), nrow=3, 
                 dimnames=list(c("id","rater","resid"),
                               c("df","ss","ms","f","p")))
  
  msb = stats["id","ms"] # Mean Square Between
  msc = stats["rater", "ms"] # Mean Square cluster
  mse = stats["resid", "ms"] # Mean Square Error
  msw = (stats["rater","ss"]+stats["resid","ss"])/(stats["rater","df"]+stats["resid","df"]) # Mean Square Within
  
  # Get ICCs #
  icc1 <- (msb - mse)/(msb + (n - 1)*mse)
  icc2 <- (msb - mse)/(msb + (n - 1)*mse + n*(msc - mse)/m)
  icc3 <- (msb - msw)/(msb + (n - 1)*msw)
  
  # Produce Output #
  output = matrix(NA, nrow=3, ncol=1,
                  dimnames=list(c("fixed.raters_rate.all",
                                  "random.raters_rate.all",
                                  "random.raters_rate.one"), 
                                c("icc")))
  output[,"icc"] = c(icc1, icc2, icc3)

  return(output)
}
#############################################################




get.rev.diag = function(m, i = 0) {
  m[ (row(m) + col(m) - 1) %% ncol(m) == i ]
}

fill.rev.diag = function(m, i = 0, value) { 
  m.ix <- matrix(seq_along(m), nrow(m))
  replace(m, get.rev.diag(m.ix, i), value)
}

offset.rev.diag <- function(x, nrow, ncol, offset=1) {
  diag(x, nrow, ncol)[rev(1:nrow), c((offset+1):ncol, 1:offset)]
}


############################
### G & Phi Coefficients ###

### Two-Facet Crossed Design ###
gphi = function(xdat, r=2, t=2, mr=r+3, mt=t+3) {
  
  p = nrow(xdat) # p for person
  id = factor(rep(1:p, r*t))
  score = stack(xdat)[,1]
  rater = factor(rep(rep(1:r, each=p), t))
  time = factor(rep(1:t, each=r*p))
  xdat = data.frame(score, id, rater, time) 
  
  mod = summary(aov(score ~ id*rater*time, data=xdat))
  modmat = matrix(unlist(mod), ncol=3)
  
  y = modmat[,3] # mean squared variance
  n = modmat[1:3,1]+1 # n contains n_person, n_rater, n_time
  
  ### Coefficient Matrix ###
  mat = matrix(0, ncol=length(y), nrow=length(y))
  # fill the first quadrant
  temp = matrix(rep(n[3:1],3), ncol=3, byrow=T) 
  temp = fill.rev.diag(temp, 0, 0)
  mat[1:3, 4:6] = temp
  # fill the second quadrant
  mat[1,1] = n[2]*n[3]
  mat[2,2] = n[1]*n[3]
  mat[3,3] = n[1]*n[2]
  # fill the fourth quadrant
  diag(mat[4:6,4:6]) = n[3:1]
  # fill the ones
  mat[,7] = 1
  
  ### Solve for the Variance ###
  var = solve(t(mat)%*%mat) %*% t(mat) %*% y
  
  ### Get Relative & Absolute Variance ###
  #varrel = sum(var[c(4,5,7)])
  #varabs = sum(var[-1])
  ### Get G & Phi Coefficients ###
  #gcoef =  var[1]/(var[1]+varrel)
  #phicoef =  var[1]/(var[1]+varabs)
  
  ### Dicision Step ###
  gmat = matrix(NA, nrow=mr, ncol=mt, dimnames=list(1:mr, 1:mt))
  phimat = matrix(NA, nrow=mr, ncol=mt, dimnames=list(1:mr, 1:mt))
  
  for (nr in 1:mr) {
    for (no in 1:mt) {
      nopt = c(p, nr, no, nr, no, nr*no, nr*no)
      
      rel = sum (var[c(4,5,7)] / nopt[c(4,5,7)])
      abs = sum (var[-1] / nopt[-1])
      g =  var[1]/(var[1]+rel)
      phi =  var[1]/(var[1]+abs)
      
      gmat[nr,no] = g
      phimat[nr,no] = phi
    }
  }
  
  ### Get G & Phi Coefficients ###
  gcoef =  gmat[r,t]
  phicoef =  phimat[r,t]
  
  list(g = gcoef, phi = phicoef, gmat=gmat, phimat=phimat)
  
}
###########################################################



