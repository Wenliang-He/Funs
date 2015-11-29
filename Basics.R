#################################
### Wenliang's Basics Backage ### Revised Email Version of Stats202.R
#################################

### Dimension ###
size = function(x) {
  if (!any(is.vector(x), is.matrix(x), is.data.frame(x))) stop ("x must be a vector, matrix, or dataframe")
  if (is.null(dim(x))) { size = c(length(x),1)  } else { size = dim(x) }
  return(size)
}

### Vector Norm ###
vecnorm = function(xvec, drop=T) { 
  out = sqrt( t(xvec) %*% xvec ); if (drop) { out = c(out) }
  return(out)
}

### Vector Inner Product ###
inner = function(xvec, yvec=NULL, drop=T, center=F) {
  if (is.null(yvec)) { yvec = xvec }
  if (center) { # if mean-centering is required
    xvec = xvec - mean(xvec, na.rm=T)
    yvec = yvec - mean(yvec, na.rm=T)
  } 
  out = t(xvec) %*% yvec; if (drop) {out = c(out)}
  return(out)
}

### Matrix Outer Product ###
outer.mat = function(xvecmat, yvecmat=NULL, fun, drop=T) {
  # perform row-wise pair-wise function-defined operation
  if (!is.character(fun)) stop('Correct function format is fun="fun.name"')
  if (is.null(yvecmat)) {yvecmat = xvecmat} # if yvecmat is not provided
  if (!any(is.vector(xvecmat), is.matrix(xvecmat))) { xvecmat = as.matrix(xvecmat) } # convert dataframe to matrix
  if (!any(is.vector(yvecmat), is.matrix(yvecmat))) { yvecmat = as.matrix(yvecmat) } # convert dataframe to matrix
  if (is.null(dim(xvecmat))) { xvecmat = t(xvecmat) } # convert column vector to a row vector
  if (is.null(dim(yvecmat))) { yvecmat = t(yvecmat) } # convert column vector to a row vector
  
  out = outer(1:nrow(xvecmat), 1:nrow(yvecmat), 
              FUN=Vectorize(function(x,y)  eval(parse(text=paste0(fun,"(xvecmat[x,], yvecmat[y,])")))))
  if (drop & min(size(out))==1) { out = c(out) } # if output is a vector
  return(out)
}


newton.Binom = function(xmat, yvec, digits=4, miter=100, small=1e-10) {
  sigmoid = function(x) {1/(1+exp(-x))} # Define sigmoid function
  
  m = nrow(xmat); n = ncol(xmat)
  # theta = rep(0, n) # initialize parameters
  y.logit = log((yvec + 0.5) / (1 - yvec + 0.5)) # initialize parameters
  theta = solve(t(xmat) %*% xmat) %*% t(xmat) %*% y.logit # starting value is important 
  L.old = 0
  
  for (iter in 1:miter) {
    z = xmat%*%theta # m_by_1
    L.vec = (yvec*log(sigmoid(z)) + (1-yvec)*log(1-sigmoid(z))) # Log-likelihood for each case
    L = sum(L.vec) # Get Log-likelihood
    
    if (abs(L-L.old) < small) break
    L.old = L
    
    S = (t(xmat)%*%(yvec-sigmoid(z))) # n_by_1
    H = -(t(xmat) %*% diag(as.vector(sigmoid(z)*(1-sigmoid(z)))) %*% xmat)
    theta = theta - solve(H)%*%S
  }
  
  Resid.Dev.vec = -2*L.vec # Residual Deviance for each case
  Dev.resid = sign(yvec - sigmoid(z)) * sqrt(Resid.Dev.vec)
  Dev.resid.sum = quantile(Dev.resid)
  
  se = sqrt(diag(solve(-H))) # -H is observed Fisher Information
  coef = matrix(NA, ncol=4, nrow=n)
  colnames(coef) = c("Estimate", "Std.Err", "z value", "Pr(>|z|)")
  coef[,"Estimate"] = theta
  coef[,"Std.Err"] = se
  coef[,"z value"] = theta/se
  coef[,"Pr(>|z|)"] = 2*(1-pnorm(abs(theta/se)))
  
  colnames = c("Resid.Dev", "Df", "AIC", "BIC", "BIC.exact")
  mat = matrix(NA, ncol=length(colnames), nrow=1)
  colnames(mat) = colnames
  mat[,"Resid.Dev"] = -2*L # sum(Resid.Dev.vec)
  mat[,"Df"] = m - n
  mat[,"AIC"] =  2*n + mat[,"Resid.Dev"] # AIC = 2*n-2*log-likelihood
  mat[,"BIC"] = log(m)*n + mat[,"Resid.Dev"] # BIC = log(m)*n - 2*log-likelihood
  mat[,"BIC.exact"] = log(m)*n - log(2*pi)*n + mat[,"Resid.Dev"] # BIC==log(m)*n - log(2*pi)*n - 2*log-likelihood
  
  list(Deviance.Residuals=round(Dev.resid.sum,digits),
       coef=round(coef,digits), mat=round(mat,digits), iter=iter)
}

# Decompose A Factor Variable into Dummy Variables #
ftodummy = function(fvec, ref=NULL) {
  # Input: A factor variable to be converted into a set of dummy variables
  # Output: All levels besides baseline become dummary variables
  
  # Change Baseline # if baseline reference level is provided
  if (!is.null(ref)) { fvec = relevel(fvec, ref=ref) }
  # Create Dummy Variables #
  lvls = levels(fvec)[-1] # remove baseline
  n = length(lvls) # number of levels to be converted
  xmat = matrix(NA, ncol=n, nrow=length(fvec))
  xmat = sapply(1:n, function(x) xmat[,x]=as.numeric(fvec==lvls[x]))
  colnames(xmat) = lvls # colnames will be replaced if function continues
  return(xmat)
}


genXmat = function(xdatvec, ref=NULL, interact=NULL, combine=NULL) {
  # xdatvec: a dataframe/vector that could contain factors to be converted into a set of dummy variables
  # ref: # need improvement; can only relevel one factor for now
  # interact: all interactions
  # combine: levels of factor(s) to be merged, e.g. combine = c("SS09/SS10/SS12/SS13", "SS12/SS13")
  # Output: Xmat ready for OLS estimation
  
  ### Construct Main Effects Matrix ###
  if ( is.null(dim(xdatvec)) ) { # if xdatvec is a vector
    if (is.factor(xdatvec)) { # if xvec is a factor 
      xmat = ftodummy(xdatvec, ref=ref)
    } else {
      xmat = as.matrix(xdatvec)
    }
  } else { # if xdatvec is a dataframe
    m = nrow(xdatvec); n = ncol(xdatvec)
    fidx = seq(n)[sapply(xdatvec, is.factor)] # factor index
    lvls = sapply(seq(fidx), function(x) levels(xdatvec[,fidx[x]]))
    fmat = matrix(NA, ncol=1, nrow=m) # matrix to accumulate dummy variables from all factors
    for (i in seq(fidx)) {
      tmat = ftodummy(xdatvec[,fidx[i]])
      fmat = cbind(fmat, tmat)
    }
    fmat = fmat[,-1] # remove the first NAs column
    xmat = cbind(fmat, xdatvec[,-fidx]) # matrix with main effects
  }
  
  ### Add More Dummy Variables Based on Merged Levels ###
  n1 = ncol(xmat) # Remember # of columns from main effects matrix
  # because added dummy variables will be removed later!
  if ( !is.null(combine) ) { # if combine is provided
    names = colnames(xmat)
    options(warn=-1)
    for (i in 1:length(combine)){
      input = unlist(strsplit(combine[i], split="/")) # levels to be combined ("lbc")
      # find out which factor do "lbc" belong to
      if ( is.list(lvls) ) { # if multiple factors are present
        wf = seq(lvls)[sapply(lvls, function(x) input %in% x)[1,]] # which factor
        ns = setdiff(lvls[[fidx[wf]]], input) # not selected; levels not to be combined ("lnbc")
      } else { # if only one factor is present
        ns = setdiff(lvls, input)
      }
      newlabels = c(rep(1, length(input)), rep(0, length(ns))) # mark "lbc" as 1 & "lnbc" as 0
      # Now generate a factor variable with "lbc" as 1 & "lnbc" as 0
      vec = factor(xdatvec[,fidx[wf]], levels=c(input,ns), labels=as.character(newlabels)) # cause warning if warning not turned off
      vec = as.numeric(as.character(vec)) # convert to numeric
      xmat = cbind(xmat, vec)
    }
    options(warn=0)
    colnames(xmat) = c(names, combine)
  }
  n2 = ncol(xmat)
  
  ### Add Interactions ###
  # Parse Interactions #
  # find all levels
  if ( is.null(combine) ) { # if combine is not provided
    alvls = unlist(lvls) # all levels
  } else { # if combine is provided
    alvls = c(unlist(lvls), combine)
  }
  # Add Interactions
  if ( !is.null(interact) ) { # if interactions are provided
    for (i in seq(interact)) {
      subcomp = unlist(strsplit(interact[i], split="*", fixed=T)) # main effects
      if ( !all(subcomp %in% alvls) ) { # if some sub-component levels are not in all-levels
        stop("Please check xdatvec & combine. They must include all main effects and their levels") }
      tvec = xmat[,subcomp[1]] * xmat[,subcomp[2]]
      names = colnames(xmat)
      xmat = cbind(xmat, tvec)
      colnames(xmat) = c(names, interact[i])
    }
  }
  
  # Add The Intercept Column #
  #xmat = xmat[,c(1:n1, (n2+1):ncol(xmat))] # remove merged dummy variables in the middle
  #Intercept = rep(1, m)
  #Xmat = cbind(Intercept, xmat)
  
  return(xmat) # return(Xmat)
}

get.interactions = function(xdat, varnames, nway=2) {
  # varnames := a list of colnames from xdat
  # can only deal with 2-way interactions for now
  # cannot deal with factors
  names = combn(varnames, nway)
  newdat = with(xdat, apply(names, 2, function(x) get(x[1])*get(x[2])))
  # apply() automatically converts objects into matrices
  colnames = paste(names[1,], names[2,], sep=".")
  colnames(newdat) = colnames
  newdat = data.frame(xdat, newdat)
  return(newdat)
}

lm.me = function(y, xmat, deci=3) {
  # Remove Missing Rows # 
  tmat = na.omit(cbind(y,xmat)) # only keep rows with complete data
  y = tmat[,1]; xmat = tmat[,-1]
  
  # Initialize #
  df = dim(xmat)[1] - dim(xmat)[2]
  qr = qr(xmat)
  R = qr.R(qr)
  xxi = solve(t(R)%*%R)
  # Prepare Output Matrix #
  outmat = matrix(NA, nrow=dim(xmat)[2], ncol=4)
  rownames(outmat) = colnames(xmat)
  colnames(outmat) = c("Estimate", "Std.Err", "t value", "Pr(>|t|)")
  # Produce Output #
  outmat[,"Estimate"] = qr.solve(qr, y)
  outmat[,"Std.Err"] = sqrt(diag(var(y)*xxi))
  outmat[,"t value"] = outmat[,"Estimate"]/outmat[,"Std.Err"]
  outmat[,"Pr(>|t|)"] = pt(abs(outmat[,"t value"]), df=df, lower.tail=F)*2
  
  return(round(outmat, deci))
}

########################################
### Data Inspection and Manipulation ### Six Sets of Functions
########################################
### 1. Overview ###
overview = function(xdat, long=TRUE, nmax=100, nlvl=5) {
  # xdat must be a dataframe 
  # nlvl is the number of levels to display for factors with numerous (i.e. nmax) levels
  dim = dim(xdat); head = head(xdat); class = sapply(xdat, class) # variable class
  fidx = c(1:ncol(xdat))[class == "factor"] # factor index
  nfidx = c(1:ncol(xdat))[class != "factor"] # non-factor index
  nf = length(fidx) # number of factors in xdat
  nnf = dim[2] - nf # number of non-factors in xdat
  values = rep(NA, length=dim[2]) # unique values in each variable
  miss = rep(NA, length=dim[2]) # missing values in each variable
  miss = sapply(xdat, function(x) sum(is.na(x)))
  
  if(nf==0) { # no factors in xdat
    levels = lapply(xdat[,fidx], levels) # fill levels as empty
    values = sapply(xdat, function(x) length(unique(x))) # find & fill values
  }
  if (nf==1) { # only one factor in xdat
    if (long) { 
      levels = table(xdat[,fidx], useNA="ifany") # find levels
      values[fidx] = paste(length(levels)-(miss[fidx]!=0)*T, "lvls") # find & fill values # discount NA as a level
    } else { 
      levels = levels(xdat[,fidx])
      values[fidx] = paste(length(levels), "lvls") # find & fill values
    }
    values[nfidx] = sapply(1:nnf, function(x) length(unique(xdat[,nfidx[x]])))
    if (length(levels)>nmax) {levels=levels[1:min(length(levels),nlvl)]} # fill levels
  }
  if (nf>1) { # multiple factors in xdat
    if (long) { # if long = TRUE
      levels = lapply(xdat[,fidx], function(x) table(x, useNA="ifany")) # find levels 
      values[fidx] = paste(sapply(levels, length)-(miss[fidx]!=0)*T, "lvls") # find & fill values # discount NA as a level
    } else { # if long = FALSE
      levels = lapply(xdat[,fidx], levels) # find levels
      values[fidx] = paste(sapply(levels, length), "lvls") # find & fill values
    }
    if (nnf>0) { # only handles the case when # of non-factors > 0
      values[nfidx] = sapply(1:nnf, function(x) length(unique(xdat[,nfidx[x]])))
    }
    if (any(sapply(levels, length)>nmax)) {levels = lapply(levels, function(x) x[1:(min(length(x),nlvl))])} # fill levels
  }
  
  class = rbind(class, values, miss)
  rownames(class) = c("class", "unique.values","missing")
  list = list(dim, head, class, levels)
  names(list) = c("dim", "head", "class", "levels")
  return(list)
}

### 2. Fix Minor Problems ###
minorfix = function(xdat, deepfix=FALSE) {
  typesofFix = c(" Convert \"\" level in factors into NA.\n",
                 "Remove empty levels.\n",
                 "Remove completely NA rows.\n",
                 "Use factors' second level names as column names.\n")
  class = sapply(xdat, class)
  fidx = c(1:ncol(xdat))[class == "factor"] # factor index
  
  if (length(fidx)==1) {
    levels = levels(xdat[,fidx])
    # Type 1 #
    type = any(levels=="") # Contain "" level in a factor
    if (type) {
      xdat[,fidx] = factor(xdat[,fidx], levels=levels[-1])
      typeidx = TRUE
    } else { typeidx = FALSE }
    # Type 2 #
    type = any(table(xdat[,fidx])==0)
    if (type) {
      xdat[,fidx] =  factor(xdat[,fidx])
      typeidx = TRUE
    } else { typeidx = FALSE }
  } else { 
    levels = lapply(xdat[,fidx], levels)
    # Type 1 #
    type = sapply(levels, function(x) any(unlist(x)=="")) # Contain "" level in factors
    if (any(type)) { 
      for (i in fidx[type]) { xdat[,i] =  factor(xdat[,i], levels=levels[[colnames(xdat)[i]]][-1]) }
      typeidx = TRUE
    } else { typeidx = FALSE }
    # Type 2 #
    type = sapply(xdat[,fidx], function(x) any(table(x)==0)) # Contain empty levels
    if (any(type)) { 
      for (i in fidx[type]) { xdat[,i] =  factor(xdat[,i]) }
      typeidx = c(typeidx, TRUE)
    } else { typeidx = c(typeidx, FALSE) }
  }
  
  # Type 3 #
  type = apply(xdat, 1, function(x) all(is.na(x))) # Contain completely NA rows
  if (any(type)) { 
    xdat = na.omit(xdat) # xdat = xdat[!type,]
    typeidx = c(typeidx, TRUE)
  } else { typeidx = c(typeidx, FALSE) }
  
  if (deepfix) {
    # Type 4 #
    tlidx = c(1:length(fidx))[sapply(levels, length)==2] # two-level factor indices
    lvlname = tolower(sapply(levels, function(x) x[2])[tlidx])
    colname = tolower(colnames(xdat)[fidx[tlidx]])
    type = any(lvlname!=colname) # The second level of a factor variable does not match its colname 
    if (type) { 
      for (i in 1:length(lvlname)) {
        xdat = rename(colname[i], lvlname[i], xdat)
      }
      typeidx = c(typeidx, TRUE)
    } else { typeidx = c(typeidx, FALSE) }
    
    # Report Types-of-Fix #
    if( any(typeidx) ) { cat(typesofFix[typeidx])
    } else { cat("There is nothing to fix; Nothing is modified.") }
  }
  
  # Report Types-of-Fix #
  n = length(typeidx)
  if ( any(typeidx) ) { cat( typesofFix[ c(1:n)[typeidx] ] )
  } else { cat("There is nothing to fix; Nothing is modified.") }
  
  return(xdat)
}

#### 3. Rename & Locate ###
rename = function(oldnames, newnames, dataframe){
  # oldnames: a single or vector of indices(integers) or characters 
  # newnames: a single or vector of characters
  n = length(oldnames)
  for (i in 1:n) {
    if (is.character(oldnames[i])) { colnames(dataframe)[colnames(dataframe)==oldnames[i]] = newnames[i]
    } else { colnames(dataframe)[oldnames[i]] = newnames[i] }
  }
  return(dataframe)
}

wordidx=function(words, xvecdat) { # called by tabmod()
  # Can handle multiple words input
  if (!is.vector(words)) stop("First input must be a vector of words to be found")
  if (!is.null(dim(xvecdat))) {xvecdat = colnames(xvecdat)} # If dataframe or matrix
  index = sapply(seq(words), function(x) seq(xvecdat)[xvecdat %in% words[x]])
  #if (is.null(dim(xvecdat))) { index = seq(xvecdat)[xvecdat %in% words] } # if vector
  #else { index = seq(ncol(xvecdat))[colnames(xvecdat) %in% words] } # if dataframe or matrix
  return(index)
}

# test unit: xvec = round(runif(20, 1, 6)) #
revlikert = function(xvec, oldorder=1:6, neworder=6:1) {
  names(neworder) = oldorder
  newxvec = neworder[xvec]
  names(newxvec) = NULL
  return(newxvec)
}

dupid = function(xvecdat, id="id"){
  # allow input to be vector or matrix/dataframe #
  if(is.numeric(dim(xvecdat))) {xvec = xvecdat[,id]} else {xvec=xvecdat}
  # find duplicate ID #
  if (length(xvec)==length(unique(xvec))) 
  {out="No Duplicate ID"} else {
    tab=table(xvec)
    tab=tab[tab>1]
    out=matrix(NA, nrow=length(tab), ncol=2)
    colnames(out)=c("id", "times")
    out[,"id"]=as.numeric(names(tab))
    out[,"times"]=tab
    out = out[order(out[,"times"]),]
  }
  return(as.data.frame(out))
}

dupitem = function(xdat, id="id"){ # Bug! Need Fix!
  idtab = dupid(xdat[,id])
  idlist = idtab[,"id"]
  
  item = data.frame()
  for (i in 1:length(idlist)) {
    temp = xdat[xdat[,id]==idlist[i],]
    item = rbind(item, temp)
  }
  #item = xdat[is.element(xdat[,id], idlist), ]
  list(idtab=idtab, item=item)
}

### 4. Combine Multiple Levels/Labels ###
### Compare Labels/Levels ###
cplabels = function(labelset, labelspace, split=";") { # Old-Set, New-Space
  xvec = labelspace
  xvec = paste(xvec, collapse=split)
  xvec = strsplit(xvec, split=split)[[1]]
  xvec = gsub("\\n", "", xvec) # remove all \n
  xvec = gsub("^ *", "", xvec) # remove all whitespaces in the front
  xvec = gsub(" *$", "", xvec) # remove all whitespaces in the back
  inOldButNotNew = setdiff(labelset, xvec)
  inNewButNotOld = setdiff(xvec, labelset)
  
  if (length(inNewButNotOld)==0 & length(inOldButNotNew)==0) {
    cat("The two labels contain the same elements! \n")
    return(output <- NULL)
  } else {
    output = list(inSetButNotSpace = sort(inOldButNotNew), inSpaceButNotSet=sort(inNewButNotOld))
    return(output)
  }
}

### Combine Labels/Levels ###
relabel = function(variable, oldlabels, newlabels, split=";") {
  
  labels = lapply(oldlabels, function(x) unlist(strsplit(x, split=split)))
  labels = lapply(labels, function(x) gsub("\n","", x))   # remove all \n 
  labels = lapply(labels, function(x) gsub("^ *", "", x)) # remove all whitespaces in the front
  labels = lapply(labels, function(x) gsub(" *$", "", x)) # remove all whitespaces in the back
  
  newvec = vector()
  for (i in 1:length(labels)) {
    newvec[ variable %in% labels[[i]] ] = newlabels[i]
  }
  newvec = factor(newvec, levels=newlabels)
  
  return(newvec)
}

### 5. Reorder Factor Levels ###
genfcode = function(xdat) { # generate "factor coding scheme" from a well-cleansed dataframe
  n = ncol(xdat); class = sapply(xdat, class); fidx = seq(n)[class == "factor"] 
  if ( length(fidx)==0 ) { stop("There are no factors in the dataframe!") }
  yesno = vector(); newlevels = vector()
  for (i in 1:length(fidx)) {
    newlvl = levels(xdat[,fidx[i]]) # target order of levels
    oldlvl = sort(newlvl) # if saved as csv and read back into R
    if ( all(newlvl==oldlvl) ) { yesno[i] = "no"; newlevels[i] = NA }
    else {
      yesno[i] = "yes"
      seq = 1:length(newlvl) # contents of a dictionary
      names(seq) = oldlvl # entries of a dictionary
      order = seq[newlvl] # get contents by calling the entries
      newlevels[i] = paste(order, collapse=";")
    }
  }
  if ( all(yesno=="no") ) { cat("Just store/read it; You do not need a coding scheme!\n") }
  else { list(yesno=yesno, newlevels=newlevels) }
}

transfcode = function(yesno, newlevels) { # translate "factor coding scheme"
  # Completely subsumed under refactor()
  yesno = gsub("\\\n", "", yesno) # remove \n
  yesno = gsub("\\\"", "", tolower(yesno)) # remove \"
  yesno = strsplit(yesno, split=" ")[[1]]
  yesno = yesno[yesno!=""]
  
  newlevels = gsub("\\\n", "", newlevels) # remove \n
  newlevels = gsub("\\\"", "", newlevels) # remove \"
  newlevels = strsplit(newlevels, split=" ")[[1]]
  newlevels = newlevels[newlevels!=""]
  
  list(yesno = yesno, newlevels=newlevels)
}


refactor = function(xdat, yesno=NULL, newlevels=NULL, byhand=FALSE) {
  if ( is.null(yesno) && !byhand ) { stop('Please either provide the coding scheme or set "byhand=T"!\n') }
  
  if ( !is.null(yesno) ) { code = transfcode(yesno, newlevels) } # if coding scheme is provided
  dim = dim(xdat)
  class = sapply(xdat, class)
  fidx = c(1:ncol(xdat))[class == "factor"] # factor index
  if ( length(fidx)==0 ) { stop("There are no factors in the dataframe!") }
  fname = colnames(xdat)[fidx]
  levels = lapply(xdat[,fidx], levels)
  tab = lapply(xdat[,fidx], table)
  
  if (byhand) { # If user wants to input by hand, espcially for the first time
    yesno = vector()
    newlvls = vector()
    for (i in 1:length(fidx)) {
      #print(tab[[i]])
      print(rbind(1:length(tab[[i]]),tab[[i]]))
      cat(paste("Do you want to refactor", fname[i]), "?", sep="")
      yesno[i] = readline('Type "yes" or "no": ')
      if (tolower(yesno[i]) == "yes") {
        newlvls[i] = readline('Input new levels (or indices of levels) in order: \nPlease use ";" as the separator, e.g. 1;4;2;3 ')
        
        input = unlist(strsplit(newlvls[i], split=";"))
        if (length(unlist(sapply(letters, function(x) grep(x, input[1:min(5,length(input))]))))==0) { # if indices/integers are provided
          input = as.numeric(input)
          xdat[,fname[i]] = factor(xdat[,fname[i]], levels=levels[[i]][input])
        } else { # if labels/characters are provided
          xdat[,fname[i]] = factor(xdat[,fname[i]], levels=input)
        }
      }
    }
    
  list(code=list(yesno=yesno, newlevels=newlvls), xdat=xdat)
    
  } else { # Already have the refactor_code
    yesno = code[[1]]; newlvls = code[[2]]
    for(i in 1:length(yesno)){
      if ( yesno[i]=="yes" ) { 
        input = unlist(strsplit(newlvls[i], split=";"))
        if (length(unlist(sapply(letters, function(x) grep(x, input[1:min(5,length(input))]))))==0) { # if indices/integers are provided
          input = as.numeric(input)
          xdat[,fname[i]] = factor(xdat[,fname[i]], levels=levels[[i]][input])
        } else { # if labels/characters are provided
          xdat[,fname[i]] = factor(xdat[,fname[i]], levels=input)
        }
      }
    }
    return(xdat)
  }
}

### 6. Impute Missing Data ###
impute.simpleOLS = function(mvec, xvecmat, likert=6, truncate=FALSE, trim=FALSE) { 
  # called by impute.OLS(); call genXmat()
  # mvec must contain missing data and can be continous or on likert scale 
  # xvecmat must NOT contain missing data
  if (!any(is.na(mvec))) stop("missing data 'mvec' does not contain any missing data!")
  if (any(is.na(xvecmat))) stop("observed data 'xvecmat' must NOT contain missing data!")
  xvecmat = genXmat(xvecmat) # convert factors into numeric variables for lsfit() to process
  n = length(mvec); midx = seq(n)[is.na(mvec)] # missing data index
  m = length(midx) # number of missing data
  mod = lsfit(xvecmat[-midx,], mvec[-midx]) # note: lsfit cannot process factors
  df = length(mod$coef)
  rse = sqrt( sum(mod$residuals^2) / (n-m-df) ) # residual standard deviation
  cmean = as.vector(as.matrix(cbind(1,xvecmat[midx,])) %*% mod$coef)
  cdraw = cmean + rse*rnorm(length(cmean)) # prediction using conditional draw
  if (trim) { # trim continuous predicted values to be on proper likert scale
    mvec[midx] = as.numeric(cut(cdraw, breaks=c(-100, seq(from=1.5, to=(likert-0.5), by=1), 100)))
  } else {
    mvec[midx] = cdraw
  }
  if (truncate) { # truncate continuous predicted values within bounds of [1, likert]
    mvec[mvec>likert & !is.na(mvec)] = likert # upper bound
    mvec[mvec<1 & !is.na(mvec)] = 1 # lower bound
  }
  return(mvec)
}

impute.OLS = function(mvecmat, xvecmat, method=NULL, likert=6, truncate=FALSE, trim=FALSE) {
  # call impute.simpleOLS()
  # mvecmat contains missing data, whose columns must be in the right order from left to right for imputation
  # xvecmat contains observed data; it's OK to have some missing data in xvecmat (casewise deletion is applied)
  # method = 1-impute separately; 2-impute sequentially; 3-impute the average
  mdat = as.data.frame(mvecmat); nm = ncol(mdat)
  xdat = as.data.frame(xvecmat); nx = ncol(xdat)
  if (!any(is.na(mdat))) stop("mvecmat does not contain any missing data!")
  if (is.null(method)) { method = as.numeric(readline('Please specify method: \n1-impute separately; 2-sequentially; 3-average. Decision: ')) }
  if (nrow(mdat)==nrow(xdat)) {m = nrow(xdat)} else { stop("mvecmat & xvecmat have different # of cases") }
  mrow = seq(m)[apply(xdat, 1, function(x) any(is.na(x)))] # missing row index of xdat
  
  ### Fill in Missing Data ###
  if (method == 1) { # Impute Items Separately
    for (i in seq(nm)) {
      if (length(mrow)==0) { # If xdat contains no missing data
        mdat[,i] = impute.simpleOLS(mvec=mdat[,i], xvecmat=xdat, likert=likert, truncate=truncate, trim=trim)
      } else { # If some rows of xdat contain missing data
        mdat[-mrow,i] = impute.simpleOLS(mvec=mdat[-mrow,i], xvecmat=xdat[-mrow,], likert=likert, truncate=truncate, trim=trim) } } }
  if (method == 2) { # Impute Items Sequentially
    for (i in seq(nm)) {
      if (length(mrow)==0) { # If xdat contains no missing data
        mdat[,i] = impute.simpleOLS(mvec=mdat[,i], xvecmat=xdat, likert=likert, truncate=truncate, trim=trim)
      } else { # If some rows of xdat contain missing data
        mdat[-mrow,i] = impute.simpleOLS(mvec=mdat[-mrow,i], xvecmat=xdat[-mrow,], likert=likert, truncate=truncate, trim=trim)
      }
      # Differ from method==1 by this one line #
      xdat = cbind(xdat, mdat[,i]) } }
  if (method == 3) { # Impute Items Average #
    mdat = simpute(mvec=rowMeans(mdat, na.rm=T), xvecmat=xdat, likert=likert, truncate=truncate, trim=trim) }
  if (nm == 1) { mdat = mdat[,1] } # if mdat is a column of a matrix, convert it back to a vector!
  return(mdat)
}
################################################################################



############################################
### Get Summary & Descriptive Statistics ### One Set of Functions
############################################
### 1. Vector Summary: Handle Two-Level Factors Only ###
summary.vec = function(xvec, confid=0.95, digits=3, long=FALSE) { # Completely subsumed under getsummary()
  if (factor <- is.factor(xvec)) {xvec = as.numeric(xvec)-1} # convert two-level factors into 0s and 1s
  alpha = 1-confid
  # Sample Statistics #
  total = length(xvec); missing = sum(is.na(xvec)); n = total-missing # non-missing cases
  t = qt(1-alpha/2, n-1) # student t distribution
  mean = mean(xvec, na.rm=T); median = median(xvec, na.rm=T)
  if (factor) {
    sampleSize = sum(xvec==1, na.rm=T) # Count
    sd = sqrt(mean*(1-mean))
  } else {
    sampleSize = n # Count
    sd = sd(xvec, na.rm=T)
  }
  min = min(xvec, na.rm=T); max = max(xvec, na.rm=T)
  pskew = sum((xvec-mean)^3, na.rm=T) /n /sd^3 /(1-1/n)^1.5 # population skewness
  skew = sqrt(n*(n-1))/(n-2) * pskew  # sample skewness 
  pkurtosis = sum((xvec-mean)^4, na.rm=T) /n /sd^4 - 3 # population excess skewness
  kurtosis = (n-1)/((n-2)*(n-3)) * ((n+1)*pkurtosis + 6) # sample excess kurtosis
  # http://www.tc3.edu/instruct/sbrown/stat/shape.htm#Skewness
  se = sd/sqrt(sampleSize) # large sample approximation by CLT
  ci = mean + c(-1,1)*t*se # large sample approximation by CLT
  
  if (long) {
    tab = matrix(c(total, missing, sampleSize, mean, sd, min, median, max, skew, kurtosis, se, ci), nrow=1, byrow=F, 
                 dimnames=list("xvec", c("Total","Missing","Count","Mean","SD","Min","Median","Max",
                                         "Skewness","Kurtosis","Std.Err.","[95% Conf.","Interval]")))
    } else {
    tab = matrix(c(total, missing, sampleSize, mean, sd, min, median, max), nrow=1, byrow=F, 
                 dimnames=list("xvec", c("Total","Missing","Count","Mean","SD","Min","Median","Max")))
    }
  return(round(tab, digits))
}

### 2. Generate Two-level Factors from a Multi-level Factor ###
genfactor = function(fvec) { # called by getsummary() & 
  # Input: A factor of n levels
  # Ouput: A dataframe with n columns; each is a two-level factor
  lvls = levels(fvec); nlvls = length(lvls)
  temp = fvec
  for (i in 1:nlvls) {
    assign( paste("var",i,sep=""), relabel(fvec, c(paste(lvls[-i], collapse=";"),lvls[i]), c("Others",lvls[i])) )
    temp = data.frame(temp, get(paste("var",i,sep="")))
  }
  dat = temp[,-1]; colnames(dat) = lvls
  return(dat)
}

### 3. Any Summary: Vectors, Matrices, DataFrames ###
getsummary = function(xvecmat, confid=0.95, digits=3, long=FALSE) {
  xvecmat = as.data.frame(xvecmat)
  # Call summary.vec() & genfactor()
  if (ncol(xvecmat)==1) { # if xvecmat is a vector
    tab = summary.vec(xvecmat[,1], confid=confid, digits=digits, long=long)
  } else { # if xvecmat is a matric or data.frame
    class = sapply(xvecmat, class)
    fidx = c(1:ncol(xvecmat))[class == "factor"] # factor index
    names = colnames(xvecmat)[class!="factor"] # names of continuous variables
    if (length(fidx)!=0) { # when factors are present
      tdat = xvecmat[,-fidx] # base cummulator with all non-factor variables
      for (i in 1:length(fidx)) { tdat = data.frame(tdat, genfactor(xvecmat[,fidx[i]])) }
      xvecmat = tdat }
    xvecmat = rename(seq(length(names)), names, xvecmat) # give continuous variables their original names
    tab=t(sapply(xvecmat, function(x) summary.vec(x, confid=confid, digits=digits, long=long)))  
  }
  if (long) { colnames(tab) = c("Total","Missing","Count","Mean","SD","Min","Median","Max",
                                "Skewness","Kurtosis","Std.Err.","[95% Conf.","Interval]")
  } else {  colnames(tab) =c("Total","Missing","Count","Mean","SD","Min","Median","Max") }
  return(tab)
}

getsummary.old = function(xvecmat, confid=0.95, digits=3, short=TRUE) {
  # Call summary.vec()
  
  if (is.null(dim(xvecmat))) { # if xvecmat is a vector
    tab = summary.vec(xvecmat)
  } else { # if xvecmat is a matric or data.frame    
    tab=t(sapply(xvecmat, function(x) summary.vec(x, confid=confid, digits=digits, short=short)))
    
    if (short) {
      colnames(tab) =c("Total","Missing","Count","Mean","SD","Min","Median","Max")
    } else {
      colnames(tab) = c("Total","Missing","Count","Mean","SD","Min","Median","Max",
                        "Skewness","Kurtosis","Std.Err.","[95% Conf.","Interval]")
    }
  }
  return(tab)
}
########################################################################################


##########################################
### Comparative Descriptive Statistics ###
##########################################

pwttest = function(xdat, paired=TRUE, deci=3) { # To improve: put pvalue & meandiff in one table
  n=ncol(xdat); m=nrow(xdat); names = colnames(xdat)
  mat = outer(1:n, 1:n, paste) # create pair-wise indices
  idx = mat[upper.tri(mat)] 
  idx = strsplit(idx, split=" ") # idx as a list
  
  pval=vector(); meandiff=vector()
  for(i in 1:length(idx)) {
    cols = as.numeric(idx[[i]])
    t = t.test(xdat[,cols[1]], xdat[,cols[2]], paired=paired) 
    pval[i] = t$p.value
    meandiff[i] = ifelse(paired, t$estimate, -diff(t$estimate)) # the way paired ttest shows mean.diff
  }
  
  pmat = matrix(1.0, n, n); rownames(pmat) = names; colnames(pmat) = names
  diffmat = matrix(0.0, n, n); rownames(diffmat) = names; colnames(diffmat) = names
  pmat[upper.tri(pmat)] = pval # fill pvalue matrix
  pmat[lower.tri(pmat)] = t(pmat)[lower.tri(pmat)]
  diffmat[upper.tri(diffmat)] = meandiff # fil mean difference matrix
  diffmat[lower.tri(diffmat)] = -1*t(diffmat)[lower.tri(diffmat)]
  
  cat(ifelse(paired, "\tPaired ttests Results:\n", "\tUnpaired ttests Results:\n"))
  cat("\tMean.diff = row - column\n")
  list(pvalue=round(pmat,deci), meandiff=round(diffmat,deci))
}

### 1. Vector Comparison ###
# Input #
# xvec = IV: continuous or nominal/ordinal variables
# cat = DV: two-level indicator variable

# Limitations #
# cat can be only two-level; "two"-sample t-tests are used for continuous variables 
# cannot handle ordinal variable

cp.vec = function(xvec, cat, digits=2) { # call summary.vec() & es() & chisqtest()
  lvls = levels(cat) # Should be "Control" & "Treat(ment)"
  if (class(xvec)=="integer" || class(xvec)=="numeric" || length(unique(xvec))>=10) { # If xvec is continuous
    tab = matrix(NA, nrow=1, ncol=4) # create an empty row #
    esvec = vector() # Effect-Size-VECtor to store statistics for computing effect size
    # Fill in Summary Statistics by Level #
    for (i in 1:length(lvls)) {
      vecsum = round(summary.vec(xvec[cat==lvls[i]]), digits)
      tab[,i] = paste0(vecsum[,"Mean"], " (", vecsum[,"SD"], ")") # "±" is "\u00b1"
      # http://www.utf8-chartable.de/unicode-utf8-table.pl
      esvec[3*(i-1)+1:3] = vecsum[,c("Mean","SD","Count")]
    }
    # Fill in Test Statistics #
    tstat = -round(t.test(xvec~cat)$stat, digits) # Add a minus sign s.t. it's Treat - Control
    pval = round(t.test(xvec~cat)$p.value, digits+1)
    tab[,3] = paste(tstat, " (", pval, ")", sep="") # t-statistic (p-value)
    tab[,4] = round(es(esvec[4], esvec[5], esvec[1], esvec[2], esvec[6], esvec[3]), 3) # effect size
  } else {# If xvec is categorical
    tab = matrix(NA, nrow=nrow(table(xvec, cat)), ncol=4) # create a empty matrix #
    # Fill in Summary Statistics #
    freq = table(xvec, cat); percent = round(100 * t(t(freq) / colSums(freq)), digits)
    tab[,1:2] = matrix(paste0(percent, "% (", freq, ")"), ncol=2)
    # Fill in Test Statistics #
    #if (is.ordered=="TRUE") {
    ##### need to figure out rowscr and colscr #########
    #  chisq = round(msqtest(freq, rowscr, colscr)$chisq, deci)
    #  pval = round(msqtest(freq, rowscr, colscr)$p.chisq, deci+1)
    #} else {
    chisqtest = chisqtest(freq)
    chisq = round(chisqtest$chisq, digits)
    pval = round(chisqtest$p.chisq, 3)
    #} 
    tab[1,3] = paste0(chisq, " (", pval, ")") # chi-squared-statistic (p-value)
  }
  return(tab)
}

### 2. Any Comparison ###
cp.table = function(xdat, cat, digits=3) {
  # Get the Statistics #
  stats = lapply(xdat, function(x) cp.vec(x, cat, digits))
  # Put Results in One Table #
  stats.tab = NULL
  for (i in 1:length(stats)) { stats.tab = rbind(stats.tab, stats[[i]]) }
  # Add Column Names #
  colnames(stats.tab) = c(rep("Mean (SD) or Percent (N)",2), "t(p) or Chi-squared(p)", "Effect Size")
  # Add Row Names #
  varnames = colnames(xdat); lvls = lapply(xdat, levels)
  rownames = vector()
  for (i in 1:length(lvls)) {
    if (is.null(lvls[[i]])) { # If the ith variable is continuous and has no levels
      rownames = c(rownames, varnames[i])
    } else {
      rownames = c(rownames, lvls[[i]])
    }
  }
  rownames(stats.tab) = rownames
  
  return(stats.tab)
}
############################################################




###################################
### Plotting Exploratory Graphs ###
###################################

hist.vec = function(xvec, freq=TRUE, label=TRUE, reverse=TRUE, varname=NULL, density=10, nmax=10, ...){
  # reverse: plot-counts-label-percentages | plot-percentages-label-counts 
  # nmax: discrete if < nmax & continuous o.w.
  name = as.character(match.call())[2]
  counts = table(xvec)
  percts = 100*counts/sum(counts)
  
  if (class(xvec)=="factor" || length(unique(xvec))<=nmax) { # if discrete, use barplot()  
    if (freq) { # plot counts
      barplot = barplot(counts, density=density, main=varname, ...)
      if (label) { # include label
        if (reverse) { # label percentages
          text(x=barplot, y=counts+max(counts)*0.05, labels=paste0(round(percts,2),"%"), xpd=TRUE) 
        } else { # label counts
          text(x=barplot, y=counts+max(counts)*0.05, labels=as.character(counts), xpd=TRUE)
        }
      }
    } else { # plot percentages
      barplot = barplot(percts, density=density, main=varname, ...)
      if (label) { # include label
        if (reverse) { # label counts
          text(x=barplot, y=percts+max(percts)*0.05, labels=as.character(counts), xpd=TRUE)
        } 
        else { # label counts
          text(x=barplot, y=percts+max(percts)*0.05, labels=paste0(round(percts,2),"%"), xpd=TRUE)
        }
      }
    }
  } else { # if continuous, use hist()
    hist(xvec, freq=freq, main=varname, ...)
  }
}

hist.me = function(xdat, figname="fig", varnames=NULL, write=FALSE, width=300, height=400, ...) {
  # To improve: handle multiple xdats together
  #if (class(xlsdat)=="list") { iter = length(xlsdat)
  #} else { iter = 1; xlsdat = list(ls=xlsdat) }
  n = ncol(xdat)
  if (is.null(varnames)) { varnames = colnames(xdat) }
  
  if (write) { # if output to hard disk
    for (i in 1:n) {
      png(paste0(figname,".", varnames[i],".png"), width=width, height=height)
      
      hist.vec(xdat[,i], varname=varnames[i], ...)
      
      dev.off()
    }
  } else { # if plot in R
    for (i in 1:n) { hist.vec(xdat[,i], varname=varnames[i], ...) }
  }
}


cp.barplot = function(xdat, beside=TRUE, density=NULL, nmax=6, ...){
  # xdat can include multiple columns, where nmax = maximum # of columns in xdat
  n = ncol(xdat); den = seq(40, 0, by=-8)
  if (n>nmax) { stop(paste0("Too many variables (>", nmax,") in dataframe. Please increase nmax to proceed!")) }
  if(is.null(density)){density=den[1:n]} else {density=density}
  perct = 100*sapply(xdat, function(x) table(x)/sum(table(x)))
  barplot(t(perct), beside=beside, density=density, ...)
}

cp.legend = function(xdat, pos="top", legend=NULL, bty="n", horiz=T, density=NULL, inset=c(0,-0.15), nmax=6, ...) {  
  par(xpd=TRUE)
  n = ncol(xdat); den = seq(40, 0, by=-8)
  if (n>nmax) { stop(paste0("Too many variables (>", nmax,") in dataframe. Please increase nmax to proceed!")) }
  if(is.null(density)){density=den[1:n]} else {density=density}
  if(is.null(legend)){legend=colnames(xdat)}
  legend(pos, legend=legend, density=density, bty=bty, horiz=horiz, inset=inset, ...)
  par(xpd=FALSE)
}

scatplot = function(xdat, pos="topleft", nf=2, xlab=NULL, ylab=NULL, main=NULL, ...){
  # Two-way or Three-way Scatterplots #
  # xdat must include y & x, where both must be continuous 
  # z is optional and z can be continuous or a factor variable of multiple levels
  # nf is the number of levels to split if z is continuous
  if(is.null(ylab)) {ylab=colnames(xdat)[1]}
  if(is.null(xlab)) {xlab=colnames(xdat)[2]}
  if(is.null(main)) {main=paste("Plotted by",colnames(xdat)[3],sep=" ")}
  
  par(xpd=FALSE) # To ensure abline() not plotting outside the plot region
  y = xdat[,1]; x = xdat[,2]
  ry = range(y, na.rm=T); rx = range(x, na.rm=T)
  
  if( ncol(xdat)==2 ) { # if no z is provided
    plot(x, y, main=main, ylab=ylab, xlab=xlab,
         ylim=ry + c(-1,1)*diff(ry)*0.05, ...)
    mod=lm(y ~ x)
    abline(mod)
  } else { # z is provided
    z = xdat[,3]; factor = is.factor(z)
    if (!factor) { # if z is a continuous
      z = cut(z, breaks = quantile(z, seq(0, 1, 1/nf)), include.lowest=T)
      z = factor(z, labels=sapply(1:nf, function(x) paste0("Level",x)))
    }
    lvls = levels(z)
    plot(x[z==lvls[1]], y[z==lvls[1]], 
         main=main, ylab=ylab, xlab=xlab, 
         xlim=rx + c(-1,1)*diff(rx)*0.05, 
         ylim=ry + c(-1,1)*diff(ry)*0.05, ...)
    abline(lm(y[z==lvls[1]] ~ x[z==lvls[1]]))
    
    for (i in 2:length(lvls)){
      points(x[z==lvls[i]], y[z==lvls[i]], col=i)
      abline(lm(y[z==lvls[i]] ~ x[z==lvls[i]]), col=i)
    }
    legend(pos, lvls, bty="n", col=1:length(lvls), lty=1)
    mod = lm(y ~ x*z)
  }
  return(summary(mod))
}


plotbygrp = function(xdat, ylab=NULL, nmin=3, nmax=6) {
  # Input: xdat has two columns: 1st is xvec (num) & 2nd yvec (num or dummy)
  names = colnames(xdat)
  xvec = xdat[,1]; yvec = xdat[,2]; factor = is.factor(yvec)
  if (factor) { if ( length(levels(yvec))!=2 ) # Stop is yvec is cat, but not dummy
  {stop("This function only allows using continuous or dummy variables as yvec")} }
  
  par( mfrow=c(2, ceiling((nmax-nmin+1)/2) ) )
  
  if ( is.null(ylab) ) {
    if (factor) { ylab = paste("Percentage of",levels(yvec)[2])
    } else { ylab = paste("Average of", names[2]) }
  }
  groupmean = list()
  for (i in nmin:nmax) {
    grp = cut(xvec, breaks = quantile(xvec, seq(0, 1, 1/i), na.rm=T), include.lowest = T)
    if (factor) { # if yvec is a factor
      tab = table(grp, yvec); perct = tab / rowSums(tab)
      mean = perct[,2] # only need the second column
    } else { # if yvec is continuous
      mean = tapply(yvec, grp, mean, na.rm=T)
    }
    r = range(mean, na.rm=T)
    plot(1:i, mean, type="o",
         xlab = paste("By", i, "Groups of", names[1]), 
         ylab = ylab,
         ylim=r + c(-1,1)*diff(r)*0.05)
    
    groupmean[[i-nmin+1]] = mean
  }
  
  par(mfrow=c(1,1))
  return(groupmean)
}
######################################################################



#################################
### Barplot for Brian's Paper ### Need to improve clarity!
## http://monkeysuncle.stanford.edu/?p=485
error.bar <- function (x, y, upper, lower=upper, length=0.1,...){
  if (length(x) != length(y) | length(y) != length(lower) | length(lower) != length(upper))
    stop ("vectors must be same length")
  arrows (x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}


barplot.me = function(y, x=NULL, group=NULL, main=NULL, ylab=NULL, ylim=NULL, legend=NULL, names.arg=NULL){
  # y is on x-axis and its mean on y-axis
  if(!is.factor(group)) {group=factor(group)}
  
  if( is.null(dim(y))==TRUE ) {
    mean = tapply(y, group, mean, na.rm=T)
    sd = tapply(y, group, sd, na.rm=T)
    n = tapply(y, group, function(x) sum(!is.na(x)))
    se = sd/sqrt(n)
    me = qt(0.975,(n-1))*se # me for marginal error
    
    figure = barplot(mean, ylim=c(0,5), main=main, names.arg=names.arg)
    arrows(figure, mean+me, figure, mean-me, angle=90, code=3, length=0.1)
    
    output = t.test(y ~ open)$p.value
  }
  
  if(is.null(dim(y))==TRUE & is.null(x)==FALSE) {
    if(is.null(legend)) {legend=levels(group)}
    if(is.null(names.arg)) {names.arg=levels(x)}
    if(is.null(ylim)) {ylim=c(0,max(y, na.rm=T)) }
    
    rn = length(levels(group)) # rn is the number of categories
    cn = length(levels(x)) # cn is the number of y to be plotted
    
    mean = matrix(NA, rn, cn); sd = matrix(NA, rn, cn)
    n = matrix(NA, rn, cn); p = vector(length=cn)
    
    for(i in 1:cn) {
      mean[,i] = tapply(y[x==levels(x)[i]], group[x==levels(x)[i]], mean, na.rm=T)
      sd[,i] = tapply(y[x==levels(x)[i]], group[x==levels(x)[i]], sd, na.rm=T)
      n[,i] = tapply(y[x==levels(x)[i]], group[x==levels(x)[i]], function(x) sum(!is.na(x)))
      p[i] = t.test(y[x==levels(x)[i]] ~ group[x==levels(x)[i]])$p.value
    }
    
    se = sd/sqrt(n)
    me = qt(0.975,(n-1))*se # me for marginal error
    figure = barplot(mean, beside=T, main=main, ylab=ylab, ylim=ylim, 
                     legend=legend, names.arg=names.arg)
    arrows(figure, mean+me, figure, mean-me, angle=90, code=3, length=0.1)
    
    output = list(mean=mean, sd=sd, n=n, se=se, me=me, p=p)
  }
  
  if(is.null(dim(y))==FALSE) {
    if(is.null(legend)) {legend=levels(group)}
    if(is.null(names.arg)) {names.arg=colnames(y)}
    
    rn = length(levels(group)) # rn is the number of categories
    cn = ncol(y) # cn is the number of y to be plotted
    
    mean = matrix(NA, rn, cn); sd = matrix(NA, rn, cn)
    n = matrix(NA, rn, cn); p = vector(length=cn)
    
    for(i in 1:cn) {
      mean[,i] = tapply(y[,i], group, mean, na.rm=T)
      sd[,i] = tapply(y[,i], group, sd, na.rm=T)
      n[,i] = tapply(y[,i], group, function(x) sum(!is.na(x)))
      p[i] = unlist(anova(lm(y[,i]~group))["Pr(>F)"])[1]
    }
    se = sd/sqrt(n)
    me = qt(0.975,(n-1))*se # me for marginal error
    figure = barplot(mean, beside=T, main=main, ylab=ylab, ylim=ylim, 
                     legend=legend, names.arg=names.arg)
    arrows(figure, mean+me, figure, mean-me, angle=90, code=3, length=0.1)
    
    output = list(mean=mean, sd=sd, n=n, se=se, me=me, p=p)
  }
  return(output)
} 
#######################################


###############################
### Plotting Summary Graphs ###
###############################



##############################################################



##################################
### Tabulate Regression Models ###
##################################
tabmod = function(list, digits=3) { # call wordidx(); subsumed under regtab()
  # Input e.g. list = list(mod1, mod2, mod3, mod4) 
  ### Collect IV Names from All Models ###
  ivnames = lapply(list, function(x) rownames(summary(x)$coef))
  ivnames = unique(unlist(ivnames))
  ### Create Regression Table ###
  tab = matrix(NA, ncol=length(list), nrow=2*length(ivnames) + 5)
  colnames(tab) = paste0("Model", seq(length(list)))
  blanks = paste0("NA", seq(length(ivnames)))
  rownames(tab) = c(rbind(ivnames, blanks),c("Cases","DF","Adj. R-squared","AIC","BIC"))
  ### Fill in Statistics ###
  for (i in seq(list)) {
    modsum = summary(list[[i]])
    coef = round(modsum$coef[,1], digits) # Must do roundup here; Will become characters!
    se = round(modsum$coef[,2], digits)
    se = paste("(", se, ")", sep="")
    p = modsum$coef[,4]
    p = ifelse(0.05<p & p<=0.1, paste("{\\super+}"), # \u2020 is hex unicode for †
               ifelse(0.01<p & p<=0.05, paste("{\\super*}"), 
                      ifelse(0.001<p & p<=0.01, paste("{\\super**}"), 
                             ifelse(p<=0.001, paste("{\\super***}"), paste("")))))
    # http://unicodelookup.com/
    coef = paste(coef, p, sep="")
    # Find Row Index of Coefficent Names from Regression Table #
    names = rownames(modsum$coef) 
    index = wordidx(names, rownames(tab))
    # Fill in Statistics #
    tab[index, i] = coef; tab[index+1, i] = se
    tab["Cases",i] = sum(modsum$df[-3])
    tab["DF",i] = modsum$df[2]
    tab["Adj. R-squared",i] = round(modsum$adj.r.squared, digits)
    tab["AIC",i] = round(AIC(list[[i]]), digits-1)
    tab["BIC",i] = round(BIC(list[[i]]), digits-1)
  }
  return(tab)
}

regtable = function(list, docname=NULL, std=FALSE, header=NULL, digits=3, NA.string="#"){ # call tabmod()
  require("rtf")
  
  regmods = tabmod(list, digits=digits)
  if ( is.null(docname) ) { docname = "RegModels.docx" }
  if ( is.null(header) ) { header = "Table X. Regression Models" }
  note.std = " All estimates are standardized beta coefficients. "
  note.se = "Standard errors are in parentheses."
  note.p = "{\\super+} < .10, {\\super*} p < .05, {\\super**} p < .01, {\\super***} p < .001"
  
  rtf = RTF(docname, font.size=12)
  addText(rtf, header)
  addNewLine(rtf)
  addTable(rtf, regmods, header.col.justify="C", col.justify="C", row.names=T, NA.string=NA.string)
  addText(rtf, "Note. ", italic=T)
  if (std) {addText(rtf, note.std)}
  addText(rtf, note.se)
  addNewLine(rtf)
  addText(rtf, note.p)
  done(rtf)
  cat("Regression table is successfully stored in the following directory \n")
  cat(getwd())
}

### Test Unit ###
#tab = matrix(c("0.12{\\super*}", "1.2{\\super**}", "0.04", "0.002"), ncol=2, nrow=2)
#colnames(tab) = c("Estimate","p value"); rownames(tab) = c("var1","var2")
#note.se = "Standard errors are in parentheses."
#note.p = "{\\super+} < .10, {\\super*} p < .05, {\\super**} p < .01, {\\super***} p < .001"

#rtf = RTF("NewDoc.docx", font.size=12)
#addText(rtf, "Table X. Regression Models")
#addNewLine(rtf)
#addTable(rtf, tab, row.names=T)
#addText(rtf, "Note. ", italic=T)
#addText(rtf, note.se)
#addNewLine(rtf)
#addText(rtf, note.p)
#done(rtf)

####################################################################




#####################################
### R Functions used in Stats 202 ### 
#####################################
### Probability ### 
p.ci = function(p, n, confid=0.95) {  
  # p for proportion
  z = qnorm((1-confid)/2, lower.tail=F)
  a = 1 + z^2/n
  b = - (2*p + z^2/n)
  c = p^2
  
  ci.standard = p+c(-1,1)*z*sqrt(p*(1-p)/n)
  ci.reverse = -(b/(2*a))+c(-1,1)*sqrt((b^2-4*a*c)/(4*a^2))
  
  x = matrix(c(ci.standard, ci.reverse), nrow=2, byrow=T)
  row.names(x) = c("CI(Standard)", "CI(Reverse)")
  
  return(x)
}

### Odds Ratio ###
or.tab = function(tab) {
  or = tab[1,1]*tab[2,2]/(tab[1,2]*tab[2,1])
  return(or)
}

or.input = function(a, b, c, d) {
  or = a*d/(b*c)
  return(or)
}

or.ci = function(tab, confid=0.95) {
  or = tab[1,1]*tab[2,2]/(tab[1,2]*tab[2,1])
  logse = sqrt(sum(1/as.vector(xmat)))
  
  alpha = 1-confid
  z = qnorm(alpha/2, lower.tail=F)
  z = c(-z, z)
  
  t = log(or)/logse
  
  logci = log(or)+z*logse
  ci = exp(logci)
  list(or = or, ci=ci)
}

es = function(tm, tsd, cm, csd, tn, cn, se=FALSE){
  if(se) {tsd = tsd * sqrt(tn); csd = csd * sqrt(cn)} # if given Standard Error instead of Standard Deviation
  diff = tm - cm
  sd = tsd*(tn-1)/(tn+cn-2) + csd*(cn-1)/(tn+cn-2)
  es = diff/sd
  return(es)
}

### Relative Risk ###
# input #
# contingency table & confidence level, e.g. 0.95
rr.ci = function(tab, confid=0.95) {
  rowtot = rowSums(tab)
  p.con = tab[,1]/rowtot
  rr = p.con[1]/p.con[2]
  
  alpha = 1-confid
  z = qnorm(alpha/2, lower.tail=F)
  z = c(-z, z)
  
  part = ((1-p.con)/(rowtot*p.con))
  logci = log(rr)+z*sqrt(sum(part))
  ci = exp(logci)
  
  list(rr=rr, ci=ci)
}

### Chi-square & likelihood ratio tests ###
# Test Unit #
#a=(1:5)*10; b=a*5 + sample(-5:5,5); xmat=cbind(a,b)
#a = (1:5)*10; b = a*5; xmat = cbind(a,b)
#chisqtest(xmat)

chisqtest = function(xmat) { 
  # xmat is a contingency table
  nrowtot = rowSums(xmat)
  ncoltot = colSums(xmat)
  prowmar = nrowtot/sum(xmat)
  pcolmar = ncoltot/sum(xmat)
  
  expected = outer(nrowtot, ncoltot) / sum(xmat)
  chisq.cell = (xmat - expected)^2 / expected
  gsq.cell = xmat * log(xmat / expected)
  sclresid = (xmat - expected) / sqrt(expected)
  stdresid = (xmat - expected) / sqrt(expected * outer(1-prowmar, 1-pcolmar))
  
  df = (ncol(xmat)-1)*(nrow(xmat)-1)
  chisq = sum(chisq.cell)
  p.chisq = pchisq(chisq, df, lower.tail=F)
  gsq = 2*sum(gsq.cell)
  p.gsq = pchisq(gsq, df, lower.tail=F)
  
  list(chisq=chisq, gsq=gsq, df=df, p.chisq=p.chisq, p.gsq=p.gsq, expected=expected, sclresid=sclresid, stdresid=stdresid)
}


### Msquare Ordinal Test ###
msqtest = function(xmat, rowscr, colscr) {
  p.cell = xmat/sum(xmat)
  p.row = rowSums(xmat)/sum(xmat)
  p.col = colSums(xmat)/sum(xmat)
  
  rowscr.mean = sum(rowscr*p.row)
  colscr.mean = sum(colscr*p.col)
  rowscr.dev = rowscr-rowscr.mean
  colscr.dev = colscr-colscr.mean
  
  p.cell.temp = p.cell*rowscr.dev
  numerator.cell = apply(p.cell.temp, 1, function(x) x*colscr.dev)
  numerator = sum(numerator.cell)
  denominator = sqrt(sum(rowscr.dev^2*p.row)*sum(colscr.dev^2*p.col))
  
  r = numerator/denominator
  msq = (sum(xmat)-1)*r^2
  p.msq = pchisq(msq, 1, lower.tail=F)
  
  list(r=r, msq=msq, p.msq=p.msq)
  
  #numerator = matrix(ncol=ncol(xmat), nrow=nrow(xmat))
  #for (r in 1:nrow(xmat)) {
  #  for (c in 1:ncol(xmat)) {
  #    numerator[r,c] = (rowscr[r]-row.mean)*(colscr[c]-col.mean)*p.cell[r,c]
  #  }
  #}
}

# test unit #
#Sex = factor(c(rep(1,6),rep(0,6)), labels=c("Female","Male"))
#Race = factor(rep(0:1, 6), labels=c("White","Black"))
#Dose = factor(c(rep(c(1,2,4,8,16,32),times=2)))
#Response = ceiling(c(1,4,9,13,18,20,0,2,6,10,12,16)/2)
#Trials = c(rep(10,12)); Trials[1:3] = ceiling(c(5,6,10)/2)
#xdat = data.frame(Trials, Response, Sex, Race, Dose) # Trials - Yes - Factors

tobinary = function(xdat) {
  # xdat[,1] := # of trials/counts
  # xdat[,2] := # of positive response/yes
  cat('Data rearranged under the assumption: \n\tFirst column is "Trials"; second column is positive "Responses"\n')
  if (any(xdat[,1] - xdat[,2] < 0)) stop("Trials must be greater than Responses")
  newdat = xdat[rep(1:nrow(xdat), xdat[,1]),-1:-2, drop=F] # in case of a single covariate
  yesno = rbind(xdat[,2], xdat[,1]-xdat[,2]) # Yes-No counts
  Response = rep(rep(1:0, ncol(yesno)), times=c(yesno))
  newdat = data.frame(Response, newdat)
  rownames(newdat) = NULL
  return(newdat)
}


tobinom <- function(xdat, dvname=NULL){
  # dvname := a string specificing name of response variable
  if(is.null(dvname)) stop("must provide name of the response variable")
  
  xdat = as.data.frame(xdat)
  index = (seq(colnames(xdat)))[colnames(xdat) == dvname]
  y = xdat[,index]
  data = xdat[,-index, drop=F] # in case of a single covariate
  newdat = aggregate(y, data, FUN=length)
  newdat = cbind( newdat, aggregate( y, data, FUN=sum)[dim(newdat)[2]] )
  n = ncol(newdat)
  newdat = cbind(newdat[,(n-1):n], newdat[,1:(n-2), drop=F])
  colnames(newdat)[1:2] = c("Trials", dvname)
  rownames(newdat) = NULL
  return(newdat)
}
##########################################################################



seecolors = function() {
  par(mar=c(5,2.5,2,2)+0.1)
  plot(c(1,10), c(-1,3), type="n", xaxt="n", cex.lab = 1.2,
       xlab="col = 1:8\nlty = 1:6", main="pch = 0:25")
  axis(1, at=1:10)
  points(x=c(1,rep(1:10, 2),1:5), y=c(-1,rep(0:1, each=10),rep(2,5)), pch=0:25, lwd=1.5)
  abline(v=seq(1,10,by=1), col=1:10, lty=1:10, lwd=1.5)
  #http://www.statmethods.net/advgraphs/parameters.html
  #pch=c(0,2,5,6,15,17)
  #col=c(1,2,4,6,7,5)
  #lty=c(1,2,4,5,6,3)
}


