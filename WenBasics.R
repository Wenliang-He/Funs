#################################
### Wenliang's Basics Backage ###
#################################

########################################
### Data Inspection and Manipulation ### Five Sets of Functions
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
    
    # Reprot Types-of-Fix #
    if( any(typeidx) ) { cat(typesofFix[typeidx])
    } else { cat("There is nothing to fix; Nothing is modified.") }
  }
  
  # Reprot Types-of-Fix #
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

wordidx=function(words, xvecdat) { # Can handle multiple words input
  # called by tabmod()
  if (is.null(dim(xvecdat))) { index = seq(xvecdat)[xvecdat==words] } # if vector
  else { index = seq(ncol(xvecdat))[colnames(xvecdat)==words] } # if dataframe or matrix
  return(index)
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

################################################################################



############################################
### Get Summary & Descriptive Statistics ### One Set of Functions
############################################
### 1. Vector Summary: Handle Two-Level Factors Only ###
summary.vec = function(xvec, confid=0.95, digits=3, short=TRUE) { # Completely subsumed under getsummary()
  if (factor <- is.factor(xvec)) {xvec = as.numeric(xvec)-1 } # convert two-level factors into 0s and 1s
  alpha = 1-confid; t = qt(1-alpha/2, length(xvec)-1) # student t distribution
  # Sample Statistics #
  total = length(xvec); missing = sum(is.na(xvec)); n = total-missing # non-missing cases
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
  
  if (short) {
    tab = matrix(c(total, missing, sampleSize, mean, sd, min, median, max), nrow=1, byrow=F, 
                 dimnames=list("xvec", c("Total","Missing","Count","Mean","SD","Min","Median","Max")))
    } else {
    tab = matrix(c(total, missing, sampleSize, mean, sd, min, median, max, skew, kurtosis, se, ci), nrow=1, byrow=F, 
                 dimnames=list("xvec", c("Total","Missing","Count","Mean","SD","Min","Median","Max",
                                         "Skewness","Kurtosis","Std.Err.","[95% Conf.","Interval]")))
    }
  return(round(tab, digits))
}

### 2. Generate Two-level Factors from a Multi-level Factor ###
genfactor = function(fvec) { 
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
getsummary = function(xvecmat, confid=0.95, digits=3, short=TRUE) {
  # Call summary.vec() & genfactor()
  if (is.null(dim(xvecmat))) { # if xvecmat is a vector
    tab = summary.vec(xvecmat, confid=confid, digits=digits, short=short)
  } else { # if xvecmat is a matric or data.frame
    class = sapply(xvecmat, class)
    fidx = c(1:ncol(xvecmat))[class == "factor"] # factor index
    if (length(fidx)!=0) { # when factors are present
      tdat = xvecmat[,-fidx] # base cummulator with all non-factor variables
      for (i in 1:length(fidx)) { tdat = data.frame(tdat, genfactor(xvecmat[,fidx[i]])) }
      xvecmat = tdat
    }
    tab=t(sapply(xvecmat, function(x) summary.vec(x, confid=confid, digits=digits, short=short)))
    
    if (short) { colnames(tab) =c("Total","Missing","Count","Mean","SD","Min","Median","Max")
    } else { colnames(tab) = c("Total","Missing","Count","Mean","SD","Min","Median","Max",
                               "Skewness","Kurtosis","Std.Err.","[95% Conf.","Interval]") 
    }
  }
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
# deci = demical point

# Limitations #
# cat can be only two-level; "two"-sample t-tests are used for continuous variables 
# cannot handle ordinal variable

cp.vec = function(xvec, cat, deci=3) { # used to be called "des.tab"
  # call summary.vec() & chisqtest()
  lvls = levels(cat)
  if (class(xvec)=="integer" || class(xvec)=="numeric" || length(unique(xvec))>=10) { # xvec is continuous
    tab = matrix(NA, nrow=1, ncol=3, byrow=F) # create an empty row #
    # fill categories #
    for (i in 1:length(lvls)) {
      mxvec = round(summary.vec(xvec[cat==lvls[i]]), digits=deci)
      tab[,i] = paste(mxvec[colnames(mxvec)=="Mean"], "±", mxvec[colnames(mxvec)=="SD"]) # "±" is "\u00b1"
      # http://www.utf8-chartable.de/unicode-utf8-table.pl
    }
    # fill statistic tests #
    tstat = round(t.test(xvec~cat)$stat, deci)
    pval = round(t.test(xvec~cat)$p.value, deci+1)
    tab[ncol(tab)] = paste(tstat, " (", pval, ")", sep="")
  }
  
  else {# If xvec is categorical
    tab = matrix(NA, nrow=length(unique(na.omit(xvec))), ncol=3, byrow=F) # create empty rows #
    # fill categories #
    freq = table(xvec, cat)
    percent = round(100*freq/colSums(t(freq)), deci)
    #percent = round(100*t(t(freq)/colSums(freq)), deci)
    tab[,-ncol(tab)] = matrix(paste(freq, " (", percent, "%)", sep=""), 
                              nrow=length(unique(na.omit(xvec))))
    # fill statistic tests #
    #if (is.ordered=="TRUE") {
    ##### need to figure out rowscr and colscr #########
    #  chisq = round(msqtest(freq, rowscr, colscr)$chisq, deci)
    #  pval = round(msqtest(freq, rowscr, colscr)$p.chisq, deci+1)
    #} else {
    chisqtest = chisqtest(freq)
    chisq = round(chisqtest$chisq, deci)
    pval = round(chisqtest$p.chisq, deci+1)
    #} 
    tab[,ncol(tab)] = paste(chisq, " (", pval, ")", sep="")
  }
  return(tab)
}

### 2. Any Comparison ###
cp.table = function(xdat, cat) {
  # Get the Statistics #
  stats = lapply(xdat, function(x) cp.vec(x, cat))
  
  # Put Results in One Table #
  stats.tab = NULL
  for (i in 1:length(stats)) {
    stats.tab = rbind(stats.tab, stats[[i]])
  }
  
  # Add Column Names #
  head = c(rep("Mean \u00b1 SD or N (%)",2), "t(p) or Chi-squared(p)")
  colnames(stats.tab) = head
  
  # Add Row Names #
  varnames = colnames(xdat)
  lvls = lapply(xdat, levels)
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
  par(mfrow=c(1,1)); n = ncol(xdat); colnames = colnames(xdat)
  if (is.null(varnames)) { varnames = colnames }
  
  if (write) { # if output to hard disk
    for (i in 1:n) {
      png(paste0(figname,".",colnames[i],".png"), width=width, height=height)
      par(mfrow=c(1,1))
      
      hist.vec(xdat[,i], varname=varnames[i], ...)
      
      dev.off()
    }
  } else { # if plot in R
    for (i in 1:n) { hist.vec(xdat[,i], varname=varnames[i], ...) }
  }
}

cp.barplot = function(xdat, beside=TRUE, density=NULL, nmax=6, ...){
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
  # xdat must include x & y, where both must be continuous 
  # z is optional and z can be continuous or a factor variable of multiple levels
  # nf is the number of levels to split if z is continuous
  if(is.null(xlab)) {xlab=colnames(xdat)[1]}
  if(is.null(ylab)) {ylab=colnames(xdat)[2]}
  if(is.null(main)) {main=paste("Plotted by",colnames(xdat)[3],sep=" ")}
  
  par(xpd=FALSE) # To ensure abline() not plotting outside the plot region
  x = xdat[,1]; y = xdat[,2]
  rx = range(x, na.rm=T); ry = range(y, na.rm=T)
  
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
    grp = cut(xvec, breaks = quantile(xvec, seq(0, 1, 1/i)), include.lowest = T)
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
  if (length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
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

