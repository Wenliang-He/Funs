#######################
### Clean Sentences ###
#######################
# test unit #
#s = "2 lists of 2nd C2 H2O 1 12 123"
#s = "  12 about I'm 2nd C2, H2O;  1230. & base\\basis; yes/no : Cl- \n \t \r \\b \\  1>34 about * (page[4]), H+? 5  "

strip = function(str, left=FALSE, right=FALSE) {
  if (left) { str = gsub("^ +", "", str) }
  if (right) { str = gsub(" +$", "", str) }
  if (!left && !right) {
    str = gsub("^ +", "", str)
    str = gsub(" +$", "", str)
    str = gsub("\\s+", " ", str)
  }
  return(str)
}


cleanst = function(c, char="'") {
  # char specifies all special characters to keep
  # convert \b into b"
  c = gsub("\\\\"," ",c)
  c = gsub("\\\b","b",c)
  # convet to lower case #
  c = tolower(c)
  # remove special characters except for "char" #
  c = gsub(paste0("[^[:alnum:]", char, "]"), " ", c)    # gsub("[^[:alnum:]'?&+/\\-]", " ", c)
  # remove numbers #
  c = gsub("\\b\\d+\\b", " ", c) # replace _*_ by _ where _ is a boundary and * is a number of any digits
  # remove whitespaces #
  c = gsub(" +", " ", c) # remove all whitespaces in between words
  c = gsub("^ +", "", c) # remove all whitespaces in the front; ^ outside [ ] means "beginning" of a line
  c = gsub(" +$", "", c) # remove all whitespaces in the back; $ outside [ ] means "end" of a line
  return(c)
}
# http://stackoverflow.com/questions/21641522/how-to-remove-specific-special-characters-in-r
# http://stackoverflow.com/questions/21653294/how-to-remove-only-actual-numbers-from-a-string-of-characters-in-r
##################################
### Find Unique Seperate Words ###
##################################
# test unit #
#c = "12 a about above I'm 2nd C2 & base 123 base base; H2O; yes/no :   1>34 with * (page[4]), H+? 5"
#c = cleanst(c)
sepwords = function(c, unique=TRUE, stopwords=NULL) { # Old_Version: unique = "yes"
  # called by count.attr()
  # split into words #
  out = unlist(strsplit(c, split=" "))
  # keep only unique words #
  if (unique) { out = unique(out) }
  # remove stop words #
  if ( !is.null(stopwords) ) { out = setdiff(out, stopwords) }
  return(out)
}

####################
### Find N-Grams ###
####################
# test unit #
#c = "12 a about above I'm 2nd C2 & base 123 base base; H2O; yes/no :   1>34 with * (page[4]), H+? 5, yes, a c h d; d"
#c = cleanst(c)

# getngram() gets all n-grams at n
# This is an advanced version of sepwords()
getngram = function(x, n, unique=TRUE) { # Old_Version: unique = "yes"
  # x is a sentence
  words = unlist(strsplit(x, split=" "))
  count = length(words)
  output = sapply(1:(count-(n-1)), function(x) paste(words[x: (x+(n-1))], collapse=" "))
  if (unique) {output = unique(output)}
  return(output)
}

# getngrams() gets all n-grams from 1 to n (i.e. cumulative ngrams)
getngrams = function(x, n, unique=TRUE) {
  count = vector(); ngrams = vector()
  for (i in 1:n) {
    temp = getngram(x, i, unique=unique)
    count[i] = length(temp)
    ngrams = c(ngrams, temp)
  }
  list(ngrams = ngrams, count = count)
}

################################################################


#######################
### Find Attributes ###
#######################
# test unit #
#c = "12 a about above I'm 2nd C2 & base 123 base base; H2O; yes/no :   1>34 with * (page[4]), H+? 5, yes, a c h d; d"
#c = cleanst(c)
# in the future write a function to count keywords from a sentence
# the following function will be greatly simplied then

count.attr = function(xvec, keywords, ngrams.count, count=TRUE) { # need to fix a bug
  # call getngram()
  # xvec is a vector, but each of its elements is a sentance
  pos = cumsum(ngrams.count) # positions of the last 1-, 2-, and 3-grams
  
  ymat = as.data.frame(matrix(0, nrow=length(xvec), ncol=length(keywords)))
  colnames(ymat) = keywords
  for (i in 1:length(keywords)) {
    n = which((pos-i)>=0)[1] # n as in n-gram is the first positive position
    ymat[,i] = as.vector(sapply(xvec, function(x) sum(keywords[i]==getngram(x, n, unique=FALSE)) ))
  }

  if ( !count ) {ymat = apply(ymat, 2, function(x) as.numeric(x>0))}
  return(ymat)
}

count.attr = function(xvec, keywords, ngrams.count, count=TRUE) { # meanwhile, use this one!
  # call sepwords()
  # although xvec is a vector, each of its elements is a sentance
  n = ngrams.count[1] # n is the number of 1-gram
  ymat = as.data.frame(matrix(0, nrow=length(xvec), ncol=length(keywords)))
  colnames(ymat) = keywords
  
  for (i in 1:n) { # all the 1-grams
    ymat[,i] = as.vector(sapply(xvec, function(x) 
      sum(keywords[i]==sepwords(x, unique=FALSE)) ))
  }
  
  for (i in (n+1):length(keywords)) { # all the 2plus-grams
    ymat[,i] = as.vector(sapply(xvec, function(x) 
      ifelse(length(grep(keywords[i], x))==1,1,0) )) 
  }
  
  if ( !count ) {ymat = apply(ymat, 2, function(x) as.numeric(x>0))}
  return(ymat)
}


#########################
### Seperate Datasets ###
#########################

sepdata = function(dat, prop=0.8) {
  m = nrow(dat) # total number of rows
  set.seed(123)
  rn = sample(m) # rn for random number
  n = ceiling(m*prop) # number of rows to be chosen
  
  chosen = dat[rn[1:n],]
  rownames(chosen)=NULL
  rest = dat[rn[-1:-n],]
  rownames(rest)=NULL
  
  list(chosen=chosen, rest=rest)
}

