has.interaction <- function(x,terms){
  out <- sapply(terms,function(i){
    sum(1-(strsplit(x,":")[[1]] %in% strsplit(i,":")[[1]]))==0
  })
  return(sum(out)>0)
}

# Function Model.select
# verbose=T gives the F-tests, dropped var and resulting model after
model.select <- function(model, keep, sig=0.05,verbose=F){
  counter=1
  # check input
  if(!is(model,"lm")) stop(paste(deparse(substitute(model)),"is not an lm object\n"))
  # calculate scope for drop1 function
  terms <- attr(model$terms,"term.labels")
  if(missing(keep)){ # set scopevars to all terms
    scopevars <- terms
  } else{            # select the scopevars if keep is used
    index <- match(keep,terms)
    # check if all is specified correctly
    if(sum(is.na(index))>0){
      novar <- keep[is.na(index)]
      warning(paste(
        c(novar,"cannot be found in the model",
          "\nThese terms are ignored in the model selection."),
        collapse=" "))
      index <- as.vector(na.omit(index))
    }
    scopevars <- terms[-index]
  }
  
  # Backward model selection : 
  
  while(T){
    # extract the test statistics from drop.
    test <- drop1(model, scope=scopevars,test="F")
    
    if(verbose){
      cat("-------------STEP ",counter,"-------------\n",
          "The drop statistics : \n")
      print(test)
    }
    
    pval <- test[,dim(test)[2]]
    
    names(pval) <- rownames(test)
    pval <- sort(pval,decreasing=T)
    
    if(sum(is.na(pval))>0) stop(paste("Model",
                                      deparse(substitute(model)),"is invalid. Check if all coefficients are estimated."))
    
    # check if all significant
    if(pval[1]<sig) break # stops the loop if all remaining vars are sign.
    
    # select var to drop
    i=1
    while(T){
      dropvar <- names(pval)[i]
      check.terms <- terms[-match(dropvar,terms)]
      x <- has.interaction(dropvar,check.terms)
      if(x){i=i+1;next} else {break}              
    } # end while(T) drop var
    
    if(pval[i]<sig) break # stops the loop if var to remove is significant
    
    if(verbose){
      cat("\n--------\nTerm dropped in step",counter,":",dropvar,"\n--------\n\n")              
    }
    
    #update terms, scopevars and model
    scopevars <- scopevars[-match(dropvar,scopevars)]
    terms <- terms[-match(dropvar,terms)]
    
    formul <- as.formula(paste(".~.-",dropvar))
    model <- update(model,formul)
    
    if(length(scopevars)==0) {
      warning("All variables are thrown out of the model.\n",
              "No model could be specified.")
      return()
    }
    counter=counter+1
  } # end while(T) main loop
  return(model)
}


##### Functions for comparing models with partly transformed data #####

# Computes the jacobian and other data related to the box-cox
BoxCoxData<-function(x,lambda){
  # Transforming the data
  xlambda <- (x^lambda-1)/lambda
  if (lambda==0){
    xlambda = log(x)
  }
  n<-length(xlambda)
  
  # Finding Jacobian in each point 
  Jacobian <- (x^(lambda-1))
  # Finding optimal mu and sigma 
  mu <- 1/n*sum(xlambda);
  sigma2 <- 1/n*sum((xlambda-mu)^2)
  
  # Returning
  returnList <- list(x = xlambda, jacobian = Jacobian, mu = mu, sigma2 = sigma2)
  return(returnList)
}

# Calculating the AIC of the model given Jacobian
AICRegModel<-function(loglik,Jacobian,n){
  AIC <- 2*n - 2*loglik - 2*Jacobian
  return(AIC)
}
# Calculating the loglikelihood given the Jacobian
logLikJacobian<-function(loglik,Jacobian){
  logLikJacobian<-loglik+Jacobian
  return(logLikJacobian)
}