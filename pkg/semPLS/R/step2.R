### Inner estimation of factor scores
step2 <-
function(Latent, innerW, model, pairwise){
  blocks <- model$blocks
  if(pairwise){
    fscores <- matrix(NA, nrow=nrow(Latent), ncol=ncol(Latent))
    colnames(fscores) <- model$latent
    for(i in model$latent){
      con <- which(innerW[,i]!=0)
      fscores[,i] <- Latent[,con, drop=FALSE] %*%
                     innerW[con,i, drop=FALSE]
    }
    Latent <- scale(fscores)
  }
  else {Latent <- scale(Latent %*% innerW)}
  # the attributes for the scale are meaningless
  attributes(Latent)[c(3,4)] <- NULL
  return(Latent)
}

### Inner estimation of factor scores by gam
step2gam <-
function(Latent, innerW, model){
  ## Latent is a matrix of component scores
  ## model is of class 'plsm'
  ## gamList is a list of fitted mgcv::gam done by smoothPathWeighting 
  Latent <- scale(as.matrix(as.data.frame(innerW)))
  # the attributes for the scale are meaningless
  attributes(Latent)[c(3,4)] <- NULL
  return(Latent)
}

