# agd.r
#
# Tools for course Analysis of Growth Data
### SvB 6may2012

z2y <- function(z   = c(-2, 0, 2), 
                x   = 1, 
                sex = "M", 
                sub = "N",
                ref = nl4.hgt,
                dist = "LMS",
                dec  = 3,
                sex.fallback = "M",
                sub.fallback = "N")
{
  
  z2y.grp <- function(z, x, ref, dist = "LMS") {
       
    if (dist=="NO") {
  	  check.names(df=ref, needed=c("x","mean","sd"))
  	  mean <- approx(x=ref[,"x"], y=ref[,"mean"], xout=x)$y
  	  sd   <- approx(x=ref[,"x"], y=ref[,"sd"],   xout=x)$y 
      return(mean + z*sd)   
    }
  	if (dist=="LMS") {
  	  check.names(df=ref, needed=c("x","L","M","S"))
  	  L <- approx(x=ref[,"x"], y=ref[,"L"], xout=x)$y
  	  M <- approx(x=ref[,"x"], y=ref[,"M"], xout=x)$y
  	  S <- approx(x=ref[,"x"], y=ref[,"S"], xout=x)$y
  	  return(ifelse(L>0.01|L<(-0.01),M*(1+L*S*z)^(1/L),M*exp(S*z)))
    }
  	if (dist=="BCCG") {
  	  check.names(df=ref, needed=c("x","nu","mu","sigma"))
  	  nu    <- approx(x=ref[,"x"], y=ref[,"nu"], xout=x)$y
  	  mu    <- approx(x=ref[,"x"], y=ref[,"mu"], xout=x)$y
  	  sigma <- approx(x=ref[,"x"], y=ref[,"sigma"], xout=x)$y
  	  return(qBCCG(pnorm(z), mu=mu, sigma=sigma, nu=nu))
    }
  	if (dist=="BCPE") {
  	  check.names(df=ref, needed=c("x","nu","mu","sigma","tau"))
  	  mu    <- approx(x=ref[,"x"], y=ref[,"mu"], xout=x)$y
  	  sigma <- approx(x=ref[,"x"], y=ref[,"sigma"], xout=x)$y
  	  nu    <- approx(x=ref[,"x"], y=ref[,"nu"], xout=x)$y
  	  tau   <- approx(x=ref[,"x"], y=ref[,"tau"], xout=x)$y
  	  return(qBCPE(pnorm(z), mu=mu, sigma=sigma, nu=nu, tau=tau))
    }
  	if (dist=="BCT") {
  	  check.names(df=ref, needed=c("x","nu","mu","sigma","tau"))
  	  mu    <- approx(x=ref[,"x"], y=ref[,"mu"], xout=x)$y
  	  sigma <- approx(x=ref[,"x"], y=ref[,"sigma"], xout=x)$y
  	  nu    <- approx(x=ref[,"x"], y=ref[,"nu"], xout=x)$y
  	  tau   <- approx(x=ref[,"x"], y=ref[,"tau"], xout=x)$y
  	  return(qBCT(pnorm(z), mu=mu, sigma=sigma, nu=nu, tau=tau))
    }
    
    stop(paste("Reference type", dist, "not implemented."))
  }

  if (!is.data.frame(ref)) stop("'ref' should be a data frame.")
  n <- length(z)
  if (n < 1)         stop("'z' must have 1 or more values")
  if(!is.vector(z))  stop("'z' must be a numeric vector")
  if(!is.numeric(z)) stop("'z' must be a numeric vector")
  
  x   <- rep(x,   length.out=length(z))
  sex <- rep(sex, length.out=length(z))
  sub <- rep(sub, length.out=length(z))
  dist <- match.arg(dist, choices=c("NO","LMS","BCCG","BCPE","BCT"))
  
  # available levels in ref: sex, sub
  lev.sex <- levels(ref$sex[, drop=TRUE])
  lev.sub <- levels(ref$sub[, drop=TRUE])

  # replace nomatching levels
  idx <- is.na(match(sub, lev.sub))
  if (any(idx)) {
    sub[idx] <- sub.fallback
    warning("Entries (n=",sum(idx),") replaced by '",sub.fallback,"'",sep="")
  }
  idx <- is.na(match(sex, lev.sex))
  if (any(idx)) {
    sex[idx] <- sex.fallback
    warning("Entries (n=",sum(idx),") replaced by '",sex.fallback,"'",sep="")
  }

  refs <- with(ref,split(ref, f=list(sub, sex)))
	xs <- split(x,list(sub, sex))
	zs <- split(z,list(sub, sex))
  ys <- vector("list",length(zs))
  names(ys) <- names(zs)

  for(i in 1:length(zs)) {
      name <- names(zs)[i]
      if(is.null(refs[[name]])) ys[[name]] <- rep(NA,length=length(zs[[name]]))
      else ys[[name]] <- z2y.grp(z=zs[[name]], x=xs[[name]], ref=refs[[name]], 
                         dist=dist)
  }

  y <- unsplit(ys,f=list(sub,sex))
  names(y) <- names(z)
  return(round(y,dec))
}
  


y2z <- function(y   = c(75, 80, 85), 
                x   = 1, 
                sex = "M", 
                sub = "N",
                ref = nl4.hgt,
                dist = "LMS",
                dec  = 3,
                sex.fallback = "M",
                sub.fallback = "N")
{

  y2z.grp <- function(y, x, ref, dist = "LMS", dec=3){
  
    if (dist=="NO") {
  	  check.names(df=ref, needed=c("x","mean","sd"))
  	  mean <- approx(x=ref[,"x"], y=ref[,"mean"], xout=x)$y
  	  sd   <- approx(x=ref[,"x"], y=ref[,"sd"],   xout=x)$y 
      return((y-mean)/sd)   
    }
  	if (dist=="LMS") {
  	  check.names(df=ref, needed=c("x","L","M","S"))
  	  L <- approx(x=ref[,"x"], y=ref[,"L"], xout=x)$y
  	  M <- approx(x=ref[,"x"], y=ref[,"M"], xout=x)$y
  	  S <- approx(x=ref[,"x"], y=ref[,"S"], xout=x)$y
  	  return(ifelse(L>0.01|L<(-0.01),(((y/M)^L)-1)/(L*S),log(y/M)/S))
    }
  	if (dist=="BCCG") {
  	  check.names(df=ref, needed=c("x","nu","mu","sigma"))
  	  nu    <- approx(x=ref[,"x"], y=ref[,"nu"], xout=x)$y
  	  mu    <- approx(x=ref[,"x"], y=ref[,"mu"], xout=x)$y
  	  sigma <- approx(x=ref[,"x"], y=ref[,"sigma"], xout=x)$y
  	  return(qnorm(pBCCG(y, mu=mu, sigma=sigma, nu=nu)))
    }
  	if (dist=="BCPE") {
  	  check.names(df=ref, needed=c("x","nu","mu","sigma","tau"))
  	  mu    <- approx(x=ref[,"x"], y=ref[,"mu"], xout=x)$y
  	  sigma <- approx(x=ref[,"x"], y=ref[,"sigma"], xout=x)$y
  	  nu    <- approx(x=ref[,"x"], y=ref[,"nu"], xout=x)$y
  	  tau   <- approx(x=ref[,"x"], y=ref[,"tau"], xout=x)$y
  	  return(qnorm(pBCPE(y, mu=mu, sigma=sigma, nu=nu, tau=tau)))
    }
  	if (dist=="BCT") {
  	  check.names(df=ref, needed=c("x","nu","mu","sigma","tau"))
  	  mu    <- approx(x=ref[,"x"], y=ref[,"mu"], xout=x)$y
  	  sigma <- approx(x=ref[,"x"], y=ref[,"sigma"], xout=x)$y
  	  nu    <- approx(x=ref[,"x"], y=ref[,"nu"], xout=x)$y
  	  tau   <- approx(x=ref[,"x"], y=ref[,"tau"], xout=x)$y
  	  return(qnorm(pBCT(y, mu=mu, sigma=sigma, nu=nu, tau=tau)))
    }
  	if (dist=="BCCG") {
  	  check.names(df=ref, needed=c("x","mu","sigma","nu"))
  	  mu    <- approx(x=ref[,"x"], y=ref[,"mu"], xout=x)$y
  	  sigma <- approx(x=ref[,"x"], y=ref[,"sigma"], xout=x)$y
  	  nu    <- approx(x=ref[,"x"], y=ref[,"nu"], xout=x)$y
  	  return()
  	  #return(ifelse(L>0.01|L<(-0.01),(((y/M)^L)-1)/(L*S),log(y/M)/S))
    }
    stop(paste("Reference type", dist, "not implemented."))
  }	

  if (!is.data.frame(ref)) stop("'ref' should be a data frame.")
  n <- length(y)
  if (n < 1)         stop("'y' must have 1 or more values")
  if(!is.vector(y))  stop("'y' must be a numeric vector")
  if(!is.numeric(y)) stop("'y' must be a numeric vector")
  
  x   <- rep(x,   length.out=length(y))
  sex <- rep(sex, length.out=length(y))
  sub <- rep(sub, length.out=length(y))
  dist <- match.arg(dist, choices=c("NO","LMS","BCCG","BCPE","BCT"))
  
  # available levels in ref: sex, sub
  lev.sex <- levels(ref$sex[, drop=TRUE])
  lev.sub <- levels(ref$sub[, drop=TRUE])
  
  # replace nomatching levels
  idx <- is.na(match(sub, lev.sub))
  if (any(idx)) {
    sub[idx] <- sub.fallback
    warning("Entries (n=",sum(idx),") replaced by '",sub.fallback,"'",sep="")
  }
  idx <- is.na(match(sex, lev.sex))
  if (any(idx)) {
    sex[idx] <- sex.fallback
    warning("Entries (n=",sum(idx),") replaced by '",sex.fallback,"'",sep="")
  }

  refs <- with(ref,split(ref, f=list(sub, sex)), drop = TRUE)
	xs <- split(x,list(sub, sex), drop = TRUE)
	ys <- split(y,list(sub, sex), drop = TRUE)
  zs <- vector("list",length(ys))
  names(zs) <- names(ys)

  for(i in 1:length(ys)) {
      name <- names(ys)[i]
      if(is.null(refs[[name]])) ys[[name]] <- rep(NA,length=length(ys[[name]]))
      else zs[[name]] <- y2z.grp(y=ys[[name]], x=xs[[name]], ref=refs[[name]], 
                         dist=dist)
  }

  z <- unsplit(zs,f=list(sub,sex))		
  names(z) <- names(y)
  return(round(z, dec))
}

check.names <- function(df, needed){
    if (missing(df)) stop("required argument 'df' not found")
    if (missing(needed)) stop("required argument 'needed' not found")
    
    notfound <- is.na(match(needed, names(df)))
 	  if (any(notfound)) stop("Not found: ",paste(needed[notfound],collapse=", "))
}

