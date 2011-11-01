# TODO: Add comment
# 
# Author: srazbash
###############################################################################
checkEigenValues<-function(eigen.values) {
	abs.eigen.values<-round(abs(eigen.values), digits=14)
	
	if(all((abs.eigen.values < 1))) {
		return(TRUE)
	} else if(all((abs.eigen.values <= 1))) {
		#print("stage one")
		if(!all((abs.eigen.values == 1))) {
			#print("stage two:")
			#print(abs.eigen.values)
			return(TRUE)
		} else {
			#print("Eigen D-2:")
			#print(abs.eigen.values)
			return(FALSE)
		}
	} else {
		#print("Eigen D:")
		#print(abs.eigen.values)
		return(FALSE)
	}
}

checkAdmissibility<-function(D, box.cox=NULL, small.phi=NULL, ar.coefs=NULL, ma.coefs=NULL) {
	#Check the range of the Box-Cox parameter
	if(!is.null(box.cox)) {
		if((box.cox < 0) | (box.cox > 1.5)) {
			return(FALSE)
		}
	}
	
	#Check the range of small.phi
	if(!is.null(small.phi)) {
		if(((small.phi < .8) | (small.phi > .98)) & (small.phi != 1)) {
			return(FALSE)
		}
	}
	
	#Check AR part for stationarity
	if(!is.null(ar.coefs)) {
		#print("as.coefs")
		arCheck <- function(ar) {
			p <- max(which(c(1, -ar) != 0)) - 1
			if (!p) 
				return(TRUE)
			all(Mod(polyroot(c(1, -ar[1L:p]))) > 1)
		}
		if(!arCheck(ar.coefs)) {
			return(FALSE)
		}
	}
	
	#Check MA part for invetibility
	if(!is.null(ma.coefs)) {
		#print("ma.coefs")
		maInvert <- function(ma) {
			q <- length(ma)
			q0 <- max(which(c(1, ma) != 0)) - 1L
			if (!q0) 
				return(ma)
			roots <- polyroot(c(1, ma[1L:q0]))
			ind <- Mod(roots) < 1
			if (all(!ind)) 
				return(ma)
			if (q0 == 1) 
				return(c(1/ma[1L], rep(0, q - q0)))
			roots[ind] <- 1/roots[ind]
			x <- 1
			for (r in roots) x <- c(x, 0) - c(0, x)/r
			c(Re(x[-1L]), rep(0, q - q0))
		}
		inverted.ma<-maInvert(ma.coefs)
		if(all(inverted.ma != ma.coefs)) {
			return(FALSE)
		}
	}
	
	#Check the eigen values of the D matrix
	D.eigen.values<-eigen(D, only.values=TRUE)$values
	#a<<-D.eigen.values
	#DD<<-D
	#print(D.eigen.values)
	answer<-checkEigenValues(D.eigen.values)
	#print(answer)
	return(answer)
}
