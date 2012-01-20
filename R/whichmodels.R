# Converts a number to a number with leading 0's specified by num.digits.
# So, e.g.,  NumConv(9,2) produces "09", and NumConv(10,2) produces "10".
NumConv <- function(n, num.digits){
    if (n == 0){
        return(paste(rep(0, num.digits), collapse=""))
    } else {
        .tmp <- paste(rep(0, (num.digits-floor(log(n,10))-1)), collapse="")
        return(paste(.tmp, n, sep=""))
    }
}

# Collapsing the max.'s into single values so we can (equivalent of) lapply over
# them. Here we make the choice that no one of the max.'s can be no larger than
# 99. Anyone who wants higher orders is clinically insane.
WhichModels <- function(max.p, max.q, max.P, max.Q, maxK){
    total.models <- (max.p+1)*(max.q+1)*(max.P+1)*(max.Q+1)*length(0:maxK)
    x <- numeric(total.models)
    i <- 1

    for (x1 in 0:max.p) for (x2 in 0:max.q){
        for (x3 in 0:max.P) for (x4 in 0:max.Q){
            for (K in 0:maxK){
                x[i] <- paste(K+1, NumConv(x1,2), NumConv(x2,2), NumConv(x3,2), NumConv(x4,2), sep="")
                i <- i+1
            }
        }
    }
    as.numeric(x)
}

# Takes a concatenated number from output of WhichModels function and converts
# it back to its component parts.  Returns vector of p, q, P, Q, K to throw into
# myarima().
UndoWhichModels <- function(n){
    K <- floor(n/10^8)
    if (K==2){
        n <- n-10^8
        K <- TRUE
    } else K <- FALSE
    
    i <- floor(n*10^(-6)) - 1e2
    j <- floor(n*10^(-4)) - 1e4 - i*1e2
    I <- floor(n*10^(-2)) - 1e6 - i*1e4 - j*1e2
    J <- n - 1e8 - i*1e6 - j*1e4 - I*1e2
    return(c(i, j, I, J, K))
}
