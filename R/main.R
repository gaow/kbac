KbacTest <- function(casectrl.dat, alpha, num.permutation, quiet = T, maf.upper.bound = 1.0, alternative = 1) {
        
    ydatIn <- as.matrix(casectrl.dat[,1])
    xmat <- as.matrix(casectrl.dat[,-1])
    mafIn <- apply(xmat, 2, function(x) sum(x[which(x > 0)])) / (length(ydatIn) * 2)
    xmat[which(xmat != 1 & xmat != 2 & xmat != 0)] <- 0
    xdatIn <- matrix(t(xmat), nrow = 1)
    xcol <- ncol(xmat)
    ylen <- nrow(xmat)
    nn <- num.permutation
    aa <- alpha
    qq <- quiet
    mafUpper <- maf.upper.bound
    pvalue <- KbacGetP(nn, qq, aa, mafUpper, xdatIn, ydatIn, mafIn, xcol, ylen, alternative)
    return(pvalue)
}
