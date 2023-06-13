# Functions ------------------------------------------------------
# 楕円フーリエ解析 (Elliptic Fourier Analysis)
efourier <- function(M, n=dim(M)[1]/2) {
    p <- dim(M)[1]
    Dx <- M[,1]-M[c(p, (1:p-1)), 1]
    Dy <- M[,2]-M[c(p, (1:p-1)), 2]
    Dt <- sqrt(Dx^2+Dy^2)
    t1 <- cumsum(Dt)
    t1m1 <- c(0, t1[-p])
    t <- sum(Dt)
    an <- bn <- cn <- dn <- numeric(n)
    for (i in 1:n) {
        an[i] <- (t/(2*pi^2*i^2))*sum((Dx/Dt)*(cos(2*i*pi*t1/t)-cos(2*pi*i*t1m1/t)))
        bn[i] <- (t/(2*pi^2*i^2))*sum((Dx/Dt)*(sin(2*i*pi*t1/t)-sin(2*pi*i*t1m1/t)))
        cn[i] <- (t/(2*pi^2*i^2))*sum((Dy/Dt)*(cos(2*i*pi*t1/t)-cos(2*pi*i*t1m1/t)))
        dn[i] <- (t/(2*pi^2*i^2))*sum((Dy/Dt)*(sin(2*i*pi*t1/t)-sin(2*pi*i*t1m1/t)))
    }
    ao <- 2*sum(M[,1]*Dt/t)
    co <- 2*sum(M[,2]*Dt/t)
    list(ao=ao, co=co, an=an, bn=bn, cn=cn, dn=dn)
}

# Normalize Elliptic Fourier
NEF <- function(M, n=dim(M)[1]/2, start=FALSE) {
    ef <- efourier(M, n)
    A1 <- ef$an[1]
    B1 <- ef$bn[1]
    C1 <- ef$cn[1]
    D1 <- ef$dn[1]
    theta <- 0.5*atan(2*(A1*B1*C1*D1)/(A1^2+C1^2-B1^2-D1^2))
    Aa <- A1*cos(theta)+B1*sin(theta)
    Cc <- C1*cos(theta)+D1*sin(theta)
    scale <- sqrt(Aa^2+Cc^2)
    psi <- atan(Cc/Aa)%%pi
    size <- (1/scale)
    rotation <- matrix(c(cos(psi), -sin(psi), sin(psi), cos(psi)), 2, 2)
    A <- B <- C <- D <- numeric(n)
    if (start) {theta <- 0}
    for (i in 1:n) {
        mat <- size*rotation%*%matrix(c(ef$an[i], ef$cn[i], ef$bn[i], ef$dn[i]), 2, 2) %*% matrix(c(cos(i*theta), sin(i*theta), -sin(i*theta), cos(i*theta)), 2, 2)
        A[i] <- mat[1,1]
        B[i] <- mat[1,2]
        C[i] <- mat[2,1]
        D[i] <- mat[2,2]
    }
    list(A=A, B=B, C=C, D=D, size=scale, theta=theta, psi=psi, ao=ef$ao, co=ef$co)
}

# 逆フーリエ変換 (Inverse Fourier Transformation)
# 形状を再構成するために使用する
iefourier <- function(an, bn, cn, dn, k, n , ao=0, co=0) {
    theta <- seq(0, 2*pi, length=n+1)[-(n+1)]
    harmx <- matrix(NA, k, n)
    harmy <- matrix(NA, k, n)
    for(i in 1:k) {
        harmx[i,] <- an[i]*cos(i*theta)+bn[i]*sin(i*theta)
        harmy[i,] <- cn[i]*cos(i*theta)+dn[i]*sin(i*theta)
    }
    x <- (ao/2) + apply(harmx, 2, sum)
    y <- (co/2) + apply(harmy, 2, sum)
    list(x=x, y=y)
}