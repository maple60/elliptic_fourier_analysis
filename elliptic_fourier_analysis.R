# Setting -----------------------------------------------------------------
rm(list=ls())
if(is.null(dev.list()) == F){
    dev.off()
}
#install.packages("colorhcplot")
#library(colorhcplot)
#library(rgl)
setwd("/Users/konrai/Library/CloudStorage/Dropbox/hawaii_leaf_personal")

#col_trans <- c("#0000FF50", "#FF000050")
#col_trans_crispula <- c("#0000FF50", NA)
#col_trans_serrata <- c(NA, "#FF000050")
#col <- c("#0000FF", "#FF0000")

#d_all <- read.csv("../Quercus_personal/EFA_check_all.csv", comment.char="#", fileEncoding="CP932")
#d_all2 <- read.csv("./Quercus_全体共有/quercus_all_data.csv", fileEncoding="CP932")
# それぞれの葉と種名を対応させる
#d_all$species <- d_all2$species[match(d_all$ID1, d_all2$ID)]
# IDをつなげる
#d_all$ID4 <- paste(d_all$ID1, d_all$ID2, sep="_")
#d_all$ID4[which(d_all$ID3!="")] <- paste(d_all$ID4[which(d_all$ID3!="")], d_all$ID3[which(d_all$ID3!="")], sep="_")

# evaluation==1のみを解析に使う
#d_all_eval1 <- d_all[! (d_all$evaluation==2|d_all$evaluation==3),]

# Functions ------------------------------------------------------
# 楕円フーリエ解析(Elliptic Fourier Analysis)
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

# File selection ----------------------------------------------------------
#site <- select.list(choices= c("鳥取大演習林", "比良山", "宮城"), preselect="鳥取", title="Select site")
#files_cnt <- list.files(path=paste0("./Quercus_画像/調査-", site, "/contour"), pattern="csv")
# files_cntの中からeval==1だけを取り出す。
#for(cnt in files_cnt){
    #cnt <- files_cnt[66]
#    ID <- strsplit(x=paste0(strsplit(x=cnt, split="_")[[1]][2], "_", strsplit(x=cnt, split="_")[[1]][3], "_", strsplit(x=cnt, split="_")[[1]][4]), split="\\.")[[1]][1]
#    if(ID %in% d_all_eval1$ID4==FALSE){
#        files_cnt <- files_cnt[-which(files_cnt==cnt)]
#    }
#}

#Fourier Coefficients & Plot ------------------------------------------------------------
#d_spp <- read.csv("./Quercus_全体共有/quercus_all_data.csv", fileEncoding="CP932")
#d <- read.csv("./contour/2014_leaf_scan/contour_H61_H62_176_4.csv")
folder <- "./contour/2014_leaf_scan"
files <- list.files(folder)
for (file in files){
    id <- sub(".csv", "", sub("contour_", "", file)[[1]][1])
    cat(file, "\n")
    d <- read.csv(paste0(folder, "/", file[1]), header=FALSE)
    plot(d, type="l" , asp=1)
    
    N <- NEF(M=d) # normalized fourier analysis
    cof <- data.frame(A=N$A, B=N$B, C=N$C, D=N$D, size=N$size, theta=N$theta, psi=N$psi)
    write.csv(x=cof, file=paste0("./coefficient/", id, "_coefficient.csv"), row.names=FALSE)
    
    # 進捗報告
    cat(paste(id, "finished : ", which(file==files), "/",length(files)), "\n")
}




for (cnt in files_cnt) {
    #cnt <- files_cnt[1]
    i <- which(cnt==files_cnt)
    ID <- strsplit(x=paste0(strsplit(x=cnt, split="_")[[1]][2], "_", strsplit(x=cnt, split="_")[[1]][3], "_", strsplit(x=cnt, split="_")[[1]][4]), split="\\.")[[1]][1]
    spp <- d_spp$species[which(d_spp$ID==strsplit(x=cnt, split="_")[[1]][2])]

    coor_cnt <- read.csv(file=paste0("./Quercus_画像/調査-", site, "/contour/", cnt), header=FALSE)
    names(coor_cnt) <- c("x", "y")
    
    # プロットしてみよう
    pdf(paste0("./Quercus_全体共有/調査-", site, "/EFD/plot_contour/", ID, "_plot_contour.pdf"))
    plot(coor_cnt$x, coor_cnt$y, type="l", asp=1, xlab="", ylab="", main=ID, frame=FALSE, axes=FALSE)
    polygon(coor_cnt$x, coor_cnt$y, col=ifelse(spp=="Serrata", "#FF000050", "#0000FF50"))
    dev.off()
    
    # Elliptic Fourier Analysis
    N <- NEF(M=coor_cnt)
    cof <- data.frame(A=N$A, B=N$B, C=N$C, D=N$D, size=N$size, theta=N$theta, psi=N$psi)
    write.csv(x=cof, file=paste0("./Quercus_全体共有/調査-", site, "/EFD/coefficient/", ID, "_coefficient.csv") ,row.names=FALSE)
    
    # 進捗報告
    cat(paste(ID, "finished : ", i, "/",length(files_cnt)), "\n")
}

# Harmonic Fourie Power ---------------------------------------------------
files_cof <- list.files(path=paste0("./Quercus_全体共有/調査-", site, "/EFD/coefficient"), pattern="csv")
# files_cntの中からeval==1だけを取り出す
for(cof in files_cof){
    #cof <- files_cof[1]
    i <- which(cof==files_cof)
    ID <- strsplit(x=paste0(strsplit(x=cof, split="_")[[1]][1], "_", strsplit(x=cof, split="_")[[1]][2]), split="\\.")[[1]][1]
    if(ID %in% d_all_eval1$ID4==FALSE){
        files_cof <- files_cof[-which(files_cof==cof)]
    }
}

for (cof in files_cof) {
    #cof <- files_cof[1]
    i <- which(cof==files_cof)
    ID <- strsplit(x=paste0(strsplit(x=cof, split="_")[[1]][1], "_", strsplit(x=cof, split="_")[[1]][2]), split="\\.")[[1]][1]
    d.cof <- read.csv(file=paste0("./Quercus_全体共有/調査-", site, "/EFD/coefficient/", cof))
    d.cof <- d.cof[,1:4] # coefficientsのみを取り出す
    
    co <- d.cof^2
    power <- apply(X=co, MARGIN=1, FUN=sum)/2 # 係数の平方和を2で割ったものがpower
    power_cum <- cumsum(power[-1])/sum(power[-1])
    d.power <- data.frame(power=power, cumulative_power=c(power[-1][1],power_cum))
    write.csv(x=d.power, file=paste0("./Quercus_全体共有/調査-", site, "/EFD/power/", ID, "_power.csv"), row.names=FALSE)
    # 進捗報告
    cat(paste(ID, "finished : ", i, "/",length(files_cof)), "\n")
}


# Coefficients Number ----------------------------------------
files_pow <- list.files(path=paste0("./Quercus_全体共有/調査-", site, "/EFD/power"), pattern="csv")
# files_powの中からeval==1だけを取り出す
for(pow in files_pow){
    #cof <- files_pow[1]
    i <- which(cof==files_pow)
    ID <- strsplit(x=paste0(strsplit(x=pow, split="_")[[1]][1], "_", strsplit(x=pow, split="_")[[1]][2]), split="\\.")[[1]][1]
    if(ID %in% d_all_eval1$ID4==FALSE){
        files_pow <- files_pow[-which(files_pow==pow)]
    }
}

v.over99 <- c()
v.ID <- c()
for (pow in files_pow) {
    #pow <- files_pow[1]
    i <- which(pow==files_pow)
    ID <- strsplit(x=paste0(strsplit(x=pow, split="_")[[1]][1], "_", strsplit(x=pow, split="_")[[1]][2]), split="\\.")[[1]][1]
    d.pow <- read.csv(file=paste0("./Quercus_全体共有/調査-", site, "/EFD/power/", pow))
    v.ID <- c(v.ID, ID)
    v.over99 <- c(v.over99, min(which(d.pow$cumulative_power >= 0.99))) # どこで初めて99%を超える?
    cat(paste(ID, "finished : ", i, "/",length(files_pow)), "\n")
}
d.over99 <- data.frame(ID=v.ID, over99=v.over99)
barplot(height=sort(x=d.over99$over99, decreasing=TRUE))
min(d.over99$over99)
max(d.over99$over99)
# 鳥取：48で十分。50番までのharmonicsを使うことにする
# 滋賀：56で十分。60番までのharmonicsを使うことにする
# 宮城：44で十分らしいので、50番までのharmonicsを使うことにする

# PCA --------------------------------------------
n_cof <- 60 # ここの数は任意。fourier powerで決める
v.ID <- c()
#files_cof <- list.files(path=paste0("./Quercus_全体共有/調査-", site, "/EFD/coefficient"), pattern="csv")
m.cof <- matrix(data=NA, nrow=length(files_cof), ncol=(n_cof+1)*4, byrow=FALSE) # 行にはIDが入り、列には係数が50コずつ入る行列

for (cof in files_cof) {
    #cof <- files_cof[790]
    i <- which(cof==files_cof)
    ID <- sub("_coefficient.csv", "", cof)
    #ID <- strsplit(x=paste0(strsplit(x=cof, split="_")[[1]][1], "_", strsplit(x=cof, split="_")[[1]][2]), split="\\.")[[1]][1]
    d.cof <- read.csv(file=paste0("./Quercus_全体共有/調査-", site, "/EFD/coefficient/", cof))
    d.cof <- d.cof[,1:4] # coefficientsのみを取り出す
    d.cof <- d.cof[1:(n_cof+1),] # 最初のn_cof個の係数をPCAに使う
    m.cof[i,] <- c(d.cof$A, d.cof$B, d.cof$C, d.cof$D)
    v.ID <- c(v.ID, ID)
    # 進捗報告
    cat(paste(ID, "finished : ", i, "/",length(files_cof)), "\n")
}
rownames(m.cof) <- v.ID
d.spp <- read.csv(file="./Quercus_全体共有/quercus_all_data.csv", fileEncoding="CP932")
ID <- substr(x=rownames(m.cof), 1, 5)
ID.spp <- unlist(lapply(ID, grep, x=d.spp$ID))
#d.spp$species[ID.spp]
# PCAをする。A1, B1, C1, D1を省くから、列の指定がすこし複雑
pc <- princomp(m.cof[,c(2:(n_cof+1), 
                        (n_cof + 1)+2:(n_cof+1), 
                        (n_cof+1)*2+2:(n_cof+1), 
                        (n_cof+1)*3+2:(n_cof+1))
                     ]
               )
# 色分け間違っている。修正しましょう。
#plot(pc$scores[,1], pc$scores[,2], asp=1, col=ifelse(d.spp$species=="Serrata", "#FF000050", "#0000FF50"), pch=16)
(pc$sdev^2/sum(pc$sdev^2))[1:5] # 各PCの説明力

#plot(pc$scores[,1], pc$scores[,2], col=ifelse(d.spp$species=="Serrata", NA, "#0000FF50"), pch=16, xlim=c(-0.2,0.1), ylim=c(-0.1,0.1))


# PCA2 --------------------------------------------------------------------
pr.out <- prcomp(m.cof[,c(2:(n_cof+1), 
                          (n_cof+1)+2:(n_cof+1), 
                          (n_cof+1)*2+2:(n_cof+1), 
                          (n_cof+1)*3+2:(n_cof+1))
                       ]
                 ) # center=TRUE, scale=TRUEにしなくても大丈夫だと思う。
#smry_pca <- summary(pr.out)
pve <- pr.out$sdev^2/sum(pr.out$sdev^2) #smry_pca$importance[2,]と同じ
pdf(paste0("./Quercus_全体共有/調査-", site, "/EFD/result/pca_contribution_rate.pdf"))
barplot(pve[1:10]*100, main="Contribution Rate (%)", ylim=c(0,50), names.arg=c(1:10))
legend("topright", legend=paste0("PC", 1:10, ": ", round(pve[1:10]*100, 1), "%"))
dev.off()
#smry_pca$importance[3,] #これは累積寄与率

pr.plot <- pr.out$x[,1:5] # Extract PC1~5
pr.plot.spp <- c()
for(i in row.names(pr.plot)){
    #print(i)
    #i <- row.names(pr.plot)[1]
    ID <-substr(i, 1, 5)
    if(ID=="G0631"){
        ID <- "G0531"
    }
    pr.plot.spp <- c(pr.plot.spp, d_spp$species[which(d_spp$ID==ID)])
}

# プロット
pdf(paste0("./Quercus_全体共有/調査-", site, "/EFD/result/PCA.pdf"))
# All
plot(pr.plot[,1], pr.plot[,2], xlab="PC1", ylab="PC2", col=col_trans[factor(pr.plot.spp)], pch=16, asp=1, main="PCA")
 # Allの左の群
#plot(pr.plot[,1], pr.plot[,2], xlab="PC1", ylab="PC2", col=col_trans[factor(pr.plot.spp)], pch=16, xlim=c(-0.1, 0))
# Allの右の群
#plot(pr.plot[,1], pr.plot[,2], xlab="PC1", ylab="PC2", col=col_trans[factor(pr.plot.spp)], pch=16, xlim=c(0.1, 0.25))

# Only Crispula
plot(pr.plot[,1], pr.plot[,2], xlab="PC1", ylab="PC2", col=col_trans_crispula[factor(pr.plot.spp)], pch=16, asp=1, main="Q. Crispula")
# Only Serrata
plot(pr.plot[,1], pr.plot[,2], xlab="PC1", ylab="PC2", col=col_trans_serrata[factor(pr.plot.spp)], pch=16, asp=1, main="Q. Serrata")
#plot(pr.plot$PC1, pr.plot$PC2, xlab="PC1", ylab="PC2", pch=16, asp=1)
dev.off()
write.csv(pr.plot, paste0("./Quercus_全体共有/調査-", site, "/EFD/result/PCA.csv"))

# Inverse Fourier Analysis(通常飛ばす) ------------------------------------------------
d.PC <- pr.out$x
PC1.mean <- mean(d.PC[,1]) # PC1の平均値。多分0になる
PC2.mean <- mean(d.PC[,2]) # PC2の平均値。多分0になる
PC1.max <- max(d.PC[,1]) # PC1の最大値
PC2.max <- max(d.PC[,2]) # PC2の最大値
PC1.min <- min(d.PC[,1]) # PC1の最大値
PC2.min <- min(d.PC[,2]) # PC2の最大値
PC1.sd <- sd(d.PC[,1]) # =pr.out2$sdev[1] #PC1の標準偏差
PC2.sd <- sd(d.PC[,2]) # =pr.out2$sdev[2] #PC2の標準偏差

v.inv <-c(PC1.min, PC1.mean, PC1.max, PC1.mean+2*PC1.sd, PC1.mean-2*PC1.sd, rep(PC1.mean,5),
          rep(PC2.mean,5), PC2.min, PC2.mean, PC2.max, PC2.mean+2*PC2.sd, PC2.mean-2*PC2.sd)
PC.mean <- apply(pr.out$x, 2, mean)[3:ncol(pr.out$x)] # PC1とPC2以外のPCスコアの平均
v.PC.mean <- rep(PC.mean, each=length(v.inv)/2)
m.inv <- matrix(c(v.inv, v.PC.mean), nrow=length(v.inv)/2, ncol=ncol(pr.out$x), byrow=FALSE)
rownames(m.inv) <- c("PC1.min", "PC1.mean", "PC1.max", "PC1.mean.plus.2SD", "PC1.mean.minus.2SD",
                     "PC2.min", "PC2.mean", "PC2max",  "PC2.mean.plus.2SD", "PC2.mean.minus.2SD") # なくても良い
colnames(m.inv) <- c(paste0(rep("PC", ncol(m.inv)), seq(ncol(m.inv)))) # なくても良い
#m.inv.cof <- m.inv %*% t(pr.out$rotation) # 固有ベクトルの転置行列をかける
m.inv.cof <- pr.out$rotation %*% t(m.inv) # 固有ベクトルの転置行列をかける(こっちが正解)
m.inv.cof <- scale(t(m.inv.cof), center=-pr.out$center, scale = FALSE)
#m.inv.cof <- t(m.inv.cof)

an.inv <- matrix(c(rep(mean(m.cof[,1]),length(v.inv)/2), m.inv.cof[,1:n_cof]), nrow=length(v.inv)/2)
bn.inv <- matrix(c(rep(mean(m.cof[,n_cof+1]),length(v.inv)/2), m.inv.cof[,n_cof+1:n_cof]), nrow=length(v.inv)/2)
cn.inv <- matrix(c(rep(mean(m.cof[,(n_cof+1)*2+1]),length(v.inv)/2), m.inv.cof[,n_cof*2+1:n_cof]), nrow=length(v.inv)/2)
dn.inv <- matrix(c(rep(mean(m.cof[,(n_cof+1)*3+1]),length(v.inv)/2), m.inv.cof[,n_cof*3+1:n_cof]), nrow=length(v.inv)/2)

n <- 300
coor.pc1.min <- iefourier(an=an.inv[1,], bn=bn.inv[1,], cn=cn.inv[1,], dn=dn.inv[1,], k=n_cof, n=n)
coor.pc1.mean <- iefourier(an=an.inv[2,], bn=bn.inv[2,], cn=cn.inv[2,], dn=dn.inv[2,], k=n_cof, n=n)
coor.pc1.max <- iefourier(an=an.inv[3,], bn=bn.inv[3,], cn=cn.inv[3,], dn=dn.inv[3,], k=n_cof, n=n)
coor.pc1.plus.2sd <- iefourier(an=an.inv[4,], bn=bn.inv[4,], cn=cn.inv[4,], dn=dn.inv[4,], k=n_cof, n=n)
coor.pc1.minus.2sd <- iefourier(an=an.inv[5,], bn=bn.inv[5,], cn=cn.inv[5,], dn=dn.inv[5,], k=n_cof, n=n)

coor.pc2.min <- iefourier(an=an.inv[6,], bn=bn.inv[6,], cn=cn.inv[6,], dn=dn.inv[6,], k=n_cof, n=n)
coor.pc2.mean <- iefourier(an=an.inv[7,], bn=bn.inv[7,], cn=cn.inv[7,], dn=dn.inv[7,], k=n_cof, n=n)
coor.pc2.max <- iefourier(an=an.inv[8,], bn=bn.inv[8,], cn=cn.inv[8,], dn=dn.inv[8,], k=n_cof, n=n)
coor.pc2.plus.2sd <- iefourier(an=an.inv[9,], bn=bn.inv[9,], cn=cn.inv[9,], dn=dn.inv[9,], k=n_cof, n=n)
coor.pc2.minus.2sd <- iefourier(an=an.inv[10,], bn=bn.inv[10,], cn=cn.inv[10,], dn=dn.inv[10,], k=n_cof, n=n)

pdf(paste0("./Quercus_全体共有/調査-", site, "/EFD/result/PC1_shift.pdf"))
plot(coor.pc1.mean$x, coor.pc1.mean$y, type="l", asp=1, frame=FALSE, axes=FALSE, xlab="", ylab="", 
     xlim=c(min(c(coor.pc1.mean$x, coor.pc1.min$x, coor.pc1.max$x, coor.pc1.minus.2sd$x, coor.pc1.plus.2sd$x)),
            max(c(coor.pc1.mean$x, coor.pc1.min$x, coor.pc1.max$x, coor.pc1.minus.2sd$x, coor.pc1.plus.2sd$x))
                ),
     ylim=c(min(c(coor.pc1.mean$y, coor.pc1.min$y, coor.pc1.max$y, coor.pc1.minus.2sd$y, coor.pc1.plus.2sd$y)),
            max(c(coor.pc1.mean$y, coor.pc1.min$y, coor.pc1.max$y, coor.pc1.minus.2sd$y, coor.pc1.plus.2sd$y))
            )
     , lwd=3, main="PC1 Shift")
#points(coor.pc1.min$x, -coor.pc1.min$y, type="l", asp=1, col="#0000FF90", lwd=3)
#points(coor.pc1.max$x, coor.pc1.max$y, type="l", asp=1, col="#FF000090", lwd=3)
points(coor.pc1.minus.2sd$x, coor.pc1.minus.2sd$y, type="l", asp=1, col="#0000FF90", lwd=3)
points(coor.pc1.plus.2sd$x, -coor.pc1.plus.2sd$y, type="l", asp=1, col="#FF000090", lwd=3)
legend("topleft", legend=c("Mean", "Mean-2SD", "Mean+2SD"), col=c("black", "#0000FF90", "#FF000090"), lwd=3)
dev.off()

pdf(paste0("./Quercus_全体共有/調査-", site, "/EFD/result/PC2_shift.pdf"))
plot(coor.pc2.mean$x, coor.pc2.mean$y, type="l", asp=1, frame=FALSE, axes=FALSE, xlab="", ylab="", 
     xlim=c(min(c(coor.pc2.mean$x, coor.pc2.min$x, coor.pc2.max$x, coor.pc2.minus.2sd$x, coor.pc2.plus.2sd$x)),
            max(c(coor.pc2.mean$x, coor.pc2.min$x, coor.pc2.max$x, coor.pc2.minus.2sd$x, coor.pc2.plus.2sd$x))
     ),
     ylim=c(min(c(coor.pc2.mean$y, coor.pc2.min$y, coor.pc2.max$y, coor.pc2.minus.2sd$y, coor.pc2.plus.2sd$y)),
            max(c(coor.pc2.mean$y, coor.pc2.min$y, coor.pc2.max$y, coor.pc2.minus.2sd$y, coor.pc2.plus.2sd$y))
     )
     , lwd=3,main="PC2 Shift")
#points(coor.pc2.min$x, coor.pc2.min$y, type="l", asp=1)
#points(coor.pc2.max$x, coor.pc2.max$y, type="l", asp=1)
points(coor.pc2.minus.2sd$x, coor.pc2.minus.2sd$y, type="l", asp=1, col="#0000FF90", lwd=3)
points(coor.pc2.plus.2sd$x, coor.pc2.plus.2sd$y, type="l", asp=1, col="#FF000090", lwd=3)
legend("topleft", legend=c("Mean", "Mean-2SD", "Mean+2SD"), col=c("black", "#0000FF90", "#FF000090"), lwd=3)
dev.off()

# IFA for perticular ID(通常飛ばす) ---------------------------------------------------
ID <- 1
eigenvector <- pr.out$rotation # 固有ベクトル
pcscore <- pr.out$x
inv.cof <- eigenvector %*% pcscore[ID,]

#inv.cof <- scale(inv.cof, center=FALSE, scale = 1/pr.out$scale) # スケールをもとに戻す
# Reverse standardization for mean value 平均値の標準化を元に戻す
inv.cof <- scale(t(inv.cof), center=-pr.out$center, scale = FALSE) # Get Shifted Harmonics!
inv.cof <- t(inv.cof)

a0.inv <- m.cof[ID, 1]
b0.inv <- m.cof[ID, (n_cof+1)*1+1]
c0.inv <- m.cof[ID, (n_cof+1)*2+1]
d0.inv <- m.cof[ID, (n_cof+1)*3+1]
an.inv <- c(a0.inv, inv.cof[1:n_cof])
bn.inv <- c(b0.inv, inv.cof[(n_cof*1+1):(n_cof*2)])
cn.inv <- c(c0.inv, inv.cof[(n_cof*2+1):(n_cof*3)])
dn.inv <- c(d0.inv, inv.cof[(n_cof*3+1):(n_cof*4)])
coor_cnt.inv <- iefourier(an.inv, bn.inv, cn.inv, dn.inv, k=20, n=300)
plot(coor_cnt.inv$x, coor_cnt.inv$y, type="l", asp=1)


# K-means Clustering --------------------------------------------------------------
set.seed(2)
d_k <- data.frame(PC1=pr.plot[,1], PC2=pr.plot[,2])
d_k_withspp <- data.frame(PC1=pr.plot[,1], PC2=pr.plot[,2], spp=pr.plot.spp)
#plot(d_k, xlab="PC1", ylab="PC2", col=col_trans[factor(pr.plot.spp)], pch=16, asp=1)
km.out <- kmeans(x=d_k, centers=2, nstart=20)
km.out
plot(d_k, col=(km.out$cluster+1), main="K-means Clustering(K=2)", xlab="", ylab="", pch=16, asp=1)

pr.plot.left <- d_k_withspp[km.out$cluster==1,]
pr.plot.right <- d_k_withspp[km.out$cluster==2,]
plot(pr.plot.left$PC1, pr.plot.left$PC2, col=col_trans[factor(pr.plot.left$spp)], pch=16)
plot(pr.plot.right$PC1, pr.plot.right$PC2, col=col_trans[factor(pr.plot.right$spp)], pch=16)
left.min <- min(pr.plot.left$PC1)
left.max <- max(pr.plot.left$PC1)
right.min <- min(pr.plot.right$PC1)
right.max <- max(pr.plot.right$PC1)

# 保存用
pdf(paste0("./Quercus_全体共有/調査-", site, "/EFD/result/k-means.pdf"))
plot(d_k, col=(km.out$cluster+1), main="K-means Clustering(K=2)\n(seed=2, nstart=20)", xlab="", ylab="", pch=16, asp=1)
# クラスターの境界
abline(v=left.min, col="#00000080", lty=1)
abline(v=left.max, col="#00000080", lty=1)
abline(v=right.min, col="#00000080", lty=2)
abline(v=right.max, col="#00000080", lty=2)
legend("bottomright", legend=c("Left", "Right"), pch=16, col=c(2,3))
dev.off()


 # Hierarchical Clustering(通常飛ばす) -------------------------------------------------
name <- c()
for(cof in files_cof){
    #cof <- files_cof[1]
    i <- which(cof==files_cof)
    name_tmp <- strsplit(x=paste0(strsplit(x=cof, split="_")[[1]][1], "_", strsplit(x=cof, split="_")[[1]][2]), split="\\.")[[1]][1]
    name <- c(name, name_tmp)
}
name <- sub("_coefficient.csv", "", files_cof)
row.names(d_k) <- name
row.names(d_k_withspp) <- name
pr.plot.left <- d_k_withspp[km.out$cluster==1,]
pr.plot.right <- d_k_withspp[km.out$cluster==2,]

hc.complete <- hclust(dist(d_k), method="complete")
hc.average <- hclust(dist(d_k), method="average")
hc.single <- hclust(dist(d_k), method="single")

hc.complete.left <- hclust(dist(pr.plot.left[,1:2]), method="complete")
hc.average.left <- hclust(dist(pr.plot.left[,1:2]), method="average")
hc.single.left <- hclust(dist(pr.plot.left[,1:2]), method="single")

hc.complete.right <- hclust(dist(pr.plot.right[,1:2]), method="complete")
hc.average.right <- hclust(dist(pr.plot.right[,1:2]), method="average")
hc.single.right <- hclust(dist(pr.plot.right[,1:2]), method="single")

#par(mfrow=c(1,1))
plot(hc.complete, main="Complete Linkage", xlab="", sub="", cex=0.1)
plot(hc.average, main="Average Linkage", xlab="", sub="")
plot(hc.single, main="Single Linkage", xlab="", sub="")

colorhcplot(hc=hc.complete, fac=as.factor(pr.plot.spp),
            hang = -1, main = "Complete Linkage",
            lab.cex = 0.01, lwd = .1, las = 1,
            color = c("red", "blue"))

# left
plot(hc.complete.left, main="Complete Linkage", xlab="", sub="", cex=0.1)
plot(hc.averagehc.complete.left, main="Average Linkage", xlab="", sub="")
plot(hc.singlehc.complete.left, main="Single Linkage", xlab="", sub="")

colorhcplot(hc=hc.complete.left, fac=as.factor(pr.plot.left$spp),
            hang = -1, main = "Complete Linkage",
            lab.cex = 0.1, lwd = .1, las = 1,
            color = c("blue", "red"))


# IDごとに集計する ---------------------------------------------------------------
d_pr.plot <- as.data.frame(pr.plot)
d_pr.plot$spp <- pr.plot.spp

ID1 <- c()
ID2 <- c()
for(i in row.names(d_pr.plot)){
    print(paste(substr(i, 1, 5), substr(i, 7, 8)))
    ID1 <- c(ID1, substr(i, 1, 5))
    ID2 <- c(ID2, substr(i, 7, 8))
}
d_pr.plot$ID <- ID1
d_pr.plot$ID2 <- ID2

pr.plot.summ <- 
    data.frame(ID=aggregate(PC1~ID, data=d_pr.plot, FUN=mean)[,1],
               spp=aggregate(spp~ID, data=d_pr.plot, FUN=unique),
               PC1=aggregate(PC1~ID, data=d_pr.plot, FUN=mean)[,2],
               PC2=aggregate(PC2~ID, data=d_pr.plot, FUN=mean)[,2])

# それぞれのIDごとに個別にプロットする
for(i in unique(d_pr.plot$ID)){
    #i <- unique(d_pr.plot$ID)[1]
    pdf(paste0("./Quercus_全体共有/調査-", site, "/EFD/result_each_pca/", i, "_pca.pdf"))
    plot(d_pr.plot$PC1, d_pr.plot$PC2,
         col=col_trans[factor(d_pr.plot$spp)], pch=16, type="n", xlab="PC1", ylab="PC2", main=i)
    d_tmp <- d_pr.plot[d_pr.plot$ID==i,]
    points(d_tmp$PC1, d_tmp$PC2, col=ifelse(unique(d_tmp$spp)=="Serrata", "#FF000090", "#0000FF90"), pch=d_tmp$ID2)
    points(mean(d_tmp$PC1), mean(d_tmp$PC2), col=ifelse(unique(d_tmp$spp)=="Serrata", "#FF000090", "#0000FF90"), pch=16)
    abline(v=left.min, col="#00000080", lty=1)
    abline(v=left.max, col="#00000080", lty=1)
    abline(v=right.min, col="#00000080", lty=2)
    abline(v=right.max, col="#00000080", lty=2)
    dev.off()
    cat(paste(i, 
              ":",
              which(unique(d_pr.plot$ID)==i), 
              "/", 
              length(unique(d_pr.plot$ID)==i)), "\n")
}

# 同じIDの葉の中でもそれぞれで2群に分かれることがわかった
# IDかつ群に分ける ---------------------------------------------------------------
# 左----
ID1 <- c()
ID2 <- c()
for(i in row.names(pr.plot.left)){
    print(paste(substr(i, 1, 5), substr(i, 7, 8)))
    ID1 <- c(ID1, substr(i, 1, 5))
    ID2 <- c(ID2, substr(i, 7, 8))
}
pr.plot.left$ID <- ID1
pr.plot.left$ID2 <- ID2

pr.plot.summ <- 
    data.frame(ID=aggregate(PC1~ID, data=pr.plot.left, FUN=mean)[,1],
               spp=aggregate(spp~ID, data=pr.plot.left, FUN=unique),
               PC1=aggregate(PC1~ID, data=pr.plot.left, FUN=mean)[,2],
               PC2=aggregate(PC2~ID, data=pr.plot.left, FUN=mean)[,2])

# それぞれのIDごとに個別にプロットする
for(i in unique(pr.plot.left$ID)){
    #i <- unique(pr.plot.left$ID)[1]
    pdf(paste0("./Quercus_全体共有/調査-", site, "/EFD/result_each_pca_left/", i, "_pca.pdf"))
    plot(pr.plot.left$PC1, pr.plot.left$PC2,
         col=col_trans[factor(pr.plot.left$spp)], pch=16, type="n", xlab="PC1", ylab="PC2", main=i)
    d_tmp <- pr.plot.left[pr.plot.left$ID==i,]
    points(d_tmp$PC1, d_tmp$PC2, col=ifelse(unique(d_tmp$spp)=="Serrata", "#FF000090", "#0000FF90"), pch=d_tmp$ID2)
    points(mean(d_tmp$PC1), mean(d_tmp$PC2), col=ifelse(unique(d_tmp$spp)=="Serrata", "#FF000090", "#0000FF90"), pch=16)
    abline(v=left.min, col="#00000080", lty=1)
    abline(v=left.max, col="#00000080", lty=1)
    abline(v=right.min, col="#00000080", lty=2)
    abline(v=right.max, col="#00000080", lty=2)
    dev.off()
    cat(paste(i, 
              ":",
              which(unique(pr.plot.left$ID)==i), 
              "/", 
              length(unique(pr.plot.left$ID)==i)), "\n")
}

# PCAの結果をプロット
pdf(paste0("./Quercus_全体共有/調査-", site, "/EFD/result/pca_left.pdf"))
plot(pr.plot.summ$PC1, pr.plot.summ$PC2, col=col_trans[factor(pr.plot.summ$spp.spp)], pch=16, xlab="PC1", ylab="PC2", xlim=c(left.min, left.max))
abline(v=left.min, col="#00000080", lty=1)
abline(v=left.max, col="#00000080", lty=1)
dev.off()

# 階層的クラスタリング
row.names(pr.plot.summ) <- pr.plot.summ$ID
hc.complete <- hclust(dist(pr.plot.summ[,4:5]), method="complete")
hc.average <- hclust(dist(pr.plot.summ[,4:5]), method="average")
hc.single <- hclust(dist(pr.plot.summ[,4:5]), method="single")

pdf(paste0("./Quercus_全体共有/調査-", site, "/EFD/result/hierarchical_culstering_left.pdf"), width=20)
#plot(hc.complete, main="Complete Linkage", xlab="", sub="", cex=0.1)
colorhcplot(hc=hc.complete, fac=as.factor(pr.plot.summ$spp.spp),
            hang = -1, main = "Complete Linkage",
            lab.cex = 1, lwd = 1, las = 1,
            color = c("blue", "red"))
colorhcplot(hc=hc.average, fac=as.factor(pr.plot.summ$spp.spp),
            hang = -1, main = "Average Linkage",
            lab.cex = 1, lwd = 1, las = 1,
            color = c("blue", "red"))
colorhcplot(hc=hc.single, fac=as.factor(pr.plot.summ$spp.spp),
            hang = -1, main = "Single Linkage",
            lab.cex = 1, lwd = 1, las = 1,
            color = c("blue", "red"))
dev.off()

 # 右----
ID1 <- c()
ID2 <- c()
for(i in row.names(pr.plot.right)){
    print(paste(substr(i, 1, 5), substr(i, 7, 8)))
    ID1 <- c(ID1, substr(i, 1, 5))
    ID2 <- c(ID2, substr(i, 7, 8))
}
pr.plot.right$ID <- ID1
pr.plot.right$ID2 <- ID2

pr.plot.summ <- 
    data.frame(ID=aggregate(PC1~ID, data=pr.plot.right, FUN=mean)[,1],
               spp=aggregate(spp~ID, data=pr.plot.right, FUN=unique),
               PC1=aggregate(PC1~ID, data=pr.plot.right, FUN=mean)[,2],
               PC2=aggregate(PC2~ID, data=pr.plot.right, FUN=mean)[,2])

# それぞれのIDごとに個別にプロットする
for(i in unique(pr.plot.right$ID)){
    #i <- unique(pr.plot.right$ID)[1]
    pdf(paste0("./Quercus_全体共有/調査-", site, "/EFD/result_each_pca_right/", i, "_pca.pdf"))
    plot(pr.plot.right$PC1, pr.plot.right$PC2,
         col=col_trans[factor(pr.plot.right$spp)], pch=16, type="n", xlab="PC1", ylab="PC2", main=i)
    d_tmp <- pr.plot.right[pr.plot.right$ID==i,]
    points(d_tmp$PC1, d_tmp$PC2, col=ifelse(unique(d_tmp$spp)=="Serrata", "#FF000090", "#0000FF90"), pch=d_tmp$ID2)
    points(mean(d_tmp$PC1), mean(d_tmp$PC2), col=ifelse(unique(d_tmp$spp)=="Serrata", "#FF000090", "#0000FF90"), pch=16)
    abline(v=right.min, col="#00000080", lty=1)
    abline(v=right.max, col="#00000080", lty=1)
    abline(v=right.min, col="#00000080", lty=2)
    abline(v=right.max, col="#00000080", lty=2)
    dev.off()
    cat(paste(i, 
              ":",
              which(unique(pr.plot.right$ID)==i), 
              "/", 
              length(unique(pr.plot.right$ID)==i)), "\n")
}

# PCAの結果をプロット
pdf(paste0("./Quercus_全体共有/調査-", site, "/EFD/result/pca_right.pdf"))
plot(pr.plot.summ$PC1, pr.plot.summ$PC2, col=col_trans[factor(pr.plot.summ$spp.spp)], pch=16, xlab="PC1", ylab="PC2", xlim=c(right.min, right.max))
abline(v=right.min, col="#00000080", lty=2)
abline(v=right.max, col="#00000080", lty=2)
dev.off()

# 階層的クラスタリング
row.names(pr.plot.summ) <- pr.plot.summ$ID
hc.complete <- hclust(dist(pr.plot.summ[,4:5]), method="complete")
hc.average <- hclust(dist(pr.plot.summ[,4:5]), method="average")
hc.single <- hclust(dist(pr.plot.summ[,4:5]), method="single")

pdf(paste0("./Quercus_全体共有/調査-", site, "/EFD/result/hierarchical_culstering_right.pdf"), width=20)
#plot(hc.complete, main="Complete Linkage", xlab="", sub="", cex=0.1)
colorhcplot(hc=hc.complete, fac=as.factor(pr.plot.summ$spp.spp),
            hang = -1, main = "Complete Linkage",
            lab.cex = 1, lwd = 1, las = 1,
            color = c("blue", "red"))
colorhcplot(hc=hc.average, fac=as.factor(pr.plot.summ$spp.spp),
            hang = -1, main = "Average Linkage",
            lab.cex = 1, lwd = 1, las = 1,
            color = c("blue", "red"))
colorhcplot(hc=hc.single, fac=as.factor(pr.plot.summ$spp.spp),
            hang = -1, main = "Single Linkage",
            lab.cex = 1, lwd = 1, las = 1,
            color = c("blue", "red"))
dev.off()


# EFD & Elevation ----------------------------------------------------------
# 左----
ID1 <- c()
ID2 <- c()
for(i in row.names(pr.plot.left)){
    print(paste(substr(i, 1, 5), substr(i, 7, 8)))
    ID1 <- c(ID1, substr(i, 1, 5))
    ID2 <- c(ID2, substr(i, 7, 8))
}
pr.plot.left$ID <- ID1
pr.plot.left$ID2 <- ID2

pr.plot.summ <- 
    data.frame(ID=aggregate(PC1~ID, data=pr.plot.left, FUN=mean)[,1],
               spp=aggregate(spp~ID, data=pr.plot.left, FUN=unique),
               PC1=aggregate(PC1~ID, data=pr.plot.left, FUN=mean)[,2],
               PC2=aggregate(PC2~ID, data=pr.plot.left, FUN=mean)[,2])

plot(pr.plot.summ$PC1, pr.plot.summ$PC2, col=col_trans[factor(pr.plot.summ$spp.spp)], pch=16, xlab="PC1", ylab="PC2", xlim=c(left.min, left.max))
abline(v=left.min, col="#00000080", lty=1)
abline(v=left.max, col="#00000080", lty=1)

names(pr.plot.summ) <- c("ID", "ID2", "spp", "PC1", "PC2")
pr.plot.summ$elevation <- d_spp$elevation[match(pr.plot.summ$ID, d_spp$ID)] # 標高データを追加する

plot3d(pr.plot.summ$PC1, pr.plot.summ$PC2, pr.plot.summ$elevation, xlab="PC1", ylab="PC2", zlab="Elevation (m)", col=col_trans[factor(pr.plot.summ$spp)])
rgl.snapshot("../Quercus_personal/3d_shiga.png", fmt="png")

# PC1
pdf("../Quercus_personal/result/pc1_elevation.pdf")
plot(pr.plot.summ$elevation, pr.plot.summ$PC1,
     col=col_trans[factor(pr.plot.summ$spp, levels=c("Crispula", "Serrata"))],
     xlab="Elevation (m)", ylab="PC1", pch=16)
reg <- lm(PC1~elevation, data=pr.plot.summ)
summary(reg)
abline(reg)
dev.off()

# PC2
pdf("../Quercus_personal/result/pc2_elevation.pdf")
plot(pr.plot.summ$elevation, pr.plot.summ$PC2,
     col=col_trans[factor(pr.plot.summ$spp, levels=c("Crispula", "Serrata"))],
     xlab="Elevation (m)", ylab="PC2", pch=16)
reg <- lm(PC2~elevation, data=pr.plot.summ)
summary(reg)
abline(reg)
dev.off()

# 右群を鏡像変換する ---------------------------------------------------------------
for(i in 1:length(row.names(pr.plot.right))){
    #i <- 1
    id <- row.names(pr.plot.right)[i]
    file <- paste0("./Quercus_画像/調査-", site, "/contour/contour_", rownames(pr.plot.right)[i], ".csv")
    
    if(file.exists(file)==TRUE){
        coor_cnt <- read.csv(file=file, header=FALSE)
    } else{
        cat(paste(id, "doesn't exist \n"))
        next
    }
    coor_cnt <- coor_cnt[2:3]
    coor_cnt <- coor_cnt[2:nrow(coor_cnt),]
    names(coor_cnt) <- c("x", "y")
    coor_cnt$x <- as.numeric(coor_cnt$x)
    coor_cnt$y <- as.numeric(coor_cnt$y)
    
    spp <- d_all$species[which(d_all$ID4==id)] # 種名を取得
    plot(coor_cnt$x, coor_cnt$y, type="l", asp=1, xlab="", ylab="", main=id, frame=FALSE, axes=FALSE)
    polygon(coor_cnt$x, coor_cnt$y, col=ifelse(spp=="Serrata", "#FF000050", "#0000FF50"))
    
    # 反転
    A <- matrix(c(-1,0,0,1), nrow=2) # y軸における変換行列
    m_cnt <- as.matrix(coor_cnt)
    trans_cnt <- t(A %*% t(m_cnt))
    trans_cnt[,1] <- trans_cnt[,1]+abs(min(trans_cnt[,1]))+abs(max(trans_cnt[,1]))
    lines(trans_cnt, type="l")
    
    # 保存
    write.table(x=trans_cnt, file=file, col.names=FALSE, row.names=FALSE, sep=",")
    
    cat(paste0(id, ": ", i, " / ", length(row.names(pr.plot.right)),  "\n"))
}
