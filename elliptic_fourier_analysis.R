# Setting -----------------------------------------------------------------
rm(list=ls())
if(is.null(dev.list()) == F){
    dev.off()
}
setwd("/Users/konrai/Library/CloudStorage/Dropbox/hawaii_leaf_personal")
source("./code/function.R") # 関数の読み込み

#Fourier Coefficients & Plot ------------------------------------------------------------
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
    flush.console()
}

# Harmonic Fourier Power ---------------------------------------------------
files_cof <- list.files(path="./coefficient")

for (file in files_cof) {
    #file <- files_cof[1] # デバッグ用
    #i <- which(cof==files_cof)
    id <- sub(".csv", "", sub("_coefficient", "", file)[[1]][1])
    
    d <- read.csv(file=paste0("./coefficient/", file))
    d <- d[,1:4] # coefficientsのみを取り出す
    
    co <- d^2
    power <- apply(X=co, MARGIN=1, FUN=sum)/2 # 係数の平方和を2で割ったものがpower
    power_cum <- cumsum(power[-1])/sum(power[-1])
    d.power <- data.frame(power=power, cumulative_power=c(power[-1][1],power_cum))
    write.csv(x=d.power, file=paste0("./harmonic_fourier_power/", id, "_power.csv"), row.names=FALSE)
    # 進捗報告
    cat(paste(id, "finished : ", which(file==files_cof), "/",length(files_cof)), "\n")
    flush.console()
}

# Coefficients Number ----------------------------------------
files_pow <- list.files("./harmonic_fourier_power", pattern="csv")

v.over99 <- c()
v.id <- c()
for (file in files_pow) {
    #file <- files_pow[2]
    id <- sub(".csv", "", sub("_power", "", file)[[1]][1])

    d <- read.csv(paste0("./harmonic_fourier_power/", file))
    v.id <- c(v.id, id)
    v.over99 <- c(v.over99, min(which(d$cumulative_power >= 0.99))) # どこで初めて99%を超える?
    cat(paste(id, "finished : ", which(file==files_pow), "/",length(files_pow)), "\n")
    flush.console()
}
d.over99 <- data.frame(ID=v.id, over99=v.over99)
barplot(height=sort(x=d.over99$over99, decreasing=TRUE))
min(d.over99$over99)
max(d.over99$over99)

# PCA --------------------------------------------
n_cof <- 50 # ここの数は任意。fourier powerで決める
v.id <- c()

m.cof <- matrix(data=NA, nrow=length(files_cof), ncol=(n_cof+1)*4, byrow=FALSE) # 行にはIDが入り、列には係数が50コずつ入る行列

for (file in files_cof) {
    #file <- files_cof[1]
    i <- which(file==files_cof)
    id <- sub("_coefficient.csv", "", file)

    d <- read.csv(file=paste0("./coefficient/", file))
    d <- d[1:(n_cof+1),1:4] # coefficientsのみを取り出す
    m.cof[i,] <- c(d$A, d$B, d$C, d$D)
    v.id <- c(v.id, id)
    # 進捗報告
    cat(paste(id, "finished : ", i, "/", length(files_cof)), "\n")
}
rownames(m.cof) <- v.id


# PCAをする。A1, B1, C1, D1を省くから、列の指定がすこし複雑
pc <- princomp(m.cof[,c(2:(n_cof+1), 
                        (n_cof + 1)+2:(n_cof+1), 
                        (n_cof+1)*2+2:(n_cof+1), 
                        (n_cof+1)*3+2:(n_cof+1))
]
)
(pc$sdev^2/sum(pc$sdev^2))[1:5] # 各PCの説明力