## Define and vectorize Input sequences
x <- "FIPFSAGPRNCIGQK"
x <- substring(x, 1:nchar(x), 1:nchar(x))
y <- "PFGFGKRSCMGRRLA"
y <- substring(y, 1:nchar(y), 1:nchar(y))

## Gap penality
gp <- 8 # Gap penalty

## Create dynamic programming matrix based on input sequences
ma <- matrix(NA, length(y)+1, length(x)+1, dimnames=list(c("gp", y), c("gp", x)))
ma[1,] <- seq(0, -(length(ma[1,])-1) * gp, -gp)
ma[,1] <- seq(0, -(length(ma[,1])-1) * gp, -gp)
ma

## If desired, write ma to tabular file 
write.table(ma, file="ma.xls", quote=FALSE, na = "", col.names = NA, sep="\t")

