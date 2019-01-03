#### 
## Create MAD-score for WES ROH
## Author : Kevin Yauy 2019-01-02
###
library(tools)
library(quantable)
args = commandArgs(trailingOnly=TRUE)
filename = tools::file_path_sans_ext(args[1])
filenametoprint = paste(filename,"MADscore.tsv", sep="_")

df <- read.delim(args[1],row.names =1)
df$X <- NULL
df$chrX <- NULL
df.MAD <- robustscale(df, dim = 2, center=TRUE, scale = TRUE, preserveScale = FALSE)
df.final <- df.MAD$data
df.final$Mean_without_1max <-apply(df.final,1,function(x) mean(x[!x %in% range(x)]))
#df.zscore <- as.data.frame(df.zscore)
write.table(df.final,filenametoprint, sep= "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)

