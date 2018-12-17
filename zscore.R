#### 
## Create z-score for WES ROH
## Author : Kevin Yauy 2018-12-14
###
library(tools)

args = commandArgs(trailingOnly=TRUE)
filename = tools::file_path_sans_ext(args[1])
filenametoprint = paste(filename,"zscore.tsv", sep="_")

df <- read.delim(args[1],row.names =1)
df$X <- NULL
df.zscore <- scale(df)
#df.zscore <- as.data.frame(df.zscore)
write.table(df.zscore,filenametoprint, sep= "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)

