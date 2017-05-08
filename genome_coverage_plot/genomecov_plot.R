#Jiawei@20160809
#Modified@20160815
library(ggplot2)
library(reshape2)

N_count <- 239850717 # N bases for hg19
args <- commandArgs(TRUE)
if (length(args) == 1){
  max_depth <- as.numeric(args[1])
}else{
  max_depth=100
}
files <- dir('.', pattern = '\\genomeCover\\.txt', full.names = TRUE)

df <- do.call(cbind, lapply(files, function(x){
    label <- unlist(strsplit(basename(x),split="_"))[1]
    df <- read.table(x, sep="\t")
    pdf <- df[df$V1=="genome",]
    pdf[1,3] <- pdf[1,3]-N_count
    pdf[,4] <- pdf[,4]-N_count
    pdf[,5] <- pdf[,3]/pdf[,4]
    pdf.rev <- pdf[order(pdf$V2, decreasing=T),]
    cdf <- cumsum(pdf.rev$V5)
    plot.df <- data.frame(label=cdf[(length(cdf)-max_depth):(length(cdf)-1)])
    names(plot.df) <- label
    return(plot.df)
}))

df <- cbind(depth=c(max_depth:1),df)
df.long <- melt(df,id.vars="depth", measure.vars = colnames(df)[2:(length(files)+1)])
ggplot(df.long, aes(x=df.long$depth, y=df.long$value*100,colour=df.long$variable))+
  geom_line(size=1)+
  scale_colour_manual(values=rainbow(length(files)))+
  xlab("Sequencing Depth")+
  ylab("Base Fraction (%)")+
  ylim(c(0,100))+
  xlim(c(1,max_depth))+
  ggtitle("Coverage Distribution")+
  theme(legend.title=element_blank(),
        text=element_text(size=14))
ggsave("genome_coverage.pdf", width=10, height=6)