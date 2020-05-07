#!/usr/local/bin/Rscript

args <- commandArgs(TRUE)
dir <- args[1]
hapfiles <- list.files(path = dir,pattern = "haplo*",full.names = T)

if(length(hapfiles)<2){
  cat("Usage: HaplotectReport.R <compute directory>\n")
  quit()
}

require(ggplot2,quietly = T)
require(dplyr,quietly = T)
require(binom,quietly = T)
require(gridExtra,quietly = T)

contam <- read.csv(hapfiles[1],skip = 1,sep="\t",header=F)[,-c(5,6,8:9)]
contam$V7 <- paste0(round(contam$V7*100,1),"%")

sites <- read.table(hapfiles[2],sep="\t",header=F,stringsAsFactors = F)

sites <- sites %>% rename(chrom = V1,
                          pos1 = V2,
                          pos2 = V3,
                          p1a1 = V4,
                          p1a2= V5,
                          p2a1 = V6,
                          p2a2= V7,
                          haps = V8,
                          dist = V9,
                          reads = V10,
                          observed = V11,
                          contam = V12)
sites$id <- do.call(paste,c(sites[,c("chrom","pos1","pos2")],sep=":"))

sites$minor <- unlist(lapply(strsplit(sites$observed,";"),function(X){ x <- as.numeric(gsub("[ACGT]+:","",X)); x<-sort(x[x!=0]); ifelse(length(x)==3,x[1],mean(x[1:2]))}))
sites <- cbind(sites,binom.confint(sites$minor*2,sites$reads,methods = "exact")[,c(5,6)])
sites$valid <- ifelse(sites$contam>0.02 & sites$minor>1,1,0)
sites$position <- factor(1:nrow(sites),levels = 1:nrow(sites),ordered=T)
point.cols <- c("light gray","red")

p1 <- ggplot() + geom_point(data=sites,aes(y=position,x=contam*100,color=factor(valid)),show.legend = T) + geom_pointrange(data=sites,aes(y=position,x=contam*100,xmin=lower*100,xmax=upper*100,color=factor(valid))) +
  geom_text(data=sites %>% filter(valid==1),aes(y = position,x=upper*100,label=paste0(minor," / ",round(contam * 100,1),"%")),nudge_x = 5,size=2.5,check_overlap = T) + 
  scale_color_manual(values = point.cols,breaks=c(0,1),labels=c("Fail","Pass (>1 read/>2%)"),name="Site evidence") +
  scale_x_continuous(limits = c(0, 50),name="Contamination estimate [%]",sec.axis = sec_axis(trans=~.,name = "Contamination estimate [%]")) +
  scale_y_discrete(name="Locus",breaks = 1:nrow(sites),labels = sites$id,limits = rev(levels(sites$position))) +
  theme_minimal() + theme(axis.text.y = element_text(size=rel(0.75),color = "black"),axis.text.x = element_text(size=rel(1),color = "black"),legend.position = c(.8,.8), legend.background = element_rect(fill="white",color = "black",size = .2))

p2 <- ggplot() + geom_bar(data=sites,aes(y=position,x=reads,fill=reads<100),stat = "identity",show.legend = T) +
  scale_y_discrete(name="Locus",breaks = 1:nrow(sites),labels = sites$id,limits = rev(levels(sites$position))) +
  scale_x_continuous(expand = expansion(mult=0),name="Site coverage",sec.axis = sec_axis(trans=~.,name = "Site coverage")) +
  scale_fill_manual(values = c("dark gray","light gray"),labels=c(">100","<100"),name="Coverage",guide = guide_legend(direction = "vertical", title.position = "top",label.position="top")) + 
  theme_minimal() + theme(axis.text.y=element_blank(),axis.title.y = element_blank())

contam <- cbind(contam[,1:2],nrow(sites),table(sites$valid)[[2]],contam[,4:5])
colnames(contam) <- c("Accession","Total sites","Informative sites","Passing sites","Mean coverage @ sites","Contamination estimate")
contam <- t(contam)

tbl <- tableGrob(contam,rows = rownames(contam),theme=ttheme_default(base_size = 10,padding = unit(c(1,.5),"lines")))
if(as.numeric(contam[4,1])>1){
  tbl$grobs[16][[1]]$label <- paste("**",tbl$grobs[16][[1]]$label,"**")
  tbl$grobs[22][[1]]$gp$fill <- alpha("red",.6)
}

if(as.numeric(gsub("%","",contam[6,1]))>2){
  tbl$grobs[18][[1]]$label <- paste("**",tbl$grobs[18][[1]]$label,"**")
  tbl$grobs[24][[1]]$gp$fill <- alpha("red",.6)
}

pdf(file=paste0("~/",contam[1,1],".contamination_report.pdf"),width = 8.5,height=nrow(sites)*.16)
grid.arrange(tbl,p1,p2,layout_matrix=matrix(c(1,1,2,3),nrow=2,byrow = T),widths=c(1,.5),heights=c(.12,nrow(sites)*.012))
dev.off()
