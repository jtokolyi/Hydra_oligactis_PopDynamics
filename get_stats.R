library(doBy)

### Read in read stats from process_radtags.pl output ##########################################################
l1 <- list.files("~/mfs/current/oligactis_RADSeq/data/lib1/samples_truncated/",full.names=T)
l1 <- l1[grep("process_radtags",l1)]
l2 <- list.files("~/mfs/current/oligactis_RADSeq/data/lib2/samples_truncated/",full.names=T)
l2 <- l2[grep("process_radtags",l2)]
l3 <- list.files("~/mfs/current/oligactis_RADSeq/data/lib3/samples_truncated/",full.names=T)
l3 <- l3[grep("process_radtags",l3)]
l4 <- list.files("~/mfs/current/oligactis_RADSeq/data/lib4/samples_truncated/",full.names=T)
l4 <- l4[grep("process_radtags",l4)]

l <- c(l1, l2, l3, l4)

res <- NULL

for(i in 1:length(l)) res <- append(res,readLines(l[i])[13:17])
res <- res[!duplicated(res)]
cat(res, file="reads.csv", sep="\n")
r <- read.csv("reads.csv",sep="\t")

all.samples <- r$Filename[c(grep("osz", r$Filename), grep("tav", r$Filename))]

r <- r[as.character(r$Filename)%in%all.samples,]

sum(r$Total)/2 ## total no. reads
sum(r$Retained)/sum(r$Total) ## prop. reads retained
mean(r$Retained/2) ## mean no. reads / sample
range(r$Retained/2)

### Read in coverage stats from denovo_map.pl output #################################################
cov1 <- readLines("~/mfs/current/oligactis_RADSeq/analysis/RAD2/results/denovo_final/denovo_map.log")
cov1 <- cov1[(grep("Depths of Coverage", cov1)+1):(grep("Depths of Coverage", cov1)+100)]
## Mean/range coverage before removal of PCR duplicates
mean(as.numeric(gsub("x$","",sapply(strsplit(cov1, split=": "), "[", 2))))
range(as.numeric(gsub("x$","",sapply(strsplit(cov1, split=": "), "[", 2))))

cov2 <- readLines("~/mfs/current/oligactis_RADSeq/analysis/RAD2/results/denovo_final/gstacks.log")
## PCR duplication rate
cov2[grep("putative PCR duplicates",cov2)]
## Coverage after romoval of PCR duplicates
cov2[grep("effective per-sample coverage",cov2)]

#cov3 <- readLines("~/mfs/current/oligactis_RADSeq/analysis/RAD2/results/denovo_final/gstacks.log.distribs")
#cov3 <- cov3[grep("^sample",cov3)[1]:(grep("END effective_coverages_per_sample",cov3)-1)]

#cov3 <- t(sapply(lapply(cov3, strsplit, split="\t"),"[[",1))
#colnames(cov3) <- cov3[1,]
#cov3 <- as.data.frame(cov3[-1,])
#cov3$mean_cov <- as.numeric(cov3$mean_cov)

#library(doBy)
#cov3 <- orderBy(~mean_cov, cov3)
