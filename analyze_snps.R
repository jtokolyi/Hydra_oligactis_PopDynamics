library(ape); library(adegenet); library(poppr); library(vcfR); library(ComplexHeatmap); library(circlize)
library(poppr); library(vcfR); library(gplots); library(stringr); library(GLMMadaptive)

v <- read.vcfR("~/mfs/current/oligactis_RADSeq/analysis/RAD2/results/denovo_final/final.vcf")

x <- read.csv("~/mfs/current/oligactis_RADSeq/analysis/RAD2/M28seasonal_full_dataset2.csv",sep="\t")
x$Strain <- gsub("/","_",x$Strain,fixed=T)

if(FALSE){
####################################################################################################
############################### COLONY #############################################################
dir.create("~/mfs/current/oligactis_RADSeq/analysis/RAD2/results/relatedness");
setwd("~/mfs/current/oligactis_RADSeq/analysis/RAD2/results/relatedness"); 
vdf <- genind2df(vcfR2genind(v), oneColPerAll=T)
vdf <- as.data.frame(lapply(vdf, as.numeric), row.names=row.names(vdf)) + 1

er1 <- 0.01; er2 <- 0.01
colony.base <- "M28seasonal_GQ30"
source("~/mfs/current/oligactis_RADSeq/analysis/RAD1/scripts/write_colony.R")
## Run Colony on the input file colony2.dat manually
}

####################################### Spectrum of genetic diversity ############################################
dist <- dist.gene(genind2df(vcfR2genind(v)), method="perc", pairwise.deletion=T)
tiff("~/mfs/current/oligactis_RADSeq/analysis/RAD2/results/spectrum.tif",res=300, height=1500, width=1500, compression="lzw")
hist(dist,breaks=100,main="Spectrum of genetic diversity",xlab="Genetic distance")
dev.off()

####################################### Clonal diversity #########################################################
col <- read.csv("~/mfs/current/oligactis_RADSeq/analysis/RAD2/results/relatedness/colony_results_GQ30_ER01/M28seasonal_GQ30.BestClone",sep="")
table(col[,2]>0.9) ## MLGs identified with high probability
col <- col[col[,2]>0.9,] ## remove MLGs that were identified with low probability
col[,3] <- gsub("tav","spring",col[,3])
col[,3] <- gsub("osz","autumn",col[,3])
samplings <- c("2018spring","2018autumn","2019spring","2019autumn")

clone.distrib <- NULL

for(i in 1:nrow(col)){
    clone.distrib <- rbind(clone.distrib, sapply(samplings, str_count, string=col[,3][i]))
}

table(apply(clone.distrib!=0,1,sum))
## Distribution of clones in seasons:
## 1  2  3  4
## 35  5  1  1
## Clones that were detected in multiple seasons:
clone.distrib[apply(clone.distrib!=0,1,sum)>1,]
## Distribution of clones that were observed in multiple seasons
##     2018spring 2018autumn 2019spring 2019autumn
##[1,]          1          3          0          0
##[2,]          0         14          7          1
##[3,]          4          2          3          1
##[4,]          4          0          7          0
##[5,]          3          0          5          0
##[6,]          3          0          1          0
##[7,]          1          0          0          1

mlg.sum <- data.frame(sampling=samplings,
                      n.strains = apply(clone.distrib,2,sum),
                      n.mlgs = apply(clone.distrib!=0,2,sum))
mlg.sum$clonal.richness <- (mlg.sum$n.mlgs - 1) / (mlg.sum$n.strains - 1)

#########################################################################################################

v <- read.vcfR("~/mfs/current/oligactis_RADSeq/analysis/RAD2/results/denovo_final/final.vcf")

## reorder samples
sample.names <- colnames(v@gt)[-1]
colnames(v@gt) <- gsub("tav","spring",colnames(v@gt))
colnames(v@gt) <- gsub("osz","autumn",colnames(v@gt))

v2 <- v[samples=unlist(strsplit(col[,3],split=","))]

vord <- v[samples=c(grep("2018tav",sample.names), grep("2018osz",sample.names), grep("2019tav",sample.names), grep("2019osz",sample.names))]



vgi <- vcfR2genind(vord)
pops <- sapply(strsplit(gsub("M28_","",row.names(vgi@tab)),split="_"),"[",1)
pop(vgi) <- factor(pops)

vgi <- repool(seppop(vgi)[c(2,1,4,3)])

###### Plot MSN
dist <- dist.gene(genind2df(vcfR2genind(vord)), method="perc", pairwise.deletion=T)
cols <- col2hex(c("darkolivegreen1","tan","springgreen4","brown")); names(cols) <- popNames(vgi)

tiff("~/mfs/current/oligactis_RADSeq/analysis/RAD2/results/msn.tif",res=300, height=1500, width=1500, compression="lzw")
poppr.msn(vgi, dist, palette=palette(cols),
          vertex.label="",wscale=F,cex=1.2)
dev.off()

######################### Heatmap ####################################
dmat <- as.matrix(dist)
relatedness <- matrix(0, ncol=ncol(dmat), nrow=nrow(dmat), dimnames=dimnames(dmat))
diag(relatedness) <- NA

clones <- read.csv("~/mfs/current/oligactis_RADSeq/analysis/RAD2/results/relatedness/colony_results_GQ30_ER01/M28seasonal_GQ30.BestClone",sep="")[,3]
clones <- strsplit(clones[grep(",",clones)],",")
clone.dyads <- unlist(lapply(clones, combn, 2, simplify=F),recursive=F)
clone.dyads <- data.frame(left=sapply(clone.dyads,"[",1), right=sapply(clone.dyads,"[",2))

clone.dyads[,1] <- gsub("osz","autumn",clone.dyads[,1])
clone.dyads[,1] <- gsub("tav","spring",clone.dyads[,1])
clone.dyads[,2] <- gsub("osz","autumn",clone.dyads[,2])
clone.dyads[,2] <- gsub("tav","spring",clone.dyads[,2])
## Next time, remember to use English in sample names...

for(i in 1:nrow(clone.dyads)){
    relatedness[match(clone.dyads[,1][i], colnames(relatedness)), match(clone.dyads[,2][i], row.names(relatedness))] <- 1
    relatedness[match(clone.dyads[,2][i], colnames(relatedness)), match(clone.dyads[,1][i], row.names(relatedness))] <- 1
}

dmat[upper.tri(dmat)] <- relatedness[upper.tri(relatedness)]

col1 <- colorRamp2(c(0,1),c("lightgrey","red"))
col2 <- colorRamp2(rev(c(0, 0.125, 0.25, 0.5)), c("#FFFFFF", "#71ACD6", "#035C96", "#03204C"))

anno <- rep("",nrow(dmat)); anno[c(20,50,80,110)] <- c("2018 spring", "2018 autumn", "2019 spring", "2019 autumn")

tiff("~/mfs/current/oligactis_RADSeq/analysis/RAD2/results/heatmap.tif",res=300, height=1500, width=1500, compression="lzw")
ht <- Heatmap(dmat, cluster_rows=F, cluster_columns=F,name="relatedness",
        show_row_names=F, show_column_names=F,
        row_order=row.names(dmat),column_order=colnames(dmat),
        show_heatmap_legend=F,
        top_annotation=HeatmapAnnotation(foo = anno_text(anno, gp=gpar(fontsize=9,font=4,col="black")),which="column"),
        left_annotation=HeatmapAnnotation(foo = anno_text(anno, gp=gpar(fontsize=9,font=4,col="black")),which="row"),
        cell_fun = function(j, i, x, y, w, h, fill) {
            if(i > j) {
                grid.rect(x, y, w, h, gp = gpar(fill = col2(dmat[i,j]),col=0))
            } else if(j == i) {
                grid.rect(x, y, w, h, gp = gpar(fill = "white", col=0))
            } else {
                grid.rect(x, y, w, h, gp = gpar(fill = col1(dmat[i,j]),col=0))
            }
        })

lgd_list=list(
    Legend(labels=c("Yes","No"), title="Clones:", legend_gp=gpar(fill=c("red", "lightgrey"))),
    Legend(labels=c("0.0","0.125","0.25","0.5"), legend_gp=gpar(fill=col2(c(0.0,0.125,0.25,0.5))),
           title="Genetic distance"))
draw(ht, heatmap_legend_list=lgd_list)

sampling.order <- factor(sapply(strsplit(colnames(dmat),split="_"), "[",2),
                         levels=c("2018spring","2018autumn","2019spring","2019autumn"))
sampling.groups <- c(0, cumsum(table(sampling.order))) / length(sampling.order)

decorate_heatmap_body("relatedness", {
    for(i in 1:length(sampling.groups)) grid.lines(c(sampling.groups[i],sampling.groups[i]),c(0,1), gp = gpar(lwd = 2, col=grey(0)))
    for(i in 1:length(sampling.groups)) grid.lines(c(0,1), 1-c(sampling.groups[i],sampling.groups[i]), gp = gpar(lwd = 2, col=grey(0)))
})
dev.off()
system("eog ~/mfs/current/oligactis_RADSeq/analysis/RAD2/results/heatmap.tif")

strains <- strsplit(col[,3], split=",");
mlg.lists <- data.frame(strain=unlist(strains), mlgID=rep(1:length(strains), times= sapply(strains, length)))
mlg.lists$sampling <- sapply(strsplit(gsub("M28_","",mlg.lists$strain),split="_"),"[",1)
strains.cc <- mlg.lists$strain[!duplicated(paste(mlg.lists$sampling, mlg.lists$mlgID, sep="_"))]

vcc <- vord[samples=strains.cc] ## Clone-corrected vcfR: one strain per MLG per sampling
vccgi <- vcfR2genind(vcc)
pops <- sapply(strsplit(gsub("M28_","",row.names(vccgi@tab)),split="_"),"[",1)
pop(vccgi) <- as.factor(pops)

vccgi <- repool(seppop(vccgi)[c(2,1,4,3)])

############## Add poppr clonal diversity stats that depend on data transformation above ########
vgc <- as.genclone(vgi)

mlg.lists <- mlg.lists[match(row.names(vgc@tab), mlg.lists$strain),]
vgc@mlg <- mlg.lists$mlgID

pp <- poppr(vgc)

mlg.sum$Evenness <- pp$E.5[1:4]
mlg.sum$Shannon <- pp$H[1:4]

################################### PopGen stats ##############################################################
setwd("/home/jtokolyi/mfs/current/oligactis_RADSeq/analysis/RAD2/results")

sep <- seppop(vgi)
sep <- sep[c("2018spring","2018autumn","2019spring","2019autumn")]

popgen.stat <- data.frame(sampling=c("2018spring","2018autumn","2019spring","2019autumn"),
                          N = NA, Hobs = NA, Hexp = NA, Fis = NA,
                          allelic.richness = NA, private.alleles = NA)

results <- list()
basic.stats <- list()
for(i in 1:4){
    results[[i]] <- summary(sep[[i]])
    basic.stats[[i]] <- basic.stats(sep[[i]])
}

AR <- allel.rich(vgi)$mean.richness
PA <- apply(private_alleles(vgi,count.alleles=F),1,sum)

popgen.stat$N <- unlist(sapply(results, "[", "n"))
popgen.stat$Hobs <- sapply((lapply(results, "[[", "Hobs")),mean)
popgen.stat$Hexp <- sapply((lapply(results, "[[", "Hexp")),mean)
popgen.stat$allelic.richness <- AR[c("2018spring","2018autumn","2019spring","2019autumn")]
popgen.stat$private.alleles <- PA[c("2018spring","2018autumn","2019spring","2019autumn")]
popgen.stat$Fis <- sapply(lapply(basic.stats,"[[","overall"),"[","Fis")

write.table(popgen.stat, file="popgen_stats.txt",col.names=T, row.names=F, quote=F, sep="\t")

sep.cc <- seppop(vccgi)
sep.cc <- sep.cc[c("2018spring","2018autumn","2019spring","2019autumn")]

popgen.stat.cc <- data.frame(sampling=c("2018spring","2018autumn","2019spring","2019autumn"),
                          N = NA, Hobs = NA, Hexp = NA, Fis = NA,
                          allelic.richness = NA, private.alleles = NA)

results.cc <- list()
basic.stats.cc <- list()
for(i in 1:4){
    results.cc[[i]] <- summary(sep.cc[[i]])
    basic.stats.cc[[i]] <- basic.stats(sep.cc[[i]])
}

AR.cc <- allel.rich(vccgi)$mean.richness
PA.cc <- apply(private_alleles(vccgi,count.alleles=F),1,sum)

popgen.stat.cc$N <- unlist(sapply(results.cc, "[", "n"))
popgen.stat.cc$Hobs <- sapply((lapply(results.cc, "[[", "Hobs")),mean)
popgen.stat.cc$Hexp <- sapply((lapply(results.cc, "[[", "Hexp")),mean)
popgen.stat.cc$allelic.richness <- AR.cc[c("2018spring","2018autumn","2019spring","2019autumn")]
popgen.stat.cc$private.alleles <- PA.cc[c("2018spring","2018autumn","2019spring","2019autumn")]
popgen.stat.cc$Fis <- sapply(lapply(basic.stats.cc,"[[","overall"),"[","Fis")

write.table(popgen.stat.cc, file="popgen_stats_cc.txt",col.names=T, row.names=F, quote=F, sep="\t")

############################################ DAPC ####################################################

dapc.vccgi <- dapc(vccgi, n.pca=150, n.da=100)
temp <- optim.a.score(dapc.vccgi)
dapc.vccgi <- dapc(vccgi, n.pca=temp$best,n.da=100)

tiff("~/mfs/current/oligactis_RADSeq/analysis/RAD2/results/dapc.tif",width=1750, height=1750, res=300,compression="lzw")
scatter(dapc.vccgi,clab=0,leg=TRUE,posi.da="bottomright",txt.leg=c("2018 Spring","2018 Autumn","2019 Spring","2019 Autumn"),
        col=col2hex(c("darkolivegreen1","tan","springgreen4","brown")),
        pch=19, cex=1.2)
dev.off()
###################################### AMOVA #########################################################
dist.cc <- dist.gene(genind2df(vcfR2genind(vcc)), method="perc", pairwise.deletion=T)

season <- as.factor(substring(gsub("M28_","",attr(dist.cc,"Labels")),1,7))
pegas::amova(dist.cc~season)

######################################## MLG-based stats ###############################################
x <- read.csv("~/mfs/current/oligactis_RADSeq/analysis/RAD2/M28seasonal_full_dataset2.csv",sep="\t")
x$Strain <- gsub("/","_",x$Strain,fixed=T)
x$MLG <- NA

for(i in 1:nrow(x)){
    x$MLG[i] <- ifelse(any(grepl(x$Strain[i], col[,3])), grep(x$Strain[i], col[,3]), NA)
}

x2 <- x[which(!is.na(x$MLG)),]
x2$ReprMode <- ifelse(x2$Sex=="ASEX",0,1)

m1 <- mixed_model(ReprMode~Season+PolypAge, random=~1|MLG, family="binomial", data=x2)

vars <- get_variance(m1)
## Proportion of variance explained by fixed effects
r2_marginal <- vars$var.fixed / (vars$var.fixed + vars$var.random + vars$var.residual) 
## Proportion of variance explained by fixed + random effects
r2_conditional <- (vars$var.fixed + vars$var.random) / (vars$var.fixed + vars$var.random + vars$var.residual)
## Proportion of variance explained by random effects
r2_random <- vars$var.random / (vars$var.fixed + vars$var.random + vars$var.residual) 

mean(prop.table(table(x2$ReprMode, x2$MLG),2)[2,]) ## average prop. of sexual individuals / MLG
prop.sex <- prop.table(table(x2$ReprMode, x2$MLG),2)[2,]

multi.season <- apply(clone.distrib!=0,1,sum)!=1
kruskal.test(prop.sex, multi.season)

