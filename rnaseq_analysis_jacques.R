library(gplots)
library(RColorBrewer)
library(DESeq)

remove_version_nb = function(x){strsplit(x,split=".",fixed=T)[[1]][1]}
### Nb total de reads mapped sur refseq, pour la normalisation 
tot=c("d0"=23058037,"d14F"=23672987,"d2"=22508391)*1e-9

makeObject <- function(d0,d2,d14,th=2)
{
rownames(d0)=sapply(rownames(d0),FUN=remove_version_nb)
rownames(d2)=sapply(rownames(d2),FUN=remove_version_nb)
rownames(d14)=sapply(rownames(d14),FUN=remove_version_nb)
common=intersect(rownames(d0),
  intersect(rownames(d14),rownames(d2)))
allnames=union(rownames(d0),
  union(rownames(d14),rownames(d2)))

common=intersect(rownames(rpkm.d0),
  intersect(rownames(rpkm.d14),rownames(rpkm.d2)))
allnames=union(rownames(rpkm.d0),
  union(rownames(rpkm.d14),rownames(rpkm.d2)))

res.data=data.frame(
  "name"=common,
  "rpkm.d0"=rpkm.d0[common,],
  "rpkm.d2"=rpkm.d2[common,],
  "rpkm.d14"=rpkm.d14[common,])
res.data$rpkm.d0[is.na(res.data$rpkm.d0)]=0
res.data$rpkm.d2[is.na(res.data$rpkm.d2)]=0
res.data$rpkm.d14[is.na(res.data$rpkm.d14)]=0
res.data=cbind(res.data,
  data.frame(ratio2_0=res.data$rpkm.d2/res.data$rpkm.d0,
             ratio14_0=res.data$rpkm.d14/res.data$rpkm.d0,ratio14_2=res.data$rpkm.d14/res.data$rpkm.d2))

return(res.data)
}

my_nbinomTest <- function(cds,condA,condB,eps=1e-04)
{
    rvf <- cds@rawVarFuncs[[cds@rawVarFuncTable[1]]]
    colA <- conditions(cds) == condA
    colB <- conditions(cds) == condB
    bmv <- getBaseMeansAndVariances(counts(cds)[, colA | colB],
        sizeFactors(cds)[colA | colB])
    rvfA <- rvf
    rvfB <- rvf
    rawScvA <- rvfA(bmv$baseMean)/bmv$baseMean^2
    rawScvB <- rvfB(bmv$baseMean)/bmv$baseMean^2
    rawScvA <- adjustScvForBias(rawScvA, attr(rawScvA, "size"))
    rawScvB <- adjustScvForBias(rawScvB, attr(rawScvB, "size"))
    pval <- nbinomTestForMatrices(counts(cds)[, colA], counts(cds)[,
        colB], sizeFactors(cds)[colA], sizeFactors(cds)[colB],
        rawScvA, rawScvB, eps)
        bmvA <- getBaseMeansAndVariances(counts(cds)[, colA],
            sizeFactors(cds)[colA])
        bmvB <- getBaseMeansAndVariances(counts(cds)[, colB],
            sizeFactors(cds)[colB])
        result=data.frame(id = rownames(counts(cds)), baseMean = bmv$baseMean,
            baseMeanA = bmvA$baseMean, baseMeanB = bmvB$baseMean,
            foldChange = bmvB$baseMean/bmvA$baseMean, log2FoldChange = log2(bmvB$baseMean/bmvA$baseMean),
            pval = pval, padj = p.adjust(pval, method = "BH"),
            resVarA = bmvA$baseVar/(bmvA$baseMean * sum(1/sizeFactors(cds)[colA])/length(condA) +
                rvf(bmv$baseMean)), resVarB = bmvB$baseVar/(bmvB$baseMean *
                sum(1/sizeFactors(cds)[colB])/length(condB) +
                rvf(bmv$baseMean)), stringsAsFactors = FALSE)
return(result) 
}   
    

##############################################################
### results vs. reprogramming.fna
len=read.table("reprogramming_cdna_len.txt",row.names=1)   #read lengths
d0.data=read.table("d0.rp.cnt",row.names=1) #counts per reprogramming_cdna
rpkm.d0=d0.data/len[row.names(d0.data),1]/tot["d0"]
d2.data=read.table("d2.rp.cnt",row.names=1) #counts per reprogramming_cdna
rpkm.d2=d2.data/len[row.names(d2.data),1]/tot["d2"]
d14.data=read.table("d14-F.rp.cnt",row.names=1) #counts per reprogramming_cdna
rpkm.d14=d14.data/len[row.names(d14.data),1]/tot["d14F"]

repro.data = makeObject(rpkm.d0,rpkm.d2,rpkm.d14)


##############################################################
len=read.table("../rs_mus-re_len.txt",row.names=1)   #read lengths
d0.data=read.table("d0.rf.cnt",row.names=1) #counts per transcript
rpkm.d0=d0.data/len[row.names(d0.data),1]/tot["d0"]
d2.data=read.table("d2.rf.cnt",row.names=1) #counts per transcript
rpkm.d2=d2.data/len[row.names(d2.data),1]/tot["d2"]
d14.data=read.table("d14-F.rf.cnt",row.names=1) #counts per transcript
rpkm.d14=d14.data/len[row.names(d14.data),1]/tot["d14F"]

refseq.data = makeObject(rpkm.d0,rpkm.d2,rpkm.d14)
venn(refseq.data[,7:9])


#refseq data
countTable <- round(refseq.data[,1:3])
cds <- newCountDataSet(countTable, as.factor(c("d0","d2","d14")))
cds <- estimateSizeFactors( cds )
cds <- estimateVarianceFunctions( cds,method="blind" )

res_refseq_d2_d0 <- my_nbinomTest(cds,"d2","d0")
res_refseq_d14_d0 <- my_nbinomTest(cds,"d14","d0")
res_refseq_d14_d2 <- my_nbinomTest(cds,"d14","d2")

plotVulcano(myData=res_refseq_d2_d0,myTitle="refseq\nd2 vs. d0")
plotVulcano(myData=res_refseq_d14_d0,myTitle="refseq\nd14 vs. d0")
plotVulcano(myData=res_refseq_d14_d2,myTitle="refseq\nd14 vs. d2")


plotVulcano <- function(myData,myTitle,myCol="blue")
{
plot(myData$log2FoldChange,-log10(myData$pval),cex=2,pch='.',col="grey",xlab="fold change (log2)",ylab="pvalue (-log10)",main=myTitle)
abline(h=-log10(0.05),lty=3,col="red")
abline(v=2,lty=3,col="red")
abline(v=-2,lty=3,col="red")
selectedRefSeq <- which(abs(myData$log2FoldChange)>2 & myData$pval<0.05)
points(myData$log2FoldChange[selectedRefSeq],-log10(myData$pval[selectedRefSeq]),cex=2,pch='.',col=myCol)
}

#Vulcano plots d0 vs. d14:
plot(res_refseq_d14_d0$log2FoldChange,-log10(res_refseq_d14_d0$pval),cex=0.8,pch=4,col="grey",xlab="fold change (log2)",ylab="pvalue (-log10)",main="refseq\nd14 vs. d0")
abline(h=-log10(0.05),lty=3,col="red")
abline(v=2,lty=3,col="red")
abline(v=-2,lty=3,col="red")
selectedRefSeq_d14_d0 <- which(abs(res_refseq_d14_d0$log2FoldChange)>2 & res_refseq_d14_d0$pval<0.05)
points(res_refseq_d14_d0$log2FoldChange[selectedRefSeq_d14_d0],-log10(res_refseq_d14_d0$pval[selectedRefSeq_d14_d0]),cex=0.8,pch=4,col="blue")


names_selectedRefSeq_d14_d0 <- res_refseq_d14_d0[selectedRefSeq_d14_d0,1]
# MA plot d0 vs. d14
ampl=log2(refseq.data$rpkm.d14*refseq.data$rpkm.d0)/2
ratios <- refseq.data$ratio14_2
plot(ampl,ratios,pch='.',main="MA plot - refseq\nd0 and d14")
ampl=log2(refseq.data[names_selectedRefSeq_d14_d0,3]*refseq.data[names_selectedRefSeq_d14_d0,1])/2
ratios <- refseq.data[names_selectedRefSeq_d14_d0,5]
points(ampl,ratios,col="red",pch='.',cex=2)

# plot d0 vs d14
plot(refseq.data[,1],refseq.data[,3],pch='.',cex=2,main="refseq\td0 vs d14",xlab="d0",ylab="d14")
points(refseq.data[names_selectedRefSeq_d14_d0,1],refseq.data[names_selectedRefSeq_d14_d0,3],col="red",pch='.',cex=2)
abline(a=0,b=1)


##############################################################
len=read.table("RetroSeq_len.txt",row.names=1)   #read lengths
d0.data=read.table("d0.rt.cnt",row.names=1) #counts per retro elts
rpkm.d0=d0.data/len[row.names(d0.data),1]/tot["d0"]
d2.data=read.table("d2.rt.cnt",row.names=1) #counts per retro elts
rpkm.d2=d2.data/len[row.names(d2.data),1]/tot["d2"]  
d14.data=read.table("d14-F.rt.cnt",row.names=1) #counts per retro elts
rpkm.d14=d14.data/len[row.names(d14.data),1]/tot["d14F"]

retro.data = makeObject(rpkm.d0,rpkm.d2,rpkm.d14)
venn(retro.data[,7:9])
myColors=brewer.pal(12, "Set3")
plot(seq(1,3,1),retro.data[1,1:3],type="l",ylim=c(0,75),col=myColors[1])
for(i in c(3:12)){lines(seq(1,3,1),retro.data[i,1:3],col=myColors[i])}

#retro data
countTable <- round(retro.data[,1:3])
cds <- newCountDataSet(countTable, as.factor(c("d0","d2","d14")))
cds <- estimateSizeFactors( cds )
cds <- estimateVarianceFunctions( cds,method="blind" )

res_retro_d2_d0 <- my_nbinomTest(cds,"d2","d0")
res_retro_d14_d0 <- my_nbinomTest(cds,"d14","d0")
res_retro_d14_d2 <- my_nbinomTest(cds,"d14","d2")



##############################################################
len=read.table("../repbase.rodent_len.txt",row.names=1)   #read lengths
d0.data=read.table("d0.rb.cnt",row.names=1) #counts per transcript
rpkm.d0=d0.data/len[row.names(d0.data),1]/tot["d0"]
d2.data=read.table("d2.rb.cnt",row.names=1) #counts per transcript
rpkm.d2=d2.data/len[row.names(d2.data),1]/tot["d2"]
d14.data=read.table("d14-F.rb.cnt",row.names=1) #counts per transcript
rpkm.d14=d14.data/len[row.names(d14.data),1]/tot["d14F"]

repbase.data = makeObject(rpkm.d0,rpkm.d2,rpkm.d14)
venn(repbase.data[,7:9])

#repbase data
countTable <- round(repbase.data[,1:3])
cds <- newCountDataSet(countTable, as.factor(c("d0","d2","d14")))
cds <- estimateSizeFactors( cds )
cds <- estimateVarianceFunctions( cds,method="blind" )

res_repbase_d2_d0 <- my_nbinomTest(cds,"d2","d0")
res_repbase_d14_d0 <- my_nbinomTest(cds,"d14","d0")
res_repbase_d14_d2 <- my_nbinomTest(cds,"d14","d2")

plotVulcano(myData=res_repbase_d2_d0,myTitle="repbase\nd2 vs. d0")
plotVulcano(myData=res_repbase_d14_d0,myTitle="repbase\nd14 vs. d0")
plotVulcano(myData=res_repbase_d14_d2,myTitle="repbase\nd14 vs. d2")


############################################################################################
#*******************************************************************************************
# !! previous version: considering only the common elts
#*******************************************************************************************
#*******************************************************************************************
############################################################################################

##############################################################
### results vs. reprogramming.fna
len=read.table("reprogramming_cdna_len.txt",row.names=1)   #read lengths
d0.data=read.table("d0.rp.cnt",row.names=1) #counts per reprogramming_cdna
rpkm.d0=d0.data/len[row.names(d0.data),1]/tot["d0"]
d2.data=read.table("d2.rp.cnt",row.names=1) #counts per reprogramming_cdna
rpkm.d2=d2.data/len[row.names(d2.data),1]/tot["d2"]
d14.data=read.table("d14-F.rp.cnt",row.names=1) #counts per reprogramming_cdna
rpkm.d14=d14.data/len[row.names(d14.data),1]/tot["d14F"]

rownames(rpkm.d0)=sapply(rownames(rpkm.d0),FUN=remove_version_nb)
rownames(rpkm.d2)=sapply(rownames(rpkm.d2),FUN=remove_version_nb)
rownames(rpkm.d14)=sapply(rownames(rpkm.d14),FUN=remove_version_nb)
common=intersect(rownames(rpkm.d0),
  intersect(rownames(rpkm.d14),rownames(rpkm.d2)))
allnames=union(rownames(rpkm.d0),
  union(rownames(rpkm.d14),rownames(rpkm.d2)))


repro.data=data.frame(
  "repro_dna"=common,
  "rpkm.d0"=rpkm.d0[common,],
  "rpkm.d2"=rpkm.d2[common,],
  "rpkm.d14"=rpkm.d14[common,])
repro.data$rpkm.d0[is.na(repro.data$rpkm.d0)]=0
repro.data$rpkm.d2[is.na(repro.data$rpkm.d2)]=0
repro.data$rpkm.d14[is.na(repro.data$rpkm.d14)]=0
repro.data=cbind(repro.data,
  data.frame(ratio2_0=repro.data$rpkm.d2/repro.data$rpkm.d0,
             ratio14_0=repro.data$rpkm.d14/repro.data$rpkm.d0,ratio14_2=repro.data$rpkm.d14/repro.data$rpkm.d2)

##############################################################
len=read.table("../rs_mus-re_len.txt",row.names=1)   #read lengths
d0.data=read.table("d0.rf.cnt",row.names=1) #counts per transcript
rpkm.d0=d0.data/len[row.names(d0.data),1]/tot["d0"]
d2.data=read.table("d2.rf.cnt",row.names=1) #counts per transcript
rpkm.d2=d2.data/len[row.names(d2.data),1]/tot["d2"]
d14.data=read.table("d14-F.rf.cnt",row.names=1) #counts per transcript
rpkm.d14=d14.data/len[row.names(d14.data),1]/tot["d14F"]
  
rownames(rpkm.d0)=sapply(rownames(rpkm.d0),FUN=remove_version_nb)
rownames(rpkm.d2)=sapply(rownames(rpkm.d2),FUN=remove_version_nb)
rownames(rpkm.d14)=sapply(rownames(rpkm.d14),FUN=remove_version_nb)
common=intersect(rownames(rpkm.d0),
  intersect(rownames(rpkm.d14),rownames(rpkm.d2)))
allnames=union(rownames(rpkm.d0),
  union(rownames(rpkm.d14),rownames(rpkm.d2)))
  
refseq.data=data.frame(
  "refseq_dna"=common,
  "rpkm.d0"=rpkm.d0[common,],
  "rpkm.d2"=rpkm.d2[common,],
  "rpkm.d14"=rpkm.d14[common,])
refseq.data$rpkm.d0[is.na(refseq.data$rpkm.d0)]=0
refseq.data$rpkm.d2[is.na(refseq.data$rpkm.d2)]=0
refseq.data$rpkm.d14[is.na(refseq.data$rpkm.d14)]=0
refseq.data=cbind(refseq.data,
  data.frame(ratio2_0=refseq.data$rpkm.d2/refseq.data$rpkm.d0,
             ratio14_0=refseq.data$rpkm.d14/refseq.data$rpkm.d0,ratio14_2=refseq.data$rpkm.d14/refseq.data$rpkm.d2))


##############################################################
len=read.table("RetroSeq_len.txt",row.names=1)   #read lengths
d0.data=read.table("d0.rt.cnt",row.names=1) #counts per retro elts
rpkm.d0=d0.data/len[row.names(d0.data),1]/tot["d0"]
d2.data=read.table("d2.rt.cnt",row.names=1) #counts per retro elts
rpkm.d2=d2.data/len[row.names(d2.data),1]/tot["d2"]
d14.data=read.table("d14-F.rt.cnt",row.names=1) #counts per retro elts
rpkm.d14=d14.data/len[row.names(d14.data),1]/tot["d14F"]

rownames(rpkm.d0)=sapply(rownames(rpkm.d0),FUN=remove_version_nb)
rownames(rpkm.d2)=sapply(rownames(rpkm.d2),FUN=remove_version_nb)
rownames(rpkm.d14)=sapply(rownames(rpkm.d14),FUN=remove_version_nb)
common=intersect(rownames(rpkm.d0),
  intersect(rownames(rpkm.d14),rownames(rpkm.d2)))
allnames=union(rownames(rpkm.d0),
  union(rownames(rpkm.d14),rownames(rpkm.d2)))

retro.data=data.frame(
  "retro_dna"=common,
  "rpkm.d0"=rpkm.d0[common,],
  "rpkm.d2"=rpkm.d2[common,],
  "rpkm.d14"=rpkm.d14[common,])
retro.data$rpkm.d0[is.na(retro.data$rpkm.d0)]=0
retro.data$rpkm.d2[is.na(retro.data$rpkm.d2)]=0
retro.data$rpkm.d14[is.na(retro.data$rpkm.d14)]=0
retro.data=cbind(retro.data,
  data.frame(ratio2_0=retro.data$rpkm.d2/retro.data$rpkm.d0,
             ratio14_0=retro.data$rpkm.d14/retro.data$rpkm.d0,ratio14_2=retro.data$rpkm.d14/retro.data$rpkm.d2))

##############################################################
len=read.table("../repbase.rodent_len.txt",row.names=1)   #read lengths
d0.data=read.table("d0.rb.cnt",row.names=1) #counts per transcript
rpkm.d0=d0.data/len[row.names(d0.data),1]/tot["d0"]
d2.data=read.table("d2.rb.cnt",row.names=1) #counts per transcript
rpkm.d2=d2.data/len[row.names(d2.data),1]/tot["d2"]
d14.data=read.table("d14-F.rb.cnt",row.names=1) #counts per transcript
rpkm.d14=d14.data/len[row.names(d14.data),1]/tot["d14F"]
common=intersect(rownames(rpkm.d0),
  intersect(rownames(rpkm.d14),rownames(rpkm.d2)))
allnames=union(rownames(rpkm.d0),
  union(rownames(rpkm.d14),rownames(rpkm.d2)))

repbase.data=data.frame(
  "repbase"=common,
  "rpkm.d0"=rpkm.d0[common,],
  "rpkm.d2"=rpkm.d2[common,],
  "rpkm.d14"=rpkm.d14[common,])
repbase.data$rpkm.d0[is.na(repbase.data$rpkm.d0)]=0
repbase.data$rpkm.d2[is.na(repbase.data$rpkm.d2)]=0
repbase.data$rpkm.d14[is.na(repbase.data$rpkm.d14)]=0
repbase.data=cbind(repbase.data,
  data.frame(ratio2_0=repbase.data$rpkm.d2/repbase.data$rpkm.d0,
             ratio14_0=repbase.data$rpkm.d14/repbase.data$rpkm.d0,ratio14_2=repbase.data$rpkm.d14/repbase.data$rpkm.d2))
##############################################################


#*************
# MA-plot d2 vs. d0
#*************
nbin=500
nlag=1
ratio=log2(refseq.data$ratio2_0)
ampl=log2(refseq.data$rpkm.d2*refseq.data$rpkm.d0)/2
retro.ratio=log2(retro.data$ratio2_0)
retro.ampl=log2(retro.data$rpkm.d2*retro.data$rpkm.d0)/2
repbase.ratio=log2(repbase.data$ratio2_0)
repbase.ampl=log2(repbase.data$rpkm.d2*repbase.data$rpkm.d0)/2
xx=quantile(ampl,0:nbin/nbin)
probs=matrix(0,ncol=8,nrow=nbin)
qtls=c(.01,.05,.25,.5,.75,.95,.99)
for (n in seq(along=xx[1:(nbin+1-nlag)])) {
   I=which(ampl >= xx[n] & ampl <= xx[n+nlag])
   probs[n,]=c(mean(xx[n+0:nlag]),quantile(ratio[I],qtls))
}
spl=list(low=smooth.spline(probs[,1],probs[,2],df=3)$fit,
 high=smooth.spline(probs[,1],probs[,8],df=3)$fit)
Iselect=c(which(repbase.ratio > predict(spl$high,repbase.ampl)$y),
 which(repbase.ratio < predict(spl$low,repbase.ampl)$y))
Iselect_retro=c(which(retro.ratio > predict(spl$high,retro.ampl)$y),
 which(retro.ratio < predict(spl$low,retro.ampl)$y))
Iselect_refseq=c(which(ratio > predict(spl$high,ampl)$y), which(ratio < predict(spl$low,ampl)$y))
rnames=repbase.data

#write selected refseq
#toWrite <- round(refseq.data[Iselect_refseq,2:7],3)
toWrite <- cbind(round(refseq.data[Iselect_refseq,2:7],3),round(log2(refseq.data[Iselect_refseq,5:7]),3))
rownames(toWrite) <- refseq.data[Iselect_refseq,1]
colnames(toWrite)[7:9] <- c("log2FC_d2_d0","log2FC_d14_d0","log2FC_d14_d2")
write.table(toWrite,file="d2_overepresented_refseq.txt",quote=FALSE,sep="\t")
refseq.de.d2_d0 <- toWrite



png("MAplot_d2vsd0.png",width=1200,height=1000)
plot(ampl,ratio,pch='.',main="d2 vs. d0")
points(repbase.ampl,repbase.ratio,col='red',pch='.',cex=3)
points(retro.ampl,retro.ratio,col='red',pch=17)

for (n in 2:ncol(probs)) {
   spl=smooth.spline(probs[,1],probs[,n],df=3)
   I=which(spl$x>-8)
   lines(spl$x[I],spl$y[I],col='blue',lwd=1)
   text(spl$x[I[1]],spl$y[I[1]],paste(qtls[n-1]*100,"%",sep=''),col='blue',pos=2)
}
labels=gsub("_M*","",unlist(strsplit(x=as.character(rnames[Iselect]),"|",fixed=T))[3*(1:length(Iselect))-1])
labels[5]="MMVL30" #name too long, take 1st column instead
xlab=repbase.ampl[Iselect]
ylab=repbase.ratio[Iselect]
pos=4
for (i in order(labels)) {
   text(x=xlab[i],y=ylab[i],lab=labels[i],col='red',pos=pos,offset=.2,cex=1)
   pos=pos+1
   if (pos>4) pos=1
}
dev.off()

#*****
ntot=length(ratio)
prob=sapply(retro.ratio,FUN=function(x) length(which(ratio>x))/ntot)

png("d2_vs_d0_expression.png",width=800,height=500)
par(lwd=1.5,cex=1.5)
hist(ratio,breaks=100,main="d2 overexpression",col='blue',xlab="log2(d2/d0)")
ypos=seq(3000,100,length.out=length(retro.ratio))
O=order(retro.ratio)
segments(x0=retro.ratio[O],y0=0,x1=retro.ratio[O],y1=ypos,col='red',lty=2)
segments(x0=retro.ratio[O],y0=ypos,x1=5.5,y1=ypos,col='red',lty=2)
text(x=5.5,y=ypos,lab=paste(retro.data$name,", p=",round(prob,3),sep='')[O],col='red',cex=.6,pos=4,offset=0)
dev.off()


#*************
# MA-plot d14 vs. d0
#*************
nbin=500
nlag=1
ratio=log2(refseq.data$ratio14_0)
ampl=log2(refseq.data$rpkm.d14*refseq.data$rpkm.d0)/2
retro.ratio=log2(retro.data$ratio14_0)
retro.ampl=log2(retro.data$rpkm.d14*retro.data$rpkm.d0)/2
repbase.ratio=log2(repbase.data$ratio14_0)
repbase.ampl=log2(repbase.data$rpkm.d14*repbase.data$rpkm.d0)/2
xx=quantile(ampl,0:nbin/nbin)
probs=matrix(0,ncol=8,nrow=nbin)
qtls=c(.01,.05,.25,.5,.75,.95,.99)
for (n in seq(along=xx[1:(nbin+1-nlag)])) {
   I=which(ampl >= xx[n] & ampl <= xx[n+nlag])
   probs[n,]=c(mean(xx[n+0:nlag]),quantile(ratio[I],qtls))
}
spl=list(low=smooth.spline(probs[,1],probs[,2],df=3)$fit,
 high=smooth.spline(probs[,1],probs[,8],df=3)$fit)
Iselect=c(which(repbase.ratio > predict(spl$high,repbase.ampl)$y),
 which(repbase.ratio < predict(spl$low,repbase.ampl)$y))
rnames=repbase.data
Iselect_retro_d14_d0=c(which(retro.ratio > predict(spl$high,retro.ampl)$y),
 which(retro.ratio < predict(spl$low,retro.ampl)$y))
Iselect_refseq=c(which(ratio > predict(spl$high,ampl)$y), which(ratio < predict(spl$low,ampl)$y))

#write selected refseq
toWrite <- cbind(round(refseq.data[Iselect_refseq,2:7],3),round(log2(refseq.data[Iselect_refseq,5:7]),3))
rownames(toWrite) <- refseq.data[Iselect_refseq,1]
colnames(toWrite)[7:9] <- c("log2FC_d2_d0","log2FC_d14_d0","log2FC_d14_d2")
write.table(toWrite,file="d14_overepresented_refseq.txt",quote=FALSE,sep="\t")
refseq.de.d14_d0 <- toWrite


png("MAplot_d14vsd0.png",width=1200,height=1000)
plot(ampl,ratio,pch='.',main="d14 vs. d0")
points(repbase.ampl,repbase.ratio,col='red',pch='.',cex=3)
points(retro.ampl,retro.ratio,col='red',pch=17)

for (n in 2:ncol(probs)) {
   spl=smooth.spline(probs[,1],probs[,n],df=3)
   I=which(spl$x>-8)
   lines(spl$x[I],spl$y[I],col='blue',lwd=1)
   text(spl$x[I[1]],spl$y[I[1]],paste(qtls[n-1]*100,"%",sep=''),col='blue',pos=2)
}
labels=gsub("_M*","",unlist(strsplit(x=as.character(rnames[Iselect]),"|",fixed=T))[3*(1:length(Iselect))-1])
xlab=repbase.ampl[Iselect]
ylab=repbase.ratio[Iselect]
pos=4
for (i in order(labels)) {
   text(x=xlab[i],y=ylab[i],lab=labels[i],col='red',pos=pos,offset=.2,cex=1)
   pos=pos+1
   if (pos>4) pos=1
}
dev.off()

#*********************
ntot=length(ratio)
prob=sapply(retro.ratio,FUN=function(x) length(which(ratio>x))/ntot)

png("d14_vs_d0_expression.png",width=800,height=500)
par(lwd=1.5,cex=1.5)
hist(ratio,breaks=100,main="d14 overexpression",col='blue',xlab="log2(d14/d0)")
ypos=seq(3000,100,length.out=length(retro.ratio))
O=order(retro.ratio)
segments(x0=retro.ratio[O],y0=0,x1=retro.ratio[O],y1=ypos,col='red',lty=2)
segments(x0=retro.ratio[O],y0=ypos,x1=5.5,y1=ypos,col='red',lty=2)
text(x=5.5,y=ypos,lab=paste(retro.data$name,", p=",round(prob,3),sep='')[O],col='red',cex=.6,pos=4,offset=0)
dev.off()

#*************
# MA-plot d14 vs. d2
#*************
nbin=500
nlag=1
ratio=log2(refseq.data$ratio14_2)
ampl=log2(refseq.data$rpkm.d14*refseq.data$rpkm.d2)/2
retro.ratio=log2(retro.data$ratio14_2)
retro.ampl=log2(retro.data$rpkm.d14*retro.data$rpkm.d2)/2
repbase.ratio=log2(repbase.data$ratio14_2)
repbase.ampl=log2(repbase.data$rpkm.d14*repbase.data$rpkm.d2)/2
xx=quantile(ampl,0:nbin/nbin)
probs=matrix(0,ncol=8,nrow=nbin)
qtls=c(.01,.05,.25,.5,.75,.95,.99)
for (n in seq(along=xx[1:(nbin+1-nlag)])) {
   I=which(ampl >= xx[n] & ampl <= xx[n+nlag])
   probs[n,]=c(mean(xx[n+0:nlag]),quantile(ratio[I],qtls))
}
spl=list(low=smooth.spline(probs[,1],probs[,2],df=3)$fit,
 high=smooth.spline(probs[,1],probs[,8],df=3)$fit)
Iselect=c(which(repbase.ratio > predict(spl$high,repbase.ampl)$y),
 which(repbase.ratio < predict(spl$low,repbase.ampl)$y))
rnames=repbase.data
Iselect_retro_d14_d2=c(which(retro.ratio > predict(spl$high,retro.ampl)$y),
 which(retro.ratio < predict(spl$low,retro.ampl)$y))
Iselect_refseq=c(which(ratio > predict(spl$high,ampl)$y), which(ratio < predict(spl$low,ampl)$y))

#write selected refseq
toWrite <- cbind(round(refseq.data[Iselect_refseq,2:7],3),round(log2(refseq.data[Iselect_refseq,5:7]),3))
rownames(toWrite) <- refseq.data[Iselect_refseq,1]
colnames(toWrite)[7:9] <- c("log2FC_d2_d0","log2FC_d14_d0","log2FC_d14_d2")
write.table(toWrite,file="d14_vs_d2_overepresented_refseq.txt",quote=FALSE,sep="\t")
refseq.de.d14_d2 <- toWrite


png("MAplot_d14vsd2.png",width=1200,height=1000)
plot(ampl,ratio,pch='.',main="d14 vs. d2")
points(repbase.ampl,repbase.ratio,col='red',pch='.',cex=3)
points(retro.ampl,retro.ratio,col='red',pch=17)

for (n in 2:ncol(probs)) {
   spl=smooth.spline(probs[,1],probs[,n],df=3)
   I=which(spl$x>-8)
   lines(spl$x[I],spl$y[I],col='blue',lwd=1)
   text(spl$x[I[1]],spl$y[I[1]],paste(qtls[n-1]*100,"%",sep=''),col='blue',pos=2)
}
labels=gsub("_M*","",unlist(strsplit(x=as.character(rnames[Iselect]),"|",fixed=T))[3*(1:length(Iselect))-1])
xlab=repbase.ampl[Iselect]
ylab=repbase.ratio[Iselect]
pos=4
for (i in order(labels)) {
   text(x=xlab[i],y=ylab[i],lab=labels[i],col='red',pos=pos,offset=.2,cex=0.6)
   pos=pos+1
   if (pos>4) pos=1
}
dev.off()


#*********************
ntot=length(ratio)
prob=sapply(retro.ratio,FUN=function(x) length(which(ratio>x))/ntot)

png("d14_vs_d2_expression.png",width=800,height=500)
par(lwd=1.5,cex=1.5)
hist(ratio,breaks=100,main="d14 overexpression\n(vs. d2)",col='blue',xlab="log2(d14/d2)")
ypos=seq(4000,100,length.out=length(retro.ratio))
O=order(retro.ratio)
segments(x0=retro.ratio[O],y0=0,x1=retro.ratio[O],y1=ypos,col='red',lty=2)
segments(x0=retro.ratio[O],y0=ypos,x1=5.5,y1=ypos,col='red',lty=2)
text(x=5.5,y=ypos,lab=paste(retro.data$name,", p=",round(prob,3),sep='')[O],col='red',cex=.6,pos=4,offset=0)
dev.off()


