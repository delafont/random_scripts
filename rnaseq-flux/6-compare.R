#! /usr/bin/env Rscript


options(width=140)
mode = "count"
#mode = "rpkm"


#---------------------------------------------------------------------------------#
# Methods comparison (COUNTS)


# Init
d = read.table("all_counts.txt", header=TRUE,
    colClasses=c("character",rep("numeric",10),"character",rep("numeric",4),"character"))

# Comparison plot
if (1) {
png("Fig2-SIM_1x4.png", height=480, width = 1600)
#pdf("Fig2-SIM_1x4.pdf", height=4, width = 13)
par(cex=1.3, cex.lab=1.5, cex.axis=1.5, mar=c(8,4,8,4)/2, oma=c(0,6,0,0), las=1)
par(mfrow=c(1,4))
par(pty="s", cex.main=1.5)
ylim=c(0,13)
xlim=c(0,13)
lwd=2
pch=20
ylab="SIM count"

comp.plot = function(alt, ref, alt.name) {
    # Correlation only on non-zero elements
    ra = data.frame(ref=ref, alt=alt)
    raz = ra[ra$alt>0,]
    xy.rho = cor(raz$ref,raz$alt, use="pairwise.complete.obs", method="spearman")
    xy.tau = cor(raz$ref,raz$alt, use="pairwise.complete.obs", method="kendall")
    xy.rho = signif(xy.rho, 2)
    xy.tau = signif(xy.tau, 2)
    ra[is.na(ra$alt),] = 0
    y = log10(ra$ref+1)
    x = log10(ra$alt+1)
    plot(x,y, xlim=xlim, ylim=ylim, xlab=paste(alt.name,"count"), pch=pch, xaxt="n", yaxt="n", ylab="",
            main=bquote(paste("N=",.(length(x)),", ",tau,"=",.(xy.tau),", ",rho,"=",.(xy.rho))),
            col=rgb(0,0,0,0.3))
    abline(0,1, col="grey", lt="dashed", lwd=lwd)
    if (alt.name=="Rnacounter") {title(ylab="SIM count")}
    # Log ticks:
    major.ticks <- axTicks(1)
    labels <- sapply(major.ticks,function(i)
                as.expression(bquote(10^ .(i)))
              )
    axis(1,at=major.ticks,labels=labels)
    axis(2,at=major.ticks,labels=labels)
}

if (mode == "rpkm") {
    ref = d$SIMrpkm
    alt = d$NNLSrpkm
    comp.plot(alt,ref, "Rnacounter")
    alt = d$RSEMfpkm
    comp.plot(alt,ref, "RSEM")
    alt = d$CUFFfpkm
    comp.plot(alt,ref, "Cuffquant")
    alt = d$SAILrpkm
    comp.plot(alt,ref, "Sailfish")
} else {
    ref = d$SIMcount
    alt = d$NNLScount
    comp.plot(alt,ref, "Rnacounter")
    alt = d$RSEMcount
    comp.plot(alt,ref, "RSEM")
    alt = d$CUFFcount
    comp.plot(alt,ref, "Cuffquant")
    alt = d$SAILcount
    comp.plot(alt,ref, "Sailfish")
}

dev.off()
}





#---------------------------------------------------------------------------------#
# Understand which transcripts are drifting off the line
# Same for all methods?

if (0) {
    out.plus = list()
    out.minus = list()
    zeroes = list()
    nas = list()
    methods = c("CUFFcount","SAILcount","RSEMcount","NNLScount")
    for (method in methods) {
        alt = d[,method]
        zero = which(alt==0)
        zeroes[[method]] = d[zero, c("GeneName","ENST")]
        ratio = ref/alt
        off = which(ratio > 10)
        out.plus[[method]] = d[off, c("GeneName","ENST")]
        off = which(ratio < 0.1)
        out.minus[[method]] = d[off, c("GeneName","ENST")]
        na = which(is.na(alt))
        nas[[method]] = d[na, c("GeneName","ENST")]
    }
    print(lapply(out.plus, FUN=nrow))
    print(lapply(out.minus, FUN=nrow))
    print(lapply(zeroes, FUN=nrow))
    print(lapply(nas, FUN=nrow))

    m = merge(out.plus[[1]], out.plus[[2]], by="ENST", all=TRUE)
    colnames(m) = c("ENST",methods[1],methods[2])
    m = merge(m, out.plus[[3]], by="ENST", all=TRUE)
    colnames(m) = c("ENST",methods[1],methods[2],methods[3])
    m = merge(m, out.plus[[4]], by="ENST", all=TRUE)
    colnames(m) = c("ENST",methods[1],methods[2],methods[3],methods[4])

    intersections = list()
    comparisons = list()
    for (k in c(2,3,4)) {
        combs = combn(methods,k)
        for (j in 1:ncol(combs)) {
            meths = combs[,j]
            print(meths)
            inter = d$ENST
            for (method in meths) {
                inter = intersect(inter, out.plus[[method]]$ENST)
            }
            comp = paste(meths,collapse='-')
            intersections[[comp]] = inter
            comparisons[[comp]] = m[m$ENST %in% inter, ]
        }
    }
    lapply(intersections, FUN=length)
    common = common[order(comparisons[["CUFFcount-SAILcount-RSEMcount-NNLScount"]]$CUFFcount),]
}


