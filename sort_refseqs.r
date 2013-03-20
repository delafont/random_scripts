
# Initialize data, remove zero counts, NAs and pseudo-genes from ref.
ref = read.table("RefSeq_MEF.txt", header=TRUE, sep="\t")
mine = read.table("MEFgenes_nosf.csv", header=TRUE, sep="\t")
deseq = read.table("DESeq_result.csv", header=TRUE, sep="\t")
deseq[,"ratio"] = deseq[,5]
mine[,"ratio"] = mine[,2]/mine[,3]
ref[,"ratio"] = ref[,9]/ref[,10]
thrs = 100
ref = ref[which(!is.na(ref$ratio) & ref$cnt.KO>=thrs & ref$cnt.WT>=thrs & !is.na(ref$mgi_symbol)),]
mine = mine[which(!is.na(mine$ratio) & mine$KO>=thrs & mine$WT>=thrs),]
deseq = deseq[which(!is.na(deseq$ratio) & deseq$baseMean>=thrs),]
max(ref$cnt.KO)
max(mine$KO)
min(ref$cnt.WT)
min(mine$WT)

dim(ref)
dim(mine)
dim(deseq)
length(intersect(intersect(ref$mgi_symbol,mine$id),deseq$id))

# Keep only genes which name appears in both files
reref = ref[pmatch(intersect(ref$mgi_symbol,mine$id),ref$mgi_symbol),]
remine = mine[pmatch(intersect(ref$mgi_symbol,mine$id),mine$id),]
redeseq = deseq[pmatch(intersect(ref$mgi_symbol,mine$id),deseq$id),]

# Order lines by gene name
orderbyname = order(reref$mgi_symbol)
reref = reref[orderbyname,]
orderbyname = order(remine$id)
remine = remine[orderbyname,]
orderbyname = order(redeseq$id)
redeseq = redeseq[orderbyname,]

# Plot ratios from both files against each other. Seems well correlated
pdf(file="refratios~mineratios_nosf.pdf")
plot(log(reref[,9]/reref[,10],2)~log(remine$ratio,2), xlab="log2(mine.ratios)", ylab="log2(ref.ratios)")
graphics.off()
pdf(file="refratios~deseq.pdf")
plot(log(reref[,9]/reref[,10],2)~log(remine$ratio,2), xlab="log2(mine.ratios)", ylab="log2(ref.ratios)")
graphics.off()

# Check coherence of values
max(reref$ratio)
max(remine$ratio)
max(redeseq$ratio)
min(reref$ratio)
min(remine$ratio)
min(redeseq$ratio)

# Write ratios from both files in front of each other in a new tab file
refnames = levels(reref$mgi_symbol)[reref$mgi_symbol]
minenames = levels(remine$id)[remine$id]
deseqnames = levels(redeseq$id)[redeseq$id]
all(refnames == minenames) & all(refnames == deseqnames)
write.table(cbind(minenames, reref$ratio, remine$ratio, redeseq$ratio),
             "MEF_sortedbyname", sep="\t", row.names=FALSE)

# Check that most expressed genes correspond across files
orderbyratio = order(reref$ratio, decreasing=TRUE)
newref = reref[orderbyratio,]
orderbyratio = order(remine$ratio, decreasing=TRUE)
newmine = remine[orderbyratio,]
orderbyratio = order(redeseq$ratio, decreasing=TRUE)
newdeseq = redeseq[orderbyratio,]
length(intersect(newref$mgi_symbol[1:50],newmine$id[1:50]))
length(intersect(newref$mgi_symbol[1:50],newdeseq$id[1:50]))
length(intersect(newref$mgi_symbol[1:500],newmine$id[1:500]))
length(intersect(newref$mgi_symbol[1:500],newdeseq$id[1:500]))

newrefnames = levels(newref$mgi_symbol)[newref$mgi_symbol]
newminenames = levels(newmine$id)[newmine$id]
write.table(cbind(newrefnames, newminenames, newref$ratio, newmine$ratio, newdeseq$ratio),
             "MEF_sortedbyratios", sep="\t", row.names=FALSE)


## DESeq
##data = read.table("MEFgenes_nosf.csv", header=TRUE, sep="\t")
#data <- read.delim("MEFgenes_nosf.csv", header=TRUE, stringsAsFactors=TRUE )
#data <- data.frame(WT=data[,3],KO=data[,2], row.names=data[,1])
#conds <- colnames(data)
#cds <- newCountDataSet(data,conds)
#cds <- estimateSizeFactors(cds)
#cds <- estimateVarianceFunctions(cds, method='blind')
#res <- nbinomTest(cds,conds[1], conds[2])
#write.table(res, "DESeq_result.csv", sep="\t", row.names=FALSE)


