library("GO.db")
library("goseq")

setwd("C:/RSYNC/LUMIER_C2H2/ChIP_seq_data/")

assayed.genes <- unique(scan("gencode.v25.protein_coding_gene.txt", what="character"))
de.genes <- unique(scan("ZNF704_SIN3A_co-binding_promoter_gene.txt", what="character"))
gene.vector=as.integer(assayed.genes%in%de.genes)
names(gene.vector)=assayed.genes
length(gene.vector)
length(assayed.genes)
supportedOrganisms()
table(gene.vector)
length(de.genes)

pwf=nullp(gene.vector,"hg19","geneSymbol")
head(pwf)

GO.wall=goseq(pwf,"hg19","geneSymbol")
head(GO.wall)
GO.samp=goseq(pwf,"hg19","geneSymbol",method="Sampling",repcnt=1000)
head(GO.samp)

plot(log10(GO.wall[,2]), log10(GO.samp[match(GO.samp[,1],GO.wall[,1]),2]),
     xlab="log10(Wallenius p-values)",ylab="log10(Sampling p-values)",
     xlim=c(-3,0))
abline(0,1,col=3,lty=2)

GO.nobias=goseq(pwf,"hg19","geneSymbol",method="Hypergeometric")
head(GO.nobias)

plot(log10(GO.wall[,2]), log10(GO.nobias[match(GO.nobias[,1],GO.wall[,1]),2]),
     xlab="log10(Wallenius p-values)", ylab="log10(Hypergeometric p-values)",
     xlim=c(-3,0), ylim=c(-3,0))
abline(0,1,col=3,lty=2)

enriched.GO=GO.nobias$category[p.adjust(GO.nobias$over_represented_pvalue, method="BH")<.05]
head(enriched.GO)


for(go in enriched.GO[1:10]){
   print(GOTERM[[go]])
   cat("--------------------------------------\n")
}

GO.nobias.BP=goseq(pwf,"hg19","geneSymbol",method="Hypergeometric", test.cats=c("GO:BP"))
head(GO.nobias.BP)
enriched.GO.BP=GO.nobias.BP$category[p.adjust(GO.nobias.BP$over_represented_pvalue, method="BH")<.05]
head(enriched.GO.BP)
GO.nobias.MF=goseq(pwf,"hg19","geneSymbol",method="Hypergeometric", test.cats=c("GO:MF"))
head(GO.nobias.MF)
enriched.GO.MF=GO.nobias.MF$category[p.adjust(GO.nobias.MF$over_represented_pvalue, method="BH")<.05]
head(enriched.GO.MF)
GO.nobias.CC=goseq(pwf,"hg19","geneSymbol",method="Hypergeometric", test.cats=c("GO:CC"))
head(GO.nobias.CC)
enriched.GO.CC=GO.nobias.CC$category[p.adjust(GO.nobias.CC$over_represented_pvalue, method="BH")<.05]
head(enriched.GO.CC)

for(go in enriched.GO.MF){
  print(GOTERM[[go]])
  cat("--------------------------------------\n")
}
