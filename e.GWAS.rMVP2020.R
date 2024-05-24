library(rMVP)
library(data.table)


args <- commandArgs(trailingOnly = TRUE)
str(args)
cat(args, sep = "\n")

args[1]
args[2]

setwd("/work/schnablelab/vladimir/NE2020_HISAT/2020_RNA/eGWAS")

phenotype = read.table("mvp.plink.phe",head=TRUE)
gene=which(colnames(phenotype) == args[1] ) #args[1] #Zm00001eb001670
colnames(phenotype[gene])
trait <- phenotype[,c(1, gene)]

colnames(trait)

index <- fread("indexGenes.csv", data.table = F)
name <- index$Name[which(index$ID == colnames(trait)[2])]

print(name)

pc <-  attach.big.matrix("mvp.plink.pc.desc")[]

genotype = attach.big.matrix("mvp.plink.geno.desc")
map <- data.table::fread("mvp.plink.geno.map", data.table = F)
Kin <- attach.big.matrix("mvp.plink.kin.desc")


imMVP <- MVP(
  phe=trait,
  geno=genotype,
  map=map,
  #CV.MLM = pc,
  nPC.MLM=3,  ##if you don't want to add PCs as covariates, please comment out the parameters instead of setting the nPC to 0.
  K = Kin,
  vc.method="EMMA",  ##only works for MLM
  ncpus=16,
  maxLoop=10,
  method=c("MLM"), 
  file.output = F, 
  p.threshold = 0.05/nrow(map)
)
imMVP <- cbind(imMVP$map, imMVP$mlm.results)
data.table::fwrite(imMVP, paste0("results/", name, ".csv.gz"))