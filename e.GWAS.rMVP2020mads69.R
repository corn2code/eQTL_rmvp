library(rMVP)
library(data.table)


args <- commandArgs(trailingOnly = TRUE)
str(args)
cat(args, sep = "\n")

args[1]
args[2]

setwd("/work/schnablelab/vladimir/NE2020_HISAT/2020_RNA/eGWAS")

phenotype = read.table("mvp.plink.phe",head=TRUE)
phenotypeRAW =fread("../boxcox_transformed_tpm.NE2020.693.csv",data.table = F)
colnames(phenotypeRAW)[1:5]
phenotypeRAW2 <- phenotypeRAW[which(phenotypeRAW$V1 =="Zm00001eb143080"),]
phenotypeRAWt <- as.data.frame(t(phenotypeRAW2))
colnames(phenotypeRAWt) <- "Zm00001eb143080"

phenotype2 <- merge(phenotypeRAWt,phenotype,by.x=0,by.y=1)

gene=which(colnames(phenotype2) == args[1] ) #args[1] "Zm00001eb143080"
colnames(phenotype2[gene])
trait <- phenotype2[,c(1, gene)]

colnames(trait)
colnames(trait)[1] <- colnames(phenotype)[1]

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