
#cit_per_ATACTF_RNAmodule
#ATAC_TF.txt  perturbation2.txt  RNA_modules.txt
#!/usr/bin/env Rscript

outfname = paste(seed, "cit.fdr.allsample.pcgene.nperm100.txt")
# permutation times
n.perm = 100

for (RNAmodule in c5olnames(RNA_modules)){
  print(RNAmodule)
  r = RNA_modules[[RNAmodule]]
  for (perturbation_count in 1:dim(perturbation2)[2]){
    print(perturbation_count)
    p = perturbation2[, perturbation_count]
    
    ATACTF_pcc_p = apply(ATAC_TF, 1, function(x) cor.test(x,r)$"p.value")
    ATACTF_pcc_fdr = p.adjust(ATACTF_pcc_p, 'BH')
    subATAC_TF = ATAC_TF[ATACTF_pcc_p < 0.05, ]
    subATAC_TF <- na.omit(subATAC_TF)
    print(dim(subATAC_TF))
    if (dim(subATAC_TF)[1] == 0) {next} # skip this because there will be no result
    
    perm.index = NA # set a new perm.index for each lifestyle
    myresults = list()
    tstcount = 1
    
    for( subATAC_TFcount in 1:dim(subATAC_TF)[1] ){
      a = as.numeric( subATAC_TF[subATAC_TFcount, ] )
      s <- cbind.data.frame(p=p, a=a, r=r)
      s = s[complete.cases(s), ]
      if (sum(is.na(perm.index)) > 0) {
        perm.index = matrix(NA, nrow=dim(s)[1], ncol=n.perm )
        for( j in 1:ncol(perm.index) ) perm.index[, j] = sample( 1:dim(s)[1] )
      }
      myresults[[tstcount]] = cit.cp(s$p, s$a, s$r, perm.index=perm.index, n.perm=n.perm, rseed = seed)
      tstcount = tstcount + 1
      print(tstcount)
    }
    fdrresult = fdr.cit( myresults )
    fdrresult = cbind.data.frame(RNAmodule = RNAmodule, perturbations = colnames(perturbation2)[perturbation_count], ATACTFs = rownames(subATAC_TF), fdrresult)
    #write.table(fdrresult, outfname, sep = "\t", col.names = !file.exists(outfname), row.names = F, quote = F, append = T)
    write.table(fdrresult, "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/GRN/integrate/cit/fdrresult_per_ATAC_RNA.txt",sep = "\t")
    
    }
}

##
#!/usr/bin/env Rscript

install.packages("/lustre2/jdhan_pkuhpc/zhusx/R/causual/cit_per_ATACTF_RNAmodule/cit_2.3.1.tar.gz", repos = NULL)

install.packages("/lustre2/jdhan_pkuhpc/zhusx/R/causual/cit_per_ATACTF_RNAmodule/gtools_3.9.4.tar.gz", repos = NULL)


library(cit)
library(data.table)
library(gtools)
seed = sample(1:1000000, 1)
print(seed)
set.seed(seed)

outfname = paste(seed, "cit.fdr.allsample.pcgene.nperm100.par.txt")

RNA_modules = read.table(file="cit_per_ATACTF_RNAmodule/RNA_modules.txt",sep = "\t")
perturbation2 = read.table(file="cit_per_ATACTF_RNAmodule/perturbation2.txt",sep = "\t")
ATAC_TF = read.table(file="cit_per_ATACTF_RNAmodule/ATAC_TF.txt",sep = "\t")

n.perm = 100

fdrresult3 = data.frame()
for (RNAmodule in colnames(RNA_modules)){
  #RNAmodule = "RNA.module.two1"
  print(RNAmodule)
  r = RNA_modules[[RNAmodule]]
  
  fdrresult2 = data.frame()
  for (perturbation_count in 1:dim(perturbation2)[2]){
    #perturbation_count = 2
    print(perturbation_count)
    p = perturbation2[, perturbation_count]
    
    ATACTF_pcc_p = apply(ATAC_TF, 1, function(x) cor.test(x,r)$"p.value")
    ATACTF_pcc_fdr = p.adjust(ATACTF_pcc_p, 'BH')
    subATAC_TF = ATAC_TF[ATACTF_pcc_p < 0.05, ]
    subATAC_TF <- na.omit(subATAC_TF)
    print(dim(subATAC_TF))
    if (dim(subATAC_TF)[1] == 0) {next} # skip this because there will be no result
    
    perm.index = NA # set a new perm.index for each lifestyle
    myresults = list()
    tstcount = 1
    
    for( subATAC_TFcount in 1:dim(subATAC_TF)[1] ){

      a = as.numeric( subATAC_TF[subATAC_TFcount, ] )
      s <- cbind.data.frame(p=p, a=a, r=r)
      s = s[complete.cases(s), ]
      if (sum(is.na(perm.index)) > 0) {
        perm.index = matrix(NA, nrow=dim(s)[1], ncol=n.perm )
        for( j in 1:ncol(perm.index) ) perm.index[, j] = sample( 1:dim(s)[1] )
      }
      myresults[[tstcount]] = cit.cp(s$p, s$a, s$r, perm.index=perm.index, n.perm=n.perm, rseed = seed)
      tstcount = tstcount + 1
      print(tstcount)
    }
    fdrresult = fdr.cit(myresults)
    fdrresult = cbind.data.frame(RNAmodule = RNAmodule, perturbations = colnames(perturbation2)[perturbation_count], ATACTFs = rownames(subATAC_TF), fdrresult)
    
    fdrresult2 = rbind(fdrresult2,fdrresult)
  }
  name = paste("fdrresult_per_ATAC_RNA",RNAmodule,".txt",sep="")
  write.table(fdrresult2,name,sep = "\t")
  fdrresult3 = rbind(fdrresult3,fdrresult2)
}
write.table(fdrresult3, "fdrresult_per_ATAC_RNAallmodule_mergeresult.txt",sep = "\t")







