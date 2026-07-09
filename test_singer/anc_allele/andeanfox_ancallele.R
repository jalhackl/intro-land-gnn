#cd /projects/psg/people/xcj768/arctic_wolves
#module load gsl/2.5 perl bcftools/1.21  
#module load openjdk/20.0.0 gcc/13.2.0 R/4.4.2
# filtering for arg: info sc + ref panel maf, ie no extra maf
#-----------------------------------------------------------------------------------------------------------------------
# 1) get ancestral state from outgroup
# read in outgroup sample from vcf - using bcftools
outgroup='AndeanFox'
in_vcf='impute/joint_impute/2905g_modern/filt/CanFam31/concat/norel/aw_infosc0.8_maf0.01.autosomes.n82.bcf'
outdir='anc_allele'
testsnps<-system(paste("bcftools view -s ",outgroup," ",in_vcf," | bcftools query -f '%CHROM %POS  %REF %ALT [%GT ]\n' ",sep=""),intern=T)
testsnps<-do.call(rbind,strsplit(testsnps,split=" "))  ## ie list merged vertically into df

#> unique(testsnps[,6])
#[1] "0|0" "1|1" "1|0" "0|1"
#-----------------------------------------------------------------------------------------------------------------------
# assign ancestral allele state | outgroup gt
proc_alleles<-function(i){
if(testsnps[i,6] == "0|0") {
allele<- testsnps[i,4]	
}
else if(testsnps[i,6] == "1|1"){
allele<- testsnps[i,5]	
}
else {
allele<-'.'
}
al<-paste("AA=",allele, sep="")
a<-cbind(t(as.data.frame(testsnps[i,])), al)
return(a)}
#-----------------------------------------------------------------------------------------------------------------------
anci<-lapply(1:nrow(testsnps), proc_alleles)
ancinfo<-do.call(rbind, anci)
#> head(ancinfo)
#                                             al    
#testsnps[i, ] "chr1" "1763" "" "T" "G" "0|0" "AA=T"
#testsnps[i, ] "chr1" "2584" "" "T" "A" "0|0" "AA=T"

write.table(ancinfo[,c(1,2,4,5,7)], file=paste(outdir,"/andeanfox_anc",sep=""), 
    quote = FALSE, sep = "\t",row.names = FALSE,col.names = FALSE)
#-----------------------------------------------------------------------------------------------------------------------