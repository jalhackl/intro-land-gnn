

#cd  /projects/psg/people/xcj768/arctic_wolves
#module load gsl/2.5 perl bcftools/1.21 openjdk/20.0.0 gcc/13.2.0 R/4.4.2
#-----------------------------------------------------------------------------------------------------------------------
library(dplyr)
vcf = "impute/joint_impute/2905g_modern/filt/CanFam31/concat/norel/aw_infosc0.8_maf0.01.autosomes.n82.nooutgroup.vcf.gz"
samples.in.vcf<-system(paste("bcftools query -l ",vcf," ",sep=""),intern=T) # -l, --list-samples: list sample names and exit
#-----------------------------------------------------------------------------------------------------------------------
inmeta = 'qc/mosdepth/aw.dataset.aut.mosdepth.metadata.23oct25.txt'
metadata<-read.table(inmeta,  fill = TRUE , header = T, sep='\t')

metadata[which(metadata$SampleNumber%in% "AlaskanMalamuteDog"),'SampleNumber']<-"SAMN03168378"
metadata[which(metadata$SampleName%in% "AlaskanMalamuteDog"),'SampleName']<-"SAMN03168378"
#-----------------------------------------------------------------------------------------------------------------------

metadata_vcf <-metadata[(metadata$SampleNumber %in% samples.in.vcf),]
metadata_vcf<-metadata_vcf[ order(match(metadata_vcf$SampleNumber, samples.in.vcf)), ]

metadata_vcf<-metadata_vcf %>% 
mutate(Type = factor(Type, levels = c( "wolf","coyote", "dog")))

# add identifier where wGreenland vs dGreenland
metadata_vcf$GeogRegion<-paste(substr(metadata_vcf$Type, 1, 1),metadata_vcf$Region, sep='')

out<- metadata_vcf[,c('SampleNumber', 'GeogRegion', 'Type')]
colnames(out)<-c('sample', 'population', 'group')
out$sex <-NA

head(out)
#        sample   population group sex
#1  D2301171079   wGreenland  wolf  NA
#2  D2301171058 wNorthBaffin  wolf  NA
#7  D2301171066 wSouthBaffin  wolf  NA
#9  D2301171068 wSouthBaffin  wolf  NA
#10 Arctic_Wolf   wGreenland  wolf  NA
#11 D2301171063      wAlaska  wolf  NA


write.table(out, file=paste("args/relate/input/aw.poplabels",sep=""), 
    quote = FALSE, sep = " ",row.names = FALSE,col.names = TRUE)

#-----------------------------------------------------------------------------------------------------------------------
#	The four columns are:
#Individual ID as specified in the .sample file [string]. Samples must be listed in the same order as the .sample file.
#Population label [string]
#Group label [string]
#Sex [integer]
#Diploid organisms:
#sample population group sex
#UNR1 PJB SAS NA
#UNR2 JPT EAS NA
#UNR3 GBR EUR NA
#UNR4 YRI AFR NA
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# make poplabels for greenland+ellesmere wolves only
#Tue May 19 10:08:21 CEST 2026
poplabels_orig<-read.table("/projects/psg/people/xcj768/arctic_wolves/args/relate/input/aw.poplabels", header = T, sep='')

grell<-poplabels_orig[poplabels_orig$population %in% c("wGreenland", "wEllesmere"),]

nrow(grell) 
#[1] 11

grell[which(grell$population=="wGreenland"),c('population')]<-"GREL"
grell[which(grell$population=="wEllesmere"),c('population')]<-"GREL"

write.table(grell, file=paste("/projects/psg/people/xcj768/arctic_wolves/args/relate/input/grellonly.poplabels",sep=""), 
    quote = FALSE, sep = " ",row.names = FALSE,col.names = TRUE)
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# make poplabels for plots groupign by arcticwolf, mainlandwolf, arcticdog, breeddog, coyote

#> unique(out$population)
# [1] "wGreenland"      "wNorthBaffin"    "wSouthBaffin"    "wAlaska"        
# [5] "wEllesmere"      "wBanksIsland"    "wVictoriaIsland" "wAtlantic"      
# [9] "wCentral"        "cAlaska"         "wAlberta"        "cAlberta"       
#[13] "dGreenland"      "dAlaska"         "dGshep"          "dChowchow"      

out[out$population %in% c("wGreenland", "wNorthBaffin","wSouthBaffin", "wEllesmere", "wBanksIsland",    "wVictoriaIsland"),'population']<-'arcticwolf'
out[out$population %in% c("wAlaska", "wAtlantic","wCentral", "wAlberta"),'population']<-'mainlandwolf'
out[out$population %in% c("cAlaska", "cAlberta"),'population']<-'coyote'
out[out$population %in% c("dGreenland", "dAlaska"),'population']<-'arcticdog'
out[out$population %in% c("dGshep", "dChowchow"),'population']<-'breeddog'


write.table(out, file=paste("args/relate/input/aw.overallgrouping.poplabels",sep=""), 
    quote = FALSE, sep = " ",row.names = FALSE,col.names = TRUE)
#-----------------------------------------------------------------------------------------------------------------------
# make poplabels for plots separating out grell wolf + gr dogs,  by arcticwolf, mainlandwolf, arcticdog, breeddog, coyote

#> unique(out$population)
# [1] "wGreenland"      "wNorthBaffin"    "wSouthBaffin"    "wAlaska"        
# [5] "wEllesmere"      "wBanksIsland"    "wVictoriaIsland" "wAtlantic"      
# [9] "wCentral"        "cAlaska"         "wAlberta"        "cAlberta"       
#[13] "dGreenland"      "dAlaska"         "dGshep"          "dChowchow"      

#out[out$population %in% c("wGreenland", "wNorthBaffin","wSouthBaffin", "wEllesmere", "wBanksIsland",    "wVictoriaIsland"),'population']<-'arcticwolf'
out[out$population %in% c("wNorthBaffin","wSouthBaffin", "wBanksIsland",    "wVictoriaIsland"),'population']<-'arcticwolf'
out[out$population %in% c("wGreenland", "wEllesmere"),'population']<-'GR_ELL_wolf'
out[out$population %in% c("wAlaska", "wAtlantic","wCentral", "wAlberta"),'population']<-'mainlandwolf'
out[out$population %in% c("cAlaska", "cAlberta"),'population']<-'coyote'
#out[out$population %in% c("dGreenland", "dAlaska"),'population']<-'arcticdog'
out[out$population %in% c("dGshep", "dChowchow"),'population']<-'breeddog'


write.table(out, file=paste("args/relate/input/aw.diffgrouping.poplabels",sep=""), 
    quote = FALSE, sep = " ",row.names = FALSE,col.names = TRUE)
#-----------------------------------------------------------------------------------------------------------------------


