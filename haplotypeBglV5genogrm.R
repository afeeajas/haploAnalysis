
makehaploblocks <- function(mapinfohap='/home/afeesa/paper3/data/haploview_dat/GENOimp_filtssa12.map',nbsize=10000){

  markermap <- read.table(paste0(mapinfohap),stringsAsFactors = F)[,-3]
  colnames(markermap) <- c('CHR','SNP','BP')
  markermap$nnum <- 1:nrow(markermap)
  gnome.chr.str <- head(markermap$BP,n=1)
  gnome.chr.end <- ceiling(tail(markermap$BP,n=1))
  slidewind_str <- seq(gnome.chr.str,gnome.chr.end,nbsize)
  slidewind_end <- c(seq(gnome.chr.str+nbsize+1,gnome.chr.end,nbsize+1),gnome.chr.end)
  haploblocks <- data.frame(blocknr=numeric(),hapstr=numeric(),
                            hapend=numeric(),stringsAsFactors = F)
  for(w in 1:length(slidewind_end)){
    nwind_snp <- markermap[which(markermap$BP>=slidewind_str[w] & 
                                   markermap$BP<=slidewind_end[w]),]
    if(nrow(nwind_snp)<1){next()}
    nwind_snp <- cbind.data.frame(blocknr=w,hapstr=nwind_snp$nnum[1],
                                  hapend=nwind_snp$nnum[nrow(nwind_snp)],
                                  stringsAsFactors=F)
    haploblocks <- rbind.data.frame(haploblocks,nwind_snp,stringsAsFactors = F)
  }
  haploblocks$nsnpinwind <- haploblocks$hapend-haploblocks$hapstr
  haploblocks <- haploblocks[which(haploblocks$nsnpinwind>1),-4]
return(haploblocks)
}
 
getblockhaploview <- function(haploviewfile){
  hapblocks <- read.table(haploviewfile,fill=T,stringsAsFactors=F,na.strings=" ")
  hapblocks <- hapblocks[which(hapblocks$V1=='BLOCK'),]
  hapblocks <- cbind.data.frame(NR=1:nrow(hapblocks),hapblocks)
  extractHAPs <- data.frame(blocknr=character(),hapstr=numeric(),hapend=numeric())
  for(i in 1:nrow(hapblocks)){
    gethapblock <- hapblocks[i,5:ncol(hapblocks)]
    gethapblock <- na.omit(t(gethapblock))
    gethapblock <- as.vector(gethapblock[which(gethapblock[,1]!=""),])
    hapblockstr <- as.numeric(as.vector(gsub("!","",gethapblock[1])))
    hapblockend <- as.numeric(as.vector(gsub("!","",gethapblock[length(gethapblock)])))
    hapsinfo <- cbind.data.frame(blocknr=hapblocks[i,3],hapstr=hapblockstr,hapend=hapblockend)
    extractHAPs <- rbind.data.frame(extractHAPs,hapsinfo)
  }
  return(extractHAPs)
}

makehaplotypes<- function(phasedbgl,mapinfohap,hapblocks,geno=T,outname='trait',hapfreqThresh=0.05){
  #phasedbgl='/home/afeesa/paper3/data/haploview_dat/GENOimphased_ssa12.vcf.gz.recode.vcf'
  hapgeno.bgl <- read.table(paste(phasedbgl),header=T,stringsAsFactors = F,comment.char="*",skip=8,check.names=F)
  #mapinfohap='/home/afeesa/paper3/data/haploview_dat/GENOimp_filtssa12.map'
  map <- read.table(paste(mapinfohap),header=F,stringsAsFactors = F)
  snpids <- hapgeno.bgl[,3]
  hapgeno.bgl <- t(hapgeno.bgl[,-c(1:9)])
  colnames(hapgeno.bgl) <- snpids
  nanim <- nrow(hapgeno.bgl)
  nhaps <- nrow(hapgeno.bgl) * 2
  animids <- rownames(hapgeno.bgl)
 
  removedose <- function(x){strsplit(x,":")}
  cat("... Removing dosage field ...\n")
  hapgeno.bgl <- matrix(unlist(apply(hapgeno.bgl,1,removedose)),ncol=2,byrow=T)
  hapgeno.bgl <- as.vector(hapgeno.bgl[,1])
  hapgeno.bgl <- matrix(hapgeno.bgl,ncol=length(snpids),byrow=T)
  colnames(hapgeno.bgl) <- snpids
  rownames(hapgeno.bgl) <- animids

  sepbypipe <- function(x){strsplit(x,"|")} 


for(j in 1:nrow(hapblocks)){
    cat("... Starting haplotype analysis for Block ",j," ...","\n")
    snps <- map[c(hapblocks[j,2]:hapblocks[j,3]),2]
    haplotype <- hapgeno.bgl[,snps]
    haplotype <-  matrix(unlist(apply(haplotype,1,sepbypipe)),ncol=3,byrow=T)
    
    #maternal haplotye
    mathaplo <- as.vector(haplotype[,1])
    mathaplo <- matrix(mathaplo,ncol=length(snps),byrow=T)
    mathaplo <- data.frame(mathaplo=apply(mathaplo,1,paste, collapse=""))

   #pathaplo
    pathaplo <- as.vector(haplotype[,3])
    pathaplo <- matrix(pathaplo,ncol=length(snps),byrow=T)
    pathaplo <- data.frame(pathaplo=apply(pathaplo,1,paste, collapse=""))

   #haps
    haps <- cbind.data.frame(mathaplo,pathaplo)
    numhaps <- unclass(c(as.factor(haps[,1]),as.factor(haps[,2])))

    haps$mathaps.num <- numhaps[1:nanim]
    haps$pathaps.num <- numhaps[(1+nanim):nhaps]
    haps$FID <- animids
    haps$IID <- animids
    haps <- haps[,c("FID","IID","mathaps.num","pathaps.num","mathaplo","pathaplo")]

    namehaplo <- paste("hapBlock-",j,sep="")
    phase.haplonum <- rbind.data.frame(data.frame(HAP=haps$mathaps.num),data.frame(HAP=haps$pathaps.num))
    phase.haplo <- rbind.data.frame(data.frame(HAP=haps$mathaplo),data.frame(HAP=haps$pathaplo))
    hapnumalleles <- unique(cbind.data.frame(HAPNUM=phase.haplonum,HAPALLELE=phase.haplo))
    hapnumalleles <- hapnumalleles[order(hapnumalleles$HAP),]

    HAPs.results <- data.frame(table(phase.haplonum))
    HAPs.results$propHAP <- round(HAPs.results$Freq/sum(HAPs.results$Freq),4)
    HAPs.results <- merge(HAPs.results,hapnumalleles,by=1)
    HAPs.results <- HAPs.results[order(-HAPs.results$Freq),]

    maphaplo <- map[map[,2] %in% snps,]

    haps <- list(haps); names(haps) <- namehaplo
    HAPs.results <- list(HAPs.results); names(HAPs.results) <- namehaplo
    maphaplo <- list(maphaplo); names(maphaplo) <- namehaplo
    
    if(j==1){
      haplos <- haps
      haplos.freq <- HAPs.results
      maphaplotype <- maphaplo
    } else {
      haplos<-c(haplos,haps)
      haplos.freq<-c(haplos.freq,HAPs.results)
      maphaplotype<-c(maphaplotype,maphaplo)
    }
  } 
    cat("... All haplotype extraction done ...\n")
    if(geno==F){
    return(list(HAP_ALLELES=haplos,HAP_FREQ=haplos.freq,MAP_info=maphaplotype))
 }

   if(geno==T){
      
     cat('... Starting conversion to Pseudo SNP ...\n')

     for(k in 1:nrow(hapblocks)){
       freqthaps <- haplos.freq[[k]][which(haplos.freq[[k]]$propHAP>=hapfreqThresh),]
       hapclasses <- as.vector(freqthaps$HAP)
       nhap <- as.numeric(length(hapclasses))
 
       haplomatrix <- matrix(0,nrow=nanim,ncol=nhap)
       names_HAP <- paste('Block',k,'HAP_',hapclasses,sep='')
       colnames(haplomatrix) <- names_HAP
        for(m in 1:nrow(haplos[[k]])){
          for(n in 1:length(hapclasses)){
          pat=ifelse(haplos[[k]][m,'pathaps.num']==hapclasses[n],1,0)
          mat=ifelse(haplos[[k]][m,'mathaps.num']==hapclasses[n],1,0)
          haplomatrix[m,n] <- pat+mat
          }
        }

      if(k==1){
        haplomat <- haplomatrix
        MAP_info= maphaplotype[[k]]
      } else {
        haplomat <- cbind.data.frame(haplomat,haplomatrix)
        MAP_info <- rbind.data.frame(MAP_info,maphaplotype[[k]])
     }
        rownames(haplomat) <- haplos[[k]]$IID
    }

     convert2ped <- function(dat){
     dat <- data.frame(dat)
     dat[which(dat[,1]==0),1] <- '11'; dat[which(dat[,1]==1),1] <- '12'; dat[which(dat[,1]==2),1] <- '22'
     return(dat)
     }

     haps2ped <- data.frame(apply(haplomat,2,convert2ped))
     colnames(haps2ped) <- colnames(haplomat)
     rownames(haps2ped) <- rownames(haplomat)
     haps2ped$FID<-rownames(haps2ped);haps2ped$IID<-rownames(haps2ped)
     haps2ped$PID<-0;haps2ped$MID<-0;haps2ped$Sex<-0;haps2ped$Phen<--9
     cat('Outputing data in plink map format\n')
     nopos <- ncol(haps2ped)
     write.table(haps2ped[,c((nopos-5):nopos,1:(nopos-6))],paste(outname,'.ped',sep=''),quote=F,col.names=F,row.names=F)
     mapR <- cbind.data.frame(CHR=unique(MAP_info[,1]),name=colnames(haplomat),GP=0,PP=1:ncol(haplomat))
     write.table(mapR,paste(outname,'.map',sep=''),quote=F,col.names=F,row.names=F)
 }
     return(list(HAP_ALLELES=haplos,HAP_FREQ=haplos.freq,MAP_info=maphaplotype,haplomatrix=haplomat))
}

