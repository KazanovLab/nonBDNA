library(data.table)
library(ggpubr)
library(gridExtra)

activity <- read.csv("/Users/mar/BIO/PROJECTS/PCAWG_APOBEC/PCAWG_enrichment_6cancers.txt",sep='\t',header=FALSE,strip.white =TRUE)
activity <- data.table(activity)
setnames(activity,c("project","sample","enrichment"))

motifs<- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/R/all_motifs_Denek.txt",header = F)
motifs <- data.table(motifs)
setnames(motifs,"V1","motif")

dataTotal <- data.table()

for(m in motifs$motif){

  i <- 1
  plots <- list()
  
  data <- read.csv(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/ALL_TRIPLETS/merg_",m,"_X0.txt"), sep='\t')
  data <- data.table(data)
  data[, inMut := in_A + in_G + in_C + in_T]
  data[, outMut := out_A + out_G + out_C + out_T]
  data <- merge(data,activity,by="sample",all.x = TRUE)

  dataTotal <- rbind(dataTotal, data)
  
  data[, logDensity := log10((inMut/in_TRG)/(outMut/out_TRG))]

    
  for(cancer in c("BLCA","BRCA","CESC","HNSC","LUAD","LUSC")){
    for(s in c("dr","mr","str","gq","ir","apr")){
      
      dt <- data[struID == s & cancID == cancer]
      
      p <- ggplot(dt,aes(x=enrichment,y=logDensity)) + geom_point(color=rgb(38,120,178,maxColorValue = 255),size=0.1) +
        scale_shape_manual(values=c(15, 16)) +
        geom_hline(yintercept=0,color=rgb(243,94,90,maxColorValue = 255),size=0.3) +
        stat_smooth(method = "lm", col = rgb(140,210,185,maxColorValue = 255),se=F,formula = y ~ log(x),size=0.5) +
        #scale_color_manual(values=c(rgb(38,120,178,maxColorValue = 255),rgb(140,210,185,maxColorValue = 255))) +
        #ggtitle(paste0("Structure=",s," , Cancer=",c," , isAPOBEC=",a)) +
        theme(panel.background = element_blank(),
              plot.title = element_text(size=8),
              axis.title = element_blank(),
              axis.line = element_line(color="black",size=0.3),
              axis.text = element_text(size=4),
              axis.ticks = element_line(size=0.3),
              legend.position = "none",
              panel.grid.major = element_line(size = rel(0.5), colour='grey92'))
      
      plots[[i]] <- p
      i <- i + 1
      
    }
  }
  
  pl <- marrangeGrob(grobs=plots,nrow=6,ncol=6)
  ggexport(pl,filename=paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/ALL_TRIPLETS/plots/",m,".pdf"))
}

data2 <- dataTotal[,.("in_TRG"=sum(in_TRG),"in_C"=sum(in_C),"in_G"=sum(in_G),"in_A"=sum(in_A),"in_T"=sum(in_T),
                      "out_TRG"=sum(out_TRG),"out_C"=sum(out_C),"out_G"=sum(out_G),"out_A"=sum(out_A),"out_T"=sum(out_T),
                      "inMut"=sum(inMut),"outMut"=sum(outMut)), by=.(struID,cancID,sample,enrichment,"motif"=substr(signat,1,2))]

data2[, logDensity := log10((inMut/in_TRG)/(outMut/out_TRG))]

for(m in unique(data2$motif)){
  
  i <- 1
  plots <- list()
  
  data2motif <- data2[motif == m]
  
  for(cancer in c("BLCA","BRCA","CESC","HNSC","LUAD","LUSC")){
    for(s in c("dr","mr","str","gq","ir","apr")){
      
      dt2 <- data2motif[struID == s & cancID == cancer]
      
      p <- ggplot(dt2,aes(x=enrichment,y=logDensity)) + geom_point(color=rgb(38,120,178,maxColorValue = 255),size=0.1) +
        scale_shape_manual(values=c(15, 16)) +
        geom_hline(yintercept=0,color=rgb(243,94,90,maxColorValue = 255),size=0.3) +
        stat_smooth(method = "lm", col = rgb(140,210,185,maxColorValue = 255),se=F,formula = y ~ log(x),size=0.5) +
        #scale_color_manual(values=c(rgb(38,120,178,maxColorValue = 255),rgb(140,210,185,maxColorValue = 255))) +
        #ggtitle(paste0("Structure=",s," , Cancer=",c," , isAPOBEC=",a)) +
        theme(panel.background = element_blank(),
              plot.title = element_text(size=8),
              axis.title = element_blank(),
              axis.line = element_line(color="black",size=0.3),
              axis.text = element_text(size=4),
              axis.ticks = element_line(size=0.3),
              legend.position = "none",
              panel.grid.major = element_line(size = rel(0.5), colour='grey92'))
      
      plots[[i]] <- p
      i <- i + 1
      
    }
  }
  
  pl <- marrangeGrob(grobs=plots,nrow=6,ncol=6)
  ggexport(pl,filename=paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/ALL_TRIPLETS/plots/",m,".pdf"))
}

data3motif <- data2[motif %in% c("AC","GC","CC")]

data3 <- data3motif[,.("in_TRG"=sum(in_TRG),"in_C"=sum(in_C),"in_G"=sum(in_G),"in_A"=sum(in_A),"in_T"=sum(in_T),
                      "out_TRG"=sum(out_TRG),"out_C"=sum(out_C),"out_G"=sum(out_G),"out_A"=sum(out_A),"out_T"=sum(out_T),
                      "inMut"=sum(inMut),"outMut"=sum(outMut)), by=.(struID,cancID,sample,enrichment)]

data3[, logDensity := log10((inMut/in_TRG)/(outMut/out_TRG))]


  i <- 1
  plots <- list()
  
  for(cancer in c("BLCA","BRCA","CESC","HNSC","LUAD","LUSC")){
    for(s in c("dr","mr","str","gq","ir","apr")){
      
      dt3 <- data3[struID == s & cancID == cancer]
      
      p <- ggplot(dt3,aes(x=enrichment,y=logDensity)) + geom_point(color=rgb(38,120,178,maxColorValue = 255),size=0.1) +
        scale_shape_manual(values=c(15, 16)) +
        geom_hline(yintercept=0,color=rgb(243,94,90,maxColorValue = 255),size=0.3) +
        stat_smooth(method = "lm", col = rgb(140,210,185,maxColorValue = 255),se=F,formula = y ~ log(x),size=0.5) +
        #scale_color_manual(values=c(rgb(38,120,178,maxColorValue = 255),rgb(140,210,185,maxColorValue = 255))) +
        #ggtitle(paste0("Structure=",s," , Cancer=",c," , isAPOBEC=",a)) +
        theme(panel.background = element_blank(),
              plot.title = element_text(size=8),
              axis.title = element_blank(),
              axis.line = element_line(color="black",size=0.3),
              axis.text = element_text(size=4),
              axis.ticks = element_line(size=0.3),
              legend.position = "none",
              panel.grid.major = element_line(size = rel(0.5), colour='grey92'))
      
      plots[[i]] <- p
      i <- i + 1
      
    }
  }
  
  pl <- marrangeGrob(grobs=plots,nrow=6,ncol=6)
  ggexport(pl,filename=paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/ALL_TRIPLETS/plots/VC.pdf"))


  ## Comparison effects in TC and VC
  dataTC <- data2[motif == "TC"]
  dataVC <- data3
    
  i <- 1
  plots <- list()
  
  for(cancer in c("BLCA","BRCA","CESC","HNSC","LUAD","LUSC")){
    for(s in c("dr","mr","str","gq","ir","apr")){
      
      dtTC <- dataTC[struID == s & cancID == cancer]
      dtVC <- data3[struID == s & cancID == cancer]
      dtVC[, motif := "VC"]
      
      dtTCVC <- rbind(dtVC,dtTC)
      dtTCVC$motif <- as.factor(dtTCVC$motif)
      
      dtTCVC <- dtTCVC[is.finite(logDensity)]
      dtTCVC[,motifBin := ifelse(motif == "TC",1,0)]
      
      reg <- lm(logDensity ~ enrichment + motifBin,dtTCVC)
      coef <- summary(reg)$coefficients[3,4]
      coefstr <- format(coef,scientific=TRUE,digits=3)
      
      maxy <- max(dtTCVC$logDensity, na.rm = TRUE)
      miny <- min(dtTCVC$logDensity, na.rm = TRUE)
      py <- miny + (maxy - miny)*0.95
      
      
      p <- ggscatter(dtTCVC, x = "enrichment", y = "logDensity", fill="motif",color="motif",
                     add = "reg.line", alpha = 0.75,size=0.1, add.params = list(size=0.2),
                     xlab = FALSE, ylab = FALSE, xticks.by = 0.5,
                     palette = c(rgb(44,0,136,maxColorValue = 255),rgb(249,204,29,maxColorValue = 255))) + annotate("text", x=2, y=py, label= paste0("p = ",coefstr),size=1.5) + 
        scale_size(range = c(0.5, 10)) +
        theme(panel.background = element_blank(),
              plot.title = element_text(size=8),
              axis.title = element_blank(),
              axis.line = element_line(color="black",size=0.3),
              axis.text = element_text(size=4),
              axis.ticks = element_line(size=0.3),
              legend.position = "none",
              panel.grid.major = element_line(size = rel(0.5), colour='grey92')) 
      
      plots[[i]] <- p
      i <- i + 1
      
    }
  }
  
  pl <- marrangeGrob(grobs=plots,nrow=6,ncol=6)
  ggexport(pl,filename=paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/ALL_TRIPLETS/plots/TCVC.pdf"))
  
      
      
