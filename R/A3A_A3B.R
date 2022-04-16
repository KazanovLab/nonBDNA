library(data.table)
library(ggplot2)

A3A_A3B_Gordenin <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/A3A_A3B/Gordenin_data/MAF_Aug31_2016_sorted_A3A_A3B_comparePlus.txt", sep='\t')
A3A_A3B_Gordenin <- data.table(A3A_A3B_Gordenin)

activity <- read.csv("/Users/mar/BIO/PROJECTS/PCAWG_APOBEC/PCAWG_enrichment_6cancers.txt",sep='\t',header=FALSE,strip.white =TRUE)
activity <- data.table(activity)
setnames(activity,c("project","sample","enrichment"))

data <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/A3A_A3B/AB_AllMu.txt", sep='\t')
data <- data.table(data)

data[, densityA := ApoA / TrgA]
data[, densityB := ApoB / TrgB]

dt <- data[, .(mutID,sample,densityA,densityB)]
dt <- merge(dt,activity,by="sample",all.x = TRUE)

dt <- merge(dt,A3A_A3B_Gordenin[,.(Tumor_Sample_Barcode,A3A_or_A3B)], by.x = "sample", by.y = "Tumor_Sample_Barcode")
dt[is.na(A3A_or_A3B),A3A_or_A3B := "ND"]
dt[, AorB_enrichment := paste0(A3A_or_A3B,"_",enrichment)]

for(cancer in c("BLCA","BRCA","HNSC","CESC","LUAD","LUSC")){

dt2 <- dt[mutID == cancer, .(AorB_enrichment,enrichment,densityA,densityB)]  
dt2 <- dt2[order(enrichment)] 
lvls <- dt2$AorB_enrichment

dt2melt <- melt(dt2, id.vars = c("AorB_enrichment","enrichment"))
dt2melt$AorB_enrichment <- factor(dt2melt$AorB_enrichment, levels=lvls)

p <- ggplot(dt2melt,aes(x=AorB_enrichment,y=value,fill=variable)) + geom_bar(stat="identity",position="dodge") +
  theme(panel.background = element_blank(),
        axis.text.x = element_text(size=8, angle = 90),
        axis.title = element_blank(),
        axis.line = element_line(color="black")#,
        #legend.position = "none"
  ) +
  scale_fill_manual(values=c(rgb(203,73,123,maxColorValue = 255), 
                             rgb(161,201,171,maxColorValue = 255), 
                             rgb(242,234,102,maxColorValue = 255),
                             rgb(139,197,229,maxColorValue = 255)))

if(cancer %in% c("BLCA","CESC"))
  wd <- 150
else if(cancer == "BRCA")
  wd <- 450
else if(cancer == "HNSC")
  wd <- 300
else if(cancer == "LUAD")
  wd <- 250
else if(cancer == "LUSC")
  wd <- 300
else
  wd <- 150

ggsave(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/A3A_A3B/plots/A3A_A3B_",cancer,".tiff"),plot=p,width=wd,height=100,units="mm",dpi=300)

}
