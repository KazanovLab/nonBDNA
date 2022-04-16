library(data.table)
library(ggplot2)

activity <- read.csv("/Users/mar/BIO/PROJECTS/PCAWG_APOBEC/PCAWG_enrichment_6cancers.txt",sep='\t',header=FALSE,strip.white =TRUE)
activity <- data.table(activity)
setnames(activity,c("project","sample","enrichment"))

for(cancer in c("BLCA","BRCA","HNSC","CESC","LUAD","LUSC")){
  
data <- read.csv(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/ReplicationStrand/GQ_rts/GQ_",cancer,"_rts.txt"),sep='\t')
data <- data.table(data)
setnames(data,"TrgStk..","TrgStkUnk")
setnames(data,"TrgLin..","TrgLinUnk")
setnames(data,"ApoStk..","ApoStkUnk")
setnames(data,"ApoLin..","ApoLinUnk")

data[, TrgStk := TrgStkFw + TrgStkBk + TrgStkUnk]
data[, TrgLin := TrgLinFw + TrgLinBk + TrgLinUnk]
data[, ApoStk := ApoStkFw + ApoStkBk + ApoStkUnk]
data[, ApoLin := ApoLinFw + ApoLinBk + ApoLinUnk]
data[, StkDensity := ApoStk / TrgStk]
data[, LinDensity := ApoLin / TrgLin]

data[, TrgFw := TrgStkFw + TrgLinFw]
data[, ApoFw := ApoStkFw + ApoLinFw]
data[, TrgBk := TrgStkBk + TrgLinBk]
data[, ApoBk := ApoStkBk + ApoLinBk]
data[, FwDensity := ApoFw / TrgFw]
data[, BkDensity := ApoBk / TrgBk]

data[, StkFwDensity := ApoStkFw / TrgStkFw]
data[, StkBkDensity := ApoStkBk / TrgStkBk]
data[, LinFwDensity := ApoLinFw / TrgLinFw]
data[, LinBkDensity := ApoLinBk / TrgLinBk]

data <- merge(data,activity,by="sample",all.x = TRUE)
data <- data[!is.na(enrichment)]

dt <- data[,.(enrichment,StkDensity,LinDensity)]
dtmelt <- melt(dt, id.vars = "enrichment")
dtmelt$enrichment <- as.factor(dtmelt$enrichment)

dt1 <- data[,.(enrichment,FwDensity,BkDensity,ApoFw,ApoBk)]
dt1melt <- melt(dt1, id = "enrichment", measure=patterns("Density","Apo"))
dt1melt[variable == 1, variable := "Fw"]
dt1melt[variable == 2, variable := "Bk"]
dt1melt$enrichment <- as.factor(dt1melt$enrichment)

dt2 <- data[,.(enrichment,StkFwDensity,ApoStkFw,StkBkDensity,ApoStkBk,LinFwDensity,ApoLinFw,LinBkDensity,ApoLinBk)]
dt2melt <- melt(dt2, id = "enrichment", measure=patterns("Density","Apo"))
dt2melt[variable == 1, variable := "StkFw"]
dt2melt[variable == 2, variable := "StkBk"]
dt2melt[variable == 3, variable := "LinFw"]
dt2melt[variable == 4, variable := "LinBk"]
dt2melt$enrichment <- as.factor(as.character(round(dt2melt$enrichment,2)))

p <- ggplot(dtmelt,aes(x=enrichment,y=value,fill=variable)) + geom_bar(stat="identity",position="dodge") +
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

ggsave(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/ReplicationStrand/plots/g4_",cancer,".tiff"),plot=p,width=wd,height=100,units="mm",dpi=300)

p <- ggplot(dt1melt,aes(x=enrichment,y=value1,fill=variable)) + geom_bar(stat="identity",position="dodge") +
  theme(panel.background = element_blank(),
        axis.text.x = element_text(size=8, angle = 90),
        axis.title = element_blank(),
        axis.line = element_line(color="black")#,
        #legend.position = "none"
  ) +
  scale_fill_manual(values=c(rgb(203,73,123,maxColorValue = 255), 
                             rgb(161,201,171,maxColorValue = 255), 
                             rgb(242,234,102,maxColorValue = 255),
                             rgb(139,197,229,maxColorValue = 255))) +
  geom_text(aes(enrichment, label = value2),position = position_dodge(width = 1),size=2)


ggsave(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/ReplicationStrand/plots/repstrand_",cancer,".tiff"),plot=p,width=wd,height=100,units="mm",dpi=300)


p <- ggplot(dt2melt,aes(x=enrichment,y=value1,fill=variable)) + geom_bar(stat="identity",position="dodge") +
  theme(panel.background = element_blank(),
        axis.text.x = element_text(size=8), #, angle = 90),
        axis.title = element_blank(),
        axis.line = element_line(color="black"),
        legend.position = "none"
  ) +
  scale_fill_manual(values=c(rgb(203,73,123,maxColorValue = 255), 
                             rgb(161,201,171,maxColorValue = 255), 
                             rgb(242,234,102,maxColorValue = 255),
                             rgb(139,197,229,maxColorValue = 255))) #+
  #geom_text(aes(enrichment, label = value2),position = position_dodge(width = 1),size=2)



ggsave(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/ReplicationStrand/plots/g4repstrand_",cancer,".tiff"),plot=p,width=wd,height=100,units="mm",dpi=300)

}

# Calculation out of tructures

for(cancer in c("BLCA","BRCA","HNSC","CESC","LUAD","LUSC")){
  data <- read.csv(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/ReplicationStrand/GQ_rts_out/GQ_",cancer,"out.txt"),sep='\t')
  data <- data.table(data)
 
  setnames(data,"Trg..","TrgUnk")
  setnames(data,"Apo..","ApoUnk")
 
  data[, FwDensity := ApoFw / TrgFw]
  data[, BkDensity := ApoBk / TrgBk]
  
  data <- merge(data,activity,by="sample",all.x = TRUE)
  
  
  dt3 <- data[,.(enrichment,FwDensity,BkDensity)]
  dt3melt <- melt(dt3, id = "enrichment", measure=patterns("Density"))
  dt3melt$enrichment <- as.factor(dt3melt$enrichment)

  p <- ggplot(dt3melt,aes(x=enrichment,y=value,fill=variable)) + geom_bar(stat="identity",position="dodge") +
    theme(panel.background = element_blank(),
          axis.text.x = element_text(size=8, angle = 90),
          axis.title = element_blank(),
          axis.line = element_line(color="black")#,
          #legend.position = "none"
    ) +
    scale_fill_manual(values=c(rgb(203,73,123,maxColorValue = 255), 
                               rgb(161,201,171,maxColorValue = 255), 
                               rgb(242,234,102,maxColorValue = 255),
                               rgb(139,197,229,maxColorValue = 255))) #+
    #geom_text(aes(enrichment, label = value2),position = position_dodge(width = 1),size=2)
  
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

  ggsave(paste0("/Users/mar/BIO/PROJECTS/APOBEC/NONBDNA/ReplicationStrand/plots/repstrandOut_",cancer,".tiff"),plot=p,width=wd,height=100,units="mm",dpi=300)
  
}
