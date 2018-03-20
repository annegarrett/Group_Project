#Load packages
library(ggpubr)
library(ggplot2)

#Load data
gene2go <- read.delim("gene2go", header=FALSE, comment.char="#")
colnames(gene2go)[3]<-"GO.Name"

ComponentOM.bvtUpAndDownRegulated <- read.csv("ComponentOM-bvtUpAndDownRegulated.tsv", sep="")
ComponentOM.bvtUpAndDownRegulated$GO_term <- "Component"
ComponentOM.bvtUpAndDownRegulated$Comparison<-"Old_trained"
ComponentYM.bv2UpAndDownRegulated <- read.csv("ComponentYM-bv2UpAndDownRegulated.tsv", sep="")
ComponentYM.bv2UpAndDownRegulated$GO_term<-"Component"
ComponentYM.bv2UpAndDownRegulated$Comparison<-"Young_2h"
ComponentYM.bvtUpAndDownRegulated <- read.csv("ComponentYM-bvtUpAndDownRegulated.tsv", sep="")
ComponentYM.bvtUpAndDownRegulated$GO_term<-"Component"
ComponentYM.bvtUpAndDownRegulated$Comparison<-"Young_trained"
FunctionOM.bvtUpAndDownRegulated <- read.csv("FunctionOM-bvtUpAndDownRegulated.tsv", sep="")
FunctionOM.bvtUpAndDownRegulated$GO_term<-"Function"
FunctionOM.bvtUpAndDownRegulated$Comparison<-"Old_trained"
FunctionYM.bvtUpAndDownRegulated <- read.csv("FunctionYM-bvtUpAndDownRegulated.tsv", sep="")
FunctionYM.bvtUpAndDownRegulated$GO_term<-"Function" 
FunctionYM.bvtUpAndDownRegulated$Comparison<-"Young_trained"
FunctionYM.bv2UpAndDownRegulated <- read.csv("FunctionYM-bv2UpAndDownRegulated.tsv", sep="")
FunctionYM.bv2UpAndDownRegulated$GO_term<-"Function"
FunctionYM.bv2UpAndDownRegulated$Comparison<-"Young_2h" 
ProcessOM.bv2UpAndDownRegulated <- read.csv("~/University/MSc Bioinformatics/Group Project/GO/ProcessOM-bv2UpAndDownRegulated.tsv", sep="")
ProcessOM.bv2UpAndDownRegulated$GO_term<-"Process"
ProcessOM.bv2UpAndDownRegulated$Comparison<-"Old 2h" 
ProcessOM.bvtUpAndDownRegulated <- read.csv("ProcessOM-bvtUpAndDownRegulated.tsv", sep="")
ProcessOM.bvtUpAndDownRegulated$GO_term<-"Process"
ProcessOM.bvtUpAndDownRegulated$Comparison<-"Old_trained"
ProcessYM.bv2UpAndDownRegulated <- read.csv("ProcessYM-bv2UpAndDownRegulated.tsv", sep="")
ProcessYM.bv2UpAndDownRegulated$GO_term<-"Process"
ProcessYM.bv2UpAndDownRegulated$Comparison<-"Young_2h"
ProcessYM.bvtUpAndDownRegulated <- read.csv("ProcessYM-bvtUpAndDownRegulated.tsv", sep="")
ProcessYM.bvtUpAndDownRegulated$GO_term<-"Process"
ProcessYM.bvtUpAndDownRegulated$Comparison<-"Young_trained"

#Join data together
all_go<-rbind(ComponentOM.bvtUpAndDownRegulated,ComponentYM.bvtUpAndDownRegulated,FunctionOM.bvtUpAndDownRegulated,FunctionYM.bvtUpAndDownRegulated,FunctionYM.bv2UpAndDownRegulated,ProcessOM.bv2UpAndDownRegulated,ProcessOM.bvtUpAndDownRegulated,ProcessYM.bv2UpAndDownRegulated,ProcessYM.bvtUpAndDownRegulated)
all_go<-left_join(all_go,gene2go[,c(3,6)])
all_go<-unique(all_go)

#Add -log 10 p-value
all_go$log10P<--log10(all_go$p.adjusted)

#Extract everything except young-trained and plot
all_go2<-all_go[-which(all_go$Comparison == "Young_trained",),]
g1<-ggplot(all_go2,aes(y=log10P, x=V6,fill=Comparison))+geom_col()+xlab("Go Term")+ylab("-log10(P-value)")+theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.text=element_text(size=8))+coord_flip()

#Extract young-trained
young_trained<-all_go[which(all_go$Comparison == "Young_trained",),]
#Filter lower p-values
young_trained<-young_trained[-which(young_trained$log10P<5),] 
#Plot young-trained
g2<-ggplot(young_trained,aes(y=log10P, x=V6,fill=Comparison))+geom_col()+xlab("Go Term")+ylab("-log10(P-value)")+scale_fill_manual(values=c("darkorchid2"))+theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.text=element_text(size=8))+coord_flip()

#Join plots
ggarrange(g1,g2,nrow = 2,labels = c("A","B"),heights=c(0.8,1))