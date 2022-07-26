proteome<-read.csv("77_cancer_proteomes_CPTAC_itraq.csv")
clinicaldata<-read.csv("clinical_data_breast_cancer.csv")
proteins<-read.csv("PAM50_proteins.csv")
################################################question1
c<-proteome[grepl("TCGA",colnames(proteome))]
rown<-c(clinicaldata[,1])
namesofproteome<-c(colnames(c))
rep<-sub('.A', "A", namesofproteome)
r<-sub("\\.","", rep)
r<-substr(r,0,6)
w<-substr(rown,6,12)
rowrep<-sub('-', "",w)
colnames(c)=r
editcolom=colnames(c)
editt=list()
for(i in r)
{
  editt<-append(editt, as.character( rowrep[grepl(i,rowrep)]))
}
editt<-as.character(editt)
install.packages("dplyr") 
library("dplyr")
newdata=data.frame()
newdata<-tibble(c[,editt])
##############################################################question2
library(tibble)
library(forcats)
library(dplyr)
library(tidyverse)
newdata$refseq=proteome$RefSeq_accession_number
newdata<-na.omit(newdata)
names(proteins)[2] <-"refseq"
newtable=newdata %>% inner_join(proteins,by="refseq")
transv<-t(newtable)
transv<-as.data.frame(transv)
colnames(transv)=newtable$refseq
clinicaldata$Complete.TCGA.ID=rowrep
transv=transv[1:80,]
transv$id=row.names(transv)
row.names(transv)=c(1:80)
names(clinicaldata)[1] <-"id"
df=transv %>% inner_join(clinicaldata,by="id")
newstatus<-sub('Positive', 1,df$HER2.Final.Status)
newstatus<-sub('Negative', -1,newstatus)
newstatus<-sub('Equivocal', 0,newstatus)
newstatus=as.numeric(newstatus)
newstatus=as.data.frame(newstatus)
df=df[,1:26]
result=list()
protein=list()
threshold=0.03
for (i in c(1:26))
{
  result<-append(result,cor(as.numeric(df[,i]), newstatus, method ="pearson"))
  proteinselect<-cor(as.numeric(df[,i]), newstatus, method ="pearson")
  if(proteinselect>threshold)
  {
    protein<-append(protein,colnames(df[i]))
  }
}
protein<-as.character(protein)
res<-as.numeric(result)
result=sort(res,decreasing = TRUE)
filterresult=list()
for(v in result)
{
  if (v>threshold)
  {
    filterresult=append(filterresult,v)
  }
  
}
filterresult=as.character(filterresult)
####################################################################question3
library("dplyr")
datafilter<-filter(newstatus,newstatus!=0)
d<-sub("-1",'Negative',datafilter$newstatus)
d<-sub("1",'postive',d)
d<-as.data.frame(d)
r<-df[1:4,]
s<-df[6:77,]
df<-rbind(r,s)
resultt=list()
proteint=list()
signficant=list()
for (i in c(1:26))
{
  resutofttest<-t.test(as.numeric(df[,i])~d$d)
  filterproteint<-resutofttest$p.value
  if (filterproteint<0.05)
  {  ####### true signficant relation
    signficant<-append(signficant,filterproteint)
    ###########p>0.05 means nosignificant no relation
    if(filterproteint>threshold)
    {
      proteint<-append(proteint,colnames(df[i]))
    }
  }
}

resultt=as.numeric(signficant)
resultoftest=sort(resultt,decreasing = TRUE)
filterresultoftest=list()
for(m in resultoftest)
{
  if (m>threshold)
  {
    filterresultoftest=append(filterresultoftest,m)
    
  }
}
filterresultoftest=as.character(filterresultoftest)
proteint<-as.character(proteint)
###############on same threshold we find 2 result of selected feature set are different #3.4