PATH_PROJECT="E:/src/r/SIRT_PTEN"
PATH="E:/data/tcga/"
DIRS=dir(path=PATH,pattern="*tcga$",full.names = TRUE,recursive = FALSE,ignore.case=TRUE)
FILE_EXPRESSION=data.frame()
data_all_tcga_mutation=data.frame(matrix(dimnames = list(c(),c("Hugo_Symbol")),ncol = 1))
data_all_tcga_mutation$Hugo_Symbol=as.character(data_all_tcga_mutation$Hugo_Symbol)
for(onedir in  DIRS){
  one=list.files(path = onedir,pattern = "data.*?RNA_Seq.*?v2.*?expression.*?median.*?txt",recursive = TRUE,ignore.case = TRUE,full.names = TRUE)
  #print(onedir)
  parent_file=stringr::str_extract(one,pattern = "^.*/")
  FILE_EXPRESSION[stringr::str_split(string = onedir,pattern = "/",simplify = TRUE)[4],"PATH"]=one
  FILE_EXPRESSION[stringr::str_split(string = onedir,pattern = "/",simplify = TRUE)[4],"PARENT_FILE"]=parent_file
  print(parent_file)
  setwd(onedir)
if(FALSE){
  ######寻找突变信息文件,合并突变信息文件
  file_mutation=list.files(path = getwd(),pattern = "^data_mutations*",all.files = TRUE,full.names = TRUE,ignore.case = TRUE,recursive = TRUE)
  print(file_mutation)
  #合并突变的信息||还是算了吧，合并的文件太大了
  for(onefile in  file_mutation){
    theonefile=read.csv(file = onefile,header = TRUE,stringsAsFactors = FALSE,sep = "\t")
    message(ncol(theonefile))
    print(colnames(theonefile))
    library(dplyr)
    ##进行变量的类型转变  
    for(thecolname in colnames(data_all_tcga_mutation)){
      data_all_tcga_mutation[,thecolname]=as.character(data_all_tcga_mutation[,thecolname])
    }
    data_mutation_selected=theonefile[which(theonefile$Entrez_Gene_Id %in% "7157"),]
    data_mutation_selected=data_mutation_selected[which(data_mutation_selected$Variant_Classification %in% c("Missense_Mutation")),]
    ##进行变量的类型转变  
    for(thecolname in colnames(data_mutation_selected)){
      data_mutation_selected[,thecolname]=as.character(data_mutation_selected[,thecolname])
    }
    data_all_tcga_mutation=full_join(data_all_tcga_mutation,data_mutation_selected)
    setwd(PATH_PROJECT)
  }
}
  }
save(file = "./data_all_tcga_mutation.save",list = c("data_all_tcga_mutation"))
 
load(file = "./data_all_tcga_mutation.save")


#########合并临床信息文件
data_clinical_all=data.frame(matrix(ncol = 1,dimnames = list(c(),c("PATIENT_ID"))))
for(onedir in  DIRS){
    setwd(onedir)
    files_clinical=list.files(path=getwd(),pattern = "^data_bcr_clinical_data_patient",all.files = TRUE,full.names = TRUE,ignore.case = TRUE,recursive = TRUE)
    data_clinical=read.csv(file=files_clinical[1],header = TRUE,skip = 4,sep="\t",stringsAsFactors = FALSE)
    for(colindex in c(1:ncol(data_clinical))){
      data_clinical[,colindex]=as.character(data_clinical[,colindex])
    }
    ######
    library(dplyr)
    for(colindex in c(1:ncol(data_clinical_all))){
      data_clinical_all[,colindex]=as.character(data_clinical_all[,colindex])
    }
    data_clinical_all=full_join(data_clinical_all,data_clinical)    
    
}
##保存合并的临床信息
setwd("E:/src/r/SIRT_PTEN");save(file = "data_clinical_all.save",data_clinical_all)



library(survival)
library(KMsurv)
data_clinical_all$OS_MONTHS=as.numeric(data_clinical_all$OS_MONTHS)
attach(data_clinical_all)
## 估计KM生存曲线
my.surv <- Surv(OS_MONTHS,OS_STATUS=='DECEASED')
kmfit1 <- survfit(my.surv~1)
summary(kmfit1)
win.graph();plot(kmfit1)


















