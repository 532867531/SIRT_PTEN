getGene="Sirt"
PATH_PROJECT="E:/src/r/SIRT_PTEN"
PATH="E:/data/tcga/"
DIRS=dir(path=PATH,pattern="*tcga$",full.names = TRUE,recursive = FALSE,ignore.case=TRUE)
FILE_EXPRESSION=data.frame()
for(onedir in  DIRS){
  one=list.files(path = onedir,pattern = "data.*?RNA_Seq.*?v2.*?expression.*?median.*?RDS",recursive = TRUE,ignore.case = TRUE,full.names = TRUE)
  print(onedir)
  parent_file=stringr::str_extract(one,pattern = "^.*/")
  FILE_EXPRESSION[stringr::str_split(string = onedir,pattern = "/",simplify = TRUE)[4],"PATH"]=one
  FILE_EXPRESSION[stringr::str_split(string = onedir,pattern = "/",simplify = TRUE)[4],"PARENT_FILE"]=parent_file
}
###加入注释情况
FILE_EXPRESSION[,"Description"]=c(
  "Adrenocortical carcinoma肾上腺皮质癌",
  "Bladder Urothelial Carcinoma膀胱尿路上皮癌",
  "Breast invasive carcinoma乳腺浸润性癌",
  "Cervical squamous cell carcinoma and endocervical adenocarcinoma宫颈鳞状细胞癌与宫颈腺癌",
  "Cholangio carcinoma胆管癌",
  "Colon adenocarcinoma结肠腺癌",
  "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma淋巴肿瘤弥漫大B细胞淋巴瘤",
  "Esophageal carcinoma食管癌",
  "Glioblastoma multiforme多形性胶质母细胞瘤",
  "Head and Neck squamous cell carcinoma头颈部鳞状细胞癌",
  "Kidney Chromophobe肾嫌色细胞",
  "Kidney renal clear cell carcinoma肾透明细胞癌",
  "Kidney renal papillary cell carcinoma肾乳头状细胞癌",
  "Acute Myeloid Leukemia急性髓系白血病",
  "Brain Lower Grade Glioma脑低级别胶质瘤",
  "Liver hepatocellular carcinoma肝细胞肝癌",
  "Lung adenocarcinoma肺腺癌",
  "Lung squamous cell carcinoma肺鳞状细胞癌",
  "Mesothelioma间皮瘤",
  "Ovarian serous cystadenocarcinoma卵巢浆液性囊腺癌",
  "Pancreatic adenocarcinoma胰腺癌",
  "Pheochromocytoma and Paraganglioma嗜铬细胞瘤与副神经节瘤",
  "Prostate adenocarcinoma前列腺腺癌",
  "Sarcoma肉瘤",
  "Skin Cutaneous Melanoma皮肤皮肤黑色素瘤",
  "Stomach adenocarcinoma胃腺癌",
  "Testicular Germ Cell Tumors睾丸生殖细胞肿瘤",
  "Thyroid carcinoma甲状腺癌",
  "Thymoma胸腺瘤",
  "Uterine Corpus Endometrial Carcinoma子宫体部子宫内膜癌",
  "Uterine Carcinosarcoma子宫肉瘤",
  "Uveal Melanoma葡萄膜黑色素瘤"
)












##开始看看表达高低情况
# apply(FILE_EXPRESSION,MARGIN = 1,FUN = function(x){
#   #fix(x)
#   x=t(x)
#   #fix(x)
#   file=x[1,"PATH"]
#   RDSfile=sub(pattern = "txt",replacement = "RDS",perl=TRUE,x=file)
#   if(file.exists(RDSfile)){
#     rna_seq_data=readRDS(file=RDSfile)
#   }else{
#     rna_data=read.csv(file = file,sep = "\t")
#     colnames(rna_data)=gsub(pattern = "\\.",replacement = "-",ignore.case = TRUE,perl = TRUE,x=colnames(rna_data))
#     saveRDS(rna_data,file=RDSfile)
#     rna_seq_data=readRDS(file=RDSfile)
#     message(paste("成功生成RDS",RDSfile))
#   }
#   setwd(PATH_PROJECT)
# })

##############读取RNAseq文件##########################
gogogo=function(cohortIndex){
  print(cohortIndex)
  thecohort=paste(rownames(FILE_EXPRESSION)[cohortIndex],FILE_EXPRESSION$Description[cohortIndex]);print(thecohort)
  setwd(FILE_EXPRESSION$PARENT_FILE[cohortIndex])
  file=FILE_EXPRESSION$PATH[cohortIndex]
  RDSfile=sub(pattern = "txt",replacement = "RDS",perl=TRUE,x=file)
  if(file.exists(RDSfile)){
    rna_seq_data=readRDS(file=RDSfile)
  }else{
    rna_data=read.csv(file = file,sep = "\t")
    colnames(rna_data)=gsub(pattern = "\\.",replacement = "-",ignore.case = TRUE,perl = TRUE,x=colnames(rna_data))
    saveRDS(rna_data,file=RDSfile)
    rna_seq_data=readRDS(file=RDSfile)
    #message(paste("成功生成RDS",RDSfile))
  }
  data_selected_gene=rna_seq_data[which(rna_seq_data$Entrez_Gene_Id=="23411"),]
  data2use=data_selected_gene[,-which(colnames(data_selected_gene) %in% c("Hugo_Symbol","Entrez_Gene_Id"))]
  data2use=t(data2use)
  library(ggpubr)
  q0=ggdensity(unlist(data_selected_gene),title = thecohort);plot(q0)
  quantile(unlist(data_selected_gene))
  high_percentile=0.70
  low_percentile=0.30
  value_low_selected_criteria=quantile(unlist(data_selected_gene),low_percentile);value_low_selected_criteria
  value_high_selected_criteria=quantile(unlist(data_selected_gene),high_percentile);value_high_selected_criteria
  
  sample_low_selected=data.frame(data2use[which(data2use<=value_low_selected_criteria),]);colnames(sample_low_selected)="value";sample_low_selected[,"group"]="low"
  sample_high_selected=data.frame(data2use[which(data2use>=value_high_selected_criteria),]);colnames(sample_high_selected)="value";sample_high_selected[,"group"]="high"
  
  names_sample_low_selected=rownames(sample_low_selected)
  names_sample_high_selected=rownames(sample_high_selected)
  ##画图表明有差异，后续需要和正常组织的对比gtex
  data_boxplot=rbind(sample_low_selected,sample_high_selected)
  library(ggplot2) #调用ggplot软件包
  p1<-ggplot(data_boxplot,aes(x=factor(group),y=value))+ggtitle(thecohort,subtitle = paste("<=",low_percentile,">=",high_percentile,sep=""))
  p1=p1+geom_boxplot(col="blue",pch=16,cex=1)+ylab(getGene)+geom_point(position="jitter",col=2,pch=16,cex=1)
  plot(p1)
  ###############################################################################################
  
  
  ##读取突变信息文件
  file_mutation=list.files(path = getwd(),pattern = "^data_mutations.*\\.txt1$",all.files = TRUE,full.names = TRUE,ignore.case = TRUE)
  print(file_mutation)
  data_mutation=data.frame()
  for(onefile in  file_mutation){
    ##测试文本有几行
    library(LaF);
    thenumberrows=LaF::determine_nlines(filename = onefile);thenumberrows
    library(readr);library(data.table);library(reader)
    theonefile=readr::read_delim(file = onefile,col_names = TRUE,delim =  "\t",trim_ws = FALSE,guess_max = 0,comment = "#",quote = "\"")
    
    #theonefile=read.delim(file = onefile,sep = "\t",header = TRUE)
    theonefile=as.data.frame(theonefile)
    if(NROW(theonefile)!=thenumberrows-1){message(paste(NROW(theonefile),"!=",thenumberrows-1));stop()}else{print(paste("行数相符，为",thenumberrows))}
    data_mutation=rbind(x=data_mutation,y=theonefile,all=T)
  }
  data_mutation=unique(data_mutation)
  
  ##查看突变的P53类型
  data_mutation_selected=data_mutation[which(data_mutation$Hugo_Symbol=="TP53"),]
  data_mutation_selected=data_mutation_selected[which(data_mutation_selected$Variant_Classification %in% c("Missense_Mutation")),]
  sample_data_mutation_selected=data_mutation_selected$Tumor_Sample_Barcode
  
  ####进行二条件组合样品:交叉高低表达（以整个cohort为参照）和突变
  names_sample_low_mutated_selected=unique(intersect(names_sample_low_selected,sample_data_mutation_selected));names_sample_low_mutated_selected
  names_sample_high_mutated_selected=unique(intersect(names_sample_high_selected,sample_data_mutation_selected));names_sample_high_mutated_selected
  
  final_patientid_low_mutated_selected=stringi::stri_extract(str=names_sample_low_mutated_selected,regex = "^.*(?=-\\d\\d$)");final_patientid_low_mutated_selected
  final_patientid_high_mutated_selected=stringi::stri_extract(str=names_sample_high_mutated_selected,regex = "^.*(?=-\\d\\d$)");final_patientid_high_mutated_selected
  ###############################################################################################
  
  
  ###############################################################################################
  ################从突变的样本中选取高表达和低表达的样本
  #data_expression_sample_data_mutation_selected=data2use[sample_data_mutation_selected,]
  ##查看突变群体的分布
  
  
  
  
  ##从突变的样本中选取高表达和低表达的样本END
  
  
  ###############################################################################################
  ##读取临床信息文件
  files_clinical=list.files(path=getwd(),pattern = "^data_bcr_clinical_data_patient",all.files = TRUE,full.names = TRUE,ignore.case = TRUE)
  data_clinical=read.csv(file=files_clinical[1],header = TRUE,skip = 4,sep="\t",stringsAsFactors = FALSE)
  ###############################################################################################
  ####先不划生存曲线，先划boxplot或者做出个四格表
  #data_clinical$OS_STATUS可以看同时期比例生存和死亡的比例而不加入时间因素
  #unique(data_clinical$OS_MONTHS)
  #data_clinical$DFS_STATUS
  #data_clinical$DFS_MONTHS
  OS_STARTUS_final_patientid_mutated_selected=list()
  OS_STARTUS_final_patientid_mutated_selected[["low"]]=data_clinical[which(data_clinical$PATIENT_ID %in%  final_patientid_low_mutated_selected),"OS_STATUS"]
  OS_STARTUS_final_patientid_mutated_selected[["high"]]=data_clinical[which(data_clinical$PATIENT_ID %in%  final_patientid_high_mutated_selected),"OS_STATUS"]
  library(DT)
  status=unique(unlist(OS_STARTUS_final_patientid_mutated_selected));print(status)
  status=sort(status,decreasing = FALSE);print(status);
  group=c("high","low")
  fourfoldTable= data.frame(matrix(nrow = length(group),ncol = length(status),dimnames = list(group,status)),stringsAsFactors = FALSE)
tryCatch({
    for(rowindex in c(1:length(OS_STARTUS_final_patientid_mutated_selected))){
        for(colindex in c(1:length(status))){
          fourfoldTable[rowindex,colindex]=length(which(OS_STARTUS_final_patientid_mutated_selected[[rownames(fourfoldTable)[rowindex]]]==colnames(fourfoldTable)[colindex]))
        }
    }
  fourfoldplot(x=as.matrix(fourfoldTable), color = c("#99CCFF", "#6699CC"),
               conf.level = 0.95,
               std = c("margins"),
               margin = c(1), space = 0.2, main = thecohort
               )
  
}        
,error=function(x){
message(x)
})
###试图绘制生存曲线
  OS_SURVIVAL_final_patientid_mutated_selected=list()
  OS_SURVIVAL_final_patientid_mutated_selected[["low"]]=data_clinical[which(data_clinical$PATIENT_ID %in%  final_patientid_low_mutated_selected),c("OS_MONTHS","OS_STATUS")]
  OS_SURVIVAL_final_patientid_mutated_selected[["high"]]=data_clinical[which(data_clinical$PATIENT_ID %in%  final_patientid_high_mutated_selected),c("OS_MONTHS","OS_STATUS")]
  library(tibble)
  OS_SURVIVAL_final_patientid_mutated_selected[["low"]]=add_column(OS_SURVIVAL_final_patientid_mutated_selected[["low"]],GROUP="LOW")
  OS_SURVIVAL_final_patientid_mutated_selected[["high"]]=add_column(OS_SURVIVAL_final_patientid_mutated_selected[["high"]],GROUP="HIGH")
  #OS_SURVIVAL_final_patientid_mutated_selected[["low"]][,"GROUP"]="low"
  #OS_SURVIVAL_final_patientid_mutated_selected[["high"]][,"GROUP"]="high"
  rbind_OS_SURVIVAL_final_patientid_mutated_selected<<-rbind(OS_SURVIVAL_final_patientid_mutated_selected[["low"]],OS_SURVIVAL_final_patientid_mutated_selected[["high"]])
  rbind_OS_SURVIVAL_final_patientid_mutated_selected$OS_MONTHS<<-as.numeric(rbind_OS_SURVIVAL_final_patientid_mutated_selected$OS_MONTHS)
  tryCatch({
      library(survival)
     # attach(rbind_OS_SURVIVAL_final_patientid_mutated_selected,warn.conflicts = FALSE)
      kmfit1 <- survfit(Surv(OS_MONTHS,OS_STATUS=='DECEASED')~GROUP,data = rbind_OS_SURVIVAL_final_patientid_mutated_selected)
      summary(kmfit1)
      plot(kmfit1)
  },error=function(x){message(x)},finally = function(x){
    #detach(rbind_OS_SURVIVAL_final_patientid_mutated_selected)
  }
  )
  
  tryCatch({
    library(survminer)
    library(ggplot2)
      q2=ggsurvplot(kmfit1,conf.int = TRUE, pval = TRUE, pval.method = TRUE,
                    risk.table = TRUE)
      plot(q2[["plot"]]);plot(q2[["table"]],ylim = 3)
      
      #print(survdiff(Surv(OS_MONTHS,OS_STATUS=='DECEASED')~GROUP,data = rbind_OS_SURVIVAL_final_patientid_mutated_selected ))
  },error=function(x){message(x)},finally = function(x){
    #detach(rbind_OS_SURVIVAL_final_patientid_mutated_selected)
  }
  )
  
  setwd(PATH_PROJECT)  
}
##载入合并的临床信息
load("data_clinical_all.save")
for(index in c(1:nrow(FILE_EXPRESSION))){
  gogogo(index)
}
