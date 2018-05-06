PATH_PROJECT="E:/src/r/SIRT_PTEN"
PATH="E:/data/tcga/"
DIRS=dir(path=PATH,pattern="*tcga$",full.names = TRUE,recursive = FALSE,ignore.case=TRUE)
FILE_EXPRESSION=data.frame()
for(onedir in  DIRS){
  one=list.files(path = onedir,pattern = "data.*?RNA_Seq.*?v2.*?expression.*?median.*?txt",recursive = TRUE,ignore.case = TRUE,full.names = TRUE)
  print(onedir)
  parent_file=stringr::str_extract(one,pattern = "^.*/")
  FILE_EXPRESSION[stringr::str_split(string = onedir,pattern = "/",simplify = TRUE)[4],"PATH"]=one
  FILE_EXPRESSION[stringr::str_split(string = onedir,pattern = "/",simplify = TRUE)[4],"PARENT_FILE"]=parent_file
}
##开始创建RDS文件
apply(FILE_EXPRESSION,MARGIN = 1,FUN = function(x){
  #fix(x)
  x=t(x)
  #fix(x)
  file=x[1,"PATH"]
  RDSfile=sub(pattern = "txt",replacement = "RDS",perl=TRUE,x=file)
  if(file.exists(RDSfile)){
#    rna_seq_data=readRDS(file=RDSfile)
  }else{
    message(paste("开始创建",RDSfile))
    rna_data=read.csv(file = file,sep = "\t")
    colnames(rna_data)=gsub(pattern = "\\.",replacement = "-",ignore.case = TRUE,perl = TRUE,x=colnames(rna_data))
    saveRDS(rna_data,file=RDSfile)
    rna_seq_data=readRDS(file=RDSfile)
    message(paste("成功生成RDS",RDSfile))
  }
  setwd(PATH_PROJECT)
})








