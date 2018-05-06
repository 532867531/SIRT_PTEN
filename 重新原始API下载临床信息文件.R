FILENAME="gdc_manifest.2018-05-05.txt"
file=read.csv(file = FILENAME,sep = "\t")
url_base="https://api.gdc.cancer.gov/data/"
library(rvest)  
library(XML) 

myHttpheader <- c(
  "User-Agent"="Mozilla/5.0 (Windows; U; Windows NT 5.1; zh-CN; rv:1.9.1.6) ",
  "Accept"="text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
  "Accept-Language"="en-us",
  "Connection"="keep-alive",
  "Accept-Charset"="GB2312,utf-8;q=0.7,*;q=0.7"
)


paste(url_base,file$id[1],sep="")
a=read_xml(x=paste(url_base,file$id[1],sep=""))
