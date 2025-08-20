library(dplyr)
library(tidyr)
library(reshape2)
library(readxl)

exp = "ST218"  

# reading in raw count data
c_ST218=read.table("SPLINTR_counts_NGS_ST218.txt",header=FALSE)
colnames(c_ST218) = c("sample", "barcode_index", "counts", "frequency", "barcode", "colours")
c_ST218= unite(c_ST218, barcode, barcode, colours, sep="_")

# add samples names
file_names=read.table("filenames.txt",header=FALSE)
file_names= separate(file_names, V2, c("path", "sample_name"), sep="a/")
file_names= separate(file_names, sample_name, c("sample_name", "file"), sep="_S")
names= select(file_names, V1, sample_name)
colnames(names)= c("sample", "sample_name")
c_ST218= left_join(c_ST218, names, by= "sample")

#I need to split c_ST218 table because it's too big for dcast
split_table <- names$sample_name
split_table.1 <- split_table[1:100]
split_table.2 <- split_table[101:200]
split_table.3 <- split_table[201:320]

c_ST218 <- filter(c_ST218, !sample_name== "Undetermined")

c_ST218.1= dcast(c_ST218, barcode~sample_name, value.var = "counts", fill = 0)
write.table(c_ST218.1, file="ST218_raw_total_counts.txt", sep="\t", row.names=TRUE, col.names=TRUE)
#c_ST218.1=read.table("ST218_raw_total_counts.txt",header=TRUE)

#import reference
GFP=read.table("/stornext/General/data/academic/lab_naik/Sara_Tomei/R_analysis/Barcode_Analysis/SPLINTR_reference/GFP_reference.txt",header=TRUE)
BFP=read.table("/stornext/General/data/academic/lab_naik/Sara_Tomei/R_analysis/Barcode_Analysis/SPLINTR_reference/BFP_reference.txt",header=TRUE)
mCh=read.table("/stornext/General/data/academic/lab_naik/Sara_Tomei/R_analysis/Barcode_Analysis/SPLINTR_reference/mCh_reference.txt",header=TRUE)

GFP$protein <- "GFP"
BFP$protein <- "BFP"
mCh$protein <- "mCh"
GFP= unite(GFP, barcode, barcode, protein, sep="_")
BFP= unite(BFP, barcode, barcode, protein, sep="_")
mCh= unite(mCh, barcode, barcode, protein, sep="_")
GFP$barcode_num <- rownames(GFP)
BFP$barcode_num <- rownames(BFP)
mCh$barcode_num <- rownames(mCh)

all <- rbind(GFP, BFP, mCh)

#filter barcodes with reference
ST218 <-  semi_join(c_ST218.1, all, by = "barcode")

ST218.1 <- ST218
ST218.1[ST218.1<2]=0
ST218.1=ST218.1[rowSums(ST218.1[,2:321])>0,]

# calculate total read per sample
#detach(package:plyr)
#detach(package:tidyr)
total.count = all_barcodes %>% group_by(protein) %>% colSums(all_barcodes[,3:482]>0)
write.table(total.count, file=paste(exp, "raw_total_counts.txt", sep="_"), sep="\t", row.names=TRUE, col.names=TRUE)

# remove sample with low total read count 
read.threshold=1000
d=all_barcodes.1[,colSums(all_barcodes.1[,3:482])>read.threshold]
names(all_barcodes.1[,colSums(all_barcodes.1[,3:482])<=read.threshold])

#Put barcode and protein as row name
rownames(ST218.1)= ST218.1$barcode
ST218.1= select(ST218.1, -barcode)

rownames(ST208.1)= ST208.1$barcode
ST208.1= select(ST208.1, -barcode)

# remove samples low pearson correlation between tech.rep.
c=ST208.1
#c=ST218.1
pearson.threshold=0.6
Samples.reads=unique(substr(names(c),1,nchar(names(c))-2))
Samples.reads.pearson=vector()

for (s in Samples.reads)
{
  e=as.data.frame(c[,grep(s,names(c),value =TRUE )])
  if(ncol(e)<2)next;
  if(cor(e,method = "pearson")[1,2]>pearson.threshold) {Samples.reads.pearson=c(Samples.reads.pearson,s)}
  else print(c(s,cor(e,method = "pearson")[1,2]))
}


# set counts to 0 when barcodes with reads in one but the other
d.clean=data.frame(row.names=row.names(c))
for (s in Samples.reads.pearson)
{
  e=as.data.frame(c[,grep(s,names(c))])
  if(ncol(e)<2)next;
  e[rowSums(e>0)<2,]=0
  d.clean=cbind(d.clean,e)
}


# average tech. rep.
d.avg=data.frame(row.names=row.names(d.clean))
new.names=vector()
for (s in Samples.reads.pearson)
{
  e=as.data.frame(d.clean[,grep(s,names(d.clean))])
  if(ncol(e)<2)next;
  new.names=c(new.names,substr(names(e)[1],1,nchar(names(e)[1])-2))
  d.avg=cbind(d.avg,rowSums(e)/2)
}

names(d.avg)=new.names

ST218.avg=d.avg

#Devide different samples per patient
#avg <- ST218.avg
avg <- ST208.avg
mCh= avg[grepl("mCh", rownames(avg)),]
GFP= avg[grepl("GFP", rownames(avg)),]
BFP= avg[grepl("BFP", rownames(avg)),]

#Normalise by cell number sorted
ST218_numbers <- read_excel("ST218_numbers.xlsx")
Cell_number <- ST218_numbers[2,]


mCh=sweep(mCh,2,colSums(mCh),`/`)
GFP=sweep(GFP,2,colSums(GFP),`/`)
BFP=sweep(BFP,2,colSums(BFP),`/`)

d.norm.cell.1= rbind(mCh, GFP, BFP)

cell_number_norm <- function (n) {
  cell_name <- unique(names(Cell_number))
  a= data.frame(row.names=row.names(n))
  for (c in cell_name) {
    x=Cell_number[[1,c]]
    b= as.data.frame(n[,grep(c, names(n))]*x)
    names(b)= c
    a=cbind(a,b)
  }
  a=a/2
  return(a)
}

d.norm.cell=cell_number_norm(d.norm.cell.1)

ST218.norm.cell <- d.norm.cell

ST218.norm.cell$barcode=rownames(ST218.norm.cell)


ST218.norm.cell_PA <- ST218.norm.cell[,grepl("barcode", names(ST218.norm.cell)) |grepl("PA", names(ST218.norm.cell))]
ST218.norm.cell_PA=ST218.norm.cell_PA[rowSums(ST218.norm.cell_PA[,1:80])>0,]
ST218.norm.cell_PB <- ST218.norm.cell[,grepl("barcode", names(ST218.norm.cell)) |grepl("PB", names(ST218.norm.cell))]
ST218.norm.cell_PB=ST218.norm.cell_PB[rowSums(ST218.norm.cell_PB[,1:80])>0,]
common_barcodes.cell <- inner_join(ST218.norm.cell_PA, ST218.norm.cell_PB, by= "barcode")
common_barcodes.cell[is.na(common_barcodes.cell)]=0
rownames(common_barcodes.cell) <- common_barcodes.cell$barcode
common_barcodes.cell= select(common_barcodes.cell, -barcode)
common_barcodes.cell=common_barcodes.cell[rowSums(common_barcodes.cell)>0,]

# normalize to 10^6
ST218.cpm.norm= d.norm.cell.1*1e6

ST218.cpm.norm$barcode=rownames(ST218.cpm.norm)

ST218.cpm.norm_PA <- ST218.cpm.norm[,grepl("barcode", names(ST218.cpm.norm)) |grepl("PA", names(ST218.cpm.norm))]
ST218.cpm.norm_PA=ST218.cpm.norm_PA[rowSums(ST218.cpm.norm_PA[,1:80])>0,]
ST218.cpm.norm_PB <- ST218.cpm.norm[,grepl("barcode", names(ST218.cpm.norm)) |grepl("PB", names(ST218.cpm.norm))]
ST218.cpm.norm_PB=ST218.cpm.norm_PB[rowSums(ST218.cpm.norm_PB[,1:80])>0,]
common_barcodes <- inner_join(ST218.cpm.norm_PA, ST218.cpm.norm_PB, by= "barcode")
common_barcodes[is.na(common_barcodes)]=0
rownames(common_barcodes) <- common_barcodes$barcode
common_barcodes= select(common_barcodes, -barcode)
common_barcodes=common_barcodes[rowSums(common_barcodes)>0,]

write.table(common_barcodes, file=paste(exp, "common_barcodes_colum_norm.txt", sep="_"), sep="\t", row.names=FALSE, col.names=TRUE)

split_sample <- function (data, m1, m2){
  m.1=data%>%select(matches(m1));m.2=data%>%select(matches(m2))
  names(m.1)=gsub(m1,"",names(m.1));names(m.2)=gsub(m2,"",names(m.2))
  row.names(m.1)=paste(m1,row.names(m.1))
  row.names(m.2)=paste(m2,row.names(m.2))
  m.1=m.1[rowSums(m.1)>0,]
  m.2=m.2[rowSums(m.2)>0,]
  w=rbind(m.1,m.2,check.names=FALSE)
  w=w[rowSums(w)>0,]
  return(w)
}


#For SIS-seq
common.barcodes.cpm= split_sample(common_barcodes, "_1", "_2")
common.barcodes.cell= split_sample(common_barcodes.cell, "_1", "_2")
common.barcodes.cpm.plate= split_sample(common.barcodes.cpm, "_PA", "_PB")
common.barcodes.cell.plate= split_sample(common.barcodes.cell, "_PA", "_PB")


