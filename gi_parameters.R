#!/usr/bin/env Rscript


setwd("path/to/your/workingdir")

### CHANGE PATHS AND TUMOR BURDEN !!!!!!

# libraries
require(plyr)
require(stringr)
library("clipr")
library(WriteXLS)


tam_gen <- 3200


## Read files and create row and col names

list_names <- list.files(pattern = "pattern")
names <- str_remove_all(list_names, "pattern")

for (i in 1:length(names)){
  
  file<-read.delim(list_names[[i]], header = TRUE, sep = "\t", dec = ".")
  
  file$length<-(file$end-file$start)/1000000
  
  file1<-file[file$length>=1, ]
  
  
  ######## annotate cns files obtained with cnvkit #######
  
  
  file2m<- file1[file1$weight>=median(file1$weight), ]
  
  
  file2m$annotation <- ifelse((file2m$cn == 1) & (file2m$cn1==1) & (file2m$cn2 ==0), "loh","undet")
  file2m$annotation[(file2m$cn == 2) & (file2m$cn1==2) & (file2m$cn2 ==0) & (file2m$baf > 0.3) & (file2m$baf < 0.7)] <- "cn-loh"
  file2m$annotation[(file2m$cn >= 3)] <- "gain"
  file2m$annotation[(file2m$cn == 1) & is.na(file2m$cn1) & is.na(file2m$cn2)] <- "loss"
  file3m<-file2m[file2m$annotation != "undet", ]
  
  file4m<- file3m[, -c(13:14)]
  
  
  nam_var1 <- paste(names[[i]],".csv", sep="")
  nam_var2 <- paste(names[[i]],".cns", sep="")
  
  
  write.table(file3m, nam_var1,
              append = FALSE, sep="\t", dec = ".", row.names = FALSE, col.names = TRUE)
  
  write.table(file4m, nam_var2,
              append = FALSE, sep="\t", dec = ".", row.names = FALSE, col.names = TRUE)
}

cp_files <- list.files(pattern = "pattern")


file.copy(from=cp_files, to="/path/to/copy>/dir", 
          overwrite = TRUE, recursive = FALSE, 
          copy.mode = TRUE)

file.remove(cp_files)


setwd("/path/to/created/files")



param_row <-c("Eventos_totales", "Media_evento", "eventos_mediana",
              "mb_totales", "pgenoma_alt", "Eventos_totales_sinCN3(0.3-0.5)",
              "eventoscn3_mediana", "Media_evento_sinCN3(0.3-0.5)", "mb_totales_sinCN3(0.3-0.5)",
              "pgenoma_alt_sinCN3(0.3-0.5)", "ganancias_totales", "mb_ganancias",
              "pgenoma_ganancias", "ganancias_totalesCN3", "mb_gananciasCN3",
              "pgananciasCN3", "perdidas_totales", "mb_perdidas",
              "pgenoma_perdidas", "loh_totales", "Mediana_loh",                    
              "Media_loh", "mb_loh", "pgenoma_loh",
              "loh>15mb_totales","mb_loh>15mb", "%genoma_LOH>15",
              "LOHmas10MB", "PLOHmas10mb", "mbLOHmas10mb"
)



my_list <- lapply(list_names, read.csv, stringsAsFactors = FALSE, sep="\t", header=T)
my_list <- lapply(my_list, function(x) { x["gene"] <- NULL; x })

## initialize matrix
param_list <- list()


for(j in 1:length(names)){
  purity <- c(
# Fill tumor burden of the samples
  )
  
  ## Correction factor to establish parameters between 0.5 and 3 of BAF
  fac_correc1 <- round(log(((1-purity[j])+purity[j]*3.5/2),2),4)
  fac_correc2 <- round(log(((1-purity[j])+purity[j]*2.5/2),2),4)
  fac_correc <- (((fac_correc1 - fac_correc2)/2) + fac_correc2)
  ## variables (n=30)
  xx<- length(my_list[[j]]$log2)
  xxx <- mean(my_list[[j]]$length)
  xxxx <- median(my_list[[j]]$length)
  a<-round(sum(my_list[[j]]$length),3)
  b<-((round(sum(my_list[[j]]$length),3))/tam_gen)*100
  c <- length(my_list[[j]]$log2[!(my_list[[j]]$annotation=="gain" & my_list[[j]]$log2<fac_correc)])
  d <- median(my_list[[j]]$length[!(my_list[[j]]$annotation=="gain" & my_list[[j]]$log2<fac_correc)])
  e <- mean(my_list[[j]]$length[!(my_list[[j]]$annotation=="gain" & my_list[[j]]$log2<fac_correc)])
  f <- sum(my_list[[j]]$length[!(my_list[[j]]$annotation=="gain" & my_list[[j]]$log2<fac_correc)])
  g <- ((round(sum(my_list[[j]]$length[!(my_list[[j]]$annotation=="gain" & my_list[[j]]$log2<fac_correc)]),3))/tam_gen)*100
  h <- length(my_list[[j]]$log2[my_list[[j]]$annotation=="gain"])
  k <- sum(my_list[[j]]$length[ my_list[[j]]$annotation=="gain"])
  l <- ((round(sum(my_list[[j]]$length[my_list[[j]]$annotation=="gain"]),3))/tam_gen)*100
  m <- length(my_list[[j]]$log2[my_list[[j]]$annotation=="gain" & my_list[[j]]$cn>3])
  n <- sum(my_list[[j]]$length[my_list[[j]]$annotation=="gain" & my_list[[j]]$cn>3])
  o <- ((round(sum(my_list[[j]]$length[my_list[[j]]$annotation=="gain" & my_list[[j]]$cn>3]),3))/tam_gen)*100
  p <- length(my_list[[j]]$log2[my_list[[j]]$annotation=="loss"])
  q <- sum(my_list[[j]]$length[ my_list[[j]]$annotation=="loss"])
  r <- ((round(sum(my_list[[j]]$length[my_list[[j]]$annotation=="loss"]),3))/tam_gen)*100
  s <- length(my_list[[j]]$log2[my_list[[j]]$annotation=="loh"])
  t <- median(my_list[[j]]$length[my_list[[j]]$annotation=="loh"])
  u <- mean(my_list[[j]]$length[my_list[[j]]$annotation=="loh"])
  v <- sum(my_list[[j]]$length[ my_list[[j]]$annotation=="loh"])
  w <- ((round(sum(my_list[[j]]$length[my_list[[j]]$annotation=="loh"]),3))/tam_gen)*100
  x <- length(my_list[[j]]$log2[my_list[[j]]$annotation=="loh" & my_list[[j]]$length>15])
  y <- sum(my_list[[j]]$length[my_list[[j]]$annotation=="loh" & my_list[[j]]$length>15])
  z <- ((round(sum(my_list[[j]]$length[my_list[[j]]$annotation=="loh"& my_list[[j]]$length>15]),3))/tam_gen)*100
  aa <- length(my_list[[j]]$log2[my_list[[j]]$annotation=="loh" & my_list[[j]]$length>10])
  ab <- ((round(sum(my_list[[j]]$length[my_list[[j]]$annotation=="loh"& my_list[[j]]$length>10]),3))/tam_gen)*100
  ac <- sum(my_list[[j]]$length[my_list[[j]]$annotation=="loh" & my_list[[j]]$length>10])
  
  ## Temporal list
  name <- paste('case',j, sep='_')
  tmp <- list(xx, xxx, xxxx,
              a,b,c,
              d,e,f,
              g,h,k,
              l,m,n,
              o,p,q,
              r,s,t,
              u,v,w,
              x,y,z,
              aa,ab,ac)
  ## List creation
  param_list[[name]] <- tmp
  
}



# Creation of df with all the data
param_df <- data.frame(matrix(unlist(param_list), nrow=length(param_list), byrow=T))
param_df <- do.call(rbind.data.frame, param_list)
param_df <- t(param_df)
colnames(param_df) <- names
rownames(param_df) <- param_row
param_df <- t(param_df)
