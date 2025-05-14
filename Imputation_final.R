library(data.table)
library(dplyr)
library(ggplot2)
library(MatchIt)
library(tidyr)
library(readxl)
library(ggrepel)
library(ggpubr)
library(seqinr)
library(cubar)
library(Biostrings)
library(cowplot)
library(ribor)
library(stringr)
library(limma)
library(edgeR)
library(treemap)
library(patchwork)
library(zCompositions)
library(foreach)
library(doParallel)



# these counts are after winsorization
ribo_count= read.csv("C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/FUS_2023/Imputation methods/ribo_raw_human_cap_995.csv")
rnaseq_count= read.csv("C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/FUS_2023/Imputation methods/rnaseq_raw_human_cap_995.csv")
polyA=read.csv("C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/FUS_2023/Imputation methods/polyA.csv")
dim(ribo_count)
dim(rnaseq_count)
dim(polyA)
ribo_count[1:5,1:5]
rnaseq_count[1:5,1:5]

# how many rows have rowSum =0

sum(rowSums(ribo_count[,2:1078])==0)
# remove thsese rows beacuse it is 0 across samples
ribozero= which(rowSums(ribo_count[,2:1078])==0) #there are 1810 rows whose rowSums are 0 in ribo counts


#same for RNA
sum(rowSums(rnaseq_count[,2:1078])==0)
rnazero=which(rowSums(rnaseq_count[,2:1078])==0) #there are 461 rows whose rowSums are 0 in rnaseq counts

# how many are zero in both
common_index= intersect(ribozero, rnazero)
length(common_index)

# 445 genes have no RNA orRibo values

ribo_count_clean= ribo_count[-common_index,]
dim(ribo_count_clean)
rna_count_clean= rnaseq_count[-common_index,]
dim(rna_count_clean)
rna_count_clean$X= sub("-\\d{3}$", "", rna_count_clean$X)


#which indexes have polyA
polyAgenes= which(rna_count_clean$X %in% polyA$GENE_true,)
length(polyAgenes)

# remove polyA genes
ribo_count_clean$X= sub("-\\d{3}$", "", ribo_count_clean$X)


#calculate the cpm to remove the genes
#has both RNA and Ribo
new_table= cbind(ribo_count_clean[,2:1078], rna_count_clean[,2:1078])
dim(new_table)

calculate_cpm <- function(expr_matrix) {
  # Calculate the column sums (total counts per sample)
  col_sums <- colSums(expr_matrix)
  
  # Calculate CPM
  cpm <- sweep(expr_matrix, 2, col_sums, FUN = "/") * 1e6
  
  return(cpm)
}

new_table_cpm= calculate_cpm(new_table) 
colSums(new_table_cpm)

# Figure ut how many rows have cpm values less than 1

threshold <- 0.2 * ncol(new_table_cpm)


# Identify rows where more than 20% of the columns have values below 1, rowsums of logical vector (TRue=1) and False=0, TRUE if it more80% samples for each gene is cpm <1
rows_result <- rowSums(new_table_cpm < 1) > threshold
sum(rows_result) #around 10830 genes have  genes who has cpm <1 for 20% samples in both RNA and Ribo togther
# add gene names
names(rows_result)=ribo_count_clean$X

sum(names(rows_result) %in% polyA$GENE_true)
rows_result[names(rows_result) %in% polyA$GENE_true]=TRUE
length(rows_result[polyA$GENE_true])

# subset these from the 
new_table_dummy= new_table[rows_result==TRUE,]

dim(new_table_dummy)

dim(new_table)

# Sum the values column-wise for the selected rows
dummy_gene <- colSums(new_table_dummy)
head(dummy_gene)


# remove these rows from the original table

common_indices <- which(rownames(new_table) %in% rownames(new_table_dummy))
length(common_indices)

# Remove rows from table1 based on common indices
Ribo_RNA <- new_table[-common_indices, ]
dim(Ribo_RNA)

# we have 8432 genes 
# Add dummy to this
Ribo_RNA_dummy= rbind(Ribo_RNA, dummy_gene)

dim(Ribo_RNA_dummy)
Ribo_RNA_dummygene =Ribo_RNA_dummy


# there are many samples whose geometric mean is screwed because of the dummy gene i.e most of the reads of the sample make the dummy gene. Hence remove these samples
hist(as.numeric(Ribo_RNA_dummygene[8433,])/colSums(Ribo_RNA_dummygene))

sum(as.numeric(Ribo_RNA_dummygene[8433,])/colSums(Ribo_RNA_dummygene) >0.3)
# there are 151 samples in which the dummy genes takes more than 70% of reads.remove those sample
dummy_cols= which(as.numeric(Ribo_RNA_dummygene[8433,])/colSums(Ribo_RNA_dummygene) >0.3)

#remove these from the table
Ribo_RNA_dummygene= Ribo_RNA_dummygene[, -dummy_cols]
dim(Ribo_RNA_dummygene)

#2003 samples instead of 2155

# add the gene_name

# Merge table1 and table2 by the common 'ID' column
# Add a column to table1 from table2 based on row names
Ribo_RNA_dummy=Ribo_RNA_dummy[,-dummy_cols]
Ribo_RNA_dummy$gene <- ifelse(rownames(Ribo_RNA_dummy) %in% rownames(ribo_count_clean),
                              ribo_count_clean[rownames(Ribo_RNA_dummy), "X"],
                              "dummy_gene")

dim(Ribo_RNA_dummy)

dim(Ribo_RNA_dummygene)
tail(Ribo_RNA_dummy$gene)
length(Ribo_RNA_dummy$gene)



Ribo_RNA_dummygene[1:5,1:5]



# transpose as mutiplicative works  within  a samples across genes 
Ribo_RNA_dummygene=t(Ribo_RNA_dummygene)
Ribo_RNA_dummygene[1:5,1:5]

# Impute  values using GBM
dim(Ribo_RNA_dummygene)



# Set up the cluster using all available cores
num_cores <- detectCores() -1  # Detect number of cores
cl <- makeCluster(num_cores)  # Create a cluster with the number of cores

# Use parApply to apply cmultRepl function in parallel

# Split the data into smaller chunks to distribute across the cores
num_rows <- nrow(Ribo_RNA_dummygene)
chunk_size <- ceiling(num_rows/num_cores)
row_indices <- split(seq_len(num_rows), ceiling(seq_len(num_rows) / chunk_size))
# Export Ribo_RNA_dummygene to the cluster
clusterExport(cl, "Ribo_RNA_dummygene")
clusterEvalQ(cl, library(zCompositions)) 


# Apply cmultRepl in parallel to each chunk
results <- parLapply(cl, row_indices, function(indices) {
  # Select the relevant subset of columns
  chunk_data <- Ribo_RNA_dummygene[indices, , drop = FALSE]
  
  # Apply the cmultRepl function
  cmultRepl(chunk_data, method = "GBM", output = "p-count")
})

# Stop the cluster after processing
stopCluster(cl)

# Combine the results back into a single table (if necessary)
final_result <- do.call(rbind, results)

# Print the result
dim(final_result)
final_result[1:5,1:5]

write.csv(final_result, "C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/FUS_2023/Imputation methods/ribo_rna_GBM_0.8.csv")

#if reading #
#final_result= read.csv("C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/FUS_2023/Imputation methods/ribo_rna_GBM_0.8.csv")
final_result[1:5,1:5]
#

# index as column, ist columninto rownames
#row.names(final_result)=final_result[,1]
#final_result=final_result[,-1]
# two columns have been removed

# Calculate clr for each sample

final_result=as.matrix(final_result)

# write a function for  calculating clr

calculate_clr <- function(df) {
  clr_result <- matrix(NA, nrow = nrow(df), ncol = ncol(df))
  colnames(clr_result) <- colnames(df)
  rownames(clr_result) <- rownames(df)
  for (i in 1:nrow(df)) {
    # Get the composition values (the row) and compute the geometric mean
    composition <- df[i, , drop = FALSE]  # Get the i-th row (a vector)
    
    # Calculate the geometric mean of the row (composition)
    geom_mean <- exp(mean(log(composition)))  # Geometric mean of the row
    
    # Calculate CLR for each component in the composition
    clr_result[i, ] <- log(composition / geom_mean)
  }
  return(as.data.frame(clr_result))
}




#clr_GBM= calculate_clr(final_result)
#write.csv(clr_GBM, "C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/FUS_2023/Imputation methods/clr_GBM.csv")
clr_GBM= read.csv("C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/FUS_2023/Imputation methods/clr_GBM.csv")
clr_GBM[1:5,1:5]
dim(clr_GBM)
row.names(clr_GBM)= clr_GBM[,1]
clr_GBM=clr_GBM[,-1]


# each row is a sample and each column is a gene
GBM_ilr= clr2ilr(clr_GBM)

GBM_ilr= read.csv("C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/FUS_2023/Imputation methods/GBM_ilr.csv")
#write.csv(GBM_ilr, "C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/FUS_2023/Imputation methods/GBM_ilr.csv")
dim(clr_GBM)
dim(GBM_ilr)


# Calculate proportional regression

write.csv(GBM_ilr, "C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/FUS_2023/Imputation methods/GBM_ilr.csv")

#split clr GBM to RNA and Ribo components
dim(clr_GBM)
clr_GBM[1:5,1:5]
clr_GBM_ribo= clr_GBM[1:995,]
clr_GBM_rna= clr_GBM[996:2001,]
clr_GBM_ribo[990:995, 1:5]
# while imputing , atotal of 2 samples was removed but both from ribo samples such that we have 1075 in Ribo and 1077
# remove these samples from RNA as well. 
#convert samples as rownames
#row.names(clr_GBM_ribo)=clr_GBM_ribo[,1]
#clr_GBM_ribo=clr_GBM_ribo[,-1]
#clr_GBM_ribo[1:5,1:5]
#row.names(clr_GBM_rna)=clr_GBM_rna[,1]
#clr_GBM_rna=clr_GBM_rna[,-1]
clr_GBM_rna[1:5,1:5]


head(row.names(clr_GBM_ribo))
head(row.names(clr_GBM_rna))
rownames(clr_GBM_ribo) <- gsub("^\\d+\\.", "", rownames(clr_GBM_ribo))
rownames(clr_GBM_rna) <- gsub("^\\d+\\.(GSM\\d+)\\..*$", "\\1",  rownames(clr_GBM_rna))
# remove the column that is not common in both 


dim(clr_GBM_rna)
dim(clr_GBM_ribo)

# Find the row names that are in dt2 but not in dt1
extra_samples <- setdiff(rownames(clr_GBM_rna), rownames(clr_GBM_ribo))
extra_samples_ribo<- setdiff(rownames(clr_GBM_ribo), rownames(clr_GBM_rna))
length(extra_samples_ribo)

# Remove the extra rows from dt2 (2)
clr_GBM_rna<- clr_GBM_rna[!rownames(clr_GBM_rna) %in% extra_samples, ]
clr_GBM_ribo<- clr_GBM_ribo[!rownames(clr_GBM_ribo) %in% extra_samples_ribo, ]

GBM_ilr[1:5,1:5]
dim(GBM_ilr)
# if red,scv
# mkae the first column as the row
row.names(GBM_ilr)= GBM_ilr[,1]
GBM_ilr=GBM_ilr[,-1]


# SImilarly split RNA and Ribo for ilr calculated 

#row.names(GBM_ilr)=clr_GBM[,1]
ilr_GBM_ribo= GBM_ilr[1:995,]
ilr_GBM_rna= GBM_ilr[996:2001,]
dim(ilr_GBM_ribo)
dim(ilr_GBM_rna)
ilr_GBM_ribo[990:995, 1:5]
ilr_GBM_rna[995:1006, 1:5]
# while imputing , atotal of 2 samples was removed during imputation ans few were removed while putting a threshold based on dummy genes
# remove these samples from RNA as well. 

tail(row.names(ilr_GBM_ribo))
tail(row.names(ilr_GBM_rna))
rownames(ilr_GBM_ribo) <- gsub("^\\d+\\.", "", rownames(ilr_GBM_ribo))
rownames(ilr_GBM_rna) <- gsub("^\\d+\\.(GSM\\d+)\\..*$", "\\1", rownames(ilr_GBM_rna))
# remove the column that is not common in both 


dim(ilr_GBM_rna)
dim(ilr_GBM_ribo)

# Find the row names that are in dt2 but not in dt1
extra_samples <- setdiff(rownames(ilr_GBM_rna), rownames(ilr_GBM_ribo))
length(extra_samples)
extra_samples_ribo<- setdiff(rownames(ilr_GBM_ribo), rownames(ilr_GBM_rna))
length(extra_samples_ribo)

# Remove the extra rows from dt2 #12
ilr_GBM_rna<- ilr_GBM_rna[!rownames(ilr_GBM_rna) %in% extra_samples, ]
ilr_GBM_ribo<- ilr_GBM_ribo[!rownames(ilr_GBM_ribo) %in% extra_samples_ribo, ]


# transpose to calculate column wise linear regression value 
ilr_GBM_rna=t(ilr_GBM_rna)
ilr_GBM_ribo=t(ilr_GBM_ribo)




# calculate the linear regresion 

library(foreach)
library(doParallel)




# Register parallel backend
num_cores <- detectCores() - 10  # Use one less core than available
cl <- makeCluster(num_cores)    # Create a cluster
registerDoParallel(cl)          # Register the cluster

# Perform parallel loop
out <- foreach(i = 1:ncol(ilr_GBM_ribo), .combine = "cbind", .packages = c("compositions")) %dopar% {
  # Fit the linear model for each column
  m <- summary(lm(ilr_GBM_ribo[, i] ~ ilr_GBM_rna[, i]))  # Linear regression
  
  # Extract residuals and transform from ILR to CLR
  data.frame(as.numeric(ilr2clr(resid(m))))  # Store transformed residuals in CLR space
}

# Stop the cluster after the computation
stopCluster(cl)

dim(out)
#write.csv(out, "C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/FUS_2023/Imputation methods/TE_GBM_clr.csv")
#TE_GBM_clr=out
TE_GBM_clr= read.csv("C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/FUS_2023/Imputation methods/TE_GBM_clr.csv")
# if read as csv ,then follwong code
TE_GBM_clr=TE_GBM_clr[,-1]
TE_GBM_clr[1:5,1:5]
dim(TE_GBM_clr)

dim(TE_GBM_clr)
TE_GBM_clr$gene= Ribo_RNA_dummy$gene
TE_GBM_clr[1:5,1:5]

TE_GBM_clr <- cbind(gene = TE_GBM_clr$gene, TE_GBM_clr[,1:994])

#row.names(TE_GBM_clr)= Ribo_RNA_dummy$gene
#TE_GBM_clr[1:5,1:5]
#TE_GBM_clr <- cbind(gene = rownames(TE_GBM_clr), TE_GBM_clr)
#row.names(TE_GBM_clr)= NULL

#Add sample names to column _# as two columns are removed after GBM. The 
# Using names() to assign column names

# Assign row names from clr_GBM_ribo
colnames(TE_GBM_clr) <- c("gene",row.names(clr_GBM_ribo))
TE_GBM_clr[1:5,1:5]
#TE_GBM_clr$gene <- sub("-2.*$", "", TE_GBM_clr$gene)
TE_GBM_clr=as.data.table(TE_GBM_clr)
colnames(TE_GBM_clr)[1] <- "transcript"

#do the same for RNA & Ribo 

clr_GBM_ribo[1:5,1:5]


Ribo_clr= t(clr_GBM_ribo)
Ribo_clr=as.data.table(Ribo_clr)
# provide genenames 
Ribo_clr[1:5,1:5]
Ribo_clr$transcript= Ribo_RNA_dummy$gene
Ribo_clr <- cbind(transcript=Ribo_clr$transcript, Ribo_clr[,1:994])
dim(Ribo_clr)
Ribo_clr[1:5,1:5]




clr_GBM_rna[1:5,1:5]

RNA_clr= t(clr_GBM_rna)
RNA_clr=as.data.table(RNA_clr)
# provide genenames 
RNA_clr[1:5,1:5]
RNA_clr$transcript= Ribo_RNA_dummy$gene
RNA_clr <- cbind(transcript=RNA_clr$transcript, RNA_clr[,1:994])
dim(RNA_clr)
RNA_clr[1:5,1:5]

# from all the samples remove samples which do not fit the line  of curve

# if use clr, transpose

clr_GBM_rna[1:5,1:5]
clr_GBM_rna_t=t(clr_GBM_rna)
clr_GBM_ribo[1:5,1:5]
clr_GBM_ribo_t=t(clr_GBM_ribo)
clr_GBM_ribo_t[1:5,1:5]
clr_GBM_rna_t[1:5,1:5]

num_cores <- detectCores() - 10  # Use one less core than available
cl <- makeCluster(num_cores)    # Create a cluster
registerDoParallel(cl)          # Register the cluster

# Perform parallel loop
r2_clr <- foreach(i = 1:ncol(clr_GBM_ribo_t), .combine = "cbind", .packages = c("compositions")) %dopar% {
  # Fit the linear model for each column
  m <- summary(lm(clr_GBM_ribo_t[, i] ~ clr_GBM_rna_t[, i]))  # Linear regression
  m$r.squared 
  
}

# Stop the cluster after the computation
stopCluster(cl)
hist(r2_clr)

r2_clr_dt= as.data.table(as.numeric(r2_clr))
r2_clr_dt$sample= colnames(clr_GBM_ribo_t)

# remove samples which have r2 =< 0.2
hist(r2_clr_dt$V1)
sum(r2_clr_dt$V1 > 0.2)
r2_clr_dt_subset= subset(r2_clr_dt, r2_clr_dt$V1 <0.2)
dim(r2_clr_dt_subset)

# remove samples which have less samples less than 0.2 R2 # this removes 30 samples
sampletoremove= r2_clr_dt_subset$sample

RNA_clr= RNA_clr[,c(sampletoremove) := NULL]
dim(RNA_clr)

# similarly for ribo

Ribo_clr= Ribo_clr[,c(sampletoremove) := NULL]
dim(Ribo_clr)

#similarly for TE
TE_GBM_clr[1:5,1:5]
TE_GBM_clr= TE_GBM_clr[,c(sampletoremove) := NULL]
dim(TE_GBM_clr)

RNA_GBM=write.csv(RNA_clr, "C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/FUS_2023/Imputation methods/RNA_clr.csv")
Ribo_GBM=write.csv(Ribo_clr, "C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/FUS_2023/Imputation methods/Ribo_clr.csv")
TE_GBM=write.csv(TE_GBM_clr, "C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/FUS_2023/Imputation methods/TE_GBM_clr.csv")


