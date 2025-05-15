ribo_count_mouse= read.csv("./ribo_raw_mouse_cap_995_input.csv")
rnaseq_count_mouse= read.csv("./rna_raw_mouse_cap_995_input.csv")
polyA_mouse=read.csv("./polyA_mouse.csv")
dim(ribo_count_mouse)
dim(rnaseq_count_mouse)
#21569 genes ans 846 samples

dim(polyA)
ribo_count_mouse[1:5,1:5]
rnaseq_count_mouse[1:5,1:5]

# how many rows have rowSum =0

sum(rowSums(ribo_count_mouse[,2:846])==0)
# remove thsese rows beacuse it is 0 across samples
ribozero= which(rowSums(ribo_count_mouse[,2:846])==0) #there are 2882 rows whose rowSums are 0 in ribo counts


#same for RNA
sum(rowSums(rnaseq_count_mouse[,2:846])==0)
rnazero=which(rowSums(rnaseq_count[,2:846])==0) #there are 590 rows whose rowSums are 0 in rnaseq counts

# how many are zero in both
common_index= intersect(ribozero, rnazero)
length(common_index)

# 77 genes have no RNA or Ribo values in common but a lot of genes have no count in Ribo. 
ribo_count_clean_mouse= ribo_count_mouse[-common_index,]
dim(ribo_count_clean_mouse)
rna_count_clean_mouse= rnaseq_count_mouse[-common_index,]
dim(rna_count_clean_mouse)


#which indexes have polyA
polyAgenes_mouse= which(rna_count_clean_mouse$transcript %in% polyA_mouse$GENE_true,)
length(polyAgenes_mouse)
#there are 61 genes that foun d in this list

##calculate the cpm to remove the genes
#has both RNA and Ribo
new_table_mouse= cbind(ribo_count_clean_mouse[,2:846], rna_count_clean_mouse[,2:846])
dim(new_table_mouse)

calculate_cpm <- function(expr_matrix) {
  # Calculate the column sums (total counts per sample)
  col_sums <- colSums(expr_matrix)
  
  # Calculate CPM
  cpm <- sweep(expr_matrix, 2, col_sums, FUN = "/") * 1e6
  
  return(cpm)
}

new_table_cpm_mouse= calculate_cpm(new_table_mouse) 
colSums(new_table_cpm_mouse)

# Figure ut how many rows have cpm values less than 1

threshold <- 0.2 * ncol(new_table_cpm_mouse)


# Identify rows where more than 20% of the columns have values below 1, rowsums of logical vector (TRue=1) and False=0, TRUE if it more 80% samples for each gene is cpm <1
rows_result_mouse <- rowSums(new_table_cpm_mouse < 1) > threshold
sum(rows_result_mouse) #around 13341 genes have  genes who has cpm <1 for 80% samples in both RNA and Ribo togther
# add gene names
names(rows_result_mouse)=ribo_count_clean_mouse$transcript

sum(names(rows_result_mouse) %in% polyA_mouse$GENE_true)
# make poy A genes as TRUE to be included as dummy genes
rows_result_mouse[names(rows_result_mouse) %in% polyA_mouse$GENE_true]=TRUE
length(rows_result_mouse[polyA_mouse$GENE_true])

# subset these from the 
new_table_dummy_mouse= new_table_mouse[rows_result_mouse==TRUE,]

dim(new_table_dummy_mouse) # there are 13358 genes that are dummy genes

dim(new_table_mouse)

# Sum the values column-wise for the selected rows
dummy_gene_mouse <- colSums(new_table_dummy_mouse)
head(dummy_gene_mouse)


# remove these rows from the original table 

common_indices_mouse <- which(rownames(new_table_mouse) %in% rownames(new_table_dummy_mouse))
length(common_indices_mouse)

# Remove rows from table1 based on common indices
Ribo_RNA_mouse <- new_table_mouse[-common_indices_mouse, ]
dim(Ribo_RNA_mouse)

# we have 8134 genes 
# Add dummy to this
Ribo_RNA_dummy_mouse= rbind(Ribo_RNA_mouse, dummy_gene_mouse)

dim(Ribo_RNA_dummy_mouse)
Ribo_RNA_dummygene_mouse =Ribo_RNA_dummy_mouse

# there are many samples whose geometric mean is screwed because of the dummy gene i.e most of the reads of the sample make the dummy gene. Hence remove these samples
hist(as.numeric(Ribo_RNA_dummygene_mouse[8135,])/colSums(Ribo_RNA_dummygene_mouse))

sum(as.numeric(Ribo_RNA_dummygene_mouse[8135,])/colSums(Ribo_RNA_dummygene_mouse) >0.3)
# there are 406 samples in which the dummy genes takes more than 80% of reads.remove those sample
dummy_cols_mouse= which(as.numeric(Ribo_RNA_dummygene_mouse[8135,])/colSums(Ribo_RNA_dummygene_mouse) >0.3)


#remove these from the table
Ribo_RNA_dummygene_mouse= Ribo_RNA_dummygene_mouse[, -dummy_cols_mouse]
dim(Ribo_RNA_dummygene_mouse)

#1284 samples instead of 1690

# add the gene_name

# Merge table1 and table2 by the common 'ID' column
# Add a column to table1 from table2 based on row names
Ribo_RNA_dummy_mouse=Ribo_RNA_dummy_mouse[,-dummy_cols_mouse]
Ribo_RNA_dummy_mouse$gene <- ifelse(rownames(Ribo_RNA_dummy_mouse) %in% rownames(ribo_count_clean_mouse),
                              ribo_count_clean_mouse[rownames(Ribo_RNA_dummy_mouse), "transcript"],
                              "dummy_gene_mouse")

dim(Ribo_RNA_dummy_mouse)

dim(Ribo_RNA_dummygene_mouse)
tail(Ribo_RNA_dummy_mouse$gene)
length(Ribo_RNA_dummy_mouse$gene)

Ribo_RNA_dummygene_mouse[1:5,1:5]
# transpose as mutiplicative works  across samples
Ribo_RNA_dummygene_mouse=t(Ribo_RNA_dummygene_mouse)
Ribo_RNA_dummygene_mouse[1:5,1:5]

# Impute  values using GBM
dim(Ribo_RNA_dummygene_mouse)


num_cores <- detectCores() -4  # Detect number of cores
cl <- makeCluster(num_cores)  # Create a cluster with the number of cores

# Use parApply to apply cmultRepl function in parallel

# Split the data into smaller chunks to distribute across the cores
num_rows <- nrow(Ribo_RNA_dummygene_mouse)
chunk_size <- ceiling(num_rows/num_cores)
row_indices <- split(seq_len(num_rows), ceiling(seq_len(num_rows) / chunk_size))
# Export Ribo_RNA_dummygene to the cluster
clusterExport(cl, "Ribo_RNA_dummygene_mouse")
clusterEvalQ(cl, library(zCompositions)) 


# Apply cmultRepl in parallel to each chunk
results <- parLapply(cl, row_indices, function(indices) {
  # Select the relevant subset of columns
  chunk_data <- Ribo_RNA_dummygene_mouse[indices, , drop = FALSE]
  
  # Apply the cmultRepl function
  cmultRepl(chunk_data, method = "GBM", output = "p-count")
})

# Stop the cluster after processing
stopCluster(cl)

# Combine the results back into a single table (if necessary)
final_result_mouse <- do.call(rbind, results)

# Print the result
dim(final_result_mouse)
final_result_mouse[1:5,1:5]

write.csv(final_result_mouse, "C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/FUS_2023/Imputation methods/ribo_rna_GBM_0.8_mouse.csv")

# Calculate clr for each sample

final_result_mouse=as.matrix(final_result_mouse)

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




clr_GBM_mouse= calculate_clr(final_result_mouse)
write.csv(clr_GBM_mouse, "./clr_GBM_mouse.csv")
clr_GBM_mouse= read.csv("./clr_GBM_mouse.csv")
dim(clr_GBM_mouse)
row.names(clr_GBM_mouse)= clr_GBM_mouse[,1]
clr_GBM_mouse=clr_GBM_mouse[,-1]
clr_GBM_mouse[1280:1284,1:5]

GBM_ilr_mouse= clr2ilr(clr_GBM_mouse)
# Print output dimensions
dim(GBM_ilr_mouse)
write.csv(GBM_ilr_mouse, "./GBM_ilr_mouse.csv")

dim(GBM_ilr_mouse)
GBM_ilr_mouse[1280:1284,1:5]


#split clr GBM to RNA and Ribo components
dim(clr_GBM_mouse)
clr_GBM_mouse[1:5,1:5]
clr_GBM_ribo_mouse= clr_GBM_mouse[1:599,]
clr_GBM_rna_mouse= clr_GBM_mouse[600:1284,]
clr_GBM_ribo_mouse[1:5, 1:5]
clr_GBM_rna_mouse[1:5, 1:5]
# while imputing , atotal of 2 samples was removed but both from ribo samples such that we have 1075 in Ribo and 1077
# remove these samples from RNA as well. 
#convert samples as rownames
#row.names(clr_GBM_ribo)=clr_GBM_ribo[,1]
#clr_GBM_ribo=clr_GBM_ribo[,-1]
#clr_GBM_ribo[1:5,1:5]
#row.names(clr_GBM_rna)=clr_GBM_rna[,1]
#clr_GBM_rna=clr_GBM_rna[,-1]



head(row.names(clr_GBM_ribo_mouse))
head(row.names(clr_GBM_rna_mouse))
rownames(clr_GBM_ribo_mouse) <- gsub("^\\d+\\.", "", rownames(clr_GBM_ribo_mouse))
rownames(clr_GBM_rna_mouse) <- gsub("^\\d+\\.(GSM\\d+)\\..*$", "\\1",  rownames(clr_GBM_rna_mouse))
# remove the column that is not common in both 


dim(clr_GBM_rna_mouse)
dim(clr_GBM_ribo_mouse)

# Find the row names that are in dt2 but not in dt1
extra_samples_mouse <- setdiff(rownames(clr_GBM_rna_mouse), rownames(clr_GBM_ribo_mouse))
extra_samples_ribo_mouse<- setdiff(rownames(clr_GBM_ribo_mouse), rownames(clr_GBM_rna_mouse))
length(extra_samples_ribo_mouse)
length(extra_samples_mouse)


# Remove the extra rows from dt2 (2)
clr_GBM_rna_mouse<- clr_GBM_rna_mouse[!rownames(clr_GBM_rna_mouse) %in% extra_samples_mouse, ]
clr_GBM_ribo_mouse<- clr_GBM_ribo_mouse[!rownames(clr_GBM_ribo_mouse) %in% extra_samples_ribo_mouse, ]

GBM_ilr_mouse[1:5,1:5]
dim(GBM_ilr_mouse)


# SImilarly split RNA and Ribo for ilr calculated 

#row.names(GBM_ilr)=clr_GBM[,1]
ilr_GBM_ribo_mouse= GBM_ilr_mouse[1:599,]
ilr_GBM_rna_mouse= GBM_ilr_mouse[600:1284,]
dim(ilr_GBM_ribo_mouse)
dim(ilr_GBM_rna_mouse)
ilr_GBM_ribo_mouse[1:5, 1:5]
ilr_GBM_rna_mouse[1:5, 1:5]
# while imputing , atotal of 2 samples was removed during imputation ans few were removed while putting a threshold based on dummy genes
# remove these samples from RNA as well. 

tail(row.names(ilr_GBM_ribo_mouse))
tail(row.names(ilr_GBM_rna_mouse))
rownames(ilr_GBM_ribo_mouse) <- gsub("^\\d+\\.", "", rownames(ilr_GBM_ribo_mouse))
rownames(ilr_GBM_rna_mouse) <- gsub("^\\d+\\.(GSM\\d+)\\..*$", "\\1", rownames(ilr_GBM_rna_mouse))
# remove the column that is not common in both 


dim(ilr_GBM_rna_mouse)
dim(ilr_GBM_ribo_mouse)

# Find the row names that are in dt2 but not in dt1
extra_samples_mouse <- setdiff(rownames(ilr_GBM_rna_mouse), rownames(ilr_GBM_ribo_mouse))
length(extra_samples_mouse)
extra_samples_ribo_mouse<- setdiff(rownames(ilr_GBM_ribo_mouse), rownames(ilr_GBM_rna_mouse))
length(extra_samples_ribo_mouse)

# Remove the extra rows from dt2 #12
ilr_GBM_rna_mouse<- ilr_GBM_rna_mouse[!rownames(ilr_GBM_rna_mouse) %in% extra_samples_mouse, ]
ilr_GBM_ribo_mouse<- ilr_GBM_ribo_mouse[!rownames(ilr_GBM_ribo_mouse) %in% extra_samples_ribo_mouse, ]


# transpose to calculate column wise linear regression value 
ilr_GBM_rna_mouse=t(ilr_GBM_rna_mouse)
ilr_GBM_ribo_mouse=t(ilr_GBM_ribo_mouse)


library(foreach)
library(doParallel)




# Register parallel backend
num_cores <- detectCores() - 10  # Use one less core than available
cl <- makeCluster(num_cores)    # Create a cluster
registerDoParallel(cl)          # Register the cluster

# Perform parallel loop
out_mouse <- foreach(i = 1:ncol(ilr_GBM_ribo_mouse), .combine = "cbind", .packages = c("compositions")) %dopar% {
  # Fit the linear model for each column
  m <- summary(lm(ilr_GBM_ribo_mouse[, i] ~ ilr_GBM_rna_mouse[, i]))  # Linear regression
  
  # Extract residuals and transform from ILR to CLR
  data.frame(as.numeric(ilr2clr(resid(m))))  # Store transformed residuals in CLR space
}

# Stop the cluster after the computation
stopCluster(cl)

dim(out_mouse)
write.csv(out_mouse, "./TE_GBM_clr_mouse.csv")
TE_GBM_clr_mouse= read.csv("./TE_GBM_clr_mouse.csv")
TE_GBM_clr_mouse[1:5,1:5]


TE_GBM_clr_mouse=TE_GBM_clr_mouse[,-1]
TE_GBM_clr_mouse[1:5,1:5]
dim(TE_GBM_clr_mouse)

dim(TE_GBM_clr_mouse)
TE_GBM_clr_mouse$gene= Ribo_RNA_dummy_mouse$gene
TE_GBM_clr_mouse[1:5,1:5]

TE_GBM_clr_mouse <- cbind(gene = TE_GBM_clr_mouse$gene, TE_GBM_clr_mouse[,1:559])

# Assign row names from clr_GBM_ribo
colnames(TE_GBM_clr_mouse) <- c("gene",row.names(clr_GBM_ribo_mouse))
TE_GBM_clr_mouse[1:5,1:5]
#TE_GBM_clr$gene <- sub("-2.*$", "", TE_GBM_clr$gene)
TE_GBM_clr_mouse=as.data.table(TE_GBM_clr_mouse)
colnames(TE_GBM_clr_mouse)[1] <- "transcript"

#do the same for RNA & Ribo 

clr_GBM_ribo_mouse[1:5,1:5]


Ribo_clr_mouse= t(clr_GBM_ribo_mouse)
Ribo_clr_mouse=as.data.table(Ribo_clr_mouse)
# provide genenames 
Ribo_clr_mouse[1:5,1:5]
Ribo_clr_mouse$transcript= Ribo_RNA_dummy_mouse$gene
Ribo_clr_mouse <- cbind(transcript=Ribo_clr_mouse$transcript, Ribo_clr_mouse[,1:559])
dim(Ribo_clr_mouse)
Ribo_clr_mouse[1:5,1:5]




clr_GBM_rna_mouse[1:5,1:5]

RNA_clr_mouse= t(clr_GBM_rna_mouse)
RNA_clr_mouse=as.data.table(RNA_clr_mouse)
# provide genenames 
RNA_clr_mouse[1:5,1:5]
RNA_clr_mouse$transcript= Ribo_RNA_dummy_mouse$gene
RNA_clr_mouse <- cbind(transcript=RNA_clr_mouse$transcript, RNA_clr_mouse[,1:559])
dim(RNA_clr_mouse)
RNA_clr_mouse[1:5,1:5]

# from all the samples remove samples which do not fit the line  of curve

# if use clr, transpose

clr_GBM_rna_mouse[1:5,1:5]
clr_GBM_rna_t_mouse=t(clr_GBM_rna_mouse)
clr_GBM_ribo_mouse[1:5,1:5]
clr_GBM_ribo_t_mouse=t(clr_GBM_ribo_mouse)
clr_GBM_ribo_t_mouse[1:5,1:5]
clr_GBM_rna_t_mouse[1:5,1:5]

num_cores <- detectCores() - 10  # Use one less core than available
cl <- makeCluster(num_cores)    # Create a cluster
registerDoParallel(cl)          # Register the cluster

# Perform parallel loop
r2_clr_mouse <- foreach(i = 1:ncol(clr_GBM_ribo_t_mouse), .combine = "cbind", .packages = c("compositions")) %dopar% {
  # Fit the linear model for each column
  m <- summary(lm(clr_GBM_ribo_t_mouse[, i] ~ clr_GBM_rna_t_mouse[, i]))  # Linear regression
  m$r.squared 
  
}

# Stop the cluster after the computation
stopCluster(cl)
hist(r2_clr_mouse)

r2_clr_dt_mouse= as.data.table(as.numeric(r2_clr_mouse))
r2_clr_dt_mouse$sample= colnames(clr_GBM_ribo_t_mouse)

# remove samples which have r2 =< 0.2
hist(r2_clr_dt_mouse$V1)
sum(r2_clr_dt_mouse$V1 > 0.2)
r2_clr_dt_subset_mouse= subset(r2_clr_dt_mouse, r2_clr_dt_mouse$V1 <0.2)
dim(r2_clr_dt_subset_mouse)

# remove samples which have less samples less than 0.2 R2 # this removes 10 samples
sampletoremove_mouse= r2_clr_dt_subset_mouse$sample

RNA_clr_mouse= RNA_clr_mouse[,c(sampletoremove_mouse) := NULL]
dim(RNA_clr_mouse)

# similarly for ribo

Ribo_clr_mouse= Ribo_clr_mouse[,c(sampletoremove_mouse) := NULL]
dim(Ribo_clr_mouse)

#similarly for TE
TE_GBM_clr_mouse[1:5,1:5]
TE_GBM_clr_mouse= TE_GBM_clr_mouse[,c(sampletoremove_mouse) := NULL]
dim(TE_GBM_clr_mouse)


RNA_mouse_GBM=write.csv(RNA_clr_mouse, "./RNA_clr_mouse.csv")
Ribo_mouse_GBM=write.csv(Ribo_clr_mouse, "./Ribo_clr_mouse.csv")
TE_mouse_GBM=write.csv(TE_GBM_clr_mouse, "./TE_GBM_clr_mouse.csv")


