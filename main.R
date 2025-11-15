
library(dplyr)         
library(tibble)        
library(ggplot2)      
library(Seurat)       
library(fastSave)       
library(COSG)      

source("~/genmark/function_nnls.R")



wilk <- readRDS("~/genmark/wilk.rds")


wilk <- NormalizeData(wilk, verbose = T)


wilk_query = subset(wilk,Donor == "H2")
wilk_ref = subset(wilk,Donor != "H2")


maker_vec = get_marker_gene(wilk_ref,'true')
# wilk_ref = wilk_ref[maker_vec,]

wilk_list = SplitObject(wilk_ref,split.by = "Donor")

prf_list = list()
for (ss in 1:length(wilk_list)) {
  ref_expr <- as.data.frame(GetAssayData(wilk_list[[ss]], slot = "data"))  # genes x ref_cells
  
  # maker_vec = get_marker_gene(wilk_list[[ss]],'true')
  # 
  ref_expr = ref_expr[maker_vec,]
  
  
  ref_meta <- wilk_list[[ss]]@meta.data
  celltype_col <- "true"  # 你的细胞类型列名
  
  celltypes <- sort(unique(ref_meta[[celltype_col]]))
  
  # 计算每个细胞类型的平均表达（pseudo-bulk）
  ref_mean_mat <- sapply(celltypes, function(ct) {
    cells_ct <- rownames(ref_meta)[ref_meta[[celltype_col]] == ct]
    if (length(cells_ct) == 0) {
      rep(0, nrow(ref_expr))
    } else {
      Matrix::rowMeans(ref_expr[, cells_ct, drop = FALSE])
    }
  })
  
  # ref_mean_mat = second_filter(ref_mean_mat)
  # ref_mean_mat: genes x K（K = 细胞类型数）
  
  query_expr <- GetAssayData(wilk_query, slot = "data")  # genes x query_cells
  
  common_genes <- intersect(rownames(ref_mean_mat), rownames(query_expr))
  ref_mean_use <- ref_mean_mat[common_genes, , drop = FALSE]
  query_use    <- query_expr[common_genes, , drop = FALSE]
  
  # 为 NNLS 准备：我们需要的是 A w ≈ b，其中：
  # A = ref_mean_use (G x K), b = x_i (G x 1)
  
  
  
  prf_list[[ss]] = nnls_prd(ref_mean_use,query_use)
}



prd_df = as.data.frame(prf_list)




# 对数据框的每一行应用函数，得到结果列
most_frequent_celltype <- apply(prd_df, 1, get_mode)

# 对数据框的每一行应用函数，得到结果列
most_frequent <- apply(prd_df, 1, get_mode2)




prd_df$anno1 = most_frequent_celltype
prd_df$anno2 = most_frequent_celltype

prd_df$most_frequent = most_frequent

prd_df[prd_df[,'most_frequent']<0.8*(ncol(prd_df) - 3),'anno2'] ='unknown'
row.names(prd_df) = colnames(wilk_query)



ref_expr <- GetAssayData(wilk_query, slot = "data") 
ref_meta <- wilk_query@meta.data
ref_meta$prd = prd_df$anno2


ref_expr = ref_expr[,row.names(prd_df[prd_df$anno2 != 'unknown',])]
ref_meta = ref_meta[row.names(prd_df[prd_df$anno2 != 'unknown',]),]
celltype_col <- "prd"  # 你的细胞类型列名

celltypes <- sort(unique(ref_meta[[celltype_col]]))

# 计算每个细胞类型的平均表达（pseudo-bulk）
ref_mean_mat <- sapply(celltypes, function(ct) {
  cells_ct <- rownames(ref_meta)[ref_meta[[celltype_col]] == ct]
  if (length(cells_ct) == 0) {
    rep(0, nrow(ref_expr))
  } else {
    Matrix::rowMeans(ref_expr[, cells_ct, drop = FALSE])
  }
})

query_expr <- GetAssayData(wilk_query, slot = "data")  # genes x query_cells

common_genes <- intersect(rownames(ref_mean_mat), rownames(query_expr))
ref_mean_use <- ref_mean_mat[common_genes, , drop = FALSE]
query_use    <- query_expr[common_genes, , drop = FALSE]


pred_celltype = nnls_prd(ref_mean_use,query_use)

names(pred_celltype) = colnames(wilk_query)

prd_df[,'anno2'] = pred_celltype


prd_df_sub = prd_df[prd_df$anno1 == prd_df$anno2,]
wilk_query2 = wilk_query@meta.data[row.names(prd_df_sub),]


table(prd_df$anno1 == wilk_query@meta.data$true)/ncol(wilk_query)

table(prd_df$anno2 == wilk_query@meta.data$true)/ncol(wilk_query)

table(prd_df_sub$anno2 == wilk_query2$true)/nrow(wilk_query2)

