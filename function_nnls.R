nnls_prd = function(ref_mean_use,query_use){
  library(nnls)
  
  K <- ncol(ref_mean_use)
  nq <- ncol(query_use)
  
  W <- matrix(NA, nrow = nq, ncol = K)  # 存放每个细胞对应的系数
  colnames(W) <- colnames(ref_mean_use)  # 细胞类型
  rownames(W) <- colnames(query_use)     # query 细胞名
  
  A <- as.matrix(ref_mean_use)  # genes x K
  
  for (i in seq_len(nq)) {
    b <- as.numeric(query_use[, i])  # 这个细胞的表达向量
    fit <- nnls(A, b)                # solve min ||A w - b|| s.t. w>=0
    W[i, ] <- fit$x
  }
  
  
  
  # 归一化系数（避免全 0 情况）
  W_norm <- W / rowSums(W + 1e-8)
  
  # 最大系数对应的 cell type
  max_ct_idx <- max.col(W_norm, ties.method = "first")
  pred_celltype <- colnames(W_norm)[max_ct_idx]
  return(pred_celltype)
}




get_mode <- function(row) {
  # 统计每个元素出现的次数
  freq <- table(row)
  # 返回出现次数最多的元素（如果有多个相同频率，返回第一个）
  names(freq)[which.max(freq)]
}

get_mode2 <- function(row) {
  # 统计每个元素出现的次数
  freq <- table(row)
  # 返回出现次数最多的元素（如果有多个相同频率，返回第一个）
  freq[which.max(freq)]
}


get_marker_gene = function(SC_ref,celltype_){

SC_ref@meta.data$celltype_label = SC_ref@meta.data[,celltype_]

cell_num =table(SC_ref@meta.data$celltype_label)
cell_num = cell_num[cell_num>10]
SC_ref = subset(SC_ref,celltype_label %in% names(cell_num))

Idents(SC_ref)<-'celltype_label'
COSG_markers <- cosg(
  SC_ref,
  groups='all',
  assay='RNA',
  slot='data',
  mu=1,
  n_genes_user=100)



genelist<-as.data.frame(COSG_markers[["names"]][,1])
colnames(genelist)<-'dge'
for(gene in 2:ncol(COSG_markers[["names"]])){
  dge2<-as.data.frame(COSG_markers[["names"]][,gene])
  colnames(dge2)<-'dge'
  genelist<-rbind(genelist,dge2)
}
genelist<-unique(genelist)
return(as.vector(genelist[,1]))
}


second_filter = function(ref_df){

secondFC <- c()
for(gene in rownames(ref_df)){
  secondFC <- c(secondFC, sort(as.numeric(ref_df[gene, ]), decreasing = T)[1]/sort(as.numeric(ref_df[gene, ]), decreasing = T)[2])
}
names(secondFC) <- row.names(ref_df)
secondFC=sort(secondFC)

secondFC_in=secondFC[secondFC>=1.5]

return(ref_df[names(secondFC_in),])
}
