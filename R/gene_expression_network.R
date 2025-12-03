library(gaschYHS)
library(dplyr)
library(corrplot)
library(ggplot2)
library(graphSI)
library(igraph)
library(ggraph)
library(scales)
library(networktools)
library(gridExtra)

# Load the gene expression dataset  
data(gaschYHS)

# group A: genes involved in carbohydrate metabolism
genesA <- c("HXT5","HXK1","GLK1","YDR516C","XKS1","PGM2","PFK26",
            "FBP26","GPM2","GLG1", "GSY2","GLC3","TPS2","TPS1",
            "TSL1","NTH1","ATH1", "GPD1")

# group B: genes involved in cellular redox reactions and defense against reactive oxygen species
genesB <- c("CCP1","TRX2","TTR1","SOD1","HYR1","GPX1",
            "CTT1","YDR453C","YBL064C", "YCL035C", "GTT1", "ECM38")

# group C: genes involved in protein folding
genesC <- c("HSP104","HSP42","HSP78", "HSP26","SSA4","SSA3","SSE2")

# group D: genes involved in protein degradation and vacuolar functions
genesD <- c("UBI4","UBC8","UBC5","APG7","UBP15","VAB31","PRC1","YBR139W","LAP4","PEP4",
            "PAI3","PRB1","PMC1","AUT7","APG1")

# group E: genes involved in DNA damage repair
genesE <- c("MAG1","MMS2")

# group F: genes involved in intracellular signaling
genesF <- c("SRA3","PKA3","PDE1", "IKS1", "YAK1", "TOR1", "MTL1", "PIG2", 
            "SDS22", "SIP2", "RIM11", "NPR1", "PTK2", "KNS1", "KIN82", 
            "YIL113W", "RRD2")

genes <- c(genesA, genesB, genesC, genesD, genesE, genesF)

expr_list <- lapply(genes, function(g) {
  if(length(gaschYHS@assayData[["exprs"]][which(gaschYHS@featureData@data[["GENE"]] == g), ]) == 0) {
    print(paste("Gene", g, "not found in dataset"))
    return(rep(NA, ncol(gaschYHS@assayData[["exprs"]])))
  } else{
    return(gaschYHS@assayData[["exprs"]][which(gaschYHS@featureData@data[["GENE"]] == g), ])
  }
})

# Match genes to their group, assign a color
gene_group <- data.frame(
  gene = genes,
  group = c(
    rep("A", length(genesA)),
    rep("B", length(genesB)),
    rep("C", length(genesC)),
    rep("D", length(genesD)),
    rep("E", length(genesE)),
    rep("F", length(genesF))
  ),
  color = c(
    rep("cornflowerblue", length(genesA)),
    rep("purple1", length(genesB)),
    rep("red2", length(genesC)),
    rep("orange", length(genesD)),
    rep("yellow2", length(genesE)),
    rep("olivedrab3", length(genesF))
  ),
  stringsAsFactors = FALSE
)

# We transpose the data to have genes in columns: we are interested in gene co-expression
df <- as.data.frame(t(data.frame(do.call(rbind, expr_list))))

colnames(df) <- genes
df <- df %>% select(where( ~!all(is.na(.))))

# Remove genes and conditions with many missing values
df <- df %>% select(which(colSums(is.na(.)) < 10))
df <- df %>% filter(rowSums(is.na(.)) == 0)

# Remove the removed genes also from gene_groups
gene_group <- gene_group %>% filter(gene %in% colnames(df))

# Compute the correlation matrix
cor_matrix <- cor(df)

# Visualize the correlation matrix
corrplot(cor_matrix, method = "color", tl.srt = 90,
         tl.col = gene_group$color, tl.cex = 0.8)

# Check normality for the first gene of each group
par(mfrow=c(2,3))
for(g in c(genesA[1], genesB[1], genesC[1], genesD[1], genesE[1], genesF[1])) {
  qqnorm(df[[g]], main = paste("Gene", g))
  qqline(df[[g]], col = "red2")
}
mtext("Normal QQ-plots on original marginals", outer = TRUE, side=3, line=1, cex=2, font=2)

# Non-paranormal transformation on the marginals
df <- huge::huge.npn(as.matrix(df), npn.func = "shrinkage")
df <- as.data.frame(df)

# Check normality again after transformation
# title for the plots
par(mfrow=c(2,3))
for(g in c(genesA[1], genesB[1], genesC[1], genesD[1], genesE[1], genesF[1])) {
  qqnorm(df[[g]], main = paste("Gene", g))
  qqline(df[[g]], col = "red2")
}
mtext("Normal QQ-plots on transformed marginals", outer = TRUE, side=3, line=1, cex=2, font=2)

# Select a graph on the entire dataset using the graphical elastic net
n <- nrow(df)
p <- ncol(df)
selected_graph <- graphSelect(df, penalty="elastic net", data.splitting = T, penalize.diagonal = F)
adjacency_matrix <- selected_graph$adjacency.matrix

# Estimate the graph in the selected model and do inference on all selected edges
estimated_graph <- graphInference(df, selected_graph, to.test="all", nullvalue=0)
W <- estimated_graph$estimated.graph

# Visualize the estimated graph with ggraph
nodes <- data.frame(id = 1:p, name = gene_group$gene, Group = gene_group$group, color = gene_group$color)
edges <- cbind(which(W != 0, arr.ind = TRUE), Weight = rescale(abs(W)[W != 0]))
edges <- as.data.frame(edges) %>% filter(row > col) %>% rename(from = row, to = col)
gene_graph <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)

# Which edges are significant?
sign_pvals <- which(estimated_graph$inference$p_value > 0.05)
i <- estimated_graph$inference$row[sign_pvals]
j <- estimated_graph$inference$col[sign_pvals]

# different color for significant and non significant edges
edges$Significance <- ifelse(paste(edges$from, edges$to) %in% paste(i, j), "Not Significant", "Significant")
gene_graph <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
ggraph(gene_graph, layout = "fr") +
  geom_edge_link(aes(alpha = Weight, color = Significance), width = 1) +
  scale_edge_color_manual(values = c("Significant" = "steelblue", "Not Significant" = "grey70")) +
  geom_node_point(aes(color = Group, size = rescale(degree(gene_graph)))) +
  scale_color_manual(values = setNames(gene_group$color, gene_group$group)) +
  geom_node_text(aes(label = name), vjust = 1.5, size = 3) +
  guides(edge_alpha = "none", size = "none", color = guide_legend(title = "Gene Group", override.aes = list(size = 3))) +
  theme_void() +
  ggtitle("Gene Co-expression Network") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold.italic"))

# Top 10 hub genes by node strength

strength_values <- igraph::strength(gene_graph)
top_hub_genes <- names(sort(strength_values, decreasing = TRUE))[1:10]
hub_table <- data.frame(
  Strength = strength_values[top_hub_genes],
  Group = V(gene_graph)$Group[match(top_hub_genes, V(gene_graph)$name)]
)
colnames(hub_table) <- c("Node Strength", "Group")

# The least connected genes are
least_connected_genes <- names(sort(strength_values, decreasing = FALSE))[1:10]
least_table <- data.frame(
  Strength = strength_values[least_connected_genes],
  Group = V(gene_graph)$Group[match(least_connected_genes, V(gene_graph)$name)]
)
colnames(least_table) <- c("Node Strength", "Group")

# Bridge centrality
bridge_values <- bridge(gene_graph, communities = V(gene_graph)$Group)
bridging_hubs <- names(sort(bridge_values$`Bridge Strength`, decreasing = TRUE))[1:10]
bridge_table <- data.frame(
  Bridge_Strength = bridge_values$`Bridge Strength`[bridging_hubs],
  Group = V(gene_graph)$Group[match(bridging_hubs, V(gene_graph)$name)]
)
colnames(bridge_table) <- c("Bridge Strength", "Group")

most_localized <- names(sort(bridge_values$`Bridge Strength`, decreasing = FALSE))[1:10]
localized_table <- data.frame(
  Bridge_Strength = bridge_values$`Bridge Strength`[most_localized],
  Group = V(gene_graph)$Group[match(most_localized, V(gene_graph)$name)]
)
colnames(localized_table) <- c("Bridge Strength", "Group")
control <- grep(
  "Heat.Shock.000.minutes|YPD.stationary.phase|steady.state|deg\\.growth|dtt.000.min..dtt.2",
  rownames(df)
)

# exclude mutants and wild-type strains
mutants <- grep("msn", "Msn", "yap", "Yap", "7286", rownames(df))
control <- setdiff(control, mutants)
stress <- setdiff(1:nrow(df), control)
stress <- setdiff(stress, mutants)

df_control <- df[control, ]
df_stress <- df[stress, ]

# Compare the co-regulation graphs

selected_graph_c <- graphSelect(df_control, penalty="elastic net", lambda=3*sqrt(log(p)/nrow(df_control)), data.splitting = F, penalize.diagonal = T)
adjacency_matrix_c <- selected_graph_c$adjacency.matrix
estimated_graph_c <- graphInference(df_control, selected_graph_c, to.test="none", nullvalue=0)
W_c <- estimated_graph_c$estimated.graph
selected_graph_s <- graphSelect(df_stress, penalty="elastic net", lambda=3*sqrt(log(p)/nrow(df_stress)), data.splitting = F, penalize.diagonal = T)
adjacency_matrix_s <- selected_graph_s$adjacency.matrix
estimated_graph_s <- graphInference(df_stress, selected_graph_s, to.test="none", nullvalue=0)
W_s <- estimated_graph_s$estimated.graph

nodes <- data.frame(id = 1:p, name = gene_group$gene, Group = gene_group$group, color = gene_group$color)

edges_c <- cbind(which(W_c != 0, arr.ind = TRUE), Weight = rescale(abs(W_c)[W_c != 0]))
edges_c <- as.data.frame(edges_c) %>% filter(row > col) %>% rename(from = row, to = col)
edges_s <- cbind(which(W_s != 0, arr.ind = TRUE), Weight = rescale(abs(W_s)[W_s != 0]))
edges_s <- as.data.frame(edges_s) %>% filter(row > col) %>% rename(from = row, to = col)

gene_graph_c <- graph_from_data_frame(d = edges_c, vertices = nodes, directed = FALSE)
gene_graph_s <- graph_from_data_frame(d = edges_s, vertices = nodes, directed = FALSE)

layout_c <- create_layout(gene_graph_c, layout = "fr")
coords <- layout_c[, c("x","y")]
coords <- as.matrix(coords)

p1 <- ggraph(layout_c) +
  geom_edge_link(aes(alpha = Weight), color = "grey70", width = 1) +
  geom_node_point(aes(color = Group,
                      size = rescale(degree(gene_graph_c)))) +
  scale_color_manual(values = setNames(gene_group$color, gene_group$group)) +
  geom_node_text(aes(label = name), vjust = 1.5, size = 3) +
  theme_void()+
  theme(legend.position = "none")
p2 <- ggraph(gene_graph_s, layout = coords) +
  geom_edge_link(aes(alpha = Weight), color = "grey70", width = 1) +
  geom_node_point(aes(color = Group,
                      size = rescale(degree(gene_graph_s)))) +
  scale_color_manual(values = setNames(gene_group$color, gene_group$group)) +
  geom_node_text(aes(label = name), vjust = 1.5, size = 3) +
  theme_void() +
  guides(edge_alpha = "none", size = "none", color = guide_legend(title = "Gene Group", override.aes = list(size = 3)))

grid.arrange(p1 + ggtitle("Control Conditions") +
               theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold.italic")),
             p2 + ggtitle("Stressed Conditions") +
               theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold.italic")),
             ncol=2)

# Histogram of edge weights in control and stress
edges_c$Condition <- "Control"
edges_s$Condition <- "Stress"
edges_all <- rbind(edges_c, edges_s)
ggplot(edges_all, aes(x = Weight, fill = Condition)) +
  geom_histogram(position = "dodge", bins = 20, alpha = 0.7) +
  scale_fill_manual(values = c("Control" = "skyblue", "Stress" = "salmon")) +
  labs(title = "Edge Weight Distribution",
       x = "Edge Weight",
       y = "Count") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold.italic"))

# Top hub genes in control
strength_values_c <- igraph::strength(gene_graph_c)
top_hub_genes_c <- names(sort(strength_values_c, decreasing = TRUE))[1:10]
hub_table_c <- data.frame(
  Strength = strength_values_c[top_hub_genes_c],
  Group = V(gene_graph_c)$Group[match(top_hub_genes_c, V(gene_graph_c)$name)]
)
colnames(hub_table_c) <- c("Node Strength", "Group")

# Top hub genes in stress
strength_values_s <- igraph::strength(gene_graph_s)
top_hub_genes_s <- names(sort(strength_values_s, decreasing = TRUE))[1:10]
hub_table_s <- data.frame(
  Strength = strength_values_s[top_hub_genes_s],
  Group = V(gene_graph_s)$Group[match(top_hub_genes_s, V(gene_graph_s)$name)]
)
colnames(hub_table_s) <- c("Node Strength", "Group")

# Less connected genes in control
least_connected_genes_c <- names(sort(strength_values_c, decreasing = FALSE))[1:10]
least_table_c <- data.frame(
  Strength = strength_values_c[least_connected_genes_c],
  Group = V(gene_graph_c)$Group[match(least_connected_genes_c, V(gene_graph_c)$name)]
)
colnames(least_table_c) <- c("Node Strength", "Group")

# Less connected genes in stress
least_connected_genes_s <- names(sort(strength_values_s, decreasing = FALSE))[1:10]
least_table_s <- data.frame(
  Strength = strength_values_s[least_connected_genes_s],
  Group = V(gene_graph_s)$Group[match(least_connected_genes_s, V(gene_graph_s)$name)]
)
colnames(least_table_s) <- c("Node Strength", "Group")

