Network analysis of gene expression data: a case study in *Saccharomyces
cerevisiae*
================

## Table of contents

- [Introduction](#introduction)
- [Exploratory Data Analysis](#exploratory-data-analysis)
- [Gene Co-expression Network Selection, Estimation and
  Inference](#gene-co-expression-network-selection-estimation-and-inference)
  - [Network Analysis: Hub and Bridging
    Genes](#network-analysis-hub-and-bridging-genes)
- [Differential analysis of Gene Co-expression Networks under Control
  and Stressed
  Conditions](#differential-analysis-of-gene-co-expression-networks-under-control-and-stressed-conditions)
  - [Differential network analysis](#differential-network-analysis)
- [References](#references)

## Introduction

In this analysis, we explore the gene expression data from the
`gaschYHS` dataset (Gasch et al. 2000), which contains measurements of
gene expression levels in *Saccharomyces cerevisiae* (baker’s yeast)
under various environmental stress conditions. Our goal is to construct
and analyze gene co-expression networks to identify key genes and their
interactions under both control and stressed conditions. The dataset is
loaded from the `gaschYHS` package available on Bioconductor.

``` r
# Load necessary libraries

# library(Biobase)
# BioCManager::install("gaschYHS")
# BioCManager::install("impute")
library(gaschYHS)
library(dplyr)
library(impute)
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
```

The `gaschYHS` dataset contains gene expression levels for various genes
under different environmental stress conditions. Note that, since each
observation is collected under a different condition, the data points
are not independent and identically distributed (i.i.d.). Our goal is to
pool these conditions and estimate an “average network” of gene
co-expression across all these conditions. In the last section, we will
also compare networks estimated separately on control and stressed
conditions.

We will focus on a subset of genes known to be involved in specific
biological processes. In particular, we will analyze genes involved in
the following functions:

- Carbohydrate metabolism (group A)
- Cellular redox reactions and defense against reactive oxygen species
  (group B)
- Protein folding (group C)
- Protein degradation and vacuolar functions (group D)
- DNA damage repair (group E)
- Intracellular signaling (group F)

These groups are identified in the original publication by Gasch et
al. (2000).

Some of the genes in these groups are excluded as they are not present
in the dataset from Bioconductor; furthermore, we remove genes with more
than missing values on more than 30% of the conditions. We impute the
remaining missing values using k-nearest neighbors imputation.

``` r
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
```

    ## [1] "Gene YDR516C not found in dataset"
    ## [1] "Gene GPX1 not found in dataset"
    ## [1] "Gene YDR453C not found in dataset"
    ## [1] "Gene YBL064C not found in dataset"
    ## [1] "Gene YCL035C not found in dataset"
    ## [1] "Gene YBR139W not found in dataset"
    ## [1] "Gene YIL113W not found in dataset"

``` r
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

# Remove genes with 30% missing values
df <- df %>% select(which(colSums(is.na(.)) < 0.3 * nrow(.)))
# Remove the removed genes also from gene_groups
gene_group <- gene_group %>% filter(gene %in% colnames(df))

# Impute other missing values with knn
df <- as.data.frame(impute.knn(as.matrix(df))$data)
```

In the final dataset, we have 173 observations (conditions) and 64 genes
with complete data for analysis.

## Exploratory Data Analysis

We begin by examining the correlation structure of the selected genes. A
correlation matrix provides insights into the pairwise relationships
between gene expression levels.

``` r
# Compute the correlation matrix
cor_matrix <- cor(df)

# Visualize the correlation matrix
corrplot(cor_matrix, method = "color", tl.srt = 90,
         tl.col = gene_group$color, tl.cex = 0.8)
```

![](README_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

From the correlation matrix, we can observe the following patterns:

- no genes are negatively correlated with each other;
- strong positive correlations exist among genes;
- genes in group A are highly correlated with each other;
- some genes in group B are less correlated with all genes, only with
  each other;
- genes in group C are highly correlated with each other and with group
  A;
- some genes in group F are highly correlated with group A and C and
  with each other;
- in particular, genes TRX2 and VAB31 have low correlation with other
  genes.

In order to apply Gaussian graphical models, we need to check the
normality assumption of the data. We will use QQ plots to assess the
normality of the gene expression levels for a representative gene from
each group. Since the data deviates from normality, we apply a
non-paranormal transformation to the marginals to better satisfy the
normality assumption. Note that to completely satisfy the assumptions of
Gaussian graphical models, we need to further assume that the
multivariate dependence structure is a Gaussian copula.

``` r
# Check normality for the first gene of each group
par(mfrow=c(2,3))
for(g in c(genesA[1], genesB[1], genesC[1], genesD[1], genesE[1], genesF[1])) {
  qqnorm(df[[g]], main = paste("Gene", g))
  qqline(df[[g]], col = "red2")
}
mtext("Normal QQ-plots on original marginals", outer = TRUE, side=3, line=1, cex=2, font=2)
```

![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
# Non-paranormal transformation on the marginals
df <- huge::huge.npn(as.matrix(df), npn.func = "shrinkage")
```

    ## Conducting the nonparanormal (npn) transformation via shrunkun ECDF....done.

``` r
df <- as.data.frame(df)

# Check normality again after transformation
# title for the plots
par(mfrow=c(2,3))
for(g in c(genesA[1], genesB[1], genesC[1], genesD[1], genesE[1], genesF[1])) {
  qqnorm(df[[g]], main = paste("Gene", g))
  qqline(df[[g]], col = "red2")
}
mtext("Normal QQ-plots on transformed marginals", outer = TRUE, side=3, line=1, cex=2, font=2)
```

![](README_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

## Gene Co-expression Network Selection, Estimation and Inference

We now construct the gene co-expression network using the `graphSI`
package. For simplicity, we use data splitting. First, select a sparse
graph using 50% of the data, and an elastic net penalty (Kovács et
al. 2021) to deal with potential multicollinearity among gene expression
levels. On the remaining 50% of the observations, estimate the edge
weights of the selected graph and perform statistical inference to
identify significant edges in the network. This ensures we are not
reusing the same data for selection and inference and thus biasing the
results. The non-significance of some edges means that there is not
enough evidence in the data to reject the hypothesis of the edge being
zero.

``` r
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
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold.italic"))
```

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

### Network Analysis: Hub and Bridging Genes

We will now analyze the network to identify hub genes, which are highly
connected nodes in the network, and bridging genes, which connect to
different gene groups. Hub genes may play crucial roles in maintaining
the structure and function of the gene co-expression network, while
bridging genes may facilitate communication between different biological
processes.

Genes with high node strength (Opsahl et al. 2010) are considered hub
genes, while genes with high bridge centrality (Jones et al. 2021)
connect different gene groups. Many genes in group A (carbohydrate
metabolism) are identified as hub genes, indicating their central role
in the network. On the other hand, genes in group B (cellular redox
reactions) tend to have lower connectivity, suggesting they may function
more independently within the network. Confirming the exploratory
analyses, genes *VAB31*, *TRX2*, *SOD1* have low values on both metrics,
indicating they are less connected and more isolated in the network.

``` r
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
```

<div style="display: inline-block; vertical-align: top; margin-right: 20px;">

<h4>

Top 10 Hub Genes
</h4>

|        | Node Strength | Group |
|:-------|--------------:|:------|
| IKS1   |            44 | F     |
| HXK1   |            38 | A     |
| GSY2   |            38 | A     |
| TSL1   |            38 | A     |
| CTT1   |            38 | B     |
| PGM2   |            37 | A     |
| NTH1   |            36 | A     |
| HSP104 |            35 | C     |
| SSE2   |            35 | C     |
| UBI4   |            35 | D     |

</div>

<div style="display: inline-block; vertical-align: top; margin-right: 20px;">

<h4>

Least Connected Genes
</h4>

|       | Node Strength | Group |
|:------|--------------:|:------|
| VAB31 |             8 | D     |
| SOD1  |             9 | B     |
| TRX2  |            12 | B     |
| NPR1  |            13 | F     |
| GLG1  |            15 | A     |
| UBC5  |            15 | D     |
| MMS2  |            16 | E     |
| CCP1  |            17 | B     |
| TTR1  |            17 | B     |
| HYR1  |            17 | B     |

</div>

<div style="display: inline-block; vertical-align: top; margin-right: 20px;">

<h4>

Top 10 Bridging Hubs
</h4>

|        | Bridge Strength | Group |
|:-------|----------------:|:------|
| CTT1   |              37 | B     |
| IKS1   |              34 | F     |
| HSP42  |              30 | C     |
| HSP104 |              29 | C     |
| SSE2   |              29 | C     |
| UBI4   |              27 | D     |
| SRA3   |              27 | F     |
| MTL1   |              27 | F     |
| GPM2   |              26 | A     |
| NTH1   |              26 | A     |

</div>

<div style="display: inline-block; vertical-align: top;">

<h4>

Most Localized Nodes
</h4>

|       | Bridge Strength | Group |
|:------|----------------:|:------|
| SOD1  |               4 | B     |
| VAB31 |               5 | D     |
| GLG1  |               6 | A     |
| TRX2  |               6 | B     |
| UBC5  |              10 | D     |
| PEP4  |              10 | D     |
| NPR1  |              10 | F     |
| GPD1  |              11 | A     |
| HYR1  |              12 | B     |
| PRC1  |              12 | D     |

</div>

## Differential analysis of Gene Co-expression Networks under Control and Stressed Conditions

Finally, we will compare the gene co-expression networks constructed
from control conditions and stressed conditions. This comparison will
help us understand how environmental stress affects gene interactions
and identify any changes in network structure or key genes.

We exclude mutants (*msn1*, *msn4*, *yap1*) and wild-type strains
(*DBY7286* strain) from the analysis, as they represent genetic
perturbations rather than environmental stress.

As the non-stressed (control) conditions, we consider:

- experiments at 0 minutes;
- stationary phase in YPD medium;
- steady state experiments;
- growth in stationary different temperatures.

As the stressed conditions, we consider all other experiments.

``` r
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
```

We obtain 30 observations for the control group and 143 observations for
the stressed group.

We first compute the correlation matrices for both groups to visualize
differences in pairwise gene relationships. All pairwise correlations
are positive in both groups, the stressed conditions show smaller
correlations overall, indicating that stress may disrupt some gene
co-expression patterns.

``` r
# Correlation matrices for control and stress conditions
cor_matrix_c <- cor(df_control)
cor_matrix_s <- cor(df_stress)
par(mfrow=c(1,2))
corrplot(cor_matrix_c, method="color", tl.srt=90,
         tl.col=gene_group$color, tl.cex=0.8)
mtext("Control Conditions", side=3, line=3)

corrplot(cor_matrix_s, method="color", tl.srt=90,
         tl.col=gene_group$color, tl.cex=0.8)
mtext("Stressed Conditions", side=3, line=3)
```

![](README_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

We compute separate gene co-expression networks for the control and
stressed conditions using the same procedure as before, but without data
splitting since inference is not performed here. We then transform the
estimated precision matrices to partial correlation matrices, to avoid
false discoveries caused by scale differences (Tu et al. 2021).

We compare the two networks by visualizing them side by side, analyzing
the top hub genes in each network, and examining the distribution of
edge weights. Due to the small sample size in the control group, we do
not perform statistical inference here, but only visualize the estimated
graphs. The penalty parameter for graph selection is adjusted according
to the sample size in each group. Due to the low sample size, and the
different sample sizes between the two groups, we should be cautious in
interpreting the results. The stress graph is denser than the control
one, with some positive interactions and some negative interactions in
both cases. This result differs from the exploratory correlation
analysis, where all correlations were positive. This highlights the
importance of considering conditional dependencies rather than just
marginal correlations when analyzing gene co-expression networks:
pairwise correlations represent the global associations between genes,
while partial correlations reflect the correlation between two genes
after conditioning on all other genes, removing the global co-regulation
effect.

``` r
# Compare the co-regulation graphs

selected_graph_c <- graphSelect(df_control, penalty="elastic net", lambda=sqrt(log(p)/nrow(df_control)), data.splitting = F, penalize.diagonal = T)
adjacency_matrix_c <- selected_graph_c$adjacency.matrix
estimated_graph_c <- graphInference(df_control, selected_graph_c, to.test="none", nullvalue=0)
selected_graph_s <- graphSelect(df_stress, penalty="elastic net", lambda=sqrt(log(p)/nrow(df_stress)), data.splitting = F, penalize.diagonal = T)
adjacency_matrix_s <- selected_graph_s$adjacency.matrix
estimated_graph_s <- graphInference(df_stress, selected_graph_s, to.test="none", nullvalue=0)

W_c <- -cov2cor(estimated_graph_c$estimated.graph)
diag(W_c) <- 1
W_s <- -cov2cor(estimated_graph_s$estimated.graph)
diag(W_s) <- 1

nodes <- data.frame(id = 1:p, name = gene_group$gene, Group = gene_group$group, color = gene_group$color)

edges_c <- cbind(which(W_c != 0, arr.ind = TRUE), Weight = rescale(abs(W_c)[W_c != 0]), Sign = sign(W_c)[W_c != 0])
edges_c <- as.data.frame(edges_c) %>% filter(row > col) %>% rename(from = row, to = col)
edges_s <- cbind(which(W_s != 0, arr.ind = TRUE), Weight = rescale(abs(W_s)[W_s != 0]), Sign = sign(W_s)[W_s != 0])
edges_s <- as.data.frame(edges_s) %>% filter(row > col) %>% rename(from = row, to = col)

gene_graph_c <- graph_from_data_frame(d = edges_c, vertices = nodes, directed = FALSE)
gene_graph_s <- graph_from_data_frame(d = edges_s, vertices = nodes, directed = FALSE)

layout_c <- create_layout(gene_graph_c, layout = "fr")
coords <- layout_c[, c("x","y")]
coords <- as.matrix(coords)

p1 <- ggraph(layout_c) +
  geom_edge_link(aes(alpha = Weight, color = factor(Sign)), width = 1) +
  geom_node_point(aes(color = Group,
                      size = rescale(degree(gene_graph_c)))) +
  scale_edge_color_manual(
    name = "Interaction",
    values = c("1" = "skyblue2", "-1" = "salmon2"),
    labels = c("Positive", "Negative")
  ) +
  scale_color_manual(values = setNames(gene_group$color, gene_group$group)) +
  geom_node_text(aes(label = name), vjust = 1.5, size = 3) +
  theme_void()+
  theme(legend.position = "none")
p2 <- ggraph(gene_graph_s, layout = coords) +
  geom_edge_link(aes(alpha = Weight, color = factor(Sign)), width = 1) +
  geom_node_point(aes(color = Group,
                      size = rescale(degree(gene_graph_s)))) +
  scale_edge_color_manual(
    name = "Interaction",
    values = c("1" = "skyblue2", "-1" = "salmon2"),
    labels = c("Positive", "Negative")
  ) +
  scale_color_manual(values = setNames(gene_group$color, gene_group$group)) +
  geom_node_text(aes(label = name), vjust = 1.5, size = 3) +
  theme_void() +
  guides(edge_alpha = "none", size = "none", color = guide_legend(title = "Gene Group", override.aes = list(size = 3)))

grid.arrange(p1 + ggtitle("Control Conditions") +
               theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold.italic")),
             p2 + ggtitle("Stressed Conditions") +
               theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold.italic")),
             ncol=2)
```

![](README_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

We analyze the top hub genes and least connected genes in each network.
In the control conditions, hubs are mostly from group A, suggesting
these genes dominate central positions under normal conditions. Under
environmental stress, group F becomes prominent, possibly reflecting the
activation of stress-response pathways.

``` r
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
```

<div style="display: inline-block; vertical-align: top; margin-right: 20px;">

<h4>

Top 10 Hub Genes in Control
</h4>

|       | Node Strength | Group |
|:------|--------------:|:------|
| KNS1  |            46 | F     |
| GLK1  |            45 | A     |
| HSP42 |            45 | C     |
| PKA3  |            45 | F     |
| PFK26 |            44 | A     |
| ATH1  |            44 | A     |
| IKS1  |            44 | F     |
| MTL1  |            43 | F     |
| TSL1  |            41 | A     |
| HSP26 |            40 | C     |

</div>

<div style="display: inline-block; vertical-align: top;">

<h4>

Top 10 Hub Genes in Stress
</h4>

|       | Node Strength | Group |
|:------|--------------:|:------|
| HSP42 |            36 | C     |
| SSE2  |            33 | C     |
| NTH1  |            32 | A     |
| CTT1  |            32 | B     |
| PKA3  |            32 | F     |
| IKS1  |            32 | F     |
| SRA3  |            31 | F     |
| PFK26 |            30 | A     |
| HSP26 |            28 | C     |
| AUT7  |            27 | D     |

</div>

<div style="display: inline-block; vertical-align: top; margin-right: 20px;">

<h4>

Least Connected Genes in Control
</h4>

|        | Node Strength | Group |
|:-------|--------------:|:------|
| SDS22  |             8 | F     |
| UBI4   |            15 | D     |
| TTR1   |            16 | B     |
| PTK2   |            16 | F     |
| TRX2   |            17 | B     |
| UBP15  |            19 | D     |
| HSP104 |            20 | C     |
| MMS2   |            22 | E     |
| LAP4   |            23 | D     |
| GPD1   |            24 | A     |

</div>

<div style="display: inline-block; vertical-align: top;">

<h4>

Least Connected Genes in Stress
</h4>

|       | Node Strength | Group |
|:------|--------------:|:------|
| TRX2  |            14 | B     |
| CCP1  |            15 | B     |
| SDS22 |            15 | F     |
| VAB31 |            16 | D     |
| SOD1  |            17 | B     |
| MMS2  |            17 | E     |
| TPS2  |            19 | A     |
| HYR1  |            19 | B     |
| SSA4  |            19 | C     |
| PEP4  |            19 | D     |

</div>

### Differential network analysis

Compute now the difference between the two estimated networks to
identify edges that change significantly between control and stressed
conditions. This differential network highlights gene interactions that
are specifically altered in response to environmental stress.

``` r
W_d <- W_s - W_c

# Visualize the differential network
edges_d <- cbind(which(W_d != 0, arr.ind = TRUE), Weight = rescale(abs(W_d)[W_d != 0]), Sign = sign(W_d)[W_d != 0])
edges_d <- as.data.frame(edges_d) %>% filter(row > col) %>% rename(from = row, to = col)
gene_graph_d <- graph_from_data_frame(d = edges_d, vertices = nodes, directed = FALSE)
ggraph(gene_graph_d, layout = coords) +
  geom_edge_link(aes(alpha = Weight, color = factor(Sign)), width = 1) +
  scale_edge_color_manual(
    name = "Interaction",
    values = c("1" = "skyblue2", "-1" = "salmon2"),
    labels = c("strengthened", "weakened")
  ) +
  geom_node_point(aes(color = Group,
                      size = rescale(degree(gene_graph_d)))) +
  scale_color_manual(values = setNames(gene_group$color, gene_group$group)) +
  geom_node_text(aes(label = name), vjust = 1.5, size = 3) +
  theme_void() +
  guides(edge_alpha = "none", size = "none", color = guide_legend(title = "Gene Group", override.aes = list(size = 3))) +
  ggtitle("Differential Gene Co-expression Network (Stress - Control)") +
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold.italic"))
```

![](README_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
strength_values_d <- igraph::strength(gene_graph_d)
top_changed_genes_d <- names(sort(strength_values_d, decreasing = TRUE))[1:10]
changed_table_d <- data.frame(
  Strength = strength_values_d[top_changed_genes_d],
  Group = V(gene_graph_d)$Group[match(top_changed_genes_d, V(gene_graph_d)$name)]
)
colnames(changed_table_d) <- c("Node Strength", "Group")
# Less changed edges in the differential network
least_changed_genes_d <- names(sort(strength_values_d, decreasing = FALSE))[1:10]
least_changed_table_d <- data.frame(
  Strength = strength_values_d[least_changed_genes_d],
  Group = V(gene_graph_d)$Group[match(least_changed_genes_d, V(gene_graph_d)$name)]
)
colnames(least_changed_table_d) <- c("Node Strength", "Group")
```

<div style="display: inline-block; vertical-align: top; margin-right: 20px;">

<h4>

Top 10 Changed Genes in Differential Network
</h4>

|       | Node Strength | Group |
|:------|--------------:|:------|
| GLK1  |            53 | A     |
| HSP42 |            52 | C     |
| PKA3  |            52 | F     |
| IKS1  |            51 | F     |
| KNS1  |            51 | F     |
| MTL1  |            49 | F     |
| ATH1  |            48 | A     |
| PFK26 |            47 | A     |
| NTH1  |            46 | A     |
| CTT1  |            46 | B     |

</div>

<div style="display: inline-block; vertical-align: top;">

<h4>

Least Changed Genes in Differential Network
</h4>

|       | Node Strength | Group |
|:------|--------------:|:------|
| SDS22 |            18 | F     |
| TRX2  |            25 | B     |
| UBP15 |            26 | D     |
| HXK1  |            29 | A     |
| TTR1  |            30 | B     |
| UBI4  |            30 | D     |
| LAP4  |            30 | D     |
| MMS2  |            30 | E     |
| PTK2  |            30 | F     |
| PGM2  |            33 | A     |

</div>

## References

Gasch, Audrey P., et al. “Genomic expression programs in the response of
yeast cells to environmental changes.” *Molecular biology of the cell*
11.12 (2000): 4241-4257.

Jones, Payton J., Ruofan Ma, and Richard J. McNally. “Bridge centrality:
a network approach to understanding comorbidity.” Multivariate
behavioral research 56.2 (2021): 353-367.

Kovács, Solt, et al. “Graphical elastic net and target matrices: Fast
algorithms and software for sparse precision matrix estimation.” *arXiv*
preprint arXiv:2101.02148 (2021).

Opsahl, Tore, Filip Agneessens, and John Skvoretz. “Node centrality in
weighted networks: Generalizing degree and shortest paths.” *Social
networks* 32.3 (2010): 245-251.

Tu, Jia-Juan, et al. “Differential network analysis by simultaneously
considering changes in gene interactions and gene expression.”
*Bioinformatics* 37.23 (2021): 4414-4423.
