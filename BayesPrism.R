library("devtools")
install_github("Danko-Lab/BayesPrism/BayesPrism")
library(BayesPrism)
library(dplyr)
library(pheatmap)


#input file in c drive GBC folder
#load single cell reference data to a dummy object
sc.data <- read.csv("Counts_matrix.csv")
#identify and calculate average of duplicate/isoform genes
sc.data <- sc.data %>%
  group_by_at(1) %>%
  summarise(across(everything(), mean, na.rm = TRUE))
#save to the actual object
sc.dat <- sc.data
sc.dat <- sc.data[, 2:24]
#set rownames for original object using dummy object
rownames(sc.dat) <- sc.data$Gene.Symbol
#transpose the sc object for test compatibility
sc.dat <- as.data.frame(t(sc.dat))

#load bulk seq dataset
bk.dat <- read.delim("OC_Khan_counts.csv")
#set rownames
rownames(bk.dat) <- bk.dat$Geneid
bk.dat <- bk.dat[,7:72]
#transpose for test compatibility
bk.dat <- as.data.frame(t(bk.dat))
metadata <- read.csv("OC_Khan_Metadata.csv")
#specify celltype and cell state lables
cell.type.lables <- rownames(sc.dat)
cell.state.lables <- rownames(sc.dat)

print(cell.type.lables) 
print(cell.state.lables)

#create the prism
myPrism <- new.prism(
  reference = sc.dat,
  mixture = bk.dat,
  input.type = "counts.matrix",
  cell.type.labels = cell.type.lables,
  cell.state.labels = cell.state.lables,
  key = NULL
  )

#quality check for celltypes
result <- myPrism@phi_cellType

#run bayesprism
bp.res <- run.prism(myPrism)
bp.res
#Save the output as a R dataset
saveRDS(bp.res, "prism_result.rds")
#load the results, avoid running the steps again and again
bp.res<- readRDS("prism_result.rds")

#view types of results available
slotNames(bp.res)

#to get posterior mean of factions of cell types in each sample(Actual result)
theta <- get.fraction (bp=bp.res,
                       which.theta = "final",
                       state.or.type = "type")
theta <- as.data.frame(theta)

#make a heatmap to visualize celltype proportions
pheatmap(theta)

#load metadata
metadata <- read.csv("OC_Khan_Metadata.csv")
#add the groups row to the result object to group the samples
theta$Group <- metadata$Sample

# Calculate the average cell proportions for each group
average_by_group <- aggregate(. ~ Group, data = theta, FUN = mean)
#take out group name to name the rows of resulting average
groups <- average_by_group$Group
# Remove the Group column to keep only the averaged cell proportions
average <- (average_by_group[, -1])
rownames(average) <- groups

#heatmap to compare celltype proportions in different groups
pheatmap(average, cluster_rows = FALSE, cluster_cols = FALSE, filename = "Heatmap_BayesPrism_SR.png")


theta$Sample <- rownames(theta)
average$Group <- rownames(average)


theta_long <- gather(theta, key = "Cell_Type", value = "theta", -Sample, -Group)
average_long <- average %>%
  pivot_longer(cols = -Group, names_to = "Cell_Type", values_to = "Mean_theta")




All_plots.bp <- ggplot() +
  geom_boxplot(data = theta_long, aes(x = Group, y = theta, fill = Group)) +
  facet_wrap(~ Cell_Type, scales = "free_y") +
  scale_fill_manual(values = c("Control" = "yellow", "Dysplasia" = "darkgreen", "HkNR" = "orange", "OSCC" = "red")) + # Custom colors
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),   # Remove major grid lines
    panel.grid.minor = element_blank(),   # Remove minor grid lines
    panel.border = element_blank(),       # Remove panel border
    plot.background = element_rect(fill = "white"),  # Set plot background to white
    panel.background = element_rect(fill = "white"), # Set panel background to white
    axis.line.x.bottom = element_line(color = "black", size = 0.5), # Add x-axis line at the bottom
    axis.line.y.left  = element_line(color = "black", size = 0.5), # Add y-axis line on the left
    axis.ticks = element_line(color = "black", size = 0.5), # Add major ticks
    axis.ticks.length = unit(0.3, "cm"), # Adjust length of major ticks
  ) +
  labs(title = "Boxplots and Means of Cell Types by disease state", y = "Cell proportion", x = "Disease state")

All_plots.bp

#to get coefficient of variance for cell fractions in each samples
theta.cv <- bp.res@posterior.theta_f@theta.cv

#identify marker genes using wilcoxon test

theta <- as.data.frame(theta)

celltypes <- colnames(theta)


expression_list <- list()

# Loop through each cell type and store the results in the list
for (i in celltypes) {
  expression_list[[i]] <- get.exp(bp = bp.res,
                          state.or.type = "type",
                          cell.name = i)
}

expression_list <- lapply(expression_list, function(x) {
  as.data.frame(x, stringsAsFactors = FALSE)
})

#Add a row to each data frame indicating the groups
expression_list <- lapply(expression_list, function(df) {
  df$Group <- metadata$Sample[match(rownames(df), metadata$Expt)]
  
  return(df)
})

#split all dataframes in the expression list dataset to evaluate individually
list2env(expression_list, envir = .GlobalEnv)


#pefrom wilcoxon test on each celltype
Adipocytes_result <- pairwise_wilcox_test(
  data = Adipocytes,
  formula = A1BG ~ Group,  # Dependent variable ~ Grouping variable
  p.adjust.method = "holm"
)

expression_list_filter <- lapply(expression_list, function(df) {
  df_filtered <- df[rowSums(df != 0) > 0, ]
  return(df)
})


