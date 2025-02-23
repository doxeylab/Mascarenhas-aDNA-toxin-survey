pid = read.csv("percent_no_gaps_matrix.csv",header=T)
cov = read.csv("coverage_matrix.csv",header=T)

removeList = c("SRR15512620","ERR4373999","ERR6096533","ERR6096535","ERR6096536","ERR6096933","ERR6096934","ERR6096935","ERR6473099","ERR6473109","ERR6473106","SRR1313909","SRR1314212","SRR1298752")
toRemove = match(removeList,pid[,1])

pid = pid[-toRemove,]
cov = cov[-toRemove,]


metadata = read.csv("metadata.csv",header=T)

metadata = metadata[match(pid[,1],metadata[,3]),]


rownames(pid) = pid[,1]
rownames(cov) = cov[,1]

rownames(pid) = paste(rownames(pid),metadata$BioSample)
rownames(cov) = paste(rownames(cov),metadata$BioSample)



library(ggplot2)

# Convert matrices/data tables into a long format data frame
df <- data.frame(
  Sample = rep(rownames(pid), times = ncol(pid)),  # Repeating row names as y-axis labels
  Gene = rep(colnames(pid), each = nrow(pid)),  # Repeating column names as x-axis labels
  PID = as.numeric(as.matrix(pid)),  # Convert pid to a vector for coloring
  Coverage = as.numeric(as.matrix(cov))  # Convert cov to a vector for bubble size
)

# Ensure numerical values and remove any NA or infinite values
df <- df[is.finite(df$PID) & is.finite(df$Coverage), ]

# Create the bubble plot
ggplot(df, aes(x = Gene, y = Sample, size = Coverage, color = PID)) +
  geom_point(alpha = 0.8) +  # Transparency for clarity
  scale_size(range = c(1, 10)) +  # Scale bubble sizes
  scale_color_gradientn(colors = c("blue", "red"), limits = c(0.95, 1.0)) +  # Color scale from 0.95 to 1.0
  theme_minimal() +  # Clean theme
  labs(x = "Gene", y = "Sample", size = "Coverage", color = "Percentage Identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels for readability


#removes subtypes to eliminate unnecessary redundancy
toRemove = c("A2","A3","A4","A5","A6","A7","A8","B2","B3","B4","B5","B6","B7","B8","CD","DC","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11","E12","F2","F3","F4","F5","F6","F7","F8","F9")

toKeep = setdiff(colnames(pid),toRemove)

pid = pid[,toKeep]

cov = cov[,toKeep]


#ERR6096937 was blast verified as TeNT not B1, so this was manually changed

pid$B1[which(pid$Dataset == "ERR6096937")] <- 0



# IF a sample has higher %ID AND COVERAGE for C than D, then identify the match as C

for (i in 1:nrow(pid))
{
	if ((pid$C[i] > pid$D[i]) & (cov$C[i] > cov$D[i]))
	{		pid$D[i] = 0
			cov$D[i] = 0
	}

}


# IF a sample has higher %ID AND COVERAGE for D than C, then identify the match as D

for (i in 1:nrow(pid))
{
	if ((pid$D[i] > pid$C[i]) & (cov$D[i] > cov$C[i]))
	{		pid$C[i] = 0
			cov$C[i] = 0
	}

}




library(ggplot2)

# Convert matrices/data tables into a long format data frame
df <- data.frame(
  Sample = rep(rownames(pid), times = ncol(pid)),  # Repeating row names as y-axis labels
  Gene = rep(colnames(pid), each = nrow(pid)),  # Repeating column names as x-axis labels
  PID = as.numeric(as.matrix(pid)),  # Convert pid to a vector for coloring
  Coverage = as.numeric(as.matrix(cov))  # Convert cov to a vector for bubble size
)


df <- df[df$PID > 0, ]


# Ensure numerical values and remove any NA or infinite values
df <- df[is.finite(df$PID) & is.finite(df$Coverage), ]

# Create the bubble plot
ggplot(df, aes(x = Gene, y = Sample, size = Coverage, color = PID)) +
  geom_point(alpha = 0.8) +  # Transparency for clarity
  scale_size(range = c(1, 10)) +  # Scale bubble sizes
  scale_color_gradientn(colors = c("blue", "red"), limits = c(0.96, 1.0)) +  # Color scale from 0.95 to 1.0
  theme_minimal() +  # Clean theme
  labs(x = "Gene", y = "Sample", size = "Coverage", color = "Percentage Identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels for readability


par(mfrow=c(2,1))

boxplot(df$PID ~ df$Gene)
boxplot(df$Coverage ~ df$Gene)

# Boxplot for PID (Percent Identity) per Gene
ggplot(df, aes(x = Gene, y = PID, fill = Gene)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  # Boxplot with transparency, hide outliers
  geom_jitter(alpha = 0.3, width = 0.2) +  # Jittered points for better visualization
  theme_minimal() +
  scale_fill_manual(values = rainbow(length(unique(df$Gene)))) +  # Unique colors per gene
  labs(x = "Gene", y = "Percent Identity (PID)", title = "Distribution of Percent Identity Across Genes") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none")  # Rotate x-axis labels

# Boxplot for Coverage per Gene
ggplot(df, aes(x = Gene, y = Coverage, fill = Gene)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(alpha = 0.3, width = 0.2) +
  theme_minimal() +
  scale_fill_manual(values = rainbow(length(unique(df$Gene)))) +
  labs(x = "Gene", y = "Coverage", title = "Distribution of Coverage Across Genes") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none")  # Rotate x-axis labels

dev.new()


library(dplyr)

# Convert the `pid` matrix to a data frame for processing
pid_df <- as.data.frame(pid)
pid_df$Dataset <- rownames(pid_df)  # Preserve dataset names



# Convert PID values to binary presence/absence (1 if match > 0, 0 otherwise)
binary_pid <- pid_df %>%
  mutate(across(-Dataset, ~ ifelse(. > 0, 1, 0)))  # Apply threshold across all genes

# Generate a column that summarizes each dataset’s gene content
binary_pid$Gene_Content <- apply(binary_pid[-1], 1, function(row) {
  paste(names(which(row == 1)), collapse = ", ")
})




# Create a summary table counting how many datasets have each unique gene or gene combo
gene_summary <- binary_pid %>%
  count(Gene_Content, name = "Dataset_Count") %>%
  arrange(desc(Dataset_Count))

# Show breakdown of datasets with their gene content
print(binary_pid[, c("Dataset", "Gene_Content")])

# Show the summary table with frequency of gene combinations
print(gene_summary)

# Create the bar plot
ggplot(gene_summary, aes(x = reorder(Gene_Content, -Dataset_Count), y = Dataset_Count, fill = Gene_Content)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = rainbow(length(unique(gene_summary$Gene_Content)))) +  # Assign unique colors
  theme_minimal() +
  labs(x = "Gene Combination", y = "Frequency (Number of Datasets)", title = "Gene Combination Frequency") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none")  # Rotate labels



dev.new()



# Convert the `pid` and `cov` matrices into a long format data frame
df <- data.frame(
  Sample = rep(rownames(pid), times = ncol(pid)),  # Sample names
  Gene = rep(colnames(pid), each = nrow(pid)),  # Gene names
  PID = as.numeric(as.matrix(pid)),  # Convert PID to numeric
  Coverage = as.numeric(as.matrix(cov))  # Convert Coverage to numeric
)

# Remove points where PID is 0 to avoid plotting unnecessary points
df <- df[df$PID > 0, ]

# Plot PID vs Coverage with colors mapped to Gene
ggplot(df, aes(x = PID, y = Coverage, color = Gene)) +
  geom_point(alpha = 0.8, size = 3) +  # Adjust transparency and point size
  scale_color_manual(values = rainbow(length(unique(df$Gene)))) +  # Assign unique colors
  theme_minimal() +  # Clean theme
  labs(x = "Percentage Identity (PID)", 
       y = "Coverage", 
       color = "Gene",
       title = "PID vs. Coverage by Gene") +
  theme(legend.position = "right")



# Convert the `pid` matrix into a binary presence/absence format (1 if PID > 0, 0 otherwise)
binary_pid <- as.data.frame(pid > 0)
binary_pid$Dataset <- rownames(pid)  # Store dataset names separately

# Remove the "Dataset" column before counting gene frequencies
gene_frequencies <- colSums(binary_pid[, !names(binary_pid) %in% "Dataset"])  # Sum occurrences per gene

# Convert to data frame for plotting
gene_freq_df <- data.frame(Gene = names(gene_frequencies), Frequency = gene_frequencies)

# Create a bar plot of gene frequency
ggplot(gene_freq_df, aes(x = reorder(Gene, -Frequency), y = Frequency, fill = Gene)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = rainbow(length(unique(gene_freq_df$Gene)))) +  # Assign colors
  theme_minimal() +
  labs(x = "Gene", y = "Frequency (Number of Datasets)", title = "Gene Frequency Across Datasets") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none")  # Rotate labels

#print out matrix
gene_freq_df


# Load necessary libraries
library(ggplot2)
library(maps)
library(dplyr)


# Merge metadata with gene content

# Convert PID values to binary presence/absence (1 if match > 0, 0 otherwise)
binary_pid <- pid_df %>%
  mutate(across(-Dataset, ~ ifelse(. > 0, 1, 0)))  # Apply threshold across all genes

# Generate a column that summarizes each dataset’s gene content
binary_pid$Gene_Content <- apply(binary_pid[-1], 1, function(row) {
  paste(names(which(row == 1)), collapse = ", ")
})


metadata$Gene_Content <- binary_pid$Gene_Content  # Ensure correct order

# Count occurrences of each latitude-longitude pair
metadata <- metadata %>%
  group_by(latitude, longitude) %>%
  mutate(count = n()) %>%
  ungroup()

# Set jittering amount (adjust as needed)
jitter_amount <- 10  # Change this if needed

# Apply jitter only to overlapping points (count > 1)
set.seed(42)  # For reproducibility
metadata <- metadata %>%
  mutate(
    jittered_latitude = ifelse(count > 1, latitude + runif(n(), -jitter_amount, jitter_amount), latitude),
    jittered_longitude = ifelse(count > 1, longitude + runif(n(), -jitter_amount, jitter_amount), longitude)
  )

# Get world map data
world_map <- map_data("world")

# Create base world map
ggplot() +
  # Plot world map
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "lightgray", color = "white") +
  # Draw lines connecting jittered points to their true locations
  geom_segment(data = metadata %>% filter(count > 1), 
               aes(x = longitude, y = latitude, 
                   xend = jittered_longitude, yend = jittered_latitude),
               color = "gray50", linewidth = 0.5, alpha = 0.7) + 
  # Plot points (jittered for overlapping, original for unique)
  geom_point(data = metadata, aes(x = jittered_longitude, y = jittered_latitude, color = Gene_Content), 
             size = 3, alpha = 0.8) +
  # Use distinct colors for each gene combination
  scale_color_manual(values = rainbow(length(unique(metadata$Gene_Content)))) +
  # Labels and theme adjustments
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude", color = "Gene Combination", 
       title = "Global Distribution of Gene Combinations (Jittered Only for Overlapping Points)") +
  theme(legend.position = "right")


#### frequency plot per toxin combo


# Load necessary libraries
library(ggplot2)
library(dplyr)

# Ensure Gene_Content is part of metadata
metadata$Gene_Content <- binary_pid$Gene_Content  # Ensure correct order

# Count occurrences of each gene combination in each study title
gene_title_summary <- metadata %>%
  count(Title, Gene_Content, name = "Frequency")  # Count occurrences

# Create the bubble plot
ggplot(gene_title_summary, aes(x = Title, y = Gene_Content, size = Frequency, color = Gene_Content)) +
  geom_point(alpha = 0.8) +  # Adjust transparency
  scale_size(range = c(2, 10)) +  # Adjust bubble sizes
  scale_color_manual(values = rainbow(length(unique(gene_title_summary$Gene_Content)))) +  # Unique colors
  theme_minimal() +
  labs(x = "Study Title", y = "Gene Combination", size = "Frequency",
       color = "Gene Combination", title = "Frequency of Gene Combinations Across Studies") +
  theme(
	
    axis.text.x = element_text(angle = 90, hjust = 1, size = 3),  # Reduce X-axis text size
    axis.text.y = element_text(size = 3),  # Reduce Y-axis text size
    axis.title.x = element_text(size = 3),  # Reduce X-axis title size
    axis.title.y = element_text(size = 3),  # Reduce Y-axis title size
    legend.text = element_text(size = 3),  # Reduce legend text size
    legend.title = element_text(size = 3),  # Reduce legend title size
    plot.title = element_text(size = 10, face = "bold")  # Reduce title size
	)  


#### frequency plot per toxin (preferred)


# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)  # For splitting gene combinations

# Ensure Gene_Content is part of metadata
metadata$Gene_Content <- binary_pid$Gene_Content  # Ensure correct order

# Split Gene_Content into individual genes
gene_study_data <- metadata %>%
  separate_rows(Gene_Content, sep = ", ") %>%  # Splits combined gene strings into rows
  count(Title, Gene_Content, name = "Frequency")  # Count occurrences of each gene in each study

# Create the bubble plot for individual gene frequencies
ggplot(gene_study_data, aes(x = Title, y = Gene_Content, size = Frequency, color = Gene_Content)) +
  geom_point(alpha = 0.8) +  # Adjust transparency
  scale_size(range = c(2, 10)) +  # Adjust bubble sizes
  scale_color_manual(values = rainbow(length(unique(gene_study_data$Gene_Content)))) +  # Unique colors
  theme_minimal() +
  labs(x = "Study Title", y = "Gene", size = "Frequency",
       color = "Gene", title = "Frequency of Individual Genes Across Studies") +
  theme(
	
    axis.text.x = element_text(angle = 90, hjust = 1, size = 3),  # Reduce X-axis text size
    axis.text.y = element_text(size = 3),  # Reduce Y-axis text size
    axis.title.x = element_text(size = 3),  # Reduce X-axis title size
    axis.title.y = element_text(size = 3),  # Reduce Y-axis title size
    legend.text = element_text(size = 3),  # Reduce legend text size
    legend.title = element_text(size = 3),  # Reduce legend title size
    plot.title = element_text(size = 10, face = "bold")  # Reduce title size
	)  

study_dataset_counts <- metadata %>%
  count(Title, name = "Num_Datasets")

# View the resulting table
print(study_dataset_counts)


#### lolliplot SNP plots


# Load required library
library(ggplot2)
library(patchwork)  # Install with install.packages("patchwork") if needed


# Type B EF028399.1

# Define sequence length and SNP positions
sequence_length <- 3876
snp_positions <- c(3724,844, 1829,2069)  # Modify this as needed

# Create data for the sequence line
sequence_df <- data.frame(x = c(0, sequence_length), y = c(0, 0))

# Create data for SNPs (lollipop annotations)
snp_df <- data.frame(
  x = snp_positions,
  y = rep(1, length(snp_positions))  # Height of lollipops
)

# Create the plot
p1 <- ggplot() +
  # Plot sequence as a horizontal line
  geom_line(data = sequence_df, aes(x, y), linewidth = 1.5, color = "black") +
  # Add SNP lollipop stems
  geom_segment(data = snp_df, aes(x = x, xend = x, y = 0, yend = y), 
               color = "red", linewidth = 1) +
  # Add SNP lollipop circles
  geom_point(data = snp_df, aes(x = x, y = y), color = "red", size = 4) +
  # Labels for SNP positions
  geom_text(data = snp_df, aes(x = x, y = y + 0.1, label = x), 
            vjust = -0.5, size = 3, color = "black") +

  # Customize theme
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) +
  # Remove y-axis and set x-axis range
  scale_y_continuous(limits = c(-0.2, 1.5)) +
  scale_x_continuous(limits = c(0, 4000)) +

  ggtitle("BoNT/B variants identified in ancient DNA")

# Display the plot
#print(p1)



# Type En

# Define sequence length and SNP positions
sequence_length <- 3840
snp_positions <- c(1299,2241,2088,537,2414,3168)  # Modify this as needed

# Create data for the sequence line
sequence_df <- data.frame(x = c(0, sequence_length), y = c(0, 0))

# Create data for SNPs (lollipop annotations)
snp_df <- data.frame(
  x = snp_positions,
  y = rep(1, length(snp_positions))  # Height of lollipops
)

# Create the plot
p2 <- ggplot() +
  # Plot sequence as a horizontal line
  geom_line(data = sequence_df, aes(x, y), linewidth = 1.5, color = "black") +
  # Add SNP lollipop stems
  geom_segment(data = snp_df, aes(x = x, xend = x, y = 0, yend = y), 
               color = "red", linewidth = 1) +
  # Add SNP lollipop circles
  geom_point(data = snp_df, aes(x = x, y = y), color = "red", size = 4) +
  # Labels for SNP positions
  geom_text(data = snp_df, aes(x = x, y = y + 0.1, label = x), 
            vjust = -0.5, size = 3, color = "black") +

  # Customize theme
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) +
  # Remove y-axis and set x-axis range
  scale_y_continuous(limits = c(-0.2, 1.5)) +
  scale_x_continuous(limits = c(0, 4000)) +
  ggtitle("BoNT/En variants identified in ancient DNA")

# Display the plot
#print(p2)



# Type C X53751.1:214-4089

# Define sequence length and SNP positions
sequence_length <- 3876
snp_positions <- c(2137, 2138, 124, 1430, 1366, 3234, 1751, 197, 542, 1989, 3426, 1886, 133, 1648, 1840, 2981, 1022, 752)  # Modify this as needed

# Create data for the sequence line
sequence_df <- data.frame(x = c(0, sequence_length), y = c(0, 0))

# Create data for SNPs (lollipop annotations)
snp_df <- data.frame(
  x = snp_positions,
  y = rep(1, length(snp_positions))  # Height of lollipops
)

# Create the plot
p3 <- ggplot() +
  # Plot sequence as a horizontal line
  geom_line(data = sequence_df, aes(x, y), linewidth = 1.5, color = "black") +
  # Add SNP lollipop stems
  geom_segment(data = snp_df, aes(x = x, xend = x, y = 0, yend = y), 
               color = "red", linewidth = 1) +
  # Add SNP lollipop circles
  geom_point(data = snp_df, aes(x = x, y = y), color = "red", size = 4) +
  # Labels for SNP positions
  geom_text(data = snp_df, aes(x = x, y = y + 0.1, label = x), 
            vjust = -0.5, size = 3, color = "black") +

  # Customize theme
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) +
  # Remove y-axis and set x-axis range
  scale_y_continuous(limits = c(-0.2, 1.5)) +
  scale_x_continuous(limits = c(0, 4000)) +
  ggtitle("BoNT/C variants identified in ancient DNA")

# Display the plot
#print(p3)



# Type D X54254.1:47-3877


# Define sequence length and SNP positions
sequence_length <- 3831
snp_positions <- c(1932, 3038, 1618, 1644, 1678, 1746, 1965, 2070, 2159, 2183, 2195, 2235, 3808, 2679)  # Modify this as needed

# Create data for the sequence line
sequence_df <- data.frame(x = c(0, sequence_length), y = c(0, 0))

# Create data for SNPs (lollipop annotations)
snp_df <- data.frame(
  x = snp_positions,
  y = rep(1, length(snp_positions))  # Height of lollipops
)

# Create the plot
p4 <- ggplot() +
  # Plot sequence as a horizontal line
  geom_line(data = sequence_df, aes(x, y), linewidth = 1.5, color = "black") +
  # Add SNP lollipop stems
  geom_segment(data = snp_df, aes(x = x, xend = x, y = 0, yend = y), 
               color = "red", linewidth = 1) +
  # Add SNP lollipop circles
  geom_point(data = snp_df, aes(x = x, y = y), color = "red", size = 4) +
  # Labels for SNP positions
  geom_text(data = snp_df, aes(x = x, y = y + 0.1, label = x), 
            vjust = -0.5, size = 3, color = "black") +

  # Customize theme
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) +
  # Remove y-axis and set x-axis range
  scale_y_continuous(limits = c(-0.2, 1.5)) +
  scale_x_continuous(limits = c(0, 4000)) +
  ggtitle("BoNT/D variants identified in ancient DNA")

# Display the plot
#print(p4)

combined_plot <- p1 / p2 / p3 / p4  # Stack plots vertically

print(combined_plot)
