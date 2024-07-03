setwd("C:/geo_data/BMP")

library("GenVisR")
library(ggplot2)
library(gplots)

file1=read.csv("ResultBMP.csv")
head(file1)

# Create a color vector for bars
file1$Color <- ifelse(file1$log2FoldChange >= 0, "Upregulated", "Downregulated")
file1

# Create the waterfall plot using ggplot2
waterfall_plot <- ggplot(file1, aes(x = Gene, y =log2FoldChange, fill = Color)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = ifelse(log2FoldChange != 0, sprintf("%.2f", log2FoldChange), "")),
            vjust = ifelse(file1$log2FoldChange >= 0, -0.5, 1),
            size = 3) +
  scale_fill_manual(values = c("Upregulated" = "green", "Downregulated" = "red"),
                    guide = guide_legend(title = "Regulation")) +
  theme_minimal() +
  labs(title = "Waterfall Plot", y = "Fold Change") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(waterfall_plot)
