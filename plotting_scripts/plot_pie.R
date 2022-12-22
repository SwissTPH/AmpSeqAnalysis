################################
# Script for plotting pie charts with plotly
#
# created 27.10.2022
#
# Contributors:
# monica.golumbeanu@unibas.ch
# sara.cantoreggi@swisstph.ch
###############################

# Install package for export function
# install.packages("webshot")
# webshot::install_phantomjs()

library(readxl)
library(ggplot2)
library(tidyverse)
library(plotly)

# Function for processing the data for plotting
process_plot_data = function(plot_data_df) {
  # Make sure that the colnames are the ones expected later for plotting
  colnames(plot_data_df) = c("Size", "Amount", "Frequency")
  
  # Remove rows with NA size
  plot_data_df = plot_data_df[!is.na(plot_data_df$Size), c("Size", "Amount", "Frequency")]
  
  # Process the data frame, define labels
  plot_data_df$Size = factor(plot_data_df$Size, levels = sort(plot_data_df$Size))
  plot_data_df$Label = paste0(plot_data_df$Amount, " (", round(plot_data_df$Amount/sum(plot_data_df$Amount), digits = 2)*100,"%)")
  
  plot_data_df = plot_data_df[order(plot_data_df$Amount, decreasing = TRUE), ]
  
  # Create sector labels
  plot_data_df$pct = as.double(round(plot_data_df$Amount/sum(plot_data_df$Amount), 2) * 100)
  plot_data_df$pct[which(plot_data_df$pct <= 1)] = 0  # Anything less than 1% should be blank
  plot_data_df$pct = paste0(plot_data_df$pct, "%")
  plot_data_df$pct[which(plot_data_df$pct == "0%")] = ""
  
  return(plot_data_df)
}

# Read the data - the path to the xlsx file needs to be modified by the user
plot_data1 = read_excel("~/genotyping/test_data/FC27_all sites_final data.xlsx")
plot_data1 = process_plot_data(plot_data1)

# Define the color palettes for different microsatellites
colfunc_FC27 <- colorRampPalette(rev(c("#f7fbff","#deebf7", "#9ecae1", "#6baed6", "#2171b5")))
colfunc_3D7 <- colorRampPalette(rev(c("#fff5f0","#fee0d2", "#fc9272", "#fb6a4a", "#ef3b2c")))
colfunc_MAD20 <- colorRampPalette(rev(c("#f7fcf5","#e5f5e0","#a1d99b","#74c476","#238b45")))
colfunc_K1 <- colorRampPalette(rev(c("#efedf5","#dadaeb","#9e9ac8","#807dba","#6a51a3")))
colfunc_Glurp <- colorRampPalette(rev(c("#e5f5f9", "#ccece6","#99d8c9","#66c2a4","#41ae76")))
colfunc_TA81 <- colorRampPalette(rev(c("#f0f0f0", "#bdbdbd", "#969696", "#525252", "#252525")))
colfunc_Polya <- colorRampPalette(rev(c("#fee6ce", "#fdd0a2", "#fdae6b", "#fd8d3c", "#d94801")))
colfunc_PfPK2 <- colorRampPalette (rev(c("#e7e1ef", "#c994c7", "#df65b0", "#e7298a", "#ce1256")))

# Plotting one pie with plotly
a = plot_ly(plot_data1, labels = ~Size, values = ~Amount, type = 'pie',
            textposition = 'outside',
            text = ~pct,
            textinfo = 'text',
            marker = list(colors = colfunc_FC27(nrow(plot_data1)),
                          line = list(color='#FFFFFF', width=0.5))) %>%
  layout(title = "<b>FC27 - All sites</b>",
         xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         uniformtext = list(minsize = 15, mode = 'hide'),
         legend = list(title = list(text = paste0('<b>Genotypes:\n', '(N=', nrow(plot_data1), ')</b>')),
                       orientation = 'h'))
# font = list(size=11))) #add this line for Glurp only (fragment sizes are longer, otherwise legend does not fit)

# Save figure to file - the path to the pdf file needs to be modified by the user
export(a, file = "/scicore/home/pothin/golmon00/genotyping/plots/FC_pie_all.pdf")
