#Load libraries
library(phyloseq)
library(microbiome)
library(readxl)
library(randomcoloR)



# Upload two Excel files 
# To graph different cohorts, the mouse2 file will change to graph the different cohorts. mouse1 will remain the same
genus_names <- read_excel('/Users/lindsayhopson/Desktop/THESIS/genus_names.xlsx')
abundance <- read_excel('/Users/lindsayhopson/Desktop/THESIS/abundance.xlsx')


OTU = otu_table(abundance, taxa_are_rows = TRUE)
TAX = tax_table(genus_names)
physeq = phyloseq(OTU, TAX)

# Define the order of the x-axis by first getting the values from mouse2.

order_values <- names(abundance)
library(formatR)
#This code revises the origanal plot_bar function so that give non default rainbow colors,and chnage other features for the bar plot to best plot our data 
custom_color_bar_plot <- function (physeq, x_pass = 'Sample', y_pass = 'Abundance', fill_variable = NULL, title = NULL, fac_grid = NULL, order_pass = NULL) 
{
     
     mdf = psmelt(physeq)
     

     if(!is.null(order_pass)) 
     {
         
         # Order the data frame first.
         
         #Turn your 'treatment' column into a character vector
         mdf$Sample <- as.character(mdf$Sample)
         
         #Then turn it back into a factor with the levels in the correct order
         mdf$Sample <- factor(mdf$Sample, levels=order_pass)
  
         
     }
     
     
     colorCount <- length(unique(mdf$ta1))
    
   
      # The number of colors specified has to be the same number of genus listed in mouse1 (the number of columns in mouse1 and mouse2 have to be the same number of specified colors )
      # labs(fill = "Genus") is added to change the legend. It's hard coded that "ta1" represents genus. So this code was added to manually change the legend
     p = ggplot(mdf, aes_string(x = x_pass, y = y_pass, fill = fill_variable)) + geom_bar(stat = 'identity', position = 'stack', color = 'black') + theme(axis.text.x = element_text(angle = -90, hjust = 0)) + scale_fill_manual(values = randomColor(colorCount, hue = "random", luminosity = " ")) + labs(fill = "Genus") + ylab("Relative Abundance") 
     if (!is.null(facet_grid)) 
     {
         
         p <- p + facet_grid(fac_grid)
         
     }
     
     if (!is.null(title)) 
     {
         
         p <- p + ggtitle(title)
         
     }
     
     return(p)
}

custom_color_bar_plot(physeq, x_pass = 'Sample', y_pass = 'Abundance', fill_variable = 'ta1', order_pass = order_values)







