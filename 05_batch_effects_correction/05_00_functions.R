# library(factoextra)
# Library(ggpubr)
# library(mixOmics)
# library(gridExtra)
# library(grid)
# library(tidyverse)
# library(gplots)

#--- Function that will calculate the variance of each row
rowVars <- function(x, na.rm = FALSE, dims = 1, unbiased = TRUE, SumSquares = FALSE, twopass = FALSE) {
  if (SumSquares) return(rowSums(x^2, na.rm, dims))
  N <- rowSums(!is.na(x), FALSE, dims)
  Nm1 <- if (unbiased) N - 1 else N
  if (twopass) {x <- if (dims==0) x - mean(x, na.rm=na.rm) else
    sweep(x, 1:dims, rowMeans(x,na.rm,dims))}
  (rowSums(x^2, na.rm, dims) - rowSums(x, na.rm, dims)^2/N) / Nm1
}

##--- Function for plotting PCA individual map with density plot

PlotPCADensity <- function(data = PCobj, data.pd = Princ.comp, batch = batch, batch.legend.title = 'Batch', 
                           density.lwd = 0.2, legend.pos = 'none',  legend.cex = 0.7, legend.title.cex = 0.75,
                           title = NULL, title.cex = 1.5){
  
  batch <- as.factor(batch)
  
  pca.plate <- fviz_pca_ind(data,
                            col.ind = batch,
                            geom = "point",
                            repel = T)
  pMain <-   ggpar(pca.plate,
                   title = "")
  
  pTop <- ggplot(data.pd, aes(x = PC1, fill = batch)) + 
                   geom_density(size = density.lwd, alpha = 0.5) + ylab('Density') +
                   theme(axis.title.x = element_blank(), axis.title.y = element_text(size = rel(0.8)), 
                         plot.title = element_text(hjust = 0.5, size = rel(title.cex)), legend.position = legend.pos,
                         axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
                         panel.background = element_blank()) + scale_fill_manual(values = color.mixo(1:51)) + labs(colour = batch)
                 
  pRight <- ggplot(data.pd, aes(x = PC2, fill = batch)) + 
                   geom_density(size = density.lwd,alpha = 0.5) +  coord_flip() + ylab('Density') +
                   theme(axis.title.x = element_text(size = rel(0.8)), 
                         axis.title.y = element_blank(), axis.line = element_blank(),
                         axis.text = element_blank(), axis.ticks = element_blank(),
                         panel.background = element_blank()) + scale_fill_manual(values = color.mixo(1:51)) + labs(colour = batch)
                 
  if (batch.legend.title == "Slide"){
    legend <- rectGrob(gp = gpar(fill = NA, col = NA))
    } else{
      g <- ggplotGrob(pTop + theme(legend.position = 'right', legend.box = 'horizontal',
                                  legend.direction = 'vertical', 
                                  legend.key.height = unit(0.2, 'cm'),
                                  legend.key.width = unit(0.1, 'cm'),
                                  legend.title = element_text(size = rel(legend.title.cex)),
                                  legend.spacing.x = unit(0.1, 'cm'),
                                  legend.spacing.y = unit(0.1, 'cm'),
                                  legend.text = element_text(size = rel(legend.cex))))$grobs
      
      legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
      }
    
  grid.arrange(pTop + ggtitle(title) + theme(legend.position = 'none', 
                                     #       legend.title = element_blank(), 
                                            plot.title = element_text(size = 10)) , 
               legend, 
               pMain + theme(legend.position = 'none'), 
               pRight + theme(legend.position = 'none'), 
               ncol = 2, nrow = 2, widths = c(3, 1), heights = c(1, 3))
}

##--- Function that creates report of PCA and ANOVA results 

GetPCAnovaReport <- function(pc.obj, prin.comp, R, pdf.fn){
  
  R <- R + 1

  ##--- ANOVA for detection of variation between PCs and batch
  
  #-- Plate 
  models.plate <- apply(prin.comp[, 2:R], 2, function(pc){
    lm(pc ~ prin.comp$Sample_Plate) })
  anova.plate.tbl <- sapply(models.plate, anova, simplify = F)
  
  #-- ANOVA for Slide
  
  models.slide <- apply(prin.comp[, 2:R], 2, function(pc){
    lm(pc ~ as.factor(as.character(prin.comp$Slide)))})
  anova.slide.tbl <- sapply(models.slide, anova, simplify = F)
  
  #-- for Array
  
  models.array <- apply(prin.comp[, 2:R], 2, function(pc){
    lm(pc ~ prin.comp$Array) })
  anova.array.tbl <- sapply(models.array, anova, simplify = F)
  
  #-- for Sex
  
  models.sex <- apply(prin.comp[, 2:R], 2, function(pc){
    lm(pc ~ prin.comp$sex) })
  anova.sex.tbl <- sapply(models.sex, anova, simplify = F)
  
  #-- for DEX
  
  models.dex <- apply(prin.comp[, 2:R], 2, function(pc){
    lm(pc ~ prin.comp$Sample_Group) })
  anova.dex.tbl <- sapply(models.dex, anova, simplify = F)
  
  #-- Combine p-values into single df to write out to pdf file
  anova.pvalues.df <- t(data.frame(rbind(
    data.frame(anova.plate.tbl)[1, c(5, 10, 15, 20, 25, 30)],
    data.frame(anova.slide.tbl)[1, c(5, 10, 15, 20, 25, 30)],
    data.frame(anova.array.tbl)[1, c(5, 10, 15, 20, 25, 30)],
    data.frame(anova.dex.tbl)[1, c(5, 10, 15, 20, 25, 30)], 
    data.frame(anova.sex.tbl)[1, c(5, 10, 15, 20, 25, 30)])))
  # anova.pvalues.df <- round(anova.pvalues.df, 5)
  colnames(anova.pvalues.df) <- c("Plate", "Slide", "Array", "DEX", "Sex")
  rownames(anova.pvalues.df) <- paste0(rep("PC", R - 1), 1:(R- 1))
  
  #-- Print out plots and table into pdf file
  pdf(pdf.fn)
  
  # PCA and densities plots
  
  title.prefix <- "PCA Individual Map and Density Plots by"
  
  PlotPCADensity(pc.obj, prin.comp, batch = as.character(prin.comp$Sample_Plate), 
                 batch.legend.title = "Plate",
                 title = paste0(title.prefix, "Plate"))
  
  PlotPCADensity(pc.obj, prin.comp, batch = as.character(prin.comp$Slide), 
                 batch.legend.title = "Slide", 
                 title = paste0(title.prefix, "Slide"))
  
  PlotPCADensity(pc.obj, prin.comp, batch = as.character(prin.comp$Array), 
                 batch.legend.title = "Array", 
                 title = paste0(title.prefix, "Array"))
  
  PlotPCADensity(pc.obj, prin.comp, batch = as.character(prin.comp$Sample_Group), 
                 batch.legend.title = "Group", legend.pos = 'right',
                 title = paste0(title.prefix, "Group (dex/veh)"))
  
  PlotPCADensity(pc.obj, prin.comp, batch = as.character(prin.comp$sex), 
                 batch.legend.title = "Sex",  legend.pos = 'none',
                 title = paste0(title.prefix, "Sex"))
  
  # ANOVA result output
  pivot.anova.tbl <- anova.pvalues.df %>% 
    as.data.frame() %>%
    rownames_to_column("PCs") %>%
    pivot_longer(-c(PCs), names_to = "Batch", values_to = "P_value") 
  
  pvalues.heatmap <- ggplot(pivot.anova.tbl, aes(x = Batch, y = PCs, fill = P_value)) + 
    geom_tile(aes(fill = P_value)) +
    geom_text(aes(label = round(P_value, 5))) +
    scale_fill_continuous(low = "red", high = "green") +
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank()) +
    ggtitle("Graphical representation of ANOVA p-values") + theme(plot.title = element_text(size = 10))
  
  print(pvalues.heatmap)
   
  textplot(capture.output(anova.pvalues.df), valign = "top", cex = 0.9)
  title("Summary table of P-values for PCs")
  
  textplot(capture.output(anova.plate.tbl), valign = "top", cex = 0.5)
  title("ANOVA results for Plate")
  
  textplot(capture.output(anova.slide.tbl), valign = "top", cex = 0.5)
  title("ANOVA results for Slide")
  
  textplot(capture.output(anova.array.tbl), valign = "top", cex = 0.5)
  title("ANOVA results for Array")
  
  textplot(capture.output(anova.dex.tbl), valign = "top", cex = 0.5)
  title("ANOVA results for Sample Group")
  
  textplot(capture.output(anova.sex.tbl), valign = "top", cex = 0.5)
  title("ANOVA results for Sex")
  dev.off()
}



##--- Function for plotting PCA individual map
PlotPCAIndMap <- function(PCobj, Princ.comp, pdf.fn){
  #-- by Plate
  pca.plate <- fviz_pca_ind(PCobj,
                            col.ind = as.factor(Prin.comp$Sample_Plate),
                            geom = "point",
                            repel = T)
  ggpar(pca.plate,
        title = "PCA: DEX-Methylation data",
        subtitle = "by Plate",
        legend.title = "Plate", legend = "bottom")
  
  #-- by Slide
  pca.plate <- fviz_pca_ind(PCobj,
                            col.ind = as.factor(as.character(Prin.comp$Slide)),
                            geom = "point",
                            repel = T)
  ggpar(pca.plate,
        title = "PCA: DEX-Methylation data",
        subtitle = "by Slide",
        legend = "none")
  
  #-- by Array
  pca.plate <- fviz_pca_ind(PCobj,
                            col.ind = as.factor(Prin.comp$Array),
                            geom = "point",
                            repel = T)
  ggpar(pca.plate,
        title = "PCA: DEX-Methylation data",
        subtitle = "by Array",
        legend.title = "Array", legend = "bottom")
  
  #-- by Group (dex, veh)
  pca.plate <- fviz_pca_ind(PCobj,
                            col.ind = as.factor(Prin.comp$Sample_Group),
                            geom = "point",
                            repel = T)
  ggpar(pca.plate,
        title = "PCA: DEX-Methylation data",
        subtitle = "by Group",
        legend.title = "Group", legend = "bottom")
  
  #-- by sex (dex, veh)
  pca.plate <- fviz_pca_ind(PCobj,
                            col.ind = as.factor(Prin.comp$sex),
                            geom = "point",
                            repel = T)
  ggpar(pca.plate,
        title = "PCA: DEX-Methylation data",
        subtitle = "by Gender",
        legend.title = "Gender", legend = "bottom")
}
