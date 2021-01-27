##--- Function that will calculate the variance of each row
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
                         panel.background = element_blank()) + scale_fill_manual(values = color.mixo(1:51))# +
                #   xlim(xlim[1], xlim[2]) + labs(title = title)
                 
  pRight <- ggplot(data.pd, aes(x = PC2, fill = batch)) + 
                   geom_density(size = density.lwd,alpha = 0.5) +  coord_flip() + ylab('Density') +
                   theme(axis.title.x = element_text(size = rel(0.8)), 
                         axis.title.y = element_blank(), axis.line = element_blank(),
                         axis.text = element_blank(), axis.ticks = element_blank(),
                         panel.background = element_blank()) + scale_fill_manual(values = color.mixo(1:51)) #+
                   #xlim(ylim[1], ylim[2])
                 
  g <- ggplotGrob(pMain + theme(legend.position = 'right', legend.box = 'horizontal',
                                  legend.direction = 'vertical', 
                                  legend.key.height = unit(0.2, 'cm'),
                                  legend.key.width = unit(0.1, 'cm'),
                                  legend.title = element_text(size = rel(legend.title.cex)),
                                  legend.spacing.x = unit(0.1, 'cm'),
                                  legend.spacing.y = unit(0.1, 'cm'),
                                  legend.text = element_text(size = rel(legend.cex))))$grobs
    
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    
  grid.arrange(pTop + ggtitle(title) + theme(legend.position = legend.pos, 
                                             legend.title = element_blank(), 
                                             plot.title = element_text(size = 10, face = "bold")) , 
               legend, 
               pMain + theme(legend.position = 'none'), 
               pRight + theme(legend.position = 'none'), 
               ncol = 2, nrow = 2, widths = c(3, 1), heights = c(1, 3))
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
