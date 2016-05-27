# mctoolsr

##############
# MISC PLOTS #
##############

#' @title Plot Venn Diagrams
#' @description Plot venn diagram showing OTU overlap between multiple sample 
#'  categories.
#' @param input The input dataset as loaded by \code{load_taxa_table()}.
#' @param category The category (header in map) to group samples by.
#' @param pres_thresh The threshold relative abundance for an OTU to be 
#'  considered present for a given sample category.
#' @concept Plots
#' @examples 
#' fv_filt = filter_data(fruits_veggies, filter_cat = 'Sample_type', 
#'                       keep_vals = c('Spinach', 'Lettuce', 'Mushrooms'))
#' plot_venn_diagram(fv_filt, category = 'Sample_type', pres_thresh = 0.005)
plot_venn_diagram = function(input, category, pres_thresh){
  if (!requireNamespace("VennDiagram", quietly = TRUE)) {
    stop(paste0("'VennDiagram' package (>= 1.6.15) needed for this function ", 
                "to work. Please install it."), call. = FALSE)
  }
  #' @keywords internal
  .venn_3cat = function(data, thresh){
    c1count = sum(data[, 2] >= thresh)
    c2count = sum(data[, 3] >= thresh)
    c3count = sum(data[, 4] >= thresh)
    n12ct = sum(rowSums(data[, c(2, 3)] > thresh) == 2)
    n13ct = sum(rowSums(data[, c(2, 4)] > thresh) == 2)
    n23ct = sum(rowSums(data[, c(3, 4)] > thresh) == 2)
    n123ct = sum(rowSums(data[, c(2, 3, 4)] > thresh) == 3)
    grid::grid.newpage()
    VennDiagram::draw.triple.venn(c1count, c2count, c3count, n12ct, n23ct, 
                                  n13ct, n123ct, category = names(data)[2:4], 
                                  fill = c('blue', 'red', 'orange'))
  }
  
  #' @keywords internal
  .venn_4cat = function(data, thresh){
    c1count = sum(data[, 2] > thresh)
    c2count = sum(data[, 3] > thresh)
    c3count = sum(data[, 4] > thresh)
    c4count = sum(data[, 5] > thresh)
    n12ct = sum(rowSums(data[, c(2, 3)] > thresh) == 2)
    n13ct = sum(rowSums(data[, c(2, 4)] > thresh) == 2)
    n14ct = sum(rowSums(data[, c(2, 5)] > thresh) == 2)
    n23ct = sum(rowSums(data[, c(3, 4)] > thresh) == 2)
    n24ct = sum(rowSums(data[, c(3, 5)] > thresh) == 2)
    n34ct = sum(rowSums(data[, c(4, 5)] > thresh) == 2)
    n123ct = sum(rowSums(data[, c(2, 3, 4)] > thresh) == 3)
    n124ct = sum(rowSums(data[, c(2, 3, 5)] > thresh) == 3)
    n134ct = sum(rowSums(data[, c(2, 4, 5)] > thresh) == 3)
    n234ct = sum(rowSums(data[, c(3, 4, 5)] > thresh) == 3)
    n1234ct = sum(rowSums(data[, c(2, 3, 4, 5)] > thresh) == 4)
    grid::grid.newpage()
    VennDiagram::draw.quad.venn(c1count, c2count, c3count, c4count, n12ct, 
                                n13ct, n14ct, n23ct, n24ct, n34ct, n123ct, 
                                n124ct, n134ct, n234ct, n1234ct, 
                                category = names(data)[2:5], 
                                fill = c('blue', 'red', 'green', 'orange'))
  }
  
  # get relative abundances
  otu_RAs = vegan::decostand(input$data_loaded, method='total', MARGIN=2)
  # get means of otus by metadata category to check for presence
  if(! category %in% names(input$map_loaded)){
    stop('Category header not found in mapping file.')
  }
  otu_RAs_t = as.data.frame(t(otu_RAs))
  otu_RAs_t$cat = input$map_loaded[, category]
  otu_RAs_melted = reshape2::melt(otu_RAs_t, id.vars = 'cat')
  otu_RAs_means = dplyr::summarize_(dplyr::group_by_(otu_RAs_melted, "variable", 
                                                   "cat"), 
                                   mean_RA = ~ mean(value))
  # plot diagram, discount OTUs with relative abundances lower than threshold
  otu_RAs_means_cast = reshape2::dcast(otu_RAs_means, variable ~ cat, 
                                       value.var = 'mean_RA')
  if(ncol(otu_RAs_means_cast) - 1 == 3){
    .venn_3cat(otu_RAs_means_cast, pres_thresh)
  } else if (ncol(otu_RAs_means_cast) - 1 == 4){
    .venn_4cat(otu_RAs_means_cast, pres_thresh)
  } else {
    stop('Can only plot Venn with 3 or 4 category levels.')
  }
}

# interactive plots
# plot_interactive = function(ordination_axes, input, color_by, name_by, out_fp){
#   to_plot = data.frame(ordination_axes, point_colors = input$map_loaded[, color_by], 
#                        point_names = input$map_loaded[, name_by])
#   clickme::clickme('points', to_plot[, 1], to_plot[, 2], 
#                    color_groups = factor(to_plot$point_colors), 
#                    names = to_plot$point_names, 
#                    opacity = 0.8,
#                    file_path = out_fp 
#   )
# }