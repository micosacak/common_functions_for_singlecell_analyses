get_fold_RectPlot <- function(
        the_df, 
        x_column = "Percentage",
        y_column = "Clusters", 
        f_column = "Samples", 
        the_colors = NULL, 
        x_limits = NULL,
        x_breaks = NULL, 
        fold_threshold = 1.25,
        the_angle = 0,
        vjs = NULL,
        hjs = NULL,
        the_line1_shape = "dashed",
        the_line1_color = "black",
        the_line1_size = 1,
        the_line2_shape = "dashed",
        the_line2_color = "red",
        the_line2_size = 1,
        showLegend = TRUE,           # Keep this for backward compatibility
        legend_position = "top",     # New parameter for legend position
        x_title = NULL,
        y_title = NULL,
        legend_title_size = 16,
        x_title_size = 16,
        y_title_size = 16,
        x_text_size = 16,
        y_text_size = 16,
        legend_text_size = 16,
        x_breaks_by = 1,
        y_breaks_by = 1
        ){

    if(is.null(x_title)) x_title = x_column
    if(is.null(y_title)) y_title = y_column
    
    # rename the columns for internal use
    colnames(the_df)[colnames(the_df) == x_column] = "x_column"
    colnames(the_df)[colnames(the_df) == y_column] = "y_column"
    colnames(the_df)[colnames(the_df) == f_column] = "f_column"
    
    #set limits and breaks
    if(is.null(x_limits)){
        x_min = abs(min(the_df$x_column, na.rm = TRUE))+0.1
        x_max = max(the_df$x_column, na.rm = TRUE)+0.1
        x_limits = c(-x_min, x_max)
    }
                           
    if(is.null(x_breaks)){
        xx_left = abs(x_limits[1]-0.1)
        yy_right = abs(x_limits[2]+0.1)
        dff = xx_left-floor(xx_left)
        if(dff > 0 & dff < 0.5){
            dff = 0.5
        }else{
            dff = 1
        }
        xx_left = floor(xx_left)+dff
        #
        dff = yy_right-floor(yy_right)
        if(dff > 0 & dff < 0.5){
            dff = 0.5
        }else{
            dff = 1
        }
        yy_right = floor(yy_right)+dff
        
        x_breaks = c(-xx_left, seq(ceiling(x_limits[1]), -1, by = y_breaks_by), seq(1, floor(x_limits[2]), by = x_breaks_by),yy_right)
        print(x_breaks)
    }

    the_df$Start <- ifelse(the_df$x_column >= 1, 1, -1)  # Bar starts at 1 or -1
    the_df$End <- the_df$x_column  # Bar ends at the actual value
    if(sum(the_df$x_column == Inf) != 0){
        idx = which(the_df$x_column == Inf)
        print(idx)
        the_df$x_column[idx] = 1
        the_df$Start[idx] = 1
        the_df$End[idx] = 1
        the_df$Deviation[idx] = 0
    }

    if(is.null(the_colors)){
        if(class(the_df$f_column) == "factor"){
            the_levels = levels(the_df$f_column)
        }else{
            the_levels = unique(the_df$f_column)
        }
        the_levels = unique(the_df$f_column)
        print(the_levels)
        the_colors = scPalette(8*length(the_levels))
        print(the_colors)
        the_colors = the_colors[seq(1, length(the_colors),8)]
        the_colors = setNames(the_colors, the_levels)
    }
    
    if(is.null(vjs) | is.null(hjs)){
            if(the_angle != 0){
                vjs = 1
                hjs = 0.5
            }else{
                vjs = 1
                hjs = 1
            }
    }
    
    # Create the plot with geom_rect
    pp <- ggplot(the_df, aes(y = y_column, fill = f_column)) +
        geom_rect(aes(xmin = Start, xmax = End, 
                  ymin = as.numeric(factor(y_column)) - 0.44, 
                  ymax = as.numeric(factor(y_column)) + 0.44), 
              position = position_dodge(width = 0.9)) +
        labs(x = x_title, y = y_title) +
        theme_classic() +
        scale_fill_manual(values = the_colors) +
        scale_x_continuous(breaks = x_breaks, limits = x_limits) +
        theme(
            axis.ticks = element_line(colour = "black"), 
            legend.title = element_text(face = "bold", colour = "black", size = legend_title_size),
            axis.text.x = element_text(angle = the_angle, vjust = vjs, hjust = hjs, face = "bold", size = x_text_size, color = "black"),
            axis.text.y = element_text(face = "bold", colour = "black", size = y_text_size),
            axis.title.x = element_text(face = "bold", size = x_title_size),
            axis.title.y = element_text(face = "bold", size = y_title_size),
            panel.background = element_rect(fill = "white", colour = "white"), 
            plot.background = element_rect(fill = "white", colour = "white"), 
            legend.text = element_text(face = "bold", size = legend_text_size)
        ) +
        geom_vline(xintercept = c(-1, 1), linetype = the_line1_shape, color = the_line1_color, size = the_line1_size) +
        geom_vline(xintercept = c(-fold_threshold, fold_threshold), linetype = the_line2_shape, color = the_line2_color, size = the_line2_size) +
        guides(fill = guide_legend(title = NULL))

    # Handle legend positioning and visibility
    if(!showLegend || legend_position == "none") {
        pp <- pp + theme(legend.position = "none")
    } else {
        # Validate legend_position
        valid_positions <- c("top", "right", "bottom", "left")
        if(!legend_position %in% valid_positions) {
            warning("Invalid legend_position. Using 'top' instead.")
            legend_position <- "top"
        }
        
        pp <- pp + theme(
            legend.position = legend_position,
            legend.justification = "center"
        )
    }
    
    print(x_limits)
    return(pp)
}


#example usage
the_counts = table(the.obj$inj, the.obj$main5)
the_counts = as.data.frame.matrix(the_counts)
the_counts = t(the_counts)
the_counts = sweep(the_counts, 2, colSums(the_counts), "/")*100
the_counts = as.data.frame(the_counts)
the_counts$Fold = (the_counts[["12hpl"]])/(the_counts$Nai)
the_counts$Fold[the_counts$Fold < 1] = -1/the_counts$Fold[the_counts$Fold < 1]
head(the_counts)
the_counts$Clusters = rownames(the_counts)
the_counts$Clusters = factor(the_counts$Clusters, levels = levels(the.obj$main5))


# add Samples for coloring!
the_counts$Samples = the_counts$Clusters # use cluster colors


p0_8 = get_fold_RectPlot(the_counts, x_column = "Fold", y_column = "Clusters", f_column = "Samples", 
            the_colors = the_main_colors, the_angle = 90, vjs = 0.5, hjs = 1, x_title = "Fold Change",
            showLegend = F,legend_position = "right", the_line1_shape = "solid", the_line1_size = 0.5,
            the_line2_shape = "dashed", the_line2_size = 0.5,
            x_title_size = 24,
            x_text_size = 22,
            y_title_size = 24,
            y_text_size = 22, x_breaks_by = 2
            )
p0_8
