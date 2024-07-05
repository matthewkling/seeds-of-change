
library(terra)
library(tidyverse)
library(patchwork)





r <- list.files("~/data/soilgrids/v1/250m_California/", full.names = T) %>%
      rast()

d <- as.data.frame(r, xy = T) %>% as_tibble() %>% na.omit()
xy <- select(d, x, y)
d <- select(d, -x, -y, 
            -contains("TEXMHT")) # this one is categorical

# deskew
boxtrans <- function(x){
      x[x==0] <- 1e-10
      bc <- MASS::boxcox(lm(x ~ 1), 
                         data = data.frame(x = x), 
                         lambda = seq(-10, 20, .2), plotit = F)
      lambda <- bc$x[which.max(bc$y)]
      if(lambda != 0) y <- (x^lambda - 1) / lambda
      if(lambda == 0) y <- log(x)
      y
}
for(i in 1:ncol(d)){
      message(i)
      d[[i]] <- boxtrans(d[[i]])
} 



# principal components
pca <- d %>% prcomp(scale. = T)

# variables
vars <- read_csv("~/data/soilgrids/v1/META_GEOTIFF_1B.csv") %>%
      mutate(var = str_remove(FileName, "_250m_ll"),
             var = str_remove(var, "M_sl"),
             var = str_remove(var, "\\.tif")) %>%
      select(var, name = VARIABLE_NAME, depth_label = DEPTH) %>%
      separate(var, c("var", "depth")) %>%
      mutate(depth = as.integer(depth)) %>%
      group_by(depth) %>%
      mutate(depth_label = depth_label[1])



# loadings

v <- tibble(comp = paste0("PC", 1:ncol(pca$rotation)),
            p_var = pca$sdev^2 / sum(pca$sdev^2) * 100)

pc <- pca$rotation %>%
      as.data.frame() %>%
      rownames_to_column("var") %>%
      gather(comp, loading, -var) %>%
      as_tibble() %>%
      mutate(var = str_remove(var, "_250m_ll"),
             var = str_remove(var, "M_sl")) %>%
      separate(var, c("var", "depth")) %>%
      mutate(depth = as.integer(depth)) %>%
      filter(comp %in% paste0("PC", 1:5),
             var != "ocs") %>%
      left_join(vars) %>%
      left_join(v) %>%
      arrange(comp, loading) %>%
      mutate(name = factor(name, levels = rev(unique(name[comp == "PC1"])))) %>%
      mutate(comp = paste0(comp, " (", round(p_var), "% var.)")) %>%
      arrange(depth)

loadings <- ggplot(pc, aes(loading, comp, 
                           color = name, alpha = depth, linewidth = depth,
                           group = paste(var, depth))) +
      geom_vline(xintercept = 0) +
      geom_path() +
      scale_y_discrete(limits=rev) +
      scale_alpha_continuous(range = c(1, .2), breaks = unique(pc$depth), labels = unique(pc$depth_label)) +
      scale_linewidth_continuous(range = c(.5, 2), breaks = unique(pc$depth), labels = unique(pc$depth_label)) +
      scale_color_manual(values = colors3d::distant_colors(length(unique(pc$var)), seed = 2)) +
      guides(color = guide_legend(order = 1), 
             alpha = guide_legend(order = 2), 
             linewidth = guide_legend(order = 2)) +
      coord_flip() +
      theme_minimal() +
      labs(color = "\nSoil variable",
           linewidth = "Depth",
           alpha = "Depth",
           x = "PC loading",
           y = NULL)



# maps

f <- bind_cols(xy, pca$x) %>% 
      select(x:PC5) %>%
      gather(comp, value, -x, -y) %>%
      mutate(comp = factor(comp, levels = rev(unique(comp)))) %>%
      group_by(comp) %>%
      mutate(value = rank(value),
             value = scales::rescale(value))

maps <- ggplot(f, aes(x, y, fill = value)) +
      facet_grid(. ~comp) +
      geom_raster() +
      scale_fill_viridis_c() +
      theme_void() +
      labs(fill = "PC quantile\n")


# combined plot
p <- maps + loadings + 
      plot_layout(ncol = 1, heights = c(1.2, 4), guides = "collect") +
      plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
ggsave("figures/soil_pca.png", p, width = 12, height = 8, units = "in")
