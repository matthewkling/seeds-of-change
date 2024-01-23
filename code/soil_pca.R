
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
      # scale_color_manual(values = c("black", "blue", "red", "darkorange", "magenta3",
      #                               "dodgerblue", "darkgreen", "darkred", "limegreen", "darkorchid4")) +
      scale_color_manual(values = colors3d::distant_colors(length(unique(pc$var)), seed = 2)) +
      guides(color = guide_legend(order = 1), 
             alpha = guide_legend(order = 2), 
             linewidth = guide_legend(order = 2)) +
      coord_flip() +
      theme_minimal() +
      # theme(legend.position = "bottom",
      #       legend.direction = "vertical") +
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
      # theme(#legend.position = "bottom",
      #       #legend.direction = "vertical",
      #       strip.text = element_blank()) +
      labs(fill = "PC quantile\n")


# combined plot
p <- maps + loadings + 
      plot_layout(ncol = 1, heights = c(1.2, 4), guides = "collect") +
      plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
ggsave("figures/soil_pca.png", p, width = 12, height = 8, units = "in")


##### deprecated below #######

# soil PCA ##################

vars <- tibble(var = c("bdod", "cec", "clay", "cfvo", "nitrogen", 
                       "ocd", "phh2o", "sand", "silt", "soc"),
               name = c("Bulk density of the fine earth fraction", "Cation Exchange Capacity of the soil", "Proportion of clay particles (< 0.002 mm) in the fine earth fraction", "Volumetric fraction of coarse fragments (> 2 mm)", "Total nitrogen (N)", 
                        "Organic carbon density", "Soil pH", "Proportion of sand particles (> 0.05 mm) in the fine earth fraction", "Proportion of silt particles (>= 0.002 mm and <= 0.05 mm) in the fine earth fraction", "Soil organic carbon content in the fine earth fraction"),
               units = c("cg/cm3", "mmol(c)/kg", "g/kg", "cm3/dm3", "cg/kg", 
                         "hg/m3", "ph*10", "g/kg", "g/kg", "dg/kg")) %>%
      expand_grid(depth = c("0-5cm", "5-15cm", "15-30cm", "30-60cm", "60-100cm", "100-200cm")) %>%
      mutate(path = paste0(var, "/", var, "_", depth, "_mean"))

template <- rast("assets/climate/ensemble 1981-2010 NA AET.tif")[[1]]
r <- list.files("~/data/soilgrids/v2/1000m_ll_california/", full.names = T) %>%
      rast() %>%
      mask(template)

d <- as.data.frame(r, xy = T) %>% as_tibble() %>% na.omit()
xy <- select(d, x, y)
d <- select(d, -x, -y, -contains("ocs"))

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


# loadings

v <- tibble(comp = paste0("PC", 1:ncol(pca$rotation)),
       p_var = pca$sdev^2 / sum(pca$sdev^2) * 100)

pc <- pca$rotation %>%
      as.data.frame() %>%
      rownames_to_column("var") %>%
      gather(comp, loading, -var) %>%
      as_tibble() %>%
      mutate(var = str_remove(var, "_1000"),
             var = str_remove(var, "_mean")) %>%
      separate(var, c("var", "depth1", "depth2")) %>%
      mutate(depth = as.integer(depth1)) %>%
      filter(comp %in% paste0("PC", 1:5),
             var != "ocs") %>%
      left_join(select(vars, var, name) %>% distinct()) %>%
      left_join(v) %>%
      mutate(comp = paste0(comp, " (", round(p_var), "% var.)"))

loadings <- ggplot(pc, aes(loading, comp, 
               color = name, alpha = depth, linewidth = depth,
               group = paste(var, depth1))) +
      geom_vline(xintercept = 0) +
      geom_path() +
      scale_y_discrete(limits=rev) +
      scale_alpha_continuous(range = c(1, .2), breaks = sort(unique(pc$depth))) +
      scale_linewidth_continuous(range = c(.5, 2), breaks = sort(unique(pc$depth))) +
      scale_color_manual(values = c("black", "blue", "red", "darkorange", "magenta3",
                                    "dodgerblue", "darkgreen", "darkred", "limegreen", "darkorchid4")) +
      guides(color = guide_legend(order = 1), 
             alpha = guide_legend(order = 2), 
             linewidth = guide_legend(order = 2)) +
      theme_minimal() +
      theme(legend.position = "bottom",
            legend.direction = "vertical") +
      labs(color = "Soil variable",
           linewidth = "Depth interval\ntop (cm)",
           alpha = "Depth interval\ntop (cm)",
           x = "PC loading",
           y = NULL)


# maps

f <- bind_cols(xy, pca$x) %>% 
      select(x:PC5) %>%
      gather(comp, value, -x, -y) %>%
      group_by(comp) %>%
      mutate(value = rank(value))

maps <- ggplot(f, aes(x, y, fill = value)) +
      facet_grid(comp ~ .) +
      geom_raster() +
      scale_fill_viridis_c() +
      theme_void() +
      theme(legend.position = "bottom",
            legend.direction = "vertical",
            strip.text = element_blank()) +
      labs(fill = "PC\nrank\n")

# combined plot
library(patchwork)
loadings + maps + plot_layout(ncol = 2, widths = c(4, 1))



cor(d[, grepl("0-5cm", colnames(d))], method = "spearman")^2 %>%
      corrplot::corrplot()
cor(d[, grepl("5-15cm", colnames(d))], method = "spearman")^2 %>%
      corrplot::corrplot()


x <- d[, grepl("5-15cm", colnames(d))]
while(ncol(x) > 5){
      message(ncol(x))
      r2 <- cor(x, method = "spearman")^2
      r2[r2 == 1] <- 0
      omit <- which.max(apply(r2, 1, max))
      x <- x[,setdiff(1:ncol(x), omit)]
}
names(x) # 5 uncorrelated variables to keep
