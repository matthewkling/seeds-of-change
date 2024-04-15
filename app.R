
library(shiny)
library(shinydashboard)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(leaflet)
library(raster)
library(terra)
library(grid)
library(purrr)
library(stringr)
select <- dplyr::select
extract <- terra::extract



#### data loading ##############################

assets <- ifelse(dir.exists("assets/all_species"), "assets/all_species/", "assets/select_species/")

# climate metadata
vars <- sort(c("PPT", "AET", "CWD", "DJF", "JJA"))
clim_files <- data.frame(path = list.files("assets/climate", full.names = T),
                         stringsAsFactors = F) %>%
      mutate(set = gsub("\\.tif", "", basename(path))) %>%
      separate(set, c("junk", "year", "ssp", "var"), remove = F, sep=" ") %>%
      select(-junk) %>%
      mutate(ssp = ifelse(ssp == "NA", "historic", toupper(ssp))) %>%
      filter(var %in% vars) %>% arrange(var)

times <- unique(clim_files$year)
ssps <- unique(clim_files$ssp)
ssps <- data.frame(ssp = ssps,
                   text = ifelse(ssps == "historic", ssps,
                                 paste0(str_sub(ssps, 1, 4), "-", str_sub(ssps, 5, 5), ".",
                                        str_sub(ssps, 6, 6))))

spps <- list.files(paste0(assets, "ranges")) %>% sub("\\.tif", "", .)
allspps <- readRDS("assets/all_species.rds")

# spps[spps == "NONE"] <- "All of California (no focal species)"
# allspps[allspps == "NONE"] <- "All of California (no focal species)"


soil_all <- rast("assets/soil_800m.tif")
names(soil_all) <- paste0("soil_PC", 1:nlyr(soil_all))

smoothed <- data.frame(path=list.files(paste0(assets, "smooth"), full.names=T),
                       stringsAsFactors=F) %>%
      mutate(set = gsub("\\.tif", "", basename(path)),
             set = str_replace(set, "NONE ", "NONE NONE ")) %>%
      separate(set, c("genus", "species", "radius"), remove=F, sep=" ") %>%
      mutate(gs = paste(genus, species),
             gs = str_replace(gs, "NONE NONE", "NONE"),
             radius = as.integer(radius)) %>%
      select(path, gs, radius)

range_summaries <- readRDS("assets/range_stats.rds")

ll <- crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

all_vars <- data.frame(abbv = c(vars, paste0("soil_PC", 1:5),
                                "prob", "clust", "sigma", paste0("d", vars)))
all_vars$desc <- c("Actual evapotranspiration", "Climatic water deficit",
                   "Winter minimum temperature", "Summer maximum temperature",
                   "Annual precipitation",
                   "Soil PC1", "Soil PC2", "Soil PC3", "Soil PC4", "Soil PC5",
                   "Environmental difference",
                   "Seed zones",
                   "Multivariate change in climate",
                   paste0("Change in ", c("actual evapotranspiration", "climatic water deficit",
                                          "winter minimum temperature", "summer maximum temperature", "annual precipitation")))
all_vars$units <- c("(mm)", "(mm)",
                    "(deg. C)", "(deg. C)",
                    "(log10 mm)",
                    rep("(standard deviations of statewide PCA)", 5),
                    rep("(dissimilarity from focal site, species range standard deviations)", 1),
                    "(environmental clusters)",
                    "(change from local baseline, species range standard deviations)",
                    "(mm)", "(mm)", "(deg. C)", "(deg. C)", "(log10 mm)")
all_vars$min <- c(rep(NA, 10), rep(0, 1), NA, 0, rep(NA, 5))
all_vars$max <- c(rep(NA, 10), rep(5, 1), NA, NA, rep(NA, 5))
all_vars$palette <- c(rep("viridis", 10), rep("sigma", 1), "cluster", rep("sigma", 1), rep("delta", 5))

sigma <- function(x, site_mean = 0, w = rep(1, nlyr(x))){
      w <- w / sum(w)
      k <- nlyr(x)
      z <- sum((x - site_mean)^2 * k * w)
      z[] <- sqrt(qchisq(pchisq(z[], k), 1))
      z
}

# default focal site location: presidio
lonlat <- "-122.52, 37.83"
s <- data.frame(species=NA, lat=38.02, lon=-122.83)
coordinates(s) <- c("lon", "lat")
crs(s) <- ll


# define a keyPressed input that activates when RETURN is pressed
js <- '
$(document).on("keyup", function(e) {
  if(e.keyCode == 13){
    Shiny.onInputChange("keyPressed", Math.random());
  }
});
'


stackBands <- function (paths, band){
      for (i in paths) {
            r <- brick(i)
            r <- raster::subset(r, band)
            if (i == paths[1]) 
                  s <- r
            else {
                  s <- stack(s, r)
            }
      }
      s
}



# UI #################
ui <- dashboardPage(
      skin = "black",
      
      dashboardHeader(title = "Seeds of Change"),
      
      dashboardSidebar(
            sidebarMenu(
                  menuItem("Analysis tool", tabName = "tool", icon = icon("map-location")),
                  menuItem("User guide", tabName = "instructions", icon = icon("circle-info")),
                  menuItem("About", tabName = "about", icon = icon("circle-user"))
            ),
            br(),
            br(),
            p(style = "margin-left: 15px;margin-right: 15px;", 
              "Seeds of Change is a decision support tool for ecological restoration efforts in California. It matches potentially suitable seed collection and planting sites by analyzing spatial data on climate change, soil characteristics, species distributions, and adaptive neighborhoods.",
              br(),
              br(),
              "University of California, Berkeley")
      ),
      
      dashboardBody(
            
            tags$head(
                  tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),
                  includeHTML("google-analytics.html")
            ),
            
            tabItems(
                  tabItem(tabName = "tool",
                          fluidRow(
                                column(
                                      width = 4,
                                      tabBox(
                                            title = "Settings", width = NULL, side = "right", selected = "Context",
                                            tabPanel("Species", 
                                                     selectizeInput("sp", 
                                                                    span("Species ", actionLink("i_species", "[?]")), 
                                                                    choices = spps, selected = "Quercus agrifolia"),
                                                     selectInput("constrain", span("Limit results to species range", actionLink("i_constrain", "[?]")), c(TRUE)),
                                                     shinyWidgets::sliderTextInput("threshold",
                                                                                   span("SDM threshold", actionLink("i_threshold", "[?]")),
                                                                                   choices = c(1, 2, 5, seq(10, 80, 10)),
                                                                                   selected = 50, post = "%", grid = T),
                                                     shinyWidgets::sliderTextInput("radius",
                                                                                   span("Smoothing radius", actionLink("i_radius", "[?]")),
                                                                                   choices=sort(unique(smoothed$radius)), selected = 2,
                                                                                   post = " km", grid = T),
                                                     sliderInput("pclim",
                                                                 span("Soil versus Climate", actionLink("i_weight", "[?]")),
                                                                 min=0, max=100, value=50, step=10, post = "% clim",
                                                                 width="100%")
                                            ),
                                            tabPanel("Context", 
                                                     selectInput("mode", 
                                                                 span("Focal site activity ", actionLink("i_mode", "[?]")), 
                                                                 c("planting", "seed collection"), "planting"),
                                                     textInput("lonlat",
                                                               span("Focal site location ", actionLink("i_location", "[?]")), 
                                                               paste(s$lon, s$lat, sep = ", ")),
                                                     selectInput("time",
                                                                 span("Time period ", actionLink("i_time", "[?]")),
                                                                 times, "2071-2100"),
                                                     selectInput("ssp",
                                                                 span("Scenario ", actionLink("i_ssp", "[?]")),
                                                                 ssps$text, "SSP5-8.5"),
                                                     # uiOutput("constraintControls"),
                                                     selectizeInput("color", 
                                                                    span("Display variable ", actionLink("i_color", "[?]")),
                                                                    all_vars$desc[c(11, 12, 1:10, 13:18)],
                                                                    all_vars$desc[11]),
                                                     uiOutput("nclust")
                                            )
                                      )
                                ),
                                column(
                                      width = 8,
                                      tabBox(
                                            title = "Results", width = NULL, side = "right", selected = "Map",
                                            tabPanel("Export", 
                                                     fluidRow(
                                                           column(12, downloadButton("download", ""), h5("Download map as raster ", actionLink("i_download", "[?]")))
                                                     )
                                            ),
                                            tabPanel("Plot", 
                                                     fluidRow(
                                                           column(12,
                                                                  span(textOutput("color_label1", inline = T), actionLink("i_variable1", "[?]")),
                                                                  plotOutput("legend2", height = "45px")
                                                           )
                                                     ),
                                                     fluidRow(
                                                           column(8, 
                                                                  plotOutput("scatter")
                                                           ),
                                                           column(4, 
                                                                  span(textOutput("target_label", inline = T), actionLink("i_arrows", "[?]"), style="color:red"),
                                                                  br(), br(),
                                                                  textOutput("points_label"),
                                                                  br(), br(),
                                                                  selectizeInput("xvar", "X variable", all_vars$desc[!str_detect(all_vars$desc, "Seed zones|Change in")], all_vars$desc[1]),
                                                                  selectizeInput("yvar", "Y variable", all_vars$desc[!str_detect(all_vars$desc, "Seed zones|Change in")], all_vars$desc[2])
                                                           )
                                                     )
                                            ),
                                            tabPanel("Map", 
                                                     span(textOutput("color_label2", inline = T), actionLink("i_variable2", "[?]")),
                                                     # span(textOutput("color_label2", inline = T), uiOutput("ivar2")),
                                                     plotOutput("legend1", height = "45px"),
                                                     leafletOutput("map", height = "calc(100vh - 210px)")
                                            )
                                      )
                                      
                                )
                          )
                  ),
                  tabItem(tabName = "instructions",
                          fluidRow(
                                column(width = 12,
                                       box(width = NULL,
                                           htmltools::includeMarkdown("assets/instructions.md")
                                       )
                                )
                          )
                  ),
                  tabItem(tabName = "about",
                          fluidRow(
                                column(width = 12,
                                       box(width = NULL,
                                           htmltools::includeMarkdown("assets/about.md")
                                       )
                                )
                          )
                  )
            )
      )
)


# server #################
server <- function(input, output, session) {
      
      
      # informational popups #############################
      
      ibox <- function(id, title, content, show = T){
            b <- modalDialog(title = title, HTML(content),
                             easyClose = TRUE, footer = modalButton("Dismiss"))
            if(show){
                  observeEvent(input[[id]], {showModal(b)} )
            }else{
                  return(b)
            }
      }
      
      modals <- read_csv("assets/modals.csv")
      
      # static elements
      pmap(modals[modals$id != "i_variable",], ibox)
      
      # dynamic popups -- modified approach to avoid triggering the modal when input$color is changed
      var_info <- reactive({ paste0(input$color, " is depicted by color in the figure. ", modals$content[modals$title == input$color]) })
      observeEvent(input$i_variable1, {showModal(ibox("i_variable1", input$color, var_info(), show = F))} )
      observeEvent(input$i_variable2, {showModal(ibox("i_variable2", input$color, var_info(), show = F))} )
      
      
      
      
      # dynamic inputs ######################
      
      output$map_title <- renderText({ paste0(input$sp, ": ",  input$color) })
      
      site <- reactiveValues(point = s)
      
      # update focal site when map is clicked
      observeEvent(input$map_click, {
            s <- data.frame(species=NA, lat=input$map_click$lat, lon=input$map_click$lng)
            updateTextInput(session, "lonlat", value = paste(round(s$lon, 2), round(s$lat, 2), sep = ", "))
            coordinates(s) <- c("lon", "lat")
            crs(s) <- ll
            site$point <- s
      })
      
      # update focal site via input text
      observeEvent(input$keyPressed, {
            crds <- as.numeric(stringr::str_split(input$lonlat, ", ")[[1]])
            s <- data.frame(species=NA, lat=crds[2], lon=crds[1])
            coordinates(s) <- c("lon", "lat")
            crs(s) <- ll
            site$point <- s
      })
      
      # show range constraint input only in collection mode
      observeEvent(input$mode, {
            if(input$mode == "planting") updateSelectInput(session, "constrain", choices = c(TRUE), selected = TRUE)
            if(input$mode == "seed collection") updateSelectInput(session, "constrain", choices = c(TRUE, FALSE))
      })
      
      # show k input only when seed zone output is selected
      output$nclust <- renderUI({
            if(input$color == "Seed zones") sliderInput("k", span("Number of seed zones", actionLink("i_nclust", "[?]")), 2, 8, 5, 1)
      })
      
      ### show only the input choices relevant to selected timeframe ###
      sel <- reactiveValues(ssp = "SSP585", color = all_vars$desc[11]) # container to record selection history
      observeEvent(input$ssp, sel$ssp <- c(sel$ssp, ssps$ssp[ssps$text == input$ssp]))
      observeEvent(input$color, sel$color <- c(sel$color, input$color))
      observeEvent(input$time, {
            freezeReactiveValue(input, "ssp")
            if(input$time == "1981-2010"){
                  updateSelectInput(session, "ssp", choices = "historic", selected = "historic")
                  updateSelectizeInput(session, "color", 
                                       choices = all_vars$desc[c(11, 12, 1:10)], 
                                       selected = tail(sel$color[!str_detect(sel$color, "Change in")], 1))
            }else{
                  updateSelectInput(session, "ssp", choices = setdiff(ssps$text, "historic"), selected = ssps$text[ssps$ssp == tail(sel$ssp[sel$ssp != "historic"], 1)])
                  updateSelectizeInput(session, "color", 
                                       choices = all_vars$desc[c(11, 12, 1:10, 13:18)], 
                                       selected = tail(sel$color, 1))
            }
      })
      
      
      
      # data processing ############################
      
      range_mask <- reactive({
            threshold <- input$threshold / 100
            p <- rast(list.files(paste0(assets, "/ranges"), pattern = input$sp, full.names = T))
            p[p < threshold] <- NA
            y <- smoothed %>% 
                  filter(gs==input$sp,
                         radius==input$radius) %>%
                  pull(path) %>%
                  terra::rast()
            crop(p, crop(y[[1]], p))
      })
      
      clip <- function(x, y) x %>% crop(y) %>% mask(y) %>% trim()
      
      # historic
      smoothed_envt <- reactive({
            req(input$sp, range_mask())
            y <- smoothed %>% 
                  filter(gs==input$sp,
                         radius==input$radius) %>%
                  pull(path) %>%
                  terra::rast() %>% 
                  clip(range_mask())
            names(y) <- sub("800m_", "PC", names(y))
            return(list(clim = y[[6:10]],
                        soil = y[[1:5]]))
      })
      
      ref_envt <- list(clim = stackBands(clim_files$path[clim_files$year=="1981-2010"], 1) %>%
                             rast() %>%
                             setNames(vars))
      
      # future 
      future_envt <- reactive({
            time <- input$time
            ssp <- ssps$ssp[ssps$text == input$ssp]
            if(time == times[1]) ssp <- ssps$ssp[1]
            if(time != times[1] & ssp == "historic") ssp <- tail(sel$ssp[sel$ssp != "historic"], 1)
            
            cf <- clim_files$path[clim_files$year==time & clim_files$ssp==ssp]
            
            clim <- stackBands(cf, 1) %>% rast() %>% setNames(vars)
            
            if(time == "1981-2010"){
                  m1 <- m2 <- m3 <- m4 <- m5 <- clim
            }else{
                  m1 <- stackBands(cf, 2) %>% rast() %>% setNames(vars)
                  m2 <- stackBands(cf, 3) %>% rast() %>% setNames(vars)
                  m3 <- stackBands(cf, 4) %>% rast() %>% setNames(vars)
                  m4 <- stackBands(cf, 5) %>% rast() %>% setNames(vars)
                  m5 <- stackBands(cf, 6) %>% rast() %>% setNames(vars)
            }
            list(clim = clim,
                 m1 = m1, m2 = m2, m3 = m3, m4 = m4, m5 = m5,
                 soil = soil_all)
      })
      
      deltas <- reactive({
            time <- input$time
            ssp <- ssps$ssp[ssps$text == input$ssp]
            var <- all_vars$abbv[all_vars$desc == input$color]
            if(! str_detect(tolower(input$color), "change in")) return(NULL)
            
            msk <- smoothed_envt()[[1]][[1]]
            files <- list.files("assets/deltas", full.names = T)
            
            if(var == "sigma"){ # compute local sigma from multivariate deltas:
                  files <- files[grepl(time, files) & 
                                       grepl(tolower(ssp), files) &
                                       grepl(paste0(names(smoothed_envt()[[1]]), collapse = "|"), files)]
                  delta <- rast(files) %>% crop(msk) %>% mask(msk) %>%
                        "/"(range_stats()$sd[1:5]) %>%
                        sigma() %>%
                        setNames(var) %>% clip(crop(range_mask(), .))
            }else{ # return univariate delta: 
                  files <- files[grepl(time, files) & 
                                       grepl(tolower(ssp), files) &
                                       grepl(str_remove(var, "d"), files)]
                  delta <- rast(files) %>% crop(msk) %>% mask(msk) %>% setNames(var) %>% clip(crop(range_mask(), .))
            }
            return(delta)
      })
      
      # reference climate across range (historic if planting mode)
      range_envt <- reactive({
            y <- switch(input$mode,
                        "planting" = smoothed_envt(),
                        "seed collection" = future_envt() )
            if(input$mode == "seed collection" & input$constrain == TRUE) y <- y %>% map(clip, y = range_mask())
            return(y)
      })
      
      
      # target environment at selected site (future if planting mode)
      target_envt <- reactive({
            future <- future_envt() %>% map(extract, y = coordinates(site$point)) %>% map(as.vector)
            historic <- smoothed_envt() %>% map(extract, y = coordinates(site$point)) %>% map(as.vector)
            reference <- ref_envt %>% map(extract, y = coordinates(site$point)) %>% map(as.vector)
            te <- list(focal = switch(input$mode, "planting" = future, "seed collection" = historic),
                       reference = switch(input$mode, "planting" = reference, "seed collection" = future))
            if(all(is.na(te$focal$clim))){
                  showModal(modalDialog(
                        title="Oops", "It seems your selected location is outside the allowed area. Please choose a site in California, and if focal site activity is set to 'seed collection', please choose a site within the species range (i.e. colored areas on the map).",
                        easyClose = TRUE, footer = modalButton("Dismiss") ))
                  updateSelectInput(session, "mode", selected = "planting")
            }
            message("target_envt")
            return(te)
      })
      
      range_stats <- reactive({
            y <- filter(range_summaries,
                        species == input$sp)
            list(mean = y$mean[c(6:10, 1:5)],
                 sd = y$sd[c(6:10, 1:5)])
      })
      
      sigmas <- reactive({
            x <- c(range_envt()$clim, range_envt()$soil)
            x <- (x - range_stats()$mean) / range_stats()$sd
            site_mean <- (unlist(c(target_envt()$focal$clim, target_envt()$focal$soil)) - range_stats()$mean) / range_stats()$sd
            w <- rep(c(input$pclim / 100, 1-(input$pclim / 100)), each = nlyr(x)/2)
            message("sigmas")
            sigma(x, site_mean, w)
      })
      
      k <- reactive({
            ifelse(is.null(input$k), 5, input$k)
      })
      
      clusters <- reactive({
            x <- c(range_envt()$clim, range_envt()$soil)
            f <- values(x)
            a <- which(is.finite(rowSums(f)))
            v <- scale(f[a,])
            w <- sqrt(rep(c(input$pclim / 100, 1-(input$pclim / 100)), each = nlyr(x)/2)*2)
            for(i in 1:ncol(v)) v[,i] <- v[,i] * w[i]
            set.seed(1)
            clust <- kmeans(v, k())
            y <- x[[1]] %>% setValues(NA) %>% setNames("clust")
            y[a] <- clust$cluster
            y
      })
      
      # overall dissimilarity combining soil and climate
      final <- reactive({
            # if(input$mode == "seed collection") browser()
            req(sigmas(), range_envt(), clusters())
            prob <- classify(sigmas(), matrix(c(5, Inf, 5), nrow = 1))
            d <- c(range_envt()$clim, range_envt()$soil, prob)
            names(d) <- all_vars$abbv[1:11]
            d <- c(d, clusters())
            if(!is.null(deltas())) d <- c(d, deltas())
            return(d)
      })
      
      
      
      
      
      # map ####################################
      
      # palettes 
      sigma_pal <- rev(c("darkred", "orange", "#FDE725FF", "#5DC863FF", "#21908CFF", "#3B528BFF", "#440154FF", "black"))
      viridis_pal <- rev(c("#FDE725FF", "#5DC863FF", "#21908CFF", "#3B528BFF", "#440154FF", "black"))
      delta_pal <- rev(c("darkred", "orange", "gray", "dodgerblue", "darkblue"))
      cluster_pal <- reactive({RColorBrewer::brewer.pal(k(), "Dark2")[1:k()]})
      
      icon <- makeAwesomeIcon("leaf", markerColor="red")
      
      output$map <- renderLeaflet({
            leaflet() %>%
                  setView(lng=-119.2, 37.2, lat=37.2, zoom=6) %>%
                  addProviderTiles(providers$Esri.WorldGrayCanvas)
      })
      
      observe({
            withProgress(
                  message = "Stand by.", 
                  value = 0, 
                  {
                        incProgress(.25, detail = "Computing environmental similarities.")
                        req(sigmas())
                        
                        v <- all_vars[all_vars$desc == input$color, ]
                        f <- final()
                        incProgress(1, detail = "Generating plots.")
                        
                        latlon <- coordinates(site$point)
                        req(final())
                        
                        r <- f[[v$abbv]]
                        if(input$color == "Seed zones"){
                              pal <- colorFactor(cluster_pal(), domain = 1:k(), na.color = "transparent")
                        }else{
                              limits <- range(c(v$min, v$max, values(r)), na.rm = T)
                              if(str_detect(input$color, "Change in")) limits <- max(abs(values(r)), na.rm = T) * c(-1, 1)
                              pal <- colorNumeric(switch(v$palette, "sigma" = sigma_pal, "viridis" = viridis_pal, "delta" = delta_pal),
                                                  limits,
                                                  na.color = "transparent")
                        }
                        leafletProxy("map") %>% clearMarkers() %>% clearImages() %>% clearControls() %>%
                              addRasterImage(r, colors = pal, opacity = 0.8,
                                             method = ifelse(v$abbv == "clust", "ngb", "bilinear")) %>%
                              addAwesomeMarkers(lng = latlon[1], lat = latlon[2], icon = icon)
                        
                  })
      })
      
      # scatter plot ####################################
      
      scat <- reactive({
            req(final())
            
            # future climate mean, and sd of either ensemble or niche
            avg <- c(target_envt()$focal$clim, target_envt()$focal$soil) %>% unlist()
            ref <- c(target_envt()$reference$clim, target_envt()$focal$soil) %>% unlist()
            std <- range_stats()$sd# %>% unlist()
            names(avg) <- names(std) <- names(ref) <- all_vars$abbv[1:length(avg)]
            
            vx <- all_vars$abbv[all_vars$desc == input$xvar]
            vy <- all_vars$abbv[all_vars$desc == input$yvar]
            vc <- all_vars$abbv[all_vars$desc == input$color]
            
            f <- final()
            d <- f[[c(vx, vy, vc)]] %>% as.data.frame() %>% na.omit()
            
            # compute range on full data set, and then subsample data
            minmax <- range(d[[vc]], na.rm = T)
            maxpoints <- 10000
            names(d) <- c("xvar", "yvar", "cvar")
            d <- d %>% sample_n(min(maxpoints, nrow(.)))
            
            fill_col <- NA
            
            x_label <- paste0(input$xvar, " ", all_vars$units[all_vars$desc == input$xvar])
            y_label <- paste(input$yvar, all_vars$units[all_vars$desc == input$yvar])
            
            # if(input$color == "Multivariate change in climate") browser()
            # plot
            req(d, vc)
            
            if(vc == "clust") d$cvar <- factor(d$cvar, levels = sort(unique(d$cvar)), labels = paste("zone", sort(unique(d$cvar))))
            # if(vc == "clust") d$cvar <- factor(d$cvar, levels = sort(unique(d$cvar)), labels = paste("zone", LETTERS[as.integer(sort(unique(d$cvar)))]))
            
            vci <- all_vars[all_vars$desc == input$color, ]
            scatter <- ggplot() +
                  geom_point(data = d, aes(xvar, yvar, color = cvar), size = .5) +
                  theme_minimal(base_size = 15) +
                  labs(x = x_label, y = y_label, color = NULL)
            
            if(vc == "clust"){
                  scatter <- scatter +
                        theme(legend.position = "none") +
                        scale_color_manual(values = cluster_pal())
            }else{
                  limits <- c(ifelse(is.na(vci$min), minmax[1], vci$min), 
                              ifelse(is.na(vci$max), minmax[2], vci$max))
                  if(str_detect(tolower(input$color), "change in") &
                     !str_detect(input$color, "Multivariate")) limits <- max(abs(minmax)) * c(-1, 1)
                  scatter <- scatter + 
                        theme(legend.position = "none") +
                        scale_color_gradientn(colours = switch(vci$palette, 
                                                               "sigma" = sigma_pal, "viridis" = viridis_pal, "delta" = delta_pal), 
                                              limits = limits)
            }
            
            
            # arrows
            if(!grepl("plot", paste(vx, vy))){
                  req(scatter)
                  if(input$mode == "planting"){
                        mods <- sapply(1:5, function(x) c(target_envt()$focal[[paste0("m", x)]], target_envt()$focal$soil) %>% unlist()) %>% t()
                        orig <- ref %>% as.matrix() %>% t()
                        dest <- mods
                  }else{
                        mods <- sapply(1:5, function(x) c(target_envt()$reference[[paste0("m", x)]], target_envt()$reference$soil) %>% unlist()) %>% t()
                        orig <- mods
                        dest <- avg %>% as.matrix() %>% t()
                  }
                  
                  scatter <- scatter +
                        annotate("point", color = "black", size = 6, x = avg[vx], y = avg[vy]) +
                        annotate("point", color = "red", size = 5, x = avg[vx], y = avg[vy])
                  
                  if(input$time != "1981-2010") scatter <- scatter +
                        annotate("segment", color="black", linewidth = 1,
                                 x=orig[,vx], y=orig[,vy], xend=dest[,vx], yend=dest[,vy],
                                 arrow=grid::arrow(type="closed", angle=15, length=unit(.1, "in"),
                                                   ends = ifelse(input$mode == "seed collection", "first", "last"))) +
                        annotate("point", color="red", size = .5, x=orig[vx], y=orig[vy]) +
                        annotate("segment", color="red", linewidth = .3,
                                 x=orig[,vx], y=orig[,vy], xend=dest[,vx], yend=dest[,vy],
                                 arrow=grid::arrow(type="open", angle=15, length=unit(.1, "in"),
                                                   ends = ifelse(input$mode == "seed collection", "first", "last")))
            }
            
            return(scatter)
      })
      
      output$scatter <- renderPlot({
            scat()
      })
      
      # color legend #############
      lgnd <- reactive({
            p <- scat() +
                  theme(legend.position = "top",
                        legend.justification = "left")
            if(input$color == "Seed zones"){
                  p <- p + guides(color = guide_legend(override.aes = list(size = 8, shape = 15)))
            }else{
                  p <- p + guides(color = guide_colorbar(barwidth = unit(.75, "npc"),
                                                         barheight = .75,
                                                         title.position = "top"))
            }
            # browser()
            # cowplot::get_plot_component(plot, "guide-box", return_all = TRUE)
            cowplot::get_legend(p)
      })
      output$legend1 <- renderPlot({ grid.draw(lgnd()) })
      output$legend2 <- renderPlot({ grid.draw(lgnd()) })
      
      
      # text elements ############
      point_info <- reactive({ifelse(input$mode == "planting", 
                                     paste0("(", times[1], ifelse(input$radius == 0, "", ", smoothed)")), 
                                     paste0("(", input$time, ", ", input$ssp, ")"))})
      color_info <- reactive({ paste(input$color,
                                     all_vars$units[all_vars$desc == input$color],
                                     ifelse(input$color %in% all_vars$desc[c(1:10, 12)], point_info(), "")) })
      output$color_label1 <- renderText({ color_info() })
      output$color_label2 <- renderText({ color_info() })
      output$target_label <- renderText({ paste0("Large marker: ", input$mode, " site (",
                                                 ifelse(input$mode == "planting", 
                                                        paste0(input$time, ", ", input$ssp), 
                                                        paste0(times[1], ifelse(input$radius == 0, "", ", smoothed"))), 
                                                 ")") })
      output$points_label <- renderText({ paste0("Points: potential ", setdiff(c("planting", "seed collection"), input$mode), 
                                                 " sites ", point_info()) })
      
      
      # downloads ########################
      
      # data download
      output$download <- downloadHandler(
            filename = function() {
                  if(input$color == "Environmental difference"){
                        ll <- coordinates(site$point)
                        f <- paste0("sigma", 
                                    "_", str_replace(input$sp, " ", "-"), 
                                    "_", input$mode, 
                                    "_pclim", input$pclim, 
                                    "_r", input$radius, "km",
                                    "_", input$time, 
                                    "_", ssps$ssp[ssps$text == input$ssp], 
                                    "_", round(ll[1], 2), "E",
                                    "_", round(ll[2], 2), "N",
                                    ".tif")
                  }else if(input$color == "Seed zones"){
                        f <- paste0("zones", 
                                    "_", str_replace(input$sp, " ", "-"), 
                                    "_k", input$k, 
                                    "_pclim", input$pclim, 
                                    "_r", input$radius, "km",
                                    ".tif")
                  }else if(str_detect(tolower(input$color), "change in")){
                        f <- paste0(input$color,
                                    "_", str_replace(input$sp, " ", "-"), "_",
                                    paste0(input$time, "_", input$ssp),
                                    ".tif")
                  }else{
                        f <- paste0(input$color,
                                    "_", str_replace(input$sp, " ", "-"), "_",
                                    ifelse(input$mode == "planting", 
                                           paste0(times[1], ifelse(input$radius == 0, "", "smoothed")), 
                                           paste0(input$time, "_", input$ssp)),
                                    ".tif")
                  }
                  f
            },
            content = function(file) {
                  y <- final()[[all_vars$abbv[all_vars$desc == input$color]]]
                  writeRaster(y, file)
            }
      )
      
      # species list download
      output$downloadSpp <- downloadHandler(
            filename = "CA_native_vascular_plants.csv",
            content = function(file) {
                  write.csv(data.frame(species = allspps), file)
            }
      )
      
      
}

# Run the application
shinyApp(ui = ui, server = server)
