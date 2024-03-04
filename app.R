
library(shiny)
library(shinydashboard)
library(shinyalert)
library(shinyBS)
library(dplyr)
library(tidyr)
library(ggplot2)
library(leaflet)
library(raster)
library(terra)
library(grid)
library(purrr)
library(stringr)
select <- dplyr::select
extract <- terra::extract


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
                                "prob", "clust", paste0("d", vars)))
all_vars$desc <- c("Actual evapotranspiration", "Climatic water deficit",
                   "Winter minimum temperature", "Summer maximum temperature",
                   "Annual precipitation",
                   "Soil PC1", "Soil PC2", "Soil PC3", "Soil PC4", "Soil PC5",
                   "Environmental difference",
                   "Seed zones",
                   paste0("Change in ", c("Actual evapotranspiration", "Climatic water deficit",
                                          "Winter minimum temperature", "Summer maximum temperature", "Annual precipitation")))
all_vars$units <- c("(mm)", "(mm)",
                    "(deg. C)", "(deg. C)",
                    "(log10 mm)",
                    rep("(standard deviations)", 5),
                    rep("(similarity to target, standard deviations)", 1),
                    "(environmental clusters)",
                    "(mm)", "(mm)", "(deg. C)", "(deg. C)", "(log10 mm)")
all_vars$min <- c(rep(NA, 10), rep(0, 1), NA, rep(NA, 5))
all_vars$max <- c(rep(NA, 10), rep(5, 1), NA, rep(NA, 5))
all_vars$palette <- c(rep("viridis", 10), rep("sigma", 1), "cluster", rep("delta", 5))

sigma <- function(x, site_mean, w){
      w <- w / sum(w)
      k <- length(site_mean)
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
                  menuItem("Tool", tabName = "tool", icon = icon("map-location")),
                  menuItem("Instructions", tabName = "instructions", icon = icon("circle-info")),
                  menuItem("About", tabName = "about", icon = icon("circle-user"))
            ),
            br(),
            br(),
            p(style = "margin-left: 15px;margin-right: 15px;", 
              "Seeds of Change is a decision support tool for ecological restoration. It identifies areas in California that may be suitable for seed collection or planting, by analyzing spatial data on climate change, soil characteristics, species distributions, and adaptive neighborhoods.")
      ),
      
      dashboardBody(
            
            tags$head(
                  tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
            ),
            
            tabItems(
                  tabItem(tabName = "tool",
                          fluidRow(
                                column(
                                      width = 4,
                                      box(
                                            title = "Model settings", width = NULL, #collapsible = T, collapsed = T,
                                            selectInput("mode", 
                                                        span("Focal site activity ", actionLink("i_mode", "[?]")), 
                                                        c("planting", "seed collection"), "planting"),
                                            uiOutput("constraintControls"),
                                            textInput("lonlat",
                                                      span("Focal site location ", actionLink("i_location", "[?]")), 
                                                      paste(s$lon, s$lat, sep = ", ")),
                                            selectizeInput("sp", 
                                                           span("Species ", actionLink("i_species", "[?]")), 
                                                           choices = spps, selected = "Quercus agrifolia"),
                                            selectInput("time",
                                                        span("Time period ", actionLink("i_time", "[?]")),
                                                        times, "2071-2100"),
                                            selectInput("ssp",
                                                        span("Scenario ", actionLink("i_ssp", "[?]")),
                                                        ssps$text, "SSP5-8.5"),
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
                                                        width="100%"),
                                            selectizeInput("color", "Display variable",
                                                           all_vars$desc[c(11, 12, 1:10, 13:17)],
                                                           all_vars$desc[11]),
                                            uiOutput("nclust")
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
                                                           column(8, 
                                                                  span(textOutput("color_label1")),
                                                                  plotOutput("legend2", height = "45px"),
                                                                  plotOutput("scatter")),
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
                                                     span(textOutput("color_label2")),
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
      
      observeEvent(input$i_species,
                   {showModal(modalDialog(
                         title="Species",
                         HTML("The tool comes pre-loaded with estimated geographic range maps for most California native plant species (only 25 species important for Bay Area restoration projects are included in this beta version, but all 5221 species listed below will ultimately be added).",
                              "Select a species to load its estimated California geographic range, which is modeled based on climate, distance to known observations, and landscape intactness.",
                              "This distribution is used to estimate the species' environmental tolerance, model gene flow among nearby populations, and hypothesize which populations may be adapted to the environment of the selected focal site.",
                              "Or for a generic species-agnostic analysis of environmental similarity between sites, select 'NONE' in the species box.<br><br>"),
                         selectizeInput("allsp", "Full species list", choices = allspps), # not used server-side -- just informational
                         easyClose = TRUE, footer = modalButton("Dismiss") )) })
      observeEvent(input$i_radius, 
                   {showModal(modalDialog(
                         title="Smoothing radius (km)",
                         "A source population's suitability as a genetic match for a planting site with the same environment depends on the population's historic balance between selection and gene flow.",
                         "The relative importance of gene swamping versus local adaptation is known to vary among species.",
                         "This parameter lets you set the size of the local neighborhood around each source population within which gene flow is expected to homogenize local adaptation.",
                         "Choose a small smoothing radius to model highly localized adaptation, or a large radius to model strong effects of widespread gene flow relative to local selection.",
                         easyClose = TRUE, footer = modalButton("Dismiss") )) })
      observeEvent(input$i_time, 
                   {showModal(modalDialog(
                         title="Time period",
                         "This tool estimates how well the historic adaptive environment at each provenance site matches the projected future environment at the planting site, or the inverse, if site activity is set to `seed collection.`",
                         "For which future time period would you like to estimate climatic similarity?",
                         "You can also select the historic period, in which case baseline climate data is used for both the focal site and species range",
                         "(Note that soil data does not change across time periods.",
                         easyClose = TRUE, footer = modalButton("Dismiss") )) })
      observeEvent(input$i_ssp, 
                   {showModal(modalDialog(
                         title = "Intergovernmental Panel on Climate Change (IPCC) emissions scenarios (Shared Socio-economic Pathways (SSP))",
                         "The IPCC has defined a range of potential future climate trajectories based on different assumptions about climate change mitigation. Select which of these SSPs you would like to compare to historic climate.",
                         br(), br(), "SSP1-2.6 - effective sustainability policies, very low greenhouse gas emissions, net zero by 2050, global average temperature increase in year 2100 ~1.8 deg. C (range 1.3-2.4 deg. C) above pre-industrial",
                         br(), "SSP3-7.0 - regional imbalance, no additional policies, high greenhouse gas emissions, doubling 2020 to 2100, global average temperature increase in year 2100 ~3.6 deg. C (range 2.8-4.6 deg. C) above pre-industrial",
                         br(), "SSP5-8.5 - fossil fuel-intensive development, very high greenhouse gas emissions, doubling 2020 to 2050, global average temperature increase in year 2100 ~4.4 deg. C (range 3.3-5.7 deg. C) above pre-industrial",
                         easyClose = TRUE, footer = modalButton("Dismiss") )) })
      observeEvent(input$i_weight, 
                   {showModal(modalDialog(
                         title="Soils versus climate",
                         "By default, soil and climate are given equal weight when calculating similarity to the focal site environment.", 
                         "Adjust this slider to incorporate species-specific knowledge about their relative importance in shaping local adaptation,",
                         "or to explore how the results change based on assumptions about their importance.",
                         easyClose = TRUE, footer = modalButton("Dismiss") )) })
      observeEvent(input$i_constrain, 
                   {showModal(modalDialog(
                         title="Limit to species range",
                         "This option is available when focal site activity is set to 'seed collection.'",
                         "Check this box to limit the map of prospective planting sites to locations within the species' current modeled range.",
                         "Uncheck it to consider all of California, which may be useful if the modeled range omits areas where the species occurs, or as a means to explore where locations outside the current range could become suitable in the future.", 
                         easyClose = TRUE, footer = modalButton("Dismiss") )) })
      observeEvent(input$i_mode, 
                   {showModal(modalDialog(
                         title = "Target site activity",
                         "Select whether the selected focal location is a 'planting' or 'seed collection' site.",
                         "If it's a planting site, its environment in the selected 'time period' will be compared to historic smoothed environments across the species range (prospective locations where seeds could be collected to plant in the focal site).", 
                         "If it's a collection site, its historic smoothed environment will be compared to range-wide environments from the selected time period (prospective locations where seeds collcted in the focal site could be planted).",
                         easyClose = TRUE, footer = modalButton("Dismiss") )) })
      observeEvent(input$i_download, 
                   {showModal(modalDialog(
                         title = "Download results",
                         "Click the download button to export a raster file of the current map.",
                         easyClose = TRUE, footer = modalButton("Dismiss") )) })
      observeEvent(input$i_location, 
                   {showModal(modalDialog(
                         title = "Select a focal location",
                         "To choose a focal site, either click on the map, or enter your 'Lon, Lat' in the box and press RETURN on your keyboard.",
                         easyClose = TRUE, footer = modalButton("Dismiss") )) })
      observeEvent(input$i_nclust, 
                   {showModal(modalDialog(
                         title = "Seed zones",
                         "Use this slider to select the number of seed zones to produce. This performs a k-means cluster analysis, grouping sites across the species range into clusters with similar environments.",
                         "Seed zones are independent of the selected focal site, but the 'focal site activity' setting determines which environmental data are used:",
                         "if the activity is set to 'planting' is then clusters are based on smoothed historic data, while if it is set to 'seed collection' then clusters will be based on environmental data for the selected future scenario.",
                         "The 'soil versus climate' setting affects the relative weight given to the two categories of variable.",
                         easyClose = TRUE, footer = modalButton("Dismiss") )) })
      observeEvent(input$i_arrows, 
                   {showModal(modalDialog(
                         title = "Focal site",
                         "The red arrows indicate projected climate change at the focal site for the selected era and emissions scenario, according to five different global climate models (GCMs).",
                         "The large red point indicates the focal site's environment for the era listed in red above the plot, which is what is compared to other sites across the species range to calculate 'environmental difference'; for future time periods this is the average of the GCMs.",
                         easyClose = TRUE, footer = modalButton("Dismiss") )) })
      observeEvent(input$i_threshold, 
                   {showModal(modalDialog(
                         title = "Species distribution model (SDM) threshold",
                         "The species range maps in this tool are based on models that estimate continuous gradients of occurrence likelihood, in which every species has suitability values ranging from 0 to 100 across the state.",
                         "These are not probabilities; they represent occurrence scores as a percentage of the highest occurrence score for the species anywhere in the state.",
                         "This tool 'thresholds' these suitability maps, displaying only those areas with a modeled occurrence score higher than the value set by the `SDM threshold` slider.",
                         "Choose a lower value to select a broader range (which will likely include more actual populations, but probably also more places the species doesn't actually exist).",
                         "Choose a higher value to select a narrower range (which will reduce the inclusion of places the species doesn't actually exist, but may also omit more true populations).",
                         easyClose = TRUE, footer = modalButton("Dismiss") )) })
      
      
      
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
            # site$ll <- coordinates(s)
      })
      
      # show range constraint input only in collection mode
      # output$constraintControls <- renderUI({
      #       if(input$mode == "seed collection") checkboxInput("constrain", span("Limit results to species range", actionLink("i_constrain", "[?]")), TRUE)
      # })
      output$constraintControls <- renderUI({
            if(input$mode == "seed collection") selectInput("constrain", span("Limit results to species range", actionLink("i_constrain", "[?]")), c(TRUE, FALSE))
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
                                       choices = all_vars$desc[c(11, 12, 1:10, 13:17)], 
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
            if(! str_detect(input$color, "Change in")) return(NULL)
            files <- list.files("assets/deltas", full.names = T)
            files <- files[grepl(time, files) & 
                                 grepl(tolower(ssp), files) &
                                 grepl(str_remove(var, "d"), files)]
            msk <- smoothed_envt()[[1]][[1]]
            rast(files) %>% crop(msk) %>% mask(msk) %>% setNames(var) %>% clip(crop(range_mask(), .)) # clip(range_mask())
      })
      
      # reference climate across range (historic if planting mode)
      range_envt <- reactive({
            y <- switch(input$mode,
                        "planting" = smoothed_envt(),
                        "seed collection" = future_envt() )
            if(input$mode == "seed collection"){
                  req(input$constrain)
                  if(!is.null(input$constrain)){
                        if(input$constrain){
                              y <- y %>% map(clip, y = range_mask())
                              # map(crop, y = smoothed_envt()[[1]][[1]]) %>%
                              # map(terra::mask, mask = smoothed_envt()[[1]][[1]])
                        }#else{
                        #browser()
                        #}
                  }
            }
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
      sigma_pal <- c("red", "#FDE725FF", "#5DC863FF", "#21908CFF", "#3B528BFF", "#440154FF", "black")
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
                        
                        #                         tag.map.title <- tags$style(HTML("
                        #   .leaflet-control.map-title { 
                        #     transform: translate(45px,-65px);
                        #     top: 0px;
                        #     left: 0px;
                        #     text-align: left;
                        #     padding-left: 10px; 
                        #     padding-right: 10px; 
                        #     background: rgba(255,255,255,0.75);
                        #     font-weight: bold;
                        #     font-size: 16px;
                        #   }
                        # "))
                        #                         title <- tags$div(tag.map.title, HTML(paste0(input$color, ":<br><i>", input$sp, "</i>")))
                        
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
                              # addControl(title, position = "topleft", className="map-title") %>%
                              addAwesomeMarkers(lng = latlon[1], lat = latlon[2], icon = icon)
                        
                  })
      })
      
      # scatter plot ####################################
      # output$scatter <- renderPlot
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
            
            
            # plot
            req(d, vc)
            
            if(vc == "clust") d$cvar <- factor(d$cvar, levels = sort(unique(d$cvar)), labels = paste("zone", sort(unique(d$cvar))))
            
            vci <- all_vars[all_vars$desc == input$color, ]
            scatter <- ggplot() +
                  geom_point(data = d, aes(xvar, yvar, color = cvar), size = .5) +
                  theme_minimal(base_size = 15) +
                  labs(x = x_label, y = y_label, color = NULL)
            
            if(vc == "clust"){
                  scatter <- scatter +
                        # guides(color = guide_legend(override.aes = list(size = 8, shape = 15))) +
                        scale_color_manual(values = cluster_pal())
            }else{
                  limits <- c(ifelse(is.na(vci$min), minmax[1], vci$min), 
                              ifelse(is.na(vci$max), minmax[2], vci$max))
                  if(str_detect(input$color, "Change in")) limits <- max(abs(minmax)) * c(-1, 1)
                  scatter <- scatter + 
                        theme(legend.position = "none") +
                        # guides(color = guide_colorbar(barheight = unit(.8, "npc"))) +
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
      
      
      output$color_label1 <- renderText({ paste(input$color,
                                                all_vars$units[all_vars$desc == input$color]) })
      output$color_label2 <- renderText({ paste(input$color,
                                                all_vars$units[all_vars$desc == input$color]) })
      
      lgnd <- reactive({
            p <- scat() +
                  theme(legend.position = "top",
                        legend.justification = "left")
            if(input$color == "Seed zones"){
                  p <- p + guides(color = guide_legend(override.aes = list(size = 8, shape = 15)))
            }else{
                  p <- p + guides(color = guide_colorbar(barwidth = unit(.5, "npc"), # unit(.8, "npc"),
                                                         barheight = .75,
                                                         title.position = "top"))
            }
            require(ggpubr)
            leg <- get_legend(p)
            as_ggplot(leg)
      })
      output$legend1 <- renderPlot({ lgnd() })
      output$legend2 <- renderPlot({ lgnd() })
      
      
      output$target_label <- renderText({ paste0("Target: ", input$mode, " site (",
                                                 ifelse(input$mode == "planting", 
                                                        paste0(input$time, ", ", input$ssp), 
                                                        paste0(times[1], ifelse(input$radius == 0, "", ", smoothed"))), 
                                                 ")") })
      output$points_label <- renderText({ paste0("Points: potential ", setdiff(c("planting", "seed collection"), input$mode), " sites (",
                                                 ifelse(input$mode == "planting", 
                                                        paste0(times[1], ifelse(input$radius == 0, "", ", smoothed")), 
                                                        paste0(input$time, ", ", input$ssp)), ")") })
      
      
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
                  }else if(str_detect(input$color, "Change in")){
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
                  # y <- final()$prob
                  # y <- final()$clust
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
