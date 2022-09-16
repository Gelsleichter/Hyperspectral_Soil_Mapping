######## Fluxogram design, Yuri Sep 2022

#' ### Auto setwd
library(rstudioapi); current_path <- getActiveDocumentContext()$path; setwd(dirname(current_path)) ### https://eranraviv.com/r-tips-and-tricks-working-directory/
getwd()
# install.packages("DiagrammeR")

library(DiagrammeR)
### Para Artigo Geoderma
flux_sp_model <- grViz("
# digraph boxes_and_circles {
digraph R {

##################### align
# rankdir=LR # left - right BL, BR, TL, TR, RB, RT, LB, LT
rankdir=TB   # top - bottom
scale = both

# a 'graph' statement
graph [overlap = true, fontsize = 12]

###########################
### several 'node' statements
node [shape = box3d, fontname = Helvetica, fillcolor = Lavender, style = filled, margin = 0.25, fontsize= 12]
'DSM \n (Digital Soil Mapping)'

node [shape = cylinder, fontname = Helvetica, fillcolor = Wheat4, margin = 0.15]
'Total Carbon soil data \n (wet chemistry lab)'

node [shape = box3d, fontname = Helvetica, fillcolor = Peru, margin = 0.25]
'Soil Map \n DSM method'

node [shape = cylinder, fontname = Helvetica, fillcolor = LemonChiffon, margin = 0.15]
'Input covariates for DSM: \n covariates set 1: \n Multispectral \n terrain \n geographic'

node [shape = cylinder, fontname = Helvetica, fillcolor = LemonChiffon, margin = 0.15]
'Input covariates for \n Hyperspectral Sub. Image: \n covariates set 2: \n Multispectral \n terrain'

node [shape = cylinder, fontname = Helvetica, style = filled, fillcolor = Goldenrod, margin= 0.15]
'Proximal Soil Sensing \n (soil lab V-SWIR spectra)'

node [shape = tab, fontname = Helvetica, fillcolor = Cyan, margin = 0.15]
'Processe to generate 130 \n Hyperspectral Sub. Image \n covariates'

node [shape = cylinder, fontname = Helvetica, fillcolor = Green3, margin = 0.12]
'Input covariates for HSM: \n covariates set 3: \n set 1 and \n Hyperspectral Sub. Image'

node [shape = cylinder, fontname = Helvetica, fillcolor = Wheat4, margin = 0.15]
'Total Carbon soil data \n (wet chemistry lab)'

node [shape = box3d, fontname = Helvetica, fillcolor = SlateBlue1, style = filled, margin = 0.25]
'HSM \n (Hyperspectral Soil Mapping)'

node [shape = box3d, fontname = Helvetica, fillcolor = Peru, margin = 0.25]
'Soil Map \n HSM method'

#############################
### several 'edge' statements
'Input covariates for \n Hyperspectral Sub. Image: \n covariates set 2: \n Multispectral \n terrain' -> 'Processe to generate 130 \n Hyperspectral Sub. Image \n covariates'

'Proximal Soil Sensing \n (soil lab V-SWIR spectra)' -> 'Processe to generate 130 \n Hyperspectral Sub. Image \n covariates'

'Processe to generate 130 \n Hyperspectral Sub. Image \n covariates' -> 'Input covariates for HSM: \n covariates set 3: \n set 1 and \n Hyperspectral Sub. Image'

'Input covariates for HSM: \n covariates set 3: \n set 1 and \n Hyperspectral Sub. Image' -> 'HSM \n (Hyperspectral Soil Mapping)'

'Total Carbon soil data \n (wet chemistry lab)' -> 'HSM \n (Hyperspectral Soil Mapping)'

'HSM \n (Hyperspectral Soil Mapping)' -> 'Soil Map \n HSM method'[label = ' ', fontcolor = white] 

'Input covariates for DSM: \n covariates set 1: \n Multispectral \n terrain \n geographic' -> 'DSM \n (Digital Soil Mapping)'[minlen = 3]

'Total Carbon soil data \n (wet chemistry lab)' -> 'DSM \n (Digital Soil Mapping)'[minlen = 1]

'DSM \n (Digital Soil Mapping)' -> 'Soil Map \n DSM method'

}
")

flux_sp_model

# the "label = ' '" is needed to "rotate" the plot
# https://epirhandbook.com/en/diagrams-and-charts.html

### Exportar: Fonte ### https://github.com/rich-iannone/DiagrammeR/issues/70 ### https://stackoverflow.com/questions/57078753/how-to-save-grviz-object-from-diagrammer-r-package
### sudo apt install libv8-dev ### para pacote "DiagrammeRsvg"
### sudo apt install librsvg2-dev ### para pacote "rsvg"
# install.packages("devtools")
# library(devtools)
# devtools::install_github("rich-iannone/DiagrammeR")
# devtools::install_github("rich-iannone/DiagrammeRsvg")
# install.packages("DiagrammeRsvg")
# install.packages("rsvg")
library(DiagrammeR)
library(DiagrammeRsvg)
library(magrittr)
library(rsvg)
library(dplyr)

# /home/gelsleichter/Dropbox/Aplicativos/Overleaf/Gelsleichter_HSM_2020_Geoderma_ATUAL_elsarticle-template/figures
# grViz(diagram = flux_sp_model[[1]]$diagram) %>% export_svg %>% charToRaw %>% rsvg_png("graph.png", width = 5000, height = 5000) ### muito pesado em alta resolução
# grViz(diagram = flux_sp_model[[1]]$diagram) %>% export_svg %>% charToRaw %>% rsvg_svg("graph.svg") ### svg quase não uso
# grViz(diagram = flux_sp_model[[1]]$diagram) %>% export_svg %>% charToRaw %>% rsvg_pdf("flux_sp_model_gd.pdf")
# grViz(diagram = flux_sp_model[[1]]$diagram) %>% export_svg %>% charToRaw %>% rsvg_png("flux_sp_model_gd.png")
grViz(diagram = flux_sp_model[[1]]$diagram) %>% export_svg %>% charToRaw %>% rsvg_pdf("P:/Artigos/2022_Artigo_Geoderma_Regional_HSM_PNI/output/graphs_plots/flux_sp_model_gd.pdf")

### salvar direto na pasta do artigo
grViz(diagram = flux_sp_model[[1]]$diagram) %>% export_svg %>% charToRaw %>% rsvg_pdf("C:/Users/Gel8695/Dropbox/Aplicativos/Overleaf/Gelsleichter_etal_HSM_2022_v2_Geoderma_Regional/figures/flux_sp_model_gd.pdf")
# grViz(diagram = flux_sp_model[[1]]$diagram) %>% export_svg %>% charToRaw %>% rsvg_pdf("/home/gelsleichter/Dropbox/Aplicativos/Overleaf/Gelsleichter_HSM_2020_Geoderma_ATUAL_elsarticle-template/figures/flux_sp_model_gd.pdf")
# grViz(diagram = flux_sp_model[[1]]$diagram) %>% export_svg %>% charToRaw %>% rsvg_png("/home/gelsleichter/Dropbox/Aplicativos/Overleaf/Gelsleichter_HSM_2020_Geoderma_ATUAL_elsarticle-template/figures/flux_sp_model_gd.png")



