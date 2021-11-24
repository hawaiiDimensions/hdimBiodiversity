
## non-looping version, just selecting Araneae to start
## filter to Araneae species sum by site
## DSG Oct to Nov 2021
# relies on the "endemic" object from "hdim-beating-analysis.Rmd" (Natalie Graham)

# spoiler alert: coverage probably was never going to be informative for metabarcode reads data, even following correction to remove 4 to create singletons. There just are so few singletons or doubletons in the dataset that coverage will be close to 1 or 100% with few individuals

## 
library(iNEXT)
library(tidyr)
library(dplyr)
otutable <- read.csv("data/hdimbeat_abund.csv")
endemic <- subset(otutable, predicted_status == "endemic")
endemic$delim_spp <- endemic$spp

################################################################
####  1
####  cycle through the Orders, calculate adjusted abundances for species, sum across plots & sites
####  stash each abundance distribution in a list, with a vector for each order
####  use iNEXT estimate(D) to calculate effective diversity at q=0,1,2 for each order by coverage level
####  inside the for loop, funky data structure for iNEXT created
####  - list with vector for each order, but descending abundance & first element is # of plots


## Orders overall pooled across plots and sites
# Create an empty list
coverlist <- c()
# function to select & calc each order
for (i in unique(endemic$orderGrp)) {
inextCoverage <- function(o){
    orderMtrx <- endemic %>%
    select(orderGrp, delim_spp, abundance) %>% # select features; remove site, sitePlot
    filter(orderGrp == o) %>% #filter by order
    group_by(delim_spp) %>% #remove site, sitePlot
    mutate(abund_adj = abundance - 4) %>% #adjust before summing
    summarise(abund_adj = sum(abund_adj)) %>% # sum abundances
    #select(-abundance) %>%
    pivot_wider(names_from = delim_spp, values_from = abund_adj)
orderMtrx
orderMtrx[is.na(orderMtrx)] <- 0  #rm NAs

# drop order column, transpose, sort
y <- orderMtrx %>%
       t()
colnames(y) <- c(o)
y <- sort(y, decreasing = TRUE)

# create list object where first element of each vector is # samples; the rest are the species abundances in descending order
coverlist[[i]] <- c(length(unique(endemic$sitePlot)) , y)
    }
}

coverlist

#inextCoverage("Orthoptera")
coverlist <- lapply(unique(endemic$orderGrp), inextCoverage)
names(coverlist) <- unique(endemic$orderGrp)

# create iNEXT object
# estimate diversity @ q=0,1,2 for each order under different coverage
D80 <- estimateD(coverlist, datatype="abundance", base="coverage", level=0.80, conf=0.95)
D90 <- estimateD(coverlist, datatype="abundance", base="coverage", level=0.90, conf=0.95)
D95 <- estimateD(coverlist, datatype="abundance", base="coverage", level=0.95, conf=0.95)

#write_csv(D80, here("data/D80cover_orders.csv"))
#write_csv(D90, here("data/D90cover_orders.csv"))
#write_csv(D95, here("data/D95cover_orders.csv"))




################################################################
####  2
####  Loop through Orders species X site matrices
####  don't know why plots are not created
####  two problems: (1) too many lines/colors to plot in ggiNEXT
####  (2) iNEXT function warns "In if (class(x) == "numeric") { ... :
####  the condition has length > 1 and only the first element will be used"
####  perhaps a warning that first element does not have #plots in each vector
####  the non-loop version (single order) ignores and plots anyway. Maybe 13 too many
####

## Orders By Site
## loop with lapply for 13 orders. Drowns in errors.
#o = "Araneae"
# function to select & calc each order
inextCoverage2 <- function(o){
    O.Mtrx <- endemic %>%
        select(orderGrp, site, delim_spp, abundance) %>%
        group_by(orderGrp, site, delim_spp) %>%
        filter(orderGrp == o) %>% # filter by order
        #select(-orderGrp) %>%  #drop column
        mutate(abund_adj = abundance - 4) %>%
        select(-abundance) %>%
        summarise(abund_adj = sum(abund_adj)) %>%
        pivot_wider(names_from = delim_spp, values_from = abund_adj)
    O.Mtrx[is.na(O.Mtrx)] <- 0  #rm NAs
    O.Mtrx

    #
    mat <- O.Mtrx %>%
        #mutate_if(is.numeric, replace_na, replace = 0) %>%
        column_to_rownames("site") %>%
        select(-orderGrp) %>%
        t() #transpose to S x N

    cover.o <- iNEXT(mat, q = 0, datatype = "abundance")
    #cover.o$DataInfo  #f1 and f2 in the output are singletons and doubletons
    #cover.o$AsyEst
    #cover.o$iNextEst

    # use ggiNEXT to plot
ggiNEXT(cover.o, type = 1, se = TRUE, facet.var = "none",
            color.var = "site", grey = FALSE)

}

lapply(unique(endemic$orderGrp), inextCoverage2)



################################################################
####   3   As the loop is not functional (too many errors?)
####   Create Species X Sample (sites) Matrix for each orderGrp
####   Plot rarefaction/interpolation curves and calc coverage by site
####

# order by order pooled by site - no loop
o <- "Araneae"
#o <- "Hemiptera"
#o <- "Malacostraca"
#o <- "Myriapoda"
#o <- "Diptera"
#o <- "Lepidoptera"
#o <- "Hymenoptera"
#o <- "Coleoptera"
#o <- "Psocoptera"
#o <- "Neuroptera"
#o <- "Orthoptera"
#o <- "Collembola"
#o <- "non-spider Arachnida"
#o <- "Orthoptera"

O.Mtrx <- endemic %>%
    select(orderGrp, site, delim_spp, abundance) %>%
    group_by(orderGrp, site, delim_spp) %>%
    filter(orderGrp == o) %>% # filter by order
    #select(-orderGrp) %>%  #drop column
    mutate(abund_adj = abundance - 4) %>%
    select(-abundance) %>%
    summarise(abund_adj = sum(abund_adj)) %>%
    pivot_wider(names_from = delim_spp, values_from = abund_adj)
O.Mtrx[is.na(O.Mtrx)] <- 0  #rm NAs
O.Mtrx

mat <- O.Mtrx %>%
    #mutate_if(is.numeric, replace_na, replace = 0) %>%
    column_to_rownames("site") %>%
    select(-orderGrp) %>%
    t() #transpose to S x N

coverout <- iNEXT(mat, q = 0, datatype = "abundance")
coverout[[1]]  #f1 and f2 in the output are singletons and doubletons
#coverout$AsyEst
#coverout$iNextEst

# use ggiNEXT to plot
ggiNEXT(coverout, type = 1, se = TRUE, facet.var = "none",
        color.var = "site", grey = FALSE)

