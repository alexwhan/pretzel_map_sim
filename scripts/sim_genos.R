library(qtl)
library(purrr)
library(magrittr)
library(tidyr)
library(genomap)

set.seed(123)
initial_map <- read.csv("data/initial_genetic_map.csv", stringsAsFactors = FALSE)
initial_map$lg <- as.factor(initial_map$lg)
lgs <- table(initial_map$lg)
marker_names <- split(initial_map$marker_name, initial_map$lg)

pmap <- sim.map(lgs * rnorm(21, 1, 0.1), lgs, include.x = FALSE, sex.sp = FALSE, eq.spacing = FALSE)
names(pmap) <- names(lgs)

pcross <- sim.cross(pmap, n.ind = 100, type = "riself")
pcross_geno <- map2(pcross$geno, marker_names, function(.x, .y) set_colnames(.x$data, .y)) %>% 
  map(as.data.frame) %>% 
  map(~mutate(.x, germplasm_id = paste("sim", 1:100, sep = "_"))) %>% 
  map(~gather(.x, marker_name, genotype, -germplasm_id)) %>% 
  bind_rows

write.csv(pcross_geno, "data/simulated_geno.csv", row.names = FALSE)
