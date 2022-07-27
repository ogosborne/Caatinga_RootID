library(RootID)
library(rgl)
library(RColorBrewer)

# load data
load("results/processed_data/best.match.md0.RData")
metadata <- read.table("data/metadata.csv", sep = ",", header = T, na.strings = "NA", stringsAsFactors = F)
# make ind.list
ind.list <- list()
for(i in sort(unique(metadata$species))){
  ind.list[[i]] <- metadata[which(metadata$species == i), "sample"]
}
rm(i)
# create output directory
dir.create("results/3d_plots/")
# make root.pos
samps <- best.matches[["species.diag.markers"]]$sample
samps <- c(samps[grepl("_50", samps, fixed = T)], samps[grepl("_20", samps, fixed = T)], samps[grepl("_10", samps, fixed = T)], samps[grepl("_05", samps, fixed = T)])
root.pos <- array(samps, dim = c(5, 5, 4))
rm(samps)
# tree data, convert coordinates into units of one square (i.e. 2 M) starting from x=0,y=0. Divide heights by 2 so they match square width
tree.dat <- data.frame(sample = metadata$sample, 
                       x = (metadata$tree.pos.x - 2) / 2,
                       y = (metadata$tree.pos.y - 2) / 2,
                       rad =  apply(metadata[, c("canopy_diameter_NS","canopy_diameter_EW")], 1, function(x) max(x)/2) / 2,
                       height = metadata$h_tree / 2,
                       crown.height = metadata$h_crown / 2)

# Fig. 2
open3d()
# manually resize rgl window
mfrow3d(2,2,sharedMouse = T)
# Commiphora_leptophloeos sp
plot_roots_3d(matches = best.matches, 
              ind.list = ind.list, 
              root.pos = root.pos, 
              taxa = "Commiphora_leptophloeos", 
              tax.type = "sp", 
              draw.root.grid = TRUE, 
              draw.trees = TRUE, 
              tree.dat = tree.dat, 
              draw.particles = FALSE, 
              alpha.scale.cubes = TRUE,
              col=c(Commiphora_leptophloeos = "forestgreen"), 
              draw.axes = TRUE,
              which.axes = c(TRUE, TRUE, FALSE),
              which.axis.sides = c("-", "+", "--"),
              root.x.names = 1:5, 
              root.y.names = c("A", "B", "C", "D", "E"), 
              root.z.names = c("0 - 5 cm", "5 - 10 cm", "15 - 20 cm", "45 - 50 cm"),
              draw.legend = FALSE)
next3d()
# Commiphora_leptophloeos ind
my.cols <- RColorBrewer::brewer.pal(2,"Set3")
my.cols <- my.cols[1:2]
names(my.cols) <- ind.list$Commiphora_leptophloeos
plot_roots_3d(matches = best.matches, 
              ind.list = ind.list, 
              root.pos = root.pos, 
              taxa = "Commiphora_leptophloeos", 
              tax.type = "sp.ind", 
              draw.root.grid = TRUE, 
              draw.trees = TRUE, 
              tree.dat = tree.dat, 
              draw.particles = TRUE, 
              alpha.scale.cubes = FALSE,
              abundance.scale.type = "max",
              col = my.cols, 
              draw.axes = TRUE,
              which.axes = c(TRUE, TRUE, FALSE),
              which.axis.sides = c("-", "+", "--"),
              root.x.names = 1:5, 
              root.y.names = c("A", "B", "C", "D", "E"), 
              draw.legend = FALSE)
next3d()
# Cenostigma_microphyllum sp
plot_roots_3d(matches = best.matches, 
              ind.list = ind.list, 
              root.pos = root.pos, 
              taxa = "Cenostigma_microphyllum", 
              tax.type = "sp", 
              draw.root.grid = TRUE, 
              draw.trees = TRUE, 
              tree.dat = tree.dat, 
              draw.particles = FALSE, 
              alpha.scale.cubes = TRUE,
              col=c(Cenostigma_microphyllum = "forestgreen"), 
              draw.axes = TRUE,
              which.axes = c(TRUE, TRUE, FALSE),
              which.axis.sides = c("-", "+", "--"),
              root.x.names = 1:5, 
              root.y.names = c("A", "B", "C", "D", "E"), 
              root.z.names = c("0 - 5 cm", "5 - 10 cm", "15 - 20 cm", "45 - 50 cm"),
              draw.legend = FALSE)
next3d()
# Cenostigma_microphyllum ind
my.cols <- RColorBrewer::brewer.pal(length(ind.list$Cenostigma_microphyllum), "Set3")
my.cols <- my.cols[1:2]
names(my.cols) <- ind.list$Cenostigma_microphyllum
plot_roots_3d(matches = best.matches, 
              ind.list = ind.list, 
              root.pos = root.pos, 
              taxa = "Cenostigma_microphyllum", 
              tax.type = "sp.ind", 
              draw.root.grid = T, 
              draw.trees = TRUE, 
              tree.dat = tree.dat, 
              draw.particles = TRUE, 
              alpha.scale.cubes = FALSE,
              abundance.scale.type = "max",
              col = my.cols, 
              draw.axes = TRUE,
              which.axes = c(TRUE, TRUE, FALSE),
              which.axis.sides = c("-","+","--"),
              root.x.names = 1:5, 
              root.y.names = c("A", "B", "C", "D", "E"), 
              root.z.names = c("0 - 5 cm", "5 - 10 cm", "15 - 20 cm", "45 - 50 cm"),
              draw.legend = FALSE)
# manually adjust zoom and viewing angle 
rgl::rgl.snapshot("./results/3d_plots/Fig2.Cl.Cm.inds.png")
rgl::rgl.close()

## Supplementary figures
# species

# Calliandra_depauperata
rgl::open3d()
plot_roots_3d(matches = best.matches, 
              ind.list = ind.list, 
              root.pos = root.pos, 
              taxa = "Calliandra_depauperata",
              col = c(Calliandra_depauperata = "forestgreen"),
              tax.type = "sp", 
              draw.legend = FALSE,
              draw.particles = FALSE, 
              alpha.scale.cubes = TRUE,
              draw.trees = TRUE, 
              tree.dat = tree.dat, 
              draw.axes = TRUE, 
              draw.root.grid = TRUE,
              root.x.names = 1:5, 
              root.y.names = c("A", "B", "C", "D", "E"), 
              which.axes = c(TRUE, TRUE, FALSE)
)
rgl::rgl.snapshot("./results/3d_plots/Calliandra_depauperata.png")
rgl::rgl.close()

# Cereus_albicaulis
rgl::open3d()
plot_roots_3d(matches = best.matches, 
              ind.list = ind.list, 
              root.pos = root.pos, 
              taxa = "Cereus_albicaulis",
              col = c(Cereus_albicaulis = "forestgreen"),
              tax.type = "sp", 
              draw.legend = FALSE, 
              draw.particles = FALSE, 
              alpha.scale.cubes = TRUE,
              which.axes = c(TRUE, TRUE, FALSE),
              draw.trees = TRUE, 
              tree.dat = tree.dat, 
              draw.axes = TRUE, 
              draw.root.grid = TRUE,
              root.x.names = 1:5, 
              root.y.names = c("A", "B", "C", "D", "E") 
              
)
rgl::rgl.snapshot("./results/3d_plots/Cereus_albicaulis.png")
rgl::rgl.close()

# Chloroleucon_foliolosum
rgl::open3d()
plot_roots_3d(matches = best.matches, 
              ind.list = ind.list, 
              root.pos = root.pos, 
              taxa = "Chloroleucon_foliolosum",
              col = c(Chloroleucon_foliolosum = "forestgreen"),
              tax.type = "sp", 
              draw.legend = FALSE, 
              draw.particles = FALSE, 
              alpha.scale.cubes = TRUE,
              which.axes = c(TRUE, TRUE, FALSE),
              draw.trees = TRUE, 
              tree.dat = tree.dat, 
              draw.axes = TRUE, 
              draw.root.grid = TRUE,
              root.x.names = 1:5, 
              root.y.names = c("A", "B", "C", "D", "E"), 
              root.z.names = c("0 - 5 cm", "5 - 10 cm", "15 - 20 cm", "45 - 50 cm")
)
rgl::rgl.snapshot("./results/3d_plots/Chloroleucon_foliolosum.png")
rgl::rgl.close()

# species
# Cnidoscolus_quercifolius
rgl::open3d()
plot_roots_3d(matches = best.matches, 
              ind.list = ind.list, 
              root.pos = root.pos, 
              taxa = "Cnidoscolus_quercifolius",
              col = c(Cnidoscolus_quercifolius = "forestgreen"),
              tax.type = "sp", 
              draw.legend = FALSE, 
              draw.particles = FALSE, 
              alpha.scale.cubes = TRUE,
              which.axes = c(TRUE, TRUE, FALSE), 
              draw.trees = TRUE, 
              tree.dat = tree.dat, 
              draw.axes = TRUE, 
              draw.root.grid = TRUE,
              root.x.names = 1:5, 
              root.y.names = c("A", "B", "C", "D", "E"), 
              root.z.names = c("0 - 5 cm", "5 - 10 cm", "15 - 20 cm", "45 - 50 cm"), 
              which.axis.sides = c("+", "-", "++"),
              which.root.grid = c("x+", "y-", "z-")
              
)
rgl::rgl.snapshot("./results/3d_plots/Cnidoscolus_quercifolius.png")
rgl::rgl.close()

# species
# Commiphora_leptophloeos
rgl::open3d()
plot_roots_3d(matches = best.matches, 
              ind.list = ind.list, 
              root.pos = root.pos, 
              taxa = "Commiphora_leptophloeos",
              col = c(Commiphora_leptophloeos = "forestgreen"),
              tax.type = "sp", 
              draw.legend = FALSE, 
              max.n.particles = 500, 
              draw.trees = TRUE, 
              tree.dat = tree.dat, 
              draw.axes = TRUE, 
              draw.root.grid = TRUE,
              root.x.names = 1:5, 
              root.y.names = c("A", "B", "C", "D", "E"), 
              root.z.names = c("0 - 5 cm", "5 - 10 cm", "15 - 20 cm", "45 - 50 cm")
)
rgl::rgl.snapshot("./results/3d_plots/Commiphora_leptophloeos.png")
rgl::rgl.close()

# species
# Croton_echioides
rgl::open3d()
plot_roots_3d(matches = best.matches, 
              ind.list = ind.list, 
              root.pos = root.pos, 
              taxa = "Croton_echioides",
              col = c(Croton_echioides = "forestgreen"),
              tax.type = "sp", 
              draw.legend = FALSE, 
              draw.particles = FALSE, 
              alpha.scale.cubes = TRUE,
              which.axes = c(TRUE, TRUE, FALSE), 
              draw.trees = TRUE, 
              tree.dat = tree.dat, 
              draw.axes = TRUE, 
              draw.root.grid = TRUE,
              root.x.names = 1:5, 
              root.y.names = c("A", "B", "C", "D", "E"), 
              root.z.names = c("0 - 5 cm", "5 - 10 cm", "15 - 20 cm", "45 - 50 cm")
)
rgl::rgl.snapshot("./results/3d_plots/Croton_echioides.png")
rgl::rgl.close()

# species
# Ditaxis_desertorum
rgl::open3d()
plot_roots_3d(matches = best.matches, 
              ind.list = ind.list, 
              root.pos = root.pos, 
              taxa = "Ditaxis_desertorum",
              col = c(Ditaxis_desertorum = "forestgreen"),
              tax.type = "sp", 
              draw.legend = FALSE, 
              draw.particles = FALSE, 
              alpha.scale.cubes = TRUE,
              which.axes = c(TRUE, TRUE, FALSE), 
              draw.trees = TRUE, 
              tree.dat = tree.dat, 
              draw.axes = TRUE, 
              draw.root.grid = TRUE,
              root.x.names = 1:5, 
              root.y.names = c("A", "B", "C", "D", "E"), 
              root.z.names = c("0 - 5 cm", "5 - 10 cm", "15 - 20 cm", "45 - 50 cm")
)
rgl::rgl.snapshot("./results/3d_plots/Ditaxis_desertorum.png")
rgl::rgl.close()

# species
# Handroanthus_spongiosus
rgl::open3d()
plot_roots_3d(matches = best.matches, 
              ind.list = ind.list, 
              root.pos = root.pos, 
              taxa = "Handroanthus_spongiosus",
              col = c(Handroanthus_spongiosus = "forestgreen"),
              tax.type = "sp", 
              draw.legend = FALSE, 
              draw.particles = FALSE, 
              alpha.scale.cubes = TRUE,
              which.axes = c(TRUE, TRUE, FALSE), 
              draw.trees = TRUE, 
              tree.dat = tree.dat, 
              draw.axes = TRUE, 
              draw.root.grid = TRUE,
              root.x.names = 1:5, 
              root.y.names = c("A", "B", "C", "D", "E"), 
              root.z.names = c("0 - 5 cm", "5 - 10 cm", "15 - 20 cm", "45 - 50 cm")
)
rgl::rgl.snapshot("./results/3d_plots/Handroanthus_spongiosus.png")
rgl::rgl.close()

# species
# Jatropha_mollissima
rgl::open3d()
plot_roots_3d(matches = best.matches, 
              ind.list = ind.list, 
              root.pos = root.pos, 
              taxa = "Jatropha_mollissima",
              col = c(Jatropha_mollissima = "forestgreen"),
              tax.type = "sp", 
              draw.legend = FALSE, 
              draw.particles = FALSE, 
              alpha.scale.cubes = TRUE,
              which.axes = c(TRUE, TRUE, FALSE), 
              draw.trees = TRUE, 
              tree.dat = tree.dat, 
              draw.axes = TRUE, 
              draw.root.grid = TRUE,
              root.x.names = 1:5, 
              root.y.names = c("A", "B", "C", "D", "E"), 
              root.z.names = c("0 - 5 cm", "5 - 10 cm", "15 - 20 cm", "45 - 50 cm")
)
rgl::rgl.snapshot("./results/3d_plots/Jatropha_mollissima.png")
rgl::rgl.close()

# species
# Manihot_carthagenensis
rgl::open3d()
plot_roots_3d(matches = best.matches, 
              ind.list = ind.list, 
              root.pos = root.pos, 
              taxa = "Manihot_carthagenensis",
              col = c(Manihot_carthagenensis = "forestgreen"),
              tax.type = "sp", 
              draw.legend = FALSE, 
              draw.particles = FALSE, 
              alpha.scale.cubes = TRUE,
              which.axes = c(TRUE, TRUE, FALSE), 
              draw.trees = TRUE, 
              tree.dat = tree.dat, 
              draw.axes = TRUE, 
              draw.root.grid = TRUE,
              root.x.names = 1:5, 
              root.y.names = c("A", "B", "C", "D", "E"), 
              root.z.names = c("0 - 5 cm", "5 - 10 cm", "15 - 20 cm", "45 - 50 cm")
)
rgl::rgl.snapshot("./results/3d_plots/Manihot_carthagenensis.png")
rgl::rgl.close()

# species
# Mimosa_arenosa
rgl::open3d()
plot_roots_3d(matches = best.matches, 
              ind.list = ind.list, 
              root.pos = root.pos, 
              taxa = "Mimosa_arenosa",
              col = c(Mimosa_arenosa="forestgreen"),
              tax.type = "sp", 
              draw.legend = FALSE, 
              draw.particles = FALSE, 
              alpha.scale.cubes = TRUE,
              which.axes = c(TRUE, TRUE, FALSE), 
              draw.trees = TRUE, 
              tree.dat = tree.dat, 
              draw.axes = TRUE, 
              draw.root.grid = TRUE,
              root.x.names = 1:5, 
              root.y.names = c("A", "B", "C", "D", "E"), 
              root.z.names = c("0 - 5 cm", "5 - 10 cm", "15 - 20 cm", "45 - 50 cm")
)
rgl::rgl.snapshot("./results/3d_plots/Mimosa_arenosa.png")
rgl::rgl.close()


# species
# Neoglaziovia_variegata
rgl::open3d()
plot_roots_3d(matches = best.matches, 
              ind.list = ind.list, 
              root.pos = root.pos, 
              taxa = "Neoglaziovia_variegata",
              col = c(Neoglaziovia_variegata = "forestgreen"),
              tax.type = "sp", 
              draw.legend = FALSE, 
              draw.particles = FALSE, 
              alpha.scale.cubes = TRUE,
              which.axes = c(TRUE, TRUE, FALSE), 
              draw.trees = TRUE, 
              tree.dat = tree.dat, 
              draw.axes = TRUE, 
              draw.root.grid = TRUE,
              root.x.names = 1:5, 
              root.y.names = c("A", "B", "C", "D", "E"), 
              root.z.names = c("0 - 5 cm", "5 - 10 cm", "15 - 20 cm", "45 - 50 cm")
)
rgl::rgl.snapshot("./results/3d_plots/Neoglaziovia_variegata.png")
rgl::rgl.close()

# species
# Pseudobombax_simplicifolium
rgl::open3d()
plot_roots_3d(matches = best.matches, 
              ind.list = ind.list, 
              root.pos = root.pos, 
              taxa = "Pseudobombax_simplicifolium",
              col = c(Pseudobombax_simplicifolium = "forestgreen"),
              tax.type = "sp", 
              draw.legend = FALSE, 
              draw.particles = FALSE, 
              alpha.scale.cubes = TRUE,
              which.axes = c(TRUE, TRUE, FALSE), 
              draw.trees = TRUE, 
              tree.dat = tree.dat, 
              draw.axes = TRUE, 
              draw.root.grid = TRUE,
              root.x.names = 1:5, 
              root.y.names = c("A", "B", "C", "D", "E"), 
              root.z.names = c("0 - 5 cm", "5 - 10 cm", "15 - 20 cm", "45 - 50 cm")
)
rgl::rgl.snapshot("./results/3d_plots/Pseudobombax_simplicifolium.png")
rgl::rgl.close()



# species
# Sapium_glandulosum
rgl::open3d()
plot_roots_3d(matches = best.matches, 
              ind.list = ind.list, 
              root.pos = root.pos, 
              taxa = "Sapium_glandulosum",
              col = c(Sapium_glandulosum = "forestgreen"),
              tax.type = "sp", 
              draw.legend = FALSE, 
              draw.particles = FALSE, 
              alpha.scale.cubes = TRUE,
              which.axes = c(TRUE, TRUE, FALSE), 
              draw.trees = TRUE, 
              tree.dat = tree.dat, 
              draw.axes = TRUE, 
              draw.root.grid = TRUE,
              root.x.names = 1:5, 
              root.y.names = c("A", "B", "C", "D", "E"), 
              root.z.names = c("0 - 5 cm", "5 - 10 cm", "15 - 20 cm", "45 - 50 cm")
)
rgl::rgl.snapshot("./results/3d_plots/Sapium_glandulosum.png")
rgl::rgl.close()

# species
# Schinopsis_brasiliensis
rgl::open3d()
plot_roots_3d(matches = best.matches, 
              ind.list = ind.list, 
              root.pos = root.pos, 
              taxa = "Schinopsis_brasiliensis",
              col = c(Schinopsis_brasiliensis = "forestgreen"),
              tax.type = "sp", 
              draw.legend = FALSE, 
              draw.particles = FALSE, 
              alpha.scale.cubes = TRUE,
              which.axes = c(TRUE, TRUE, FALSE), 
              draw.trees = TRUE, 
              tree.dat = tree.dat, 
              draw.axes = TRUE, 
              draw.root.grid = TRUE,
              root.x.names = 1:5, 
              root.y.names = c("A", "B", "C", "D", "E"), 
              root.z.names = c("0 - 5 cm", "5 - 10 cm", "15 - 20 cm", "45 - 50 cm")
)
rgl::rgl.snapshot("./results/3d_plots/Schinopsis_brasiliensis.png")
rgl::rgl.close()

# species
# Senna_macranthera
rgl::open3d()
plot_roots_3d(matches = best.matches, 
              ind.list = ind.list, 
              root.pos = root.pos, 
              taxa = "Senna_macranthera",
              col = c(Senna_macranthera = "forestgreen"),
              tax.type = "sp", 
              draw.legend = FALSE, 
              draw.particles = FALSE, 
              alpha.scale.cubes = TRUE,
              which.axes = c(TRUE, TRUE, FALSE), 
              draw.trees = TRUE, 
              tree.dat = tree.dat, 
              draw.axes = TRUE, 
              draw.root.grid = TRUE,
              root.x.names = 1:5, 
              root.y.names = c("A", "B", "C", "D", "E"), 
              root.z.names = c("0 - 5 cm", "5 - 10 cm", "15 - 20 cm", "45 - 50 cm")
)
rgl::rgl.snapshot("./results/3d_plots/Senna_macranthera.png")
rgl::rgl.close()

# species
# Tacinga_inamoena
rgl::open3d()
plot_roots_3d(matches = best.matches, 
              ind.list = ind.list, 
              root.pos = root.pos, 
              taxa = "Tacinga_inamoena",
              col = c(Tacinga_inamoena = "forestgreen"),
              tax.type = "sp", 
              draw.legend = FALSE, 
              draw.particles = FALSE, 
              alpha.scale.cubes = TRUE,
              which.axes = c(TRUE, TRUE, FALSE), 
              draw.trees = TRUE, 
              tree.dat = tree.dat, 
              draw.axes = TRUE, 
              draw.root.grid = TRUE,
              root.x.names = 1:5, 
              root.y.names = c("A", "B", "C", "D", "E"), 
              root.z.names = c("0 - 5 cm", "5 - 10 cm", "15 - 20 cm", "45 - 50 cm")
)
rgl::rgl.snapshot("./results/3d_plots/Tacinga_inamoena.png")
rgl::rgl.close()

# species
# Varronia_leucocephala
rgl::open3d()
plot_roots_3d(matches = best.matches, 
              ind.list = ind.list, 
              root.pos = root.pos, 
              taxa = "Varronia_leucocephala",
              col = c(Varronia_leucocephala = "forestgreen"),
              tax.type = "sp", 
              draw.legend = FALSE, 
              draw.particles = FALSE, 
              alpha.scale.cubes = TRUE,
              which.axes = c(TRUE, TRUE, FALSE), 
              draw.trees = TRUE, 
              tree.dat = tree.dat, 
              draw.axes = TRUE, 
              draw.root.grid = TRUE,
              root.x.names = 1:5, 
              root.y.names = c("A", "B", "C", "D", "E"), 
              root.z.names = c("0 - 5 cm", "5 - 10 cm", "15 - 20 cm", "45 - 50 cm")
)
rgl::rgl.snapshot("./results/3d_plots/Varronia_leucocephala.png")
rgl::rgl.close()

  # individuals
  # Cenostigma_microphyllum
  my.cols <- brewer.pal(2, "Set3")
  my.cols <- my.cols[1:2]
  names(my.cols) <- ind.list$Cenostigma_microphyllum
  rgl::open3d()
  plot_roots_3d(matches = best.matches, 
                ind.list = ind.list, 
                root.pos = root.pos, 
                taxa = "Cenostigma_microphyllum",
                col = my.cols,
                tax.type = "sp.ind", 
                draw.legend = FALSE, 
                max.n.particles = 500, 
                draw.trees = TRUE, 
                tree.dat = tree.dat, 
                draw.axes = TRUE, 
                draw.root.grid = TRUE,
                root.x.names = 1:5, 
                root.y.names = c("A", "B", "C", "D", "E"), 
                root.z.names = c("0 - 5 cm", "5 - 10 cm", "15 - 20 cm", "45 - 50 cm")
  )
  rgl::rgl.snapshot("./results/3d_plots/Cenostigma_microphyllum.inds.png")
  rgl::rgl.close()

# individuals

open3d()
my.cols <- RColorBrewer::brewer.pal(length(ind.list$Cnidoscolus_quercifolius),"Set3")
my.cols <- my.cols[1:length(ind.list$Cnidoscolus_quercifolius)]
names(my.cols) <- ind.list$Cnidoscolus_quercifolius
plot_roots_3d(matches = best.matches, 
              ind.list = ind.list, 
              root.pos = root.pos, 
              taxa = "Cnidoscolus_quercifolius", 
              tax.type = "sp.ind", 
              draw.root.grid = TRUE, 
              draw.trees = TRUE, 
              tree.dat = tree.dat, 
              draw.particles = TRUE,
              which.axes = c(TRUE, TRUE, FALSE),
              alpha.scale.cubes = FALSE,
              abundance.scale.type = "max",
              col = my.cols, 
              draw.axes = TRUE,
              which.axis.sides = c("+", "-", "++"),
              root.x.names = 1:5, 
              root.y.names = c("A", "B", "C", "D", "E"), 
              root.z.names = c("0 - 5 cm", "5 - 10 cm", "15 - 20 cm", "45 - 50 cm"),
              draw.legend = FALSE)
rgl::rgl.snapshot("./results/3d_plots/Cnidoscolus_quercifolius.ind.png")
rgl::rgl.close()
#Croton_echioides
open3d()
my.cols <- RColorBrewer::brewer.pal(length(ind.list$Croton_echioides), "Set3")
my.cols <- my.cols[1:length(ind.list$Croton_echioides)]
names(my.cols) <- ind.list$Croton_echioides
plot_roots_3d(matches = best.matches, 
              ind.list = ind.list, 
              root.pos = root.pos, 
              taxa = "Croton_echioides", 
              tax.type = "sp.ind", 
              draw.root.grid = TRUE, 
              draw.trees = TRUE, 
              tree.dat = tree.dat, 
              draw.particles = TRUE,
              which.axes = c(TRUE, TRUE, FALSE),
              alpha.scale.cubes = FALSE,
              abundance.scale.type = "max",
              col = my.cols, 
              draw.axes = TRUE,
              which.axis.sides = c("-", "+", "--"),
              root.x.names = 1:5, 
              root.y.names = c("A", "B", "C", "D", "E"), 
              root.z.names = c("0 - 5 cm", "5 - 10 cm", "15 - 20 cm", "45 - 50 cm"),
              draw.legend = FALSE)
rgl::rgl.snapshot("./results/3d_plots/Croton_echioides.ind.png")
rgl::rgl.close()

# Handroanthus_spongiosus
open3d()
my.cols <- RColorBrewer::brewer.pal(length(ind.list$Handroanthus_spongiosus),"Set3")
my.cols <- my.cols[1:length(ind.list$Handroanthus_spongiosus)]
names(my.cols) <- ind.list$Handroanthus_spongiosus
plot_roots_3d(matches = best.matches, 
              ind.list = ind.list, 
              root.pos = root.pos, 
              taxa = "Handroanthus_spongiosus", 
              tax.type = "sp.ind", 
              draw.root.grid = TRUE, 
              draw.trees = TRUE, 
              tree.dat = tree.dat, 
              draw.particles = TRUE,
              which.axes = c(TRUE, TRUE, FALSE),
              alpha.scale.cubes = FALSE,
              abundance.scale.type = "max",
              col = my.cols, 
              draw.axes = TRUE,
              which.axis.sides = c("-", "+", "--"),
              root.x.names = 1:5, 
              root.y.names = c("A", "B", "C", "D", "E"), 
              root.z.names = c("0 - 5 cm", "5 - 10 cm", "15 - 20 cm", "45 - 50 cm"),
              draw.legend = FALSE)
rgl::rgl.snapshot("./results/3d_plots/Handroanthus_spongiosus.ind.png")
rgl::rgl.close()

# Jatropha_mollissima
open3d()
my.cols <- RColorBrewer::brewer.pal(length(ind.list$Jatropha_mollissima),"Set3")
my.cols <- my.cols[1:length(ind.list$Jatropha_mollissima)]
names(my.cols) <- ind.list$Jatropha_mollissima
# remove those with no roots detected
my.cols <- my.cols[c("L_1","L_10","L_19","L_23","L_41")]
plot_roots_3d(matches = best.matches, 
              ind.list = ind.list, 
              root.pos = root.pos, 
              taxa = "Jatropha_mollissima", 
              tax.type = "sp.ind", 
              draw.root.grid = TRUE, 
              draw.trees = TRUE, 
              tree.dat = tree.dat, 
              draw.particles = TRUE,
              which.axes = c(TRUE, TRUE, FALSE),
              alpha.scale.cubes = F,
              abundance.scale.type = "max",
              col = my.cols, 
              draw.axes = TRUE,
              which.axis.sides = c("-", "+", "--"),
              root.x.names = 1:5, 
              root.y.names = c("A", "B", "C", "D", "E"), 
              root.z.names = c("0 - 5 cm", "5 - 10 cm", "15 - 20 cm", "45 - 50 cm"),
              draw.legend = FALSE)
rgl::rgl.snapshot("./results/3d_plots/Jatropha_mollissima.ind.png")
rgl::rgl.close()

# Pseudobombax_simplicifolium
open3d()
my.cols <- RColorBrewer::brewer.pal(length(ind.list$Pseudobombax_simplicifolium),"Set3")
my.cols <- my.cols[1:length(ind.list$Pseudobombax_simplicifolium)]
names(my.cols) <- ind.list$Pseudobombax_simplicifolium
plot_roots_3d(matches = best.matches, 
              ind.list = ind.list, 
              root.pos = root.pos, 
              taxa = "Pseudobombax_simplicifolium", 
              tax.type = "sp.ind", 
              draw.root.grid = TRUE, 
              draw.trees = TRUE, 
              tree.dat = tree.dat, 
              draw.particles = TRUE,
              which.axes = c(TRUE, TRUE, FALSE),
              alpha.scale.cubes = FALSE,
              abundance.scale.type = "max",
              col = my.cols, 
              draw.axes = TRUE,
              which.axis.sides = c("-", "+", "--"),
              root.x.names = 1:5, 
              root.y.names = c("A", "B", "C", "D", "E"), 
              root.z.names = c("0 - 5 cm", "5 - 10 cm", "15 - 20 cm", "45 - 50 cm"),
              draw.legend = FALSE)
rgl::rgl.snapshot("./results/3d_plots/Pseudobombax_simplicifolium.ind.png")
rgl::rgl.close()

# Sapium_glandulosum
open3d()
my.cols <- RColorBrewer::brewer.pal(length(ind.list$Sapium_glandulosum),"Set3")
my.cols <- my.cols[1:length(ind.list$Sapium_glandulosum)]
names(my.cols) <- ind.list$Sapium_glandulosum
plot_roots_3d(matches = best.matches, 
              ind.list = ind.list, 
              root.pos = root.pos, 
              taxa = "Sapium_glandulosum", 
              tax.type = "sp.ind", 
              draw.root.grid = TRUE, 
              draw.trees = TRUE, 
              tree.dat = tree.dat, 
              draw.particles = TRUE,
              which.axes = c(TRUE, TRUE, FALSE),
              alpha.scale.cubes = FALSE,
              abundance.scale.type = "max",
              col = my.cols, 
              draw.axes = TRUE,
              which.axis.sides = c("-", "+", "--"),
              root.x.names = 1:5, 
              root.y.names = c("A", "B", "C", "D", "E"), 
              root.z.names = c("0 - 5 cm", "5 - 10 cm", "15 - 20 cm", "45 - 50 cm"),
              draw.legend = FALSE)
rgl::rgl.snapshot("./results/3d_plots/Sapium_glandulosum.ind.png")
rgl::rgl.close()

# Schinopsis_brasiliensis
open3d()
my.cols <- RColorBrewer::brewer.pal(length(ind.list$Schinopsis_brasiliensis),"Set3")
my.cols <- my.cols[1:length(ind.list$Schinopsis_brasiliensis)]
names(my.cols) <- ind.list$Schinopsis_brasiliensis
plot_roots_3d(matches = best.matches, 
              ind.list = ind.list, 
              root.pos = root.pos, 
              taxa = "Schinopsis_brasiliensis", 
              tax.type = "sp.ind", 
              draw.root.grid = TRUE, 
              draw.trees = TRUE, 
              tree.dat = tree.dat, 
              draw.particles = TRUE,
              which.axes = c(TRUE, TRUE, FALSE),
              alpha.scale.cubes = FALSE,
              abundance.scale.type = "max",
              col = my.cols, 
              draw.axes = TRUE,
              which.axis.sides = c("-", "+", "--"),
              root.x.names = 1:5, 
              root.y.names = c("A", "B", "C", "D", "E"), 
              root.z.names = c("0 - 5 cm", "5 - 10 cm", "15 - 20 cm", "45 - 50 cm"),
              draw.legend = FALSE)
rgl::rgl.snapshot("./results/3d_plots/Schinopsis_brasiliensis.ind.png")
rgl::rgl.close()
