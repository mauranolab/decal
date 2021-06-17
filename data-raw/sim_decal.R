## code to prepare `sc_simulated` dataset goes here
set.seed(919564)
devtools::load_all()

nclone <- 25
ncells <- 1000
ngenes <- 1500
lfc <- rep(c(0, 1, 2, 3, -1, -2, -3), c(80, 5, 5, 5, 5, 5, 5))

depth <- floor(rlnorm(ncells, 9, .5))
names(depth) <- sprintf("cell%03d", seq_along(depth))
ratio <- 10 ** seq(-4, -2, length.out = ngenes)
names(ratio) <- sprintf("gene%03d", seq_along(ratio))

sim_decal <- sim_experiment(
  depth, ratio, lfc, nclones = nclone, min_n = 5, max_n = 30, theta = 100L
)
## replace indexes by names
clone_name <- sprintf("clone-%03d", seq_len(nclone))
names(sim_decal$clone) <- clone_name
sim_decal$perturbations$gene  <- names(ratio)[sim_decal$perturbations$gene]
sim_decal$perturbations$clone <- clone_name[sim_decal$perturbations$clone]

sim_decal$clone <- lapply(sim_decal$clone, function(x) names(depth)[x])

usethis::use_data(sim_decal, overwrite = TRUE)
