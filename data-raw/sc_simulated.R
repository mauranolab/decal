## code to prepare `sc_simulated` dataset goes here
set.seed(199155)
dp <- ceiling(rlnorm(500, 9, .5))
ps <- seq_log(0.0001, 0.01, 1500)

sc_simulated <- as(sim_count(dp, ps, theta = 80), "dgTMatrix")
rownames(sc_simulated) <- sprintf("GENE%05d", seq_along(ps))
colnames(sc_simulated) <- sprintf("CELL%03d", seq_along(dp))

usethis::use_data(sc_simulated, overwrite = TRUE)
