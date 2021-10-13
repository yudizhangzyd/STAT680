grid <- readRDS("grid.rds")
k_result <- lapply(1:nrow(grid), function(idx) {
  SynClustR::kcdf(grid$data_result[[idx]]$tt$n_rs, xgrid = grid$nrstbl_result[[idx]]$n_rs_cluster)
})
# k_result <- SynClustR::kcdf(all_container[[1]]$data$n_rs, xgrid = all_container[[1]]$misclass$n_rs_cluster)
saveRDS(k_result, "k_result.rds")
