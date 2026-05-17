### the following are bad functions, i'm aware
### the main idea is that this will find an arithmetic solution to a mean + sd
### closer to what we think is the true mean + sd of the core
### found using a solution to a set of linear equations

arithmetic_mean_adjust <- function(coremean, coresd, top_z, water_val) {
  beta_x <- (((coremean ** 3) - (coremean ** 2 * coresd * top_z) + (coremean * water_val * coresd * top_z)
  ) /
    ((coremean ** 2 * coresd * top_z) + (coresd ** 3 * top_z ** 3))) + 1 - (coremean / (coresd * top_z))
  adjust_mean <- coremean * beta_x
}

arithmetic_sd_adjust <- function(coremean, coresd, top_z, water_val) {
  beta_s <- ((coremean ** 2) - (coremean * coresd * top_z) + (water_val * coresd * top_z)) /
    ((coremean ** 2) + (coresd ** 2 * top_z ** 2))
  adjust_sd <- coresd * beta_s
}

naive_mean_sd_adjust <- function(coremean, coresd, err) {
  corecv <- coresd / coremean
  naive_mean <- coremean + err
  naive_sd <- corecv * naive_mean
  naive_mean_sd <- c(naive_mean, naive_sd)
}
