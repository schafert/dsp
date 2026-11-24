# Stress test for the time-varying regression sampler (JSS review, item 2).
# Fits the CAPM example of the article with the reviewer's MCMC settings
# for one seed and fails if the fit errors or does not finish in time.
# Usage: Rscript capm_stress.R <seed>
library(dsp)
seed <- as.integer(commandArgs(TRUE)[1])
df_capm <- readRDS(file.path(".github", "stress", "df_capm.rds"))
spec_reg <- dsp_spec(family = "gaussian", model = "regression",
                     X = as.matrix(df_capm$x))
set.seed(seed)
t0 <- Sys.time()
fit <- dsp_fit(y = df_capm$y, model_spec = spec_reg,
               nsave = 5000, nburn = 20000, verbose = FALSE)
cat(sprintf("seed %d: %.1f min, QR fallbacks = %d\n", seed,
            as.numeric(difftime(Sys.time(), t0, units = "mins")),
            get(".dsp_state", asNamespace("dsp"))$n_sqrt))
stopifnot(all(is.finite(as.mcmc(fit, pars = "beta"))))
