devtools::clean_dll()
devtools::load_all()


saveRDS(
    spatialEnhance(
        readRDS("__tools/data.RDS"),
        q = 4,
        platform = "ST",
        d = 7,
        model = "t",
        gamma = 2,
        jitter_prior = 0.3,
        jitter_scale = 3.5,
        nrep = 1000,
        burn.in = 100,
        save.chain = TRUE
    ),
    "__tools/results.RDS"
)
