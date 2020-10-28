# ~ ~ ~ ~#
# Part C #
# ~ ~ ~ ~#

# Load R-packages
library(BacArena)
library(data.table)
sybil::SYBIL_SETTINGS("SOLVER","cplexAPI"); ok <- 1

# Load reconstructed models
er <- readRDS("eure.RDS") # E. rectale
bl <- readRDS("bilo.RDS") # B. longum

# Small fix to D/L-Lactate secretion *
bl <- rmReact(bl, react = "EX_cpd00221_e0")

# Construct the organism objects for BacArena simulations
eure <- Bac(er)
bilo <- Bac(bl)

# Construct the arena size 10x10 grid cells
arena <- Arena(n = 10, m = 10)

# For each organism, populate randomly 2 grid cells in the Arena as 
# 'starter culture'
arena <- addOrg(arena, eure, amount = 2)
arena <- addOrg(arena, bilo, amount = 2)

# add substrates to arena
arena_subs <- fread("gf_medium.csv") # same as gapfill medium
arena_subs[, ex.rxn := paste0("EX_", compounds, "_e0")]

arena <- addSubs(arena, smax = arena_subs$maxFlux, 
                 mediac = arena_subs$ex.rxn, unit = "mM", addAnyway = T)
# Remove acetate from initial substrate list to see effect of Cross-Feeding
arena <- rmSubs(arena, mediac = "EX_cpd00029_e0") 

# ~ ~ ~ ~#
# Part B #
# ~ ~ ~ ~#

# Simulation for 20 time steps
CF_sim <- simEnv(arena,time=20, sec_obj = "mtf")

# Plot levels of Acetate, Buyrate, and Lactate as well as growth
#png("CF_eure_bilo.png", width = 600, height = 320)
par(mfrow=c(1,2))
plotCurves2(CF_sim,legendpos = "topleft",
            subs = c("cpd00211_e0","cpd00029_e0","cpd00159_e0"),
            dict = list(cpd00211_e0 = "Butyrate", 
                        cpd00029_e0 = "Acetate", 
                        cpd00159_e0 = "Lactate"))
#dev.off()

# ~ ~ ~ ~#
# Part C #
# ~ ~ ~ ~#

# Lets get the exchange fluxs for each grid cell at time step 20
dt.cf <- CF_sim@exchangeslist[[20]]

# Limit output to Acetate, Butyrate, and Lactate
dt.cf <- as.data.table(dt.cf[,c("species",
                                "EX_cpd00029_e0",
                                "EX_cpd00211_e0",
                                "EX_cpd00159_e0")])

# Rename column names (this is just aestetics)
dt.cf <- dt.cf[,.(species, 
                  Acetate = EX_cpd00029_e0, 
                  Butyrate = EX_cpd00211_e0,
                  Lactate = EX_cpd00159_e0)]

# Wide-To-Long table transformation
dt.cf <- melt(dt.cf, 
              id.vars = "species", 
              variable.name = "Metabolite", 
              value.name = "Flux")
dt.cf <- dt.cf[!is.na(Flux)] # rm NA entries (no exchange reaction in grid cell)

# Sum exhanges for each species and metabolite over all 100 grid cells
dt.cf[, .(summed.Flux = sum(Flux)), by = c("species","Metabolite")]

