# Cross-feeding of two gastrointestinal bacteria

##### Background

The intestinal bacterium *Eubacterium rectale* is known to be able to use acetate as energy source under anaerobic conditions and thereby forms butyrate as end product ([Rivère *et al.* (2015) Appl Envrion Microbiol](https://pubmed.ncbi.nlm.nih.gov/26319874/)). Acetate is a common fermentation end product in a number of different other intestinal bacteria, including Bifidobacteria (e.g. *Bifidobacterium longum*). In this tutorial, genome-scale models for *E. rectale* and *B. longum* are reconstructed using **gapseq**. Subsequently, the two models are simulated in co-growth and their interaction is investigated.

*NOTE: All intermediate files produced by the commands below are stored at the github repository (https://github.com/Waschina/gapseq.tutorial.data), which you could download/clone if you wish to start not at the beginning but at a later step of this tutorial.*

##### Input

- Genome assemblies:

  - *Eubacterium rectale* ATCC 33656

    RefSeq: `GCF_000020605.1`

  - *Bifidobacterium longum* NCC2705: 

    RefSeq: `GCF_000007525.1`

- Growth media file: `gf_medium.csv` 

  This is basically a glucose and acetate minimal medium. No amino acids are added since both organisms are likely prototrophic for all proteinogenic amino acids, based on predictions using [GapMind](http://papers.genomics.lbl.gov/cgi-bin/gapView.cgi) ([Price et al. (2019) mSystems](https://doi.org/10.1101/741918 )).  

  E. rectale: [View Gapmind results](http://papers.genomics.lbl.gov/cgi-bin/gapView.cgi?orgs=NCBI__GCF_000020605.1&set=aa); B. longum: [View Gapmind results](http://papers.genomics.lbl.gov/cgi-bin/gapView.cgi?orgs=NCBI__GCF_000007525.1&set=aa)



##### Preparations

Download genome assemblies and gapfill medium. Renaming files.

```sh
#!/bin/bash

# Download genome assemblies 
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/020/605/GCF_000020605.1_ASM2060v1/GCF_000020605.1_ASM2060v1_genomic.fna.gz .
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/525/GCF_000007525.1_ASM752v1/GCF_000007525.1_ASM752v1_genomic.fna.gz .

# Download gapfill-medium file
wget https://github.com/Waschina/gapseq.tutorial.data/raw/master/CF_eure_bilo/gf_medium.csv .

# Rename genomes to "eure" (E. rectale) and "bilo" (B. longum) 
mv GCF_000020605.1_ASM2060v1_genomic.fna.gz eure.fna.gz
mv GCF_000007525.1_ASM752v1_genomic.fna.gz bilo.fna.gz
```



##### gapseq reconstruction 

Now we have the genome sequences and a gapfill medium. That is all we need. Lets reconstruct models:

```sh
#!/bin/bash

modelA="eure"
modelB="bilo"

# If not set already (e.g via .bashrc): set the path to gapseq
# There are different ways to do this. One example:
gapseq=~/workspace/2018/gapseq/./gapseq

# Reaction & Pathway prediction
$gapseq find -p all -b 200 $modelA.fna.gz
$gapseq find -p all -b 200 $modelB.fna.gz

# Transporter prediction
$gapseq find-transport -b 200 $modelA.fna.gz 
$gapseq find-transport -b 200 $modelB.fna.gz

# Building Draft Model - based on Reaction-, Pathway-, and Transporter prediction
$gapseq draft -r $modelA-all-Reactions.tbl -t $modelA-Transporter.tbl -p $modelA-all-Pathways.tbl -c $modelA.fna.gz -u 200 -l 100
$gapseq draft -r $modelB-all-Reactions.tbl -t $modelB-Transporter.tbl -p $modelB-all-Pathways.tbl -c $modelB.fna.gz -u 200 -l 100

# Gapfilling
$gapseq fill -m $modelA-draft.RDS -n gf_medium.csv -c $modelA-rxnWeights.RDS -g $modelA-rxnXgenes.RDS -b 100
$gapseq fill -m $modelB-draft.RDS -n gf_medium.csv -c $modelB-rxnWeights.RDS -g $modelB-rxnXgenes.RDS -b 100
```

The final models are stored as R-Object files: `eure.RDS` and `bilo.RDS`, which can be loaded in R using the `readRDS()` command (After the *sybil* package has bee loaded). 

NOTE: If you repeat the gapfilling step, gapseq will <u>not</u> overwrite the final model files. Instead it saves them as `eure-gapfilled.RDS` and/or `bilo-gapfilled.RDS`. However, if you repeat the gapfilling step for the third time or even more, it will overwrite the `...-gapfilled.RDS` files (a behavior that we will fix).



##### Community simulation

Here, we will use the R-Package `BacArena` to perform an agent-based simulation for the co-growth of *B. longum* and *E. rectale*. The following code-block shows the R-source code for a simple community metabolism simulation.

```R
# Load R-pakages
library(BacArena)
library(data.table)

# Load reconstructed models
er <- readRDS("glycan_cf/eure.RDS") # E. rectale
bl <- readRDS("glycan_cf/bilo.RDS") # B. longum

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

```

( * *gapseq frequently predicts, that if the organism is producing lactate as fermentation end product,  the optimal solution could involve both  enantiomers: D-/L-Lactate. For plotting & analysis reasons we prohibit the production of D-Lactate to ensure that we see the level of produced Lactate only as one metabolite: L-Lactate. This has otherwise no effect on simulation results.*)

Now we are ready to perform the actual community simulation and plot results:

```R
# Simulation for 20 time steps
CF_sim <- simEnv(arena,time=20, sec_obj = "mtf")

# Plot levels of Acetate, Buyrate, and Lactate as well as growth
par(mfrow=c(1,2))
plotCurves2(CF_sim,legendpos = "topleft",
            subs = c("cpd00211_e0","cpd00029_e0","cpd00159_e0"),
            dict = list(cpd00211_e0 = "Butyrate", 
                        cpd00029_e0 = "Acetate", 
                        cpd00159_e0 = "Lactate"))
```

![](https://github.com/Waschina/gapseq.tutorial.data/raw/master/CF_eure_bilo/CF_eure_bilo.svg)

The simulations predicted, that acetate, butyrate, and lactate are produced during co-growth of *E. rectale* and *B. longum*.

Next, let's see, if some of the fermentation products are partially consumed by one of the organisms. This involves a few lines of data wrangling:

```R
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
```

The output:

```
   species Metabolite   summed.Flux
1:    bilo    Acetate  1.621512e+03
2:    eure    Acetate -9.676227e+02
3:    bilo   Butyrate -8.231749e-13
4:    eure   Butyrate  1.900784e+03
5:    bilo    Lactate  5.858259e+02
```

We can see, that *B. longum* (bilo) secretes acetate as main end product.  Approximately 60 % of the acetate is consumed by *E. rectale* (eure). *B. longum* produces in addition lactate and *E. rectale* secretes butyrate. The predicted consumption of butyrate by *B. longum* with a rate of `-8.231749e-13` should be considered as zero flux, as this number might be due to numeric instabilities in R or the LP-solver (here CPLEX).
