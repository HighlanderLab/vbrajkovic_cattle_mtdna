# Title: "Quantifying the effects of the mitochondrial genome on milk production traits in dairy cows: empirical results and modelling challenges"
# Authors: Vladimir Brajkovic, Ivan Pocrnic, Miroslav Kaps, Marija Špehar, Vlatka Cubric-Curik, Strahil Ristov, Dinko Novosel, Gregor Gorjanc, Ino Curik
# DOI: 10.3168/jds.2024-25203
# Please cite: Brajkovic et al., 2024 

# Last Update: May 2024

# Tested on Intel macOS Sonoma 14.2.1 and Windows 10 64-bit

# R version 4.4.0 and RStudio 2024.04.0 Build 735

# Due to the new version of R update for some packages like INLA, it is necessary to install a test version

# ---- Setup ----

# Cleanup
rm(list = ls())

# Set working directory

setwd("/path-to-your-working-directory/")
getwd()

# The order of the packages is as they are required in the code

install.packages("tidyverse")
install.packages("pedigreemm")
install.packages("nadiv")
install.packages("terra")
install.packages("raster")
install.packages("INLA")

library(package = "tidyverse")
library(package = "pedigreemm")
library(package = "nadiv") # installed from Package Archive file (.tgz; .tar.gz) https://cran.r-project.org/src/contrib/Archive/nadiv/nadiv_2.17.2.tar.gz
library(package = "terra")    
library(package = "raster")   
library(package = "INLA")

## Set and get global options for INLA

## check it
#inla.getOption("num.threads")

## set number of threads
#inla.setOption("num.threads", 3)

## check it after setting it up
#inla.getOption("num.threads")


# Set source for the functions_v5.R

source("/path-to-your-working-directory/functions_v5.R")



###################################
# ---- Data: pedigree import ---- #
###################################

Peds = read_csv(file = "Pedigree.csv", col_names = TRUE)
head(Peds)
nrow(Peds) # 6336

Ped = arrange(Peds, animalID)

# Check if animalID, sire, dam  is sorted - animals after the parents

Tmp = Ped %>%
  filter(!is.na(sire))
Test = Tmp$sire > Tmp$animalID
if (any(Test)) {
  stop("Wrong sire-animal coding")
}
Tmp = Ped %>%
  filter(!is.na(dam))
Test = Tmp$dam > Tmp$animalID
if (any(Test)) {
  stop("Wrong dam-animal coding")
}

sum(is.na(Ped$animalID))
sum(is.na(Ped$sire))
sum(is.na(Ped$dam))

# ---- Build pedigree precision matrix ---- #

pedX = pedigree(label = Ped$animalID,
                sire  = Ped$sire,
                dam   = Ped$dam)
#AAA = getA(ped = pedX)
TInv = as(pedX, "sparseMatrix")
DInv = Diagonal(x = 1 / Dmat(pedX))
QPedigree = crossprod(sqrt(DInv) %*% TInv)

# ---- Build Sex precision matrix ---- #

# Order of ped for this function must be: animal, dam , sire, sex !!!

pedisex = dplyr::select(Ped, animalID, dam, sire, sex)
table(pedisex$sex)
#   1    2
# 993 5343

# Function needs different format for the table, e.g. matrix
pedisex2 = as.matrix(pedisex)

QSex = makeS(pedisex2, heterogametic = "1", returnS = TRUE)

QSexsp = QSex$Sinv



#########################################
# ---- Data: phenotype data import ---- #
#########################################

Phe = read_csv(file = "Phenotype_data.csv", col_names = TRUE,
               col_types = cols(.default = col_character(),
                                lab_id = col_character(),
                                animalID = col_double(),
                                sire = col_double(),
                                dam = col_double(),
                                date_of_birth = col_date(format = "%Y-%m-%d"), 
                                number_of_calving = col_double(),
                                date_of_calving = col_date(format = "%Y-%m-%d"), 
                                age_at_calving = col_double(),
                                county = col_double(),
                                milk_yield = col_double(),
                                fat_yield = col_double(),
                                protein_yield = col_double(),
                                herd = col_double(),
                                herd_c = col_double(),
                                hys = col_double(), #herd-year-season_of_calving
                                cyto = col_double(), # founders
                                haplo = col_double(), # haplotpyes
                                amino = col_double(), # aminohaplotypes
                                evol = col_double())) # evolutionary haplotypes

head(Phe)
nrow(Phe) # 7576

Test = !Phe$animalID %in% Ped$animalID
if (any(Test)) {
  stop("Animals in phenotype data not in pedigree")
}

#######################################

Phe$number_of_calvingC = as.factor(Phe$number_of_calving)
Phe$countyC = as.factor(Phe$county)

# Permanent effect
Phe$animal_perm = Phe$animalID

# Animal for Sex matrix
Phe$animal_sex = Phe$animalID




######################################
# ---- Data: mtDNA polymorphisms ----#
######################################
# MEGA7 was used to converted CHF_set_109_al_2.fas to CCHF_set_109_al_2.meg. 
# DnaSP was used to open CCHF_set_109_al_2.meg and create Variable_sites_new.hap
# DnaSp procedure: Generate Haplotype Data File; Parameters: Site with gaps/missing: not considered, 
# Invariable sites: removed, Generate: Arlequin Haplotype list › save.as Variable_sites_new.hap

Snp = read_table(file = "Variable_sites_new.hap", col_names = FALSE) 
# head(Snp)
nrow(Snp) # 96 haplotypes

# Get sample haplotype codes
Snp$X1 = Snp$X1 %>%
  sub(pattern = "Hap_", replacement = "") %>%
  as.numeric()
head(Snp)

# Rename haplotypes from 1-96 to 2-97 so it is the same name of haplotypes as in MJ network and suplemental table
Snp$X1 = Snp$X1 + 1

# Convert SNP nucleotides to 0/1 allele dosages
Snp2 = Snp$X2 %>%
  strsplit(split = "")
names(Snp2) = 1:length(Snp2)
Snp2 = as_tibble(Snp2) %>%
  t() %>% as_tibble()
dim(Snp2) # 96 358
Snp2
Snp2 = cbind(2:(nrow(Snp2) + 1), Snp2) %>% as_tibble() 
dim(Snp2) # 96 359
nSnp = ncol(Snp2) - 1
colnames(Snp2) = c("Hap", paste("Snp_", 1:nSnp, sep = ""))
Snp2
dim(Snp2) # 96 359

# Check for multi-allelism
Tmp = vector(mode = "list", length = nSnp)
for (Snp in 1:nSnp) {
  # Snp = 1
  Tmp[[Snp]] = table(Snp2[[paste("Snp_", Snp, sep = "")]])
  TestMulti = length(Tmp[[Snp]]) > 2
  if (TestMulti) {
    cat("Snp:", Snp)
    print(Tmp[[Snp]])
    print("Multi-allelism")
  }
}

## Snp 316 is tri-allelic!
# Recode and ensure proper coding of multi-allelism
tmp3 = data.frame(matrix(NA, nrow = 96, ncol = nSnp+1))
dim(tmp3) # 96 359
Snp3 = cbind(Snp2$Hap, tmp3) %>% as_tibble()
colnames(Snp3) = c("Hap",
                   paste("Snp_", 1:315, sep = ""),
                   "Snp_316_GA",
                   "Snp_316_GC",
                   paste("Snp_", 317:358, sep = ""))
Tmp = vector(mode = "list", length = nSnp)
Cov = 0
for (Snp in 1:nSnp) {
  # Snp = 1
  Cov = Cov + 1
  Nucleotides = Snp2[[paste("Snp_", Snp, sep = "")]]
  Tmp[[Snp]] = table(Nucleotides)
  RefAllele = names(which.max(Tmp[[Snp]]))
  AltAllele = names(Tmp[[Snp]])[!names(Tmp[[Snp]]) %in% RefAllele]
  TestMulti = length(Tmp[[Snp]]) > 2
  if (TestMulti) {
    TestRefAllele = Nucleotides == RefAllele
    for (Alt in AltAllele) {
      # Alt = "A"
      # Alt = "C"
      TestAltAllele = Nucleotides == Alt
      Snp3[             , Cov + 1] = 0
      Snp3[TestAltAllele, Cov + 1] = 1
      Cov = Cov + 1
    }
    Cov = Cov - 1
  } else {
    TestRefAllele = Nucleotides == RefAllele
    Snp3[              ,  Cov + 1] = 0
    Snp3[!TestRefAllele, Cov + 1] = 1
  }
}

dim(Snp3) # 96 360

# Merge with phenotype data

Test = Snp3$Hap %in% Phe$haplo
sum(Test) # 96
if (sum(Test) != 96) {
  Snp3$Hap[!Test]
  stop("Something went wrong with haplotype coding")
}
Phe = merge(x = Phe, y = Snp3, by.x = "haplo", by.y = "Hap")
nrow(Phe) # 7576

# ---- Data: recalculate age at first calving ---- #

# hist(Phe$age_at_calving)
Tmp = as.numeric((Phe$date_of_calving - Phe$date_of_birth) / 30)
# plot(Tmp ~ Phe$age_at_calving)
sd(Tmp)                  # 22.86
sd(Phe$age_at_calving) # 22.49
Phe$age_at_calving = Tmp

# ---- Data: scale phenotypes and covariates ---- #

(milk_yieldMean = mean(Phe$milk_yield, na.rm = TRUE)) # 7158.711
(milk_yieldSDev = sd(Phe$milk_yield, na.rm = TRUE)) # 1843.216
Phe$milk_yieldS = (Phe$milk_yield - milk_yieldMean) / milk_yieldSDev

(fat_yieldMean = mean(Phe$fat_yield, na.rm = TRUE)) # 278.1104
(fat_yieldSDev = sd(Phe$fat_yield, na.rm = TRUE)) # 79.11727
Phe$fat_yieldS = (Phe$fat_yield - fat_yieldMean) / fat_yieldSDev

(protein_yieldMean = mean(Phe$protein_yield, na.rm = TRUE)) # 235.1427
(protein_yieldSDev = sd(Phe$protein_yield, na.rm = TRUE)) # 60.49171
Phe$protein_yieldS = (Phe$protein_yield - protein_yieldMean) / protein_yieldSDev

(age_at_calvingMean = Phe %>%
    group_by(number_of_calvingC) %>%
    summarise(age_at_calvingMean = mean(age_at_calving, na.rm = TRUE)))
(age_at_calvingSDev = Phe %>%
    group_by(number_of_calvingC) %>%
    summarise(age_at_calvingSDev = sd(age_at_calving, na.rm = TRUE)))
age_at_calvingStats = merge(x = age_at_calvingMean, y = age_at_calvingSDev)
Phe = merge(x = Phe, y = age_at_calvingStats)
Phe$age_at_calvingS = (Phe$age_at_calving - Phe$age_at_calvingMean) / Phe$age_at_calvingSDev




##################################################
# ---- Repeatability model with Barrier SPDE ----# 
##################################################

length(unique(Phe$cyto))
# 109
length(unique(Phe$haplo))
# 96
length(unique(Phe$amino))
# 48
length(unique(Phe$evol))
# 10

PheLact15 = Phe %>%
  filter(number_of_calving < 6) %>%
  droplevels()
nrow(PheLact15) # 7115

PheLact15$herd15_c = PheLact15$herd %>% factor() %>% as.numeric()

###########################################
# ---- Data: coordinates of each herd ----#
###########################################

herd_coordinates = read.csv("Herd_coordinates.csv", header = TRUE)

# Add coordinates to each animal based on the herd:
PheLact15loc = merge(x = PheLact15, y = herd_coordinates, by.x = "herd", by.y = "herd")

# Add map of Croatia:
mapaHR = getData('GADM', country="HRV", level= 1) # gadm36_HRV_1_sp.rds file is made available in the working directory, as the web is sometimes not available

library(sf) #activate the package, it is pre-installed in the INLA package

# Extract the border from the map
mapaHR <- st_as_sf(mapaHR)
hr.border <- mapaHR %>% st_set_precision(1e8) %>% summarize()

# Formatting boundary for INLA
hr.bdry = inla.sp2segment(hr.border) 

#1. Mesh:
locations = cbind(herd_coordinates$long, herd_coordinates$lat)

mesh3 = inla.mesh.2d(loc = locations, boundary = hr.bdry, max.edge = c(0.4, 0.8), cutoff = 0.08, offset = c(0.7, 0.7))

plot(mesh3, main="3rd mesh"); points(herd_coordinates$long, herd_coordinates$lat, col = "red") 

# Defining the spatial barriers:
# Pick out which triangles in the mesh belong to the barrier area:
tl = length(mesh3$graph$tv[,1])
# - the number of triangles in the mesh
posTri = matrix(0, tl, 2)
for (t in 1:tl){
  temp = mesh3$loc[mesh3$graph$tv[t, ], ]
  posTri[t,] = colMeans(temp)[c(1,2)] 
}

posTriCent = SpatialPoints(posTri)
# The positions of the triangle centres

#Add projection attribute:
raster::crs(posTriCent) = "+proj=longlat +datum=WGS84 +no_defs"

# barrier = over(hr.border, posTriCent, returnList = T) 

posTriCent = st_as_sf(posTriCent)
barrier = st_intersects(posTriCent, hr.border)

#pointsSP <- SpatialPoints(posTri, proj4string = CRS(proj4string(hr.border))) #################### https://stackoverflow.com/questions/64898152/how-to-fix-error-in-over-identicalcrsx-y-is-not-true
#barrier = over(hr.border, pointsSP , returnList = T)   #################### solution if problem occur https://stackoverflow.com/questions/64898152/how-to-fix-error-in-over-identicalcrsx-y-is-not-true  

# - checking which mesh triangles are inside the barrier area
barrier = unlist(barrier)
barrier.triangles = setdiff(1:tl, barrier)

poly.barrier = inla.barrier.polygon(mesh3, barrier.triangles)
# - the Barrier model's polygon
# - in most cases this should be the same as poly.original
# Ignore the warning messages -> same as https://haakonbakkagit.github.io/btopic107.html 

plot(mesh3) #################### 
plot(poly.barrier, add=T, col='lightblue')
plot(mesh3, add=T); points(locations, col='red', pch=20, cex=0.5)

#Define barrier model:
barrier.model = inla.barrier.pcmatern(mesh3, barrier.triangles = barrier.triangles, prior.range = c(0.3, 0.5), prior.sigma = c(1, 0.1) )


# 2.Projector matrix (A-Matrix)
# To connect the measurement locations to the mesh representation

A3 = inla.spde.make.A(mesh = mesh3, loc = cbind(PheLact15loc$long, PheLact15loc$lat))

dim(A3)
# 7115 807
# No. Obs. x No. mesh param.

# Create a Stack and model for INLA

# Model matrix, and have 1 included to indicate that we have an intercept in the model
MM = model.matrix(~1 + number_of_calvingC*age_at_calvingS, PheLact15loc)

# Intercept is together with the SPDE field in the stack, so it is removed from the factors 
FactorNames =  tail(colnames(MM), ncol(MM)-1)

# Prepare for adding the factors to the stack
FactorStack = paste0( "list(", FactorNames, " = MM[,", 2:(ncol(MM)), "]),", collapse = "" )
FactorStack = ( gsub('.{1}$', '', FactorStack ) )
FactorStack = ( gsub(':', 'xx', FactorStack ) )
FactorStack = paste0("c(", FactorStack, ")")

stackbarrier = inla.stack(data = list(y = PheLact15loc$milk_yieldS),
                          A = list(A3, 1, 1), tag = "Estimate",
                          effects=list(c(list(Intercept = 1),
                                         list(i = 1:mesh3$n)),
                                       c(list(hys = PheLact15loc$hys),
                                         list(animal = PheLact15loc$animalID),
                                         list(animal_sex = PheLact15loc$animal_sex),
                                         list(animal_perm = PheLact15loc$animal_perm),
                                         list(haplo = PheLact15loc$haplo)),
                                       c(eval(parse(text = FactorStack)))
                          ))

FactorNames2 = ( gsub(':', 'xx', FactorNames ) )

# Basic haplotype-based model (Mito haplotype effect is IID) 
ModelHapBarrier = c("0", "Intercept", 
                        "f(hys, model = \"iid\")", 
                        "f(animal, model = \"generic0\", Cmatrix = QPedigree)",
                        "f(haplo, model = \"iid\")", 
                        "f(animal_sex, model = \"generic0\", Cmatrix = QSexsp)", 
                        "f(animal_perm, model = \"iid\")", 
                        "f(i, model = barrier.model)")

FormulaFactor = c(ModelHapBarrier, FactorNames2)  

y = "y"

ModelMilkHapBarrier = as.formula(paste(y, paste(FormulaFactor, collapse=" + "), sep=" ~ "))

# Put all together and run INLA:
FitMilkLact15HapBarrier = inla(formula = ModelMilkHapBarrier, data = inla.stack.data(stackbarrier),
                            control.predictor = list(A = inla.stack.A(stackbarrier), compute = T),
                            control.compute = list(mlik = T, dic = T, waic = T))


# 4. Summarise the results
summary(FitMilkLact15HapBarrier)
SummarizeInlaVars(FitMilkLact15HapBarrier)

capture.output(summary(FitMilkLact15HapBarrier), file = "Summary_FitMilkLact15HapBarrier_SF.txt")
write.csv(SummarizeInlaVars(FitMilkLact15HapBarrier), file = "Inla_var_FitMilkLact15HapBarrier_SF.csv")

# Variance for spatial component:
teta1 = FitMilkLact15HapBarrier$internal.marginals.hyperpar[[7]]
# Theta 1 = log(SD)
# To get VAR:
te1 = inla.tmarginal(function(x) (exp(x))^2, teta1)
rm(Samples)
Samples = matrix(data = numeric(), nrow = 10000, ncol = 1)

Samples[, 1] = inla.rmarginal(n = 10000, te1)
mean(Samples) #################### 1.764439

# Range for spatial component:
teta2 = FitMilkLact15HapBarrier$internal.marginals.hyperpar[[8]]
# Theta 2 = log(range)
# To get range:
te2 = inla.tmarginal(function(x) (exp(x)), teta2)
rm(Samples)
Samples = matrix(data = numeric(), nrow = 10000, ncol = 1)
Samples[, 1] = inla.rmarginal(n = 10000, te2)
mean(Samples) #################### 0.8078734




###########################################################
# ---- Repeatability model with SNP Markers and SPDE ---- #  
###########################################################


# Uses barrier SPDE model defined previously

# Add marker necessities: 
# Create Z matrix with unique haplotypes only:

PheLact15un = PheLact15 %>%
  distinct(haplo, .keep_all = TRUE)
nrow(PheLact15un) #96
# hcode = PheLact15un$haplo
# sort(hcode)
# hcode 2-97 (1 was reference and was removed)
SnpTmpun = PheLact15un %>%
  dplyr::select(starts_with("SNP_"), "haplo")
SnpZun = as.matrix(SnpTmpun)
rownames(SnpZun) = SnpTmpun$haplo
SnpZun = SnpZun[order(as.numeric(rownames(SnpZun))), -ncol(SnpZun)]
dim(SnpZun) #96 359

###### SNP/Z-model by using stack

# Design matrix to connect those unique haplotypes with (repeated) measurements
# Done by haplo
# So haplo #2 is 1. column in M, #3 is 2., etc...

M = sparse.model.matrix( ~ 0 + (as.factor(PheLact15$haplo)))
dim(M) # 7115 x 96

stacksnpspde = inla.stack(data = list(y = PheLact15$milk_yieldS),
                          A = list(A3, M, 1, 1), tag = "Estimate",
                          effects=list(c(list(Intercept = 1),
                                         list(i = 1:mesh3$n)), 
                                       c(list(snpovi = 1:ncol(M))), 
                                       c(list(hys = PheLact15$hys),
                                         list(animal = PheLact15$animalID),
                                         list(animal_sex = PheLact15$animal_sex),
                                         list(animal_perm = PheLact15$animal_perm)),
                                       c(eval(parse(text = FactorStack)))
                          ))


# Basic SNP (Z) model:
ModelSNP = c("0", "Intercept", 
             "f(hys, model = \"iid\")", 
             "f(animal, model = \"generic0\", Cmatrix = QPedigree)",
             "f(animal_sex, model = \"generic0\", Cmatrix = QSexsp)", 
             "f(animal_perm, model = \"iid\")", 
             "f(snpovi, model =\"z\", Z = SnpZun)", 
             "f(i, model = barrier.model)")

FormulaFactor = c(ModelSNP, FactorNames2)  

y = "y"

ModelMilkZ = as.formula(paste(y, paste( FormulaFactor, collapse=" + "), sep=" ~ "))

# Put all together and run INLA:
FitMilkLact15ZSstack = inla(formula = ModelMilkZ, data = inla.stack.data(stacksnpspde),
                            control.predictor=list(A = inla.stack.A(stacksnpspde), compute = T),
                            control.compute = list(mlik = T, dic = T, waic = T, config = TRUE))

summary(FitMilkLact15ZSstack)
SummarizeInlaVars(FitMilkLact15ZSstack)

capture.output(summary(FitMilkLact15ZSstack), file = "Summary_FitMilkLact15ZSstack_SF.txt")
write.csv(SummarizeInlaVars(FitMilkLact15ZSstack), file = "Inla_var_FitMilkLact15ZSstack_SF.csv")

# Variance for spatial component:
teta1 = FitMilkLact15ZSstack$internal.marginals.hyperpar[[7]]
# Theta 1 = log(SD)
# To get VAR:
te1 = inla.tmarginal(function(x) (exp(x))^2, teta1)
rm(Samples)
Samples = matrix(data = numeric(), nrow = 10000, ncol = 1)
Samples[, 1] = inla.rmarginal(n = 10000, te1)
mean(Samples) #################### 1.023772

# Range for spatial component:
teta2 = FitMilkLact15ZSstack$internal.marginals.hyperpar[[8]]
# Theta 2 = log(range)
# To get range:
te2 = inla.tmarginal(function(x) (exp(x)), teta2)
rm(Samples)
Samples = matrix(data = numeric(), nrow = 10000, ncol = 1)
Samples[, 1] = inla.rmarginal(n = 10000, te2)
mean(Samples) #################### 0.3972663




############################ 
# ---- Sample Effects ---- # 
############################ 


# Sample
sample1 <- inla.posterior.sample(1000, FitMilkLact15ZSstack)

# save.image("aftermarkersample_new.RData")
# load("aftermarkersample_new.RData")


# Pull out correct part from samples: 
# https://haakonbakka.bitbucket.io/btopic112.html 

contents = FitMilkLact15ZSstack$misc$configs$contents
effect = "snpovi"
id.effect = which(contents$tag==effect)
ind.effect = contents$start[id.effect]-1 + (1:contents$length[id.effect])

samples.effect = lapply(sample1, function(x) x$latent[ind.effect])
s.eff = matrix(unlist(samples.effect), byrow = T, nrow = length(samples.effect))
colnames(s.eff) = rownames(sample1[[1]]$latent)[ind.effect]

# Compare the variance from sampling to the variances from the INLA model summary:
# For the SNP solutions:
sol_m = s.eff[,97:455]
# Calculate variance:
mean(matrixStats::rowVars(sol_m)) #################### 0.02122246
# Equal to var from hyperparam: 
mean(1 / inla.rmarginal(n = 1000, marginal = FitMilkLact15ZSstack$marginals.hyperpar[["Precision for snpovi"]])) #################### 0.02069168

# For the haplotype solutions:
sol_haps = s.eff[,1:96]
# Calculate Var(hap)
mean(matrixStats::rowVars(sol_haps)) #################### 0.08430112
# The haplotype solutions have to be from sampling

# Residual has to go from the hyperparameter:
mean(1 / inla.rmarginal(n = 1000, marginal = FitMilkLact15ZSstack$marginals.hyperpar[["Precision for the Gaussian observations"]])) #################### 0.3508088
mean(unlist(lapply(sample1, function(x) 1/(x$hyperpar["Precision for the Gaussian observations"])))) #################### 0.3505232


#########################
# Samples by region  ####
#########################


#####################################################
# ---- Data: SNP position in mitogenome regions ----#
#####################################################

snploc = read.csv("SNP_position_region_new.csv", header = TRUE)

# 358 not 359 --> because we had one SNP triallelic
# table(snploc$Region)
# length(unique(snploc$Region))
# 29 regions

# Construct 1000 covariance matrices, based on the 1000 SNP-effect samples 
# Hmat = Var(Z*m_i)
rm(Hmat)
Hmat = list()
for(r in 1:1000){
  sol_mm = sol_m[r,]
  SNPmatrix = SnpZun
  for(i in 1:359){
    SNPmatrix[,i] = SNPmatrix[,i] * sol_mm[i]
  }
  colnames(SNPmatrix) = paste0("Snp_", 1:ncol(SNPmatrix))
  # Create covariance matrix:
  SNPmatrixVar = var(SNPmatrix)
  # Calculate statistics:
  totalGeneticVar = sum(SNPmatrixVar)
  totalGenicVar = sum(diag(SNPmatrixVar))
  totalLD = sum(SNPmatrixVar[upper.tri(SNPmatrixVar, diag = FALSE)])
  # Extract Genic variance of each region:  
  Rmat = matrix(, ncol = 1, nrow = 29) 
  Rmat[1,] = sum(diag(SNPmatrixVar[1:12,1:12])) 
  Rmat[2,] = sum(diag(SNPmatrixVar[13:25,13:25])) 
  Rmat[3,] = sum(SNPmatrixVar[26:26,26:26])
  Rmat[4,] = sum(diag(SNPmatrixVar[27:44,27:44])) 
  Rmat[5,] = sum(SNPmatrixVar[45:45,45:45])
  Rmat[6,] = sum(diag(SNPmatrixVar[46:66,46:66])) 
  Rmat[7,] = sum(SNPmatrixVar[67:67,67:67])
  Rmat[8,] = sum(SNPmatrixVar[68:68,68:68])
  Rmat[9,] = sum(diag(SNPmatrixVar[69:90,69:90])) 
  Rmat[10,] = sum(SNPmatrixVar[91:91,91:91]) 
  Rmat[11,] = sum(SNPmatrixVar[92:92,92:92]) 
  Rmat[12,] = sum(diag(SNPmatrixVar[93:117,93:117])) 
  Rmat[13,] = sum(SNPmatrixVar[118:118,118:118]) 
  Rmat[14,] = sum(diag(SNPmatrixVar[119:128,119:128])) 
  Rmat[15,] = sum(diag(SNPmatrixVar[129:134,129:134])) 
  Rmat[16,] = sum(diag(SNPmatrixVar[135:146,135:146])) 
  Rmat[17,] = sum(diag(SNPmatrixVar[147:162,147:162])) 
  Rmat[18,] = sum(diag(SNPmatrixVar[163:169,163:169])) 
  Rmat[19,] = sum(SNPmatrixVar[170:170,170:170]) 
  Rmat[20,] = sum(diag(SNPmatrixVar[171:174,171:174])) 
  Rmat[21,] = sum(diag(SNPmatrixVar[175:209,175:209])) 
  Rmat[22,] = sum(diag(SNPmatrixVar[210:211,210:211])) 
  Rmat[23,] = sum(diag(SNPmatrixVar[212:254,212:254])) 
  Rmat[24,] = sum(diag(SNPmatrixVar[255:270,255:270])) 
  Rmat[25,] = sum(SNPmatrixVar[271:271,271:271])
  Rmat[26,] = sum(diag(SNPmatrixVar[272:293,272:293])) 
  Rmat[27,] = sum(SNPmatrixVar[294:294,294:294])
  Rmat[28,] = sum(diag(SNPmatrixVar[295:296,295:296])) 
  Rmat[29,] = sum(diag(SNPmatrixVar[297:358,297:358])) 
  #
  Hmat[[r]] = list(SNPmatrixVar, totalGeneticVar, totalGenicVar, totalLD, Rmat)
  names(Hmat[[r]]) = c("CoV", "totalGeneticVar", "totalGenicVar", "TotalLD", "RegionGenic")
}

# Check dimensions:
dim(Hmat[[1]]$CoV)
# 359x359
dim(Hmat[[1]]$RegionGenic)
# 29x1

# Take the mean over 1000 samples:
Genetic = lapply(Hmat, function(x) x$totalGeneticVar)
mean(unlist(Genetic)) #################### 0.08430144
# This corresponds to Var(hap) 

Genic = lapply(Hmat, function(x) x$totalGenicVar)
mean(unlist(Genic)) #################### 0.1096756

TotLD = lapply(Hmat, function(x) x$TotalLD) 
mean(unlist(TotLD)) #################### -0.01268706

GenicRegion = lapply(Hmat, function(x) x$RegionGenic)
GenicRegion_mean = apply(simplify2array(GenicRegion), c(1,2), mean)

###### Construct covariance matrix per region (Compress presented above by region) ####

rm(RegMat)
RegMat = list()
for(r in 1:1000){
  sol_mm = sol_m[r,]
  SNPmatrix = SnpZun
  for(i in 1:359){
    SNPmatrix[,i] = SNPmatrix[,i] * sol_mm[i]
  }
  Rmat = matrix(, ncol = length(unique(snploc$Region)), nrow = nrow(SNPmatrix)) 
  Rmat[,1] = as.matrix(rowSums(SNPmatrix[,1:12])) 
  Rmat[,2] = as.matrix(rowSums(SNPmatrix[,13:25])) 
  Rmat[,3] = as.matrix(SNPmatrix[,26:26])
  Rmat[,4] = as.matrix(rowSums(SNPmatrix[,27:44])) 
  Rmat[,5] = as.matrix(SNPmatrix[,45:45])
  Rmat[,6] = as.matrix(rowSums(SNPmatrix[,46:66])) 
  Rmat[,7] = as.matrix(SNPmatrix[,67:67])
  Rmat[,8] = as.matrix(SNPmatrix[,68:68])
  Rmat[,9] = as.matrix(rowSums(SNPmatrix[,69:90])) 
  Rmat[,10] = as.matrix(SNPmatrix[,91:91]) 
  Rmat[,11] = as.matrix(SNPmatrix[,92:92]) 
  Rmat[,12] = as.matrix(rowSums(SNPmatrix[,93:117])) 
  Rmat[,13] = as.matrix(SNPmatrix[,118:118]) 
  Rmat[,14] = as.matrix(rowSums(SNPmatrix[,119:128])) 
  Rmat[,15] = as.matrix(rowSums(SNPmatrix[,129:134])) 
  Rmat[,16] = as.matrix(rowSums(SNPmatrix[,135:146])) 
  Rmat[,17] = as.matrix(rowSums(SNPmatrix[,147:162])) 
  Rmat[,18] = as.matrix(rowSums(SNPmatrix[,163:169])) 
  Rmat[,19] = as.matrix(SNPmatrix[,170:170]) 
  Rmat[,20] = as.matrix(rowSums(SNPmatrix[,171:174])) 
  Rmat[,21] = as.matrix(rowSums(SNPmatrix[,175:209])) 
  Rmat[,22] = as.matrix(rowSums(SNPmatrix[,210:211])) 
  Rmat[,23] = as.matrix(rowSums(SNPmatrix[,212:254])) 
  Rmat[,24] = as.matrix(rowSums(SNPmatrix[,255:270])) 
  Rmat[,25] = as.matrix(SNPmatrix[,271:271])
  Rmat[,26] = as.matrix(rowSums(SNPmatrix[,272:293])) 
  Rmat[,27] = as.matrix(SNPmatrix[,294:294])
  Rmat[,28] = as.matrix(rowSums(SNPmatrix[,295:296])) 
  Rmat[,29] = as.matrix(rowSums(SNPmatrix[,297:358])) 
  # dim(Rmat)
  # 96x29
  colnames(Rmat) = paste0("Region_", 1:ncol(Rmat))
  # Create covariance matrix:
  RegionMatrixVar = var(Rmat)
  # Calculate statistics:
  totalGeneticVar = sum(RegionMatrixVar)
  # totalGenicVar = sum(diag(RegionMatrixVar))
  totalLD = sum(RegionMatrixVar[upper.tri(RegionMatrixVar, diag = FALSE)])
  RegMat[[r]] = list(RegionMatrixVar, totalGeneticVar, totalLD)
  names(RegMat[[r]]) = c("CoV", "totalGeneticVar", "TotalLD")
  # GeneticVar = GenicVar + 2*totalLD
}

# Check dimensions:
dim(RegMat[[1]]$CoV)
# 29 29

# Take the means over 1000 samples:
TotGenetic = lapply(RegMat, function(x) x$totalGeneticVar)
mean(unlist(TotGenetic)) #################### 0.08323125
# Results comparable to Hmat and Var(hap)

TotLD2 = lapply(RegMat, function(x) x$TotalLD) 
mean(unlist(TotLD)) #################### -0.01268706
# Results comparable to Hmat

# Genetic variance for each Region:
GeneticVarRegions = NULL
for(i in 1:29){
  ri = lapply(RegMat, function(x) x$CoV[[i,i]])
  rim = mean(unlist(ri))
  GeneticVarRegions = rbind(GeneticVarRegions, rim)
}

# LD between the regions:
rm(LDBetweenRegions)
LDBetweenRegions = list()
for(r in 1:1000){
  ri3 = NULL
  for(i in 1:29){ 
    ri = as.matrix(RegMat[[r]]$CoV[,i])
    ri2 = sum(ri[-i,])
    ri3 = rbind(ri3, ri2)
  }
  LDBetweenRegions[[r]] = ri3
  # This number is already total LD between the regions - no need to 2xLD later
}

LDb_mean = apply(simplify2array(LDBetweenRegions), 1, mean) 

# Prepare table for regions:
# Genetic Var, Genic Var, LD between, LD within:
# Note: LDWithinRegion was calculated from the means not actually sampled!
rm(table1)
table1 = round(cbind(GeneticVarRegions, GenicRegion_mean, GeneticVarRegions-GenicRegion_mean, LDb_mean), 5)
colnames(table1) = c("GeneticVarRegion", "GenicVarRegion", "LDWithinRegion", "LDBetweenRegions") 
row.names(table1) = unique(snploc$Region)

table2 = rbind(table1, colSums(table1))
row.names(table2)[30] = c("SUM")

write.csv(table2, "~/Desktop/table2_milk.csv") 

# Take samples from posterior for all the effects
# Calculate covariances and correlations (Effects and EffectxRegions)

# Calculate proper variances and ratios for SNP model

# Based on INLA posterior sample (as all calculations before)
# sample1 <- inla.posterior.sample(1000, FitMilkLact15ZSstack)


SampleEfektPost = function(modelfit, efekt){
  contents = modelfit$misc$configs$contents
  effect = efekt
  id.effect = which(contents$tag==effect)
  ind.effect = contents$start[id.effect]-1 + (1:contents$length[id.effect])
  samples.effect = lapply(sample1, function(x) x$latent[ind.effect])
  s.eff = matrix(unlist(samples.effect), byrow = T, nrow = length(samples.effect))
  colnames(s.eff) = rownames(sample1[[1]]$latent)[ind.effect]
  return(s.eff)
}

# Take posterior samples of random effects
# Remove column name "effect_name:effect_ID" and take name from model summary
# Z-model: Has to be separated to hap and marker part
# SPDE: Has to be multiplied by "A3" matrix to get solution for each animal

sol_a = SampleEfektPost(FitMilkLact15ZSstack, "animal")
colnames(sol_a) = colnames(sol_a) %>% sub(pattern = "animal:", replacement = "", x = ., fixed = TRUE)

sol_s=SampleEfektPost(FitMilkLact15ZSstack, "animal_sex")
colnames(sol_s) = colnames(sol_s) %>% sub(pattern = "animal_sex:", replacement = "", x = ., fixed = TRUE)

sol_p=SampleEfektPost(FitMilkLact15ZSstack, "animal_perm")
# colnames(sol_p) = colnames(sol_p) %>% sub(pattern = "animal_perm:", replacement = "", x = ., fixed = TRUE)
colnames(sol_p) = FitMilkLact15ZSstack$summary.random[4]$animal_perm$ID

sol_hys=SampleEfektPost(FitMilkLact15ZSstack, "hys")
# colnames(sol_hys) = colnames(sol_hys) %>% sub(pattern = "hys:", replacement = "", x = ., fixed = TRUE)
colnames(sol_hys) = FitMilkLact15ZSstack$summary.random[1]$hys$ID

sol_Z = SampleEfektPost(FitMilkLact15ZSstack, "snpovi")
# colnames(sol_Z) = colnames(sol_Z) %>% sub(pattern = "snpovi:", replacement = "", x = ., fixed = TRUE)

sol_haps =  sol_Z[,1:96]
colnames(sol_haps) = seq(2,97)

sol_m = sol_Z[,97:455]  
colnames(sol_m) = seq(1, 359)

sol_l = SampleEfektPost(FitMilkLact15ZSstack, "i")
# colnames(sol_l) = colnames(sol_l) %>% sub(pattern = "i:", replacement = "", x = ., fixed = TRUE)
sol_ll=NULL
for(r in 1:1000){
  AA = A3
  spatialAllAnim = as.vector(AA %*% sol_l[r,])
  sol_ll = rbind(sol_ll, spatialAllAnim)
}
colnames(sol_ll) = PheLact15loc$animal
# This effects is following order of dataset used to create A and mesh for SPDE models
# That dataset has different order than the other one used in analysis  

# They all have different dimensions
# Create design matrices for each of them (for a,s, and p are the same - I created anyways to be clear)
# Those design matrices will multiply solutions in loop for each replicate

Za = sparse.model.matrix( ~ 0 + (as.factor(PheLact15$animalID)))
colnames(Za) = colnames(Za) %>% sub(pattern = "as.factor(PheLact15$animalID)", replacement = "", x = ., fixed = TRUE)

Zs = sparse.model.matrix( ~ 0 + (as.factor(PheLact15$animal_sex)))
colnames(Zs) = colnames(Zs) %>% sub(pattern = "as.factor(PheLact15$animal_sex)", replacement = "", x = ., fixed = TRUE)

Zp = sparse.model.matrix( ~ 0 + (as.factor(PheLact15$animal_perm)))
colnames(Zp) = colnames(Zp) %>% sub(pattern = "as.factor(PheLact15$animal_perm)", replacement = "", x = ., fixed = TRUE)

Zhys = sparse.model.matrix( ~ 0 + (as.factor(PheLact15$hys)))
colnames(Zhys) = colnames(Zhys) %>% sub(pattern = "as.factor(PheLact15$hys)", replacement = "", x = ., fixed = TRUE)

Zhap = sparse.model.matrix( ~ 0 + (as.factor(PheLact15$haplo)))
colnames(Zhap) = colnames(Zhap) %>% sub(pattern = "as.factor(PheLact15$haplo)", replacement = "", x = ., fixed = TRUE)

Zloc = Za

# Sample variances and ratios:

# Residual has to go from the hyperparameter:
resid = unlist(lapply(sample1, function(x) 1/(x$hyperpar["Precision for the Gaussian observations"])))

# Main Loop: 
rm(OverSamples)
OverSamples = list()
for(r in 1:1000){
  r2 = resid[r]
  a2 = var(sol_a[r,])
  s2 = var(sol_s[r,])
  p2 = var(sol_p[r,])
  hys2 = var(sol_hys[r,])
  hap2 = var(sol_haps[r,])
  loc2 = var(unique(sol_ll[r,]))
  
  sum_all = sum(r2, a2, s2, p2, hys2, hap2, loc2)
  
  r_2 = r2/sum_all
  h_2 = a2/sum_all
  s_2 = s2/sum_all
  p_2 = p2/sum_all
  hys_2 = hys2/sum_all
  m_2 = hap2/sum_all
  loc_2 = loc2/sum_all
  
  OverSamples[[r]] = list(r2, a2, s2, p2, hys2, hap2, loc2, r_2, h_2, s_2, p_2, hys_2, m_2, loc_2)
  names(OverSamples[[r]]) = c("Ve", "Va", "Vs", "Vp", "Vhys", "Vhap", "Vloc", "r2", "h2", "s2", "p2", "hys2", "m2", "loc2")
}

#### correct variances for SNP model #####

# Variances:
# Equivalent to: SummarizeInlaVars(FitFatLact15ZSstack)
round(mean(unlist(lapply(OverSamples, function(x) x$Ve))), 4) #################### 0.3505
round(mean(unlist(lapply(OverSamples, function(x) x$Va))), 4) #################### 0.3168
round(mean(unlist(lapply(OverSamples, function(x) x$Vs))), 4) #################### 3e-04
round(mean(unlist(lapply(OverSamples, function(x) x$Vp))), 4) #################### 0.0789
round(mean(unlist(lapply(OverSamples, function(x) x$Vhys))), 4) #################### 0.1002
round(mean(unlist(lapply(OverSamples, function(x) x$Vhap))), 4) #################### 0.0843
round(mean(unlist(lapply(OverSamples, function(x) x$Vloc))), 4) #################### 0.0439
 
round(sd(unlist(lapply(OverSamples, function(x) x$Ve))), 4) #################### 0.0094
round(sd(unlist(lapply(OverSamples, function(x) x$Va))), 4) #################### 0.028
round(sd(unlist(lapply(OverSamples, function(x) x$Vs))), 4) #################### 7e-04
round(sd(unlist(lapply(OverSamples, function(x) x$Vp))), 4) #################### 0.0182
round(sd(unlist(lapply(OverSamples, function(x) x$Vhys))), 4) #################### 0.0081
round(sd(unlist(lapply(OverSamples, function(x) x$Vhap))), 4) #################### 0.0142
round(sd(unlist(lapply(OverSamples, function(x) x$Vloc))), 4) #################### 0.0069

# Ratios:
round(mean(unlist(lapply(OverSamples, function(x) x$r2))), 4) #################### 0.3597
round(mean(unlist(lapply(OverSamples, function(x) x$h2))), 4) #################### 0.3248
round(mean(unlist(lapply(OverSamples, function(x) x$s2))), 4) #################### 3e-04
round(mean(unlist(lapply(OverSamples, function(x) x$p2))), 4) #################### 0.081
round(mean(unlist(lapply(OverSamples, function(x) x$hys2))), 4) #################### 0.1029
round(mean(unlist(lapply(OverSamples, function(x) x$m2))), 4) #################### 0.0864
round(mean(unlist(lapply(OverSamples, function(x) x$loc2))), 4) #################### 0.045

round(sd(unlist(lapply(OverSamples, function(x) x$r2))), 4) #################### 0.0114
round(sd(unlist(lapply(OverSamples, function(x) x$h2))), 4) #################### 0.0256
round(sd(unlist(lapply(OverSamples, function(x) x$s2))), 4) #################### 7e-04
round(sd(unlist(lapply(OverSamples, function(x) x$p2))), 4) #################### 0.019
round(sd(unlist(lapply(OverSamples, function(x) x$hys2))), 4) #################### 0.0085
round(sd(unlist(lapply(OverSamples, function(x) x$m2))), 4) #################### 0.0138
round(sd(unlist(lapply(OverSamples, function(x) x$loc2))), 4) #################### 0.0068


