# Title: "XXXX"
# Please cite: Brajkovic et al., XXXX
# Last Update: September 2020


# ---- Setup ----

# Cleanup
rm(list = ls())

# working directory
setwd(dir = "/Users/ipocrnic/Documents/Roslin Posao/Projekti/2019 cattle_mtDNA/Data_original/2019 VladoDairyCattle")
# setwd("~/Desktop/Repeatability model with Barrier SPDE SNP without region")

library(package = "tidyverse")
library(package = "pedigreemm")
library(package = "nadiv")
library(package = "INLA")
library(package = "ggplot2")
library(package = "raster")
library(package = "maptools")
library(package = "rgeos")

# there were 2 more geo packages but I thnik not needed, we nned to check so the script runs without them...


if (.Platform$OS.type == "unix") {
  inla.setOption(mkl = TRUE)
}

inla.getOption("mkl")
inla.getOption("num.threads")

source("/Users/ipocrnic/Documents/Roslin Posao/Projekti/2019 cattle_mtDNA/ipocrnic_cattle_mtdna/funkcija_nova.R")
# source("/Users/vladimir/Desktop/6.prociscenikod_Pocrnic_3.5.2020/updatekoda/funkcija_nova.R")

# ---- Data: pedigree import ----

Peds = read_csv(file = "6.Export_import_csv_txt/Ped1_sex.csv", col_names = TRUE)
# Peds = read_csv(file = "Ped1_sex.csv", col_names = TRUE)
head(Peds)
nrow(Peds) # 6336

Ped = arrange(Peds, acode)

#Check if acode, dcode, scode is sorted - animals after the parents

Tmp = Ped %>%
  filter(!is.na(scode))
Test = Tmp$scode > Tmp$acode
if (any(Test)) {
  stop("Wrong sire-animal coding")
}
Tmp = Ped %>%
  filter(!is.na(dcode))
Test = Tmp$dcode > Tmp$acode
if (any(Test)) {
  stop("Wrong dam-animal coding")
}

sum(is.na(Ped$acode))
sum(is.na(Ped$scode))
sum(is.na(Ped$dcode))

# ---- Build pedigree precision matrix ----

pedX = pedigree(label = Ped$acode,
                sire  = Ped$scode,
                dam   = Ped$dcode)
#AAA = getA(ped = pedX)
TInv = as(pedX, "sparseMatrix")
DInv = Diagonal(x = 1 / Dmat(pedX))
QPedigree = crossprod(sqrt(DInv) %*% TInv)


# ---- Build Sex precision matrix ----

#Order of ped for this function must be: A,D,S,Sex !!!

pedisex = dplyr::select(Ped, acode, dcode, scode, sex)
table(pedisex$sex)
#   1    2
# 993 5343

#Function needs diferent format for the table, e.g. matrix
pedisex2 = as.matrix(pedisex)

QSex = makeS(pedisex2, heterogametic = "1", returnS = TRUE)

QSexsp = QSex$Sinv



# ---- Data: phenotype import ----
# file = "Hap_merge_3.csv"
Phe = read_csv(file = "6.Export_import_csv_txt/koor_feno_merge.csv", col_names = TRUE,
               col_types = cols(.default = col_character(),
                                id = col_character(),
                                x = col_double(),
                                lab_id = col_character(),
                                evohap = col_double(),
                                godina_i_sezona_teljenja = col_character(),
                                godina_i_sezona_teljenja_c = col_double(),
                                stado_godina_teljenja = col_character(),
                                stado_godina_teljenja_c = col_double(),
                                godina_teljenja = col_double(),
                                godina_teljenja_c = col_double(),
                                dob_prvog_teljenja = col_double(),
                                dob_prvog_teljenja_c = col_double(),
                                founder = col_character(),
                                founder_c = col_double(),
                                stado_godina_sezona_teljenja = col_character(),
                                stado_godina_sezona_teljenja_c = col_double(),
                                stado = col_double(),
                                stado_c = col_double(),
                                amino_hap = col_character(),
                                freq = col_double(),
                                asNL = col_double(),
                                redni = col_double(),
                                broj_teljenja = col_double(),
                                datum_rodjenja = col_date(format = "%d/%m/%Y"),
                                datum_teljenja = col_date(format = "%d/%m/%Y"),
                                godina_rodjenja = col_double(),
                                dob_kod_teljenja = col_double(),
                                mjesec_teljenja = col_double(),
                                zupanija = col_double(),
                                datum_zasusenja = col_date(format = "%d/%m/%Y"),
                                servis_period = col_double(),
                                datum_izlucenja = col_date(format = "%d/%m/%Y"),
                                vrsta_izlucenja = col_double(),
                                trajanje_laktacije = col_double(),
                                kolicina_mlijeka = col_double(),
                                kolicina_m_m = col_double(),
                                sadrzaj_m_m = col_double(),
                                kolicina_proteina = col_double(),
                                sadrzaj_proteina = col_double(),
                                stalna_okolina = col_character(),
                                sezona_teljenja = col_double(),
                                r_year = col_double(),
                                id_nu = col_double(),
                                paritya = col_double(),
                                regiona = col_double(),
                                scalva = col_double(),
                                herda = col_double(),
                                hyeara = col_double(),
                                penva = col_double(),
                                linea = col_double(),
                                idanima = col_double(),
                                COX1_1 = col_character(),
                                COX1_2 = col_character(),
                                COX2_1 = col_character(),
                                COX3_1 = col_character(),
                                COX3_2 = col_character(),
                                COX3_3 = col_character(),
                                COX3_4 = col_character(),
                                COX3_5 = col_character(),
                                COX3_6 = col_character(),
                                COX3_7 = col_character(),
                                ATP6_1 = col_character(),
                                ATP6_2 = col_character(),
                                ATP6_3 = col_character(),
                                ATP6_4 = col_character(),
                                ATP8_1 = col_character(),
                                ATP8_2 = col_character(),
                                ND1_1 = col_character(),
                                ND1_2 = col_character(),
                                ND1_3 = col_character(),
                                ND1_4 = col_character(),
                                ND1_5 = col_character(),
                                ND1_6 = col_character(),
                                ND1_7 = col_character(),
                                ND2_1 = col_character(),
                                ND2_2 = col_character(),
                                ND2_3 = col_character(),
                                ND2_4 = col_character(),
                                ND2_5 = col_character(),
                                ND3_1 = col_character(),
                                ND3_2 = col_character(),
                                ND4_1 = col_character(),
                                ND4_2 = col_character(),
                                ND4_3 = col_character(),
                                ND4_4 = col_character(),
                                ND4_5 = col_character(),
                                ND4_6 = col_character(),
                                ND4L_1 = col_character(),
                                ND5_1 = col_character(),
                                ND5_2 = col_character(),
                                ND5_3 = col_character(),
                                ND5_4 = col_character(),
                                ND5_5 = col_character(),
                                ND5_6 = col_character(),
                                ND5_7 = col_character(),
                                ND5_8 = col_character(),
                                ND5_9 = col_character(),
                                ND5_10 = col_character(),
                                ND5_11 = col_character(),
                                ND6_1 = col_character(),
                                CYTB_1 = col_character(),
                                CYTB_2 = col_character(),
                                CYTB_3 = col_character(),
                                CYTB_4 = col_character(),
                                CYTB_5 = col_character(),
                                CYTB_6 = col_character(),
                                CYTB_7 = col_character(),
                                CYTB_8 = col_character(),
                                CYTB_9 = col_character(),
                                CYTB_10 = col_character(),
                                S12_1 = col_character(),
                                S12_2 = col_character(),
                                S12_3 = col_character(),
                                S12_4 = col_character(),
                                S12_5_I = col_character(),
                                S12_6 = col_character(),
                                S12_7 = col_character(),
                                S12_8 = col_character(),
                                S12_9 = col_character(),
                                S12_10 = col_character(),
                                S12_11_I = col_character(),
                                S12_12_I = col_character(),
                                S12_13 = col_character(),
                                S12_14 = col_character(),
                                S12_15 = col_character(),
                                S12_16 = col_character(),
                                S16_1 = col_character(),
                                S16_2 = col_character(),
                                S16_3_D = col_character(),
                                S16_4 = col_character(),
                                S16_5 = col_character(),
                                S16_6 = col_character(),
                                S16_7 = col_character(),
                                S16_8 = col_character(),
                                S16_9 = col_character(),
                                S16_10 = col_character(),
                                S16_11 = col_character(),
                                S16_12 = col_character(),
                                S16_13 = col_character(),
                                S16_14 = col_character(),
                                S16_15 = col_character(),
                                S16_16 = col_character(),
                                S16_17 = col_character(),
                                S16_18 = col_character(),
                                S16_19 = col_character(),
                                S16_20 = col_character(),
                                tRNA_Val = col_character(),
                                tRNA_Leu = col_character(),
                                tRNA_Gln = col_character(),
                                tRNA_Met = col_character(),
                                tRNA_Asn = col_character(),
                                tRNA_Cys = col_character(),
                                tRNA_Arg = col_character(),
                                tRNA_Ser_1 = col_character(),
                                tRNA_Ser_2 = col_character(),
                                tRNA_Glu = col_character(),
                                tRNA_Thr_1 = col_character(),
                                tRNA_Thr_2 = col_character(),
                                amino_hap_code = col_double(),
                                hap_code = col_double(),
                                srn = col_character(),
                                drn = col_character(),
                                dtbirth = col_character(),
                                a1 = col_double(),
                                s1 = col_double(),
                                d1 = col_double(),
                                acode = col_double(),
                                scode = col_double(),
                                dcode = col_double()))
head(Phe)
nrow(Phe) # 7576

Test = !Phe$acode %in% Ped$acode
if (any(Test)) {
  stop("Animals in phenotype data not in pedigree")
}

###########################################################
###########################################################

Phe$broj_teljenjaC = as.factor(Phe$broj_teljenja)
Phe$zupanijaC = as.factor(Phe$zupanija)

# Permanent effect
Phe$acode_perm = Phe$acode

# Acode for Sex matrix
Phe$acode_sex = Phe$acode

# ---- Data: mtDNA polymorphisms ----

Snp = read_table(file = "9.nex/arle.hap", col_names = FALSE)
# head(Snp)
nrow(Snp) # 97 (First one is the reference and others are sample haplotypes)

# Remove the reference
Snp = Snp %>%
  filter(!X1 == "Hap_1")
nrow(Snp) # 96

# Get sample haplotype codes
Snp$X1 = Snp$X1 %>%
  sub(pattern = "Hap_", replacement = "") %>%
  as.numeric()
head(Snp)

# Convert SNP nucleotides to 0/1 allele dosages
Snp2 = Snp$X2 %>%
  strsplit(split = "")
names(Snp2) = 1:length(Snp2)
Snp2 = as_tibble(Snp2) %>%
  t() %>% as_tibble()
dim(Snp2) # 96 361
Snp2
Snp2 = cbind(2:(nrow(Snp2) + 1), Snp2) %>% as_tibble() # NOTE: Hap names are from 2 to 97
dim(Snp2) # 96 362
nSnp = ncol(Snp2) - 1
colnames(Snp2) = c("Hap", paste("Snp_", 1:nSnp, sep = ""))
Snp2

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
# Snp 319 is tri-allelic!
# Recode and ensure proper coding of multi-allelism
tmp3 = data.frame(matrix(NA, nrow = 96, ncol = nSnp+1))
dim(tmp3) # 96 362
Snp3 = cbind(Snp2$Hap, tmp3) %>% as_tibble()
colnames(Snp3) = c("Hap",
                   paste("Snp_", 1:318, sep = ""),
                   "Snp_319_GA",
                   "Snp_319_GC",
                   paste("Snp_", 320:361, sep = ""))
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

dim(Snp3) # 96 363

# Merge with phenotype data
Test = Snp3$Hap %in% Phe$hap_code
sum(Test) # 96
if (sum(Test) != 96) {
  Snp3$Hap[!Test]
  stop("Something went wrong with haplotype coding")
}
Phe = merge(x = Phe, y = Snp3, by.x = "hap_code", by.y = "Hap")
nrow(Phe) # 7576

# ---- Data: recalculate age at first calving ----

# hist(Phe$dob_kod_teljenja)
Tmp = as.numeric((Phe$datum_teljenja - Phe$datum_rodjenja) / 30)
# plot(Tmp ~ Phe$dob_kod_teljenja)
sd(Tmp)                  # 22.86
sd(Phe$dob_kod_teljenja) # 22.49
Phe$dob_kod_teljenja = Tmp

# ---- Data: scale phenotypes and covariates ----

(kolicina_miljekaMean = mean(Phe$kolicina_mlijeka, na.rm = TRUE)) # 7158.711
(kolicina_miljekaSDev = sd(Phe$kolicina_mlijeka, na.rm = TRUE)) # 1843.216
Phe$kolicina_mlijekaS = (Phe$kolicina_mlijeka - kolicina_miljekaMean) / kolicina_miljekaSDev

(kolicina_m_mMean = mean(Phe$kolicina_m_m, na.rm = TRUE)) # 278.1104
(kolicina_m_mSDev = sd(Phe$kolicina_m_m, na.rm = TRUE)) # 79.11727
Phe$kolicina_m_mS = (Phe$kolicina_m_m - kolicina_m_mMean) / kolicina_m_mSDev

(kolicina_proteinaMean = mean(Phe$kolicina_proteina, na.rm = TRUE)) # 235.1427
(kolicina_proteinaSDev = sd(Phe$kolicina_proteina, na.rm = TRUE)) # 60.49171
Phe$kolicina_proteinaS = (Phe$kolicina_proteina - kolicina_proteinaMean) / kolicina_proteinaSDev

(dob_kod_teljenjaMean = Phe %>%
    group_by(broj_teljenjaC) %>%
    summarise(dob_kod_teljenjaMean = mean(dob_kod_teljenja, na.rm = TRUE)))
(dob_kod_teljenjaSDev = Phe %>%
    group_by(broj_teljenjaC) %>%
    summarise(dob_kod_teljenjaSDev = sd(dob_kod_teljenja, na.rm = TRUE)))
dob_kod_teljenjaStats = merge(x = dob_kod_teljenjaMean, y = dob_kod_teljenjaSDev)
Phe = merge(x = Phe, y = dob_kod_teljenjaStats)
Phe$dob_kod_teljenjaS = (Phe$dob_kod_teljenja - Phe$dob_kod_teljenjaMean) / Phe$dob_kod_teljenjaSDev


#############################################
# Repeatability model with Barrier SPDE ##### 
#############################################

length(unique(Phe$founder))
# 109
length(unique(Phe$hap_code))
# 96
length(unique(Phe$amino_hap_code))
# 48
length(unique(Phe$evohap))
# 10

PheLact15 = Phe %>%
  filter(broj_teljenja < 6) %>%
  droplevels()
nrow(PheLact15) # 7115

PheLact15$stado15_c = PheLact15$stado %>% factor() %>% as.numeric()

# Read coordiates of each herd:
koordinate = read.csv("/Users/ipocrnic/Documents/Roslin Posao/Projekti/2019 cattle_mtDNA/Data_original/2019 VladoDairyCattle/8.Koordinate/2.Dodane_koordinate_Marijinim_kojih_nema/stado_lat_lon_3.csv", header = TRUE)
# koordinate = read.csv("stado_lat_lon_3.csv", header = TRUE)

# Add coordinates to each animal based on the herd:
PheLact15loc = merge(x = PheLact15, y = koordinate, by.x = "stado_c", by.y = "stado_c")

# Add map of Croatia:
mapaHR = getData('GADM', country="HRV", level= 1)


#####   TABLE 1 - paper  #######

# Extract the border from the map
hr.border <- unionSpatialPolygons(mapaHR, rep(1, nrow(mapaHR)))

# Formatting boundary for INLA
hr.bdry <- inla.sp2segment(hr.border)

# 1. Mesh:
locations = cbind(koordinate$long, koordinate$lat)

mesh3 = inla.mesh.2d(loc = locations, boundary = hr.bdry, max.edge = c(0.4, 0.8), cutoff = 0.08, offset = c(0.7, 0.7))
# plot(mesh3, main="3rd mesh"); points(koordinate$long, koordinate$lat, col = "red")

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

barrier = over(hr.border, posTriCent, returnList = T)
# - checking which mesh triangles are inside the barrier area
barrier = unlist(barrier)
barrier.triangles = setdiff(1:tl, barrier)

poly.barrier = inla.barrier.polygon(mesh3, barrier.triangles)
# - the Barrier model's polygon
# - in most cases this should be the same as poly.original
# Ignore the warning messages -> same as https://haakonbakkagit.github.io/btopic107.html 

# plot(poly.barrier, add=T, col='lightblue')
# plot(mesh3, add=T); points(lokacije, col='red', pch=20, cex=0.5)

#Define barrier model:
barrier.model = inla.barrier.pcmatern(mesh3, barrier.triangles = barrier.triangles, prior.range = c(0.3, 0.5), prior.sigma = c(1, 0.1) )

# 2.Projector matrix (A-Martix)
# To connect the measurement locations to the mesh representation
A3 = inla.spde.make.A(mesh = mesh3, loc = cbind(PheLact15loc$long.y, PheLact15loc$lat.y))

# dim(A3)
# 7115 771
# No. Obs. x No. mesh param.

# Create a Stack and model for INLA

# Model matrix, and have 1 included to indicate that we have an intercept in the model
MM = model.matrix(~1 + broj_teljenjaC*dob_kod_teljenjaS, PheLact15loc)

# Intercept is together with the SPDE field in the stack, so it is removed from the factors 
FactorNames =  tail(colnames(MM), ncol(MM)-1)

# Prepare for adding the factors to the stack
FactorStack = paste0( "list(", FactorNames, " = MM[,", 2:(ncol(MM)), "]),", collapse = "" )
FactorStack = ( gsub('.{1}$', '', FactorStack ) )
FactorStack = ( gsub(':', 'xx', FactorStack ) )
FactorStack = paste0("c(", FactorStack, ")")

stackbarrier = inla.stack(data = list(y = PheLact15loc$kolicina_mlijekaS),
                          A = list(A3, 1, 1), tag = "Estimate",
                          effects=list(c(list(Intercept = 1),
                                         list(i = 1:mesh3$n)),
                                       c(list(stado_godina_sezona_teljenja_c = PheLact15loc$stado_godina_sezona_teljenja_c),
                                         list(acode = PheLact15loc$acode),
                                         list(acode_sex = PheLact15loc$acode_sex),
                                         list(acode_perm = PheLact15loc$acode_perm),
                                         list(hap_code = PheLact15loc$hap_code)),
                                       c(eval(parse(text = FactorStack)))
                          ))

FactorNames2 = ( gsub(':', 'xx', FactorNames ) )

# Basic haplotype-based model (Mito haplotype effect is IID) 
ModelHapBarrier = c("0", "Intercept", 
                        "f(stado_godina_sezona_teljenja_c, model = \"iid\")", 
                        "f(acode, model = \"generic0\", Cmatrix = QPedigree)",
                        "f(hap_code, model = \"iid\")", 
                        "f(acode_sex, model = \"generic0\", Cmatrix = QSexsp)", 
                        "f(acode_perm, model = \"iid\")", 
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

# capture.output(summary(FitMilkLact15HapBarrier), file = "Summary_FitMilkLact15HapBarrier.txt")
# write.csv(SummarizeInlaVars(FitMilkLact15HapBarrier), file = "Inla_var_FitMilkLact15HapBarrier.csv")

# Variance for spatial component:
teta1 = FitMilkLact15HapBarrier$internal.marginals.hyperpar[[7]]
# Theta 1 = log(SD)
# To get VAR:
te1 = inla.tmarginal(function(x) (exp(x))^2, teta1)
rm(Samples)
Samples = matrix(data = numeric(), nrow = 10000, ncol = 1)
Samples[, 1] = inla.rmarginal(n = 10000, te1)
mean(Samples)

# Range for spatial component:
teta2 = FitMilkLact15HapBarrier$internal.marginals.hyperpar[[8]]
# Theta 2 = log(range)
# To get range:
te2 = inla.tmarginal(function(x) (exp(x)), teta2)
rm(Samples)
Samples = matrix(data = numeric(), nrow = 10000, ncol = 1)
Samples[, 1] = inla.rmarginal(n = 10000, te2)
mean(Samples)


#########################################################
# Repeatability model with SNP Markers and SPDE #########  
#########################################################

#####   TABLE 1 - paper  #######

# Uses barrier SPDE model defined previously

# Add marker necessities: 
# Create Z matrix with unique haplotypes only:

PheLact15un = PheLact15 %>%
  distinct(hap_code, .keep_all = TRUE)
nrow(PheLact15un) #96
# hcode = PheLact15un$hap_code
# sort(hcode)
# hcode 2-97 (1 was reference and was removed)
SnpTmpun = PheLact15un %>%
  dplyr::select(starts_with("SNP_"), "hap_code")
SnpZun = as.matrix(SnpTmpun)
rownames(SnpZun) = SnpTmpun$hap_code
SnpZun = SnpZun[order(as.numeric(rownames(SnpZun))), -ncol(SnpZun)]
dim(SnpZun) #96 362

###### SNP/Z-model by using stack

# Design matrix to connect those unique haplotypes with (repeated) measurements
# Done by hap_code
# So hap_code #2 is 1. column in M, #3 is 2., etc...

M = sparse.model.matrix( ~ 0 + (as.factor(PheLact15$hap_code)))
dim(M) # 7115 x 96

stacksnpspde = inla.stack(data = list(y = PheLact15$kolicina_mlijekaS),
                          A = list(A3, M, 1, 1), tag = "Estimate",
                          effects=list(c(list(Intercept = 1),
                                         list(i = 1:mesh3$n)), 
                                       c(list(snpovi = 1:ncol(M))), 
                                       c(list(stado_godina_sezona_teljenja_c = PheLact15$stado_godina_sezona_teljenja_c),
                                         list(acode = PheLact15$acode),
                                         list(acode_sex = PheLact15$acode_sex),
                                         list(acode_perm = PheLact15$acode_perm)),
                                       c(eval(parse(text = FactorStack)))
                          ))


# Basic SNP (Z) model:
ModelSNP = c("0", "Intercept", 
             "f(stado_godina_sezona_teljenja_c, model = \"iid\")", 
             "f(acode, model = \"generic0\", Cmatrix = QPedigree)",
             "f(acode_sex, model = \"generic0\", Cmatrix = QSexsp)", 
             "f(acode_perm, model = \"iid\")", 
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

# capture.output(summary(FitMilkLact15ZSstack), file = "Summary_FitMilkLact15ZSstack.txt")
# write.csv(SummarizeInlaVars(FitMilkLact15ZSstack), file = "~/Desktop/Inla_var_FitMilkLact15ZSstack.csv")

# Variance for spatial component:
teta1 = FitMilkLact15ZSstack$internal.marginals.hyperpar[[7]]
# Theta 1 = log(SD)
# To get VAR:
te1 = inla.tmarginal(function(x) (exp(x))^2, teta1)
rm(Samples)
Samples = matrix(data = numeric(), nrow = 10000, ncol = 1)
Samples[, 1] = inla.rmarginal(n = 10000, te1)
mean(Samples) 

# Range for spatial component:
teta2 = FitMilkLact15ZSstack$internal.marginals.hyperpar[[8]]
# Theta 2 = log(range)
# To get range:
te2 = inla.tmarginal(function(x) (exp(x)), teta2)
rm(Samples)
Samples = matrix(data = numeric(), nrow = 10000, ncol = 1)
Samples[, 1] = inla.rmarginal(n = 10000, te2)
mean(Samples) 


######################## 
#### Sample Effects #### 
######################## 


# Sample
sample1 <- inla.posterior.sample(1000, FitMilkLact15ZSstack)

# save.image("/Users/ipocrnic/Desktop/aftermarkersample_new.RData")
# load("/Users/ipocrnic/Desktop/aftermarkersample_new.RData")


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
sol_m = s.eff[,97:458]
# Calculate variance:
mean(matrixStats::rowVars(sol_m)) 
# Equal to var from hyperparam: 
mean(1 / inla.rmarginal(n = 1000, marginal = FitMilkLact15ZSstack$marginals.hyperpar[["Precision for snpovi"]]))

# For the haplotype solutions:
sol_haps = s.eff[,1:96]
# Calculate Var(hap)
mean(matrixStats::rowVars(sol_haps))
# The haplotype solutions have to be from sampling

# Residual has to go from the hyperparameter:
mean(1 / inla.rmarginal(n = 1000, marginal = FitMilkLact15ZSstack$marginals.hyperpar[["Precision for the Gaussian observations"]]))
mean(unlist(lapply(sample1, function(x) 1/(x$hyperpar["Precision for the Gaussian observations"]))))


#########################
# Samples by region  ####
#########################

# SNP-Map:
# snploc = read.csv("SNP_location_region.csv", header = TRUE)
snploc = read.csv("/Users/ipocrnic/Documents/Roslin Posao/Projekti/2019 cattle_mtDNA/Data_original/2019 VladoDairyCattle/10.SNP_location_region/SNP_location_region.csv", header = TRUE)

# 361 not 362 --> because we had one SNP triallelic
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
   for(i in 1:362){
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
   Rmat[4,] = sum(diag(SNPmatrixVar[27:45,27:45])) 
   Rmat[5,] = sum(SNPmatrixVar[46:46,46:46])
   Rmat[6,] = sum(diag(SNPmatrixVar[47:67,47:67])) 
   Rmat[7,] = sum(SNPmatrixVar[68:68,68:68])
   Rmat[8,] = sum(SNPmatrixVar[69:69,69:69])
   Rmat[9,] = sum(diag(SNPmatrixVar[70:91,70:91])) 
   Rmat[10,] = sum(SNPmatrixVar[92:92,92:92]) 
   Rmat[11,] = sum(SNPmatrixVar[93:93,93:93]) 
   Rmat[12,] = sum(diag(SNPmatrixVar[94:118,94:118])) 
   Rmat[13,] = sum(SNPmatrixVar[119:119,119:119]) 
   Rmat[14,] = sum(diag(SNPmatrixVar[120:129,120:129])) 
   Rmat[15,] = sum(diag(SNPmatrixVar[130:135,130:135])) 
   Rmat[16,] = sum(diag(SNPmatrixVar[136:147,136:147])) 
   Rmat[17,] = sum(diag(SNPmatrixVar[148:164,148:164])) 
   Rmat[18,] = sum(diag(SNPmatrixVar[165:171,165:171])) 
   Rmat[19,] = sum(SNPmatrixVar[172:172,172:172]) 
   Rmat[20,] = sum(diag(SNPmatrixVar[173:176,173:176])) 
   Rmat[21,] = sum(diag(SNPmatrixVar[177:211,177:211])) 
   Rmat[22,] = sum(diag(SNPmatrixVar[212:213,212:213])) 
   Rmat[23,] = sum(diag(SNPmatrixVar[214:257,214:257])) 
   Rmat[24,] = sum(diag(SNPmatrixVar[258:273,258:273])) 
   Rmat[25,] = sum(SNPmatrixVar[274:274,274:274])
   Rmat[26,] = sum(diag(SNPmatrixVar[275:296,275:296])) 
   Rmat[27,] = sum(SNPmatrixVar[297:297,297:297])
   Rmat[28,] = sum(diag(SNPmatrixVar[298:299,298:299])) 
   Rmat[29,] = sum(diag(SNPmatrixVar[300:362,300:362])) 
   #
   Hmat[[r]] = list(SNPmatrixVar, totalGeneticVar, totalGenicVar, totalLD, Rmat)
   names(Hmat[[r]]) = c("CoV", "totalGeneticVar", "totalGenicVar", "TotalLD", "RegionGenic")
}

# Check dimensions:
dim(Hmat[[1]]$CoV)
# 362x362
dim(Hmat[[1]]$RegionGenic)
# 29x1

# Take the mean over 1000 samples:
Genetic = lapply(Hmat, function(x) x$totalGeneticVar)
mean(unlist(Genetic)) 
# This corresponds to Var(hap) 

Genic = lapply(Hmat, function(x) x$totalGenicVar)
mean(unlist(Genic)) 

TotLD = lapply(Hmat, function(x) x$TotalLD)
mean(unlist(TotLD)) 

GenicRegion = lapply(Hmat, function(x) x$RegionGenic)
GenicRegion_mean = apply(simplify2array(GenicRegion), c(1,2), mean)

###### Construct covariance matrix per region (Compress presented above by region) ####

rm(RegMat)
RegMat = list()
for(r in 1:1000){
  sol_mm = sol_m[r,]
  SNPmatrix = SnpZun
  for(i in 1:362){
    SNPmatrix[,i] = SNPmatrix[,i] * sol_mm[i]
  }
  Rmat = matrix(, ncol = length(unique(snploc$Region)), nrow = nrow(SNPmatrix)) 
  Rmat[,1] = as.matrix(rowSums(SNPmatrix[,1:12])) 
  Rmat[,2] = as.matrix(rowSums(SNPmatrix[,13:25])) 
  Rmat[,3] = as.matrix(SNPmatrix[,26:26])
  Rmat[,4] = as.matrix(rowSums(SNPmatrix[,27:45])) 
  Rmat[,5] = as.matrix(SNPmatrix[,46:46])
  Rmat[,6] = as.matrix(rowSums(SNPmatrix[,47:67])) 
  Rmat[,7] = as.matrix(SNPmatrix[,68:68])
  Rmat[,8] = as.matrix(SNPmatrix[,69:69])
  Rmat[,9] = as.matrix(rowSums(SNPmatrix[,70:91])) 
  Rmat[,10] = as.matrix(SNPmatrix[,92:92]) 
  Rmat[,11] = as.matrix(SNPmatrix[,93:93]) 
  Rmat[,12] = as.matrix(rowSums(SNPmatrix[,94:118])) 
  Rmat[,13] = as.matrix(SNPmatrix[,119:119]) 
  Rmat[,14] = as.matrix(rowSums(SNPmatrix[,120:129])) 
  Rmat[,15] = as.matrix(rowSums(SNPmatrix[,130:135])) 
  Rmat[,16] = as.matrix(rowSums(SNPmatrix[,136:147])) 
  Rmat[,17] = as.matrix(rowSums(SNPmatrix[,148:164])) 
  Rmat[,18] = as.matrix(rowSums(SNPmatrix[,165:171])) 
  Rmat[,19] = as.matrix(SNPmatrix[,172:172]) 
  Rmat[,20] = as.matrix(rowSums(SNPmatrix[,173:176])) 
  Rmat[,21] = as.matrix(rowSums(SNPmatrix[,177:211])) 
  Rmat[,22] = as.matrix(rowSums(SNPmatrix[,212:213])) 
  Rmat[,23] = as.matrix(rowSums(SNPmatrix[,214:257])) 
  Rmat[,24] = as.matrix(rowSums(SNPmatrix[,258:273])) 
  Rmat[,25] = as.matrix(SNPmatrix[,274:274])
  Rmat[,26] = as.matrix(rowSums(SNPmatrix[,275:296])) 
  Rmat[,27] = as.matrix(SNPmatrix[,297:297])
  Rmat[,28] = as.matrix(rowSums(SNPmatrix[,298:299])) 
  Rmat[,29] = as.matrix(rowSums(SNPmatrix[,300:362])) 
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
mean(unlist(TotGenetic))
# Results comparable to Hmat and Var(hap)

TotLD2 = lapply(RegMat, function(x) x$TotalLD)
mean(unlist(TotLD))
# Results comparable to Hmat

# Genetic variance for each Region:
GeneticVarRegije = NULL
for(i in 1:29){
  ri = lapply(RegMat, function(x) x$CoV[[i,i]])
  rim = mean(unlist(ri))
  GeneticVarRegije = rbind(GeneticVarRegije, rim)
}

# LD between the regions:
rm(LDBetweenRegije)
LDBetweenRegije = list()
for(r in 1:1000){
  ri3 = NULL
  for(i in 1:29){ 
    ri = as.matrix(RegMat[[r]]$CoV[,i])
    ri2 = sum(ri[-i,])
    ri3 = rbind(ri3, ri2)
  }
  LDBetweenRegije[[r]] = ri3
  # This number is already total LD between the regions - no need to 2xLD later
}

LDb_mean = apply(simplify2array(LDBetweenRegije), 1, mean) 

# Prepare table for regions:
# Genetic Var, Genic Var, LD between, LD within:
# Note: LDWithinRegion was calculated from the means not actually sampled!
rm(tablica1)
tablica1 = round(cbind(GeneticVarRegije, GenicRegion_mean, GeneticVarRegije-GenicRegion_mean, LDb_mean), 5)
colnames(tablica1) = c("GeneticVarRegion", "GenicVarRegion", "LDWithinRegion", "LDBetweenRegions") 
row.names(tablica1) = unique(snploc$Region)

tablica2 = rbind(tablica1, colSums(tablica1))
row.names(tablica2)[30] = c("SUM")

#####   TABLE 2 - paper  #######
# write.csv(tablica2, "~/Desktop/Tablica2_milk.csv") 

# Heath-map:
# covM = lapply(RegMat, function(x) x$CoV)
# covM_mean = apply(simplify2array(covM), c(1,2), mean)
# heatmap(covM_mean, scale = "none", Rowv = NA, Colv = NA, symm = TRUE)

# Test if cummulative LD-between regions is different than zero ####

#####   TABLE 2 - paper  ####### 
LDb_P = apply(X = simplify2array(LDBetweenRegije), MARGIN = 1, FUN = ProbDiffFromZero)

write.csv(round(LDb_P,2),"~/Desktop/LDb_MLP.csv", row.names = T)


# Take samples from posterior for all the effects
# Calculate covariances and correlations (Effects and EffectxRegions)

# Calculate proper variances and ratios for SNP model

# Based on INLA posterior sample (as all calculations before)
# sample1 <- inla.posterior.sample(1000, FitMilkLact15ZSstack)
# save.image("/Users/ipocrnic/Desktop/aftermarkersample.RData")
# load("/Users/ipocrnic/Desktop/aftermarkersample.RData")


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

sol_a = SampleEfektPost(FitMilkLact15ZSstack, "acode")
colnames(sol_a) = colnames(sol_a) %>% sub(pattern = "acode:", replacement = "", x = ., fixed = TRUE)

sol_s=SampleEfektPost(FitMilkLact15ZSstack, "acode_sex")
colnames(sol_s) = colnames(sol_s) %>% sub(pattern = "acode_sex:", replacement = "", x = ., fixed = TRUE)

sol_p=SampleEfektPost(FitMilkLact15ZSstack, "acode_perm")
# colnames(sol_p) = colnames(sol_p) %>% sub(pattern = "acode_perm:", replacement = "", x = ., fixed = TRUE)
colnames(sol_p) = FitMilkLact15ZSstack$summary.random[4]$acode_perm$ID

sol_hys=SampleEfektPost(FitMilkLact15ZSstack, "stado_godina_sezona_teljenja_c")
# colnames(sol_hys) = colnames(sol_hys) %>% sub(pattern = "stado_godina_sezona_teljenja_c:", replacement = "", x = ., fixed = TRUE)
colnames(sol_hys) = FitMilkLact15ZSstack$summary.random[1]$stado_godina_sezona_teljenja_c$ID

sol_Z = SampleEfektPost(FitMilkLact15ZSstack, "snpovi")
# colnames(sol_Z) = colnames(sol_Z) %>% sub(pattern = "snpovi:", replacement = "", x = ., fixed = TRUE)

sol_haps =  sol_Z[,1:96]
colnames(sol_haps) = seq(2,97)

sol_m = sol_Z[,97:458]  
colnames(sol_m) = seq(1, 362)

sol_l = SampleEfektPost(FitMilkLact15ZSstack, "i")
# colnames(sol_l) = colnames(sol_l) %>% sub(pattern = "i:", replacement = "", x = ., fixed = TRUE)
sol_ll=NULL
for(r in 1:1000){
  AA = A3
  spatialAllAnim = as.vector(AA %*% sol_l[r,])
  sol_ll = rbind(sol_ll, spatialAllAnim)
}
colnames(sol_ll) = PheLact15loc$acode
# This effects is following order of dataset used to create A and mesh for SPDE models
# That dataset has different order than the other one used in analysis  

# They all have different dimensions
# Create design matrices for each of them (for a,s, and p are the same - I created anyways to be clear)
# Those design matrices will multiply solutions in loop for each replicate

Za = sparse.model.matrix( ~ 0 + (as.factor(PheLact15$acode)))
colnames(Za) = colnames(Za) %>% sub(pattern = "as.factor(PheLact15$acode)", replacement = "", x = ., fixed = TRUE)

Zs = sparse.model.matrix( ~ 0 + (as.factor(PheLact15$acode_sex)))
colnames(Zs) = colnames(Zs) %>% sub(pattern = "as.factor(PheLact15$acode_sex)", replacement = "", x = ., fixed = TRUE)

Zp = sparse.model.matrix( ~ 0 + (as.factor(PheLact15$acode_perm)))
colnames(Zp) = colnames(Zp) %>% sub(pattern = "as.factor(PheLact15$acode_perm)", replacement = "", x = ., fixed = TRUE)

Zhys = sparse.model.matrix( ~ 0 + (as.factor(PheLact15$stado_godina_sezona_teljenja_c)))
colnames(Zhys) = colnames(Zhys) %>% sub(pattern = "as.factor(PheLact15$stado_godina_sezona_teljenja_c)", replacement = "", x = ., fixed = TRUE)

Zhap = sparse.model.matrix( ~ 0 + (as.factor(PheLact15$hap_code)))
colnames(Zhap) = colnames(Zhap) %>% sub(pattern = "as.factor(PheLact15$hap_code)", replacement = "", x = ., fixed = TRUE)

Zloc = Za

# Main Loop - Calculate 1000 5x5 Cov Matrices
rm(OverSamples)
OverSamples = list()
for(r in 1:1000){
  a2 = Za %*% sol_a[r,colnames(Za)]
  s2 = Zs %*% sol_s[r,colnames(Zs)]
  p2 = Zp %*% sol_p[r,colnames(Zp)]
  hys2 = Zhys %*% sol_hys[r,colnames(Zhys)]
  hap2 = Zhap %*% sol_haps[r,colnames(Zhap)]
  loc2 = Zloc %*% sol_ll[r,colnames(Zloc)]
  
  S = cbind(a2@x, s2@x, p2@x, hys2@x, hap2@x, loc2@x)
  Scov = cov(S)
  Scor = cor(S)
  OverSamples[[r]] = list(Scov,Scor)
  names(OverSamples[[r]]) = c("Cov", "Cor")
}

# Covariance between the effects:
Cov_samples = lapply(OverSamples, function(x) x$Cov)
cs_mean = round(apply(simplify2array(Cov_samples), c(1,2), mean), 4)
colnames(cs_mean) = c("a", "s", "pe", "hys", "hap", "loc")
rownames(cs_mean) = c("a", "s", "pe", "hys", "hap", "loc")
cs_mean

cs_sd = round(apply(simplify2array(Cov_samples), c(1,2), sd), 4)
colnames(cs_sd) = c("a", "s", "pe", "hys", "hap", "loc")
rownames(cs_sd) = c("a", "s", "pe", "hys", "hap", "loc")
cs_sd

# Correlation between the effects:
Cor_samples = lapply(OverSamples, function(x) x$Cor)
cr_mean = round(apply(simplify2array(Cor_samples), c(1,2), mean), 4)
colnames(cr_mean) = c("a", "s", "pe", "hys", "hap", "loc")
rownames(cr_mean) = c("a", "s", "pe", "hys", "hap", "loc")
cr_mean

cr_sd = round(apply(simplify2array(Cor_samples), c(1,2), sd), 4)
colnames(cr_sd) = c("a", "s", "pe", "hys", "hap", "loc")
rownames(cr_sd) = c("a", "s", "pe", "hys", "hap", "loc")
cr_sd


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

#####   TABLE 1 - paper  ####### correct variances for SNP model

# Variances:
# Equivalent to: SummarizeInlaVars(FitFatLact15ZSstack)
round(mean(unlist(lapply(OverSamples, function(x) x$Ve))), 4)
round(mean(unlist(lapply(OverSamples, function(x) x$Va))), 4)
round(mean(unlist(lapply(OverSamples, function(x) x$Vs))), 4)
round(mean(unlist(lapply(OverSamples, function(x) x$Vp))), 4)
round(mean(unlist(lapply(OverSamples, function(x) x$Vhys))), 4)
round(mean(unlist(lapply(OverSamples, function(x) x$Vhap))), 4)
round(mean(unlist(lapply(OverSamples, function(x) x$Vloc))), 4)

round(sd(unlist(lapply(OverSamples, function(x) x$Ve))), 4)
round(sd(unlist(lapply(OverSamples, function(x) x$Va))), 4)
round(sd(unlist(lapply(OverSamples, function(x) x$Vs))), 4)
round(sd(unlist(lapply(OverSamples, function(x) x$Vp))), 4)
round(sd(unlist(lapply(OverSamples, function(x) x$Vhys))), 4)
round(sd(unlist(lapply(OverSamples, function(x) x$Vhap))), 4)
round(sd(unlist(lapply(OverSamples, function(x) x$Vloc))), 4)


# Ratios:
round(mean(unlist(lapply(OverSamples, function(x) x$r2))), 4)
round(mean(unlist(lapply(OverSamples, function(x) x$h2))), 4)
round(mean(unlist(lapply(OverSamples, function(x) x$s2))), 4)
round(mean(unlist(lapply(OverSamples, function(x) x$p2))), 4)
round(mean(unlist(lapply(OverSamples, function(x) x$hys2))), 4)
round(mean(unlist(lapply(OverSamples, function(x) x$m2))), 4)
round(mean(unlist(lapply(OverSamples, function(x) x$loc2))), 4)

round(sd(unlist(lapply(OverSamples, function(x) x$r2))), 4)
round(sd(unlist(lapply(OverSamples, function(x) x$h2))), 4)
round(sd(unlist(lapply(OverSamples, function(x) x$s2))), 4)
round(sd(unlist(lapply(OverSamples, function(x) x$p2))), 4)
round(sd(unlist(lapply(OverSamples, function(x) x$hys2))), 4)
round(sd(unlist(lapply(OverSamples, function(x) x$m2))), 4)
round(sd(unlist(lapply(OverSamples, function(x) x$loc2))), 4)



