install.packages("devtools")
devtools::install_github("ribailey/gghybrid")
library(gghybrid)
library(coda)
library(data.table)
library(MCMCpack)
install.packages("https://cran.r-project.org/src/contrib/Archive/hzar/hzar_0.2-5.tar.gz",repos = NULL,type = "source")
library(hzar)
library(ggplot2)
library(geosphere)

setwd("PATH/05_aDNA/02_results/gghybrid")

# (0) Prepare input data
PREPARE <- FALSE #if you need to process plink's structure input for gghybrid
# (1) Estimate HI
GGHYBRID <- FALSE #if you want to generate HI using gghybrid
NEUTRAL <- TRUE #else outliers
# (2) Prepare input for geographical cline
# (3) Geogrpahical cline
# (4) Plot

my_colors <- c(
  "E"  = "#7F3B08",
  "GB" = "#B35806",
  "F"  = "#E08214",
  "D"  = "#FDB863",
  "BE" = "#FEE0B6",
  "NL"=  "#FFDD00",
  "HZ"  = "#E78AC3", #uli's sample: austria, eastgermany,northitaly
  "IRE"= "#D6266C",
  "S" =  "#A6CEE3",
  "PL"=  "#2899B2",
  "RUS"= "#084081",
  "COR" = "#66C2A5",
  "IT" = "#66C2A5",
  "B"  = "#117733",
  "ISR"= "#7A2A7E", 
  "IRQ"= "#54278F")

##############################################################################
# (0) Prepare input data #

if(PREPARE==TRUE) {
popdata <- read.table("PATH/05_aDNA/01_angsd_enrichedallfresh_rescaled/pop.txt",header=FALSE) 
colnames(popdata) <- c("country", "sampleID", "age","group")
pop <- popdata %>% mutate(country = ifelse(sampleID == "LSS003_TE","F_grey",country))

#for outliers
dat <- read.table("enrichedallfresh_rescaled_outlier_1111_geno_maxmis_q20_dp3_plink.recode.strct_in", skip = 2)
first_row <- readLines("enrichedallfresh_rescaled_outlier_1111_geno_maxmis_q20_dp3_plink.recode.strct_in",n = 1)

fields <- strsplit(first_row, "\\s+")[[1]]
loci <- fields[3:length(fields)]
INDLABEL    <- pop[,2]
POPID <- paste(pop[,1], pop[,4], sep = "_")
paste(pop[,1], pop[,4], sep = "_")
geno <- dat[,-c(1,2)]
geno[geno == 0] <- NA
geno[geno == 2] <- 0
geno[geno == 1] <- 1
odd  <- geno[, seq(1, ncol(geno), 2)]
even <- geno[, seq(2, ncol(geno), 2)]
colnames(odd)  <- loci
colnames(even) <- loci

#keep only the peak region on scaffold 78 used to make pca on fig 3b = 83 loci (from total 112 sites in fig 3 PCA)
targets <- read.table("SNPIDs_pca.txt",header=FALSE,stringsAsFactors = FALSE)[[1]]
keep <- colnames(odd) %in% targets
odd  <- odd[, keep, drop = FALSE]
even <- even[, keep, drop = FALSE]
out1 <- cbind(INDLABEL, POPID, odd)
out2 <- cbind(INDLABEL, POPID, even)
out <- do.call(rbind, lapply(seq_len(nrow(out1)), function(i) {
    rbind(out1[i, ], out2[i, ])}))
  
#remove anomaly ancient samples
bad <- c("BDG002_TE","KCZ001_TE","DSZ007_TE","DVT017_TE","HOC007_TE")
out <- out[!(out$INDLABEL %in% bad), ]
# write.table(out,file = "gghybrid_scaffold78outliers_input.txt",quote = FALSE,row.names = FALSE,col.names = TRUE,sep = "\t")

#for neutral
dat <- read.table("enrichedallfresh_rescaled_neutral_geno_maxmis_q20_dp3_plink.recode.strct_in", skip = 2)
first_row <- readLines("enrichedallfresh_rescaled_neutral_geno_maxmis_q20_dp3_plink.recode.strct_in",n = 1)

fields <- strsplit(first_row, "\\s+")[[1]]
loci <- fields[3:length(fields)]
INDLABEL    <- pop[,2]
POPID <- paste(pop[,1], pop[,4], sep = "_")
paste(pop[,1], pop[,4], sep = "_")
geno <- dat[,-c(1,2)]
geno[geno == 0] <- NA
geno[geno == 2] <- 0
geno[geno == 1] <- 1
odd  <- geno[, seq(1, ncol(geno), 2)]
even <- geno[, seq(2, ncol(geno), 2)]
colnames(odd)  <- loci
colnames(even) <- loci

out1 <- cbind(INDLABEL, POPID, odd)
out2 <- cbind(INDLABEL, POPID, even)
out <- do.call(rbind, lapply(seq_len(nrow(out1)), function(i) {
  rbind(out1[i, ], out2[i, ])}))

#remove anomaly ancient samples
bad <- c("BDG002_TE","KCZ001_TE","DSZ007_TE","DVT017_TE","HOC007_TE")
out <- out[!(out$INDLABEL %in% bad), ]
# write.table(out,file = "gghybrid_neutral_input.txt",quote = FALSE,row.names = FALSE,col.names = TRUE,sep = "\t")
}

##############################################################################
# (1) Estimate HI #
if (GGHYBRID==TRUE) {
if(NEUTRAL==TRUE) {
dat2 <- read.data("gghybrid_scaffold78outliers_input.txt",nprecol=2,MISSINGVAL=NA, NUMINDS=159, PLOIDY = 2)
prepdata=data.prep(data=dat2$data,loci=dat2$loci,alleles=dat2$alleles,precols=dat2$precols,
                   S0="E_present", S1="IRQ_present", 
                   max.S.MAF = 0.1,	#Filtering by parental minor allele frequency# 
                   return.genotype.table=T, return.locus.table=T)
} else {
dat2 <- read.data("gghybrid_neutral_input.txt",nprecol=2,MISSINGVAL=NA, NUMINDS=159, PLOIDY = 2)
prepdata=data.prep(data=dat2$data,loci=dat2$loci,alleles=dat2$alleles,precols=dat2$precols,
                   S0="E_present",S1="IRQ_present", 
                   max.S.MAF = 0.1,	return.genotype.table=T,return.locus.table=T) 
}

hindlabel=esth(data.prep.object=prepdata$data.prep, read.data.precols=dat2$precols,
  include.Source=TRUE,	                 #Leave at default TRUE if you want hybrid indices for the parental reference individuals, which is often useful#
  plot.col = c("blue","green","cyan","purple","magenta","red"), nitt=3000,burnin=1000)
setkey(hindlabel$hi,beta_mean)

#HI for all samples (in supplementary))
abc = plot_h(data=hindlabel$hi, test.subject=hindlabel$test.subject,
              mean.h.by="POPID",			             #Calculate the mean hybrid index for each value of the "POPID" column#
              sort.by=c("mean_h","POPID","h_posterior_mode"), col.group="POPID",
              group.sep=NULL,fill.source=FALSE,basic.lines=FALSE,
              source.col=c("black","grey"),cex=1,pch=16,cex.lab=1,cex.main=1,ylim=c(0,1))
setkey(abc,rn);
par(xpd = TRUE)
par(mar = c(4, 4, 2, 2)) 
legend("right",inset = c(0, 0),legend = abc[, POPID],
  horiz = FALSE,pch = 16,col = abc[, col.Dark2],
  pt.bg = abc[, col.Dark2],ncol = ceiling(1),
  bty = "n",cex = 0.7) #6x8 landscape

#HI for ancient sampled from 1.5-3k
ancient2k <- c("S0", "S1","E_present","IRQ_present","F_1500-3000", "D_1500-3000", "GB_1500-3000", "RUS_1500-3000", "F_grey_1500-3000", "B_1500-3000", "PL_1500-3000")
present <- c("S0", "S1","E_present","IRQ_present", "D_present", "F_present", "IRE_present", "COR_present","PL_present","S_present","IT_present","B_present","RUS_present","ISR_present")
hi_plot1 <- hindlabel$hi[POPID %in% ancient2k]
abc1 = plot_h(data=hi_plot1,
             test.subject=hindlabel$test.subject,
             mean.h.by="POPID",			             #Calculate the mean hybrid index for each value of the "POPID" column#
             sort.by=c("mean_h","POPID","h_posterior_mode"),  #Order test subjects along the x axis by the mean hybrid index calculated above and also by individual hybrid index ("POPID" is included as some population pairs may have identical mean hi).
             col.group="POPID",
             group.sep=NULL,
             fill.source=FALSE,
             basic.lines=FALSE,
             source.col=c("black","grey"),
             cex=1,pch=16,
             cex.lab=1,cex.main=1,ylim=c(0,1))
setkey(abc1,rn);
par(xpd = TRUE)
par(mar = c(4, 4, 2, 2)) 
legend("right", inset = c(0, 0), legend = abc1[, POPID],
       horiz = FALSE, pch = 16, col = abc1[, col.Dark2], pt.bg = abc1[, col.Dark2],
       ncol = ceiling(1), bty = "n", cex = 0.7)

hi_plot2 <- hindlabel$hi[POPID %in% present]
abc2 = plot_h(data=hi_plot2,
              test.subject=hindlabel$test.subject,
              mean.h.by="POPID",
              sort.by=c("mean_h","POPID","h_posterior_mode"), 
              col.group="POPID",
              group.sep=NULL,
              fill.source=FALSE,
              basic.lines=FALSE,
              source.col=c("black","grey"),
              cex=1,pch=16,
              cex.lab=1,cex.main=1,ylim=c(0,1))
setkey(abc2,rn); par(xpd = TRUE)
par(mar = c(4, 4, 2, 2)) 
legend("right", inset = c(0, 0), legend = abc2[, POPID],
       horiz = FALSE, pch = 16, col = abc2[, col.Dark2],
       pt.bg = abc2[, col.Dark2], ncol = ceiling(1), bty = "n", cex = 0.7)

if(NEUTRAL==TRUE) {
  # fwrite(hindlabel$hi, "HI_table_neutral.tsv")
  } else { 
  #fwrite(hindlabel$hi, "HI_table_outlier.tsv") 
  }
}

##############################################################################
# (2) Prepare input for geographical cline #
if(NEUTRAL==TRUE) {
  hi <- fread("HI_table_neutral.tsv") 
  } else { 
  hi <- fread("HI_table_outlier.tsv") 
}

# Now I add the geo distance info 
loc <- as.data.table(read.table("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/pop_map_more.txt",header=TRUE))
setDT(loc)
loc_unique <- loc[, .SD[1], by = "sample"]
setnames(loc_unique, "sample", "INDLABEL")
loc_unique <- loc_unique[, .(INDLABEL, longitude, latitude)]
hi_merged <- merge(hi,loc_unique,by = "INDLABEL",all.x = TRUE,sort = FALSE)

# Now I add pc1 as clinal axis but this will not be used
pc1_cline <- as.data.table(read.table("/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/01_angsd_enrichedallfresh_rescaled/enrichedallfresh_rescaled_outlier_1111_hapconsensus_maxmis_q20_miss100_smartsnp_withCOR_pc1.txt",header=TRUE))
setnames(pc1_cline, "sample", "INDLABEL")
hi_merged_final <- merge(hi_merged,pc1_cline,by = "INDLABEL",all.x = TRUE)

# European equal-area projection - does not work for the curved transect
## previous projection onto a stright transect ##
#spain <- c(lon = -5.477000 , lat = 42.58800) #endpoint1
#iraq  <- c(lon = 44.361400, lat =  33.94000) #endpoint2
#endpts <- data.frame( lon = c(spain["lon"], iraq["lon"]), lat = c(spain["lat"], iraq["lat"]))
#endpts_sf <- st_as_sf(endpts,coords = c("lon","lat"),crs = 4326)
#endpts_sf <- st_transform(endpts_sf, 3035)
#samples_sf <- st_as_sf(hi_merged_final,coords = c("longitude", "latitude"),crs = 4326)
#samples_sf <- st_transform(samples_sf, 3035)
#A <- st_coordinates(endpts_sf)[1, ]
#B <- st_coordinates(endpts_sf)[2, ]
#P <- st_coordinates(samples_sf)
#AB <- B - A
#t <- ((P[,1]-A[1])*AB[1] + (P[,2]-A[2])*AB[2]) / sum(AB^2)
#transect_length_km <- sqrt(sum(AB^2))/1000
#hi_merged_final$dist_km <- t * transect_length_km
#hi_merged_final$country <- sub("_.*", "", hi_merged_final$POPID) 
#hi_merged_final <- hi_merged_final %>% mutate(country = factor(country, levels = c("E","GB","F","D","BE","NL","IRE","S","PL","RUS","COR","IT","B","ISR","IRQ")))
transect <- st_read("crow_hybrid_zone.kml", quiet = TRUE)
line_coords <- st_coordinates(transect)[, c("X", "Y")]
samples <- hi_merged_final[, c("longitude", "latitude")]
d <- dist2Line(p = samples,line = line_coords)
hi_merged_final$dist_km <- d[, "distance"] / 1000
hi_merged_final$country <- sub("_.*", "", hi_merged_final$POPID)
hi_merged_final <- hi_merged_final %>%
  mutate(country = factor(country,levels = c("E","GB","F","D","BE","NL","IRE","S","PL","RUS","COR","IT","B","ISR","IRQ")))
west <- c("E","GB","F","D","BE","NL")
hi_merged_final$signed_dist_km <-ifelse(hi_merged_final$country %in% west,-hi_merged_final$dist_km,hi_merged_final$dist_km)

#HI vs dist
hy1 <- ggplot(hi_merged_final) +
  geom_point(aes(signed_dist_km,h_posterior_mode,fill = country,shape = period),size = 2,colour = "black") +
  geom_text(data=hi_merged_final%>%filter(POPID %in% c("E_present","IRQ_present")),aes(signed_dist_km,h_posterior_mode,label = POPID),size = 2.5, nudge_y=-0.03)+
  scale_shape_manual(values = c("100-1500" = 22,"1500-3000" = 22,"5000-6000" = 23,"more_than_10k" = 24,"more_than_20k" = 25,"present" = 21)) +
  scale_fill_manual(values= my_colors) +
  guides(fill = guide_legend(override.aes = list(shape = 21,colour = "black"))) +
  theme_minimal() #7x7

# remove Ireland and UK samples as they likely belong to a different hybrid zone transect
# remove non-comparable populations between time 
# separate modern and ancient (1500-3000 ya)

hi_clean <- hi_merged_final[!(POPID %in% c("GB_1500-3000", "IRE_present", "F_grey_1500-3000"))]
hi_present <- hi_clean[period == "present"] %>% filter(!POPID %in% c("COR_present","IT_present","S_present","ISR_present")) #7pops

#exclude RUS
# hi_clean <- hi_merged_final[!(POPID %in% c("GB_1500-3000", "IRE_present", "F_grey_1500-3000", "RUS_1500-3000"))]
# hi_present <- hi_clean[period == "present"] %>% filter(!POPID %in% c("COR_present","IT_present","S_present","ISR_present", "RUS_present")) 
hi_ancient <- hi_clean[period == "1500-3000" | POPID %in% c("E_present", "IRQ_present")] #8 pop as hooded in France as a separate pop =36 inds

#ensure similar sample size per pop and similar locality!!
hi_ancient[, .N, by = POPID]
hi_present[, .N, by = POPID] #RUS to 3, F to 2, D to 7, B to 2, Pl to 1
hi_subsampled <- hi_present %>% filter(INDLABEL %in% c("A07","A08","A09","FPa01","FVa01","D01","D02","D03","D04","D05","D06","D07","B01","B02-b","P01") | POPID %in% c("E_present","IRQ_present"))


if(NEUTRAL==TRUE) {
  #fwrite(rbind(hi_ancient, hi_subsampled), "HI_distance_neutral_AUG2026.tsv")
  #fwrite(rbind(hi_ancient, hi_subsampled), "HI_distance_neutral_AUG2026_noRUS.tsv")
} else { 
  #fwrite(rbind(hi_ancient, hi_subsampled), "HI_distance_outlier_AUG2026.tsv") 
  #fwrite(rbind(hi_ancient, hi_subsampled), "HI_distance_outlier_AUG2026_noRUS.tsv") 
}


##############################################################################
# (3) Geogrpahical cline #

#for pop
hi_present_pop <- hi_present[,.(HI = mean(h_posterior_mode),n = .N, dist = mean(signed_dist_km)),by = POPID]
hi_ancient_pop <-  hi_ancient[,.(HI = mean(h_posterior_mode),n = .N, dist = mean(signed_dist_km)),by = POPID]
hi_subsampled_pop <-  hi_subsampled[,.(HI = mean(h_posterior_mode),n = .N, dist = mean(signed_dist_km)),by = POPID]
setorder(hi_present_pop, dist)
setorder(hi_ancient_pop, dist)
setorder(hi_subsampled_pop, dist)

obsData_pre <- hzar.doMolecularData1DPops(distance = hi_present_pop$dist,
  pObs = hi_present_pop$HI, nEff = hi_present_pop$n)
obsData_anc <- hzar.doMolecularData1DPops(distance = hi_ancient_pop$dist,
                                              pObs = hi_ancient_pop$HI, nEff = hi_ancient_pop$n)
obsData_subsampled <- hzar.doMolecularData1DPops(distance = hi_subsampled_pop$dist,
                                          pObs = hi_subsampled_pop$HI, nEff = hi_subsampled_pop$n)

### analysis for ancient ###
anc <- list()
anc$obs <- obsData_anc
anc$models <- list()

loadModel <- function(obj, scaling, tails, id){
  obj$models[[id]] <- hzar.makeCline1DFreq(
    obj$obs,
    scaling = scaling,
    tails = tails)
  obj
}

anc <- loadModel(anc, "fixed", "none", "modelI")
anc <- loadModel(anc, "free",  "none", "modelII")
anc <- loadModel(anc, "free",  "both", "modelIII")

anc$fitRs <- list()
anc$fitRs$init <- sapply(anc$models, hzar.first.fitRequest.old.ML,
  obsData = anc$obs, verbose = FALSE, simplify = FALSE)

chainLength <- 1000000 #changed from 100000 in outliers as the neutral did not converge
for(m in names(anc$fitRs$init)) {
  anc$fitRs$init[[m]]$mcmcParam$chainLength <- chainLength
  anc$fitRs$init[[m]]$mcmcParam$burnin <-
    chainLength %/% 10
}

## run intial fit for ancient obs
anc$runs$init <- list()
anc$runs$init$modelI <- hzar.doFit(anc$fitRs$init$modelI)
#plot(hzar.mcmc.bindLL(anc$runs$init$modelI))
anc$runs$init$modelII <- hzar.doFit(anc$fitRs$init$modelII)
#plot(hzar.mcmc.bindLL(anc$runs$init$modelII))
anc$runs$init$modelIII <- hzar.doFit(anc$fitRs$init$modelIII)
#plot(hzar.mcmc.bindLL(anc$runs$init$modelIII))

#3 replicates each
anc$fitRs$chains <- lapply(anc$runs$init, hzar.next.fitRequest)
anc$fitRs$chains <- hzar.multiFitRequest(anc$fitRs$chains, each=3, baseSeed=NULL)
anc$runs$chains <- hzar.doChain.multi(anc$fitRs$chains, doPar=TRUE, inOrder=FALSE, count=3)

## Did models converge?
summary(do.call(mcmc.list, lapply(anc$runs$chains[1:3], function(x) hzar.mcmc.bindLL(x[[3]]) )) )
summary(do.call(mcmc.list, lapply(anc$runs$chains[4:6], function(x) hzar.mcmc.bindLL(x[[3]]) )) )
summary(do.call(mcmc.list, lapply(anc$runs$chains[7:9], function(x) hzar.mcmc.bindLL(x[[3]]) )) )
###Gelman-Rubin diagnostic (<1 would be excellent)
mcmc_obj_1 <- do.call(mcmc.list,lapply(anc$runs$chains[1:3],function(x) hzar.mcmc.bindLL(x[[3]])))
mcmc_obj_2 <- do.call(mcmc.list,lapply(anc$runs$chains[4:6],function(x) hzar.mcmc.bindLL(x[[3]])))
mcmc_obj_3 <- do.call(mcmc.list,lapply(anc$runs$chains[7:9],function(x) hzar.mcmc.bindLL(x[[3]])))
gelman.diag(mcmc_obj_1) #1 (out) / 4.1 (neu)
gelman.diag(mcmc_obj_2) #1 (out) / 3.05 (neu)
gelman.diag(mcmc_obj_3) #3.7 (out) / 6.7 (neu)

#anc$fitRs$extra.modelIII <- lapply(lapply(anc$runs$chains[7:9], function(x) x[[3]]), hzar.next.fitRequest)
#anc$runs$extra.modelIII <- hzar.doFit.multi(anc$fitRs$extra.modelIII, doPar=TRUE, inOrder=FALSE)
#summary(do.call(mcmc.list, lapply(anc$runs$extra.modelI, hzar.mcmc.bindLL )))

## Create a model data group for the null model 
anc$analysis$initDGs <- list(nullModel = hzar.dataGroup.null(anc$obs))
## Create a model data group for each initial runs
anc$analysis$initDGs$modelI <- hzar.dataGroup.add(anc$runs$init$modelI)
anc$analysis$initDGs$modelII <- hzar.dataGroup.add(anc$runs$init$modelII)
anc$analysis$initDGs$modelIII <- hzar.dataGroup.add(anc$runs$init$modelIII)
anc$analysis$oDG <- hzar.make.obsDataGroup(anc$analysis$initDGs)
anc$analysis$oDG <- hzar.copyModelLabels(anc$analysis$initDGs, anc$analysis$oDG)
anc$analysis$oDG <- hzar.make.obsDataGroup(lapply(anc$runs$chains, hzar.dataGroup.add), anc$analysis$oDG)
print(summary(anc$analysis$oDG$data.groups))
hzar.plot.cline(anc$analysis$oDG)

## Do model selection based on the AICc scores
print(anc$analysis$AICcTable <- hzar.AICc.hzar.obsDataGroup(anc$analysis$oDG));
print(anc$analysis$model.name <- rownames(anc$analysis$AICcTable)[[ which.min(anc$analysis$AICcTable$AICc )]])
## Extract the hzar.dataGroup object for the selected model: modelI (for both outlier and neutral)
anc$analysis$model.selected <- anc$analysis$oDG$data.groups[[anc$analysis$model.name]]
## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(anc$analysis$model.selected, names(anc$analysis$model.selected$data.param)));
## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(anc$analysis$model.selected))
## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(anc$analysis$model.selected)
## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(anc$analysis$model.selected) #6x6
outancplot <- hzar.plot.fzCline(anc$analysis$model.selected) #6x6

### analysis for subsampled present populations ###
subsampled <- list()
subsampled$obs <- obsData_subsampled
subsampled$models <- list()

loadModel <- function(obj, scaling, tails, id){
  obj$models[[id]] <- hzar.makeCline1DFreq(
    obj$obs,
    scaling = scaling,
    tails = tails)
  obj
}

subsampled <- loadModel(subsampled, "fixed", "none", "modelI")
subsampled <- loadModel(subsampled, "free",  "none", "modelII")
subsampled <- loadModel(subsampled, "free",  "both", "modelIII")

subsampled$fitRs <- list()
subsampled$fitRs$init <- sapply(subsampled$models, hzar.first.fitRequest.old.ML,
                         obsData = subsampled$obs, verbose = FALSE, simplify = FALSE)

chainLength <- 1000000
for(m in names(subsampled$fitRs$init)) {
  subsampled$fitRs$init[[m]]$mcmcParam$chainLength <- chainLength
  subsampled$fitRs$init[[m]]$mcmcParam$burnin <-
    chainLength %/% 10
}

## run intial fit for subsampledient obs
subsampled$runs$init <- list()
subsampled$runs$init$modelI <- hzar.doFit(subsampled$fitRs$init$modelI)
#plot(hzar.mcmc.bindLL(subsampled$runs$init$modelI))
subsampled$runs$init$modelII <- hzar.doFit(subsampled$fitRs$init$modelII)
#plot(hzar.mcmc.bindLL(subsampled$runs$init$modelII))
subsampled$runs$init$modelIII <- hzar.doFit(subsampled$fitRs$init$modelIII)
#plot(hzar.mcmc.bindLL(subsampled$runs$init$modelIII))

#3 replicates each
subsampled$fitRs$chains <- lapply(subsampled$runs$init, hzar.next.fitRequest)
subsampled$fitRs$chains <- hzar.multiFitRequest(subsampled$fitRs$chains, each=3, baseSeed=NULL)
subsampled$runs$chains <- hzar.doChain.multi(subsampled$fitRs$chains, doPar=TRUE, inOrder=FALSE, count=3)

## Did models converge?
summary(do.call(mcmc.list, lapply(subsampled$runs$chains[1:3], function(x) hzar.mcmc.bindLL(x[[3]]) )) )
summary(do.call(mcmc.list, lapply(subsampled$runs$chains[4:6], function(x) hzar.mcmc.bindLL(x[[3]]) )) )
summary(do.call(mcmc.list, lapply(subsampled$runs$chains[7:9], function(x) hzar.mcmc.bindLL(x[[3]]) )) )
###Gelman-Rubin diagnostic (<1 would be excellent)
mcmc_obj_1 <- do.call(mcmc.list,lapply(subsampled$runs$chains[1:3],function(x) hzar.mcmc.bindLL(x[[3]])))
mcmc_obj_2 <- do.call(mcmc.list,lapply(subsampled$runs$chains[4:6],function(x) hzar.mcmc.bindLL(x[[3]])))
mcmc_obj_3 <- do.call(mcmc.list,lapply(subsampled$runs$chains[7:9],function(x) hzar.mcmc.bindLL(x[[3]])))
gelman.diag(mcmc_obj_1) #1 (out) / 3,98 (neu)
gelman.diag(mcmc_obj_2) #1 (out) / 3.84 (neu)
gelman.diag(mcmc_obj_3) #15.3 (out) / 6.61 (neu)

## Create a model data group for the null model 
subsampled$analysis$initDGs <- list(nullModel = hzar.dataGroup.null(subsampled$obs))
## Create a model data group for each initial runs
subsampled$analysis$initDGs$modelI <- hzar.dataGroup.add(subsampled$runs$init$modelI)
subsampled$analysis$initDGs$modelII <- hzar.dataGroup.add(subsampled$runs$init$modelII)
subsampled$analysis$initDGs$modelIII <- hzar.dataGroup.add(subsampled$runs$init$modelIII)
subsampled$analysis$oDG <- hzar.make.obsDataGroup(subsampled$analysis$initDGs)
subsampled$analysis$oDG <- hzar.copyModelLabels(subsampled$analysis$initDGs, subsampled$analysis$oDG)
subsampled$analysis$oDG <- hzar.make.obsDataGroup(lapply(subsampled$runs$chains, hzar.dataGroup.add), subsampled$analysis$oDG)
print(summary(subsampled$analysis$oDG$data.groups))
hzar.plot.cline(subsampled$analysis$oDG)

## Do model selection based on the AICc scores
print(subsampled$analysis$AICcTable <- hzar.AICc.hzar.obsDataGroup(subsampled$analysis$oDG))
print(subsampled$analysis$model.name <- rownames(subsampled$analysis$AICcTable)[[ which.min(subsampled$analysis$AICcTable$AICc )]])
## Extract the hzar.dataGroup object for the selected model: modelI
subsampled$analysis$model.selected <- subsampled$analysis$oDG$data.groups[[subsampled$analysis$model.name]]
## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(subsampled$analysis$model.selected, names(subsampled$analysis$model.selected$data.param)))
## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(subsampled$analysis$model.selected))
## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(subsampled$analysis$model.selected)
## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(subsampled$analysis$model.selected)
outsubsampledplot <- hzar.plot.fzCline(subsampled$analysis$model.selected) #6x6

### analysis for present populations (not needed) ###
pre <- list()
pre$obs <- obsData_pre
pre$models <- list()

loadModel <- function(obj, scaling, tails, id){
  obj$models[[id]] <- hzar.makeCline1DFreq(
    obj$obs,
    scaling = scaling,
    tails = tails)
  obj
}

pre <- loadModel(pre, "fixed", "none", "modelI")
pre <- loadModel(pre, "free",  "none", "modelII")
pre <- loadModel(pre, "free",  "both", "modelIII")

pre$fitRs <- list()
pre$fitRs$init <- sapply(pre$models, hzar.first.fitRequest.old.ML,
                         obsData = pre$obs, verbose = FALSE, simplify = FALSE)

chainLength <- 1000000
for(m in names(pre$fitRs$init)) {
  pre$fitRs$init[[m]]$mcmcParam$chainLength <- chainLength
  pre$fitRs$init[[m]]$mcmcParam$burnin <-
    chainLength %/% 10
}

## run intial fit for preient obs
pre$runs$init <- list()
pre$runs$init$modelI <- hzar.doFit(pre$fitRs$init$modelI)
#plot(hzar.mcmc.bindLL(pre$runs$init$modelI))
pre$runs$init$modelII <- hzar.doFit(pre$fitRs$init$modelII)
#plot(hzar.mcmc.bindLL(pre$runs$init$modelII))
pre$runs$init$modelIII <- hzar.doFit(pre$fitRs$init$modelIII)
#plot(hzar.mcmc.bindLL(pre$runs$init$modelIII))

#3 replicates each
pre$fitRs$chains <- lapply(pre$runs$init, hzar.next.fitRequest)
pre$fitRs$chains <- hzar.multiFitRequest(pre$fitRs$chains, each=3, baseSeed=NULL)
pre$runs$chains <- hzar.doChain.multi(pre$fitRs$chains, doPar=TRUE, inOrder=FALSE, count=3)

## Did models converge?
summary(do.call(mcmc.list, lapply(pre$runs$chains[1:3], function(x) hzar.mcmc.bindLL(x[[3]]) )) )
summary(do.call(mcmc.list, lapply(pre$runs$chains[4:6], function(x) hzar.mcmc.bindLL(x[[3]]) )) )
summary(do.call(mcmc.list, lapply(pre$runs$chains[7:9], function(x) hzar.mcmc.bindLL(x[[3]]) )) )
###Gelman-Rubin diagnostic (<1 would be excellent)
mcmc_obj_1 <- do.call(mcmc.list,lapply(pre$runs$chains[1:3],function(x) hzar.mcmc.bindLL(x[[3]])))
mcmc_obj_2 <- do.call(mcmc.list,lapply(pre$runs$chains[4:6],function(x) hzar.mcmc.bindLL(x[[3]])))
mcmc_obj_3 <- do.call(mcmc.list,lapply(pre$runs$chains[7:9],function(x) hzar.mcmc.bindLL(x[[3]])))
gelman.diag(mcmc_obj_1)
gelman.diag(mcmc_obj_2)
gelman.diag(mcmc_obj_3) 

## Create a model data group for the null model 
pre$analysis$initDGs <- list(nullModel = hzar.dataGroup.null(pre$obs))
## Create a model data group for each initial runs
pre$analysis$initDGs$modelI <- hzar.dataGroup.add(pre$runs$init$modelI)
pre$analysis$initDGs$modelII <- hzar.dataGroup.add(pre$runs$init$modelII)
pre$analysis$initDGs$modelIII <- hzar.dataGroup.add(pre$runs$init$modelIII)
pre$analysis$oDG <- hzar.make.obsDataGroup(pre$analysis$initDGs)
pre$analysis$oDG <- hzar.copyModelLabels(pre$analysis$initDGs, pre$analysis$oDG)
pre$analysis$oDG <- hzar.make.obsDataGroup(lapply(pre$runs$chains, hzar.dataGroup.add), pre$analysis$oDG)
print(summary(pre$analysis$oDG$data.groups))
hzar.plot.cline(pre$analysis$oDG)

## Do model selection based on the AICc scores
print(pre$analysis$AICcTable <- hzar.AICc.hzar.obsDataGroup(pre$analysis$oDG))
print(pre$analysis$model.name <- rownames(pre$analysis$AICcTable)[[ which.min(pre$analysis$AICcTable$AICc )]])
## Extract the hzar.dataGroup object for the selected model: modelII
pre$analysis$model.selected <- pre$analysis$oDG$data.groups[[pre$analysis$model.name]]
## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(pre$analysis$model.selected, names(pre$analysis$model.selected$data.param)))
## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(pre$analysis$model.selected))
## Plot the maximum likelihood cline for the selected model
hzar.plot.cline(pre$analysis$model.selected)
## Plot the 95% credible cline region for the selected model
hzar.plot.fzCline(pre$analysis$model.selected)
outpreplot <- hzar.plot.fzCline(pre$analysis$model.selected) #6x6

# saveRDS(subsampled$analysis$oDG, "hzr_results_outlier_subsampled_cleaned_updatedAUG2026.rds")
# saveRDS(anc$analysis$oDG, "hzr_results_outlier_ancient_cleaned_updatedAUG2026.rds")
# saveRDS(subsampled$analysis$oDG, "hzr_results_neutral_subsampled_cleaned_updatedAUG2026.rds")
# saveRDS(anc$analysis$oDG, "hzr_results_neutral_ancient_cleaned_updatedAUG2026.rds")

# saveRDS(subsampled$analysis$oDG, "hzr_results_outlier_subsampled_cleaned_updatedAUG2026_noRUS.rds")
# saveRDS(anc$analysis$oDG, "hzr_results_outlier_ancient_cleaned_updatedAUG2026_noRUS.rds")
# saveRDS(subsampled$analysis$oDG, "hzr_results_neutral_subsampled_cleaned_updatedAUG2026_noRUS.rds")
# saveRDS(anc$analysis$oDG, "hzr_results_neutral_ancient_cleaned_updatedAUG2026_noRUS.rds")

##############################################################################

# (4) Plot (from previous runs)

if(NEUTRAL==TRUE){
hzr_results_ancient <- readRDS("hzr_results_neutral_ancient_cleaned_updatedAUG2026.rds")
hzr_results_ancient$analysis$oDG <- readRDS("hzr_results_neutral_ancient_cleaned_updatedAUG2026.rds")
hzr_results_subsampled <- readRDS("hzr_results_neutral_subsampled_cleaned_updatedAUG2026.rds")
hzr_results_subsampled$analysis$oDG <- readRDS("hzr_results_neutral_subsampled_cleaned_updatedAUG2026.rds") 
} else{
  hzr_results_ancient <- readRDS("hzr_results_outlier_ancient_cleaned_updatedAUG2026.rds")
hzr_results_ancient$analysis$oDG <- readRDS("hzr_results_outlier_ancient_cleaned_updatedAUG2026.rds")
hzr_results_subsampled <- readRDS("hzr_results_outlier_subsampled_cleaned_updatedAUG2026.rds")
hzr_results_subsampled$analysis$oDG <- readRDS("hzr_results_outlier_subsampled_cleaned_updatedAUG2026.rds") 
}

print(hzr_results_ancient$analysis$AICcTable <- hzar.AICc.hzar.obsDataGroup(hzr_results_ancient$analysis$oDG))
print(hzr_results_ancient$analysis$model.name <- rownames(hzr_results_ancient$analysis$AICcTable)[[ which.min(hzr_results_ancient$analysis$AICcTable$AICc )]])
hzr_results_ancient$analysis$model.selected <- hzr_results_ancient$analysis$oDG$data.groups[[hzr_results_ancient$analysis$model.name]]
center1 <- hzr_results_ancient$analysis$model.selected$ML.cline$param.all$center
width1 <- hzr_results_ancient$analysis$model.selected$ML.cline$param.all$width
left1  <- center1 - width1/2
right1 <- center1 + width1/2

print(hzr_results_subsampled$analysis$AICcTable <- hzar.AICc.hzar.obsDataGroup(hzr_results_subsampled$analysis$oDG))
print(hzr_results_subsampled$analysis$model.name <- rownames(hzr_results_subsampled$analysis$AICcTable)[[ which.min(hzr_results_subsampled$analysis$AICcTable$AICc )]])
hzr_results_subsampled$analysis$model.selected <- hzr_results_subsampled$analysis$oDG$data.groups[[hzr_results_subsampled$analysis$model.name]]
center2 <- hzr_results_subsampled$analysis$model.selected$ML.cline$param.all$center
width2 <- hzr_results_subsampled$analysis$model.selected$ML.cline$param.all$width
left2  <- center2 - width2/2
right2 <- center2 + width2/2

hzar.plot.fzCline(hzr_results_ancient$analysis$model.selected, col="darkorange", alpha=0.1)
hzar.plot.fzCline(hzr_results_subsampled$analysis$model.selected, col="blue",add=TRUE,alpha=0.1)
usr <- par("usr")
segments(x0 = left1,x1 = right1,y0 = usr[3], y1 = usr[3],col = "darkorange",lwd = 5)
abline(v = center1, col = "darkorange", lwd = 2, lty=3)
segments(x0 = left2,x1 = right2,y0 = usr[3], y1 = usr[3],col = "blue",lwd = 5)
abline(v = center2, col = "blue", lwd = 2, lty=3) #5x5
