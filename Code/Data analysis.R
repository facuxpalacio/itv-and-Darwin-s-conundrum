###########################################################################
### TREE ASSEMBLAGES ######################################################
# Packages
library(here) # path reading
library(lme4) # linear mixed models
library(dplyr) # data handling
library(tidyr) # data handling
library(lmerTest) # mixed models p-values
library(ggplot2) # plotting
library(funspace) # trait spaces
library(hypervolume) # hypervolumes
library(ggnewscale) # color palettes

### Tree functional space
trees <- read.table(here("Data", "tree_leaves_talares.txt"), head = T)
trees$LDMC <- 1000*trees$dry_mass_g/trees$wet_mass_g

# Remove sites with mixed or Eucalyptus forests
sites_rem <- c("L17", "L44", "L45", "L46", "L47", "L49", "L50", "L53",
               "T14", "T15", "T16")
trees <- trees[!trees$pointID %in% sites_rem, ]
trees <- trees[-which.max(trees$height), ] # outliers
trees <- trees[-which.min(trees$LDMC), ] # outliers

# Mean number of individuals per species and trait
height_data <- trees[!is.na(trees$height), ]
mean(table(height_data$plant_species))

ldmc_data <- trees[!is.na(trees$LDMC),]
mean(table(ldmc_data$plant_species))

plants <- read.delim(here("Data", "fruit_talares.txt"))
fruit <- plants[!is.na(plants$Peso_fruto_g),]
fruit$inv <- ifelse(fruit$Sp_planta == "Ligustrum_lucidum", "invasive", "native")
nrow(fruit)
fruit$xPeso_fruto <- fruit$Peso_fruto_g/fruit$N_frutos
fruit$xPeso_sem <- fruit$Peso_sem_g/fruit$N_sem
table(fruit$Sp_planta)

fruit %>% group_by(inv) %>% summarise(sd_fruto = sd(xPeso_fruto),
                                      sd_sem = sd(xPeso_sem,na.rm=T))

t.test(fruit$xPeso_fruto ~ fruit$inv)
t.test(fruit$xPeso_sem ~ fruit$inv)

# Collinearity
cor.test(trees$height, trees$LDMC, use = "complete.obs")

# Summary at the plot level
trees %>%
  distinct(plant_species) %>%
  summarise(regional_richness = n())

native_species <- c(
  "Celtis_tala",
  "Scutia_buxifolia",
  "Jodina_rhombifolia",
  "Schinus_longifolius",
  "Phytolacca_dioica",
  "Sapium_haematospermum",
  "Sambucus_australis",
  "Neltuma_alba"
)

trees <- trees %>%
  mutate(origin = ifelse(
    plant_species %in% native_species,
    "Native", "Non-native"
  ))

richness_plot <- trees %>%
  filter(!grepl("^NA", plant_species)) %>%
  distinct(pointID, habitat, plant_species, origin) %>%
  group_by(pointID, habitat, origin) %>%
  summarise(richness = n_distinct(plant_species), .groups = "drop") %>%
  pivot_wider(
    names_from = origin,
    values_from = richness,
    values_fill = 0
  ) %>%
  mutate(total_richness = Native + `Non-native`)

richness_plot %>% group_by(habitat) %>%
  summarize(mean_rich = mean(total_richness),
            sd_rich = sd(total_richness),
            min_rich = min(total_richness),
            max_rich = max(total_richness))

# Exotic species different from Ligustrum lucidum
trees %>%
  filter(
    !plant_species %in% native_species,
    plant_species != "Ligustrum lucidum"
  ) %>%
  distinct(plant_species) %>%
  arrange(plant_species) %>% filter(-c(10:14))

trees %>%
  filter(
    !plant_species %in% native_species,
    plant_species != "Ligustrum lucidum"
  ) %>%
  distinct(habitat, plant_species) %>%
  count(habitat, name = "n_exotic_species")

trees %>%
  filter(
    !plant_species %in% native_species,
    plant_species != "Ligustrum lucidum"
  ) %>%
  distinct(habitat, pointID, plant_species) %>%
  group_by(habitat, plant_species) %>%
  summarise(n_plots = n_distinct(pointID), .groups = "drop")

table(trees$habitat)

trees %>%
  distinct(treeID, .keep_all = TRUE) %>%
  count(habitat, plant_species, name = "n_individuals") %>%
  arrange(habitat, desc(n_individuals))

# one observation per individual
individuals <- trees %>%
  filter(!is.na(treeID), treeID != "") %>%
  distinct(habitat, pointID, plant_species, treeID, .keep_all = TRUE)

# Variability within each plot
variability_plot <- individuals %>%
  group_by(plant_species, pointID) %>%
  summarise(
    across(
      c(height, LDMC),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        SD = ~sd(.x, na.rm = TRUE),
        CV = ~100 * sd(.x, na.rm = TRUE) / mean(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    n = n(),
    .groups = "drop"
  ) 

# Summary between plots for each species
variability_species <- variability_plot %>%
  pivot_longer(
    cols = -c(plant_species, pointID, n),
    names_to = c("trait", ".value"),
    names_pattern = "(height|LDMC)_(mean|SD|CV)"
  ) %>%
  group_by(plant_species, trait) %>%
  summarise(
    n_plots = n(),
    mean = mean(CV, na.rm = TRUE),
    SD = sd(CV, na.rm = TRUE),
    mean_trait = mean(mean, na.rm = TRUE),
    SD_trait = sd(mean, na.rm = TRUE),
    mean_CV = mean(CV, na.rm = TRUE),
    SD_CV = sd(CV, na.rm = TRUE),
    .groups = "drop"
  )

# Boxplots of traits
individuals <- trees %>%
  filter(
    !is.na(treeID),
    treeID != "",
    !is.na(plant_species),
    !grepl("^NA", plant_species)
  ) %>%
  distinct(treeID, .keep_all = TRUE) %>%
  group_by(plant_species) %>%
  filter(n_distinct(treeID) > 1) %>%
  ungroup()

# Order according to median height
species_order_height <- individuals %>%
  group_by(plant_species) %>%
  summarise(mediana = median(height, na.rm = TRUE)) %>%
  arrange(desc(mediana)) %>%
  pull(plant_species)

# Order according to median LDMC
species_order_LDMC <- individuals %>%
  filter(is.finite(LDMC)) %>%
  group_by(plant_species) %>%
  summarise(mediana = median(LDMC)) %>%
  arrange(desc(mediana)) %>%
  pull(plant_species)

p_height <- ggplot(individuals,
                   aes(x = plant_species, y = height,
                       fill = habitat)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    alpha = 0.6,
    outlier.shape = NA
  ) +
  geom_point(
    aes(color = habitat),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.8
    ),
    alpha = 0.35,
    size = 1.5
  ) +
  labs(x = NULL, y = "Plant height (m)",
       fill = "Community type",
       color = "Community type") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  ) + scale_x_discrete(limits = species_order_height)

p_height

p_LDMC <- ggplot(
  filter(individuals, is.finite(LDMC)),
  aes(x = plant_species, y = LDMC, fill = habitat)
)  +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    alpha = 0.6,
    outlier.shape = NA
  ) +
  geom_point(
    aes(color = habitat),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.8
    ),
    alpha = 0.35,
    size = 1.5
  ) +
  labs(x = NULL, y = "Leaf dry matter content (mg/g)",
       fill = "Community type",
       color = "Community type") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )+ scale_x_discrete(limits = species_order_LDMC)

p_LDMC

# Sources of variation
trees$pointID <- as.factor(trees$pointID)
native <- subset(trees, habitat == "native")
invaded <- subset(trees, habitat == "invaded") 

mean(native$height); sd(native$height)
mean(invaded$height); sd(invaded$height)

mean(native$LDMC,na.rm=T); sd(native$LDMC,na.rm=T)
mean(invaded$LDMC,na.rm=T); sd(invaded$LDMC,na.rm=T)

## 1) Count by combination pointID x plant_species
counts <- native %>%
  group_by(pointID, plant_species) %>%
  summarise(n_ind = n(), .groups = "drop")

# 2) Count summary
summary(counts$n_ind)

# 3) % of combinations with >1 individuals
pct_mult <- mean(counts$n_ind > 1)*100
pct_mult # If pct_mult is high (e.g.> 50%), nestedeness can be estimated in general

# 4) check combinations with 1 obs and with multiple obs
counts %>% filter(n_ind == 1) %>% slice_head(n = 10)   # first 10 combs with 1 obs
counts %>% arrange(desc(n_ind)) %>% slice_head(n = 10) # top 10 combs with more reps

# Check native spcies into invaded plots and invasive species into native plots
native_species <- c(
  "Celtis_tala",
  "Scutia_buxifolia",
  "Jodina_rhombifolia",
  "Schinus_longifolius"
)

trees <- trees %>%
  mutate(
    species_origin = case_when(
      plant_species %in% native_species ~ "Native",
      plant_species == "Ligustrum_lucidum" ~ "Ligustrum",
      TRUE ~ "Other_exotic"
    )
  )

composition_habitat <- trees %>%
  count(habitat, species_origin, name = "n_individuals") %>%
  group_by(habitat) %>%
  mutate(
    proportion = n_individuals / sum(n_individuals) * 100
  ) %>%
  ungroup()

composition_habitat

# Variance components
hist(native$height)
hist(invaded$height)
trees$habitat <- relevel(as.factor(trees$habitat), ref = "native")
var_comps_h <- lmer(height ~ habitat + (1|pointID/plant_species), 
                     REML = TRUE, data = trees)
summary(var_comps_h)
vc_h <- VarCorr(var_comps_h)
varcomps_h <- c(unlist(lapply(vc_h, diag)), attr(vc_h, "sc")^2)
varcomps_h <- as.data.frame(t(100*(varcomps_h/sum(varcomps_h))))
colnames(varcomps_h) <- c("Between-species variation (%)", 
                         "Between-plot variation (%)", 
                         "Within-species variation (%)")
varcomps_h

var_comps_l <- lmer(LDMC ~ habitat + (1|pointID/plant_species), 
                    REML = TRUE, data = trees)
summary(var_comps_l)
vc_l <- VarCorr(var_comps_l)
varcomps_l <- c(unlist(lapply(vc_l, diag)), attr(vc_l, "sc")^2)
varcomps_l <- as.data.frame(t(100*(varcomps_l/sum(varcomps_l))))
colnames(varcomps_l) <- c("Between-species variation (%)", 
                          "Between-plot variation (%)", 
                          "Within-species variation (%)")
varcomps_l

# Barplots variance components
varcomps_h$trait <- "Plant height"
varcomps_l$trait <- "LDMC"

varcomps_all <- bind_rows(varcomps_h, varcomps_l)

varcomps_long <- varcomps_all %>%
  pivot_longer(
    cols = -trait,
    names_to = "component",
    values_to = "variance"
  )

ggplot(varcomps_long, aes(x = trait, y = variance, fill = component)) +
  geom_bar(stat = "identity", width = 0.5,
           color = "black",
           linewidth = 0.5) +
  ylab("Variance (%)") +
  xlab("Trait") +
  scale_fill_grey(
    start = 1,   
    end = 0     
  ) +
  theme_minimal()

# Species contribution to ITV
species_contribution <- trees %>%
  filter(
    !is.na(treeID),
    !is.na(plant_species),
    !grepl("^NA", plant_species)
  ) %>%
  distinct(treeID, .keep_all = TRUE) %>%
  pivot_longer(
    cols = c(height, LDMC),
    names_to = "trait",
    values_to = "value"
  ) %>%
  filter(is.finite(value)) %>%
  group_by(habitat, trait, plant_species) %>%
  summarise(
    n = n(),
    variance = var(value,na.rm=T),
    SS_within = (n - 1)*variance,
    .groups = "drop"
  ) %>%
  group_by(habitat, trait) %>%
  mutate(
    contribution_percent = 100*SS_within/ sum(SS_within,na.rm=T)
  ) %>%
  ungroup()

species_contribution

# Trait spaces based on means
# 1. Standardize traits
df_traits <- trees %>%
  mutate(
    height_z = as.numeric(scale(height)),
    LDMC_z   = as.numeric(scale(LDMC))
  )

mean_sp <- trees %>%
  filter(
    !plant_species %in% c(
      "Acacia_sp",
      "Acer_sp",
      "Catalpa_speciosa",
      "Eucalyptus_sp",
      "Fraxinus_sp",
      "Ligustrum_sinense",
      "Celtis_australis",
      "NA1","NA2","NA4","NA5","NA6",
      "Populus_alba",
      "Melia_azedarach",
      "Quercus_sp",
      "Ulmus_populus"
    )
  ) %>%
  group_by(plant_species) %>%
  summarise(
    height = mean(height, na.rm = TRUE),
    LDMC   = mean(LDMC, na.rm = TRUE),
    .groups = "drop"
  )

mean_sp <- mean_sp %>%
  mutate(
    group = ifelse(
      plant_species == "Ligustrum_lucidum",
      "Ligustrum lucidum",
      "Native species"
    ))

ggplot(mean_sp, aes(x = height, y = LDMC)) +
  geom_point(aes(shape = group),size = 4) +
  scale_shape_manual(
    values = c(
      "Native species" = 16,
      "Ligustrum lucidum" = 17
    )
  ) +
  labs(
    x = "Tree height (m)",
    y = "Leaf dry matter content (g/mg)",
    shape = NULL
  ) +
  theme_bw() +
  geom_text(
    aes(label = plant_species),
    vjust = -0.8,
    size = 3.5
  )

# 2. Trait means (natives)
native_species <- df_traits %>%
  filter(
    habitat == "native",
    !plant_species %in% c(
      "Ligustrum_lucidum",
      "Acacia_sp",
      "Celtis_australis",
      "NA1", "NA2", "NA6"
    )) %>%
  group_by(plant_species) %>%
  summarise(
    height_z = mean(height_z, na.rm = TRUE),
    LDMC_z   = mean(LDMC_z, na.rm = TRUE),
    .groups = "drop"
  )

# 3. Ligustrum lucidum mean
ligustrum_mean <- df_traits %>%
  filter(plant_species == "Ligustrum_lucidum") %>%
  summarise(
    height_z = mean(height_z, na.rm = TRUE),
    LDMC_z   = mean(LDMC_z, na.rm = TRUE),
    .groups = "drop"
  )

# 4. Funcional distance between Ligustrum and each native
native_dist <- native_species %>%
  mutate(
    distance_Ligustrum = sqrt(
      (height_z - ligustrum_mean$height_z)^2 +
        (LDMC_z - ligustrum_mean$LDMC_z)^2
    )
  ) %>%
  select(plant_species, distance_Ligustrum)
native_dist

# Null model
species_pool <- bind_rows(
  native_species %>%
    select(plant_species, height_z, LDMC_z),
  ligustrum_mean %>%
    mutate(plant_species = "Ligustrum_lucidum")
)

MPD_observed <- mean(native_dist$distance_Ligustrum,na.rm=T)
MPD_observed

set.seed(123)

n_sim <- 999

MPD_null <- numeric(n_sim)

for (i in 1:n_sim) {
  
  # Randomly select one species as an "invader"
  invader_random <- sample(
    species_pool$plant_species,
    size = 1
  )
  
  # Traits of selected species
  invader_traits <- species_pool %>%
    filter(plant_species == invader_random)
  
  # Pool comparison: remaining species
  native_random <- species_pool %>%
    filter(plant_species != invader_random)
  
  # Distances
  distances <- sqrt(
    (native_random$height_z - invader_traits$height_z)^2 +
      (native_random$LDMC_z - invader_traits$LDMC_z)^2
  )
  
  # MPD
  MPD_null[i] <- mean(distances,na.rm=T)
}

SES_MPD <- (
  MPD_observed - mean(MPD_null,na.rm=T)
) / sd(MPD_null,na.rm=T)
SES_MPD

p_value <- (
  sum(abs(MPD_null - mean(MPD_null,na.rm=T)) >=
        abs(MPD_observed - mean(MPD_null,na.rm=T)),na.rm=T) + 1
) / (n_sim + 1)

p_value

hist(
  MPD_null,
  breaks = 30,
  xlab = "Mean functional distance",
  main = "Null distribution of MPD"
)

abline(
  v = MPD_observed,
  lwd = 3
)

# Trait space accounting for ITV
plot(trees$dry_mass_g, trees$wet_mass_g)
plot(trees$height, trees$LDMC)
cor(trees$height, trees$LDMC, method = "spearman", use = "complete.obs")

trees_na <- na.omit(trees)
nrow(trees_na)
hist(log(trees_na$height))
hist(trees_na$LDMC)

ggplot(data = trees_na, aes(x = height, y = LDMC)) +
  geom_point(aes(col = as.factor(habitat)), size = 2.5)

common_trees <- trees_na %>% 
  filter(plant_species %in% c("Celtis_tala", "Jodina_rhombifolia", "Schinus_longifolius",
                              "Ligustrum_lucidum", "Scutia_buxifolia"))
ggplot(data = common_trees, aes(x = height, y = LDMC)) +
  geom_point(aes(shape = as.factor(plant_species), 
                 col = as.factor(plant_species)), size = 2.5)

zheight <- scale(trees_na$height)
zLDMC <- scale(trees_na$LDMC)

trait_space_habitat <- funspace(data.frame(zheight, zLDMC),
                                n_divisions = 300, 
                                group.vec = as.factor(trees_na$habitat),
                                threshold = 0.999)

summary(trait_space_habitat)


# Plots
plot(x = trait_space_habitat, # funspace object
     type = "global", # plot the global TPD
     quant.plot = TRUE) # add quantile lines

plot(x = trait_space_habitat,
     type = "groups", # a plot for each group (habitat)
     quant.plot = T, # Plot quantile lines
     globalContour = T, # The contour of the global TPD
     pnt = F, # add points for species of each habitat
     pnt.cex = 0.5, # points should be small
     pnt.col = rgb(0, 0, 0, alpha = 0.2), # color for points
     axis.title.line = 1,
     xlim = c(-2,2), ylim = c(-4,4), 
     axis.title.x = "Tree height",
     axis.title.y = "Leaf dry matter content")

# The two hypervolumes in the same plot
grid <- trait_space_habitat$parameters$evaluation_grid

df_invaded <- data.frame(
  grid,
  density = as.vector(
    trait_space_habitat$groups$invaded$images$noThreshold.quantiles
  ),
  quantile = as.vector(
    trait_space_habitat$groups$invaded$images$TPD.quantiles
  ),
  habitat = "Invaded"
)

df_native <- data.frame(
  grid,
  density = as.vector(
    trait_space_habitat$groups$native$images$noThreshold.quantiles
  ),
  quantile = as.vector(
    trait_space_habitat$groups$native$images$TPD.quantiles
  ),
  habitat = "Native"
)

df <- bind_rows(df_invaded, df_native)

ggplot() +
  
  # Native
  geom_raster(
    data = filter(df, habitat == "Native", !is.na(quantile)),
    aes(x = trait1, y = trait2, fill = density),
    alpha = 0.8
  ) +
  scale_fill_gradientn(
    colours = c("#0000FF", "#C6DBEF"),
    name = "Native"
  ) +
  
  ggnewscale::new_scale_fill() +
  
  # Invaded
  geom_raster(
    data = filter(df, habitat == "Invaded", !is.na(quantile)),
    aes(x = trait1, y = trait2, fill = density),
    alpha = 0.5
  ) +
  scale_fill_gradientn(
    colours = c("#FF0000","#FFFF00"),
    name = "Invaded"
  ) +
  
  # Contornos de ambos hábitats
  geom_contour(
    data = df,
    aes(
      x = trait1,
      y = trait2,
      z = quantile,
      color = habitat
    ),
    breaks = c(0.25, 0.50, 0.99),
    linewidth = 0.65,
    na.rm = TRUE
  ) +
  
  scale_color_manual(
    values = c(
      "Invaded" = "orange",
      "Native" = "blue"
    )
  ) +
  
  scale_linetype_manual(
    values = c("Native" = "solid", "Invaded" = "dashed"),
    guide = "none"
  ) +
  
  coord_cartesian(
    xlim = c(-4, 6),
    ylim = c(-4, 4)
  ) +
  
  labs(
    x = "Tree height",
    y = "Leaf dry matter content",
    color = "Habitat"
  ) +
  
  theme_classic(base_size = 13) +
  theme(
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    ),
    axis.line = element_blank()
  ) +
  geom_vline(xintercept=0,linetype="dashed",col="gray40")+
  geom_hline(yintercept=0,linetype="dashed",col="gray40")

# Overlap between hypervolumes
native_traits <- na.omit(trees[trees$habitat == "native", c("height", "LDMC")])
invaded_traits <- na.omit(trees[trees$habitat == "invaded", c("height", "LDMC")])
set.seed(123)
hv_native <- hypervolume(scale(native_traits))
hv_invaded <- hypervolume(scale(invaded_traits))

# Observed statistic
hvlist <- hypervolume_set(hv1 = hv_native, hv2 = hv_invaded,
                          check.memory = F)
hypervolume_overlap_statistics(hvlist)
Jaccard <- hypervolume_overlap_statistics(hvlist)[1]

# Hypothesis testing (permutation test)
# Create a folder named "Objects" and within this, a folder named "trees_perm"
perm_path <- hypervolume_permute("trees_perm", 
                                 hv1 = hv_native, hv2 = hv_invaded, n = 999, cores = 20)
HV_overlap <- hypervolume_overlap_test(hv1 = hv_native, hv2 = hv_invaded, 
                                       path = "Objects/trees_perm", alternative = "one-sided")

HV_overlap$p_values
HV_overlap$plots

#5x4 in
ggplot(data = data.frame(HV_overlap$distribution), aes(x = jaccard)) + 
  geom_histogram(aes(y = after_stat(density)), 
                 color = "black", fill = "white", bins = 30) + 
 #geom_density(alpha = .2, fill = "#FF6666") +
  xlab("Jaccard index") + 
  ylab("Frequency") +
  #xlim(0.6, 0.86) +
  geom_vline(xintercept = Jaccard,
             color = "blue", linetype = "dashed", linewidth = 1) +
  theme_bw()
  
#### Functional richness
trees_na <- na.omit(trees)
zheight <- scale(trees_na$height)
zLDMC <- scale(trees_na$LDMC)
trait_space_habitat <- funspace(data.frame(zheight, zLDMC),
                                n_divisions = 300, 
                                group.vec = as.factor(trees_na$habitat),
                                threshold = 0.999)
FRich <- trait_space_habitat$FD$groups$FRich
set.seed(123)
# Number of iterations
n_iter <- 999

# Store results
FRich_boot <- matrix(NA, nrow = n_iter, ncol = length(unique(trees_na$habitat)))
colnames(FRich_boot) <- levels(as.factor(trees_na$habitat))

for(i in 1:n_iter){
  
  # Subsampling: 100 individuals per habitat type
  samp <- trees_na %>%
    group_by(habitat) %>%
    slice_sample(n = 100) %>%
    ungroup()
  
  # Rescale traits
  zheight <- scale(samp$height)
  zLDMC   <- scale(samp$LDMC)
  
  # Compute functional space
  fs <- funspace(
    data.frame(zheight, zLDMC),
    n_divisions = 300,
    group.vec = as.factor(samp$habitat),
    threshold = 0.999
  )
  
  # Store FRich
  FRich_boot[i, ] <- fs$FD$groups$FRich
}

# Summary
FRich_boot <- read.csv("C:/RD/FRich_boot.csv")
apply(FRich_boot, 2, mean)
apply(FRich_boot, 2, sd)

# Distributions
boxplot(FRich_boot)
t.test(FRich_boot$invaded,FRich_boot$native,var.equal=F)

##########################################################
### SENSITIVITY ANALYSIS
##########################################################
# Prepare data and standardize traits globally
dat <- trees %>%
  filter(complete.cases(height, LDMC)) %>%
  mutate(
    zheight = as.numeric(scale(height)),
    zLDMC = as.numeric(scale(LDMC))
  )

native <- c("Celtis_tala", "Scutia_buxifolia",
            "Jodina_rhombifolia", "Schinus_longifolius")

# Three scenarios
scenarios <- list(
  Original = dat,
  Invaded_without_natives = dat %>%
    filter(!(habitat == "invaded" & plant_species %in% native)),
  Native_without_Ligustrum = dat %>%
    filter(!(habitat == "native" & plant_species == "Ligustrum_lucidum"))
)

# Check sample sizes
sapply(scenarios, function(x) table(x$habitat))

# Calculate Jaccard overlap
set.seed(123)
n_iter <- 999

results <- lapply(scenarios, function(d) {
  
  replicate(n_iter, {
    
    samp <- d %>%
      group_by(habitat) %>%
      slice_sample(n = 100) %>%
      ungroup()
    
    hv_native <- hypervolume(
      samp[samp$habitat == "native", c("zheight", "zLDMC")]
    )
    
    hv_invaded <- hypervolume(
      samp[samp$habitat == "invaded", c("zheight", "zLDMC")]
    )
    
    hv <- hypervolume_set(hv_native, hv_invaded,
                          check.memory = FALSE)
    
    as.numeric(hypervolume_overlap_statistics(hv)[1])
  })
})

# Summary
do.call(rbind, lapply(results, function(x)
  c(mean = mean(x), SD = sd(x),
    lower = quantile(x, 0.025),
    upper = quantile(x, 0.975))
))
