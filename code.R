# Resolving pytorch issues ------------------
Sys.setenv(RETICULATE_USE_UV = "FALSE")
Sys.setenv(RETICULATE_PYTHON = "/Users/joaobiosmac/miniforge3/envs/r-sjsdm/bin/python")
library(reticulate)
use_python("/Users/joaobiosmac/miniforge3/envs/r-sjsdm/bin/python", required = TRUE)
py_config()
library(sjSDM)

# Packages ----------------------------
remove.packages("ggtern")
remove.packages("ggplot2")

install.packages("remotes")
remotes::install_version(
  "ggplot2",
  version = "3.5.1",
  repos = "https://cloud.r-project.org"
)

install.packages("ggtern")
library(ggplot2)
library(ggtern)
library(reticulate)
library(tidyverse)
library(readxl)
library(vegan)
library(betapart)
library(adespatial)
library(betareg)
library(cowplot)
library(gstat)
library(reshape2)
library(tidyr)
library(stringr)
library(forcats)
library(ggtext)
library(ggrepel) 
library(ggcorrplot)

# Getting data ready ----------------------------------
data = read_excel("raw_data.xlsx")
com = data[,16:34]
data$`maximum_length (m)` = as.numeric(data$`maximum_length (m)`)
data$`maximum_width (m)` = as.numeric(data$`maximum_width (m)` )
data$`maximum_depth (cm)` = as.numeric(data$`maximum_depth (cm)`)
data$depth_m = data$`maximum_depth (cm)` / 100
V = (2/3) * pi * data$`maximum_length (m)`* data$`maximum_width (m)`  * data$depth_m 
data$Volume <- V
oc = decostand(com, method = "pa")
Occ = oc
environ = data[,c(11:15, 36)]
environ <- as.data.frame(sapply(environ, as.numeric))
environ <- scale(environ)  
Env = environ
coord = data[,9:10]
coord <- as.data.frame(sapply(coord, as.numeric))
SP = coord
Occ = as.matrix(Occ)
Env = as.matrix(Env)
SP = as.matrix(SP)
keep <- colSums(Occ, na.rm = TRUE) > 0 & colSums(Occ, na.rm = TRUE) < nrow(Occ)
Occ2 <- Occ[, keep, drop = FALSE]
apply(Env, 1, function(x) var(x, na.rm = TRUE))
cors <- cor(Env, use = "pairwise.complete.obs")
high <- which(abs(cors) > 0.95 & abs(cors) < 1, arr.ind = TRUE)
cors[high] 
Env_s  <- as.data.frame(scale(Env))
SP_s   <- as.data.frame(scale(SP[, c("lat","long")]))
spatial = sjSDM::linear(data = SP_s, formula = ~ lat + long)
com_m = as.matrix(com)
var_ok <- apply(com_m, 2, function(x) length(unique(x)) > 1)
site_occ <- colSums(com_m > 0)
site_ok  <- site_occ > 1
keep <- var_ok & site_ok
com_m_filt <- com_m[, keep, drop = FALSE]
species_removed <- colnames(com_m)[!keep]
species_kept    <- colnames(com_m_filt)

if ("sf" %in% class(SP)) {
  SP2 <- as.data.frame(st_coordinates(SP))
  colnames(SP2) <- c("long", "lat")
} else {
  SP2 <- as.data.frame(SP)
}

SP2$lat_c  <- as.numeric(scale(SP2$lat))
SP2$long_c <- as.numeric(scale(SP2$long))
SP2$latlong <- scale(SP2$lat_c * SP2$long_c)[,1]

# sJSDM Model ----------------------------------
set.seed(123)
m <- sjSDM(
  Y = com_m_filt,
  env = linear(data = Env_s,  formula = ~ Volume + stream_distance + Temp + DO + pH + porc_forest),
  spatial = linear(data = SP2, formula = ~ lat_c + long_c + latlong),
  family = "nbinom",
  se = TRUE,
  sampling = 5000L,
  verbose = FALSE
)

summary(m)
an = anova(m, verbose = FALSE)
summary(an)
results = internalStructure(an)
plot(results)

# Figure 2 --------------------------------
df_sites   <- results$raws$Sites   |> mutate(group = "Sites")
df_species <- results$raws$Species |> mutate(group = "Species")

df <- bind_rows(df_sites, df_species) |>
  mutate(
    env    = pmax(env, 0),
    spa    = pmax(spa, 0),
    codist = pmax(codist, 0)
  ) |>
  mutate(total = env + spa + codist) |>
  mutate(
    env    = env    / total,
    spa    = spa    / total,
    codist = codist / total
  )

df$category <- NA_character_
df$category[df$group == "Sites"] <- data$category  

sp_names <- colnames(com_m_filt) 

n_sp <- sum(df$group == "Species")

df$species <- NA_character_
df$species[df$group == "Species"] <- sp_names[seq_len(n_sp)]

amphibious_vec <- c(
  "Atlantirivulus peruibensis",
  "Leptopanchax itanhaensis",
  "Callichthys callichthys",
  "Phalloceros reisi",
  "Poecilia reticulata",
  "Synbranchus marmoratus"
)

df <- df |>
  mutate(
    fill_group = case_when(
      group == "Sites" ~ category,  # Natural / Artificial
      group == "Species" & species %in% amphibious_vec ~ "Amphibious lifestyle",
      group == "Species" ~ "No amphibious lifestyle"
    )
  )


fig2 = ggtern(
  data = df,
  aes(x = env, z = spa, y = codist)
) +
  geom_point(
    aes(size = r2 * 100, fill = fill_group),
    shape = 21, color = "black", alpha = 0.8
  ) +
  scale_fill_manual(
    values = c(
      "Natural"                 = "#66c2a5",
      "Artificial"              = "#fc8d62",
      "Amphibious lifestyle"    = "#4575b4",
      "No amphibious lifestyle" = "#d73027"
    ),
    breaks = c(
      "Natural",
      "Artificial",
      "Amphibious lifestyle",
      "No amphibious lifestyle"
    ),
    name = NULL
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(size = 3),
      nrow = 1                      
    ),
    size = guide_legend(order = 2)   
  )+
  scale_size_continuous(
    name   = expression(R^2~"(%)"),
    range  = c(2, 6),
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100)
  ) +
  facet_wrap(~ group) +
  theme_bw() +
  theme_showarrows() +
  labs(
    x = "Environment",
    y = "Space",
    z = "Species associations"
  ) +
  theme(
    legend.position = "bottom",
    panel.background = element_rect(fill = "white"),
    tern.axis.title.T = element_blank(),
    tern.axis.title.L = element_blank(),
    tern.axis.title.R = element_blank()
  )


fig2

ggsave("figure_2.tiff", fig2, width = 9, height = 7, dpi = 600, compression = "lzw")

# r2 values ----------------------------
summary_sites_all <- df |>
  filter(group == "Sites") |>
  summarise(
    r2_mean = mean(r2, na.rm = TRUE),
    r2_sd   = sd(r2, na.rm = TRUE),
    n       = sum(!is.na(r2))
  )

summary_sites_all

summary_sites_by_cat <- df |>
  filter(group == "Sites") |>
  group_by(category) |>
  summarise(
    r2_mean = mean(r2, na.rm = TRUE),
    r2_sd   = sd(r2, na.rm = TRUE),
    n       = sum(!is.na(r2)),
    .groups = "drop"
  )

summary_sites_by_cat

df$period <- NA_character_
df$period[df$group == "Sites"] <- as.character(data$period)

summary_sites_by_period <- df |>
  filter(group == "Sites") |>
  group_by(period) |>
  summarise(
    r2_mean = mean(r2, na.rm = TRUE),
    r2_sd   = sd(r2, na.rm = TRUE),
    n       = sum(!is.na(r2)),
    .groups = "drop"
  )

summary_sites_by_period

summary_species_all <- df |>
  filter(group == "Species") |>
  summarise(
    r2_mean = mean(r2, na.rm = TRUE),
    r2_sd   = sd(r2, na.rm = TRUE),
    n       = sum(!is.na(r2))
  )

summary_species_all

df <- df |>
  mutate(
    lifestyle = case_when(
      group == "Species" & species %in% amphibious_vec      ~ "Amphibious",
      group == "Species" & !species %in% amphibious_vec     ~ "Non-amphibious",
      TRUE                                                   ~ NA_character_
    )
  )

summary_species_by_lifestyle <- df |>
  filter(group == "Species") |>
  group_by(lifestyle) |>
  summarise(
    r2_mean = mean(r2, na.rm = TRUE),
    r2_sd   = sd(r2, na.rm = TRUE),
    n       = sum(!is.na(r2)),
    .groups = "drop"
  )

summary_species_by_lifestyle

# Figure 3 -----------------------------------------
sp_names_true <- colnames(com_m_filt)
betas <- coef(m)$env[[1]]     
ses_raw <- m$se               
ses     <- t(ses_raw)        

if (!is.null(sp_names_true)) {
  sp_names <- sp_names_true
} else {
  sp_names <- rownames(betas)
  if (is.null(sp_names)) {
    sp_names <- paste0("sp", seq_len(nrow(betas)))
  }
}

var_names <- colnames(betas)
if (is.null(var_names)) {
  
  var_names <- c("(Intercept)", "Volume", "stream_distance", "Temp", "DO", "pH",
                 "porc_forest")
}

rownames(betas) <- sp_names
colnames(betas) <- var_names

rownames(ses)   <- sp_names
colnames(ses)   <- var_names

df_b <- as.data.frame(betas) |>
  mutate(species = sp_names) |>
  pivot_longer(
    cols = all_of(var_names),
    names_to = "term",
    values_to = "estimate"
  )

df_se <- as.data.frame(ses) |>
  mutate(species = sp_names) |>
  pivot_longer(
    cols = all_of(var_names),
    names_to = "term",
    values_to = "se"
  )

coefs <- left_join(df_b, df_se,
                   by = c("species", "term"))

coefs <- coefs |>
  filter(term != "(Intercept)") |>
  mutate(
    ci_low  = estimate - 1.96 * se,
    ci_high = estimate + 1.96 * se,
    z       = estimate / se,
    p       = 2 * pnorm(abs(z), lower.tail = FALSE),
    sig     = p < 0.05
  )

ord_terms <- c("DO", "pH", "stream_distance", "Temp", "Volume", "porc_forest")

coefs <- coefs |>
  mutate(
    term = factor(term, levels = ord_terms),
    species = fct_rev(factor(species, levels = sp_names_true))
  )

fig3 = ggplot(coefs, aes(x = estimate, y = species)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                 height = 0, color = "black") +
  geom_point(aes(fill = sig), size = 2, shape = 21, show.legend = FALSE) +
  scale_fill_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey40")) +
  facet_wrap(
    ~ term,
 #    scales = "free_x",
    labeller = as_labeller(c(
      stream_distance = "Stream distance",
      Temp = "Temperature",
      DO = "DO",
      pH = "pH",
      Volume = "Volume",
      porc_forest = "Forest cover"
    ))
  ) +
  labs(
    x = "Coefficient estimate",
    y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.background = element_rect(fill = "white"),
    axis.text.y = element_text(size = 8, face = "italic")
  )

fig3

ggsave("figure_3.tiff", plot = fig3,
       width = 18, height = 14, units = "cm", dpi = 600, compression = "lzw")

# Figure 4 --------------------------------------------------
## Figure 4A ---------------
cooc <- getCor(m)
dim(cooc)       
head(cooc[, 1:5]) 
rownames(cooc) <- colnames(com_m_filt)
colnames(cooc) <- colnames(com_m_filt)


fig4_A = ggcorrplot(
  cooc,
  hc.order = TRUE,       
  type = "lower",      
  lab = FALSE,  
  outline.col = "black",
  tl.cex = 10,
  tl.col = "black",
  tl.srt = 45,
  colors = c("red", "white", "blue"),
  ggtheme = theme_minimal,
  legend.title = NULL
) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.x = element_markdown(face = "italic", size = 8),
    axis.text.y = element_markdown(face = "italic", size = 8)
  )

fig4_A

## Figure 4B -----------------------------------------
cooc <- getCor(m)
rownames(cooc) <- colnames(com_m_filt)
colnames(cooc) <- colnames(com_m_filt)
dist_mat <- as.dist((1 - cooc) / 2)
hc       <- hclust(dist_mat, method = "complete")
ord      <- hc$order
cooc_ord <- cooc[ord, ord]
vals    <- cooc_ord[upper.tri(cooc_ord) | lower.tri(cooc_ord)]
pos_cut <- quantile(vals[vals > 0], probs = 0.975, na.rm = TRUE)
neg_cut <- quantile(vals[vals < 0], probs = 0.025, na.rm = TRUE)
cooc_extreme <- cooc_ord
keep <- (cooc_ord >= pos_cut | cooc_ord <= neg_cut) &
  (row(cooc_ord) != col(cooc_ord))

cooc_extreme[!keep] <- 0 

fig4_B <- ggcorrplot(
  cooc_extreme,
  hc.order   = FALSE,    
  type       = "lower",
  lab        = FALSE,
  outline.col = "black",
  tl.cex     = 10,
  tl.col     = "black",
  tl.srt     = 45,
  colors     = c("red", "white", "blue"),
  ggtheme    = theme_minimal,
  legend.title = NULL
) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.x = element_markdown(face = "italic", size = 8),
    axis.text.y = element_markdown(face = "italic", size = 8)
  )

fig4_B

fig4_complete = plot_grid(fig4_A, fig4_B, labels = "AUTO", nrow = 2)


ggsave("figure_4.tiff", fig4_complete, width = 6, 
       height = 11, dpi = 600, compression = "lzw")

# Figure S1 ---------------------------------------
df_sites <- results$raws$Sites |>
  mutate(
    env    = pmax(env, 0),
    spa    = pmax(spa, 0),
    codist = pmax(codist, 0)
  ) |>
  mutate(total = env + spa + codist) |>
  mutate(
    env    = env    / total,
    spa    = spa    / total,
    codist = codist / total
  )

df_sites$category <- data$category
df_sites$period   <- data$period

df_sites$category <- factor(df_sites$category, 
                            levels = c("Natural", "Artificial"))

df_sites$period <- factor(df_sites$period,
                          levels = c("Wet Period", "Dry Period"))


fig_pools_period <- ggtern(
  data = df_sites,
  aes(x = env, z = spa, y = codist)
) +
  geom_point(
    aes(
      size  = r2 * 100,
      fill  = category,  
      shape = period    
    ),
    color = "black",
    alpha = 0.8
  ) +
  scale_fill_manual(
    values = c(
      "Natural"    = "#66c2a5",
      "Artificial" = "#fc8d62"
    ),
    name = "Habitat"
  ) +
  scale_shape_manual(
    values = c(
      "Wet Period" = 21,  
      "Dry Period" = 24  
    ),
    name = "Hydrological period"
  ) +
  scale_size_continuous(
    name   = expression(R^2~"(%)"),
    range  = c(2, 6),
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100)
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(
        shape = 21,     
        size  = 3,
        color = "black" 
      ),
      order = 1
    ),
    shape = guide_legend(order = 2),
    size  = guide_legend(order = 3)
  ) +
  theme_bw() +
  theme_showarrows() +
  labs(
    x = "Environment",
    y = "Space",
    z = "Species associations"
  ) +
  theme(
    legend.position   = "bottom",
    panel.background  = element_rect(fill = "white"),
    tern.axis.title.T = element_blank(),
    tern.axis.title.L = element_blank(),
    tern.axis.title.R = element_blank()
  )

fig_pools_period

ggsave("Figure_S1.tiff", fig_pools_period, width = 10, height = 8, dpi = 600, compression = "lzw")

# Figure S2 ----------------------------
amphibious_vec <- c(
  "Atlantirivulus peruibensis",
  "Leptopanchax itanhaensis",
  "Callichthys callichthys",
  "Phalloceros reisi",
  "Poecilia reticulata",
  "Synbranchus marmoratus"
)

df_sp <- results$raws$Species |>
  mutate(group = "Species") |>
  mutate(
    env    = pmax(env, 0),
    spa    = pmax(spa, 0),
    codist = pmax(codist, 0)
  ) |>
  mutate(total = env + spa + codist) |>
  mutate(
    env    = env    / total,
    spa    = spa    / total,
    codist = codist / total
  )


df_sp$species <- colnames(com_m_filt)

df_sp <- df_sp |>
  mutate(
    lifestyle = if_else(
      species %in% amphibious_vec,
      "Amphibious lifestyle" ,
      "No amphibious lifestyle" 
    )
  )

fig_sp <- ggtern(
  data = df_sp,
  aes(x = env, z = spa, y = codist)
) +
  geom_point(
    aes(size = r2 * 100, fill = lifestyle),
    shape = 21,
    color = "black",
    alpha = 0.9
  ) +
  geom_text(
    aes(label = species),
    size = 2,
    vjust = -0.5
  ) +
  scale_fill_manual(
    values = c(
      "Amphibious lifestyle"      = "#4575b4",
      "No amphibious lifestyle"  = "#d73027"
    ),
    name = "Lifestyle"
  ) +
  scale_size_continuous(
    name   = expression(R^2~"(%)"),
    range  = c(2, 6),
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100)
  ) +
  facet_wrap(~ lifestyle) +
  theme_bw() +
  theme_showarrows() +
  labs(
    x = "Environment",
    y = "Space",
    z = "Species associations"
  ) +
  theme(
    legend.position = "bottom",
    panel.background = element_rect(fill = "white"),
    tern.axis.title.T = element_blank(),
    tern.axis.title.L = element_blank(),
    tern.axis.title.R = element_blank()
  )

fig_sp

ggsave("Figure_S2.tiff", fig_sp, width = 9, height = 7, dpi = 600, compression = "lzw")
