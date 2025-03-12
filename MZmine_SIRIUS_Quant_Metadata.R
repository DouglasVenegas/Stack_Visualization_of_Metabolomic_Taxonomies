library(fuzzyjoin)
library(tidyverse)
library(readxl)
library(openxlsx)
library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(stringr)
library(writexl)

# Quant_File
Quant <- read_csv("Quant_base.csv") # change quant file name
Quant$`row ID` <- as.character(Quant$`row ID`)

## Load CANOPUS predictions
Canopus <- read_tsv("canopus_structure_summary.tsv") # change name
Canopus$mappingFeatureId <- as.character(Canopus$mappingFeatureId)


#############
# Merge CANOPUS data with experimental data and filter out unwanted columns
# Note that you are working with superclass. In case you want to work with another taxonomy, you must keep the desired taxonomy (remove the -NPC#xxx in the command) and remove the superclass (add -`NPC#superclass`). 
Quant_canopus_superclass <- Quant %>%
  left_join(Canopus, by = c("row ID" = "mappingFeatureId")) %>%
  select(   -`identified by n=`,-`auto MS2 verify`,-`neutral M mass`, -`partners`, 
-`annotation network number`,  -`best ion`,-`correlation group ID`, -`row CCS`,-`row m/z`,
-`row retention time`, -`row ion mobility`, -`row ion mobility unit` ,-ionMass, -formulaRank, 
-formulaId, -alignedFeatureId, -overallFeatureQuality, -retentionTimeInSeconds,
-retentionTimeInMinutes, -molecularFormula, -adduct, -precursorFormula, -`NPC#pathway`,
-`NPC#class`, -`NPC#pathway Probability`, -`NPC#superclass Probability`, 
-`ClassyFire#most specific class`, -`ClassyFire#most specific class Probability`,
-`ClassyFire#level 5`, -`ClassyFire#level 5 Probability`, -`ClassyFire#superclass probability`, 
-`ClassyFire#all classifications`, -`ClassyFire#class Probability`, -`ClassyFire#class`, 
-`ClassyFire#superclass`, -`ClassyFire#subclass`, -`ClassyFire#subclass Probability`, -`NPC#class Probability`)

# Function to obtain the columns in which a feature is present
get_present_samples <- function(feature_row) {
  present_samples <- names(feature_row)[!is.na(feature_row) & feature_row != 0]
  paste(unique(present_samples), collapse = ", ")
}

# Select relevant columns and prepare the dataframe
Quant_precurssor <- Quant_canopus_superclass %>%
  select(-`NPC#superclass`)

# Create dataframe with results
Feature_per_sample <- data.frame(
  Precursor = Quant_precurssor$`row ID`,
  samples = sapply(apply(Quant_precurssor[, 2:ncol(Quant_precurssor)], 1, get_present_samples), 
                   function(x) paste(sub("\\..*", "", unlist(strsplit(x, ","))), collapse = ",")),
  NPC_superclass = Quant_canopus_superclass$`NPC#superclass`
)

# Separate the `samples` column into individual rows
Feature_per_sample_long <- Feature_per_sample %>%
  separate_rows(samples, sep = ",")%>%
  mutate(samples = str_trim(samples)) 


# Count the amount of precursors per sample.
Feature_sample <- Feature_per_sample_long %>%
  group_by(samples) %>%
  summarise(count = n())


# plot
plot=ggplot(Feature_sample, aes(x = samples, y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = count), vjust = -0.5, color = "black") +
  labs(x = "Sample",
       y = "Number of features") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  annotate("text", x = Inf, y = Inf, label = "A)", 
           hjust = 1.1, vjust = 1.1, size = 6, fontface = "bold", color = "black")

# Create the 'plot' folder if it does not exist
if (!dir.exists("plot")) {
  dir.create("plot")
}

# Save the plot in the 'plot' folder
ggsave("plot/features_sample.png", plot = plot, width = 10, height = 6)

plot


########### Proportion pathway per sample


Quant_canopus_pathway <- Quant %>%
  left_join(Canopus, by = c("row ID" = "mappingFeatureId")) %>%
  select(   -`identified by n=`,-`auto MS2 verify`,-`neutral M mass`, -`partners`, 
            -`annotation network number`,  -`best ion`,-`correlation group ID`, -`row CCS`,-`row m/z`,
            -`row retention time`, -`row ion mobility`, -`row ion mobility unit` ,-ionMass, -formulaRank, 
            -formulaId, -alignedFeatureId, -overallFeatureQuality, -retentionTimeInSeconds,
            -retentionTimeInMinutes, -molecularFormula, -adduct, -precursorFormula, -`NPC#superclass`,
            -`NPC#class`, -`NPC#pathway Probability`, -`NPC#superclass Probability`, 
            -`ClassyFire#most specific class`, -`ClassyFire#most specific class Probability`,
            -`ClassyFire#level 5`, -`ClassyFire#level 5 Probability`, -`ClassyFire#superclass probability`, 
            -`ClassyFire#all classifications`, -`ClassyFire#class Probability`, -`ClassyFire#class`, 
            -`ClassyFire#superclass`, -`ClassyFire#subclass`, -`ClassyFire#subclass Probability`, -`NPC#class Probability`)

# Create dataframe with results
Feature_per_sample_pathway <- data.frame(
  Precursor = Quant_precurssor$`row ID`,
  samples = sapply(apply(Quant_precurssor[, 2:ncol(Quant_precurssor)], 1, get_present_samples), 
                   function(x) paste(sub("\\..*", "", unlist(strsplit(x, ","))), collapse = ",")),  # Limpia cada muestra
  NPC_pathway = Quant_canopus_pathway$`NPC#pathway`
)


# Transforming data to long format for visualization
Feature_per_sample_long_pathway  <- Feature_per_sample_pathway %>%
  separate_rows(samples, sep = ",") %>%
  mutate(samples = trimws(samples)) %>%
  mutate(NPC_pathway = ifelse(NPC_pathway %in% (Feature_per_sample_pathway %>%
                                                        count(NPC_pathway) %>%
                                                        top_n(15, wt = n) %>%
                                                        pull(NPC_pathway)), NPC_pathway, "Others")) %>%
  group_by(samples, NPC_pathway) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(samples) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup() %>%
  mutate(samples = factor(samples, levels = sort(unique(samples))))

# plot
pathway_sample = ggplot(Feature_per_sample_long_pathway, aes(x = samples, y = prop, fill = NPC_pathway)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = custom_color_palette) +
  labs(fill = "Feature Pathway", x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  annotate("text", x = Inf, y = Inf, label = "", 
           hjust = 1.1, vjust = 1.1, size = 6, fontface = "bold", color = "black")

# Save the plot in the 'plot' folder'
ggsave("plot/proportion_pathway_sample.png", plot = pathway_sample, width = 10, height = 6)

pathway_sample



########### Proportion superclass per sample
# Custom color palette
custom_color_palette <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666",
  "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5",
  "#e41a1c", "#8B8878", "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c",
  "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928"
)


# Transforming data to long format for visualization
Feature_per_sample_long_superclass  <- Feature_per_sample_long %>%
  separate_rows(samples, sep = ",") %>%
  mutate(samples = trimws(samples)) %>%
  mutate(NPC_superclass = ifelse(NPC_superclass %in% (Feature_per_sample_long %>%
                                                        count(NPC_superclass) %>%
                                                        top_n(25, wt = n) %>%
                                                        pull(NPC_superclass)), NPC_superclass, "Others")) %>%
  group_by(samples, NPC_superclass) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(samples) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup() %>%
  mutate(samples = factor(samples, levels = sort(unique(samples))))

# plot
superclass_sample = ggplot(Feature_per_sample_long_superclass, aes(x = samples, y = prop, fill = NPC_superclass)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = custom_color_palette) +
  labs(fill = "Feature Superclass", x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  annotate("text", x = Inf, y = Inf, label = "", 
           hjust = 1.1, vjust = 1.1, size = 6, fontface = "bold", color = "black")

# Save the plot in the 'plot' folder'
ggsave("plot/proportion_superclass_sample.png", plot = superclass_sample, width = 10, height = 6)

superclass_sample


  
########################## class
  # custom color pallette
custom_color_palette <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666",
  "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5",
  "#e41a1c", "#8B8878", "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c",
  "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928"
)
  # Note that you are working with superclass. In case you want to work with another taxonomy, you must keep the desired taxonomy (remove the -NPC#xxx in the command) and remove the superclass (add -`NPC#superclass`). 
  Quant_canopus_class <- Quant %>%
    left_join(Canopus, by = c("row ID" = "mappingFeatureId")) %>%
    select(   -`identified by n=`,-`auto MS2 verify`,-`neutral M mass`, -`partners`, 
              -`annotation network number`,  -`best ion`,-`correlation group ID`, -`row CCS`,-`row m/z`,
              -`row retention time`, -`row ion mobility`, -`row ion mobility unit` ,-ionMass, -formulaRank, 
              -formulaId, -alignedFeatureId, -overallFeatureQuality, -retentionTimeInSeconds,
              -retentionTimeInMinutes, -molecularFormula, -adduct, -precursorFormula, -`NPC#superclass`,
              -`NPC#pathway`, -`NPC#pathway Probability`, -`NPC#superclass Probability`, 
              -`ClassyFire#most specific class`, -`ClassyFire#most specific class Probability`,
              -`ClassyFire#level 5`, -`ClassyFire#level 5 Probability`, -`ClassyFire#superclass probability`, 
              -`ClassyFire#all classifications`, -`ClassyFire#class Probability`, -`ClassyFire#class`, 
              -`ClassyFire#superclass`, -`ClassyFire#subclass`, -`ClassyFire#subclass Probability`, -`NPC#class Probability`)
  
  Feature_per_sample_class <- data.frame(
    Precursor = Quant_precurssor$`row ID`,
    samples = sapply(apply(Quant_precurssor[, 2:ncol(Quant_precurssor)], 1, get_present_samples), 
                     function(x) paste(sub("\\..*", "", unlist(strsplit(x, ","))), collapse = ",")),  # Limpia cada muestra
    NPC_class = Quant_canopus_class$`NPC#class`
  )
  
  
  # Transforming data to long format for visualization
  Feature_per_sample_long_class  <- Feature_per_sample_class %>%
    separate_rows(samples, sep = ",") %>%
    mutate(samples = trimws(samples)) %>%
    mutate(NPC_class = ifelse(NPC_class %in% (Feature_per_sample_class %>%
                                                    count(NPC_class) %>%
                                                    top_n(30, wt = n) %>%
                                                    pull(NPC_class)), NPC_class, "Others")) %>%
    group_by(samples, NPC_class) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(samples) %>%
    mutate(prop = count / sum(count)) %>%
    ungroup() %>%
    mutate(samples = factor(samples, levels = sort(unique(samples))))
  
  # plot
  class_sample = ggplot(Feature_per_sample_long_class, aes(x = samples, y = prop, fill = NPC_class)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_manual(values = custom_color_palette) +
    labs(fill = "Feature Class", x = NULL, y = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    annotate("text", x = Inf, y = Inf, label = "", 
             hjust = 1.1, vjust = 1.1, size = 6, fontface = "bold", color = "black")
  
  # Save the plot in the 'plot' folder'
  ggsave("plot/proportion_class_sample.png", plot = class_sample, width = 10, height = 6)
  
  class_sample
