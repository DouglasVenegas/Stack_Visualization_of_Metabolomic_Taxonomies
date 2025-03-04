# Load packages
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

# Metaboanalyst
Metaboanalyst <- read_csv("Metaboanalyst_base.csv") # change metaboanalyst name
Metaboanalyst_WT_Metadata<- Metaboanalyst[-1, ]
Metaboanalyst_Expanded <- Metaboanalyst_WT_Metadata %>%
  separate(Filename, into = c("FeatureID", "Precursor", "Tr"), sep = "/")

# SIRIUS identifications
Structure_identifications <- read_tsv("structure_identifications_all.tsv")
Structure_identifications <- Structure_identifications %>%
  distinct(mappingFeatureId, .keep_all = TRUE)
Structure_identifications$mappingFeatureId <- as.character(Structure_identifications$mappingFeatureId)

## Load CANOPUS predictions
Canopus <- read_tsv("canopus_structure_summary.tsv") # change name
Canopus$mappingFeatureId <- as.character(Canopus$mappingFeatureId)

# Filter by experimental annotations
Experimental_annotations <- Metaboanalyst_Expanded %>%
  filter(!str_ends(Precursor, "mz"))

# Remove experimental annotations from SIRIUS to keep the experimental annotation as priority
Structure_identifications_AIS <- Structure_identifications%>%
  anti_join(Experimental_annotations, by = c("mappingFeatureId" = "FeatureID"))

# Obtain the list of annotable features with priority for the experimental features
Features_annotated <- c(Experimental_annotations$FeatureID, Structure_identifications_AIS$mappingFeatureId)

# Fiter and merge
Features_annotated  <- Metaboanalyst_Expanded %>%
  filter(FeatureID %in% Features_annotated) %>%
  mutate(smiles_annotation = "smiles") %>%
  left_join(Structure_identifications_AIS, by = c("FeatureID" = "mappingFeatureId")) %>%
  select(-name, -overallFeatureQuality,-alignedFeatureId,-formulaId,-retentionTimeInMinutes, -ConfidenceScoreExact,-ConfidenceScoreApproximate, -SiriusScoreNormalized, -precursorFormula,-structurePerIdRank, -formulaRank, -ZodiacScore, -molecularFormula, -adduct, -InChIkey2D, -InChI, -xlogp, -pubchemids, -`CSI:FingerIDScore`, -links, -dbflags, -ionMass, -retentionTimeInSeconds, -SiriusScore, -smiles_annotation, -Tr)

# Dataframe for SIRIUS annotations with SMILES
Metaboanalyst_smiles <- Metaboanalyst_Expanded %>%
  left_join(Structure_identifications_AIS, by = c("FeatureID" = "mappingFeatureId")) %>%
  mutate(smiles_annotation = "smiles") %>%
  select(-overallFeatureQuality,-alignedFeatureId,-formulaId,-retentionTimeInMinutes, -ConfidenceScoreExact,-ConfidenceScoreApproximate, -SiriusScoreNormalized, -precursorFormula,-structurePerIdRank, -formulaRank, -ZodiacScore, -molecularFormula, -adduct, -InChIkey2D, -InChI, -name, -xlogp, -pubchemids, -`CSI:FingerIDScore`, -links, -dbflags, -ionMass, -retentionTimeInSeconds, -SiriusScore, -smiles_annotation, -Tr)

# Create the 'CSV' folder if it does not exist
if (!dir.exists("CSV")) {
  dir.create("CSV")
}

# Exportar archivos CSV
write.csv(rbind(Metaboanalyst[1, ], Experimental_annotations %>% unite(Filename, FeatureID, Precursor, Tr, sep = "/")), "CSV/Metaboanalyst_Experimental_annotations.csv", row.names = FALSE)
write.csv(rbind(Metaboanalyst[1, ], Features_annotated %>% unite(Filename, FeatureID, smiles, Precursor, sep = "/")), "CSV/Metaboanalyst_features_annotated.csv", row.names = FALSE)
write.csv(rbind(Metaboanalyst[1, ], Metaboanalyst_smiles %>% unite(Filename, FeatureID, smiles, Precursor, sep = "/")), "CSV/Metaboanalyst_exprimental_insilico_annotations.csv", row.names = FALSE)



########

# Merge CANOPUS data with experimental data and filter out unwanted columns
# Note that you are working with superclass. In case you want to work with another taxonomy, you must keep the desired taxonomy (remove the -NPC#xxx in the command) and remove the superclass (add -`NPC#superclass`). 
Metaboanalyst_canopus_superclass <- Metaboanalyst_Expanded %>%
  left_join(Canopus, by = c("FeatureID" = "mappingFeatureId")) %>%
  select(-ionMass, -formulaRank, -formulaId, -alignedFeatureId, -overallFeatureQuality, -retentionTimeInSeconds,-retentionTimeInMinutes, -molecularFormula, -adduct, -precursorFormula, -`NPC#pathway`, -`NPC#class`, -`NPC#pathway Probability`, -`NPC#superclass Probability`, -`ClassyFire#most specific class`, -`ClassyFire#most specific class Probability`, -`ClassyFire#level 5`, -`ClassyFire#level 5 Probability`, -`ClassyFire#superclass probability`, -`ClassyFire#all classifications`, -`ClassyFire#class Probability`, -`ClassyFire#class`, -`ClassyFire#superclass`, -`ClassyFire#subclass`, -`ClassyFire#subclass Probability`, -`NPC#class Probability`, -Tr)

# Function to obtain the columns in which a feature is present
get_present_samples <- function(feature_row) {
  present_samples <- names(feature_row)[!is.na(feature_row) & feature_row != 0]# names obtiene nombres de las columnas
  paste(unique(present_samples), collapse = ", ")
}

# Select relevant columns and prepare the dataframe
Metaboanalyst_precurssor <- Metaboanalyst_canopus_superclass %>%
  select(-FeatureID, -`NPC#superclass`)

# Create dataframe with results
Feature_per_sample <- data.frame(
  Precursor = Metaboanalyst_precurssor$Precursor,
  samples = apply(Metaboanalyst_precurssor[, 2:ncol(Metaboanalyst_precurssor)], 1, get_present_samples),
  NPC_superclass = Metaboanalyst_canopus_superclass$`NPC#superclass`
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


###########
# Custom color palette
custom_color_palette <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666",
  "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5",
  "#e41a1c"
)


# Transforming data to long format for visualization
Feature_per_sample_long_superclass  <- Feature_per_sample_long %>%
  separate_rows(samples, sep = ",") %>%
  mutate(samples = trimws(samples)) %>%
  mutate(NPC_superclass = ifelse(NPC_superclass %in% (Feature_per_sample_long %>%
                                                        count(NPC_superclass) %>%
                                                        top_n(15, wt = n) %>%
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

##########

# custom color pallette
custom_color_palette <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666",
  "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5",
  "#e41a1c", "blue"
)


# Merge CANOPUS data with experimental data and filter out unwanted columns
Metaboanalyst_canopus_class <- Metaboanalyst_Expanded %>%
  left_join(Canopus, by = c("FeatureID" = "mappingFeatureId")) %>%
  select(-ionMass, -formulaRank, -formulaId, -alignedFeatureId, -overallFeatureQuality, -retentionTimeInSeconds,-retentionTimeInMinutes, -molecularFormula, -adduct, -precursorFormula, -`NPC#pathway`, -`NPC#superclass`, -`NPC#pathway Probability`, -`NPC#superclass Probability`, -`ClassyFire#most specific class`, -`ClassyFire#most specific class Probability`, -`ClassyFire#level 5`, -`ClassyFire#level 5 Probability`, -`ClassyFire#superclass probability`, -`ClassyFire#all classifications`, -`ClassyFire#class Probability`, -`ClassyFire#class`, -`ClassyFire#superclass`, -`ClassyFire#subclass`, -`ClassyFire#subclass Probability`, -`NPC#class Probability`, -Tr)

# Select relevant columns and prepare the dataframe
Metaboanalyst_precursor_class <- Metaboanalyst_canopus_class %>%
  select(-FeatureID, -`NPC#class`)

# Create dataframe with results
Feature_per_sample_class <- data.frame(
  Precursor = Metaboanalyst_precursor_class$Precursor,
  samples = apply(Metaboanalyst_precursor_class[, 2:ncol(Metaboanalyst_precursor_class)], 1, get_present_samples),
  NPC_class = Metaboanalyst_canopus_class$`NPC#class`
)


# Transform data to long format for visualization
Feature_per_sample_class_long <- Feature_per_sample_class %>%
  separate_rows(samples, sep = ",") %>%
  mutate(samples = trimws(samples)) %>%
  mutate(NPC_class = ifelse(NPC_class %in% (Feature_per_sample_class %>%
                                              count(NPC_class) %>%
                                              top_n(15, wt = n) %>%
                                              pull(NPC_class)), NPC_class, "Others")) %>%
  group_by(samples, NPC_class) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(samples) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup() %>%
  mutate(samples = factor(samples, levels = sort(unique(samples))))

# Plot

class_sample=ggplot(Feature_per_sample_class_long, aes(x = samples, y = prop, fill = NPC_class)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = custom_color_palette) +
  labs(fill = "Feature class", x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  annotate("text", x = Inf, y = Inf, label = "", 
           hjust = 1.1, vjust = 1.1, size = 6, fontface = "bold", color = "black")

class_sample


####

# Samples names

prefixes <- c("202","203","204","205", "206","207","208","209", "210","211","212","213")
nreplicates<- 2 # features in 2/3 réplicates

Metaboanalyst_long <- Metaboanalyst %>%
  pivot_longer(
    cols = -Filename,  # Excluir la columna Filename
    names_to = "variable",  # Nombre de la columna que contendrá los nombres de las variables originales
    values_to = "Value"    # Nombre de la columna que contendrá los valores
  )
# Create a list for storing clean dataframes
FPT <- list()

# Iterate over each prefix
for (prefix in prefixes) {
  FPT_filtered <- Metaboanalyst_long %>%
    filter(str_starts(variable, prefix))
  
  FPT_filtered_w <- FPT_filtered %>%
    pivot_wider(
      names_from = variable,  
      values_from = Value     
    )
  
  FPT_filtered_w_c <- FPT_filtered_w %>%
    rowwise() %>%
    mutate(
      replicas_presentes = sum(!is.na(c_across(starts_with(prefix))))
    ) %>%
    filter(replicas_presentes >= nreplicates) %>%  # Filters out only those variables present in at least 2 replicates
    ungroup() %>%
    select(-replicas_presentes)  
  
  FPT[[prefix]] <- FPT_filtered_w_c 
}

# Extract the Filename column from each dataframe and combine into a single dataframe.
Filename_FPT <- FPT %>%
  lapply(function(df) df %>% select(Filename)) %>% 
  bind_rows() %>% 
  distinct()

# Join each clean dataframe with the combined dataframe based on Filename
combined_FPT <- Filename_FPT

for (prefix in names(FPT)) {
  df <- FPT[[prefix]]
  combined_FPT <- combined_FPT %>%
    left_join(df, by = "Filename")
}

# Save the combined dataframe in a CSV file
write_csv(combined_FPT, "CSV/2_Metaboanalyts_FPT.csv")

#Prepare an excel sheet showing the samples in which each feature appeared.

df_long <- combined_FPT%>%
  pivot_longer(
    cols = -Filename,  
    names_to = "variable", 
    values_to = "Value"    
  )

df_condensado <- df_long%>%
  filter(!is.na(Value)) %>%
  group_by(Filename) %>%
  summarise(variables_presentes = paste(variable, collapse = "/")) %>%
  ungroup()

write_xlsx(df_condensado, "CSV/2_features_en_muestras_FPT.xlsx")


# Obtain in silico annotations per sample
expanded_data <-combined_FPT%>%
  separate(Filename, into = c("FeatureID", "Precursor", "Tr"), sep = "/")

experimental_data <- expanded_data %>%
  filter(!str_ends(Precursor, "mz"))


# Remove experimental annotations from SIRIUS to keep the experimental annotation as priority
sirius_filtered <- Structure_identifications%>%
  anti_join(Experimental_annotations, by = c("mappingFeatureId" = "FeatureID"))

# Obtain the list of annotable features with priority for the experimental features
Features_annotated <- c(experimental_data$FeatureID, sirius_filtered$mappingFeatureId)

# Filter and merge
annotated_features <- expanded_data %>%
  filter(FeatureID %in% Features_annotated) %>%
  left_join(sirius_filtered, by = c("FeatureID" = "mappingFeatureId")) %>%
  mutate(smiles_annotation = "smiles") %>%
  select(-ionMass, -formulaRank, -formulaId, -alignedFeatureId, 
         -overallFeatureQuality, -retentionTimeInSeconds, -retentionTimeInMinutes, 
         -molecularFormula, -adduct, -precursorFormula, -Tr, 
         -ConfidenceScoreExact, -ConfidenceScoreApproximate, -ZodiacScore,
         -links, -dbflags, -pubchemids, -structurePerIdRank, -SiriusScoreNormalized,-SiriusScore,
         -smiles_annotation, -InChI, -InChIkey2D, -xlogp, -`CSI:FingerIDScore`)


write.csv(rbind(combined_FPT[1, ], annotated_features %>% unite(Filename, FeatureID, smiles, name,Precursor, sep = "/")), "CSV/2_Metaboanalyst_features_anotados_FPT.csv", row.names = FALSE)


Metaboanalyst_features_annoatated_FPT <- read_csv("CSV/2_Metaboanalyst_features_anotados_FPT.csv")

df_long <- Metaboanalyst_features_annoatated_FPT %>%
  pivot_longer(
    cols = -Filename,  
    names_to = "variable",  
    values_to = "Value" )   

df_condensed <- df_long%>%
  filter(!is.na(Value)) %>%
  group_by(Filename) %>%
  summarise(variables_presentes = paste(variable, collapse = "/")) %>%
  ungroup()

df_condensado <- df_condensado %>%
  separate(Filename, into = c("FeatureID", "smiles","Anotación_in_silico", "Anotación_Experimental" ), sep = "/")

write_xlsx(df_condensado, "CSV/2_Feutures_samples_Ainsilico_FPT.xlsx")


result_df <- data.frame(
  Precursor = combined_FPT$Filename,
  samples = apply(combined_FPT[, 2:ncol(combined_FPT)], 1, get_present_samples)
)

df_long <- result_df %>%
  separate_rows(samples, sep = ",")%>%
  mutate(samples = str_trim(samples))


# Precurssor per sample
count_precursors <- df_long %>%
  group_by(samples) %>%
  summarise(count = n())

# plot
plot <- ggplot(count_precursors, aes(x = samples, y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = count), vjust = -0.5, color = "black") +
  labs( x = "Sample",
        y = "Number of features") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate("text", x = Inf, y = Inf, label = "B)", 
           hjust = 1.1, vjust = 1.1, size = 6, fontface = "bold", color = "Black")

# Crear la carpeta 'plot' si no existe
if (!dir.exists("plot")) {
  dir.create("plot")
}

# Guardar el gráfico en la carpeta 'plot'
ggsave("plot/features_sample_FPT.png", plot = plot, width = 10, height = 6)

plot

########
Metaboanalyts_FPT <- read_csv("CSV/2_Metaboanalyts_FPT.csv")
filtered_data <- Metaboanalyts_FPT[-1, ]  
canopus_data <- Canopus
canopus_data$featureId <- as.character(canopus_data$mappingFeatureId)
expanded_data <- filtered_data %>%
  separate(Filename, into = c("FeatureID", "Precursor", "Tr"), sep = "/")
combined_data <- expanded_data %>%
  left_join(canopus_data, by = c("FeatureID" = "featureId")) %>%
  select(-alignedFeatureId, -formulaId, -formulaRank, -ionMass, -mappingFeatureId, -retentionTimeInMinutes, -retentionTimeInSeconds,-molecularFormula, -adduct, -precursorFormula, -`NPC#pathway`, -`NPC#class`, -`NPC#pathway Probability`, -`NPC#superclass Probability`, -`ClassyFire#most specific class`, -`ClassyFire#most specific class Probability`, -`ClassyFire#level 5`, -`ClassyFire#level 5 Probability`, -`ClassyFire#superclass probability`, -`ClassyFire#all classifications`, -`ClassyFire#class Probability`, -`ClassyFire#class`, -`ClassyFire#superclass`, -`ClassyFire#subclass`, -`ClassyFire#subclass Probability`, -`NPC#class Probability`, -Tr)
get_present_samples <- function(feature_row) {
  present_samples <- names(feature_row)[!is.na(feature_row) & feature_row != 0]
  paste(unique(present_samples), collapse = ", ")
}
processed_data <- combined_data %>%
  select(-FeatureID, -`NPC#superclass`)
result_df <- data.frame(
  Precursor = processed_data$Precursor,
  samples = apply(processed_data[, 2:ncol(processed_data)], 1, get_present_samples),
  NPC_superclass = combined_data$`NPC#superclass`
) 
long_format_df <- result_df %>%
  separate_rows(samples, sep = ",") %>%
  mutate(samples = trimws(samples)) %>%
  mutate(NPC_superclass = ifelse(NPC_superclass %in% (result_df %>% count(NPC_superclass) %>%
                                                        top_n(15, wt = n) %>%
                                                        pull(NPC_superclass)), NPC_superclass, "Others")) %>%
  group_by(samples, NPC_superclass) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(samples) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup() %>%
  mutate(samples = factor(samples, levels = sort(unique(samples))))
custom_color_palette <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666",
  "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5",
  "#e41a1c"
)
# Plot
superclass_sample=ggplot(long_format_df, aes(x = samples, y = prop, fill = NPC_superclass)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = custom_color_palette) +
  labs(fill = "Feature Superclass", x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  annotate("text", x = Inf, y = Inf, label = "", 
           hjust = 1.1, vjust = 1.1, size = 6, fontface = "bold", color = "black")
ggsave("plot/Superclass_sample_FPT.png", plot = superclass_sample, width = 10, height = 6)

superclass_sample



############
combined_data <- expanded_data %>%
  left_join(canopus_data, by = c("FeatureID" = "mappingFeatureId")) %>%
  select(-alignedFeatureId, -formulaId, -formulaRank, -ionMass, -retentionTimeInMinutes, -retentionTimeInSeconds,-molecularFormula, -adduct, -precursorFormula, -`NPC#pathway`, -`NPC#superclass`, -`NPC#pathway Probability`, -`NPC#superclass Probability`, -`ClassyFire#most specific class`, -`ClassyFire#most specific class Probability`, -`ClassyFire#level 5`, -`ClassyFire#level 5 Probability`, -`ClassyFire#superclass probability`, -`ClassyFire#all classifications`, -`ClassyFire#class Probability`, -`ClassyFire#class`, -`ClassyFire#superclass`, -`ClassyFire#subclass`, -`ClassyFire#subclass Probability`, -`NPC#class Probability`, -Tr)
get_present_samples <- function(feature_row) {
  present_samples <- names(feature_row)[!is.na(feature_row) & feature_row != 0]
  paste(unique(present_samples), collapse = ", ")
}
processed_data <- combined_data %>%
  select(-FeatureID, -`NPC#class`)
result_df <- data.frame(
  Precursor = processed_data$Precursor,
  samples = apply(processed_data[, 2:ncol(processed_data)], 1, get_present_samples),
  NPC_class = combined_data$`NPC#class`
)
long_format_df <- result_df %>%
  separate_rows(samples, sep = ",") %>%
  mutate(samples = trimws(samples)) %>%
  mutate(NPC_class = ifelse(NPC_class %in% (result_df %>%
                                              count(NPC_class) %>%
                                              top_n(15, wt = n) %>%
                                              pull(NPC_class)), NPC_class, "Others")) %>%
  group_by(samples, NPC_class) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(samples) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup() %>%
  mutate(samples = factor(samples, levels = sort(unique(samples)))) 
custom_color_palette <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666",
  "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5",
  "#e41a1c"
)

# Graficar

class_sample=ggplot(long_format_df, aes(x = samples, y = prop, fill = NPC_class)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = custom_color_palette) +
  labs(fill = "Feature class", x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  annotate("text", x = Inf, y = Inf, label = "B)", 
           hjust = 1.1, vjust = 1.1, size = 6, fontface = "bold", color = "black")
ggsave("plot/Class_sample_FPT.png", plot = class_sample, width = 10, height = 6)


class_sample


