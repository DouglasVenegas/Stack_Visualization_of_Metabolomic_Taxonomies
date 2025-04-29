if (!require(pacman)) install.packages("pacman")
pacman::p_load(svglite, readr, fuzzyjoin, tidyverse, readxl, openxlsx, dplyr, stringr, tidyr, ggplot2, writexl, UpSetR, Cairo)

# Part 1: Load data
## Quant_File
Quant <- read_csv("Quant_base.csv") %>%
  mutate(`row ID` = as.character(`row ID`)) %>%
  select(-c(`row m/z`, `row retention time`, `row ion mobility`, `row ion mobility unit`, 
            `row CCS`, `correlation group ID`, `annotation network number`, `best ion`, 
            `auto MS2 verify`, `identified by n=`, `partners`, `neutral M mass`))

## Load CANOPUS predictions
Canopus <- read_tsv("canopus_structure_summary.tsv") # Change name
Canopus$mappingFeatureId <- as.character(Canopus$mappingFeatureId)

## Load metadata 
Metadata <- read_delim("Metadata.txt")%>%
  select("Filename","ATTRIBUTE_Sample" )

#Part 2:taxonomy selection
Taxonomy <- "ClassyFire#" #"NPC#"

taxonomy_columns <- list(
  "NPC#" = c("pathway", "superclass", "class"),
  "ClassyFire#" = c("superclass", "class", "subclass")
)
categories <- taxonomy_columns[[Taxonomy]]


#Part 3: data description
## Function to obtain the columns in which a feature is present
get_present_samples <- function(feature_row) {
  present_samples <- names(feature_row)[!is.na(feature_row) & feature_row != 0]
  paste(unique(present_samples), collapse = ",")
}

## Create dataframe with results
Feature_per_sample <- data.frame(
  Precursor = Quant$`row ID`,
  samples = sapply(apply(Quant[, 2:ncol(Quant)], 1, get_present_samples), 
                   function(x) paste(sub("\\..*", "", unlist(strsplit(x, ","))), collapse = ",")))

## Separate the `samples` column into individual rows
Feature_per_sample_long <- Feature_per_sample %>%
  separate_rows(samples, sep = ",")%>%
  mutate(samples = str_trim(samples)) 


## Count the amount of precursors per sample.
Feature_sample <- Feature_per_sample_long %>%
  group_by(samples) %>%
  summarise(count = n())


## plot
plot=ggplot(Feature_sample, aes(x = samples, y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = count), vjust = -0.5, color = "black") +
  labs(x = "Sample",
       y = "Number of features") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  annotate("text", x = Inf, y = Inf, label = "", 
           hjust = 1.1, vjust = 1.1, size = 6, fontface = "bold", color = "black")

## Create the 'plot' folder if it does not exist
if (!dir.exists("plot")) {
  dir.create("plot")
}

## Save the plot in the 'plot' folder
ggsave("plot/features_sample.png", plot = plot, width = 10, height = 6)
ggsave("plot/features_sample.svg", plot = plot, width = 10, height = 6, device = "svg")


## Considering replicates 

Quant_metadata <- Quant
colnames(Quant_metadata)[2:ncol(Quant_metadata)] <- sapply(colnames(Quant_metadata)[2:ncol(Quant_metadata)], function(x) 
  paste(sub(" .*", "", x), collapse = ","))
col_names <- colnames(Quant_metadata)[-1]
name_map <- setNames(Metadata$ATTRIBUTE_Sample, Metadata$Filename)
new_col_names <- ifelse(col_names %in% names(name_map), name_map[col_names], col_names)
colnames(Quant_metadata)[-1] <- new_col_names

# Create dataframe with results
Feature_per_sample_Metadata <- data.frame(
  Precursor = Quant_metadata$`row ID`,  
  samples = sapply(apply(Quant_metadata[, 2:ncol(Quant_metadata)], 1, get_present_samples), function(x) paste(x, collapse = ","))
)

# Separate the `samples` column into individual rows
Feature_per_sample_Metadata_long <- Feature_per_sample_Metadata %>%
  separate_rows(samples, sep = ",")%>%
  mutate(samples = str_trim(samples)) 

# Count the amount of precursors per sample.
Feature_sample_Metadata <- Feature_per_sample_Metadata_long %>%
  group_by(samples) %>%
  summarise(count = n())

# plot
Features_per_sample_metadata_plot=ggplot(Feature_sample_Metadata, aes(x = samples, y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = count), vjust = -0.5, color = "black") +
  labs(x = "Sample",
       y = "Number of features") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  annotate("text", x = Inf, y = Inf, label = "", 
           hjust = 1.1, vjust = 1.1, size = 6, fontface = "bold", color = "black")


# Save the plot in the 'plot' folder
ggsave("plot/features_sample_Metadata.png", plot = Features_per_sample_metadata_plot, width = 10, height = 6)
ggsave("plot/features_sample_Metadata.svg", plot = Features_per_sample_metadata_plot, width = 10, height = 6, device = "svg")


# Upsetplot

df_upsetplot <- strsplit(as.character(Feature_per_sample_Metadata$samples), ",")
names(df_upsetplot) <- Feature_per_sample_Metadata$Precursor
unique_samples <- unique(unlist(df_upsetplot))
presence_matrix <- do.call(rbind, lapply(df_upsetplot, function(x) {
  as.numeric(unique_samples %in% x)
}))
colnames(presence_matrix) <- unique_samples
presence_df <- as.data.frame(presence_matrix)


png("plot/upset_plot.png", width = 3000, height = 800, res = 150)
upset(
  presence_df, 
  sets = colnames(presence_df), 
  order.by = "freq", 
  main.bar.color = "#2F4F4F",          
  sets.bar.color = "#68838B",          
  matrix.color = "#8B5742",           
  shade.color = "white",               
  set_size.show = TRUE,               
)
dev.off()

svg("plot/upset_plot.svg", width = 18, height = 8)
upset(
  presence_df, 
  sets = colnames(presence_df), 
  order.by = "freq", 
  main.bar.color = "#2F4F4F",          
  sets.bar.color = "#68838B",          
  matrix.color = "#8B5742",           
  shade.color = "white",               
  set_size.show = TRUE,
  mb.ratio = c(0.7, 0.3),
  set_size.angles = 0
)
dev.off()

print("Now you can check the folder 'plot' for your saved graphics")

#############################################################

# Taxonomy

# Proportion pathway per sample
## Custom color palette
custom_color_palette <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666",
  "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5",
  "#e41a1c", "#8B8878", "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c",
  "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928"
)


##################

Quant_canopus <- Quant %>%
  left_join(Canopus, by = c("row ID" = "mappingFeatureId")) %>% 
  select("row ID", contains("Peak Area"), matches(paste0("^", Taxonomy, "(", paste(taxonomy_columns[[Taxonomy]], collapse = "|"), ")$")))

#samples processing
process_sample_names <- function(x) {
  parts <- unlist(strsplit(x, ","))
  base_names <- sub("\\..*", "", parts)
  paste(base_names, collapse = ",")
}

# category loop
for (category in taxonomy_columns[[Taxonomy]]) {
  current_tax_col <- paste0(Taxonomy, category)
  sample_data <- apply(Quant_canopus %>% select(contains("Peak Area")), 1, get_present_samples)
  processed_samples <- sapply(sample_data, process_sample_names)
  Feature_per_sample <- data.frame(
    Precursor = Quant_canopus$`row ID`,
    samples = processed_samples,
    Category = Quant_canopus[[current_tax_col]]
  ) %>% 
    replace(is.na(.), "Unassigned") %>%
    filter(!is.na(Category))
  if(nrow(Feature_per_sample) == 0) {
    warning(paste("Taxonomy column not found:", category))
    next
  }
  Feature_per_sample_long <- Feature_per_sample %>%
    separate_rows(samples, sep = ",") %>%
    mutate(samples = trimws(samples)) %>%
    mutate(Category = ifelse(Category %in% (Feature_per_sample %>%
                                              count(Category) %>%
                                              top_n(15, wt = n) %>%
                                              pull(Category)), 
                             Category, "Others")) %>%
    group_by(samples, Category) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(samples) %>%
    mutate(prop = count / sum(count)) %>%
    ungroup() %>%
    mutate(samples = factor(samples, levels = sort(unique(samples))))
  p <- ggplot(Feature_per_sample_long, aes(x = samples, y = prop, fill = Category)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_manual(values = custom_color_palette) +
    labs(title = paste("Distribution by", str_to_title(category)),
         subtitle = paste("Taxonomy system:", sub("#", "", Taxonomy)),
         x = NULL, y = "Proportion") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  file_prefix <- paste0("plot/", sub("#", "", Taxonomy), "_", category)
  ggsave(paste0(file_prefix, "_sample.png"), plot = p, width = 10, height = 6)
  ggsave(paste0(file_prefix, "_sample.svg"), plot = p, width = 10, height = 6, device = "svg")
  assign(paste0(sub("#", "", Taxonomy), "_", category, "_sample"), p, envir = .GlobalEnv)
}

##### replicates

Quant_canopus_metadata <- Feature_per_sample_Metadata %>%
  left_join(Canopus, by = c("Precursor" = "mappingFeatureId")) %>% 
  select("Precursor","samples", contains("Peak Area"), matches(paste0("^", Taxonomy, "(", paste(taxonomy_columns[[Taxonomy]], collapse = "|"), ")$")))
process_sample_column <- function(sample_string) {
  samples <- unlist(strsplit(sample_string, " ,"))
  trimws(sub("\\..*", "", samples))
}
for (category in taxonomy_columns[[Taxonomy]]) {
  # Get the full taxonomy column name
  current_tax_col <- paste0(Taxonomy, str_remove(category, "\\$"))  
  # Verify the column exists in the data
  if (!current_tax_col %in% names(Quant_canopus_metadata)) {
    warning(paste("Taxonomy column not found:", current_tax_col))
    next  # Skip to next category if column doesn't exist
  }
  Feature_per_sample <- Quant_canopus_metadata %>%
    select(Precursor, samples, Category = all_of(current_tax_col)) %>%
    mutate(
      samples = map(samples, ~ trimws(unlist(strsplit(.x, split = "\\s*,\\s*")))),
      samples = map(samples, ~ sub("\\..*", "", .x))
    ) %>%
    unnest(samples) %>% 
    filter(!is.na(samples) & samples != "") %>%  
    mutate(Category = replace_na(Category, "Unassigned"))  
  # Calculate top categories (15 most frequent)
  top_categories <- Feature_per_sample %>%
    count(Category, sort = TRUE) %>%
    slice_head(n = 15) %>%
    pull(Category)
  Feature_per_sample_long <- Feature_per_sample %>%
    mutate(
      Category = if_else(Category %in% top_categories, Category, "Others"),
      samples = factor(samples, levels = sort(unique(samples)))  
    ) %>%
    group_by(samples, Category) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(samples) %>%
    mutate(prop = count / sum(count)) %>%
    ungroup()

  p <- ggplot(Feature_per_sample_long, aes(x = samples, y = prop, fill = Category)) +
    geom_col(position = "fill") + 
    scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_manual(values = custom_color_palette) +
    labs(
      title = paste("Distribution by", str_to_title(str_remove(category, "\\$"))),
      subtitle = paste("Taxonomy system:", sub("#", "", Taxonomy)),
      x = NULL, 
      y = "Proportion",
      fill = NULL
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
      legend.position = "right",
      plot.title = element_text(face = "bold")
    )
  base_name <- paste0(sub("#", "", Taxonomy), "_", str_remove(category, "\\$"), "_sample_metadata")
  ggsave(
    paste0("plot/", base_name, ".png"), 
    plot = p, 
    width = 10, 
    height = 6,
    dpi = 300 
  )
  ggsave(
    paste0("plot/", base_name, ".svg"), 
    plot = p, 
    width = 10, 
    height = 6
  )
  assign(paste0(base_name, "_plot"), p, envir = .GlobalEnv)
}




