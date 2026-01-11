# Parkinson's Disease Drug Delivery Text Mining Project
# This code searches PubMed for articles about drug delivery in Parkinson's disease
# and creates visualizations to understand the research landscape


# Set where outputs should go (change this to your folder)
OUTPUT_DIR <- "C:/Rprojects/parkinsons_textmining/Outputs"

# Should we search PubMed or just load existing data?
RUN_SEARCH <- TRUE

# How many articles to get from each search
MAX_ARTICLES <- 800

# What years to search
YEAR_MIN <- 1995
YEAR_MAX <- 2025

# Make the output folder if it doesn't exist
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Install and load packages
# List all the packages we need
# I learned about these from various tutorials and Stack Overflow
packages_needed <- c(
  "rentrez",     # for searching PubMed
  "XML",         # for reading XML data from PubMed
  "dplyr",       # for data manipulation
  "tidyr",       # for reshaping data
  "stringr",     # for working with text
  "forcats",     # for working with categories
  "purrr",       # for applying functions
  "tibble",      # for making data frames
  "ggplot2",     # for making plots
  "scales",      # for plot scales
  "tidytext",    # for text mining
  "ggalluvial",  # for flow diagrams
  "igraph",      # for network graphs
  "ggraph",      # for network graphs
  "RColorBrewer", # for colors
  "textmineR",   # for text mining
  "topicmodels", # for topic modeling
  "ggrepel"      # for labels on plots
)

# Install packages if not already installed, then load them
# This loop goes through each package one by one
for (package in packages_needed) {
  if (!requireNamespace(package, quietly=TRUE)) {
    install.packages(package)
  }
  library(package, character.only=TRUE)
}

# Set up colors and theme for plots
# Make a custom theme for our plots
# I based this on examples from the ggplot2 documentation
theme_pd <- function() {
  theme_minimal(base_size=14) +
    theme(
      plot.title=element_text(face="bold", size=18, hjust=0.5),
      plot.subtitle=element_text(size=12, hjust=0.5),
      legend.position="bottom",
      panel.grid.minor=element_blank()
    )
}

# Colors for different study types
study_colors <- c(
  "Clinical"="#00A676",
  "In vitro"="#FF8C42",
  "In vivo"="#3F88C5"
)

# Colors for different categories
category_colors <- c(
  "drug"="#E4572E",
  "delivery"="#F6AE2D",
  "formulation"="#4E79A7"
)

# Function to save plots - makes it easier
save_plot <- function(plot, filename, width=10, height=6) {
  ggsave(
    file.path(OUTPUT_DIR, filename), 
    plot, 
    width=width, 
    height=height, 
    dpi=300
  )
}

# Define the terms we're looking for
# Sometimes different words mean the same thing, so we normalize them
# For example "nanoparticles" and "nanocarriers" both become "nanoparticle"
normalization_dict <- c(
  "nanoparticles"="nanoparticle",
  "nanocarriers"="nanoparticle",
  "nanocapsules"="nanocapsule",
  "nanogels"="nanoparticle",
  "nanoemulsions"="nanoparticle",
  "nanovesicles"="nanoparticle",
  "liposomes"="liposome",
  "micelles"="micelle",
  "microspheres"="microsphere",
  "hydrogels"="hydrogel",
  "patches"="patch",
  "implants"="implant",
  "encapsulated"="encapsulation",
  "encapsulating"="encapsulation",
  "injection"="injectable",
  "injections"="injectable"
)

# Function to normalize a word
normalize_term <- function(word) {
  # Check if the word is in our dictionary
  if(word %in% names(normalization_dict)) {
    return(normalization_dict[word])
  } else {
    return(word)
  }
}

# List of Parkinson's drugs we want to track
drug_terms <- c(
  "levodopa", "l-dopa", "carbidopa", "benserazide",
  "entacapone", "tolcapone", "opicapone",
  "rasagiline", "selegiline", "safinamide",
  "pramipexole", "ropinirole", "rotigotine",
  "amantadine", "apomorphine"
)

# List of drug delivery methods
delivery_terms <- c(
  "nanoparticle", "liposome", "micelle", "microsphere",
  "hydrogel", "patch", "transdermal", "intranasal", "nasal",
  "inhalation", "oral", "implant", "injectable",
  "buccal", "sublingual", "gel", "spray"
)

# List of formulation-related terms
formulation_terms <- c(
  "controlled", "extended", "sustained", "modified",
  "matrix", "polymeric", "plga", "peg", "chitosan",
  "formulation", "release", "encapsulated", "encapsulation"
)

# Function to assign a category to each word
# This helps us organize our findings
assign_category <- function(word) {
  if(word %in% drug_terms) {
    return("drug")
  } else if(word %in% delivery_terms) {
    return("delivery")
  } else if(word %in% formulation_terms) {
    return("formulation")
  } else {
    return("other")
  }
}

# Function to get articles from PubMed
# This function searches PubMed and downloads article information
# I learned how to use the rentrez package from their documentation
fetch_corpus <- function(search_term, max_results, label) {
  
  # Search PubMed
  # The tryCatch helps handle errors if the search fails
  search_result <- tryCatch(
    rentrez::entrez_search(db="pubmed", term=search_term, retmax=max_results),
    error=function(e) {
      return(NULL)
    }
  )
  
  # If search failed or no results, return empty data frame
  if(is.null(search_result) || length(search_result$ids)==0) {
    return(tibble())
  }
  
  # PubMed only lets us fetch 200 articles at a time
  # So we need to split the IDs into batches
  id_batches <- split(search_result$ids, ceiling(seq_along(search_result$ids)/200))
  
  # Now we fetch each batch
  all_records <- list()  # Create empty list to store results
  
  for(i in 1:length(id_batches)) {
    
    batch <- id_batches[[i]]
    
    # Fetch the XML data for this batch
    xml_data <- tryCatch(
      rentrez::entrez_fetch(db="pubmed", id=batch, rettype="xml", parsed=TRUE),
      error=function(e) {
        return(NULL)
      }
    )
    
    if(is.null(xml_data)) {
      next  # Skip this batch if it failed
    }
    
    # Convert XML to a list format we can work with
    parsed_list <- tryCatch(
      xmlToList(xml_data),
      error=function(e) {
        return(NULL)
      }
    )
    
    if(is.null(parsed_list)) {
      next
    }
    
    # Go through each article in the batch
    for(j in 1:length(parsed_list)) {
      record <- parsed_list[[j]]
      
      # Skip if this isn't a valid record
      if(is.null(record$MedlineCitation)) {
        next
      }
      
      article <- record$MedlineCitation$Article
      
      # Extract the title
      title <- tryCatch(
        as.character(article$ArticleTitle),
        error=function(e) ""
      )
      
      # Extract the abstract
      abstract <- tryCatch({
        abstract_text <- article$Abstract$AbstractText
        if(is.null(abstract_text)) {
          ""
        } else {
          # Sometimes abstract is in multiple parts, so combine them
          paste(unlist(abstract_text), collapse=" ")
        }
      }, error=function(e) "")
      
      # Extract the year
      year <- tryCatch(
        as.integer(article$Journal$JournalIssue$PubDate$Year),
        error=function(e) NA_integer_
      )
      
      # Create a data frame for this article
      article_data <- tibble(
        pmid=as.character(record$MedlineCitation$PMID[[1]]),
        year=year,
        title=title,
        abstract=abstract,
        study_type=label
      )
      
      # Add to our list
      all_records[[length(all_records) + 1]] <- article_data
    }
  }
  
  # Combine all records into one data frame
  if(length(all_records)==0) {
    return(tibble())
  }
  
  combined_records <- bind_rows(all_records)
  
  # Remove duplicates (same article might appear in multiple searches)
  combined_records <- combined_records %>%
    distinct(pmid, .keep_all=TRUE) %>%
    filter(!is.na(abstract), abstract!="")
  
  return(combined_records)
}

# Build the search queries for PubMed
# Create the year range filter
year_range <- paste0(YEAR_MIN, ":", YEAR_MAX, "[DP]")

# Base query that all searches will use
# This searches for Parkinson's AND drug delivery related terms
base_query <- paste(
  '(Parkinson Disease[MeSH] OR parkinson*[tiab])',
  'AND ("drug delivery"[tiab] OR nanoparticle*[tiab] OR liposome*[tiab] OR formulation[tiab])'
)

# Query for in vitro studies
query_invitro <- paste(
  base_query,
  'AND ("in vitro"[tiab] OR cell culture[tiab])',
  year_range
)

# Query for in vivo studies (animal studies)
query_invivo <- paste(
  base_query,
  'AND ("in vivo"[tiab] OR rat[tiab] OR mouse[tiab])',
  year_range
)

# Query for clinical studies (human studies)
query_clinical <- paste(
  base_query,
  'AND ("clinical trial"[pt] OR phase[tiab])',
  year_range
)

# Get the data from PubMed
if(RUN_SEARCH) {
  
  # Search for each type of study
  df_invitro <- fetch_corpus(query_invitro, MAX_ARTICLES, "In vitro")
  df_invivo <- fetch_corpus(query_invivo, MAX_ARTICLES, "In vivo")
  df_clinical <- fetch_corpus(query_clinical, MAX_ARTICLES, "Clinical")
  
  # Combine all the results
  article_df <- bind_rows(df_invitro, df_invivo, df_clinical)
  
  # Remove rows with missing years
  article_df <- article_df %>% filter(!is.na(year))
  
  # Add a document ID to each article
  article_df <- article_df %>% mutate(doc_id=paste0("doc", row_number()))
  
  # Save the data so we don't have to search again
  saveRDS(article_df, file.path(OUTPUT_DIR, "article_df_clean.rds"))
  
} else {
  # Load previously saved data
  article_df <- readRDS(file.path(OUTPUT_DIR, "article_df_clean.rds"))
}

# Process the text - break into words and categorize
# Load standard stop words (common words like "the", "and", etc.)
data("stop_words")

# Add custom stop words specific to our analysis
custom_stop_words <- tibble(word=c(
  "study", "result", "results", "conclusion", "analysis",
  "method", "patients", "clinical", "parkinson",
  "parkinsons", "disease"
))

# Combine all stop words
all_stop_words <- bind_rows(stop_words, custom_stop_words)
all_stop_words <- distinct(all_stop_words)  # Remove duplicates

# Now process the text
# Combine title and abstract for each article
article_df <- article_df %>%
  mutate(text=paste(title, abstract))

# Break text into individual words (tokenization)
df_tokens <- article_df %>%
  select(doc_id, study_type, year, text) %>%
  unnest_tokens(word, text)  # This breaks text into words

# Remove stop words
df_tokens <- df_tokens %>%
  anti_join(all_stop_words, by="word")

# Normalize terms and assign categories
df_tokens <- df_tokens %>%
  mutate(year=as.integer(year))

# Normalize each word one at a time
df_tokens$term <- sapply(df_tokens$word, normalize_term)

# Add category for each word
df_tokens$category <- sapply(df_tokens$word, assign_category)

# Normalize each word one at a time
df_tokens$term <- sapply(df_tokens$word, normalize_term)

# Add category for each word
df_tokens$category <- sapply(df_tokens$word, assign_category)

# Make study_type a factor with specific order
df_tokens$study_type <- factor(
  df_tokens$study_type,
  levels=c("Clinical", "In vitro", "In vivo")
)

article_df$study_type <- factor(
  article_df$study_type,
  levels=c("Clinical", "In vitro", "In vivo")
)

# Save the processed tokens
saveRDS(df_tokens, file.path(OUTPUT_DIR, "df_tokens.rds"))

# Helper functions for making visualizations
# Function to get top terms for each category and study type
get_top_terms <- function(tokens, category_name, n=9) {
  # Filter for specific category
  filtered <- tokens %>% filter(category==category_name)
  
  # Count how many times each term appears in each study type
  counted <- filtered %>% count(study_type, term, name="n")
  
  # For each study type, get the top N terms
  top_terms <- counted %>%
    group_by(study_type) %>%
    slice_max(n, n=n, with_ties=FALSE) %>%
    ungroup()
  
  return(top_terms)
}

# Function to plot how terms evolve over time
plot_term_evolution <- function(tokens, category_name, n=8, filename) {
  
  # First, find the overall top N terms for this category
  top_terms <- tokens %>%
    filter(category==category_name) %>%
    count(term, name="total") %>%
    slice_max(total, n=n, with_ties=FALSE)
  
  top_term_list <- top_terms$term
  
  # Now count mentions by year for these top terms
  evolution_data <- tokens %>%
    filter(
      category==category_name,
      term %in% top_term_list,
      !is.na(year)
    ) %>%
    count(year, study_type, term, name="n")
  
  # Create the plot
  p <- ggplot(evolution_data, aes(x=year, y=n, colour=study_type)) +
    geom_line(linewidth=0.7) +
    facet_wrap(~term, scales="free_y") +  # Separate plot for each term
    scale_colour_manual(values=study_colors) +
    labs(
      title=paste0("Evolution of Top ", n, " ", 
                   tools::toTitleCase(category_name), " Terms Over Time"),
      x="Year",
      y="Mentions per Year"
    ) +
    theme_pd() +
    theme(strip.text=element_text(face="bold"))
  
  # Save the plot
  save_plot(p, filename, width=12, height=8)
}

# Figure 1: Publications over time by study type
# Count publications per year for each study type
pub_counts <- article_df %>% count(year, study_type)

# Create line plot
p_publications <- ggplot(pub_counts, aes(x=year, y=n, color=study_type)) +
  geom_line(linewidth=1.2) +
  geom_point(size=2) +
  scale_color_manual(values=study_colors) +
  labs(
    title="Publications Over Time by Study Type",
    x="Year",
    y="Number of Publications"
  ) +
  theme_pd()

# Save the plot
save_plot(p_publications, "01_publications_by_study_type.png", 11, 6)

# Figure 2: Top 9 drug terms by study type
# Get top drug terms
top_drugs <- get_top_terms(df_tokens, "drug", n=9)

# Reorder terms by frequency for better visualization
top_drugs <- top_drugs %>% mutate(term=fct_reorder(term, n))

# Create bar plot
p_drugs <- ggplot(top_drugs, aes(x=term, y=n)) +
  geom_col(fill=category_colors["drug"]) +
  geom_text(aes(label=n), hjust=-0.1, size=3) +  # Add count labels
  coord_flip() +  # Horizontal bars
  facet_wrap(~study_type, ncol=1, scales="free_y") +
  labs(
    title="Top 9 Drug Terms by Study Type",
    x=NULL,
    y="Frequency"
  ) +
  theme_pd()

save_plot(p_drugs, "02_top_terms_drug_by_studytype.png", 9, 9)

# Figure 3: Top 9 delivery terms by study type
top_delivery <- get_top_terms(df_tokens, "delivery", n=9)
top_delivery <- top_delivery %>% mutate(term=fct_reorder(term, n))

p_delivery <- ggplot(top_delivery, aes(x=term, y=n)) +
  geom_col(fill=category_colors["delivery"]) +
  geom_text(aes(label=n), hjust=-0.1, size=3) +
  coord_flip() +
  facet_wrap(~study_type, ncol=1, scales="free_y") +
  labs(
    title="Top 9 Delivery Terms by Study Type",
    x=NULL,
    y="Frequency"
  ) +
  theme_pd()

save_plot(p_delivery, "03_top_terms_delivery_by_studytype.png", 9, 9)

# Figure 4: Top 9 formulation terms by study type
top_formulation <- get_top_terms(df_tokens, "formulation", n=9)
top_formulation <- top_formulation %>% mutate(term=fct_reorder(term, n))

p_formulation <- ggplot(top_formulation, aes(x=term, y=n)) +
  geom_col(fill=category_colors["formulation"]) +
  geom_text(aes(label=n), hjust=-0.1, size=3) +
  coord_flip() +
  facet_wrap(~study_type, ncol=1, scales="free_y") +
  labs(
    title="Top 9 Formulation Terms by Study Type",
    x=NULL,
    y="Frequency"
  ) +
  theme_pd()

save_plot(p_formulation, "04_top_terms_formulation_by_studytype.png", 9, 9)

# Figures 5-7: Term evolution over time
plot_term_evolution(df_tokens, "drug", n=8, filename="05_term_evolution_drug.png")
plot_term_evolution(df_tokens, "delivery", n=8, filename="06_term_evolution_delivery.png")
plot_term_evolution(df_tokens, "formulation", n=8, filename="07_term_evolution_formulation.png")

# Figure 8: Translational gap analysis
# Count delivery terms by study type
trans_counts <- df_tokens %>%
  filter(category=="delivery") %>%
  count(study_type, term, name="n")

# Reshape data to have columns for each study type
trans_wide <- trans_counts %>%
  pivot_wider(names_from=study_type, values_from=n, values_fill=0)

# Calculate total mentions and get top 20
trans_wide <- trans_wide %>%
  mutate(total=Clinical + `In vitro` + `In vivo`) %>%
  slice_max(total, n=20, with_ties=FALSE)

# Create scatter plot
p_translational <- ggplot(
  trans_wide,
  aes(x=`In vitro`, y=Clinical, size=`In vivo`, label=term)
) +
  geom_point(alpha=0.7, colour="#FF8C42") +
  geom_text_repel(size=3) +  # Labels that don't overlap
  scale_size(range=c(2, 10)) +
  labs(
    title="Translational Gap in Parkinson's Drug Delivery",
    subtitle="In vitro vs Clinical; size = In vivo",
    x="In Vitro Mentions",
    y="Clinical Mentions"
  ) +
  theme_pd()

save_plot(p_translational, "08_translational_gap.png", 10, 7)

# Figure 9: Drug x Delivery co-occurrence heatmap (global)
# Find which drugs appear in each document
doc_drugs <- df_tokens %>%
  filter(category=="drug") %>%
  distinct(doc_id, term) %>%
  rename(drug=term)

# Find which delivery systems appear in each document
doc_delivery <- df_tokens %>%
  filter(category=="delivery") %>%
  distinct(doc_id, term) %>%
  rename(delivery=term)

# Join to find which drugs and delivery systems appear together
pairs <- inner_join(doc_drugs, doc_delivery, by="doc_id")

# Count co-occurrences
cooccurrence <- pairs %>% count(drug, delivery, name="n")

# Get top 10 drugs
top_drugs_list <- cooccurrence %>%
  count(drug, wt=n, name="total") %>%
  slice_max(total, n=10) %>%
  pull(drug)

# Get top 10 delivery systems
top_delivery_list <- cooccurrence %>%
  count(delivery, wt=n, name="total") %>%
  slice_max(total, n=10) %>%
  pull(delivery)

# Filter to only top terms
cooccurrence_filtered <- cooccurrence %>%
  filter(drug %in% top_drugs_list, delivery %in% top_delivery_list)

# Create heatmap
p_heatmap <- ggplot(cooccurrence_filtered, aes(x=delivery, y=drug, fill=n)) +
  geom_tile(colour="white", linewidth=0.6) +
  scale_fill_gradient(
    low="#FFF5EB",
    high="#E4572E",
    name="Co-occurrence"
  ) +
  labs(
    title="Drug × Delivery System Co-occurrence (Same Article)",
    x="Delivery System",
    y="Drug"
  ) +
  theme_pd() +
  theme(axis.text.x=element_text(angle=45, hjust=1))

save_plot(p_heatmap, "09_heatmap_drug_delivery_global.png", 10, 7)

# Figure 10: Drug x Delivery co-occurrence by study type
# Similar to above but keeping study type
doc_drugs_by_study <- df_tokens %>%
  filter(category=="drug") %>%
  distinct(doc_id, study_type, term) %>%
  rename(drug=term)

doc_delivery_by_study <- df_tokens %>%
  filter(category=="delivery") %>%
  distinct(doc_id, study_type, term) %>%
  rename(delivery=term)

pairs_by_study <- inner_join(
  doc_drugs_by_study,
  doc_delivery_by_study,
  by=c("doc_id", "study_type")
)

cooccurrence_by_study <- pairs_by_study %>%
  count(study_type, drug, delivery, name="n")

# Get top 8 drugs across all study types
top_drugs_study <- cooccurrence_by_study %>%
  group_by(drug) %>%
  summarise(total=sum(n), .groups="drop") %>%
  slice_max(total, n=8) %>%
  pull(drug)

# Get top 8 delivery systems
top_delivery_study <- cooccurrence_by_study %>%
  group_by(delivery) %>%
  summarise(total=sum(n), .groups="drop") %>%
  slice_max(total, n=8) %>%
  pull(delivery)

# Filter to top terms
cooccurrence_study_filtered <- cooccurrence_by_study %>%
  filter(drug %in% top_drugs_study, delivery %in% top_delivery_study)

# Create heatmap with facets
p_heatmap_study <- ggplot(
  cooccurrence_study_filtered,
  aes(x=delivery, y=drug, fill=n)
) +
  geom_tile(colour="white", linewidth=0.5) +
  scale_fill_gradient(
    low="#FFF5EB",
    high="#E4572E",
    name="Co-occurrence"
  ) +
  facet_wrap(~study_type, nrow=1) +
  labs(
    title="Drug × Delivery System Co-occurrence by Study Type",
    x="Delivery System",
    y="Drug"
  ) +
  theme_pd() +
  theme(
    axis.text.x=element_text(angle=45, hjust=1),
    panel.border=element_rect(colour="grey30", fill=NA, linewidth=1)
  )

save_plot(p_heatmap_study, "10_heatmap_drug_delivery_by_studytype.png", 13, 5)

# Figure 11: Topic modeling with LDA
# Prepare data for topic modeling
# Count how often each term appears in each document
dtm_data <- df_tokens %>%
  filter(term!="") %>%
  count(doc_id, term) %>%
  rename(freq=n)

# Convert to document-term matrix format
dtm <- cast_dtm(dtm_data, document=doc_id, term=term, value=freq)

# Remove documents with no terms (empty rows)
dtm <- dtm[rowSums(as.matrix(dtm))>0, ]

# Run LDA topic model with 4 topics
# I learned about LDA from online tutorials
set.seed(123)  # Set seed for reproducibility
lda_model <- LDA(dtm, k=4, control=list(seed=123))

# Extract top terms for each topic
lda_terms <- tidy(lda_model, matrix="beta")

# Get top 12 terms per topic
lda_top_terms <- lda_terms %>%
  group_by(topic) %>%
  slice_max(beta, n=12) %>%
  arrange(topic, desc(beta))

# Save results to CSV
write.csv(
  lda_top_terms,
  file.path(OUTPUT_DIR, "11_LDA_top_terms_by_topic.csv"),
  row.names=FALSE
)

# Create visualization
p_lda <- lda_top_terms %>%
  mutate(term=reorder_within(term, beta, topic)) %>%
  ggplot(aes(x=term, y=beta, fill=factor(topic))) +
  geom_col(show.legend=FALSE) +
  facet_wrap(~topic, scales="free") +
  scale_x_reordered() +
  coord_flip() +
  labs(
    title="Top Terms for Each LDA Topic",
    x="Term",
    y="Weight (β)"
  ) +
  theme_minimal(base_size=14)

save_plot(p_lda, "11_LDA_top_terms_plot.png", 12, 8)

# Figure 12: ABC Hypothesis Generation
# This section finds untested drug delivery combinations
# We use the "ABC method" - if Drug A works with Delivery B (papers exist),
# and Delivery B works with Formulation C (papers exist),
# then maybe Drug A + Delivery B + Formulation C is worth testing!
# But we only want combinations where all three together haven't been tested yet

# Load extra package we need for making tables
library(gt)

# Step 1: Find which drugs appear with which delivery methods

cat("\nFinding Drug + Delivery Combinations...\n")

# Find documents that mention both a drug AND a delivery method
drug_delivery_docs <- df_tokens %>%
  filter(category %in% c("drug", "delivery")) %>%
  distinct(doc_id, category, term) %>%
  # Reshape so we can see both categories at once
  pivot_wider(names_from = category, values_from = term, values_fn = list) %>%
  # Expand the lists
  unnest(cols = drug) %>%
  unnest(cols = delivery) %>%
  distinct(doc_id, drug, delivery)

# Count how many papers mention each drug-delivery pair
drug_delivery_support <- drug_delivery_docs %>%
  count(drug, delivery, name = "drug_delivery_papers")

cat("Found", nrow(drug_delivery_support), "drug-delivery combinations\n")

# Step 2: Find which delivery methods appear with which formulations

cat("\nFinding Delivery + Formulation Combinations...\n")

# Find documents that mention both a delivery method AND a formulation
delivery_formulation_docs <- df_tokens %>%
  filter(category %in% c("delivery", "formulation")) %>%
  distinct(doc_id, category, term) %>%
  pivot_wider(names_from = category, values_from = term, values_fn = list) %>%
  unnest(cols = delivery) %>%
  unnest(cols = formulation) %>%
  distinct(doc_id, delivery, formulation)

# Count how many papers mention each delivery-formulation pair
delivery_formulation_support <- delivery_formulation_docs %>%
  count(delivery, formulation, name = "delivery_formulation_papers")

cat("Found", nrow(delivery_formulation_support), "delivery-formulation combinations\n")

# Step 3: Connect the paths (Drug to Delivery to Formulation)

cat("\nBuilding Complete Paths (Drug -> Delivery -> Formulation)...\n")

# Join drug-delivery pairs with delivery-formulation pairs
# The "delivery" term is what connects them
all_paths <- drug_delivery_docs %>%
  distinct(drug, delivery) %>%
  inner_join(
    delivery_formulation_docs %>% distinct(delivery, formulation), 
    by = "delivery"
  ) %>%
  distinct(drug, delivery, formulation)

cat("Found", nrow(all_paths), "possible complete paths\n")

# Step 4: Find which combinations are ALREADY tested (to exclude them)

cat("\nFinding Already Tested Combinations...\n")

# Find documents where drug AND formulation appear together
# (These are combinations that have been directly tested)
known_combinations <- df_tokens %>%
  filter(category %in% c("drug", "formulation")) %>%
  distinct(doc_id, category, term) %>%
  pivot_wider(names_from = category, values_from = term, values_fn = list) %>%
  unnest(cols = drug) %>%
  unnest(cols = formulation) %>%
  distinct(drug, formulation)

cat("Found", nrow(known_combinations), "already tested drug-formulation combinations\n")

# Step 5: Keep only UNTESTED combinations

cat("\nFinding Novel (Untested) Hypotheses...\n")

# Remove paths where the drug-formulation combination already exists
untested_hypotheses <- all_paths %>%
  # Remove combinations that are already known
  anti_join(known_combinations, by = c("drug", "formulation")) %>%
  # Add the supporting paper counts
  left_join(drug_delivery_support, by = c("drug", "delivery")) %>%
  left_join(delivery_formulation_support, by = c("delivery", "formulation")) %>%
  # Replace missing values with 0
  mutate(
    drug_delivery_papers = coalesce(drug_delivery_papers, 0L),
    delivery_formulation_papers = coalesce(delivery_formulation_papers, 0L),
    total_support = drug_delivery_papers + delivery_formulation_papers
  ) %>%
  # Sort by total support (best hypotheses first)
  arrange(desc(total_support), desc(drug_delivery_papers), desc(delivery_formulation_papers))

cat("Found", nrow(untested_hypotheses), "UNTESTED hypotheses!\n")

# Save all hypotheses to a file
write.csv(
  untested_hypotheses,
  file.path(OUTPUT_DIR, "12_all_untested_hypotheses.csv"),
  row.names = FALSE
)

cat("Saved complete list to: 12_all_untested_hypotheses.csv\n")

# Step 6: Select the Top 10 most promising hypotheses

cat("\nSelecting Top 10 Most Promising Hypotheses...\n")

# Get top 10 with variety (avoid showing same drug multiple times)
# First, for each drug, find its best hypothesis
best_per_drug <- untested_hypotheses %>%
  filter(total_support >= 5) %>%
  group_by(drug) %>%
  slice_max(total_support, n = 1, with_ties = FALSE) %>%
  ungroup()

# Now take the top 10 from these best ones, sorted by total support
top10 <- best_per_drug %>%
  arrange(desc(total_support), desc(drug_delivery_papers), desc(delivery_formulation_papers)) %>%
  slice(1:10) %>%
  # Make the names look nice (capitalize first letter)
  mutate(
    drug_clean = str_to_title(drug),
    delivery_clean = str_to_title(delivery),
    formulation_clean = str_to_title(formulation)
  )

cat("Selected top 10 diverse hypotheses\n")

# Print the best one
cat("\nBest Hypothesis:\n")
cat("Drug:", top10$drug_clean[1], "\n")
cat("Delivery:", top10$delivery_clean[1], "\n")
cat("Formulation:", top10$formulation_clean[1], "\n")
cat("Total Support:", top10$total_support[1], "papers\n")
cat("(", top10$drug_delivery_papers[1], "drug-delivery papers +", 
    top10$delivery_formulation_papers[1], "delivery-formulation papers)\n")

# Figure 13: Table of Top 10 Hypotheses

cat("\nCreating Table...\n")

# Prepare data for the table
table_data <- top10 %>%
  transmute(
    Drug = drug_clean,
    Delivery = delivery_clean,
    Formulation = formulation_clean,
    `Drug Delivery Papers` = drug_delivery_papers,
    `Delivery Formulation Papers` = delivery_formulation_papers,
    `Total Support` = total_support
  )

# Create the table using gt
# I learned about gt from online tutorials
top10_table <- gt(table_data) %>%
  tab_header(
    title = md("**Top 10 Untested Drug Delivery Hypotheses**"),
    subtitle = "Hypotheses inferred from indirect literature evidence (ABC method)"
  ) %>%
  cols_label(
    `Drug Delivery Papers` = md("**Drug-Delivery**<br>papers"),
    `Delivery Formulation Papers` = md("**Delivery-Formulation**<br>papers"),
    `Total Support` = md("**Total**<br>support")
  ) %>%
  cols_align(align = "left", columns = c(Drug, Delivery, Formulation)) %>%
  cols_align(align = "center", columns = c(`Drug Delivery Papers`, `Delivery Formulation Papers`, `Total Support`)) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  ) %>%
  tab_style(
    style = cell_fill(color = "#F0F0F0"),
    locations = cells_body(rows = seq(1, 10, 2))  # Zebra striping for every other row
  ) %>%
  tab_style(
    style = cell_text(size = px(14)),
    locations = cells_body()
  ) %>%
  tab_options(
    table.font.size = px(14),
    heading.title.font.size = px(18),
    heading.subtitle.font.size = px(14)
  ) %>%
  tab_footnote(
    footnote = "Drug-Delivery papers: papers mentioning both drug and delivery method",
    locations = cells_column_labels(columns = `Drug Delivery Papers`)
  ) %>%
  tab_footnote(
    footnote = "Delivery-Formulation papers: papers mentioning both delivery method and formulation",
    locations = cells_column_labels(columns = `Delivery Formulation Papers`)
  )

# Save table to HTML and PDF
gtsave(
  top10_table,
  file.path(OUTPUT_DIR, "13_Table_Top10_Hypotheses.html")
)

gtsave(
  top10_table,
  file.path(OUTPUT_DIR, "13_Table_Top10_Hypotheses.png")
)

cat("Table saved to: 13_Table_Top10_Hypotheses.html and .pdf\n")

# Figure 14: Bar chart showing evidence for each hypothesis

cat("\nCreating Bar Chart...\n")

# Prepare data for plotting
# We need to make it so each hypothesis has two bars (one for each type of evidence)
plot_data <- top10 %>%
  mutate(
    rank = row_number(),
    # Create a single line label with all three components
    hypothesis_label = paste0(rank, ". ", drug_clean, " + ", delivery_clean, " + ", formulation_clean)
  ) %>%
  select(rank, hypothesis_label, drug_delivery_papers, delivery_formulation_papers) %>%
  # Reshape the data so we can have two bars per hypothesis
  pivot_longer(
    cols = c(drug_delivery_papers, delivery_formulation_papers),
    names_to = "evidence_type",
    values_to = "papers"
  ) %>%
  mutate(
    # Make the legend labels nicer
    evidence_type = case_when(
      evidence_type == "drug_delivery_papers" ~ "Drug-Delivery",
      evidence_type == "delivery_formulation_papers" ~ "Delivery-Formulation"
    ),
    # Reverse the order so rank 1 is at the top
    hypothesis_label = factor(hypothesis_label, levels = rev(unique(hypothesis_label)))
  )

# Create the bar chart
p_hypotheses <- ggplot(
  plot_data, 
  aes(x = papers, y = hypothesis_label, fill = evidence_type)
) +
  # Make the bars
  geom_col(position = position_dodge(width = 0.8), width = 0.75, alpha = 0.9) +
  # Add numbers on the bars
  geom_text(
    aes(label = papers),
    position = position_dodge(width = 0.8),
    hjust = -0.2,
    size = 4,
    fontface = "bold"
  ) +
  # Set the colors
  scale_fill_manual(
    values = c(
      "Drug-Delivery" = "#3F88C5",
      "Delivery-Formulation" = "#E4572E"
    ),
    name = NULL
  ) +
  # Add labels
  labs(
    title = "Top 10 Novel Drug Delivery Hypotheses for Parkinson's Disease",
    subtitle = "Ranked by strength of indirect evidence from the literature",
    x = "Number of Supporting Papers",
    y = NULL
  ) +
  # Extend x-axis so numbers don't get cut off
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  # Apply our custom theme
  theme_pd() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "#555555"),
    legend.position = "bottom",
    legend.text = element_text(size = 11),
    axis.title.x = element_text(face = "bold", size = 12, margin = margin(t = 10)),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

# Save the plot
save_plot(p_hypotheses, "14_Top10_Hypotheses.png", width = 14, height = 8)

cat("Bar chart saved to: 14_Top10_Hypotheses.png\n")

cat("\nAll ABC hypothesis visualizations complete!\n")
cat("These represent new research directions that haven't been tried yet!\n")


# PRISMA-STYLE FLOWCHART (Figure 15)
# Adds a clean PRISMA-style flowchart using available statistics.

# Install/load flowchart export packages (kept local so nothing else in your script changes)
for (pkg in c("DiagrammeR", "DiagrammeRsvg", "rsvg")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
library(DiagrammeR)
library(DiagrammeRsvg)

# Try to use saved search statistics if available (preferred)
stats_path <- file.path(OUTPUT_DIR, "search_statistics.rds")

if (file.exists(stats_path)) {
  search_stats <- readRDS(stats_path)

  total_found <- as.integer(search_stats$invitro$found + search_stats$invivo$found + search_stats$clinical$found)

  # "Retrieved" here corresponds to how many valid records were fetched before deduplication in the tracked script
  retrieved_total <- as.integer(search_stats$total_before_dedup)

  excluded_missing_abstract <- as.integer(search_stats$invitro$excluded_no_abstract +
                                            search_stats$invivo$excluded_no_abstract +
                                            search_stats$clinical$excluded_no_abstract)

  # This script does not separately track "downloaded in MEDLINE"; keep it consistent with the tracked script logic:
  downloaded_total <- as.integer(retrieved_total - excluded_missing_abstract)

  duplicates_removed <- as.integer(search_stats$duplicates_removed)
  unique_records <- as.integer(search_stats$final_corpus)

  clinical_count <- as.integer(search_stats$final_by_type$clinical)
  invitro_count  <- as.integer(search_stats$final_by_type$invitro)
  invivo_count   <- as.integer(search_stats$final_by_type$invivo)

  final_corpus <- as.integer(search_stats$final_corpus)

} else {
  # Fallback: compute what we can from the objects in this script
  # Note: this script does not track "records identified" or "excluded: missing abstract" during fetch,
  # so those fields are set to NA and not displayed as exclusions.
  retrieved_total <- as.integer(nrow(article_df))

  if ("pmid" %in% names(article_df)) {
    unique_df <- article_df %>% dplyr::distinct(pmid, .keep_all = TRUE)
  } else {
    unique_df <- article_df
  }

  unique_records <- as.integer(nrow(unique_df))

  clinical_count <- as.integer(sum(unique_df$study_type == "Clinical", na.rm = TRUE))
  invitro_count  <- as.integer(sum(unique_df$study_type == "In vitro", na.rm = TRUE))
  invivo_count   <- as.integer(sum(unique_df$study_type == "In vivo", na.rm = TRUE))

  duplicates_removed <- as.integer(retrieved_total - unique_records)

  total_found <- NA_integer_
  excluded_missing_abstract <- NA_integer_
  downloaded_total <- NA_integer_
  final_corpus <- unique_records
}

# Build a clean, aligned PRISMA-style flowchart
flowchart <- DiagrammeR::grViz(sprintf("
digraph PRISMA {

  graph [
    layout = dot,
    rankdir = TB,
    bgcolor = \"white\",
    pad = 0.30,
    nodesep = 0.65,
    ranksep = 0.55,
    splines = curved
  ]

  node [
    shape = box,
    style = \"rounded,filled\",
    fontname = \"Helvetica\",
    fontsize = 11,
    margin = 0.20,
    penwidth = 1.2,
    color = \"#D1D5DB\",
    fillcolor = \"#F3F4F6\"
  ]

  edge [
    color = \"#4B5563\",
    penwidth = 2.4,
    arrowsize = 0.85
  ]

  // Main column (fixed widths help alignment)
  n1 [label = \"Records identified\\nNLM search\\n%s\", width=3.6, height=0.95]
  n2 [label = \"Records retrieved from PubMed\\n(n = %d)\", width=3.6, height=0.85]
  n3 [label = \"Unique records after deduplication\\n(n = %d)\", width=3.6, height=0.85]
  n4 [label = \"Manual validation\\nClassified by study type\\n(n = %d)\", width=3.6, height=0.95]

  // Right-side notes (subtle callouts)
  ex2 [label = \"Removed\\nDuplicates (PMIDs)\\n(n = %d)\",
       fillcolor=\"#FFF3E6\", color=\"#FFD7AE\", width=2.55, height=0.75, fontsize=10]

  // Study type row
  node [fillcolor=\"#F3E8FF\", color=\"#E2D2FF\"]
  c1 [label = \"Clinical\\n(n = %d)\", width=2.15, height=0.75]
  c2 [label = \"In vitro\\n(n = %d)\", width=2.15, height=0.75]
  c3 [label = \"In vivo\\n(n = %d)\", width=2.15, height=0.75]

  // Final
  node [fillcolor=\"#6D28D9\", color=\"#6D28D9\", fontcolor=white, fontsize=12]
  final [label = \"FINAL CORPUS\\n(n = %d)\", width=2.75, height=0.90]

  // Flow
  n1 -> n2
  n2 -> n3
  n3 -> n4

  // Side callout
  n2 -> ex2 [style=dashed, color=\"#F59E0B\", penwidth=2, constraint=false]

  // Split
  n4 -> c1 [color=\"#7C3AED\"]
  n4 -> c2 [color=\"#7C3AED\"]
  n4 -> c3 [color=\"#7C3AED\"]

  // Keep study boxes on same row
  { rank = same; c1; c2; c3 }

  // Merge to final
  c1 -> final [color=\"#7C3AED\"]
  c2 -> final [color=\"#7C3AED\"]
  c3 -> final [color=\"#7C3AED\"]
}
",
if (is.na(total_found)) "(n = Not tracked in this script)" else sprintf("(n = %d)", total_found),
retrieved_total, unique_records, unique_records, duplicates_removed,
clinical_count, invitro_count, invivo_count, final_corpus))

flowchart

# Export as high-resolution PNG and PDF
flowchart %>%
  DiagrammeRsvg::export_svg() %>%
  charToRaw() %>%
  rsvg::rsvg_png(file.path(OUTPUT_DIR, "15_PRISMA_Clean.png"), width = 2600, height = 2000)

flowchart %>%
  DiagrammeRsvg::export_svg() %>%
  charToRaw() %>%
  rsvg::rsvg_pdf(file.path(OUTPUT_DIR, "15_PRISMA_Clean.pdf"))
