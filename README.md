# Parkinson’s Drug Delivery Text Mining (R)

This repository contains the full analysis pipeline used in the accompanying article to map the **drug delivery landscape in Parkinson’s disease** using PubMed abstracts and text mining.

The script is designed to be **single-file and reproducible**: it searches PubMed (optional), builds a clean corpus, tokenises text, extracts drug/delivery/formulation signals, and exports publication-ready figures and tables.

## What the code does

1. **Fetches PubMed records** for Parkinson’s + drug delivery terms  
   - Splits the search into **in vitro**, **in vivo**, and **clinical** queries  
   - Tracks exclusions (missing abstracts, missing year, fetch/parsing errors) for PRISMA-style reporting  
2. **Cleans and tokenises** titles + abstracts  
   - Removes stop words + a small custom stop word list  
   - Normalises common variants (for example `nanoparticles` → `nanoparticle`)  
   - Assigns each token into one of: **drug**, **delivery**, **formulation**, or **other**
3. **Exports figures** describing the landscape  
   - Publications per year (by study type)  
   - Top terms per category (drug/delivery/formulation)  
   - Term evolution plots (top terms over time)  
   - Translational gap scatter (in vitro vs clinical; size = in vivo)  
   - Drug × delivery co-occurrence heatmaps (overall and by study type)  
4. **Topic modelling (LDA)**  
   - Fits an LDA model with 4 topics  
   - Exports top terms per topic and a faceted bar plot  
5. **ABC hypothesis generation**  
   - Generates candidate **Drug → Delivery → Formulation** paths  
   - Removes combinations that are already observed as **Drug–Formulation** within the same articles  
   - Exports a full hypothesis list plus a top-10 summary table and bar chart

## Files

- `parkinsons_drug_delivery_textmining.R` — main analysis script (run this)
- Outputs are written to `config$output_dir` (you choose this path)

## How to run

1. Install R (4.2+ recommended) and RStudio (optional)
2. Download or clone this repository
3. Open `parkinsons_drug_delivery_textmining.R`
4. Edit the configuration block at the top:

```r
config <- list(
  output_dir   = "Insert Path",
  run_search   = TRUE,
  max_articles = 800,
  year_min     = 1995,
  year_max     = 2025,
  seed         = 123
)
```

5. Run the script (it will install missing packages the first time)

### If you want to skip PubMed searching

If you already ran it once (and saved the cleaned corpus), set:

```r
config$run_search <- FALSE
```

The script will load `article_df_clean.rds` and `search_statistics.rds` from your output directory.

## Outputs you should expect

Figures (PNG):

- `01_publications_by_study_type.png`
- `02_top_terms_drug_by_studytype.png`
- `03_top_terms_delivery_by_studytype.png`
- `04_top_terms_formulation_by_studytype.png`
- `05_term_evolution_drug.png`
- `06_term_evolution_delivery.png`
- `07_term_evolution_formulation.png`
- `08_translational_gap.png`
- `09_heatmap_drug_delivery_global.png`
- `10_heatmap_drug_delivery_by_studytype.png`
- `11_LDA_top_terms_plot.png`
- `14_Top10_Hypotheses.png`
- 'Prisma_Flowchart'
Tables / data:

- `search_statistics.rds` (exclusion counts and totals)
- `article_df_clean.rds` (clean corpus)
- `df_tokens.rds` (token-level dataset)
- `11_LDA_top_terms_by_topic.csv`
- `12_all_untested_hypotheses.csv`
- `13_Table_Top10_Hypotheses.html`
- `13_Table_Top10_Hypotheses.png`

## Customising dictionaries

Three character vectors control classification:

- `drug_terms`
- `delivery_terms`
- `formulation_terms`

To extend the analysis (for example adding new delivery technologies), edit the relevant vector and re-run.

## Notes and limitations

- PubMed metadata can be inconsistent (some records have no usable year or abstract). The script excludes these and reports counts.
- Terms are dictionary-based and intentionally simple. This is transparent and reproducible, but it may miss synonyms not in the dictionary.

## Citation

If you use this code, please cite the associated article and link back to this repository.
