## Notebook templates

This folder includes `Rmd` files for step `02_clustering`.

```bash
cd /path/to/OpenScPCA-analysis
cd analyses/cell-type-wilms-tumor-14

# required for rendering Rmd files.
# This path exists on OpenScPCA Lightsail for Research instances, but you may need to specify a different path if running on a different platform, or run `rmarkdown::render` inside your rstudio session.
export RSTUDIO_PANDOC="/usr/lib/rstudio/resources/app/bin/quarto/bin/tools"

metadata="../../data/current/SCPCP000014/single_cell_metadata.tsv"
for library_id in `sed '1d' ${metadata} | cut -f3`; do
    level="celltype"
    Rscript -e "rmarkdown::render('./notebook_templates/02_clustering/02_clustering_explore.Rmd', \
        clean = TRUE, \
        output_file = '02_clustering_explore_${library_id}_${level}.html', \
            params = list(library_id = '${library_id}', \
                          level = '${level}') \
        )"
    level="compartment"
    Rscript -e "rmarkdown::render('./notebook_templates/02_clustering/02_clustering_explore.Rmd', \
        clean = TRUE, \
        output_file = '02_clustering_explore_${library_id}_${level}.html', \
            params = list(library_id = '${library_id}', \
                          level = '${level}') \
        )"
done
```