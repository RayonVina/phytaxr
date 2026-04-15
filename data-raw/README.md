# data-raw/

This directory contains the scripts that **generate** the internal data
objects shipped with the package. They are kept here for reproducibility
but **do not need to be re-run** unless the source data changes.

---

## Files

### `vernacular.R`

Builds the internal vernacular-to-taxon dictionary (`R/sysdata.rda`).

The dictionary maps common phytoplankton vernacular names to accepted
scientific names (e.g. `"diatom"` → `"Bacillariophyceae"`). It is used
automatically by `process_generic_taxa()` during the cleaning pipeline.

**When to re-run:**

- You want to add, update, or remove entries from the dictionary.

**How to re-run:**

```r
# From the package root
source("data-raw/vernacular.R")
```

This overwrites `R/sysdata.rda`. Rebuild the package afterwards:

```r
devtools::build()
```

> The current dictionary covers ~80 vernacular terms across diatoms,
> dinoflagellates, coccolithophores, cyanobacteria, haptophytes,
> chrysophytes, and other microplankton groups.
