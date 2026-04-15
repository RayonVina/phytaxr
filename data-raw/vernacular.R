# data-raw/vernacular.R
# Run once to regenerate R/sysdata.rda after editing the dictionary.

.VERNACULAR_TO_TAXON <- c(
  # Chlorophytes
  "chlorophytes" = "Chlorophyta",
  "chlorophyte" = "Chlorophyta",
  # Chrysophytes
  "chrysophytes" = "Chrysophyceae",
  "chrysophyte" = "Chrysophyceae",
  "crysophytes" = "Chrysophyceae",
  # Cryptophytes
  "cryptophytes" = "Cryptophyceae",
  "cryptophyte" = "Cryptophyceae",
  # Cyanobacteria
  "cyanophytes" = "Cyanobacteria",
  "cyanophyte" = "Cyanobacteria",
  "cyanobactereal" = "Cyanobacteria",
  "cyanobacteral" = "Cyanobacteria",
  # Dinoflagellates
  "dinoflagellatae" = "Dinoflagellata",
  "dinoflagellates" = "Dinoflagellata",
  "dinoflagellate" = "Dinoflagellata",
  "dinophyceae" = "Dinophyceae",
  # Euglenoids
  "euglenophytes" = "Euglenozoa",
  "euglenophyte" = "Euglenozoa",
  "euglenoid" = "Euglenozoa",
  "euglenophycota" = "Euglenozoa",
  "euglenophyceae" = "Euglenophyceae",
  # Haptophytes
  "haptophytes" = "Haptophyta",
  "haptophyte" = "Haptophyta",
  "prymnesiophytes" = "Prymnesiophyceae",
  "prymnesiophyte" = "Prymnesiophyceae",
  "prasinophytes" = "Prasinophyceae",
  "prasinophyte" = "Prasinophyceae",
  "praesinophyte" = "Prasinophyceae",
  # Raphidophytes
  "raphidophytes" = "Raphidophyceae",
  "raphidophyte" = "Raphidophyceae",
  # Pyrrophytes
  "pyrrophycophyta" = "Dinoflagellata",
  # Silicoflagellates
  "silicoflagellatae" = "Dictyochophyceae",
  "silicoflagellates" = "Dictyochophyceae",
  "silicoflagellate" = "Dictyochophyceae",
  # Coccolithophores
  "coccolithophores" = "Coccolithophyceae",
  "coccolithophore" = "Coccolithophyceae",
  "coccolithophorid" = "Coccolithophyceae",
  "coccolitophorids" = "Coccolithophyceae",
  "coccolithophorids" = "Coccolithophyceae",
  "heterococcolithophorid" = "Coccolithophyceae",
  "heterococcoliths" = "Coccolithophyceae",
  "holococcoliths" = "Coccolithophyceae",
  "holococcolithophores" = "Coccolithophyceae",
  "coccoliths" = "Coccolithophyceae",
  "coccolith" = "Coccolithophyceae",
  "coccospheres" = "Coccolithophyceae",
  "coccolithinae" = "Coccolithophyceae",
  # Diatoms
  "diatomaceae" = "Bacillariophyceae",
  "diatoms" = "Bacillariophyceae",
  "diatom" = "Bacillariophyceae",
  "naviculoids" = "Naviculales",
  "centric diatom" = "Coscinodiscophyceae",
  "pennate diatom" = "Bacillariophyceae",
  # Other groups
  "ebridians" = "Ebriida",
  "radiolarians" = "Radiolaria",
  "radiolarian" = "Radiolaria",
  "tintinnids" = "Tintinnida",
  "tintinnid" = "Tintinnida",
  "tintinnoideae" = "Tintinnida",
  "scuticociliates" = "Scuticociliatida",
  "peridinians" = "Peridiniales",
  "peridinioid" = "Peridiniales",
  "diplopsalioideae" = "Diplopsalidaceae",
  "infusoria" = "Ciliophora",
  "heliozoon" = "Heliozoa",
  "nanoflagellates" = "Flagellate",
  "nanoflagellate" = "Flagellate",
  "microflagellates" = "Flagellate",
  "microflagellate" = "Flagellate",
  "flagellates" = "Flagellate",
  "flagellate" = "Flagellate",
  "phytoflagellate" = "Flagellate",
  "gymnodimaceae" = "Gymnodiniaceae",
  "gynmodiniaceae" = "Gymnodiniaceae",
  "gynodiniaceae" = "Gymnodiniaceae",
  "halopappinae" = "Halopappidae",
  "warnowiaceae" = "Warnowiaceae"
)

usethis::use_data(.VERNACULAR_TO_TAXON, internal = TRUE, overwrite = TRUE)
