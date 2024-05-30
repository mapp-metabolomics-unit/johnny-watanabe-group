
## Sirius

```bash
ssh commons-server
tmux a -n sirius
cd git_repos/mapp-metabolomics-unit/johnny-watanabe-group/docs/mapp_project_00002/mapp_batch_00003/results
```

### Login

``````bash
sirius login --user-env SIRIUS_USERNAME --password-env SIRIUS_PASSWORD
```


### Run sirius

(from the results of mapp_batch_00003)

sirius -i ./mzmine/mapp_batch_00003_sirius.mgf --output ./sirius/mapp_batch_00003 --maxmz 2500 config --IsotopeSettings.filter=true --FormulaSearchDB=BIO --Timeout.secondsPerTree=0 --FormulaSettings.enforced=HCNOP --Timeout.secondsPerInstance=0 --AdductSettings.detectable=[[M+H3N+H]+,[M+Na]+,[M-H4O2+H]+,[M+K]+,[M-H2O+H]+,[M+H]+] --UseHeuristic.mzToUseHeuristicOnly=650 --AlgorithmProfile=orbitrap --IsotopeMs2Settings=IGNORE --MS2MassDeviation.allowedMassDeviation=5.0ppm --NumberOfCandidatesPerIon=1 --UseHeuristic.mzToUseHeuristic=300 --FormulaSettings.detectable=Cl,Br,S --NumberOfCandidates=10 --ZodiacNumberOfConsideredCandidatesAt300Mz=10 --ZodiacRunInTwoSteps=true --ZodiacEdgeFilterThresholds.minLocalConnections=10 --ZodiacEdgeFilterThresholds.thresholdFilter=0.95 --ZodiacEpochs.burnInPeriod=2000 --ZodiacEpochs.numberOfMarkovChains=10 --ZodiacNumberOfConsideredCandidatesAt800Mz=50 --ZodiacEpochs.iterations=20000 --AdductSettings.enforced=, --AdductSettings.fallback=[[M+Na]+,[M+K]+,[M-H2O+H]+,[M+H]+] --FormulaResultThreshold=true --InjectElGordoCompounds=true --StructureSearchDB=BIO --RecomputeResults=false formula zodiac fingerprint structure canopus write-summaries

### Rerun sirius

sirius -i ./sirius/mapp_batch_00003 --output ./sirius/mapp_batch_00003 --maxmz 2500 config --IsotopeSettings.filter=true --FormulaSearchDB=BIO --Timeout.secondsPerTree=0 --FormulaSettings.enforced=HCNOP --Timeout.secondsPerInstance=0 --AdductSettings.detectable=[[M+H3N+H]+,[M+Na]+,[M-H4O2+H]+,[M+K]+,[M-H2O+H]+,[M+H]+] --UseHeuristic.mzToUseHeuristicOnly=650 --AlgorithmProfile=orbitrap --IsotopeMs2Settings=IGNORE --MS2MassDeviation.allowedMassDeviation=5.0ppm --NumberOfCandidatesPerIon=1 --UseHeuristic.mzToUseHeuristic=300 --FormulaSettings.detectable=Cl,Br,S --NumberOfCandidates=10 --ZodiacNumberOfConsideredCandidatesAt300Mz=10 --ZodiacRunInTwoSteps=true --ZodiacEdgeFilterThresholds.minLocalConnections=10 --ZodiacEdgeFilterThresholds.thresholdFilter=0.95 --ZodiacEpochs.burnInPeriod=2000 --ZodiacEpochs.numberOfMarkovChains=10 --ZodiacNumberOfConsideredCandidatesAt800Mz=50 --ZodiacEpochs.iterations=20000 --AdductSettings.enforced=, --AdductSettings.fallback=[[M+Na]+,[M+K]+,[M-H2O+H]+,[M+H]+] --FormulaResultThreshold=true --InjectElGordoCompounds=true --StructureSearchDB=BIO --RecomputeResults=false formula zodiac fingerprint structure canopus write-summaries


### Taxonomy handling

Install taxonomical-utils

```bash
pip install taxonomical_utils
```

Run the following command to get the taxonomy for the results:

1. resolve

```bash
taxonomical-utils resolve --input-file docs/fibl-pilot/pos/metadata/original/fibl_pilot_pos_metadata.tsv --output-file docs/fibl-pilot/pos/metadata/original/fibl_pilot_pos_metadata_resolved.csv --org-column-header source_taxon
```

2. retrieve upper taxa lineage
    
```bash
taxonomical-utils append-taxonomy --input-file docs/fibl-pilot/pos/metadata/original/fibl_pilot_pos_metadata_resolved.csv --output-file docs/fibl-pilot/pos/metadata/original/metadata_upper_taxa_lineage.csv
```

3. merge the taxonomy with the results

```bash
taxonomical-utils merge --input-file docs/fibl-pilot/pos/metadata/original/fibl_pilot_pos_metadata.tsv --resolved-taxa-file docs/fibl-pilot/pos/metadata/original/fibl_pilot_pos_metadata_resolved.csv --upper-taxa-lineage-file docs/fibl-pilot/pos/metadata/original/metadata_upper_taxa_lineage.csv --output-file docs/fibl-pilot/pos/metadata/original/fibl_pilot_pos_metadata.csv --org-column-header source_taxon --delimiter '\t'
```


### Met-annot-unifer




```bash
met-annot-unifier-cli align-horizontally --canopus-file /Users/pma/Dropbox/git_repos/COMMONS_Lab/EMI/fibl-metabolomics/docs/fibl-pilot/fibl_pilot_pos/results/sirius/canopus_compound_summary.tsv --gnps-file /Users/pma/Dropbox/git_repos/COMMONS_Lab/EMI/fibl-metabolomics/docs/fibl-pilot/fibl_pilot_pos/results/met_annot_enhancer/71042c319fa444d088e3704141a96354/nf_output/library/merged_results_with_gnps.tsv --gnps-mn-file /Users/pma/Dropbox/git_repos/COMMONS_Lab/EMI/fibl-metabolomics/docs/fibl-pilot/fibl_pilot_pos/results/met_annot_enhancer/71042c319fa444d088e3704141a96354/nf_output/networking/clustersummary_with_network.tsv --sirius-file /Users/pma/Dropbox/git_repos/COMMONS_Lab/EMI/fibl-metabolomics/docs/fibl-pilot/fibl_pilot_pos/results/sirius/compound_identifications.tsv --isdb-file /Users/pma/Dropbox/git_repos/COMMONS_Lab/EMI/fibl-metabolomics/docs/fibl-pilot/fibl_pilot_pos/results/met_annot_enhancer/fibl_pilot_pos_source_taxon/fibl_pilot_pos_source_taxon_spectral_match_results_repond_flat.tsv --output /Users/pma/Dropbox/git_repos/COMMONS_Lab/EMI/fibl-metabolomics/docs/fibl-pilot/fibl_pilot_pos/results/tmp/fibl_pilot_pos__met_annot_unified.tsv
```