# JZ_Lucy

Project Lucy analyses.

## BL-251 Scar-seq AUC and 40 ng protein comparison

The BL-251 comparison workflow builds pairwise scatter plots comparing
Scar-seq time AUCs across tissues and cell lines with protein expression
measured at 40 ng and 24 h.

### Inputs

The default workflow expects the official Project Lucy SharePoint result
exports in `data/raw`:

| Dataset | SharePoint export |
| --- | --- |
| BL-251 in vivo | `LUCY-IVM-Results-BL251Scarseq-2026-11-20-JZ.csv` |
| BL-251 in vitro Scar-seq | `LUCY-IVM-Results-invitroBL251Scarseq-2026-11-20-JZ.csv` |
| BL-251 in vitro protein expression | `LUCY-IVB-TEST-001.json` |

Normalize those result exports before plotting:

```powershell
python .\scripts\prepare_bl251_sharepoint_inputs.py
python .\scripts\prepare_bl251_lucystudio_inputs.py
```

That writes the tidy `cyto_abund` tables configured in
`config/bl251_auc_inputs.yml`. If your Scar-seq export is already a tidy table
with `dataset`, `context`, `lnp`, `time_h`, `rep`, `embedded_barcode`,
`cytosol_reads`, `total_reads`, and `qpcr_norm_gl`, replace a `cyto_abund`
configuration path with your CSV:

```yaml
cyto_abund: path/to/tidy_cyto_abundance.csv
```

Protein-expression measurements are normalized from the LucyStudio 15 LNP
Pilot export. To use another protein-expression source, save a table at:

```text
data/processed/bl251_invitro_protein_expression.csv
```

with these columns:

```text
dataset,context,lnp,dose_ng,rep,protein_expression
```

The BL-251 LNP palette lives in `config/bl251_lnp_colors.csv`. Update that file
if a newer Lucy deck changes the canonical colors. The SharePoint normalization
keeps the 20-plex `86/56` Scar-seq rows for downstream AUCs and ignores the
125/250/500/1000-plex `86/56` rows. The ignored rows still contribute to the
per-sample cytosolic read denominator before they are filtered away.

### Run

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -e .[dev]
python .\scripts\prepare_bl251_sharepoint_inputs.py
python .\scripts\prepare_bl251_lucystudio_inputs.py
python .\scripts\plot_bl251_auc_pairs.py --config .\config\bl251_auc_inputs.yml
Rscript .\scripts\plot_bl251_auc_pairs.R
```

Outputs are written to:

| Output | Path |
| --- | --- |
| Scar-seq AUC summary | `outputs/bl251_auc_summary.csv` |
| Plotted Scar-seq AUC and 40 ng protein values | `outputs/bl251_comparison_values_40ng.csv` |
| In vitro Scar-seq AUC vs 40 ng protein PDF | `outputs/figures/bl251_auc_invitro_scarseq_vs_invitro_protein.pdf` |
| In vivo Scar-seq AUC vs 40 ng protein PDF | `outputs/figures/bl251_auc_invivo_scarseq_vs_invitro_protein.pdf` |
| 600 dpi PNGs | `outputs/figures/*.png` |

### Formula

For each sample, qPCR-normalized GL abundance is allocated across internal
barcodes by cytosolic read share:

```text
cyto_abund = norm_GL * cytosol_reads / sum(cytosol_reads in sample)
```

The workflow then sums the two internal barcodes for each LNP, averages
replicates at each time point, and computes trapezoid AUC per
`dataset / context / LNP`.

The in vitro protein-expression input uses the 24 h LucyStudio Alexa 488 mean
intensity measurement at 40 ng, averaged across replicate wells for each
cell context and LNP. It is a single-dose measurement, not an AUC. The
default plots compare the 40 ng protein expression with in vitro and in vivo
Scar-seq AUCs separately.
Each figure is rendered as a context grid with shared, log-scaled axes. In the
protein comparisons, only the Scar-seq axis represents an AUC.
