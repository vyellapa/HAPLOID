# Reference resources

Nothing here is distributed with the code; build it once. Paths are relative to the repository root.

| Path | Source | How |
|---|---|---|
| `hg38-blacklist.v2.bed` | ENCODE blacklist v2 (Boyle lab, github.com/Boyle-Lab/Blacklist) | download `lists/hg38-blacklist.v2.bed.gz`, gunzip |
| `cytoband_b38.bed` | UCSC hg38 `cytoBand.txt.gz` | download, gunzip, keep as BED (chrom, start, end, band, stain) |
| `chain/hg19ToHg38.over.chain.gz` | UCSC liftOver chains | download from hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/ |
| `panel_hg38/1kGP_phase3.chr<N>.hg38.biSNP_maf0.01.bcf` | 1000 Genomes phase 3 phased GRCh37 VCFs | `PANEL_B37=<dir> REF=<hg38.fa> scripts/liftover_panel.sh <N>` then `scripts/panel_add_acan.sh <N>`, for N = 1..22 |
| `sites_hg38/chr<N>.sites.{vcf,tsv}.gz` | the lifted panel | `scripts/make_sites.sh <N>`, for N = 1..22 |
| `maps_b38/chr<N>.b38.gmap.gz` | GLIMPSE repository, `maps/genetic_maps.b38/` | download |
| `exclusion/` | GIAB GRCh38 genome stratifications (`GRCh38_segdups.bed.gz`, `GRCh38_AllTandemRepeats_ge101bp_slop5.bed.gz`), GIAB HG002 SV v0.6 `vdj_loci.GRCh38.bed`, Illumina DRAGEN `human.hg38.excl.bed` | download, then `REF=<hg38.fa> python3 scripts/build_exclusion.py <dir>` (prints every masked region for review) |
| `glimpse_bin/` | GLIMPSE2 v2.0.0 release | `scripts/fetch_glimpse.sh` |

Optional:

| Path | Made by | Purpose |
|---|---|---|
| `control_bin_mask_500kb.tsv`, `pipeline_bias_500kb.tsv` | `scripts/control_bin_mask.py` from healthy controls run through the same pipeline | mask of bins that misbehave in controls; bias subtraction only with `HAPLOID_CONTROL_BIAS=1` |
| `cohort_bias_500kb.tsv` | `scripts/cohort_bias.py` (second pass) | residual depth bias from the cohort's tumor-negative samples |
| `fit_overrides.tsv` | you (template: `fit_overrides.example.tsv`) | per-sample purity / ploidy / diploid constraints, applied to every joint fit |

Files derived from your own samples (control masks, cohort bias, constraints) describe your data; keep them out of
version control (`.gitignore` does).
