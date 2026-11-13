# HAPLOID container

The Linux amd64 image includes the HAPLOID scripts and their Python/R dependencies,
ASCAT 3.2.0, DNAcopy, samtools, bcftools, htslib, liftOver, and the four GLIMPSE2
2.0.0 tools. GLIMPSE runs directly inside this image. Docker is not needed inside
the container.

Docker Hub image: [`vyellapa/haploid:latest`](https://hub.docker.com/r/vyellapa/haploid).

Two images:

- `vyellapa/haploid:latest` - code and tools (including ichorCNA, pinned commit, with its panels of normals removed).
  Reference data and sequencing data are supplied as mounts, as described below.
- `vyellapa/haploid:full` - the same plus the public reference data (~15 GB): hg38, the prepared GLIMPSE2 1000
  Genomes reference, panel sites, genetic maps and region exclusions. `haploid run` goes from a BAM to final
  results with one `/out` mount; see the main README's quick start. It holds nothing derived from samples (no
  healthy-control mask, cohort bias, constraints or panel of normals). Users build their own from their healthy
  controls with `haploid controls` and pass the folder to `run --controls` (main README).

Neither image includes BAMs, existing results or credentials. The fragment-count alternative to ichorCNA depth
requires a prepared panel of normals (see below).

## Build the full image

From a machine with the prepared resources (`resources/`, a complete GLIMPSE2 reference with its `.done.chr*`
markers, hg38.fa + .fai) and the base image built:

```bash
container/bundle_resources.sh resources/ work/glimpse/ref /path/to/hg38.fa /tmp/haploid_bundle   # hard links
docker build --platform linux/amd64 -f container/Dockerfile.full \
  --build-arg BASE=vyellapa/haploid:latest -t vyellapa/haploid:full /tmp/haploid_bundle
```

`bundle_resources.sh` copies an allow-list only and refuses an exclusion set built with the DRAGEN file.

## Build and publish

From the repository root, with Docker installed:

```bash
export HAPLOID_IMAGE=YOUR_DOCKERHUB_USERNAME/haploid:latest
docker build --platform linux/amd64 -t "$HAPLOID_IMAGE" .
docker run --rm --platform linux/amd64 "$HAPLOID_IMAGE" check
docker login
docker push "$HAPLOID_IMAGE"
```

The build runs dependency checks and the synthetic regression tests. The static
GLIMPSE release supports x86-64, so the image is built for `linux/amd64`.

## Mount data and run

Prepare `resources/` as described in [the reference guide](../resources/README.md).
Include `hg38.fa` and `hg38.fa.fai` for the same hg38 reference used to align the
BAMs. Each input BAM must be indexed. The expected input layout is:

```text
/path/to/bams/SAMPLE_500k/SAMPLE.bam
/path/to/bams/SAMPLE_500k/SAMPLE.bam.bai
/path/to/bams/SAMPLE_500k/SAMPLE_Y.correctedDepth.txt
/path/to/bams/SAMPLE_500k/SAMPLE_Y.params.txt
/path/to/bams/SAMPLE_500k/SAMPLE_Y.seg.txt
```

Define a reusable command in Bash, from the repository root:

```bash
export HAPLOID_IMAGE=vyellapa/haploid:latest
mkdir -p work logs
HAPLOID_RUN=(docker run --rm --init --platform linux/amd64
  --user "$(id -u):$(id -g)"
  -v "$PWD/resources:/opt/haploid/resources"
  -v "$PWD/work:/opt/haploid/work"
  -v "$PWD/logs:/opt/haploid/logs"
  -v /path/to/bams:/data/bams:ro
  -e CORES=4 -e JOBS=1 -e PREP_JOBS=1
  "$HAPLOID_IMAGE")

# All samples: shared reference preparation, duplicate marking, likelihoods,
# phasing, haplotype counts, depth/bins, segmentation, joint fits, and reports.
"${HAPLOID_RUN[@]}" cohort

# Or one sample, including shared reference preparation and reports:
"${HAPLOID_RUN[@]}" sample /data/bams/SAMPLE_500k/SAMPLE.bam

# Generate genome-wide comparison plots after a successful run:
"${HAPLOID_RUN[@]}" plot SAMPLE
```

The `resources` mount must be writable: shared preparation can add panel AC/AN
tags, and cohort calibration writes bias tables there. `work` and `logs` must also
be writable. Using your UID/GID keeps generated files owned by you. Results persist
in `work/` and runs resume from its completion markers. A failed or incomplete
sample makes the container exit with a nonzero status.

The defaults are one sample at a time, four sample cores, and one reference
preparation job. Each reference preparation job itself uses four threads; increase
these settings according to available memory and CPU. Other existing `HAPLOID_*`,
`ICHOR_*`, `MODES`, and `KEEP_BAM` settings can be passed with Docker's `-e` option.
All file paths in environment variables must be container paths.

Without ichorCNA corrected-depth files, put `ichor_pon_hg38_500kb.tsv` in
`work/depth/` and add `-e ALLOW_FRAGMENT_DEPTH=1` before the image name. The runner
then builds the shared GC table and uses its fragment-count depth method.

## Second pass and constraints

Cohort calibration is a separate command so its masked regions can be reviewed
before refitting, as in the existing workflow. It needs at least 30 suitable
tumor-negative samples by default.

```bash
"${HAPLOID_RUN[@]}" bias
# Review resources/cohort_bias_review.tsv and work/cohort_bias.png.
"${HAPLOID_RUN[@]}" refit SAMPLE1 SAMPLE2
"${HAPLOID_RUN[@]}" resolve

# Known sample values, or rebuild reports without rerunning the fits:
"${HAPLOID_RUN[@]}" constrain SAMPLE1 purity=0.4-0.6 ploidy=2.2-2.7
"${HAPLOID_RUN[@]}" report
```

All original scripts are available through `exec`, including panel preparation:

```bash
"${HAPLOID_RUN[@]}" exec bash scripts/make_sites.sh 1
"${HAPLOID_RUN[@]}" exec python3 -m unittest discover -s tests -v
```

## Apptainer

The same Docker Hub image can be used on a cluster:

```bash
export HAPLOID_IMAGE=vyellapa/haploid:latest
apptainer pull haploid.sif "docker://$HAPLOID_IMAGE"
apptainer run --cleanenv \
  --bind "$PWD/resources:/opt/haploid/resources" \
  --bind "$PWD/work:/opt/haploid/work" \
  --bind "$PWD/logs:/opt/haploid/logs" \
  --bind /path/to/bams:/data/bams:ro \
  --env CORES=4,JOBS=1,PREP_JOBS=1 \
  haploid.sif cohort
```
