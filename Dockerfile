FROM mambaorg/micromamba:2.10.0@sha256:e1476182cd40eb1f7c0cfae4adec9bfd78f3d7f05c61368306abfdc8939c1198

ARG ASCAT_VERSION=3.2.0
ARG GLIMPSE_VERSION=2.0.0
ARG ICHORCNA_COMMIT=5bfc03ed854f0e93fe5b624c97c1290fa0053837
ARG VCS_REF=unknown

LABEL org.opencontainers.image.title="HAPLOID" \
      org.opencontainers.image.description="Haplotype-aware LOH and copy-number analysis of low-pass cfDNA sequencing" \
      org.opencontainers.image.source="https://github.com/vyellapa/HAPLOID" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.revision="${VCS_REF}"

USER root
COPY environment.yml /tmp/environment.yml
RUN micromamba config set use_sharded_repodata false \
    && micromamba install --yes --strict-channel-priority --name base --file /tmp/environment.yml \
        python=3.12 r-base=4.4 curl r-data.table r-doparallel r-foreach \
        r-rcolorbrewer bioconductor-genomicranges \
        bioconductor-hmmcopy bioconductor-genomeinfodb r-plyr r-optparse hmmcopy \
    && micromamba clean --all --yes

ENV PATH="/opt/conda/bin:${PATH}" \
    GLIMPSE_RUNNER=native \
    GLIMPSE_BIN=/opt/glimpse/bin \
    HAPLOID_RES=/opt/haploid/resources \
    REF=/opt/haploid/resources/hg38.fa \
    BAM_DIR=/data/bams \
    PYTHONUNBUFFERED=1 \
    MPLBACKEND=Agg \
    MPLCONFIGDIR=/tmp/matplotlib \
    XDG_CACHE_HOME=/tmp/cache \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1

SHELL ["/bin/bash", "-o", "pipefail", "-c"]
RUN curl --fail --location --retry 3 \
        "https://github.com/VanLoo-lab/ascat/archive/refs/tags/v${ASCAT_VERSION}.tar.gz" \
        --output /tmp/ascat.tar.gz \
    && tar -xzf /tmp/ascat.tar.gz -C /tmp \
    && R CMD INSTALL "/tmp/ascat-${ASCAT_VERSION}/ASCAT" \
    && mkdir -p /opt/src && mv /tmp/ascat.tar.gz "/opt/src/ascat-${ASCAT_VERSION}.tar.gz" \
    && rm -rf "/tmp/ascat-${ASCAT_VERSION}"

# ichorCNA (broadinstitute/ichorCNA, pinned commit) for BAMs without ichorCNA results; its panels of normals are
# removed (the pipeline runs ichorCNA without one), leaving the hg38 GC / mappability / centromere files.
# GPL-3 sources of ASCAT and ichorCNA are kept in /opt/src.
RUN curl --fail --location --retry 3 \
        "https://github.com/broadinstitute/ichorCNA/archive/${ICHORCNA_COMMIT}.tar.gz" --output /tmp/ichorCNA.tar.gz \
    && tar -xzf /tmp/ichorCNA.tar.gz -C /tmp \
    && R CMD INSTALL "/tmp/ichorCNA-${ICHORCNA_COMMIT}" \
    && install -m 755 "/tmp/ichorCNA-${ICHORCNA_COMMIT}/scripts/runIchorCNA.R" /opt/conda/bin/runIchorCNA.R \
    && rm -f "$(Rscript -e 'cat(system.file("extdata", package = "ichorCNA"))')"/*PoN* \
    && mv /tmp/ichorCNA.tar.gz "/opt/src/ichorCNA-${ICHORCNA_COMMIT}.tar.gz" \
    && rm -rf "/tmp/ichorCNA-${ICHORCNA_COMMIT}"

# Behind a TLS-inspecting proxy, pass its root certificate at build time (used for this download only, not stored):
#   docker build --secret id=extra_ca,src=proxy_root.pem ...
RUN --mount=type=secret,id=extra_ca,required=false \
    test "$(uname -m)" = x86_64 \
    && CA=() && if [ -s /run/secrets/extra_ca ]; then \
        cat /opt/conda/ssl/cacert.pem /run/secrets/extra_ca > /tmp/ca.pem && CA=(--cacert /tmp/ca.pem); fi \
    && mkdir -p /opt/glimpse/bin \
    && for tool in chunk split_reference phase ligate; do \
        curl "${CA[@]}" --fail --location --retry 3 \
            "https://github.com/odelaneau/GLIMPSE/releases/download/v${GLIMPSE_VERSION}/GLIMPSE2_${tool}_static" \
            --output "/opt/glimpse/bin/GLIMPSE2_${tool}_static" \
        && chmod 755 "/opt/glimpse/bin/GLIMPSE2_${tool}_static" || exit 1; \
    done \
    && rm -f /tmp/ca.pem

WORKDIR /opt/haploid
COPY scripts/ scripts/
COPY tests/ tests/
COPY container/ container/
COPY resources/README.md resources/fit_overrides.example.tsv resources/
COPY README.md LICENSE environment.yml ./
RUN chmod 755 scripts/*.sh container/*.sh \
    && mkdir -p work logs /data/bams \
    && chmod 777 work logs resources

USER $MAMBA_USER
RUN bash container/entrypoint.sh check \
    && python3 -m unittest discover -s tests -v

ENTRYPOINT ["/bin/bash", "/opt/haploid/container/entrypoint.sh"]
CMD ["--help"]
