FROM mambaorg/micromamba:2.10.0@sha256:e1476182cd40eb1f7c0cfae4adec9bfd78f3d7f05c61368306abfdc8939c1198

ARG ASCAT_VERSION=3.2.0
ARG GLIMPSE_VERSION=2.0.0
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
    && micromamba clean --all --yes

ENV PATH="/opt/conda/bin:${PATH}" \
    GLIMPSE_RUNNER=native \
    GLIMPSE_BIN=/opt/glimpse/bin \
    HAPLOID_RES=/opt/haploid/resources \
    REF=/opt/haploid/resources/hg38.fa \
    BAM_DIR=/data/bams \
    PYTHONUNBUFFERED=1 \
    MPLBACKEND=Agg \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1

SHELL ["/bin/bash", "-o", "pipefail", "-c"]
RUN curl --fail --location --retry 3 \
        "https://github.com/VanLoo-lab/ascat/archive/refs/tags/v${ASCAT_VERSION}.tar.gz" \
        --output /tmp/ascat.tar.gz \
    && tar -xzf /tmp/ascat.tar.gz -C /tmp \
    && R CMD INSTALL "/tmp/ascat-${ASCAT_VERSION}/ASCAT" \
    && rm -rf /tmp/ascat.tar.gz "/tmp/ascat-${ASCAT_VERSION}"

RUN test "$(uname -m)" = x86_64 \
    && mkdir -p /opt/glimpse/bin \
    && for tool in chunk split_reference phase ligate; do \
        curl --fail --location --retry 3 \
            "https://github.com/odelaneau/GLIMPSE/releases/download/v${GLIMPSE_VERSION}/GLIMPSE2_${tool}_static" \
            --output "/opt/glimpse/bin/GLIMPSE2_${tool}_static" \
        && chmod 755 "/opt/glimpse/bin/GLIMPSE2_${tool}_static" || exit 1; \
    done

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
