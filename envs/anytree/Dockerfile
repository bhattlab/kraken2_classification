FROM mambaorg/micromamba:1.4.3
COPY --chown=$MAMBA_USER:$MAMBA_USER anytree_docker.yaml /tmp/anytree_docker.yaml
RUN micromamba install -y -n base -f /tmp/anytree_docker.yaml && \
    micromamba clean --all --yes
ENV PATH="/opt/conda/bin:/opt/conda/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"