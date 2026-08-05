# containers/Dockerfile.py
#
# Python environment: single-cell GPU path (with Scanpy fallback), clustering
# agreement metrics, version collection, provenance helpers.
#
# Used by: RAPIDS_SINGLECELL, CLUSTERING_AGREEMENT, SOFTWARE_VERSIONS

FROM python:3.11-slim

LABEL org.opencontainers.image.title="kazrna-py" \
      org.opencontainers.image.description="Python environment for KazRNA-Pipe" \
      org.opencontainers.image.source="https://github.com/ORG/kazrna-pipe" \
      org.opencontainers.image.licenses="MIT"

RUN apt-get update && apt-get install -y --no-install-recommends \
        procps git build-essential \
    && rm -rf /var/lib/apt/lists/*

# Versions pinned so the image is reproducible from the Dockerfile alone.
RUN pip install --no-cache-dir \
        numpy==1.26.4 \
        pandas==2.2.2 \
        scipy==1.13.1 \
        scikit-learn==1.5.1 \
        matplotlib==3.9.1 \
        anndata==0.10.8 \
        scanpy==1.10.2 \
        harmonypy==0.0.10 \
        igraph==0.11.6 \
        leidenalg==0.10.2 \
        h5py==3.11.0 \
        PyYAML==6.0.1

RUN python -c "import numpy, pandas, scipy, sklearn, matplotlib, anndata, scanpy, yaml, h5py; \
    print('All packages present')"

CMD ["python"]
