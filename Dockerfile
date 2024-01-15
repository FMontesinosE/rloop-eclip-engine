# Multi-stage Dockerfile for rloop-eclip-engine
# Uses uv for fast dependency resolution and installation

# Stage 1: Build dependencies with uv
FROM python:3.11-slim as builder

# Install uv package manager
RUN pip install --no-cache-dir uv==0.1.18

WORKDIR /build

# Copy dependency specifications
COPY pyproject.toml ./

# Install dependencies to virtual environment using uv
RUN uv venv /opt/venv && \
    . /opt/venv/bin/activate && \
    uv pip install --no-cache -e .

# Stage 2: Runtime image
FROM python:3.11-slim

# Install runtime dependencies (minimal)
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Copy virtual environment from builder
COPY --from=builder /opt/venv /opt/venv

# Copy application source
WORKDIR /app
COPY src/ ./src/
COPY pyproject.toml ./

# Activate virtual environment
ENV PATH="/opt/venv/bin:$PATH"
ENV PYTHONPATH="/app:$PYTHONPATH"

# Install package in editable mode
RUN pip install --no-cache-dir -e .

# Create data mount points
RUN mkdir -p /data /output

# Default command shows help
ENTRYPOINT ["rloop-eclip"]
CMD ["--help"]

# Metadata
LABEL maintainer="your.email@institution.edu"
LABEL description="High-performance genomic interval enrichment analysis for eCLIP and DRIPc-seq data"
LABEL version="0.2.1"
