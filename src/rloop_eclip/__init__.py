"""
rloop-eclip-engine: High-performance genomic interval enrichment analysis.
"""

__version__ = "0.2.1"

from .engine import EnrichmentEngine
from .schemas import EnrichmentConfig, EnrichmentResult

__all__ = ["EnrichmentEngine", "EnrichmentConfig", "EnrichmentResult", "__version__"]
