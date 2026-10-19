"""Typed configuration objects for the pipeline."""
from dataclasses import dataclass, field
from typing import Literal, Optional


@dataclass(frozen=True)
class DataConfig:
    source: Literal["synthetic", "hmp", "local", "biom"] = "synthetic"
    path: Optional[str] = None
    seed: int = 42
    min_depth: int = 1000
    min_prevalence: float = 0.10
    taxon_rank: Literal["genus", "phylum", "species"] = "genus"


@dataclass(frozen=True)
class TransformConfig:
    method: Literal["relative", "clr", "rarefy"] = "clr"
    pseudocount: float = 0.5
    rarefy_depth: int = 5000


@dataclass(frozen=True)
class DAConfig:
    method: Literal["wilcoxon", "ancom_bc"] = "wilcoxon"
    group_col: str = "body_site"
    group_a: Optional[str] = None
    group_b: Optional[str] = None
    fdr_alpha: float = 0.05


@dataclass(frozen=True)
class MLConfig:
    target_col: str = "body_site"
    n_features: int = 50
    folds: int = 5
    seed: int = 42


@dataclass(frozen=True)
class PipelineConfig:
    data: DataConfig = field(default_factory=DataConfig)
    transform: TransformConfig = field(default_factory=TransformConfig)
    da: DAConfig = field(default_factory=DAConfig)
    ml: MLConfig = field(default_factory=MLConfig)
    out_dir: str = "results"
