# SVrefiner/config.py

import os
from pathlib import Path
from typing import Optional
from dataclasses import dataclass

@dataclass
class Config:
    """Configuration for the pipeline"""
    output_dir: str
    threads: int = 1
    window_size: int = 4
    vcf_file: str = None
    ref_fasta: str = None
    tr_vcf: str = None 
    grouped_variants_file: str = None
    alignments_dir: Optional[str] = None
    
    def __post_init__(self):
        """Validate inputs and create output directory"""
        if self.vcf_file and not os.path.exists(self.vcf_file):
            raise FileNotFoundError(f"VCF file not found: {self.vcf_file}")
        if self.ref_fasta and not os.path.exists(self.ref_fasta):
            raise FileNotFoundError(f"Reference FASTA not found: {self.ref_fasta}")
        if self.tr_vcf and not os.path.exists(self.tr_vcf):
            raise FileNotFoundError(f"TR Database VCF not found: {self.tr_vcf}")
        if self.grouped_variants_file and not os.path.exists(self.grouped_variants_file):
            raise FileNotFoundError(f"Grouped variants file not found: {self.grouped_variants_file}")
        os.makedirs(self.output_dir, exist_ok=True)