import os
import re
import pysam
from pathlib import Path
from tqdm import tqdm
from collections import defaultdict
from .motif_utils import VcfIndexManager, MotifOptimizer
from .alignment_engine import align_lavs_strict_left
from ..utils.logging_utils import get_logger

logger = get_logger(__name__)

class TRAlignmentManager:
    def __init__(self, tr_db_path, output_dir):
        self.db = VcfIndexManager(tr_db_path)
        self.output_dir = Path(output_dir) / "alignment_results"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.scores_cache = {} 

    def process_groups(self, groups, config):
        chrom_counters = defaultdict(int)

        with tqdm(total=len(groups), desc="TR-Align Processing", unit="grp") as pbar:
            for grp in groups:
                chrom_counters[grp.chrom] += 1
                grp_index = chrom_counters[grp.chrom]
                
                # Group_{chrom}_{index}_{pos}
                grp_name = f"Group_{grp.chrom}_{grp_index}_{grp.variants[0].start}"
                
                try:
                    res = self._process_single_group(grp, grp_name)
                    if res:
                        self.scores_cache[grp_name] = res
                except Exception as e:
                    logger.warning(f"TR Align failed for {grp_name}: {e}")
                
                pbar.update(1)
        
        return self.scores_cache

    def _process_single_group(self, group, group_name):
        chrom = group.chrom
        start = group.start
        end = group.end
        
        match = self.db.search_region(chrom, start, end)
        if not match: return None
        _, tr_start, tr_end, motifs = match
        if not motifs: return None

        symbol_map = {} 
        for i, m in enumerate(motifs):
            symbol_map[m] = i 
        
        next_unk_id = 100 
        GAP_ID = -1
        
        id_to_seq_str = {} 
        for m, mid in symbol_map.items():
            id_to_seq_str[mid] = m
        id_to_seq_str[GAP_ID] = "-" 
        
        total_ratio = 0
        valid_vars = 0

        raw_seqs = {} 
        
        base_ref = group.variants[0].ref
        global_padding = base_ref[0]
        
        raw_seqs["seq0"] = base_ref

        seq_idx = 1
        for var in group.variants:
            for alt in var.alt:
                if "<" in alt or "[" in alt or "]" in alt: continue 
                seq_id = f"seq{seq_idx}"
                raw_seqs[seq_id] = alt
                seq_idx += 1

        seq_vectors = {}  
        
        for sid, full_seq in raw_seqs.items():
            if not full_seq:
                seq_vectors[sid] = []
                continue

            body_seq = full_seq[1:]
            
            vec = []
            if len(body_seq) > 0:
                segments, score = MotifOptimizer.decompose_structure(body_seq, motifs)
                if sid != "seq0": 
                    total_ratio += (score / len(body_seq))
                    valid_vars += 1
                
                for (seg_str, type_id) in segments:
                    if type_id >= 0:
                        vec.append(type_id)
                    else:
                        if seg_str not in symbol_map:
                            symbol_map[seg_str] = next_unk_id
                            id_to_seq_str[next_unk_id] = seg_str
                            next_unk_id += 1
                        vec.append(symbol_map[seg_str])
            
            seq_vectors[sid] = vec

        aligned_vectors = align_lavs_strict_left(seq_vectors)

        final_seqs = {}
        for sid, vec_body in aligned_vectors.items():
            str_parts = [global_padding] 
            for token in vec_body:
                if token == GAP_ID:
                    str_parts.append("-") 
                else:
                    str_parts.append(id_to_seq_str.get(token, "N"))
            final_seqs[sid] = str_parts

        seq_ids = list(final_seqs.keys())
        if not seq_ids: return None
        
        ali_len = len(final_seqs[seq_ids[0]])
        aligned_strings = {sid: "" for sid in seq_ids}
        
        for col in range(ali_len):
            max_w = 0
            for sid in seq_ids:
                token_str = final_seqs[sid][col]
                if token_str != "-":
                    max_w = max(max_w, len(token_str))
            
            if max_w == 0: max_w = 1
            
            for sid in seq_ids:
                token_str = final_seqs[sid][col]
                if token_str == "-":
                    aligned_strings[sid] += ("-" * max_w)
                else:
                    padding = max_w - len(token_str)
                    aligned_strings[sid] += (token_str + "-" * padding)

        out_name = f"{group_name}_aligned.fasta"
        out_path = self.output_dir / out_name

        sorted_ids = sorted(seq_ids, key=lambda x: int(x[3:]) if x[3:].isdigit() else -1)
        
        with open(out_path, 'w') as f:
            for sid in sorted_ids:
                f.write(f">{sid}\n{aligned_strings[sid]}\n")
        
        avg_score = total_ratio / valid_vars if valid_vars > 0 else 0
        return {"TR_SCORE": avg_score}