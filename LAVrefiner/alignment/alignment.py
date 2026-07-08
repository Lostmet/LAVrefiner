import re
import os
import time
import click
import shutil
import subprocess
import pysam
from tqdm import tqdm
from Bio import SeqIO
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from ..exceptions import AlignmentError
from ..utils.logging_utils import get_logger, log_tqdm_summary
from concurrent.futures import ProcessPoolExecutor, as_completed
from ..TR_align.motif_utils import VcfIndexManager, MotifOptimizer
from ..TR_align.alignment_engine import align_lavs_strict_left

logger = get_logger(__name__)

# === Global Cache for Worker Processes ===
_GLOBAL_TR_INDEX = None

class AlignmentResult:
    def __init__(self, group_id: str, sequences: Dict[str, str], 
                 tr_div: float = None, tr_type: str = None, rl: int = None):
        self.group_id = group_id
        self.sequences = sequences
        self.reference = sequences.get('reference', '')
        self.tr_div = tr_div
        self.tr_type = tr_type
        self.rl = rl

    @property
    def variant_count(self) -> int:
        return len(self.sequences) - 1  # Exclude reference

def read_fasta(fasta_file: str) -> Dict[str, List[str]]:
    """Read a FASTA file and group sequences by group header."""
    sequences = {}
    current_group = None
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id.startswith("Group"):
            current_group = record.id
            sequences[current_group] = [str(record.seq)]
        elif record.id.startswith("Variant") and current_group:
            sequences[current_group].append(str(record.seq))
    return sequences


def parse_fasta(fasta_file: str) -> Dict[str, str]:
    """Parse a FASTA file and return a dictionary of id → sequence."""
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences


def run_tr_core(sequences: List[str], motifs: List[str]) -> Tuple[Dict[str, str], Dict[str, str], float]:
    """
    Perform TR-aware alignment in memory.
    Output: (aligned_strings, aligned_tr_strings, avg_score)
    """
    # 1. Prepare Symbol Map
    symbol_map = {}
    for i, m in enumerate(motifs):
        symbol_map[m] = i

    next_unk_id = -2
    GAP_ID = -1
    id_to_seq_str = {mid: m for m, mid in symbol_map.items()}
    id_to_seq_str[GAP_ID] = "-"

    seq_vectors = {}
    total_ratio = 0.0
    valid_vars = 0
    
    # 2. Decompose Sequences
    for i, seq_str in enumerate(sequences):
        seq_id = f"seq{i}"
        vec = []
        if len(seq_str) > 0:
            segments, score = MotifOptimizer.decompose_structure(seq_str, motifs)

            if i > 0:
                eff_len = max(1, len(seq_str))
                total_ratio += (score / eff_len)
                valid_vars += 1

            for (seg_str, type_id) in segments:
                if type_id >= 0:
                    vec.append(type_id)
                else:
                    current_id = next_unk_id
                    next_unk_id -= 1 
                    id_to_seq_str[current_id] = seg_str
                    vec.append(current_id)
        
        seq_vectors[seq_id] = vec
    
    # 3. Align Vectors
    aligned_vectors = align_lavs_strict_left(seq_vectors)

    # 4. Reconstruct & Pad
    if not aligned_vectors: return {}, {}, 0.0
    
    seq_ids = list(aligned_vectors.keys())
    ali_len = len(aligned_vectors[seq_ids[0]])
    
    aligned_strings = {sid: "" for sid in seq_ids}
    aligned_tr_strings = {sid: "" for sid in seq_ids}

    for col in range(ali_len):
        # Pass 1: Calc max width for Standard Sequence (ACTG)
        max_w_seq = 0
        for sid in seq_ids:
            token_id = aligned_vectors[sid][col]
            s = "-" if token_id == GAP_ID else id_to_seq_str.get(token_id, "N")
            if s != "-": max_w_seq = max(max_w_seq, len(s))
        if max_w_seq == 0: max_w_seq = 1

        # Pass 2: Calc max width for TR Annotation
        max_w_tr = 0
        col_tr_tokens = {}
        for sid in seq_ids:
            token_id = aligned_vectors[sid][col]
            if token_id == GAP_ID:
                anno = "(-)"
            elif token_id < -1:
                u_num = abs(token_id) - 1
                anno = f"(u{u_num})"
            else:
                anno = f"({token_id})"
            col_tr_tokens[sid] = anno
            max_w_tr = max(max_w_tr, len(anno))

        # Pass 3: Build Strings
        for sid in seq_ids:
            token_id = aligned_vectors[sid][col]
            
            # Standard Output
            s = "-" if token_id == GAP_ID else id_to_seq_str.get(token_id, "N")
            pad_len = max_w_seq - len(s)
            if s == "-":
                aligned_strings[sid] += ("-" * max_w_seq)
            else:
                aligned_strings[sid] += (s + "-" * pad_len)

            # TR Output (IDs only)
            anno = col_tr_tokens[sid]
            pad_len_tr = max_w_tr - len(anno)
            aligned_tr_strings[sid] += (anno + " " * pad_len_tr)

    avg_score = total_ratio / valid_vars if valid_vars > 0 else 0.0
                
    return aligned_strings, aligned_tr_strings, avg_score

# ==========================================
# MAFFT Runner
# ==========================================
def run_mafft(threads: int, input_fasta: Path, output_fasta: Path, log_dir: Path = None, config=None) -> None:
    tmp_dir = Path(config.output_dir) / "mafft_tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    env = os.environ.copy()
    env["TMPDIR"] = str(tmp_dir)

    command = ["mafft", "--thread", "1", str(input_fasta)]
    max_retries = 5
    last_error_msg = ""

    for attempt in range(max_retries):
        try:
            result = subprocess.run(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True,
                env=env
            )

            with open(output_fasta, "w") as output_file:
                for line in result.stdout.splitlines():
                    if line.startswith(">"):
                        output_file.write(line + "\n")
                    else:
                        output_file.write(line.upper() + "\n")
            return

        except subprocess.CalledProcessError as e:
            last_error_msg = e.stderr
            logger.warning(f"MAFFT attempt {attempt + 1} failed for {input_fasta.name}. Retrying in 1s...")
            time.sleep(1)
        except FileNotFoundError as e:
            last_error_msg = f"Executable not found: {e}"
            logger.warning(f"MAFFT attempt {attempt + 1} failed for {input_fasta.name} (Not found). Retrying in 1s...")
            time.sleep(1)
        except Exception as e:
            last_error_msg = str(e)
            logger.warning(f"MAFFT attempt {attempt + 1} failed for {input_fasta.name}. Retrying in 1s...")
            time.sleep(1)

    if log_dir:
        log_dir.mkdir(parents=True, exist_ok=True)
        error_log = log_dir / f"{input_fasta.stem}_mafft_error.log"
        with open(error_log, "w") as log_file:
            log_file.write(f"MAFFT alignment failed after {max_retries} attempts for {input_fasta}:\n")
            log_file.write(last_error_msg)
            
    raise AlignmentError(
        f"MAFFT alignment failed for {input_fasta} after {max_retries} retries. Error log saved at: {error_log if log_dir else 'Not logged'}"
    )


def align_group(threads: str, group_name: str, sequences: List[str], align_dir: Path, log_dir: Path, has_insertion: bool, poly_ins_list, config, tr_vcf: str = None) -> AlignmentResult:
    """
    Align a single group.
    """
    output_fasta = align_dir / f"{group_name}_aligned.fasta"
    output_tr_file = align_dir / f"{group_name}_aligned.tr" 
    
    accumulated_scores = []

    try:
        parts = group_name.split('_')
        chrom = parts[1]
        pos = int(parts[3])
    except:
        pos = 0
        chrom = "unknown"

    global _GLOBAL_TR_INDEX
    tr_match = None
    
    tr_type = None
    rl = None
    
    if tr_vcf and _GLOBAL_TR_INDEX:
        tr_match = _GLOBAL_TR_INDEX.search_region(chrom, pos, pos + 1)
        if tr_match:
            _, _, _, motifs = tr_match
            if motifs and len(motifs) > 0:
                first_len = len(motifs[0])
                rl = first_len
                tr_type = "VNTR" if first_len > 6 else "STR"

    final_tr_div = None

    if has_insertion:
        input_fasta = align_dir / f"{group_name}_input_origin.fasta"
        origin_fasta = {f"seq{i}": seq for i, seq in enumerate(sequences)}
        origin_tr_fasta = {f"seq{i}": seq for i, seq in enumerate(sequences)}

        relevant_ins = [item for item in poly_ins_list if item['pos'] == pos]

        for i, ins in enumerate(relevant_ins):
            slice_seqs = []
            for j in range(len(sequences)):
                s = sequences[j][ins['start']:ins['end']].replace("-", "")
                slice_seqs.append(s)

            use_tr_algo = False
            aligned_slice_dict = {}
            tr_slice_dict = {}
            slice_score = 0.0

            if tr_match:
                _, _, _, motifs = tr_match
                if motifs:
                    try:
                        aligned_slice_dict, tr_slice_dict, slice_score = run_tr_core(slice_seqs, motifs)
                        accumulated_scores.append(slice_score)
                        use_tr_algo = True
                    except Exception as e:
                        logger.warning(f"TR Align failed for {group_name}, fallback to MAFFT: {e}")
                        use_tr_algo = False

            if not use_tr_algo:
                input_sliced = align_dir / f"{group_name}_input_sliced_{i+1}.fasta"
                output_sliced = align_dir / f"{group_name}_aligned_sliced_{i+1}.fasta"

                with open(input_sliced, "w") as f:
                    for j, s in enumerate(slice_seqs):
                        f.write(f">seq{j}\n{s}\n")

                run_mafft(threads, input_sliced, output_sliced, log_dir=log_dir, config=config)
                aligned_slice_dict = parse_fasta(output_sliced)
                tr_slice_dict = aligned_slice_dict

            try:
                for key in origin_fasta:
                    if key in aligned_slice_dict:
                        ori_list = list(origin_fasta[key])
                        ori_list[ins['start']:ins['end']] = list(aligned_slice_dict[key])
                        origin_fasta[key] = "".join(ori_list)
                        
                        ori_tr_list = list(origin_tr_fasta[key])
                        tr_segment = tr_slice_dict.get(key, "")
                        ori_tr_list[ins['start']:ins['end']] = list(tr_segment)
                        origin_tr_fasta[key] = "".join(ori_tr_list)
                        
            except Exception as e:
                click.echo(f"Warning: Stitching failed in {group_name}: {e}")

        if accumulated_scores:
            final_tr_div = sum(accumulated_scores) / len(accumulated_scores)

        with open(output_fasta, "w") as f:
            for idx in range(len(sequences)):
                sid = f"seq{idx}"
                if sid in origin_fasta:
                    f.write(f">{sid}\n{origin_fasta[sid]}\n")
        
        if tr_type:
            with open(output_tr_file, "w") as f:
                for idx in range(len(sequences)):
                    sid = f"seq{idx}"
                    if sid in origin_tr_fasta:
                        f.write(f">{sid}\n{origin_tr_fasta[sid]}\n")
    else:
        with open(output_fasta, "w") as f:
            for i, seq in enumerate(sequences):
                f.write(f">seq{i}\n{seq}\n")
        

    aligned_sequences = {}
    if output_fasta.exists():
        with open(output_fasta) as f:
            current_id = None
            current_seq = []
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_id:
                        aligned_sequences[current_id] = "".join(current_seq)
                    current_id = line[1:].strip()
                    current_seq = []
                else:
                    current_seq.append(line)
            if current_id:
                aligned_sequences[current_id] = "".join(current_seq)

    return AlignmentResult(
        group_id=group_name,
        sequences=aligned_sequences,
        tr_div=final_tr_div,
        tr_type=tr_type, 
        rl=rl            
    )

def run_alignments(config, fasta_file: str, has_insertion_dict: Dict[str, bool], poly_ins_list, tr_vcf: str = None) -> List[AlignmentResult]:
    """
    Run alignments for all groups using MAFFT or TR-Motif aligner.
    """
    group_sequences = read_fasta(fasta_file)

    align_dir = Path(config.output_dir) / "alignment_results"
    log_dir = Path(config.output_dir) / "alignment_error_logs"
    align_dir.mkdir(parents=True, exist_ok=True)

    ordered_items = list(group_sequences.items())
    ordered_items.sort(key=lambda x: len(x[1]) * len(x[1][0]), reverse=True)
    total = len(ordered_items)
    results = [None] * total

    global _GLOBAL_TR_INDEX
    if tr_vcf:
        if _GLOBAL_TR_INDEX is None:
            try:
                _GLOBAL_TR_INDEX = VcfIndexManager(tr_vcf)
            except Exception as e:
                logger.error(f"Failed to pre-load TR index: {e}")

    with ProcessPoolExecutor(max_workers=config.threads) as executor, tqdm(
        total=total, desc="Processed Groups", unit='groups'
    ) as pbar:
        futures = {
            executor.submit(
                align_group,
                config.threads,
                ordered_items[idx][0],
                ordered_items[idx][1],
                align_dir,
                log_dir,
                has_insertion_dict.get(ordered_items[idx][0], False),
                poly_ins_list,
                config,
                tr_vcf 
            ): idx
            for idx in range(total)
        }

        for future in as_completed(futures):
            idx = futures[future]
            group_name = ordered_items[idx][0]
            try:
                res = future.result()
                results[idx] = res
            except Exception as e:
                if not log_dir.exists():
                    log_dir.mkdir(parents=True, exist_ok=True)
                with open(log_dir / f"{group_name}_error.log", "a") as log_file:
                    log_file.write(f"Alignment failed for {group_name}: {str(e)}\n")
            finally:
                pbar.update(1)

    pbar.close()
    log_tqdm_summary(pbar, logger)
    
    fasta_path = Path(config.output_dir) / "variants_pre_aligned.fasta"
    fasta_path.unlink(missing_ok=True)

    tmp_dir = Path(config.output_dir) / "mafft_tmp"
    if tmp_dir.exists() and tmp_dir.is_dir():
        shutil.rmtree(tmp_dir)

    return results