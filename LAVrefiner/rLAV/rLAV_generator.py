import os
import re
import click
import numpy as np
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm
from pathlib import Path
from collections import OrderedDict
from ..utils.logging_utils import get_logger, log_tqdm_summary
from concurrent.futures import ThreadPoolExecutor, as_completed

logger = get_logger(__name__)

def save_D_matrix(D, seq_ids, lav_ids_full, out_path):
    lav_ids = [x.split("_")[-1] for x in lav_ids_full]
    df = pd.DataFrame(D, columns=lav_ids)
    df.insert(0, "lav_id", seq_ids) 
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    df.to_csv(out_path, index=False)

# === Encoding Map ===
base_map = {"-": 0, "A": 1, "a": 1, "T": 2, "t": 2, "C": 3, "c": 3, "G": 4, "g": 4, "N": 5, "n": 5}
reverse_map = {0: "-", 1: "A", 2: "T", 3: "C", 4: "G", 5: "N"}

def encode_sequence(seq):
    return np.array([base_map.get(base.upper(), 5) for base in seq], dtype=int)

def num_to_base(n):
    return reverse_map.get(n, "-")

# === FASTA Input Utilities ===
def parse_filename_metadata(filename):
    pattern = r"Group_(\d+)_(\d+)_(\d+)_aligned\.fasta"
    match = re.search(pattern, filename)
    if match:
        return {
            "chrom": int(match.group(1)),
            "group_num": int(match.group(2)),
            "pos": int(match.group(3))
        }
    return None

def list_aligned_fasta_files(folder):
    aligned_files = []
    if not os.path.exists(folder):
        return []
    for fname in os.listdir(folder):
        if fname.endswith("_aligned.fasta"):
            full_path = os.path.join(folder, fname)
            meta = parse_filename_metadata(fname)
            if meta:
                aligned_files.append({"path": full_path, "meta": meta})
    return aligned_files

def preload_fasta_files(file_infos):
    in_memory_data = []
    for file_info in file_infos:
        records = list(SeqIO.parse(file_info["path"], "fasta"))
        info_copy = file_info.copy()
        info_copy["records"] = records
        in_memory_data.append(info_copy)
    return in_memory_data

# === Block Construction Utilities ===
def is_rank1_block(mat):
    if np.all(mat == 0): return False
    return np.linalg.matrix_rank(mat) == 1

def is_rank1_pair(col1, col2):
    return np.linalg.matrix_rank(np.stack([col1, col2], axis=1)) == 1

def find_blocks_with_rank(matrix, window=10):
    blocks = []
    i = 0
    while i < matrix.shape[1]:
        if np.all(matrix[:, i] == 0):
            blocks.append((i, i))
            i += 1
            continue
        end = min(i + window, matrix.shape[1])
        submat = matrix[:, i:end]
        if is_rank1_block(submat):
            blocks.append((i, end - 1))
            i = end
        else:
            start = i
            i += 1
            while i < matrix.shape[1]:
                if np.all(matrix[:, i] == 0):
                    blocks.append((start, i - 1))
                    blocks.append((i, i))
                    i += 1
                    break
                if not is_rank1_pair(matrix[:, i - 1], matrix[:, i]):
                    blocks.append((start, i - 1))
                    start = i
                i += 1
            else:
                if start < matrix.shape[1]:
                    blocks.append((start, matrix.shape[1] - 1))
    return blocks

def find_blocks_mask_based(matrix):
    mask = matrix != 0  
    is_rank1 = np.all(mask[:, :-1] == mask[:, 1:], axis=0) 
    blocks = []
    start = 0
    for i, related in enumerate(is_rank1):
        if not related:
            blocks.append((start, i))
            start = i + 1
    blocks.append((start, matrix.shape[1] - 1)) 
    return blocks

def smart_split_by_row_patterns(mat):
    unique_patterns, inverse_indices = np.unique(mat, axis=0, return_inverse=True)
    sub_matrices = []
    for idx in range(unique_patterns.shape[0]):
        mask = (inverse_indices == idx)
        sub_mat = np.zeros_like(mat)
        sub_mat[mask] = mat[mask]
        if np.any(sub_mat != 0):
            sub_matrices.append(sub_mat)
    return sub_matrices

def extract_blocks_with_split(matrix, blocks, start_index_shift=True):
    expanded, new_blocks = [], []
    for (start, end) in blocks:
        mat = matrix[:, start:end + 1]
        if np.unique(mat, axis=0).shape[0] <= 2:
            expanded.append(mat)
            new_blocks.append((start, end))
        else:
            sub_mats = smart_split_by_row_patterns(mat)
            for i, submat in enumerate(sub_mats):
                new_start = start + i * mat.shape[1] if start_index_shift else 0
                new_end = new_start + mat.shape[1] - 1
                expanded.append(submat)
                new_blocks.append((new_start, new_end))
    return expanded, new_blocks

# === Matrix and Metadata Construction ===
def convert(row):
    return ''.join([num_to_base(x) for x in row if x != 0])

def build_D_and_meta(expanded_mats):
    D_cols, metas = [], []
    ref_prefix, pos, pos_buffer = "", 0, 0
    for i, mat in enumerate(expanded_mats):
        first_col = mat[:, 0]
        last_col = mat[:, -1]
        if i == 0:
            ref_prefix = num_to_base(last_col[0])
            continue
        d_col = (first_col != 0).astype(int) if first_col[0] == 0 else (first_col == 0).astype(int)
        D_cols.append(d_col[:, None])
        row_idx = np.where(d_col == 1)[0]
        if len(row_idx) == 0: continue
        ref = ref_prefix + convert(mat[0])
        alt = ref_prefix + convert(mat[row_idx[0]])
        pos_buffer = pos
        if last_col[0] != 0:
            ref_prefix = num_to_base(last_col[0])
            pos += mat.shape[1]
        metas.append({"pos": pos_buffer, "ref": ref, "alt": alt})
    
    if not D_cols: return np.array([]), []
    D = np.hstack(D_cols)[1:, :]
    return D, metas


def merge_rlav_by_pos_and_mask(D: np.ndarray, metas: list):
    """
    Merge adjacent rLAVs that share the same genotype mask.
    Crucial for TRs: CAG -> merges C, A, G columns into one.
    """
    if D.size == 0 or len(metas) == 0:
        return D, metas

    # Group by POS first (usually adjacent columns from one block share a pos buffer in build_D_and_meta)
    pos_to_cols = OrderedDict()
    for j, m in enumerate(metas):
        p = int(m["pos"])
        pos_to_cols.setdefault(p, []).append(j)

    keep_cols = []
    metas_new = []

    for p in sorted(pos_to_cols.keys()):
        cols = pos_to_cols[p] 
        # Group by Mask within the same pos
        mask_groups = OrderedDict() 
        for j in cols:
            mask = tuple(int(x) for x in D[:, j].tolist())
            mask_groups.setdefault(mask, []).append(j)

        for mask, js in mask_groups.items():
            js_sorted = sorted(js)            
            j0 = js_sorted[0]                 
            merged_ref = metas[j0]["ref"]
            merged_alt = metas[j0]["alt"]
            # Concatenate ALTs
            for j in js_sorted[1:]:
                aj = metas[j]["alt"]
                merged_alt += (aj[1:] if len(aj) > 0 else "")
            keep_cols.append(j0)
            metas_new.append({"pos": p, "ref": merged_ref, "alt": merged_alt})

    D_new = D[:, keep_cols] if len(keep_cols) > 0 else D
    return D_new, metas_new


def process_fasta_in_memory(file_info, has_insertion_dict, tr_group_set):
    records = file_info["records"]
    pos_prefix = file_info["meta"]["pos"]
    fasta_path = file_info["path"]
    group_name = re.sub(r'_aligned\.fasta$', '', os.path.basename(fasta_path))

    is_tr = group_name in tr_group_set
    has_insertion = has_insertion_dict.get(group_name, False)
    
    should_merge = has_insertion or is_tr

    seqs = [str(rec.seq) for rec in records]
    matrix = np.vstack([encode_sequence(seq) for seq in seqs])

    blocks = find_blocks_with_rank(matrix) if should_merge else find_blocks_mask_based(matrix)
    
    expanded_mats, _ = extract_blocks_with_split(matrix, blocks)

    D, metas = build_D_and_meta(expanded_mats)

    for m in metas:
        m["pos"] += pos_prefix

    if should_merge:
        D, metas = merge_rlav_by_pos_and_mask(D, metas)

    return os.path.basename(fasta_path), metas, D

# === Output & Main Pipeline ===

def save_meta_csv(meta_list, out_path):
    meta_tuples = [(item["group_name"], eval(item["meta_array"])) for item in meta_list]
    
    def sort_key(x):
        gname = x[0]
        m = re.search(r'Group_(\d+)_(\d+)_(\d+)_rLAV(\d+)', gname)
        if m:
            chrom, idx, pos, rlav_i = map(int, m.groups())
            return (chrom, idx, pos, rlav_i)
        return (999, 999, 999, 999)

    meta_tuples.sort(key=sort_key)
    sorted_meta = [{"group_name": g, "meta_array": str(m)} for g, m in meta_tuples]
    df = pd.DataFrame(sorted_meta)
    df.to_csv(out_path, index=False)

def rLAV_generator(out: str, has_insertion_dict: dict, threads: int = 16, tr_metadata: dict = None):
    matrix_dir = Path(out) / "matrix_results"
    matrix_dir.mkdir(parents=True, exist_ok=True)

    alignments_dir = Path(out) / "alignment_results"
    all_files = list_aligned_fasta_files(alignments_dir)
    in_memory_data = preload_fasta_files(all_files)

    tr_group_set = set(tr_metadata.keys()) if tr_metadata else set()

    results = []
    
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(process_fasta_in_memory, item, has_insertion_dict, tr_group_set): item for item in in_memory_data}
        pbar = tqdm(as_completed(futures), total=len(futures), desc="rLAV generating", unit="groups")
        for future in pbar:
            try:
                fname, metas, D = future.result()
                meta_info = futures[future]["meta"]
                results.append({
                    "source": fname,
                    "chrom": meta_info["chrom"],
                    "group": meta_info["group_num"],
                    "original_pos": meta_info["pos"],
                    "meta": metas,
                    "D": D
                })
            except Exception as e:
                logger.error(f"Failed to process file: {futures[future]['path']} — {repr(e)}")

    log_tqdm_summary(pbar, logger)

    meta_records = []
    for item in results:
        fname = item["source"]
        metas = item["meta"]
        D = item["D"]
        group_name = fname.replace("_aligned.fasta", "")

        tr_info = None
        if tr_metadata and group_name in tr_metadata:
            tr_info = tr_metadata[group_name] 
        
        lav_ids_full = [f"{group_name}_rLAV{i+1}" for i in range(D.shape[1])]
        seq_ids = [f"seq{i+1}" for i in range(D.shape[0])]
        
        save_D_matrix(D, seq_ids, lav_ids_full, matrix_dir / f"{group_name}_D_matrix.csv")
        
        for i, meta in enumerate(metas, start=1):
            if tr_info:
                meta.update(tr_info)
            
            meta_records.append({
                "group_name": f"{group_name}_rLAV{i}",
                "meta_array": str(meta)
            })

    try:
        save_meta_csv(meta_records, Path(out) / "rLAV_meta.csv")
    except Exception as e:
        logger.error(f"Failed to save rLAV_meta.csv: {repr(e)}")