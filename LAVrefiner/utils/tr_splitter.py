import os
import sys
import gzip
import pysam
from collections import defaultdict


def smart_open(filename, mode='rt'):
    if filename.endswith(".gz"):
        return gzip.open(filename, mode, encoding='utf-8', errors='replace')
    else:
        return open(filename, mode, encoding='utf-8', errors='replace')

def normalize_chrom(chrom):
    if chrom.startswith("chr"):
        return chrom[3:]
    return chrom

def parse_info(info_str):
    d = {}
    for part in info_str.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            d[k] = v
    return d

def load_tr_database(db_path):
    print(f"[INFO] Loading TRCompDB: {db_path}")
    tr_index = defaultdict(lambda: defaultdict(list))
    count = 0
    BIN_SIZE = 10000
    
    try:
        with smart_open(db_path, 'rt') as f:
            for line in f:
                if line.startswith("#"): continue
                
                parts = line.strip().split("\t")
                if len(parts) < 8: continue
                
                chrom = normalize_chrom(parts[0])
                try:
                    start = int(parts[1])
                    info = parse_info(parts[7])
                    
                    if "END" not in info: continue
                    end = int(info["END"])
                    
                    start_bin = start // BIN_SIZE
                    end_bin = end // BIN_SIZE
                    record = (start, end)
                    
                    for b in range(start_bin, end_bin + 1):
                        tr_index[chrom][b].append(record)
                    count += 1
                except ValueError:
                    continue
    except FileNotFoundError:
        print(f"[ERROR] Database file not found: {db_path}")
        sys.exit(1)
            
    print(f"[INFO] DB Loaded. {count} regions.")
    return tr_index

def compress_and_index(vcf_path):
    abs_vcf_path = os.path.abspath(vcf_path)
    
    if not os.path.exists(abs_vcf_path):
        print(f"[ERROR] File missing before indexing: {abs_vcf_path}")
        return None
    
    print(f"       Indexing {vcf_path} ...", end="\r")
    
    try:
        pysam.tabix_index(abs_vcf_path, preset="vcf", force=True)
        
        gz_path = abs_vcf_path + ".gz"
        tbi_path = gz_path + ".tbi"
        
        if os.path.exists(gz_path) and os.path.exists(tbi_path):
            if os.path.exists(abs_vcf_path):
                os.remove(abs_vcf_path)
            return gz_path
        else:
            return vcf_path
            
    except Exception as e:
        gz_path = abs_vcf_path + ".gz"
        tbi_path = gz_path + ".tbi"
        if os.path.exists(gz_path) and os.path.exists(tbi_path):
            if os.path.exists(abs_vcf_path):
                os.remove(abs_vcf_path)
            return gz_path
            
        print(f"\n[ERROR] Tabix indexing failed: {e}")
        return vcf_path

def split_vcf_by_tr(vcf_path, tr_db_path, output_dir, prefix):
    os.makedirs(output_dir, exist_ok=True)
    BIN_SIZE = 10000
    
    f_tr_tmp = os.path.join(output_dir, f"{prefix}.TR.vcf")
    f_nontr_tmp = os.path.join(output_dir, f"{prefix}.NonTR.vcf")
    
    print(f"[INFO] Splitting VCF: {vcf_path}")
    
    tr_index = load_tr_database(tr_db_path)
    
    tr_count = 0
    nontr_count = 0
    total_count = 0

    try:
        with smart_open(vcf_path, 'rt') as fin, \
             open(f_tr_tmp, 'w') as out_tr, \
             open(f_nontr_tmp, 'w') as out_nontr:
            
            for line in fin:
                if line.startswith("#"):
                    out_tr.write(line)
                    out_nontr.write(line)
                    continue
                
                parts = line.strip().split("\t")
                if len(parts) < 4: continue
                
                chrom = normalize_chrom(parts[0])
                try:
                    pos = int(parts[1])
                    ref = parts[3]
                    var_start = pos
                    var_end = pos + len(ref) - 1
                except ValueError:
                    out_nontr.write(line)
                    continue

                bin_id = pos // BIN_SIZE
                candidates = tr_index[chrom].get(bin_id, [])
                is_pure_tr = False
                
                for tr_start, tr_end in candidates:
                    if tr_start <= var_start and var_end <= tr_end:
                        is_pure_tr = True
                        break
                
                if is_pure_tr:
                    out_tr.write(line)
                    tr_count += 1
                else:
                    out_nontr.write(line)
                    nontr_count += 1
                
                total_count += 1
                if total_count % 10000 == 0:
                    print(f"       Processed {total_count} variants...", end="\r")
        
        print(f"\n[INFO] Writing complete. Starting compression and indexing...")

        final_tr_path = compress_and_index(f_tr_tmp)
        final_nontr_path = compress_and_index(f_nontr_tmp)

    except FileNotFoundError:
        print(f"[ERROR] Input VCF not found: {vcf_path}")
        sys.exit(1)
    except Exception as e:
        print(f"[ERROR] An error occurred: {e}")
        sys.exit(1)

    tr_pct = (tr_count / total_count * 100) if total_count > 0 else 0
    nontr_pct = (nontr_count / total_count * 100) if total_count > 0 else 0

    print(f"\n[DONE] Finished.")
    print(f"  Total Variants: {total_count}")
    print(f"  TR Regions:     {tr_count} ({tr_pct:.2f}%) -> {final_tr_path}")
    print(f"  Non-TR Regions: {nontr_count} ({nontr_pct:.2f}%) -> {final_nontr_path}")