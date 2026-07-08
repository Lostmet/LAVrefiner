import os
import sys
import pickle
import gzip
from collections import defaultdict

try:
    import Levenshtein
except ImportError:
    pass

def _nested_dict_factory():
    return defaultdict(list)

class VcfIndexManager:
    def __init__(self, vcf_path):
        self.vcf_path = vcf_path
        self.index_path = vcf_path + ".memory_index.v2.pkl"
        
        self.tr_data = defaultdict(_nested_dict_factory) 
        self.BIN_SIZE = 10000

        if self._check_index_valid():
            self.load_index()
        else:
            self.build_index()

    def _check_index_valid(self):
        if not os.path.exists(self.index_path):
            return False
        if os.path.getmtime(self.vcf_path) > os.path.getmtime(self.index_path):
            return False
        return True

    def _smart_open(self, path):
        if path.endswith(".gz"):
            return gzip.open(path, "rt", encoding="utf-8", errors="replace")
        return open(path, "r", encoding="utf-8", errors="replace")

    def build_index(self):
        print(f"[TR-Align] Building in-memory index for {self.vcf_path}...")
        print("           (This may take 1-2 mins for the first time, but will speed up alignment significantly)")
        
        count = 0
        self.tr_data = defaultdict(_nested_dict_factory)
        
        try:
            with self._smart_open(self.vcf_path) as f:
                for line in f:
                    if line.startswith("#"): continue
                    
                    parts = line.strip().split("\t")
                    if len(parts) < 8: continue
                    
                    chrom = parts[0]
                    if chrom.startswith("chr"): chrom = chrom[3:]
                    
                    try:
                        start = int(parts[1])
                        info_str = parts[7]
                        end = start + len(parts[3]) - 1 # Fallback
                        
                        motifs = []
                        for field in info_str.split(";"):
                            if field.startswith("END="):
                                end = int(field[4:])
                            elif field.startswith("RU="):
                                motifs = field[3:].split(",")
                        
                        if not motifs: continue

                        record = (start, end, motifs)
                        
                        start_bin = start // self.BIN_SIZE
                        end_bin = end // self.BIN_SIZE
                        
                        for b in range(start_bin, end_bin + 1):
                            self.tr_data[chrom][b].append(record)
                        
                        count += 1
                        if count % 100000 == 0:
                            print(f"           Indexed {count} records...", end="\r")
                            
                    except ValueError:
                        continue
                        
            print(f"\n[TR-Align] Index built. {count} records loaded into memory.")
            
            with open(self.index_path, 'wb') as f:
                pickle.dump(self.tr_data, f)
            print(f"[TR-Align] Index saved to {self.index_path}")
                
        except Exception as e:
            import traceback
            traceback.print_exc()
            print(f"\n[Error] Failed to build index: {e}")
            sys.exit(1)

    def load_index(self):
        print(f"[TR-Align] Loading cached index: {self.index_path}")
        try:
            with open(self.index_path, 'rb') as f:
                self.tr_data = pickle.load(f)
        except Exception as e:
            print(f"[Warning] Failed to load cache, rebuilding... ({e})")
            self.build_index()

    def search_region(self, chrom, start, end):
        if chrom.startswith("chr"):
            chrom_key = chrom[3:]
        else:
            chrom_key = chrom

        start_bin = start // self.BIN_SIZE
        end_bin = end // self.BIN_SIZE

        best_match = None
        max_overlap_len = 0
        
        candidates = []
        for b in range(start_bin, end_bin + 1):
            if b in self.tr_data[chrom_key]:
                candidates.extend(self.tr_data[chrom_key][b])

        if not candidates and chrom_key != chrom:
             for b in range(start_bin, end_bin + 1):
                if b in self.tr_data[chrom]:
                    candidates.extend(self.tr_data[chrom][b])

        if not candidates: return None

        for (db_start, db_end, motifs) in candidates:
            o_s = max(start, db_start)
            o_e = min(end, db_end)
            
            if o_s <= o_e:
                overlap_len = o_e - o_s + 1
                
                if overlap_len >= max_overlap_len:
                    max_overlap_len = overlap_len
                    best_match = (db_start, db_end, motifs)
        
        if best_match:
            return (None, best_match[0], best_match[1], best_match[2])
            
        return None

class MotifOptimizer:
    @staticmethod
    def get_seq_change(ref, alt):
        if len(ref) > 0 and len(alt) > 0 and ref[0] == alt[0]:
            if len(ref) == 1: return alt[1:]
            if len(alt) == 1: return ref[1:]
            return alt[len(ref):] if len(alt) > len(ref) else ref[len(alt):] 
        return alt

    @staticmethod
    def decompose_structure(sequence, raw_motifs):
        if 'Levenshtein' not in sys.modules:
            raise ImportError("Levenshtein library required for TR mode.")

        if not sequence: return [], 0
        
        N = len(sequence)
        UNK_PENALTY = 1
        
        dp = [float('inf')] * (N + 1)
        dp[0] = 0
        parent = [None] * (N + 1)
        
        motifs_data = []
        for i, m in enumerate(raw_motifs):
            m_len = len(m)
            threshold = max(2, int(m_len * 0.2))
            motifs_data.append((i, m, m_len, threshold))

        # --- DP ---
        for i in range(1, N + 1):
            # Option 1: UNK
            cost_unk = dp[i-1] + UNK_PENALTY
            if cost_unk < dp[i]:
                dp[i] = cost_unk
                parent[i] = (i-1, "UNK", None)
            
            # Option 2: Motif
            for m_id, m_seq, m_len, limit in motifs_data:
                min_sub = max(1, m_len - limit)
                max_sub = m_len + limit
                
                start_j = max(0, i - max_sub)
                end_j = i - min_sub
                if end_j < start_j: continue
                
                for j in range(start_j, end_j + 1):
                    if dp[j] == float('inf'): continue
                    sub = sequence[j:i]
                    
                    # Levenshtein distance check
                    dist = Levenshtein.distance(sub, m_seq)
                    if dist <= limit:
                        new_cost = dp[j] + dist
                        if new_cost < dp[i]:
                            dp[i] = new_cost
                            parent[i] = (j, "MOTIF", m_id)

        # --- Backtrack ---
        segments = [] 
        curr = N
        unk_buffer = []
        
        while curr > 0 and parent[curr] is not None:
            prev, desc, m_id = parent[curr]
            seq_slice = sequence[prev:curr]
            
            if desc == "UNK":
                unk_buffer.append(seq_slice)
            else:
                if unk_buffer:
                    full_unk = "".join(reversed(unk_buffer))
                    segments.append((full_unk, -1))
                    unk_buffer = []
                segments.append((seq_slice, m_id))
            curr = prev
            
        if unk_buffer:
            full_unk = "".join(reversed(unk_buffer))
            segments.append((full_unk, -1))
            
        return list(reversed(segments)), dp[N]