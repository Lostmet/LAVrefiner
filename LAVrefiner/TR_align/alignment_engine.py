import numpy as np

def align_lavs_strict_left(lav_dict):
    seqs = list(lav_dict.values())
    names = list(lav_dict.keys())

    indexed_seqs = sorted(enumerate(seqs), key=lambda x: len(x[1]), reverse=True)
    anchor_idx, anchor_seq = indexed_seqs[0]

    msa = [list(anchor_seq)]
    msa_row_indices = [anchor_idx]
    
    GAP_VAL = -1

    for original_idx, query_seq in indexed_seqs[1:]:
        current_anchor = msa[0]
        n, m = len(current_anchor), len(query_seq)
        
        # DP Initialization
        score = [[0] * (m+1) for _ in range(n+1)]
        ptr = [[0] * (m+1) for _ in range(n+1)]
        
        GAP_SCORE = -1
        MATCH_SCORE = 10
        MISMATCH_SCORE = -999999
        
        for i in range(n+1):
            score[i][0] = i * GAP_SCORE
            ptr[i][0] = 1 # Up (Ref Char, Query Gap)
        for j in range(m+1):
            score[0][j] = j * GAP_SCORE
            ptr[0][j] = 2 # Left (Ref Gap, Query Char)
            
        # DP Fill
        for i in range(1, n+1):
            for j in range(1, m+1):
                ref_item = current_anchor[i-1]
                query_item = query_seq[j-1]

                if ref_item == GAP_VAL or query_item == GAP_VAL:
                    is_match = False
                else:
                    is_match = (ref_item == query_item)
                
                s_diag = score[i-1][j-1] + (MATCH_SCORE if is_match else MISMATCH_SCORE)
                s_up = score[i-1][j] + GAP_SCORE
                s_left = score[i][j-1] + GAP_SCORE
                
                best = max(s_diag, s_up, s_left)
                score[i][j] = best

                if best == s_up:
                    ptr[i][j] = 1
                elif best == s_diag:
                    ptr[i][j] = 0
                else:
                    ptr[i][j] = 2
                    
        # Traceback
        aligned_anchor = []
        aligned_query = []
        ops = []
        
        i, j = n, m
        while i > 0 or j > 0:
            p = ptr[i][j]
            if i > 0 and j > 0 and p == 0: # Match
                aligned_anchor.append(current_anchor[i-1])
                aligned_query.append(query_seq[j-1])
                ops.append('M')
                i -= 1; j -= 1
            elif i > 0 and (j == 0 or p == 1): # Up (Query Gap)
                aligned_anchor.append(current_anchor[i-1])
                aligned_query.append(GAP_VAL)
                ops.append('U')
                i -= 1
            else: # Left (Insertion in Query)
                aligned_anchor.append(GAP_VAL)
                aligned_query.append(query_seq[j-1])
                ops.append('L')
                j -= 1
        
        ops.reverse()
        aligned_query.reverse()
        
        # Update MSA
        new_msa = [[] for _ in range(len(msa))]
        old_col_idx = 0
        
        for op in ops:
            if op == 'M' or op == 'U':
                for r in range(len(msa)):
                    new_msa[r].append(msa[r][old_col_idx])
                old_col_idx += 1
            else: # L -> Insert Gap into all existing rows
                for r in range(len(msa)):
                    new_msa[r].append(GAP_VAL)
        
        new_msa.append(aligned_query)
        msa = new_msa
        msa_row_indices.append(original_idx)
        
    # Reconstruct result dict
    result = {}
    idx_to_row = {orig_idx: msa[r_idx] for r_idx, orig_idx in enumerate(msa_row_indices)}
    
    for i, name in enumerate(names):
        result[name] = idx_to_row[i]
        
    return result