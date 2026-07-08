import os
import gzip
import pysam
from collections import defaultdict
from tqdm import tqdm
import concurrent.futures
from ..TR_align.motif_utils import VcfIndexManager, MotifOptimizer
from ..utils.logging_utils import get_logger

logger = get_logger(__name__)

def gt_to_dosage(sample_rec):
    gt = sample_rec.get('GT', (None, None))
    if gt is None or any(a is None for a in gt): return None
    return sum(1 for a in gt if a == 1)

def exact_to_gt(ds, d):
    if ds is None: return './.'
    if ds == 0:    return '0/0'
    if ds == d:    return '0/1'
    if d != 0 and ds == 2 * d: return '1/1'
    return '0/0'

def count_motifs(segments):
    counts = defaultdict(int)
    for _, mid in segments:
        if mid >= 0: counts[mid] += 1
    return counts

def read_vcf_meta_header(vcf_path):
    lines = []
    if str(vcf_path).endswith('.gz'):
        file_obj = gzip.open(vcf_path, 'rt', encoding='utf-8')
    else:
        file_obj = open(vcf_path, 'r', encoding='utf-8')
        
    with file_obj as f:
        for line in f:
            if line.startswith('##'): 
                lines.append(line.rstrip('\n'))
            else: 
                break
                
    return lines

def write_vcf_header(fout, meta_lines, new_info_lines, samples):
    for line in meta_lines: fout.write(line + '\n')
    for line in new_info_lines: fout.write(line + '\n')
    cols = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
    fout.write('\t'.join(cols + samples) + '\n')

NEW_INFO_LINES = [
    '##INFO=<ID=RL,Number=1,Type=Integer,Description="Repeat unit length (bp)">',
    '##INFO=<ID=CN_DELTA,Number=1,Type=Integer,Description="Net copy-number change at locus (exact, aggregated across all variants)">',
    '##INFO=<ID=MOTIF_ID,Number=1,Type=Integer,Description="Motif index in TRCompDB RU list for this locus">',
    '##INFO=<ID=MOTIF_SEQ,Number=1,Type=String,Description="Motif sequence">',
    '##INFO=<ID=MOTIF_DELTA,Number=1,Type=Integer,Description="Net motif copy change at locus (exact, aggregated across all variants)">',
    '##INFO=<ID=TR_LOCUS,Number=1,Type=String,Description="TRCompDB locus identifier chrom:start-end">',
]

def process_locus_task(locus_key, tr_start, tr_end, motifs, rl_locus, n_samp, slices):
    chrom = locus_key[0]
    locus_name = f"{chrom}:{tr_start}-{tr_end}"
    
    cn_bi_lines, cn_ex_lines = [], []
    md_bi_lines, md_ex_lines = [], []

    cn_exact = [0] * n_samp
    motif_exact = defaultdict(lambda: [0] * n_samp)

    for s in slices:
        ref_slice = s['ref']
        alt_slice = s['alt']
        dos = s['dos_arr']

        is_ins = len(alt_slice) > len(ref_slice)
        is_del = len(alt_slice) < len(ref_slice)
        delta_sign = +1 if is_ins else (-1 if is_del else 0)
        
        changed_seq = alt_slice if is_ins else ref_slice
        motif_delta = {}
        
        if changed_seq and motifs:
            segments, _ = MotifOptimizer.decompose_structure(changed_seq, motifs)
            counts = count_motifs(segments)
            motif_delta = {mid: cnt * delta_sign for mid, cnt in counts.items()}
        
        cn_delta = round((len(alt_slice) - len(ref_slice)) / rl_locus) if rl_locus > 0 else 0

        for i, d in enumerate(dos):
            if d is None: cn_exact[i] = None
            elif cn_exact[i] is not None: cn_exact[i] += d * cn_delta
        
        for mid, mdelta in motif_delta.items():
            for i, d in enumerate(dos):
                if d is None: motif_exact[mid][i] = None
                elif motif_exact[mid][i] is not None: motif_exact[mid][i] += d * mdelta

    cn_distinct = sorted(set(v for v in cn_exact if v is not None and v != 0))
    for d in cn_distinct:
        row_id = f"{locus_name}_copy{d}"
        gt_fields = [exact_to_gt(ds, d) for ds in cn_exact]
        cn_bi_lines.append('\t'.join([str(chrom), str(tr_start), row_id, 'N', '<TR>', '.', '.', f"RL={rl_locus};CN_DELTA={d};TR_LOCUS={locus_name}", 'GT'] + gt_fields) + '\n')

    if any(v is not None and v != 0 for v in cn_exact):
        repeat_unit = ','.join(motifs[:4]) if motifs else "."
        cn_ex_lines.append('\t'.join([str(chrom), str(tr_start), locus_name, repeat_unit, str(rl_locus)] + ['NA' if v is None else str(v) for v in cn_exact]) + '\n')

    for mid in sorted(motif_exact.keys()):
        mex = motif_exact[mid]
        if all(v is None or v == 0 for v in mex): continue
        motif_seq = motifs[mid] if mid < len(motifs) else f"motif{mid}"
        motif_id = f"{locus_name}_motif{mid}({motif_seq})"
        
        m_distinct = sorted(set(v for v in mex if v is not None and v != 0))
        for d in m_distinct:
            row_id = f"{motif_id}_delta{d}"
            gt_fields = [exact_to_gt(ds, d) for ds in mex]
            md_bi_lines.append('\t'.join([str(chrom), str(tr_start), row_id, 'N', '<TR>', '.', '.', f"RL={rl_locus};MOTIF_ID={mid};MOTIF_SEQ={motif_seq};MOTIF_DELTA={d};TR_LOCUS={locus_name}", 'GT'] + gt_fields) + '\n')
        md_ex_lines.append('\t'.join([str(chrom), str(tr_start), motif_id, motif_seq, str(len(motif_seq))] + ['NA' if v is None else str(v) for v in mex]) + '\n')

    return (cn_bi_lines, cn_ex_lines, md_bi_lines, md_ex_lines)

def process_locus_task_wrapper(kwargs):
    return process_locus_task(**kwargs)

def run_dosage(vcf_path, tr_db_path, out_dir, prefix, threads=1):
    os.makedirs(out_dir, exist_ok=True)
    tr_index = VcfIndexManager(tr_db_path)
    
    vcf = pysam.VariantFile(vcf_path, threads=threads)
    samples = list(vcf.header.samples)
    n_samp = len(samples)
    meta_lines = read_vcf_meta_header(vcf_path)
    EXACT_HEADER = ['#CHROM', 'POS', 'ID', 'REPEAT_UNIT', 'REPEAT_LENGTH'] + samples

    paths = {
        'cn_bi': os.path.join(out_dir, f"{prefix}.copy_number_biallelic.vcf"),
        'cn_ex': os.path.join(out_dir, f"{prefix}.copy_number_exact.tsv"),
        'md_bi': os.path.join(out_dir, f"{prefix}.motif_dosage_biallelic.vcf"),
        'md_ex': os.path.join(out_dir, f"{prefix}.motif_dosage_exact.tsv"),
    }

    fcn_bi = open(paths['cn_bi'], 'w')
    fcn_ex = open(paths['cn_ex'], 'w')
    fmd_bi = open(paths['md_bi'], 'w')
    fmd_ex = open(paths['md_ex'], 'w')

    try:
        write_vcf_header(fcn_bi, meta_lines, NEW_INFO_LINES, samples)
        write_vcf_header(fmd_bi, meta_lines, NEW_INFO_LINES, samples)
        fcn_ex.write('\t'.join(EXACT_HEADER) + '\n')
        fmd_ex.write('\t'.join(EXACT_HEADER) + '\n')

        active_loci = {}
        current_chrom = None
        processed_loci = 0
        batch = []
        BATCH_SIZE = 50 * threads

        executor = concurrent.futures.ProcessPoolExecutor(max_workers=threads) if threads > 1 else None

        def flush_batch():
            nonlocal batch, processed_loci
            if not batch: return
            if executor:
                results = list(executor.map(process_locus_task_wrapper, batch))
            else:
                results = [process_locus_task(**kwargs) for kwargs in batch]
            
            for cn_bi, cn_ex, md_bi, md_ex in results:
                fcn_bi.writelines(cn_bi)
                fcn_ex.writelines(cn_ex)
                fmd_bi.writelines(md_bi)
                fmd_ex.writelines(md_ex)
            processed_loci += len(batch)
            batch.clear()

        for rec in tqdm(vcf.fetch(), desc="Processing VCF variants", unit=" vars", dynamic_ncols=True):
            chrom, pos, ref, alt = rec.chrom, rec.pos, rec.ref, rec.alts[0] if rec.alts else None
            if alt is None: continue
            
            if current_chrom != chrom:
                batch.extend(active_loci.values())
                active_loci.clear()
                flush_batch()
                current_chrom = chrom
            else:
                to_flush = [k for k, v in active_loci.items() if v['tr_end'] < pos - 5000]
                for k in to_flush:
                    batch.append(active_loci.pop(k))
                if len(batch) >= BATCH_SIZE:
                    flush_batch()

            var_start = pos
            var_end = pos + len(ref) - 1
            
            overlapping_trs = []
            p = var_start
            while p <= var_end:
                result = tr_index.search_region(chrom, p, var_end)
                if not result: break
                _, t_start, t_end, motifs = result
                if t_start > var_end: break
                overlapping_trs.append((t_start, t_end, motifs))
                p = max(p + 1, t_end + 1)
            
            if not overlapping_trs: continue

            is_ins = len(alt) > len(ref)
            is_del = len(alt) < len(ref)
            dos_arr = [gt_to_dosage(rec.samples[s]) for s in samples]

            for t_start, t_end, motifs in overlapping_trs:
                if not motifs: continue
                locus_key = (chrom, t_start, t_end)
                
                slice_ref, slice_alt = "", ""
                
                if is_ins:
                    if t_start <= var_start <= t_end:
                        slice_alt = MotifOptimizer.get_seq_change(ref, alt) or alt[1:]
                elif is_del:
                    o_start = max(var_start + 1, t_start)
                    o_end = min(var_end, t_end)
                    if o_start <= o_end:
                        slice_ref = ref[o_start - var_start : o_end - var_start + 1]
                else:
                    if t_start <= var_start <= t_end:
                        slice_alt = MotifOptimizer.get_seq_change(ref, alt) or alt[1:]

                if slice_ref or slice_alt:
                    if locus_key not in active_loci:
                        rl_locus = len(motifs[0]) if motifs else 1
                        active_loci[locus_key] = {
                            'locus_key': locus_key, 'tr_start': t_start, 'tr_end': t_end,
                            'motifs': motifs, 'rl_locus': rl_locus, 'n_samp': n_samp, 'slices': []
                        }
                    active_loci[locus_key]['slices'].append({
                        'ref': slice_ref, 'alt': slice_alt, 'dos_arr': dos_arr
                    })

        batch.extend(active_loci.values())
        active_loci.clear()
        flush_batch()

    finally:
        fcn_bi.close(); fcn_ex.close()
        fmd_bi.close(); fmd_ex.close()
        vcf.close()
        if threads > 1 and executor: executor.shutdown()

    print()
    logger.info(f"Dosage calculation complete. Successfully processed {processed_loci} TR loci.")
    return paths