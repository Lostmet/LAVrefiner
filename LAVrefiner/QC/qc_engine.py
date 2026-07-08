import os
import pysam
import time
from pathlib import Path
from ..utils.logging_utils import get_logger

logger = get_logger(__name__)

def run_qc(vcf_path: str, 
           out_dir: str, 
           prefix: str, 
           tr_div_threshold: float = None, 
           exclude_type: str = None, 
           extract_type: str = None, 
           min_rl: int = None, 
           max_rl: int = None,
           tr_vcf_path: str = None,
           exclude_homopoly: bool = False): 
    """
    Filter rLAV VCF based on TR metrics (TR_DIV, TR_TYPE, RL).
    Supports 'TR' as a wildcard for both STR and VNTR.
    If tr_vcf_path is provided, it supplements missing TR_TYPE/RL info by searching the DB.
    """

    if not os.path.exists(vcf_path):
        raise FileNotFoundError(f"Input VCF not found: {vcf_path}")
    
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    out_name = f"{prefix}.QC.vcf"
    out_path = os.path.join(out_dir, out_name)

    tr_index = None
    if tr_vcf_path:
        if os.path.exists(tr_vcf_path):
            logger.info(f"Loading TR database for annotation: {tr_vcf_path}")
            from ..TR_align.motif_utils import VcfIndexManager
            tr_index = VcfIndexManager(tr_vcf_path)
        else:
            logger.warning(f"Provided --tr-vcf not found: {tr_vcf_path}. Skipping annotation.")

    vcf_in = pysam.VariantFile(vcf_path, "r")
    header = vcf_in.header

    required_info = [
        ("TR_TYPE", 1, "String", "TR Type: STR (<=6bp) or VNTR (>6bp)"),
        ("RL", 1, "Integer", "Repeat Unit Length"),
        ("TYPE", ".", "String", "Variant type"),
        ("SVTYPE", "1", "String", "Type of structural variant"),
        ("SVLEN", ".", "Integer", "Difference in length"),
        ("END", "1", "Integer", "End position"),
        ("MOTIF", "1", "String", "Repeat motif"),
        ("RU", "1", "String", "Repeat unit")
    ]

    for id_tag, num, v_type, desc in required_info:
        if id_tag not in header.info:
            header.info.add(id_tag, num, v_type, desc)

    if 'PASS' not in header.filters:
        header.filters.add('PASS', None, None, 'All filters passed')

    cmd_record = f'##LAVrefiner_QC_Command="TR_DIV={tr_div_threshold}, ..."'
    header.add_line(cmd_record)

    vcf_out = pysam.VariantFile(out_path, "w", header=header)

    stats = {
        "Total": 0,
        "Passed": 0,
        "Annotated_by_DB": 0,
        "Filtered_TR_DIV": 0,
        "Filtered_Exclude": 0,
        "Filtered_Extract": 0,
        "Filtered_MinRL": 0,
        "Filtered_MaxRL": 0,
        "Filtered_Homopoly": 0 
    }
    
    logger.info(f"Starting QC on {vcf_path}...")
    start_time = time.time()

    for record in vcf_in:
        stats["Total"] += 1

        if exclude_homopoly:
            ref = record.ref
            alt = record.alts[0] if record.alts else ""

            if not alt.startswith("<") and not alt.startswith("[") and not alt.endswith("]"):
                combined_seq = (str(ref) + str(alt)).upper()
                if len(set(combined_seq)) == 1:
                    stats["Filtered_Homopoly"] += 1
                    continue   
        # === Step 0: Supplement TR Info (Universal QC Logic) ===
        current_tr_type = record.info.get("TR_TYPE")
        current_rl = record.info.get("RL")
        
        if (current_tr_type is None or current_rl is None) and tr_index is not None:
            var_start = record.pos
            var_end = record.pos + len(record.ref) - 1
            
            if var_start == var_end: var_end += 1

            match = tr_index.search_region(record.chrom, var_start, var_end)
            
            if match:
                motifs = match[-1]
                
                if motifs:
                    rl_val = len(motifs[0])
                    type_val = "VNTR" if rl_val > 6 else "STR"

                    record.info["RL"] = rl_val
                    record.info["TR_TYPE"] = type_val

                    current_rl = rl_val
                    current_tr_type = type_val
                    
                    stats["Annotated_by_DB"] += 1

        # --- Filter 1: TR_DIV (Divergence Score) ---
        if tr_div_threshold is not None:
            div = record.info.get("TR_DIV")
            if div is not None and div >= tr_div_threshold:
                stats["Filtered_TR_DIV"] += 1
                continue

        # --- Filter 2 & 3: TR_TYPE (STR / VNTR / TR) ---
        # Exclude Logic
        if exclude_type is not None:
            should_exclude = False
            if exclude_type == 'TR':
                if current_tr_type is not None: should_exclude = True
            elif current_tr_type == exclude_type:
                should_exclude = True
            
            if should_exclude:
                stats["Filtered_Exclude"] += 1
                continue
        
        # Extract Logic
        if extract_type is not None:
            should_keep = False
            if extract_type == 'TR':
                if current_tr_type is not None: should_keep = True
            elif current_tr_type == extract_type:
                should_keep = True
            
            if not should_keep:
                stats["Filtered_Extract"] += 1
                continue

        # --- Filter 4: RL (Repeat Length) ---
        if min_rl is not None or max_rl is not None:
            if current_rl is not None:
                if min_rl is not None and current_rl < min_rl:
                    stats["Filtered_MinRL"] += 1
                    continue
                if max_rl is not None and current_rl > max_rl:
                    stats["Filtered_MaxRL"] += 1
                    continue
            else:
                pass

        # --- Passed ---
        vcf_out.write(record)
        stats["Passed"] += 1

    vcf_in.close()
    vcf_out.close()

    total = stats['Total']
    
    def _pct(count):
        if total == 0: return "0.00%"
        return f"{count/total*100:.2f}%"

    logger.info(f"QC Completed. Output: {out_path}")
    logger.info("-" * 60)
    logger.info(f"{'Metric':<30} {'Count':<10} {'Percentage':<10}")
    logger.info("-" * 60)
    logger.info(f"{'Total Variants':<30} {total:<10} {'100.00%':<10}")
    
    if tr_index:
        logger.info(f"{'Annotated (Missing info)':<30} {stats['Annotated_by_DB']:<10} {_pct(stats['Annotated_by_DB']):<10}")
        
    logger.info(f"{'Passed QC':<30} {stats['Passed']:<10} {_pct(stats['Passed']):<10}")
    logger.info("-" * 60)
    
    if tr_div_threshold: 
        logger.info(f"{'Filtered (High TR_DIV)':<30} {stats['Filtered_TR_DIV']:<10} {_pct(stats['Filtered_TR_DIV']):<10}")
    if exclude_homopoly:
        logger.info(f"{'Filtered (Homopoly)':<30} {stats['Filtered_Homopoly']:<10} {_pct(stats['Filtered_Homopoly']):<10}")
    if exclude_type:     
        logger.info(f"{'Filtered (Exclude)':<30} {stats['Filtered_Exclude']:<10} {_pct(stats['Filtered_Exclude']):<10}")
    if extract_type:     
        logger.info(f"{'Filtered (Extract)':<30} {stats['Filtered_Extract']:<10} {_pct(stats['Filtered_Extract']):<10}")
    if min_rl:           
        logger.info(f"{'Filtered (Short RL)':<30} {stats['Filtered_MinRL']:<10} {_pct(stats['Filtered_MinRL']):<10}")
    if max_rl:           
        logger.info(f"{'Filtered (Long RL)':<30} {stats['Filtered_MaxRL']:<10} {_pct(stats['Filtered_MaxRL']):<10}")
    logger.info("-" * 60)

    return out_path