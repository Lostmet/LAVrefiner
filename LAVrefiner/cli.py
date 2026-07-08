import os
import time
import click
import shutil
import logging
import datetime
from pathlib import Path
from .config import Config
from .utils.logging_utils import *
from .genotype.genotype_mapper import *
from .genotype.matrix_generator import *
from .rLAV.rLAV_generator import rLAV_generator
from .alignment.alignment import run_alignments
from .utils.file_utils import check_vcf_vs_fasta
from .variant_processing.fasta_generator import generate_fasta_sequences
from .variant_processing.vcf_parser import process_variants, filter_vcf, make_oLAV
from .TR_align.manager import TRAlignmentManager
from .TR_align.motif_utils import VcfIndexManager
from .utils.tr_splitter import split_vcf_by_tr
from .QC.qc_engine import run_qc
from .QC.dosage_calculator import run_dosage
logger = get_logger(__name__)


def common_options(require_ref: bool = True):
    """
    Decorator factory: add common CLI options.
    If require_ref=False, do not add --ref option.
    """
    def decorator(func):
        options = [
            click.option('--vcf', required=True, help='Input VCF file'),
        ]

        if require_ref:
            options.append(
                click.option('--ref', required=True, help='Reference FASTA file')
            )

        options.extend([
            click.option('--out', required=True, help='Output directory'),
            click.option('--threads', default=10, type=int, help='Number of threads'),
            click.option(
                '--write-matrix',
                default='NO',
                type=str,
                callback=lambda ctx, param, value: value.upper() == "YES",
                help='Write X & T matrices (YES/NO)'
            ),
            click.option('--tr-vcf', default=None, help='Path to TRCompDB VCF for TR mode (Optional)'),
        ])

        for opt in options:
            func = opt(func)
        return func

    return decorator


def check_tools(*tools):
    """
    Check if required tools are available in the user's PATH.
    Raise an error if any tool is missing.
    """
    missing_tools = [tool for tool in tools if shutil.which(tool) is None]
    if missing_tools:
        raise RuntimeError(
            f"The following tools are missing from your environment: {', '.join(missing_tools)}. "
            f"Please ensure they are installed and accessible from your PATH."
        )


@click.group()
@click.version_option("0.3", prog_name="LAVrefiner", message="""
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
LAVrefiner - A Python tool for generating refined Length-Altering Variants (rLAVs)  
Version 0.3 - Contact: fenglostmet@tju.edu.cn
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
""")
def cli():
    """A Python tool for generating refined Length-Altering Variants (rLAVs)."""
    pass


@cli.command("stats")
@common_options(require_ref=False)
def stats(vcf: str, out: str, threads: int, write_matrix: bool, tr_vcf):
    """Only get overlapping information of input VCF."""
    setup_logging(Path(out) / "stats.log")
    try:
        logger.info(f"{'Command:':<5}{get_clean_command()}")
        start_time = time.time()
        logger.info(f"Start time: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        log_thread_info(threads)

        log_step("Processing VCF")

        config = Config(
            vcf_file=vcf,
            ref_fasta=None,
            output_dir=out,
            threads=threads
        )

        grouped_variants_list, single_lav_count, \
        multi_bp, percentage_lav_overlapped, \
        single_group, inv_count, variant_count, snp_groups = process_variants(config)

        stats_dict = {
            "Total variants": f"{variant_count:,}",
            "INV count": f"{inv_count:,}",
            "SNP count": f"{len(snp_groups):,}",
            "nLAV (non-overlap) count": f"{single_lav_count:,}",
            "Overlapping LAVs": f"{variant_count - single_lav_count:,}",
            "Overlapping LAVs percentage": f"{percentage_lav_overlapped:.2f}%",
            "Total variant groups": f"{len(grouped_variants_list):,}",
        }

        # --- [Enhanced] TR Analysis ---
        if tr_vcf:
            logger.info(f"Analyzing TR variants using database: {tr_vcf}")
            if os.path.exists(tr_vcf):
                tr_index = VcfIndexManager(tr_vcf)
                
                # Helper function to count TRs in a list of groups
                def get_tr_counts(groups):
                    str_c, vntr_c, total_c = 0, 0, 0
                    for grp in groups:
                        for var in grp.variants:
                            total_c += 1
                            match = tr_index.search_region(var.chrom, var.start, var.end)
                            if match:
                                _, _, _, motifs = match
                                if motifs:
                                    if len(motifs[0]) > 6:
                                        vntr_c += 1
                                    else:
                                        str_c += 1
                    return str_c, vntr_c, total_c

                # 1. Calculate for Overlapping LAVs
                ov_str, ov_vntr, ov_total = get_tr_counts(grouped_variants_list)
                ov_tr = ov_str + ov_vntr

                # 2. Calculate for Non-overlapping LAVs (nLAVs)
                nov_str, nov_vntr, nov_total = get_tr_counts(single_group)
                nov_tr = nov_str + nov_vntr

                # 3. Calculate for SNPs (to complete the Total count)
                snp_str, snp_vntr, snp_total = get_tr_counts(snp_groups)
                snp_tr = snp_str + snp_vntr

                # 4. Grand Totals
                total_str = ov_str + nov_str + snp_str
                total_vntr = ov_vntr + nov_vntr + snp_vntr
                total_tr = total_str + total_vntr
                grand_total_vars = ov_total + nov_total + snp_total

                # Helper for formatted string with percentage
                def fmt(cnt, base, label):
                    pct = (cnt / base * 100) if base > 0 else 0.0
                    return f"{cnt:,} ({pct:.2f}% of {label})"

                stats_dict["----------------"] = "----------------"
                stats_dict["TR Analysis DB"] = os.path.basename(tr_vcf)
                
                # --- Total Statistics ---
                stats_dict["Total TR (STR+VNTR)"] = fmt(total_tr, grand_total_vars, "Total Variants")
                stats_dict["  - Total STR"] = fmt(total_str, grand_total_vars, "Total Variants")
                stats_dict["  - Total VNTR"] = fmt(total_vntr, grand_total_vars, "Total Variants")

                # --- Overlapping Statistics ---
                stats_dict["---- Overlapping ----"] = "----------------"
                stats_dict["Overlapping TR"] = fmt(ov_tr, ov_total, "Overlapping LAVs")
                stats_dict["  - Overlap STR"] = fmt(ov_str, ov_total, "Overlapping LAVs")
                stats_dict["  - Overlap VNTR"] = fmt(ov_vntr, ov_total, "Overlapping LAVs")

                # --- Non-Overlapping Statistics ---
                stats_dict["---- Non-Overlapping ----"] = "----------------"
                stats_dict["nLAV TR"] = fmt(nov_tr, nov_total, "nLAVs")
                stats_dict["  - nLAV STR"] = fmt(nov_str, nov_total, "nLAVs")
                stats_dict["  - nLAV VNTR"] = fmt(nov_vntr, nov_total, "nLAVs")
                
            else:
                logger.warning(f"TR VCF file not found: {tr_vcf}")
        # ---------------------------------

        end_time = time.time()
        runtime = end_time - start_time
        log_step("Summary")
        log_summary_block(
            cmd=get_clean_command(),
            start=start_time,
            duration=runtime,
            stats=stats_dict,
        )
        log_all_warnings_and_errors()

    except Exception as e:
        logger.error(f"Error in stats: {str(e)}")
        logging.error(f"Error in stats: {str(e)}")
        raise click.Abort()


@cli.command("align")
@common_options(require_ref=True)
def align_time(vcf: str, ref: str, out: str, threads: int, write_matrix: bool, tr_vcf: str = None):
    """Do alignment of overlapping LAVs"""
    import json
    setup_logging(Path(out) / "align.log")
    start_time = time.time()
    
    logger.info(f"{'Command:':<5}{get_clean_command()}")
    logger.info(f"Start time: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log_thread_info(threads)
    
    log_step("Step 1 Processing VCF and generating FASTA")
    config = Config(vcf_file=vcf, ref_fasta=ref, output_dir=out, threads=threads, tr_vcf=tr_vcf)
    
    grouped_variants_list, _, _, _, \
    _, _, _, _ = process_variants(config)

    grouped_variants_dict = {}
    for group in grouped_variants_list:
        grouped_variants_dict.setdefault(group.chrom, []).append(group)
    total_groups = sum(len(groups) for groups in grouped_variants_dict.values())
    logger.info(f"Overlapping variants grouped: {total_groups:,} groups")

    fasta_path, has_insertion_dict, poly_ins_list = generate_fasta_sequences(config, grouped_variants_dict, total_groups)
    logger.info(f"FASTA file created: {fasta_path}")

    log_step("Step 2 Running alignments")
    if tr_vcf:
        logger.info(f"TR Mode Enabled within Alignment Engine. Using DB: {tr_vcf}")
        
    alignments_config = Config(grouped_variants_file=fasta_path, ref_fasta=ref, output_dir=out, threads=threads)
    
    alignment_results = run_alignments(alignments_config, fasta_path, has_insertion_dict, poly_ins_list, tr_vcf=tr_vcf)
    
    tr_metadata = {}
    if alignment_results:
        for res in alignment_results:
            if res:
                entry = {}
                if res.tr_div is not None:
                    entry["TR_DIV"] = res.tr_div
                if res.tr_type is not None:
                    entry["TR_TYPE"] = res.tr_type
                if res.rl is not None:
                    entry["RL"] = res.rl
                
                if entry:
                    tr_metadata[res.group_id] = entry
                    
    if tr_metadata:
        metadata_path = Path(out) / "tr_metadata.json"
        with open(metadata_path, 'w') as f:
            json.dump(tr_metadata, f)
        logger.info(f"TR metadata cached to {metadata_path}")
    
    logger.info("Alignments completed. Results saved in 'alignment_results' directory.")    
    end_time = time.time()
    total_time = end_time - start_time
    hours, remainder = divmod(total_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    log_step(f"Total runtime: {int(hours)}:{int(minutes):02d}:{int(seconds):02d}")


@cli.command("process-vcf")
@common_options(require_ref=False)
@click.option('--ref', default=None, help='Reference FASTA file (optional, only used to check VCF vs FASTA compatibility)')
def process_vcf(vcf: str, out: str, threads: int, write_matrix: bool, tr_vcf: str = None, ref: str = None):
    """Process overlapping variants, group them, and generate nLAVs, oLAVs. 
    This step produces 'variants_pre_aligned.fasta', 'oLAV.vcf' and 'nLAV.vcf'."""
    setup_logging(Path(out) / "process_vcf.log")
    
    if ref:
        check_vcf_vs_fasta(vcf, ref)
        
    try:
        logger.info(f"{'Command:':<5}{get_clean_command()}")
        start_time = time.time()
        logger.info(f"Start time: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        log_thread_info(threads)

        log_step("Processing VCF")
        config = Config(vcf_file=vcf, ref_fasta=None, output_dir=out, threads=threads, tr_vcf=tr_vcf)
        grouped_variants_list, single_lav_count, \
        multi_bp, percentage_lav_overlapped, \
        single_group, inv_count, variant_count, snp_groups = process_variants(config)

        est_time = multi_bp / 3000
        hours, remainder = divmod(est_time, 3600)  
        minutes, _ = divmod(remainder, 60)  
        click.echo(f"Estimated runtime: {int(hours)}:{int(minutes):02d}:00")
        click.echo(click.style("Note: Actual runtime depends on CPU performance.", fg="yellow"))

        click.echo("Processing non-overlapping LAVs (nLAVs) and generating nLAV.vcf...")
        _ = filter_vcf(config, single_group, snp_groups)

        click.echo("Processing overlapping LAVs (oLAVs) and generating oLAV.vcf...")
        make_oLAV(config, grouped_variants_list) # make oLAV.vcf
        
        stats = {
                    "Total variants": f"{variant_count:,}",
                    "INV count": f"{inv_count:,}",
                    "nLAV count": f"{single_lav_count:,}",
                    "Excluded SNP count": f"{len(snp_groups)}",
                    "Overlapping LAVs": f"{variant_count - single_lav_count:,}",
                    "Overlapping LAVs percentage": f"{percentage_lav_overlapped:.2f}%",
                    "Total variant groups": f"{len(grouped_variants_list):,}",
                }
  
        end_time = time.time()
        runtime = end_time - start_time
        log_step("Summary")
        log_summary_block(
            cmd=get_clean_command(),
            start=start_time,
            duration=runtime,
            stats=stats)
        log_all_warnings_and_errors()

    except Exception as e:
        logger.info(f"Error in process-vcf: {str(e)}", err=True)
        logging.error(f"Error in process-vcf: {str(e)}")
        raise click.Abort()


@cli.command("run-all")
@common_options(require_ref=True)
def run_all(vcf: str, ref: str, out: str, threads: int, write_matrix: bool, tr_vcf: str = None):
    """Execute the complete variant processing pipeline."""
    setup_logging(Path(out) / "run-all.log")
    check_tools("mafft")
    check_vcf_vs_fasta(vcf, ref)
    writedown = write_matrix

    try:
        logger.info(f"{'Command:':<5}{get_clean_command()}")
        start_time = time.time()
        logger.info(f"Start time: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        log_thread_info(threads)

        log_step("Step 1 Processing VCF")
        config = Config(vcf_file=vcf, ref_fasta=ref, output_dir=out, threads=threads, tr_vcf=tr_vcf)
        
        grouped_variants_list, single_lav_count, \
        multi_bp, percentage_lav_overlapped, \
        single_group, inv_count, variant_count, snp_groups = process_variants(config)

        total_groups = len(grouped_variants_list)

        est_time = multi_bp / 3000
        hours, remainder = divmod(est_time, 3600)  
        minutes, _ = divmod(remainder, 60)  
        click.echo(f"Estimated runtime: {int(hours)}:{int(minutes):02d}:00")
        click.echo(click.style("Note: Actual runtime depends on CPU performance.", fg="yellow"))

        click.echo("Processing non-overlapping LAVs (nLAVs) and generating nLAV.vcf...")
        nLAV_name = filter_vcf(config, single_group, snp_groups)

        click.echo("Processing overlapping LAVs (oLAVs) and generating oLAV.vcf...")
        make_oLAV(config, grouped_variants_list) 
        
        grouped_variants_dict = {}
        for group in grouped_variants_list:
            grouped_variants_dict.setdefault(group.chrom, []).append(group)
        
        logger.info(f"Overlapping variants grouped: {total_groups:,} groups")

        has_insertion_dict = {}
        poly_ins_list = []

        log_step("Step 1b Generating FASTA")
        fasta_path, has_insertion_dict, poly_ins_list = generate_fasta_sequences(config, grouped_variants_dict, total_groups)
        logger.info(f"FASTA file created: {fasta_path}")

        log_step("Step 2 Running Alignments (Hybrid Mode)")
        if tr_vcf:
            logger.info(f"TR Mode Enabled within Alignment Engine. Using DB: {tr_vcf}")
        
        alignments_config = Config(grouped_variants_file=fasta_path, ref_fasta=ref, output_dir=out, threads=threads)
        alignment_results = run_alignments(alignments_config, fasta_path, has_insertion_dict, poly_ins_list, tr_vcf=tr_vcf)
        
        tr_metadata = {}
        if alignment_results:
            for res in alignment_results:
                if res:
                    entry = {}
                    if res.tr_div is not None:
                        entry["TR_DIV"] = res.tr_div
                    if res.tr_type is not None:
                        entry["TR_TYPE"] = res.tr_type
                    if res.rl is not None:
                        entry["RL"] = res.rl
                    
                    if entry:
                        tr_metadata[res.group_id] = entry
        
        logger.info("Alignments completed.")

        log_step("Step 3 Generating rLAVs")
        matrix_dir = Path(out) / "matrix_results"
        
        rLAV_generator(out, has_insertion_dict, threads, tr_metadata=tr_metadata)
        logger.info("rLAV matrix generation completed.")

        log_step("Step 4 Genotype Mapping and VCF Generation")
        
        rLAV_meta_csv = os.path.join(out, "rLAV_meta.csv")
        rLAV_vcf = os.path.join(out, "rLAV.vcf")

        if writedown:
            sample_names, _ = process_vcf_to_x_matrix(out, matrix_dir, writedown=True)
            compute_t_matrix(matrix_dir)
            click.echo("D, X, and T matrices successfully generated.")
            rLAV_count, sample_names2, gt_buffer = load_encoded_gt_from_matrix_dir(out)
            logging.info(f"rLAV count: {rLAV_count:,}")
            vcf_generate_from_gt(rLAV_meta_csv, rLAV_vcf, sample_names2, gt_buffer)
        else:
            click.echo("Generating rLAV.vcf in streamed mode...")
            rLAV_count, sample_names = build_rlav_vcf_streamed_nobuf(
                vcf_dir=out,
                output_dir=matrix_dir,
                meta_csv_file=rLAV_meta_csv,
                output_vcf_file=rLAV_vcf,
                strict_meta_group_check=True,
                SAMPLE_CHUNK=4096 
            )
            logging.info(f"rLAV count: {rLAV_count:,}")

        logging.info(f"{rLAV_vcf} successfully generated!")
        click.echo("Pipeline execution completed!")

        stats = {
            "Total variants": f"{variant_count:,}",
            "INV count": f"{inv_count:,}",
            "nLAV count": f"{single_lav_count:,}",
            "Excluded SNP count": f"{len(snp_groups)}",
            "Overlapping LAVs": f"{variant_count - single_lav_count:,}",
            "Overlapping LAVs percentage": f"{percentage_lav_overlapped:.2f}%",
            "Total variant groups": f"{total_groups:,}",
            "Final rLAV count": f"{rLAV_count:,}",
            "Mode": "TR-Aware" if tr_vcf else "Standard"
        }  
        end_time = time.time()
        runtime = end_time - start_time
        log_step("Summary")
        log_summary_block(
            cmd=get_clean_command(),
            start=start_time,
            duration=runtime,
            stats=stats)
        log_all_warnings_and_errors()

    except Exception as e:
        logger.error(f"Critical error in run-all: {str(e)}", exc_info=True)
        raise click.Abort()


@cli.command("make-rlav")
@common_options(require_ref=True)
def make_rLAV(vcf: str, ref: str, out: str, threads: int, write_matrix: bool, tr_vcf: str = None):
    '''If you already have alignment results and oLAV.vcf, you can make rLAVs by this tag. 
    Ensure the output path remains unchanged.'''
    import json
    writedown = write_matrix
    setup_logging(Path(out) / "make-rlav.log")
    start_time = time.time()
    log_step("Step 1 Extracting LAV information")
    
    config = Config(vcf_file=vcf, ref_fasta=ref, output_dir=out, threads=threads, tr_vcf=tr_vcf)

    grouped_variants_list, single_lav_count, \
    _, percentage_lav_overlapped, \
    _, inv_count, variant_count, snp_groups = process_variants(config)

    log_step("Step 2 Converting to VCF format and generating matrices")

    grouped_variants_dict = {}
    for group in grouped_variants_list:
        grouped_variants_dict.setdefault(group.chrom, []).append(group)
    total_groups = sum(len(groups) for groups in grouped_variants_dict.values())
    logger.info(f"Overlapping variants grouped: {total_groups:,} groups")

    fasta_path, has_insertion_dict, poly_ins_list = generate_fasta_sequences(config, grouped_variants_dict, total_groups)

    matrix_dir = Path(out) / "matrix_results"
    sample_names = sample_name_contract(vcf)
    
    tr_metadata = None
    metadata_path = Path(out) / "tr_metadata.json"
    if metadata_path.exists():
        with open(metadata_path, 'r') as f:
            tr_metadata = json.load(f)
        logger.info("Successfully loaded cached TR metadata.")
    
    click.echo("1. Generating rLAV_meta.csv with metadata (positions, reference, and alternate alleles)...")
    
    rLAV_generator(out, has_insertion_dict, threads, tr_metadata=tr_metadata)

    rLAV_meta_csv = os.path.join(out, "rLAV_meta.csv")
    rLAV_vcf = os.path.join(out, "rLAV.vcf")

    if writedown:
        sample_names, _ = process_vcf_to_x_matrix(out, matrix_dir, writedown=True)
        compute_t_matrix(matrix_dir)
        click.echo("D, X, and T matrices successfully generated in 'matrix_results' directory.")
        rLAV_count, sample_names2, gt_buffer = load_encoded_gt_from_matrix_dir(out)
        logging.info(f"rLAV count: {rLAV_count:,}")
        vcf_generate_from_gt(rLAV_meta_csv, rLAV_vcf, sample_names2, gt_buffer)
    else:
        click.echo("Generating rLAV.vcf in streamed mode (no X/T)...")
        rLAV_count, sample_names = build_rlav_vcf_streamed_nobuf(
            vcf_dir=out,
            output_dir=matrix_dir,
            meta_csv_file=rLAV_meta_csv,
            output_vcf_file=rLAV_vcf,
            strict_meta_group_check=True,
            SAMPLE_CHUNK=4096
        )
        logging.info(f"rLAV count: {rLAV_count:,}")

    logger.info(f"{rLAV_vcf} successfully generated!")
    click.echo("Pipeline execution completed!")
    click.echo("All steps completed successfully!")

    stats = {
        "Total variants": f"{variant_count:,}",
        "INV count": f"{inv_count:,}",
        "nLAV count": f"{single_lav_count:,}",
        "Excluded SNP count": f"{len(snp_groups)}",
        "Overlapping LAVs": f"{variant_count - single_lav_count:,}",
        "Overlapping LAVs percentage": f"{percentage_lav_overlapped:.2f}%",
        "Total variant groups": f"{total_groups:,}",
        "Final rLAV count": f"{rLAV_count:,}"
    }

    end_time = time.time()
    runtime = end_time - start_time
    log_step("Summary")
    log_summary_block(
        cmd=get_clean_command(),
        start=start_time,
        duration=runtime,
        stats=stats)
    log_all_warnings_and_errors()



@cli.command("QC")
@click.option('--vcf', required=True, help='Input VCF file (supports .vcf.gz)')
@click.option('--out', required=True, help='Output directory')
@click.option('--prefix', required=True, help='Output prefix for the filtered VCF')
@click.option('--tr-div', type=float, default=None, help='Exclude variants with TR_DIV >= threshold (e.g. 0.1)')
@click.option('--exclude', type=click.Choice(['STR', 'VNTR', 'TR']), default=None, help='Exclude variants of a specific TR type. TR = STR + VNTR.')
@click.option('--extract', type=click.Choice(['STR', 'VNTR', 'TR']), default=None, help='Extract ONLY variants of a specific TR type. TR = STR + VNTR.')
@click.option('--min-rl', type=int, default=None, help='Minimum Repeat Length (RL) to keep')
@click.option('--max-rl', type=int, default=None, help='Maximum Repeat Length (RL) to keep')
@click.option('--tr-vcf', default=None, help='Optional: Path to TRCompDB VCF. Used to annotate/filter variants that lack TR_TYPE/RL info.')
@click.option('--exclude-homopoly', is_flag=True, default=False, help='Exclude homopolymer variants (where REF and ALT consist of the same single nucleotide)')
def qc_command(vcf, out, prefix, tr_div, exclude, extract, min_rl, max_rl, tr_vcf, exclude_homopoly): 
    """
    Quality Control and Filtering for rLAV VCF.
    ...
    """
    setup_logging(Path(out) / "qc.log")
    
    if exclude and extract and exclude == extract:
        logger.error("Cannot exclude and extract the same TR type simultaneously.")
        raise click.Abort()

    try:
        logger.info(f"{'Command:':<5}{get_clean_command()}")
        
        run_qc(
            vcf_path=vcf,
            out_dir=out,
            prefix=prefix,
            tr_div_threshold=tr_div,
            exclude_type=exclude,
            extract_type=extract,
            min_rl=min_rl,
            max_rl=max_rl,
            tr_vcf_path=tr_vcf,
            exclude_homopoly=exclude_homopoly
        )
        
    except Exception as e:
        logger.error(f"Error in QC: {str(e)}", exc_info=True)
        raise click.Abort()

 
@cli.command("dosage")
@click.option('--vcf', required=True, help='Input VCF file (e.g., oLAV.vcf or rLAV.vcf)')
@click.option('--tr-vcf', required=True, help='Path to TRCompDB VCF index')
@click.option('--out', required=True, help='Output directory')
@click.option('--prefix', required=True, help='Prefix for the output dosage files')
@click.option('--threads', '-t', default=1, type=int, help='Number of threads to use') 
def dosage_command(vcf, tr_vcf, out, prefix, threads):
    """
    Calculate motif dosage and copy number variations for TR loci.
    Outputs exact TSV matrices and bi-allelic VCFs suitable for GWAS/eQTL.
    """
    setup_logging(Path(out) / f"{prefix}_dosage.log")
    
    try:
        logger.info(f"{'Command:':<5}{get_clean_command()}")
        start_time = time.time()
        
        log_step("Calculating TR Motif Dosage & Copy Number")
        
        output_paths = run_dosage(
            vcf_path=vcf,
            tr_db_path=tr_vcf,
            out_dir=out,
            prefix=prefix,
            threads=threads  
        )
        
        end_time = time.time()
        runtime = end_time - start_time
        
        stats = {
            "Threads": threads,
            "Input VCF": os.path.basename(vcf),
            "TR Database": os.path.basename(tr_vcf),
        }
        
        log_step("Summary")
        log_summary_block(
            cmd=get_clean_command(),
            start=start_time,
            duration=runtime,
            stats=stats
        )
        
    except Exception as e:
        logger.error(f"Error in dosage module: {str(e)}", exc_info=True)
        raise click.Abort()
    