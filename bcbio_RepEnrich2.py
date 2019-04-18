import logging
import os
from argparse import ArgumentParser
from bcbio.distributed.transaction import file_transaction
from bcbio.utils import file_exists, safe_makedir
from bcbio.bam import sam_to_bam, is_paired
import yaml
from subprocess import check_call

python_path = "/home/rdk4/local/share/bcbio/anaconda/envs/python2/bin/python"

logging.basicConfig(level=logging.DEBUG)

def get_files(bcbio_config):
    with open(args.bcbio_config) as in_handle:
        bcbio_config = yaml.load(in_handle)
    filelist = [x["files"] for x in bcbio_config["details"]]
    samplenames = [x["description"] for x in bcbio_config["details"]]
    return {samplename: files for samplename, files in zip(samplenames, filelist)}

def run_bowtie2(files, samplename, bowtie_index, threads, out_dir):
    out_stem = os.path.join(out_dir, samplename)
    multimap_file = out_stem + "-multimap.fastq"
    out_sam = out_stem + ".sam"
    if file_exists(out_sam):
        logging.info(f"{out_sam} already exists, skipping bowtie1 alignment.")
        return out_sam
    base_cmd = "bowtie2 -q -p {threads} -x {bowtie_index} "
    if len(files) == 2:
        fq1, fq2 = files
        cmd = base_cmd + "-1 {fq1} -2 {fq2} -S {tx_out_sam} "
    else:
        fq1 = files[0]
        cmd = base_cmd + "-U {fq1} -S {tx_out_sam}"
    with file_transaction([multimap_file, out_sam]) as tx_out_files:
        tx_multimap_file = tx_out_files[0]
        tx_out_sam = tx_out_files[1]
        cmd = cmd.format(**locals())
        logging.info(f"Running bowtie1 on {files}.")
        check_call(cmd, shell=True)
    return out_sam

def subset_reads(bam_file, samplename, rep2_path):
    repenrich = os.path.join(rep2_path, "RepEnrich2_subset.py")
    out_prefix = os.path.join("subset", samplename)
    unique = out_prefix + "_unique.bam"
    fq1 = out_prefix + "_multimap_R1.fastq"
    fq2 = out_prefix + "_multimap_R2.fastq"
    if is_paired(bam_file):
        paired = "TRUE"
    else:
        paired = "FALSE"
    if not file_exists(unique):
        cmd = f"{python_path} {repenrich} {bam_file} 30 {out_prefix} --pairedend {paired}"
        check_call(cmd, shell=True)
    if file_exists(fq2):
        return [unique, fq1, fq2]
    else:
        return [unique, fq1]

def run_repenrich2(files, samplename, annotation, rep2_path, rep2_setup, threads):
    if len(files) == 3:
        bam_file, fq1, fq2 = files
    else:
        bam_file, fq1 = files
        fq2 = None
    out_dir = os.path.join("rep2")
    repenrich = os.path.join(rep2_path, "RepEnrich2.py")
    out_dir = os.path.join("rep2", samplename)
    cmd = f"{python_path} {repenrich} {annotation} {out_dir} {samplename} {rep2_setup} {fq1} {bam_file} --cpus {threads} "
    if fq2:
        cmd += f"--fastqfile2 {fq2} --pairedend TRUE"
    check_call(cmd, shell=True)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("bcbio_config", help="final bcbio configuration file.")
    parser.add_argument("bowtie_index", help="location of bowtie2 index")
    parser.add_argument("annotation", help="repeat annotation")
    parser.add_argument("rep2_setup", help="repenrich2 setup directory")
    parser.add_argument("--threads", default=1, help="Number of threads to use.")
    parser.add_argument("--outdir", default="align", help="output directory")
    parser.add_argument("--rep2_path", default="RepEnrich2", help="path to RepEnrich2 code")
    args = parser.parse_args()
    filedict = get_files(args.bcbio_config)
    safe_makedir(args.outdir)
    config = config["algorithm"] = {"num_cores": args.threads}
    repenrich = "RepEnrich2/RepEnrich.py"
    annotation = "metadata/hg38_repeatmasker_clean.txt"
    repenrich_setup = "metadata/RepEnrich2_setup_hg38"
    for samplename, files in filedict.items():
        out_dir = os.path.join(args.outdir, "repenrich", samplename)
        logging.info(f"Aligning {samplename} to {args.bowtie_index}.")
        out_sam = run_bowtie2(files, samplename, args.bowtie_index, args.threads, args.outdir)
        logging.info(f"Converting {out_sam} to BAM format.")
        out_bam = sam_to_bam(out_sam, config)
        logging.info(f"Subsetting {out_bam} into unique and multimapped reads.")
        subset_files = subset_reads(out_bam, samplename, args.rep2_path)
        logging.info(f"Running RepEnrich2 on {samplename}.")
        run_repenrich2(subset_files, samplename, args.annotation, args.rep2_path, args.rep2_setup, args.threads)
        logging.info(f"Finished {samplename}.")
