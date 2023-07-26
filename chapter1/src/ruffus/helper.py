import ruffus.cmdline as cmdline
import yaml, os, glob
from ruffus import *


def common_args(parser):
    parser.add_argument('--config', default="config/config.yaml", help="Configuration file")
    parser.add_argument('--genome', help="Reference genome")
    parser.add_argument('--aln_dir', help="Directory to store bam files")
    return parser


def add_bwa_args(parser):
    parser.add_argument('--idx_dir', help="Directory to store genome indices")
    return parser


def override_args(args, arg_type="all"):
    """
    Overrides config file arguments with run time parameters when provided.
    arg_type can be bwa, bcftools or all
    """

    # Override common arguments to bwa and bcftools
    config = yaml.safe_load(open(args.config))
    args.genome = args.genome if args.genome else config.get('genome')
    args.aln_dir = args.aln_dir if args.aln_dir else config.get('aln_dir')
    args.samples = config.get('data').get('sample_sheet')

    # Override bwa specific arguments
    if arg_type == "bwa" or arg_type == "all":
        # Override config file parameters if provided at runtime.
        args.idx_dir = args.idx_dir if args.idx_dir else config.get('idx_dir')

        # An index file
        args.idx_file = os.path.join(args.idx_dir, os.path.basename(args.genome) + ".bwt")

        # Genome index
        args.index = os.path.splitext(args.idx_file)[0]

        # Get alignment mode and flags
        args.lib = config.get('data').get('library')
        args.sam_flags = config.get('sam').get('sam_flags')
        args.aln_flags = config.get('aligner').get('aln_flags')

    # Override bcftools specific arguments
    if arg_type == "bcftools" or arg_type == "all":
        # Override config file parameters if provided at runtime.
        args.vcf_dir = args.vcf_dir if args.vcf_dir else config.get('vcf_dir')
        args.vcf_type = args.vcf_type if args.vcf_type else config.get('vcf_type')

        # Get variant call flags
        args.pflags = config.get('vcf').get('pileup_flags')
        args.cflags = config.get('vcf').get('call_flags')

    return args


def get_bwa_args():
    """
    Get run time parameters for bwa aligner.
    """
    parser = cmdline.get_argparse(description='Align data to reference using bwa-mem')
    parser = common_args(parser)
    parser = add_bwa_args(parser)
    args = parser.parse_args()
    args = override_args(args, arg_type="bwa")
    return args


def add_bcftools_args(parser):
    parser.add_argument('--vcf_dir', help="Directory to store vcf files")
    parser.add_argument('--vcf_type', help="Type of variant calling, either multi-sample or sample-vcf")
    return parser


def get_bcftools_args():
    """
    Get run time parameters
    """
    parser = cmdline.get_argparse(description='Call variants from data using bcftools')
    parser = common_args(parser)
    parser = add_bcftools_args(parser)
    args = parser.parse_args()
    args = override_args(args, arg_type="bcftools")
    return args


def get_variant_args():
    parser = cmdline.get_argparse(description='Variant calling pipeline')
    parser = common_args(parser)
    parser = add_bwa_args(parser)
    parser = add_bcftools_args(parser)
    args = parser.parse_args()
    # print(args)
    # print("****")
    args = override_args(args, "all")

    # print(args)
    # print("@@@@@")
    return args


def pipeline_dry_run():
    print("\nThe different tasks in the pipeline are \n")
    print("index_genome: Index reference genome")
    print("align_reads: Align data to reference genome")
    print("index_bam: Index bam files")
    print("call_variants: Call variants")
    print("\n")


def read_data(df, mode):
    """
    Reads data as a list of lists
    """
    store = list()
    r1, r2 = "", ""
    for i in df.index:
        r1 = df.loc[i, 'read1']
        if mode == "PE":
            r2 = df.loc[i, 'read2']
        store.append([r1, r2]) if mode == "PE" else store.append([r1])
    return store


def get_prefix_suffix(fpath, sample):
    # suffix is everything from the last underscore.
    path, name = os.path.split(fpath)
    a, ext = os.path.splitext(name)
    vals = a.split('_')
    suffix = '_' + vals[-1] + ext
    return path, suffix


def read_bam(aln_dir, samples=[]):
    store = []
    if not samples:
        return glob.glob(f"{aln_dir}/*.bam")
    for s in samples:
        store.extend(glob.glob(f"{aln_dir}/{s}*.bam"))
    return store


def run(args):
    if args.just_print:
        pipeline_dry_run()
        pipeline_printout(verbose=3)
    else:
        cmdline.run(args, checksum_level=3)
        # pipeline_run()
