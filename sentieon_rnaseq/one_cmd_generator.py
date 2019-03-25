# coding = utf-8


def score(**kwargs):
    cmd = "{} driver ".format(kwargs['sentieon'])
    cmd += "-t {} ".format(kwargs['threads'])
    cmd += "-i {} ".format(kwargs['bam'])
    cmd += "--algo LocusCollector "
    cmd += "--fun score_info "
    cmd += "{} ".format(kwargs['score_txt'])
    return cmd


def dedup(**kwargs):
    cmd = "{} driver ".format(kwargs['sentieon'])
    cmd += "-t {} ".format(kwargs['threads'])
    cmd += "-i {} ".format(kwargs['bam'])
    cmd += "--algo Dedup "
    cmd += "--rmdup "
    cmd += "--score_info {} ".format(kwargs['score_txt'])
    cmd += "--metrics {} ".format(kwargs['dedup_metric_txt'])
    cmd += "{} ".format(kwargs['dedup_bam'])
    return cmd


def split(**kwargs):
    "sentieon driver -t NUMBER_THREADS -r REFERENCE -i DEDUP_BAM \
--algo RNASplitReadsAtJunction --reassign_mapq 255:60 SPLIT_BAM"
    cmd = "{} driver ".format(kwargs['sentieon'])
    cmd += "-t {} ".format(kwargs['threads'])
    cmd += "-i {} ".format(kwargs['bam'])
    cmd += "-r {} ".format(kwargs['ref_fasta'])
    cmd += "--algo RNASplitReadsAtJunction "
    cmd += "--reassign_mapq {} ".format(kwargs['reassign_mapq'])
    cmd += "{} ".format(kwargs['split_bam'])
    return cmd


def realign(**kwargs):
    "sentieon driver -t NUMBER_THREADS -r REFERENCE \
-i DEDUP_BAM --algo Realigner [-k KNOWN_SITES] REALIGNED_BAM"
    cmd = "{} driver ".format(kwargs['sentieon'])
    cmd += "-t {} ".format(kwargs['threads'])
    cmd += "-i {} ".format(kwargs['bam'])
    cmd += "-r {} ".format(kwargs['ref_fasta'])
    cmd += "--algo Realigner "
    if kwargs['known_sites']:
        cmd += "-k {} ".format(kwargs['known_sites'])
    cmd += "{} ".format(kwargs['realigned_bam'])
    return cmd


def recalibrate(**kwargs):
    "sentieon driver -t NUMBER_THREADS -r REFERENCE \
-i REALIGNED_BAM --algo QualCal [-k KNOWN_SITES] RECAL_DATA.TABLE"
    cmd = "{} driver ".format(kwargs['sentieon'])
    cmd += "-t {} ".format(kwargs['threads'])
    cmd += "-i {} ".format(kwargs['bam'])
    cmd += "-r {} ".format(kwargs['ref_fasta'])
    cmd += "--algo QualCal "
    if kwargs['known_sites']:
        cmd += "-k {} ".format(kwargs['known_sites'])
    cmd += "{} ".format(kwargs['recal_data_table'])
    return cmd


def calling(**kwargs):
    "sentieon driver -t NUMBER_THREADS -r REFERENCE -i RECALIBRATED_BAM \
--algo Genotyper [-d dbSNP] VARIANT_VCF"
    "sentieon driver -t NUMBER_THREADS -r REFERENCE -i REALIGNED_BAM \
-q RECAL_DATA.TABLE --algo Haplotyper --trim_soft_clip \
--call_conf 20 --emit_conf 20 [-d dbSNP] VARIANT_VCF"

    cmd = "{} driver ".format(kwargs['sentieon'])
    cmd += "-t {} ".format(kwargs['threads'])
    cmd += "-i {} ".format(kwargs['bam'])
    cmd += "-r {} ".format(kwargs['ref_fasta'])
    cmd += "-q {} ".format(kwargs['recal_data_table'])
    cmd += "--algo Haplotyper "
    cmd += "--trim_soft_clip "
    cmd += "--call_conf 20 "
    cmd += "--emit_conf 20 "
    if kwargs['dbSNP']:
        cmd += "-d {} ".format(kwargs['dbSNP'])
    cmd += "{} ".format(kwargs['variant_vcf'])
    return cmd
