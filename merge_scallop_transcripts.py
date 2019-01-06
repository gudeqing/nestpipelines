# coding=utf-8
import os
from subprocess import check_call
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-gtf', help='gtf file path', nargs='+', required=True)
parser.add_argument('-refgtf', required=True, help='reference gtf')
parser.add_argument('-refgenome', required=True, help='reference genome fasta')
parser.add_argument('-outdir', help='output dir', default=os.getcwd())
parser.add_argument('-reftrans', default=None, 
                    help='reference transcripts fasta. If None, use refgtf and refgenome to extract it')
parser.add_argument('-gtfmerge', help='gtfmerge path, https://github.com/Kingsford-Group/rnaseqtools',
                    default='gtfmerge')
parser.add_argument('-gffcompare', help='gffcompare path', default='gffcompare')
parser.add_argument('-gtfcuff', help='gtfcuff path', default = 'gtfcuff')
parser.add_argument('-gffread', help='gffread path', default = 'gffread')
args = parser.parse_args()


def write_gtf_list(gtf_list, out_file):
    with open(out_file, 'w') as f:
        for gtf in gtf_list:
            f.write(gtf+'\n')
    return out_file


cmds = list()

gtf_list_file = write_gtf_list(args.gtf, os.path.join(args.outdir, 'gtf.list'))
merged_file = os.path.join(args.outdir, 'merged_assembled.gtf')
cmd = '{} union {} {} -n '.format(args.gtfmerge, gtf_list_file, merged_file)
cmds.append(cmd)

out_prefix = os.path.join(args.outdir, 'gffall')
cmd = '{} -o {} -r {} {}'.format(args.gffcompare, out_prefix, args.refgtf, merged_file)
cmds.append(cmd)

novel_gtf = os.path.join(args.outdir, 'novel.gtf')
cmd = '{} puniq {}.{}.tmap {} {} {}'.format(args.gtfcuff, out_prefix, merged_file, merged_file, args.refgtf, novel_gtf)
cmds.append(cmd)

transcriptome = args.reftrans
if not transcriptome or not os.path.exists(transcriptome):
    transcriptome = os.path.join(args.outdir, 'reference.transcriptome.fa')
    cmd = '{} {} -g {} -w {}'.format(args.gffread, args.refgtf, args.refgenome, transcriptome)
    cmds.append(cmd)

new_transcripts = os.path.join(args.outdir, 'novel.transcripts.fa')
cmd = '{} {} -g {} -w {}'.format(args.gffread, novel_gtf, args.refgenome, new_transcripts)
cmds.append(cmd)

all_transcripts = os.path.join(args.outdir, 'all.transcripts.fa')
cmd = 'cat {} {} > {}'.format(new_transcripts, transcriptome, all_transcripts)
cmds.append(cmd)

all_gtf = os.path.join(args.outdir, 'all.transcripts.gtf')
cmd = 'cat {} {} > {}'.format(args.refgtf, novel_gtf, all_gtf)
cmds.append(cmd)

with open('process.cmds.list', 'w') as f:
    for cmd in cmds:
        f.write(cmd+"\n")

for cmd in cmds:
    check_call(cmd, shell=True)

