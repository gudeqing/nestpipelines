import os
from subprocess import check_call


def annovar_annotation(vcf):
    """
    调用annovar注释vcf而已
    :return:
    """
    if os.path.exists(f'{vcf}.hg19_multianno.vcf'):
        print('annotated result already exists, skip now!')
        return f'{vcf}.hg19_multianno.vcf', f'{vcf}.hg19_multianno.txt'
    annovar = "/nfs2/software/annovar/ANNOVAR_2019.10.24/table_annovar.pl"
    annovardb = "/nfs2/database/annovar/humandb/"
    cmd = f"{annovar} {vcf} {annovardb} -buildver hg19 -out {vcf} -remove "
    g_r_pro = ['refGeneWithVer', 'cytoBand']
    f_pro = [
        'avsnp150',
        'cosmic90', 'CGI',
        'IntOGen', 'TCGA', 'icgc21', 'CancerHotspots', 'nci60', 'clinvar_20200316',
        'esp6500siv2_all', '1000g2015aug_all', 'exac03', 'gnomad_exome', 'gme',
        'intervar_20180118', 'dbnsfp33a'
    ]
    protocols = ','.join(g_r_pro + f_pro)
    operations = 'g,r,' + ','.join(['f'] * len(f_pro))
    args = ','.join(['-hgvs'] * len(g_r_pro + f_pro))
    cmd += f'-protocol {protocols} '
    cmd += f'-operation {operations} '
    cmd += f"--argument '{args}' "
    cmd += "-nastring . -vcfinput --dot2underline --thread 12 --maxgenethread 20 "
    check_call(cmd, shell=True)
    return f'{vcf}.hg19_multianno.vcf', f'{vcf}.hg19_multianno.txt'



if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['annovar_annotation'])
