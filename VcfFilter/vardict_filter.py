"""
测序错误过滤突变的思路：
假设事先知道测序过程中发生的碱基转换的概率，例如，知道A被错误测成T的概率为0.001，
针对一个突变位点A-T，其AF=0.008，是否可以将当前位点判定为假阳性呢?
假设还知道当前位点深度为1000，根据这1000 Reads，在95%的置信水平下，估计出来的测序错误概率与真实的错误概率偏差可以使用如下公式推算：
E = 1.96/(N/(P*(1-P)))**0.5
99的置信水平：E = 2.58/(N/(P*(1-P)))**0.5
代入 N=1000, P=0.001得 E=0.00196，也就是，如果A>T是测序错误导致的，那么这个错误频率的置信区间为
[P-E, P+E] = [0, 0.00296], 由于现在测得的AF=0.008 已经远在置信区间水平以外，所以不应该过滤掉。

思考2：
当有对照样本时，且对照样本中该位点也存在突变，设其突变频率为ctrP, 如果使用上述方法判定该突变为假阳，
则可以使用该P作为真实测序错误率对测试样本进行上述类似的统计.

思考3：
当没有对照样本时，而只有阴性设计样本如NA12878时，我们可以通过统计阴性样本中碱基测序错误率作为上述过滤思路的输入

思考4：
当既没有对照样本，又没有阴性样本，我们假设肿瘤样本中call出来的低频突变中绝大部分为假阳性突变，那么我们也可以根据肿瘤样本粗略估计测序错误。
最后可以用肿瘤样本自己作为输入并实现过滤。
"""
import json
# import pandas as pd
import pysam
import scipy.stats as stats
import statistics
"""
要求vcf每一行只包含一个突变，这个可以通过bcftools norm 快速实现
变异类型的分类
https://www.ebi.ac.uk/training-beta/online/courses/human-genetic-variation-introduction/what-is-genetic-variation/types-of-genetic-variation/
"""


class VardictFilter():
    def __init__(self, vcf_path):
        self.vcf = pysam.VariantFile(vcf_path)
        self.vcf_path = vcf_path

    def chrom_name_is_numeric(self):
        with open(self.vcf_path) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                else:
                    if line.split()[0].isnumeric():
                        return True
                    else:
                        return False

    def add_contig_header(self, ref='hg19'):
        contig_info = [
            f'##assembly={ref}',
            "##contig=<ID=chr1,length=249250621>",
            "##contig=<ID=chr2,length=243199373>",
            "##contig=<ID=chr3,length=198022430>",
            "##contig=<ID=chr4,length=191154276>",
            "##contig=<ID=chr5,length=180915260>",
            "##contig=<ID=chr6,length=171115067>",
            "##contig=<ID=chr7,length=159138663>",
            "##contig=<ID=chr8,length=146364022>",
            "##contig=<ID=chr9,length=141213431>",
            "##contig=<ID=chr10,length=135534747>",
            "##contig=<ID=chr11,length=135006516>",
            "##contig=<ID=chr12,length=133851895>",
            "##contig=<ID=chr13,length=115169878>",
            "##contig=<ID=chr14,length=107349540>",
            "##contig=<ID=chr15,length=102531392>",
            "##contig=<ID=chr16,length=90354753>",
            "##contig=<ID=chr17,length=81195210>",
            "##contig=<ID=chr18,length=78077248>",
            "##contig=<ID=chr19,length=59128983>",
            "##contig=<ID=chr20,length=63025520>",
            "##contig=<ID=chr21,length=48129895>",
            "##contig=<ID=chr22,length=51304566>",
            "##contig=<ID=chrMT,length=16569>",
            "##contig=<ID=chrX,length=155270560>",
            "##contig=<ID=chrY,length=59373566>",
        ]
        if self.chrom_name_is_numeric():
            contig_info = [x.replace('=chr', '=') for x in contig_info]
        for line in contig_info:
            self.vcf.header.add_line(line)

    def poll_error_conf(self, error_rate, depth, z=2.58):
        e = z/(depth/(error_rate*(1-error_rate)))**0.5
        lower = 0 if (error_rate - e <= 0) else (error_rate - e)
        upper = error_rate + e
        return lower, upper

    def qual_to_error_rate(self, base_qual):
        # Phred = -10 * log(error_p)
        error_prob = 10**(base_qual*(-0.1))
        return error_prob

    def pass_seq_error(self, record, sample, seq_error:float=None, z:float=2.58, factor:float=1.0, read_len=150):
        dp = record.samples[sample]['DP']
        if type(dp) != int:
            dp = sum(dp)
        if 'HIAF' in record.info and 'AF' in record.info:
            # vardict style
            af = min([record.info['HIAF'], record.info['AF'][0]])
        elif 'AF' in record.samples[sample]:
            af = record.samples[sample]['AF'][0]
        elif 'FREQ' in record.samples[sample]:
            # for varscan2 style
            af = record.samples[sample]['FREQ']
            if '%' in af:
                af = float(af.strip('%'))*0.01
        else:
            raise Exception('No AF or FREQ field found !')
        if seq_error is None:
            # 认为提供的是phred base quality
            error_rate = 1e-6
            for each in ['NM', 'QUAL', 'MBQ']:
            # for each in ['QUAL', 'MBQ']:
                # QUAL for vardict
                # MBQ for mutect2
                if each in record.info:
                    if each == 'NM':
                        if type(record.info[each]) == float:
                            # vardict style, 'NM':"Mean mismatches in reads"
                            # 根据平均错配数量推算测序错误率的上限，这个错误率肯定偏大，因为其包含了真实突变导致的错配
                            # 假设包含真实突变dna的含量占比0.5，假设一个碱基发生突变，理论mismatch为0.5
                            error_rate = (record.info['NM'])/read_len * 0.3
                            if error_rate <= 1e-6:
                                error_rate = 1e-6
                            break
                    else:
                        error_rate = self.qual_to_error_rate(record.info[each])
                        break
            print(f'error rate for {record.start}:{record.ref}>{record.alts[0]}:', error_rate)
        else:
            error_rate = seq_error

        lower, upper = self.poll_error_conf(error_rate, dp, z)
        # print(dp, r.qual, error_rate, lower, upper)
        if af >= upper*factor:
            return True, lower, upper
        else:
            return False, lower, upper

    def pass_depth_bias(self, record, cutoff=0.05, info_field='', normal=None, tumor=None):
        """
        tumor vs normal depth bias
                    tumour  normal
        ref_depth   2475    2269
        var_depth   28    14
        odds, pvalue = stats.fisher_exact([[2475, 28], [2269, 14]], alternative='less')
        :param record: pysam.variant.record
        :param cutoff: 0.05
        :param info_field:
        :param fmt_filed:
        :param sample:
        :return:
        """
        passed = True
        pvalue = None
        possible_info_fields = [info_field] + ['SPV']
        for field in possible_info_fields:
            if field in record.info:
                if type(record.info[field]) == float:
                    pvalue = record.info[field]
                    break

        def get_depth(sample, field='DP'):
            dp = record.samples[sample][field]
            if type(dp) != int:
                if field == 'AD':
                    dp = dp[1]
                else:
                    dp = sum(dp)
            return dp

        if (normal and tumor) and (pvalue is None):
            normal_dp = get_depth(normal, 'DP')
            normal_ad = get_depth(normal, 'AD')
            tumor_dp = get_depth(tumor, 'DP')
            tumor_ad = get_depth(tumor, 'AD')
            odds, pvalue = stats.fisher_exact(
                [[tumor_dp-tumor_ad, tumor_ad], [normal_dp-normal_ad, normal_ad]],
                alternative='less'
            )
        # judge
        if pvalue is not None:
            if pvalue >= cutoff:
                passed = False
        else:
            print('tumor vs normal depth bias filter is not applied for no valid field is found !')
            pass

        return passed, pvalue

    def pass_strand_bias(self, record, info_field='', fmt_field='', cutoff=0.005, sample=None):
        passed = True
        pvalue = None
        possible_info_fields = [info_field] + ['SBF']
        for field in possible_info_fields:
            if field in record.info:
                if type(record.info[field]) == float:
                    pvalue = record.info[field]
                    break
        # find info in "format"
        possible_fmt_fields = [fmt_field] + ['SB', 'DP4']
        if (pvalue is None) and (sample is not None):
            ref_fwd, ref_bwd, alt_fwd, alt_bwd = 0, 0, 0, 0
            for field in possible_fmt_fields:
                if field in record.samples[sample]:
                    if type(record.samples[sample][field]) == tuple:
                        ref_fwd, ref_bwd, alt_fwd, alt_bwd = record.samples[sample][field]
                        break
            if sum([ref_fwd, ref_bwd, alt_fwd, alt_bwd]) > 0:
                odds_ratio, pvalue = stats.fisher_exact([[ref_fwd, ref_bwd], [alt_fwd, alt_bwd]])
                if ref_fwd == 583:
                    print(ref_fwd, ref_bwd, alt_fwd, alt_bwd)
                    print(pvalue)

        # judge
        if pvalue is not None:
            if pvalue <= cutoff:
                passed = False
        else:
            print('strand bias filter is not applied for no valid field is found !')
        #
        return passed, pvalue

    def pass_pstd(self, record, cutoff=1e-4):
        """
        vardict中PSTD表示突变在reads中位置的标准差，如果突变在reads中出现的位置一样，那么该值为0，假阳的可能性极高
        通常，alt reads序列完全一致会导致这种情况，这意味着支持该突变的uniq read只有1条
        对于QIA的数据预处理，由于会统一去除primer序列，这使得最终同一个primer捕获的插入片段的测序结果中read1的起始位置都一致
        对于测通的情况，read1和read2包含的序列信息还会完全一致
        所以针对上述过滤可能不大合适，因此改为仅在支持的reads数目<=2的情况使用该过滤条件。
        """
        passed = True
        pstd = 0
        if 'PSTD' in record.info and 'VD' in record.info:
            pstd = record.info['PSTD']
            support_reads = record.info['VD']
            if pstd <= cutoff and support_reads <=2:
                passed = False
        return passed, pstd

    def filtering(self, genome, out_prefix=None, seq_error=None, tumor_index=1, center_size:tuple=None, basic_error_dict=None):
        samples = list(self.vcf.header.samples)
        if len(samples) == 2:
            if tumor_index == 1:
                normal, tumor = samples
            else:
                tumor, normal = samples
        else:
            tumor = samples[0]
            normal = None

        # 先给vcf的header添加新字段定义才能往添加新的字段信息
        self.vcf.header.info.add('LOD', number=2, type='Float',
                                 description='The first value is input error rate which will be used '
                                             'as theoretical frequency to calculate the second value. '
                                             'The second value is the upper value of estimated error rate '
                                             'base on assumption of simple sampling')
        self.vcf.header.filters.add('SeqErrorOrGermline', number=None, type=None, description='likely SeqErrorOrGermline')
        self.vcf.header.filters.add('SeqError', number=None, type=None, description='likely seq error')
        if 'StrandBias' not in self.vcf.header.filters:
            self.vcf.header.filters.add('StrandBias', number=None, type=None, description="severe strand bias")
        if normal is not None:
            if 'DepthBias' not in self.vcf.header.filters:
                self.vcf.header.filters.add('DepthBias', number=None, type=None,
                                            description="not enough depth bias between norm and tumor")
        if self.vcf.header.contigs.__len__() < 1:
            self.add_contig_header()
        if out_prefix is None:
            out_prefix = self.vcf_path.rsplit('.', 1)[0]
        vcf_out = pysam.VariantFile(out_prefix+'.filtered.vcf', "w", header=self.vcf.header)
        vcf_discard = pysam.VariantFile(out_prefix+'.discarded.vcf', "w", header=self.vcf.header)

        # 读取seq error信息，解析成字典，记录每个序列突变成ATCG和''的频率
        if seq_error.replace(".", "", 1).isdigit():
            seq_error = float(seq_error)
            tmp_dict = dict(zip('ATCGID', [seq_error]*6))
            seq_error_dict = dict()
            for i in 'ATCGID':
                seq_error_dict[i] = tmp_dict
            key_left = key_right = 0
        else:
            seq_error_dict = json.load(open(seq_error))
            if center_size:
                key_left, key_right = map(int, center_size)
            else:
                key_len = max(len(x) for x in seq_error_dict['A'].keys())//2
                key_left = key_right = key_len
        if basic_error_dict:
            basic_error_dict = json.load(open(basic_error_dict))

        gn = pysam.FastaFile(genome)
        # print(seq_error_dict.keys())
        lod_list = []
        discard = 0
        total = 0
        for r in self.vcf:
            total += 1
            if '<' in r.alts[0]:
                discard += 1
                # 跳过vardict中输出的特殊突变
                print('skip', r.contig, r.ref, list(r.alts))
                continue
            reasons = []
            # 1.根据测序错误率或germline突变频率过滤
            ctrl_af_as_error_rate = False
            # r.pos正好是1-based, r.start 是0-based
            if len(r.alts[0]) == len(r.ref) == 1:
                # snv
                key = gn.fetch(r.contig, r.start - key_left, r.start + 1 + key_right).upper()
                # error_rate = seq_error_dict[key][r.alts[0]]
                if key in seq_error_dict[r.alts[0]]:
                    error_rate = seq_error_dict[r.alts[0]][key][r.alts[0]]
                else:
                    if basic_error_dict:
                        error_rate = basic_error_dict[r.alts[0]][r.ref][r.alts[0]]
                    else:
                        error_rate = 1e-6
            elif r.alts[0].startswith(r.ref):
                # insertion
                key = gn.fetch(r.contig, r.start - key_left, r.start + 1 + key_right).upper()
                if key in seq_error_dict['I']:
                    error_rate = seq_error_dict['I'][key]['I']
                    if len(r.alts[0])-len(r.ref) >= 3:
                        error_rate = error_rate**2
                else:
                    if basic_error_dict:
                        error_rate = basic_error_dict['I'][r.ref]['I']
                        if len(r.alts[0]) - len(r.ref) >= 3:
                            error_rate = error_rate ** 2
                    else:
                        # 模型中没有该类突变的参考频率
                        error_rate = 1e-6
            elif r.ref.startswith(r.alts[0]):
                # deletion
                key = gn.fetch(r.contig, r.pos - key_left, r.pos + 1 + key_right).upper()
                if 'D' in seq_error_dict['D']:
                    # error_rate = seq_error_dict[key]['']
                    if key in seq_error_dict['D']:
                        error_rate = seq_error_dict['D'][key]['D']
                        if len(r.ref) - len(r.alts[0]) >= 3:
                            error_rate = error_rate ** 2
                    else:
                        if basic_error_dict:
                            error_rate = basic_error_dict['D'][r.ref[1]]['D']
                            if len(r.ref) - len(r.alts[0]) >= 3:
                                error_rate = error_rate ** 2
                        else:
                            # 模型中没有该类突变的参考频率
                            error_rate = 1e-6
                else:
                    error_rate = 1e-6
                # print('del', key, error_rate, r.ref, list(r.alts))
            else:
                # complex
                # error_rate = seq_error_dict[key][r.alts[0][0]]
                key = gn.fetch(r.contig, r.start - key_left, r.start + 1 + key_right).upper()
                # print(r.pos, key, r.ref, r.alts)
                if key in seq_error_dict[r.alts[0][0]]:
                    error_rate = seq_error_dict[r.alts[0][0]][key][r.alts[0][0]]
                else:
                    error_rate = 1e-6

            if normal:
                # 当存在对照样本时，如果某个位点在对照样本也存在突变，且突变频率大于seq_error时，可以把对照样本中的突变频率作为测序错误率进行过滤
                # 这样既可以过滤掉germline突变，也可以过滤掉由于测序错误而导致的突变
                # 由于字段里存储的af通常是四舍五入的结果，所以为了保证精度，根据深度重新计算
                normal_dp = r.samples[normal]['DP']
                if type(normal_dp) == tuple:
                    normal_dp = sum(normal_dp)
                normal_ad = r.samples[normal]['AD']
                if type(normal_ad) == tuple:
                    normal_ad = normal_ad[1]
                normal_af = normal_ad/normal_dp
                if normal_af >= 1:
                    normal_af = 0.999
                # print(normal_af, seq_error)
                if normal_af > (seq_error or 0):
                    ctrl_af_as_error_rate = True
                    error_rate = normal_af
            # judge = self.pass_seq_error(r, tumor, error_rate, z=1.96, factor=1.2)
            judge = self.pass_seq_error(r, tumor, error_rate, z=2.58, factor=1.2)
            if not judge[0]:
                if ctrl_af_as_error_rate:
                    reasons.append('SeqErrorOrGermline')
                else:
                    # print(key, error_rate, r.ref, list(r.alts))
                    reasons.append('SeqError')
                pass

            r.info['LOD'] = (round(error_rate, 5), round(judge[2], 5))
            lod_list.append(judge[2])
            # 2.根据strand bias进行过滤
            judge2 = self.pass_strand_bias(r, cutoff=0.003, sample=tumor)
            if not judge2[0]:
                reasons.append('StrandBias')
                pass

            # 3. depth bias filtering
            if normal is not None:
                judge3 = self.pass_depth_bias(r, cutoff=0.05, normal=normal, tumor=tumor)
                if not judge3[0]:
                    reasons.append('DepthBias')
                    pass

            # 4. position std filtering, 如果突变在reads中出现的位置基本不变,且支持的read<=2，需要过滤掉
            judge4 = self.pass_pstd(r, cutoff=0.00001)
            if not judge4[0]:
                reasons.append('pSTD')

            if reasons:
                discard += 1
                for each in reasons:
                    r.filter.add(each)
                vcf_discard.write(r)
            else:
                vcf_out.write(r)
        vcf_out.close()
        vcf_discard.close()
        print('median LOD:', statistics.median(lod_list))
        print('min LOD:', min(lod_list))
        print('max LOD:', max(lod_list))
        print(f'discard {discard} variants while keep {total-discard} ones!')


def filterVcf(vcf, genome, out_prefix=None,  seq_error=None, tumor_index:int=0,
              center_size:tuple=None, basic_error_dict=None):
    VardictFilter(vcf).filtering(
        out_prefix=out_prefix, seq_error=seq_error, tumor_index=tumor_index, genome=genome,
        center_size=center_size, basic_error_dict=basic_error_dict
    )




if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())








