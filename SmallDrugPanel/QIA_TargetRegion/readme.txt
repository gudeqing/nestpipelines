从cosmic的all coding mutation: CosmicCodingMuts.normal.vcf.gz  提取目标基因的突变
target.mutation.list 是彬彬这边整理的目标检测位点，由于对应到蛋白变化时，cosmic可以找到蛋白变化相同但碱基变化不同的突变
根据target.mutation.list里的AA或Coding变化可以在cosmic找到71条突变target.genes.target_mutation.vcf
我们的目标突变列表中有四个AA变化出现多次：
      4 p.G776delinsVC
      4 p.Q61H
      2 p.E746_A750del （cosmic中出现3次）
      2 p.G12D  （cosmic中出现4次）

发现我们关注的突变位点仅仅使用小部分的genePrimer就可以捕获，因此对于其他在捕获区域的cosmic突变也可以纳入进来
target_genes.cosmic_coding_mutation.in_ROI.vcf 共包含2799个突变，这些突变都在捕获区域， 这个文件来源如下：
 >  bedtools intersect -a target.genes.vcf -b QIAseq_DNA_panel.CDHS-33088Z-140.rangeOfInterest.bed -wa -header > tmp
 >  bedtools intersect -a tmp -b ../PrimerInfo/primer_targeted_250bp.bed -wa -header -f 1 > tmp2
 >  grep -v '#' tmp2 | sort |uniq > tmp3
 >  grep '#' tmp2 >> target_genes.cosmic_coding_mutation.in_ROI.vcf
 >  grep -v '#' tmp2 | sort |uniq >> target_genes.cosmic_coding_mutation.in_ROI.vcf

