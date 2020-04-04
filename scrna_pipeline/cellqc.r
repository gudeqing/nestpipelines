library(optparse)
# make options
opt_list = list(
    make_option(c('-c', '--count'), help='count matrix file'),
    make_option(c('-e', '--expr'), help='expr matrix file'),
    make_option(c('-m', '--control'), default='/nfs2/database/gencode_v29/chrM.gene.list',
        help='a file containing control gene list')
)

# make options
opt_parser = OptionParser(option_list=opt_list)
opt = parse_args(opt_parser)
if (is.null(opt$count)){
    stop("count file must be provided!")
}
if (is.null(opt$expr)){
    stop("expr file must be provided!")
}
# load library
library(ggplot2)
library(SingleCellExperiment)
library(scater)
library(scran)

# 构建SCE对象
mt_genes = read.table(opt$control)  # 读取线粒体基因
counts <- read.table(opt$count, sep = "\t", row.name=1, header=1)
rownames(counts) <- gsub('\\..*', '', rownames(counts))
exprs = read.table(opt$expr, sep='\t', row.names=1, header=1)
rownames(exprs) <- gsub('\\..*', '', rownames(exprs))
sce <- SingleCellExperiment(assays=list(counts=as.matrix(counts), exprs=as.matrix(exprs)))
isSpike(sce, "MT") <- rownames(sce) %in% mt_genes[,1]

# 计算QC metics，主要是统计count数，文库大小，MT基因的count占比，没有count的基因数目
sce <- calculateQCMetrics(sce, exprs_values = "counts", feature_controls = list(MT=isSpike(sce, "MT")))
write.table(colData(sce),'qc.metrics.xls', sep='\t', col.names=NA, quote=FALSE)

# Classification of cell cycle phase
hsa.pairs = readRDS(system.file("exdata", "human_cycle_markers.rds", package='scran'))
assignments <- cyclone(sce, hsa.pairs, gene.names=rownames(sce), assay.type="exprs")
sce$phases <- assignments$phases
score_table = as.data.frame(assignments$normalized.scores, colnames(sce))
phase_table = as.data.frame(assignments$phases, colnames(sce))
write.table(score_table, 'cycle_score.xls', sep='\t', col.names=NA, quote=FALSE)
write.table(phase_table, 'cycle_phases.xls', sep='\t', col.names=NA, quote=FALSE)
