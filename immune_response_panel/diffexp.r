library(optparse)
# make options
opt_list = list(
    make_option(c('-e', '--expr'), help='expr matrix file'),
    make_option(c('-g', '--group'), help='sample group file'),
    make_option(c('-c', '--contrast'), help='comparison info file'),
    make_option(c('-p', '--proportion'), default=0.1, help='comparison info file', type='numeric'),
    make_option(c('-f', '--fold_change'), default=2.0, help='fold change cutoff', type='numeric'),
    make_option(c('-s', '--stat_cutoff'), default=0.05, help='pvalue or adjust pvalue cutoff', type='numeric'),
    make_option(c('-t', '--type'), default='padjust', help='use uncorrected pvalue if set to pvalue')
)

# make options
opt_parser = OptionParser(option_list=opt_list)
opt = parse_args(opt_parser)
if (is.null(opt$e)){
    stop("expr file must be provided!")
}
if (is.null(opt$g)){
    stop("group file must be provided!")
}
if (is.null(opt$c)){
    stop('contrast info file must be provided!')
}

# -----------main--------------
library(limma)

diff_test <-function(expr_matrix, group_list, contrast_exp,
    proportion=0.1, trend=TRUE, robust=TRUE){
    design <- model.matrix(~0+factor(group_list))
    colnames(design) = levels(factor(group_list))
    rownames(design) = colnames(expr_matrix)
    contrast.matrix <- makeContrasts(contrast_exp, levels=design)
    fit <- lmFit(expr_matrix, design)
    ##step2 根据对比模型进行差值计算
    fit2 <- contrasts.fit(fit, contrast.matrix)
    ##step3 贝叶斯检验
    fit2 <- eBayes(fit2, proportion=proportion, trend=trend, robust=robust)
    ##step4 生成所有基因的检验结果报告
    result = topTable(fit2, coef=1, n=Inf)
}

group_matrix = read.table(opt$g, header=T)
contrast = read.table(opt$c, header=T)
all_expr_table = read.table(opt$e, header=T, row.names=1)
for (i in 1:nrow(contrast)){
    ctrl = as.character(contrast[i, 1])
    test = as.character(contrast[i, 2])
    test_samples = group_matrix[which(group_matrix[, 2] == test), 1]
    ctrl_samples = group_matrix[which(group_matrix[, 2] == ctrl), 1]
    print(paste(ctrl, '(', length(ctrl_samples), ')', '_vs_', test, '(', length(test_samples), ')', sep=''))
    group_list = append(rep(ctrl, length(ctrl_samples)), rep(test, length(test_samples)))
    expr_matrix = all_expr_table[, append(as.vector(ctrl_samples), as.vector(test_samples))]
    contrast_exp = paste(test, '-', ctrl)
    result = diff_test(expr_matrix, group_list, contrast_exp, proportion=opt$p)
    if (opt$t == 'padjust'){
        result$significant = (abs(result$'logFC') >= log2(opt$f)) & (result[, 'adj.P.Val'] <= opt$s)
        result[order(result$adj.P.Val), ]
    }else{
        result$significant = (abs(result$'logFC') >= log2(opt$f)) & (result[, 'P.Value'] <= opt$s)
        result[order(result[, "P.Value"]), ]
    }
    out_name = paste(ctrl, '_vs_', test, '.xls', sep='')
    write.table(result, out_name, sep='\t', col.names=NA, quote=FALSE)
}
