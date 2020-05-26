library(optparse)
# make options
opt_list = list(
    make_option(c('-e', '--expr'), help='expr matrix file'),
    make_option(c('-g', '--group'), help='sample group file'),
    make_option(c('-c', '--contrast'), help='comparison info file'),
    make_option(c('-p', '--proportion'), default=0.1, help='numeric value between 0 and 1, assumed proportion of genes which are differentially expressed', type='numeric'),
    make_option(c('-f', '--fold_change'), default=2.0, help='fold change cutoff, if = -1, log2fc_cutoff = mean(abs(result$logFC)) + 2*sd(abs(result$logFC))', type='numeric'),
    make_option(c('-s', '--stat_cutoff'), default=0.05, help='pvalue or adjust pvalue cutoff', type='numeric'),
    make_option(c('-t', '--type'), default='padjust', help='use uncorrected pvalue if set to pvalue'),
    make_option(c('-x', '--x_paired'), default=0, help='1 for paired t-test, 0 for unpaired t-test', type='integer'),
    make_option(c('-o', '--out_prefix'), default='', help='out prefix')
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
    ##step4 生成所有基因的检验结果报告, 不排序，方便后面添加其他列
    result = topTable(fit2, coef=1, n=Inf, adjust.method='BH', sort.by="none")
}

group_matrix = read.table(opt$g, header=T)
contrast = read.table(opt$c, header=T)
all_expr_table = read.table(opt$e, header=T, row.names=1)
for (i in 1:nrow(contrast)){
    ctrl = as.character(contrast[i, 1])
    test = as.character(contrast[i, 2])
    test_samples = group_matrix[which(group_matrix[, 2] == test), 1]
    ctrl_samples = group_matrix[which(group_matrix[, 2] == ctrl), 1]
    if (length(ctrl_samples) <= 1 | length(test_samples) <= 1){
        next
    }
    print(paste(ctrl, '(', length(ctrl_samples), ')', '_vs_', test, '(', length(test_samples), ')', sep=''))
    group_list = append(rep(ctrl, length(ctrl_samples)), rep(test, length(test_samples)))

    # ctrl stat
    ctrl_mean_value = apply(all_expr_table[, as.vector(ctrl_samples)], 1, mean)
    ctrl_median_value = apply(all_expr_table[, as.vector(ctrl_samples)], 1, median)
    ctrl_sd_value = apply(all_expr_table[, as.vector(ctrl_samples)], 1, sd)

    # test stat
    test_mean_value = apply(all_expr_table[, as.vector(test_samples)], 1, mean)
    test_median_value = apply(all_expr_table[, as.vector(test_samples)], 1, median)
    test_sd_value = apply(all_expr_table[, as.vector(test_samples)], 1, sd)
    # print(head(test_mean_value, 30))
    # do ebayes test
    expr_matrix = all_expr_table[, append(as.vector(ctrl_samples), as.vector(test_samples))]
    contrast_exp = paste(test, '-', ctrl)
    result = diff_test(expr_matrix, group_list, contrast_exp, proportion=opt$p)
    # print(head(result, 30))
    ctrl_num = length(ctrl_samples) + 1
    # 如果输入的是log2转换的数据，那么t-test的假设或许存在问题？
    ttest_pvalue = apply(expr_matrix, 1,
        function(x){
            t.test(as.numeric(x[1:length(ctrl_samples)]),
                as.numeric(x[ctrl_num:dim(expr_matrix)[2]]),
                var.equal=F,
                paired=opt$x)$p.value
        }
    )

    result$ttest_pvalue = ttest_pvalue
    result$adj.ttest_pvalue = as.vector(p.adjust(ttest_pvalue, method='BH'))

    # add mean and median and sd columns
    result[, paste('mean.', ctrl, sep='')] = ctrl_mean_value
    result[, paste('mean.', test, sep='')] = test_mean_value
    result[, paste('median.', ctrl, sep='')] = ctrl_median_value
    result[, paste('median.', test, sep='')] = test_median_value
    result[, paste('std.', ctrl, sep='')] = ctrl_sd_value
    result[, paste('std.', test, sep='')] = test_sd_value

    if (opt$f == -1){
        fc_cutoff = mean(abs(result$logFC)) + 2*sd(abs(result$logFC))
    }else{
        fc_cutoff = abs(log2(opt$f))
    }
    if (opt$t == 'padjust'){
        result$significant = (abs(result$'logFC') >= fc_cutoff) & (result[, 'adj.P.Val'] <= opt$s)
        result = result[order(result$adj.P.Val), ]
        print(paste('cutoff: |log2fc| >= ', fc_cutoff, " & pajust", ' <= ', opt$s, sep=''))
    }else if (opt$t == 'ttest') {
        result$significant = (abs(result$'logFC') >= fc_cutoff) & (result[, 'ttest_pvalue'] <= opt$s)
        result = result[order(result[, "ttest_pvalue"]), ]
        print(paste('cutoff: |log2fc| >= ', fc_cutoff, " & ttest_pvalue", ' <= ', opt$s, sep=''))
    }
    else{
        result$significant = (abs(result$'logFC') >= fc_cutoff) & (result[, 'P.Value'] <= opt$s)
        result = result[order(result[, "P.Value"]), ]
        print(paste('cutoff: |log2fc| >= ', fc_cutoff, " & pvalue", ' <= ', opt$s, sep=''))
    }
    result$regulate = 'up'
    result[result$logFC < 0, 'regulate'] = 'down'
    # output
    out_name = paste(opt$o, ctrl, '_vs_', test, '.xls', sep='')
    write.table(result, out_name, sep='\t', col.names=NA, quote=FALSE)
    out_name = paste(opt$o, ctrl, '_vs_', test, '.DE.list', sep='')
    deg = rownames(result[result$significant==T, ])
    print(paste('DEG number: ', length(deg), sep=''))
    deg_reg = cbind(deg, result[result$significant==T, 'regulate'])
    write.table(deg_reg, out_name, sep='\t', col.names=F, quote=FALSE, row.names=F)
}
