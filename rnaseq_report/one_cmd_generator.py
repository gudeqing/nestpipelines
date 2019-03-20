def gene_body_coverage(**kwargs):
    cmd = '{} '.format(kwargs['python'])
    cmd += '{} '.format(kwargs['script'])
    cmd += 'gene_body_coverage '
    cmd += '-files {} '.format(kwargs['files'])
    cmd += '-outdir {} '.format(kwargs['outdir'])
    cmd += '-formats {} '.format(kwargs['formats'])
    cmd += '-scale {} '.format(kwargs['scale'])
    cmd += '-height {} '.format(kwargs['height'])
    cmd += '-width {} '.format(kwargs['width'])
    cmd += '-file_from {} '.format(kwargs['file_from'])
    return cmd


def fragment_length(**kwargs):
    cmd = '{} '.format(kwargs['python'])
    cmd += '{} '.format(kwargs['script'])
    cmd += 'fragment_length '
    cmd += '-files {} '.format(kwargs['files'])
    cmd += '-outdir {} '.format(kwargs['outdir'])
    cmd += '-formats {} '.format(kwargs['formats'])
    cmd += '-scale {} '.format(kwargs['scale'])
    cmd += '-height {} '.format(kwargs['height'])
    cmd += '-width {} '.format(kwargs['width'])
    cmd += '-min_len {} '.format(kwargs['min_len'])
    cmd += '-max_len {} '.format(kwargs['max_len'])
    return cmd


def inner_distance(**kwargs):
    cmd = '{} '.format(kwargs['python'])
    cmd += '{} '.format(kwargs['script'])
    cmd += 'inner_distance '
    cmd += '-files {} '.format(kwargs['files'])
    cmd += '-outdir {} '.format(kwargs['outdir'])
    cmd += '-formats {} '.format(kwargs['formats'])
    cmd += '-scale {} '.format(kwargs['scale'])
    cmd += '-height {} '.format(kwargs['height'])
    cmd += '-width {} '.format(kwargs['width'])
    cmd += '-min_dist {} '.format(kwargs['min_dist'])
    cmd += '-max_dist {} '.format(kwargs['max_dist'])
    return cmd


def read_distribution(**kwargs):
    cmd = '{} '.format(kwargs['python'])
    cmd += '{} '.format(kwargs['script'])
    cmd += 'read_distribution '
    cmd += '-files {} '.format(kwargs['files'])
    cmd += '-outdir {} '.format(kwargs['outdir'])
    cmd += '-formats {} '.format(kwargs['formats'])
    cmd += '-scale {} '.format(kwargs['scale'])
    cmd += '-height {} '.format(kwargs['height'])
    cmd += '-width {} '.format(kwargs['width'])
    return cmd


def read_duplication(**kwargs):
    cmd = '{} '.format(kwargs['python'])
    cmd += '{} '.format(kwargs['script'])
    cmd += 'read_duplication '
    cmd += '-files {} '.format(kwargs['files'])
    cmd += '-outdir {} '.format(kwargs['outdir'])
    cmd += '-formats {} '.format(kwargs['formats'])
    cmd += '-scale {} '.format(kwargs['scale'])
    cmd += '-height {} '.format(kwargs['height'])
    cmd += '-width {} '.format(kwargs['width'])
    cmd += '-max_dup {} '.format(kwargs['max_dup'])
    return cmd


def chromosome_read_distribution(**kwargs):
    cmd = '{} '.format(kwargs['python'])
    cmd += '{} '.format(kwargs['script'])
    cmd += 'chromosome_read_distribution '
    cmd += '-files {} '.format(kwargs['files'])
    cmd += '-outdir {} '.format(kwargs['outdir'])
    cmd += '-formats {} '.format(kwargs['formats'])
    cmd += '-scale {} '.format(kwargs['scale'])
    cmd += '-height {} '.format(kwargs['height'])
    cmd += '-width {} '.format(kwargs['width'])
    cmd += '-top {} '.format(kwargs['top'])
    return cmd


def exp_saturation(**kwargs):
    cmd = '{} '.format(kwargs['python'])
    cmd += '{} '.format(kwargs['script'])
    cmd += 'exp_saturation '
    cmd += '-files {} '.format(kwargs['files'])
    cmd += '-outdir {} '.format(kwargs['outdir'])
    cmd += '-formats {} '.format(kwargs['formats'])
    cmd += '-scale {} '.format(kwargs['scale'])
    cmd += '-height {} '.format(kwargs['height'])
    cmd += '-width {} '.format(kwargs['width'])
    cmd += '-outlier_limit {} '.format(kwargs['outlier_limit'])
    return cmd


def exp_pca(**kwargs):
    cmd = '{} '.format(kwargs['python'])
    cmd += '{} '.format(kwargs['script'])
    cmd += 'exp_pca '
    cmd += '-exp_table {} '.format(kwargs['exp_table'])
    cmd += '-row_sum_cutoff {} '.format(kwargs['row_sum_cutoff'])
    cmd += '-exp_cutoff {} '.format(kwargs['exp_cutoff'])
    cmd += '-cv_cutoff {} '.format(kwargs['cv_cutoff'])
    cmd += '-explained_ratio {} '.format(kwargs['explained_ratio'])
    if kwargs['group_dict'] is not None:
        cmd += '-group_dict {} '.format(kwargs['group_dict'])
    cmd += '-outdir {} '.format(kwargs['outdir'])
    cmd += '-formats {} '.format(kwargs['formats'])
    cmd += '-scale {} '.format(kwargs['scale'])
    cmd += '-height {} '.format(kwargs['height'])
    cmd += '-width {} '.format(kwargs['width'])
    return cmd


def exp_density(**kwargs):
    cmd = '{} '.format(kwargs['python'])
    cmd += '{} '.format(kwargs['script'])
    cmd += 'exp_density '
    cmd += '-exp_table {} '.format(kwargs['exp_table'])
    cmd += '-row_sum_cutoff {} '.format(kwargs['row_sum_cutoff'])
    cmd += '-exp_cutoff {} '.format(kwargs['exp_cutoff'])
    cmd += '-outdir {} '.format(kwargs['outdir'])
    cmd += '-formats {} '.format(kwargs['formats'])
    cmd += '-scale {} '.format(kwargs['scale'])
    cmd += '-height {} '.format(kwargs['height'])
    cmd += '-width {} '.format(kwargs['width'])
    return cmd


def CollectAlignmentSummaryMetrics(**kwargs):
    cmd = '{} '.format(kwargs['python'])
    cmd += '{} '.format(kwargs['script'])
    cmd += 'CollectAlignmentSummaryMetrics '
    cmd += '-files {} '.format(kwargs['files'])
    cmd += '-outdir {} '.format(kwargs['outdir'])
    cmd += '-formats {} '.format(kwargs['formats'])
    cmd += '-scale {} '.format(kwargs['scale'])
    cmd += '-height {} '.format(kwargs['height'])
    cmd += '-width {} '.format(kwargs['width'])
    return cmd


def CollectInsertSizeMetrics(**kwargs):
    cmd = '{} '.format(kwargs['python'])
    cmd += '{} '.format(kwargs['script'])
    cmd += 'CollectInsertSizeMetrics '
    cmd += '-files {} '.format(kwargs['files'])
    cmd += '-outdir {} '.format(kwargs['outdir'])
    cmd += '-formats {} '.format(kwargs['formats'])
    cmd += '-scale {} '.format(kwargs['scale'])
    cmd += '-height {} '.format(kwargs['height'])
    cmd += '-width {} '.format(kwargs['width'])
    return cmd


def CollectRnaSeqMetrics(**kwargs):
    cmd = '{} '.format(kwargs['python'])
    cmd += '{} '.format(kwargs['script'])
    cmd += 'CollectRnaSeqMetrics '
    cmd += '-files {} '.format(kwargs['files'])
    cmd += '-outdir {} '.format(kwargs['outdir'])
    cmd += '-formats {} '.format(kwargs['formats'])
    cmd += '-scale {} '.format(kwargs['scale'])
    cmd += '-height {} '.format(kwargs['height'])
    cmd += '-width {} '.format(kwargs['width'])
    return cmd


def CollectTargetedPcrMetrics(**kwargs):
    cmd = '{} '.format(kwargs['python'])
    cmd += '{} '.format(kwargs['script'])
    cmd += 'CollectTargetedPcrMetrics '
    cmd += '-files {} '.format(kwargs['files'])
    cmd += '-outdir {} '.format(kwargs['outdir'])
    cmd += '-formats {} '.format(kwargs['formats'])
    cmd += '-scale {} '.format(kwargs['scale'])
    cmd += '-height {} '.format(kwargs['height'])
    cmd += '-width {} '.format(kwargs['width'])
    return cmd


def sample_correlation(**kwargs):
    cmd = '{} '.format(kwargs['python'])
    cmd += '{} '.format(kwargs['script'])
    cmd += '-data_file {} '.format(kwargs['data_file'])
    cmd += '-out_name {} '.format(kwargs['out_name'])
    cmd += '-corr_method {} '.format(kwargs['corr_method'])
    cmd += '-log_base {} '.format(kwargs['log_base'])
    cmd += '-log_additive {} '.format(kwargs['log_additive'])
    cmd += '-lower_exp_cutoff {} '.format(kwargs['lower_exp_cutoff'])
    cmd += '-row_sum_cutoff {} '.format(kwargs['row_sum_cutoff'])
    cmd += '-cv_cutoff {} '.format(kwargs['cv_cutoff'])
    cmd += '-height {} '.format(kwargs['height'])
    cmd += '-width {} '.format(kwargs['width'])
    cmd += '-sample_label_size {} '.format(kwargs['sample_label_size'])
    cmd += '-sample_label_angle {} '.format(kwargs['sample_label_angle'])
    cmd += '-color_scale {} '.format(kwargs['color_scale'])
    cmd += '--do_correlation_cluster '
    if kwargs['sample_group']:
        cmd += '-sample_group {} '.format(kwargs['sample_group'])
    return cmd


def make_slider(**kwargs):
    cmd = '{} '.format(kwargs['python'])
    cmd += '{} '.format(kwargs['script'])
    cmd += 'make_slider '
    cmd += '-images {} '.format(kwargs['images'])
    if kwargs['image_ids']:
        cmd += '-image_ids {} '.format(kwargs['image_ids'])
    cmd += '-image_desc "{}" '.format(kwargs['image_desc'])
    cmd += '-template {} '.format(kwargs['template'])
    cmd += '-out {} '.format(kwargs['out'])
    if kwargs['link_images'] != 'no' and kwargs['link_images']:
        cmd += '--link_images {} '
    return cmd
