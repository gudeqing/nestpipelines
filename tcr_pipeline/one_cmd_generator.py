# coding = utf-8
import os


def CalcBasicStats(**kwargs):
    cmd = '{} '.format(kwargs['vdjtools'])
    cmd += '{} '.format(kwargs['sub_cmd'])
    cmd += '-m {} '.format(kwargs['metadata'])
    cmd += '{} '.format(kwargs['prefix'])
    return cmd


def CalcSegmentUsage(**kwargs):
    cmd = '{} '.format(kwargs['vdjtools'])
    cmd += '{} '.format(kwargs['sub_cmd'])
    cmd += '-m {} '.format(kwargs['metadata'])
    if kwargs['plot'] == 'yes':
        cmd += '-p '
        cmd += '-f {} '.format(kwargs['factor_field'])
        if kwargs['factor_is_numeric'] == 'yes':
            cmd += '-n '.format(kwargs['factor_is_numeric'])
        cmd += '-l {} '.format(kwargs['label_field'])
    cmd += '{} '.format(kwargs['prefix'])
    return cmd


def CalcSpectratype(**kwargs):
    cmd = '{} '.format(kwargs['vdjtools'])
    cmd += '{} '.format(kwargs['sub_cmd'])
    cmd += '-m {} '.format(kwargs['metadata'])
    cmd += '{} '.format(kwargs['prefix'])
    return cmd


def PlotFancySpectratype(**kwargs):
    cmd = '{} '.format(kwargs['vdjtools'])
    cmd += '{} '.format(kwargs['sub_cmd'])
    cmd += '-t {} '.format(kwargs['top'])
    cmd += '{} '.format(kwargs['sample'])
    cmd += '{} '.format(kwargs['prefix'])
    return cmd


def PlotFancyVJUsage(**kwargs):
    cmd = '{} '.format(kwargs['vdjtools'])
    cmd += '{} '.format(kwargs['sub_cmd'])
    cmd += '{} '.format(kwargs['sample'])
    cmd += '{} '.format(kwargs['prefix'])
    return cmd


def PlotSpectratypeV(**kwargs):
    cmd = '{} '.format(kwargs['vdjtools'])
    cmd += '{} '.format(kwargs['sub_cmd'])
    cmd += '-t {} '.format(kwargs['top'])
    cmd += '{} '.format(kwargs['sample'])
    cmd += '{} '.format(kwargs['prefix'])
    return cmd


def PlotQuantileStats(**kwargs):
    cmd = '{} '.format(kwargs['vdjtools'])
    cmd += '{} '.format(kwargs['sub_cmd'])
    cmd += '-t {} '.format(kwargs['top'])
    cmd += '{} '.format(kwargs['sample'])
    cmd += '{} '.format(kwargs['prefix'])
    return cmd


def RarefactionPlot(**kwargs):
    cmd = '{} '.format(kwargs['vdjtools'])
    cmd += '{} '.format(kwargs['sub_cmd'])
    cmd += '-i {} '.format(kwargs['intersect-type'])
    cmd += '-s {} '.format(kwargs['steps'])
    cmd += '-f {} '.format(kwargs['factor_field'])
    if kwargs['factor_is_numeric'] == 'yes':
        cmd += '-n '.format(kwargs['factor_is_numeric'])
    cmd += '-l {} '.format(kwargs['label_field'])
    if kwargs['label-exact'] == 'yes':
        cmd += '--label-exact '
    cmd += '-m {} '.format(kwargs['metadata'])
    cmd += '{} '.format(kwargs['prefix'])
    return cmd


def CalcDiversityStats(**kwargs):
    cmd = '{} '.format(kwargs['vdjtools'])
    cmd += '{} '.format(kwargs['sub_cmd'])
    cmd += '-m {} '.format(kwargs['metadata'])
    cmd += '-i {} '.format(kwargs['intersect-type'])
    if kwargs['downsample_to'] != 'smallest':
        cmd += '-x {} '.format(kwargs['downsample_to'])
    if kwargs['extrapolate_to'] != 'largest':
        cmd += '-X {} '.format(kwargs['extrapolate_to'])
    cmd += '{} '.format(kwargs['prefix'])
    return cmd


def CalcPairwiseDistances(**kwargs):
    cmd = '{} '.format(kwargs['vdjtools'])
    cmd += '{} '.format(kwargs['sub_cmd'])
    cmd += '-m {} '.format(kwargs['metadata'])
    cmd += '-i {} '.format(kwargs['intersect-type'])
    if kwargs['plot'] == 'yes':
        cmd += '-p '
    cmd += '{} '.format(kwargs['prefix'])
    return cmd


def ClusterSamples(**kwargs):
    cmd = '{} '.format(kwargs['vdjtools'])
    cmd += '{} '.format(kwargs['sub_cmd'])
    cmd += '-i {} '.format(kwargs['intersect-type'])
    cmd += '-e {} '.format(kwargs['measure'])
    if kwargs['plot'] == 'yes':
        cmd += '-p '
        cmd += '-f {} '.format(kwargs['factor_field'])
        if kwargs['factor_is_numeric'] == 'yes':
            cmd += '-n '.format(kwargs['factor_is_numeric'])
        cmd += '-l {} '.format(kwargs['label_field'])
    cmd += '{} '.format(kwargs['input_prefix'])
    cmd += '{} '.format(kwargs['out_prefix'])
    return cmd

