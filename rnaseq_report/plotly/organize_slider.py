import os
import os.path as path
from os.path import join
import shutil
from jinja2 import Template


def make_report(report_dir):
    index_html_path = join(report_dir, 'index.html')
    html_utils_dir = os.path.dirname(os.path.dirname(__file__))
    shutil.copytree(join(html_utils_dir, 'html.utils'), join(report_dir, 'html.utils'))
    index_template = join(html_utils_dir, 'templates', 'index.jinja2')
    result_dict = dict()
    result_dict['质量评估'] = {
        "GeneBodyCoverage": join(report_dir, 'Summary', 'GeneBodyCoverage.html'),
        "ReadDistribution": join(report_dir, 'Summary', 'ReadDistribution.html'),
        "ChrReadDistribution": join(report_dir, 'Summary', 'ChrReadDistribution.html'),
        "FragmentSize": join(report_dir, 'Summary', 'FragmentSize.html'),
        "ReadDuplication": join(report_dir, 'Summary', 'ReadDuplication.html'),
        "InnerDistance": join(report_dir, 'Summary', 'InnerDistance.html'),
        "ExpSaturation": join(report_dir, 'Summary', 'ExpSaturation.html'),
        "InsertSize": join(report_dir, 'Summary', 'CollectInsertSizeMetrics.html'),
    }
    result_dict['表达量分析'] = {
        "ExpressionDistribution": join(report_dir, 'Summary', 'ExpressionDistribution.html'),
        "SampleCorrelationCluster": join(report_dir, 'Summary', 'CorrelationCluster.html'),
        "PrincipalComponentAnalysis": join(report_dir, 'Summary', 'PCA.html'),
    }
    result_dict['比对结果统计表'] = {
        "AlignmentSummary": join(report_dir, 'Summary', 'CollectAlignmentSummaryMetrics.html'),
        "RnaSeqMetrics": join(report_dir, 'Summary', 'CollectRnaSeqMetrics.html'),
        "TargetedSummary": join(report_dir, 'Summary', 'TargetedSummaryMetrics.html'),
    }


    # use relative path, this is important
    for section in result_dict:
        for each in result_dict[section]:
            # print(section)
            tmp_path = path.abspath(result_dict[section][each])
            result_dict[section][each] = path.relpath(tmp_path, start=path.abspath(report_dir))
    # format html
    template = Template(open(index_template).read())
    content = template.render(result=result_dict)
    with open(index_html_path, 'w', encoding='utf-8') as f:
        f.write(content)


def make_slider(images:list, image_ids:list=None, image_desc:list=None, template="templates/slide.jinja2",
                out='slider.html', link_images=False):
    out_dir = os.path.dirname(out)
    if not out_dir:
        out_dir = os.getcwd()
    else:
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
    if not images:
        return

    if link_images:
        out_dir = os.path.join(out_dir, os.path.basename(out)[:-5]+'.images')
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        tmp_list = list()
        for each in images:
            target = os.path.join(out_dir, os.path.basename(each))
            if os.path.exists(target):
                continue
            os.symlink(each, target)
            tmp_list.append(target)
        images = tmp_list

    if not image_ids:
        image_ids = [os.path.basename(x).split('.')[0] for x in images]
    if len(image_ids) < len(images):
        raise Exception('number of image ids is less than the number of images')

    if not image_desc:
        image_desc = ['']*len(images)
    if len(image_desc) == 1:
        image_desc = image_desc*len(images)
    else:
        if len(image_desc) < len(images):
            raise Exception('number of image desc is less than the number of images')

    # env = Environment(loader=FileSystemLoader("templates"))
    # template = env.get_template("slide.jinja2")
    template = Template(open(template).read())
    img_info_dict = dict()
    img_cls = ["slide"]*len(images)
    img_cls[0] = "slide slide--current"
    for img, img_id, desc, tag in zip(images, image_ids, image_desc, img_cls):
        img_info_dict[img_id] = dict()
        img_info_dict[img_id]['path'] = os.path.relpath(os.path.abspath(img), start=os.path.abspath(out_dir))
        img_info_dict[img_id]['cls'] = tag
        img_info_dict[img_id]['description'] = desc
    # print(img_info_dict)
    content = template.render(img_info_dict=img_info_dict)
    with open(out, 'w', encoding='utf-8') as f:
        f.write(content)
    return out


if __name__ == '__main__':
    class Func2Command(object):
        def __init__(self, callable_dict):
            self.call(callable_dict)

        @staticmethod
        def introduce_command(func):
            import argparse
            import inspect
            import json
            import time
            if isinstance(func, type):
                description = func.__init__.__doc__
            else:
                description = func.__doc__
            if description:
                _ = [print(x.strip()) for x in description.split('\n') if x.strip()]
                parser = argparse.ArgumentParser(add_help=False)
            else:
                parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description=description)
            func_args = inspect.getfullargspec(func)
            arg_names = func_args.args
            if not arg_names:
                func()
                return
            arg_defaults = func_args.defaults
            if not arg_defaults:
                arg_defaults = []
            arg_defaults = ['None']*(len(arg_names) - len(arg_defaults)) + list(arg_defaults)
            sig = inspect.signature(func)
            for arg, value in zip(arg_names, arg_defaults):
                arg_type = sig.parameters[arg].annotation
                if arg == 'self':
                    continue
                if value == 'None':
                    if arg_type in [list, tuple, set]:
                        parser.add_argument('-' + arg, nargs='+', required=True, metavar=arg)
                    elif arg_type in [int, str, float]:
                        parser.add_argument('-' + arg, type=arg_type, required=True, metavar=arg)
                    else:
                        parser.add_argument('-'+arg, required=True, metavar=arg)
                elif type(value) == bool:
                    if value:
                        parser.add_argument('--'+arg, action="store_false", help='default: True')
                    else:
                        parser.add_argument('--'+arg, action="store_true", help='default: False')
                elif value is None:
                    if arg_type in [list, tuple, set]:
                        parser.add_argument('-' + arg, nargs='+', required=False, metavar=arg)
                    elif arg_type in [int, str, float]:
                        parser.add_argument('-' + arg, type=arg_type, required=False, metavar=arg)
                    else:
                        parser.add_argument('-'+arg, required=False, metavar='Default:' + str(value))
                else:
                    if arg_type in [list, tuple, set] or (type(value) in [list, tuple, set]):
                        default_value = ' '.join(str(x) for x in value)
                        if type(value) in [list, tuple]:
                            one_value = value[0]
                        else:
                            one_value = value.pop()
                        parser.add_argument('-' + arg, default=value, nargs='+', type=type(one_value),
                                            metavar='Default:'+default_value, )
                    else:
                        parser.add_argument('-' + arg, default=value, type=type(value),
                                            metavar='Default:' + str(value), )
            if func_args.varargs is not None:
                print("warning: *varargs is not supported, and will be neglected! ")
            if func_args.varkw is not None:
                print("warning: **keywords args is not supported, and will be neglected! ")
            args = parser.parse_args().__dict__
            try:
                with open("Argument_detail.json", 'w') as f:
                    json.dump(args, f, indent=2, sort_keys=True)
            except IOError:
                print('Current Directory is not writable, thus argument log is not written !')
            start = time.time()
            func(**args)
            print("total time: {}s".format(time.time() - start))

        def call(self, callable_dict):
            import sys
            excludes = ['introduce_command', 'Func2Command']
            _ = [callable_dict.pop(x) for x in excludes if x in callable_dict]
            if len(callable_dict) >= 2:
                if len(sys.argv) <= 1:
                    print("The tool has the following sub-commands: ")
                    _ = [print(x) for x in callable_dict]
                    return
                sub_cmd = sys.argv[1]
                sys.argv.remove(sub_cmd)

                if sub_cmd in callable_dict:
                    self.introduce_command(callable_dict[sub_cmd])
                else:
                    print('sub-command: {} is not defined'.format(sub_cmd))
            elif len(callable_dict) == 1:
                self.introduce_command(callable_dict.pop(list(callable_dict.keys())[0]))

    callable_dict = {x: y for x, y in locals().items() if callable(y)}
    _ = [callable_dict.pop(x) for x in {'Func2Command'} if x in callable_dict]
    Func2Command(callable_dict)
