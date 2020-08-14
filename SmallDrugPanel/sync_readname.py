import gzip
import sys
import os


def sync_name(fq, fq2):
    """assign read id of fq to fq2"""
    if fq.endswith('.gz'):
        file = gzip.open
    else:
        file = open
    tmp_out = '__tmp.' + fq2.rsplit('.', 1)[1]
    with file(fq) as fr, file(fq2) as fr2, file(tmp_out, 'w') as fw:
        for i, (line, line2) in enumerate(zip(fr, fr2), start=1):
            if i % 4 == 1:
                fw.write(line.split()[0]+' '+line2.split(' ')[1])
            else:
                fw.write(line2)
    os.rename(tmp_out, fq2)


if __name__ == '__main__':
    args = sys.argv
    sync_name(args[1], args[2])
