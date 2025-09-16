from cns.pipelines import exists
from cns.utils.files import exists


import time
from os.path import exists


def save_time(action, out_file, runtime, start, threads):
        timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(start))
        filepath = "../docs/times.tsv"
        write_mode = "a" if exists(filepath) else "w"
        with open(filepath, write_mode) as f:
            f.write(f"{timestamp}\t{action}\t{threads}\t{out_file}\t{runtime}\n")


def parse_cncols(cncols):
    if cncols != None:
        cncols = cncols.split(",")
        if len(cncols) > 2:
            raise ValueError("Only one or two columns can be specified.")
    return cncols