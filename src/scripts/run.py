#!/usr/bin/env python3
import os

from cmi_pb_server.run import run


if __name__ == '__main__':
    #os.chdir("../..")
    run(
        "build/mro-tables.db",
        "src/table.tsv",
        cgi_path="/MRO/branches/demo-site/views/src/scripts/run.py",
        log_file="mro.log",
        title="MRO",
    )
