#!/usr/bin/env python3.9
import os

from cmi_pb_server.run import run


if __name__ == '__main__':
    os.chdir("../..")
    run(
        "build/mro-tables.db",
        "src/table.tsv",
        cgi_path="/MRO/branches/demo-site/views/src/scripts/run.py",
        log_file="mro.log",
      	hide_index=True,
        max_children=100,
        term_index="index",
        title="MRO",
    )
