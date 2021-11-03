#!/usr/bin/env python3

import sqlite3
import sys

if len(sys.argv) != 2:
    sys.exit("Need provide sqlite3 DB path")

db = sys.argv[1] 

con = sqlite3.connect(db)
cur = con.cursor()

cur.execute("select * from aml_clearance_somatic_variants")
for row in cur.fetchall():
    print(*row, sep='\t')

con.close()
