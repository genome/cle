#!/usr/bin/env python3

import sqlite3

db = 'test.db'

con = sqlite3.connect(db)
cur = con.cursor()

cur.execute("select * from aml_clearance_somatic_variants")
for row in cur.fetchall():
    print(*row, sep='\t')

con.close()
