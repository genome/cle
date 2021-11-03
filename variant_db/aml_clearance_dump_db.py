#!/usr/bin/env python3

import sqlite3

db = '/storage1/fs1/gtac-mgi/Active/CLE/validation/cle_validation/CLE_variant_database/sqlite_variant_DB/aml_clearance_variant.db'

con = sqlite3.connect(db)
cur = con.cursor()

cur.execute("select * from aml_clearance_somatic_variants")
for row in cur.fetchall():
    print(*row, sep='\t')

con.close()
