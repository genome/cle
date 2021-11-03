#!/usr/bin/env python3

from pyliftover import LiftOver
import sqlite3
import warnings
import fnmatch
import glob
import sys
import os
import re

if len(sys.argv) != 3:
    sys.exit("Need two arguments: upload report directory path and sqlite DB path")

tsv_dir = sys.argv[1]
db = sys.argv[2]

con = sqlite3.connect(db)
cur = con.cursor()

if os.path.exists(db) and os.path.getsize(db) > 0:
    print("AML_clearance varaint DB existing")
else:
    print("Create aml_clearance_somatic_variants table")
    cur.execute("create table aml_clearance_somatic_variants (id integer primary key, assay_type text, case_name text, chromosome text, position integer, reference text, variant text, variant_type text, transcript_name text, consequence text, symbol text, c_position text, amino_acid_change text, variant_callers text)")

for f in sorted(glob.glob(tsv_dir + "/*.tsv", recursive=True)):
    if fnmatch.fnmatch(os.path.basename(f), '*full_variant_report*'):
        p = re.compile("(.*).full_variant_report")
        result = p.search(os.path.basename(f))
        case_name = result.group(1)
        with open(f, "r") as new_fh:
            for line in new_fh:
                if not line.startswith("\"Manual"):
                    columns = line.rstrip().replace('"', "").split("\t")
                    if len(columns[4]) == len(columns[5]):
                        var_type = 'snp'
                    elif len(columns[4]) > len(columns[5]):
                        var_type = 'del'
                    elif len(columns[4]) < len(columns[5]):
                        var_type = 'ins'
                    else:
                        warnings.warn(case_name + ": " +columns[4] + " " + columns[5])
                    
                    if len(columns) == 35:
                        cur.execute("insert into aml_clearance_somatic_variants (assay_type, case_name, chromosome, position, reference, variant, variant_type, transcript_name, consequence, symbol, c_position, amino_acid_change, variant_callers) values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", ('CWL', case_name, columns[1], columns[2], columns[4], columns[5], var_type, columns[20], columns[17], columns[18], columns[21], columns[22], columns[6].replace('-', ',')))
                    elif len(columns) == 33:
                        cur.execute("insert into aml_clearance_somatic_variants (assay_type, case_name, chromosome, position, reference, variant, variant_type, transcript_name, consequence, symbol, c_position, amino_acid_change, variant_callers) values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", ('CWL', case_name, columns[1], columns[2], columns[4], columns[5], var_type, columns[18], columns[15], columns[16], columns[19], columns[20], columns[6].replace('-', ',')))
                            
        new_fh.close()
        print(case_name + " CWL done")
    else:  
        p = re.compile("(.*)_.*.xlsx.tsv")
        result = p.search(os.path.basename(f))
        case_name = result.group(1)
        with open(f, "r") as old_fh:
            for line in old_fh:
                if not line.startswith("\"Manual"):
                    columns = line.rstrip().replace('"', "").split("\t")
                    chrom = 'chr' + columns[1]
                    lo = LiftOver('hg19', 'hg38')
                    new_pos = lo.convert_coordinate(chrom, int(columns[2]))[0][1]
                    cur.execute("insert into aml_clearance_somatic_variants (assay_type, case_name, chromosome, position, reference, variant, variant_type, transcript_name, consequence, symbol, c_position, amino_acid_change, variant_callers) values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", ('GMS', case_name, chrom, new_pos, columns[4], columns[5], columns[6], columns[7], columns[8], columns[12], columns[10], columns[11], columns[33].replace(' ', '')))
        old_fh.close()
        print(case_name + " GMS done")

con.commit()
con.close()
