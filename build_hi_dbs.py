#!/usr/bin/env python3
### Import Modules ###
import sys
from Bio import SeqIO
import os
import re
import operator
import csv
import pprint as pp
import locale
import argparse
import urllib3
import json
from subprocess import *
#call[("module", "load", "Python/3.4")]
#call[("module", "load", "ncbi-blast+/2.2.29")]
#call[("export", "BLASTDB=/blast/db")]
encoding = locale.getdefaultlocale()[1]
SCRIPT_PATH = os.path.realpath(__file__)
DIR_PATH = os.path.dirname(SCRIPT_PATH)
CUSTOM_DB = os.path.join(DIR_PATH,"custom_allele_sets")
OUTPUT_DIR = ""
genes_to_get = ["HAEM1876","HAEM1877","HAEM1878","HAEM1879","HAEM1880","HAEM1881","HAEM1882","HAEM1883","HAEM1873","HAEM1874","HAEM1875",
"HAEM1147","HAEM1148","HAEM1149","HAEM1150","HAEM1153","HAEM1156","HAEM1151","HAEM1154","HAEM1155","HAEM1157","HAEM1158",
"HAEM1160","HAEM1161","HAEM1163","HAEM1164","HAEM1165"]
def set_output(output):
	global OUTPUT_DIR
	OUTPUT_DIR = os.path.join(DIR_PATH,output)
	if os.path.isdir(output):
		print("Output Directory",output,"aleady exists, not creating")
	else:
		os.system("mkdir {}".format(OUTPUT_DIR))
		os.system("cp -r {}/hinfluenzae {}".format(CUSTOM_DB,OUTPUT_DIR))
		print("Created Output Directory",output)

def make_blast_db(allele_fasta,allele_name,mol_type,allele_db_path):
	call(["makeblastdb","-in", allele_fasta,"-input_type","fasta","-title",allele_name,"-dbtype",mol_type,"-out",os.path.join(allele_db_path,allele_name)],shell=False)
	
def main():		
	output = "hinfluenzae_capsule_DB"
	set_output(output)	
	http = urllib3.PoolManager()
	dna = ["A","G","C","T"]
	db_names = ["hinfluenzae"] #PubMLST Species DB names
	for db in db_names:
		all_loci_dna = ""
		for locus in genes_to_get:
			locus_url = "http://rest.pubmlst.org/db/pubmlst_hinfluenzae_seqdef/loci/{}".format(locus)
			fasta_request = http.request('GET',locus_url+"/alleles_fasta")
			fasta_data = fasta_request.data.decode('utf-8')		
			allele_name = locus
			if "{" in fasta_data:
				continue # no results found
			allele_db_path = os.path.join(OUTPUT_DIR,"{}".format(db),allele_name.replace("'","#"))
			if os.path.exists(allele_db_path):
				new_allele_count = int(fasta_data.count(">"))
				with open(os.path.join(allele_db_path,"{}.fasta".format(allele_name)),"r") as f:
					existing_file = f.read()
				existing_file_count = int(existing_file.count(">"))
				if new_allele_count == existing_file_count:
					print("no changes detected for {}, skipping".format(allele_name))
					continue
				else:
					os.system("rm -r {}".format(allele_db_path)) #overwrite system
					os.system("mkdir {}".format(allele_db_path))
			else:
				os.system("mkdir {}".format(allele_db_path))
			allele_fasta = os.path.join(allele_db_path,"{}.fasta".format(allele_name))		
			first = True
			with open(os.path.join(OUTPUT_DIR,"hinfluenzae",allele_fasta),"w") as f:
				f.write(fasta_data)
			first_rec = True
			for seq_record in SeqIO.parse(allele_fasta,"fasta"):
				if first_rec:
					mol_type = "nucl"
					seq = seq_record.seq
					first_rec = False
					for letter in seq:
						if letter not in dna:
							mol_type = "prot"
				else:
					break	
			call(["makeblastdb","-in", allele_fasta,"-input_type","fasta","-title",allele_name,"-dbtype",mol_type,"-out",os.path.join(allele_db_path,allele_name)],shell=False)					
	
		
if __name__ == "__main__":
	main()
				
