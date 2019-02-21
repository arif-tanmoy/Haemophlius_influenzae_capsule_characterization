#!/usr/bin/env python3.4
### Import Modules ###
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import os
import re
import json
from collections import OrderedDict
import operator
import csv
import pprint as pp
import string
import locale
import time
import tempfile
import argparse
import urllib3
import math
from multiprocessing import Pool, Process, Queue
from subprocess import *
#call[("module", "load", "Python/3.4")]
#call[("module", "load", "ncbi-blast+/2.2.29")]
#call[("export", "BLASTDB=/blast/db")]

### SET UP DB connections and dictionaries ###
encoding = locale.getdefaultlocale()[1]
SCRIPT_PATH = os.path.realpath(__file__)
DIR_PATH = os.path.dirname(os.path.dirname(SCRIPT_PATH))
ABS_PATH = os.path.dirname(SCRIPT_PATH)
OUTPUT_DIR = ""
http = urllib3.PoolManager()
PUBMLST_DB = os.path.join(ABS_PATH,"hinfluenzae_capsule_DB")
fasta_extensions = [".fasta",".fna",".fas"]
DNA = ["A","T","C","G"]
gene_names ={"HAEM1876":"ecs1","HAEM1877":"ecs2","HAEM1878":"ecs3","HAEM1879":"ecs4","HAEM1880":"ecs5","HAEM1881":"ecs6","HAEM1882":"ecs7","HAEM1883":"ecs8","HAEM1873":"fcs1","HAEM1874":"fcs2","HAEM1875":"fcs3",
"HAEM1147":"hcsB","HAEM1148":"hcsA","HAEM1149":"bcs4","HAEM1150":"bcs3","HAEM1152":"bcs1","HAEM1153":"bexD","HAEM1156":"bexA","HAEM1151":"bcs2","HAEM1154":"bexC","HAEM1155":"bexB","HAEM1157":"hcsB","HAEM1158":"hcsA",
"HAEM1160":"bcs3","HAEM1161":"bcs2","HAEM1162":"bcs1","HAEM1163":"bexD","HAEM1164":"bexC","HAEM1165":"bexB"}
#allele_exceptions = ["Insertion_Element","Insertion_Elements"]

serotypes = {
    "f":{"essential":["fcs1","fcs2","fcs3","bexA","bexB","bexC","bexD","hcsA","hcsB"],"nonessential":[]},
    "a":{"essential":["acs1","acs2","acs3","acs4","bexA","bexB","bexC","bexD","hcsA","hcsB"],"nonessential":[]},
    "c":{"essential":["ccs1","ccs2","ccs3","ccs4","bexA","bexB","bexC","bexD","hcsA","hcsB"],"nonessential":[]},
    "b":{"essential":["bcs1","bcs2","bcs3","bcs4","bexA","bexB","bexC","bexD","hcsA","hcsB"],"nonessential":[]},
    "d":{"essential":["dcs1","dcs2","dcs3","dcs4","dcs5","bexA","bexB","bexC","bexD","hcsA","hcsB"],"nonessential":[]},
    "e":{"essential":["ecs1","ecs2","ecs3","ecs4","ecs5","ecs6","ecs7","ecs8","bexA","bexB","bexC","bexD","hcsA","hcsB"],"nonessential":[]},    
}
ST_unique = {"fcs1":"f","fcs2":"f","fcs3":"f",
            "acs1":"a","acs2":"a","acs3":"a","acs4":"a",
            "ccs1":"c","ccs2":"c","ccs3":"c","ccs4":"c",
            "bcs1":"b","bcs2":"b","bcs3":"b","bcs4":"b",
            "dcs1":"d","dcs2":"d","dcs3":"d","dcs4":"d","dcs5":"d",
            "ecs1":"e","ecs2":"e","ecs3":"e","ecs4":"e","ecs5":"e","ecs6":"e","ecs7":"e","ecs8":"e"}
			

species_dict = {"Haemophilus influenzae":"hinfluenzae"}

### Check if output directory exists, if not, create, if yes, set global output to dir ###
def set_output(output):
	global OUTPUT_DIR
	OUTPUT_DIR = output
	if os.path.isdir(output):
		print("Output Directory",output,"already exists, not creating")
	else:
		os.system("mkdir {}".format(output))
		print("Created Output Directory",output)
		
		
def create_gff(results_dict,seq_dict):
	for in_file in results_dict:
		text = ""
		header_info = ""		
		fasta_text = ""
		out_path = os.path.join(OUTPUT_DIR,"gff")
		with open(os.path.join(out_path,"{}.gff".format(in_file)),"w") as f:		
			for contig in results_dict[in_file]["contigs"]:
				if contig == "species":
					continue					
				count=0
				for allele in results_dict[in_file]["contigs"][contig]:
					if results_dict[in_file]["contigs"][contig][allele]:
						for hit in results_dict[in_file]["contigs"][contig][allele]:
							seq_dict[in_file]["contigs"][contig]["alleles"][count] = hit 
							count+=1
			for contig,_ in sorted(seq_dict[in_file]["contigs"].items(),key=lambda x: int(x[1]["length"]),reverse=True):
				sequence = seq_dict[in_file]["contigs"][contig]["seq"]
				seq_length = seq_dict[in_file]["contigs"][contig]["length"]				
				if seq_dict[in_file]["contigs"][contig]["alleles"]:
					for hit,_ in sorted(seq_dict[in_file]["contigs"][contig]["alleles"].items(),key= lambda x: int(x[1]["qstart"])):			
						allele_name = seq_dict[in_file]["contigs"][contig]["alleles"][hit]["allele_name"]
						if "Insertion" in allele_name:
							scheme_type = "Insertion_Element"
						else:
							scheme_type = "Capsule Locus"
						allele = seq_dict[in_file]["contigs"][contig]["alleles"][hit]["allele"]
						allele_id = str(seq_dict[in_file]["contigs"][contig]["alleles"][hit]["allele_id"])
						start = str(seq_dict[in_file]["contigs"][contig]["alleles"][hit]["qstart"])						
						stop = str(seq_dict[in_file]["contigs"][contig]["alleles"][hit]["qend"])						
						score = str(seq_dict[in_file]["contigs"][contig]["alleles"][hit]["score"])						
						strand = str(seq_dict[in_file]["contigs"][contig]["alleles"][hit]["strand"])						
						raw_flags = seq_dict[in_file]["contigs"][contig]["alleles"][hit]["flags"]
						if len(raw_flags) == 1:
							flags = raw_flags[0]
						elif len(raw_flags) > 1:
							flags = ",".join(raw_flags)
						else:
							flags = "N/A"
						if "N/A" in flags and "internal" in flags:
							flags.replace("N/A","")
						annotations = seq_dict[in_file]["contigs"][contig]["alleles"][hit]["annotations"]
						new = seq_dict[in_file]["contigs"][contig]["alleles"][hit]["new"]
						edge_match = seq_dict[in_file]["contigs"][contig]["alleles"][hit]["edge"]
						identity = seq_dict[in_file]["contigs"][contig]["alleles"][hit]["identity"]
						a_length = seq_dict[in_file]["contigs"][contig]["alleles"][hit]["align_length"]
						s_length = seq_dict[in_file]["contigs"][contig]["alleles"][hit]["subject_length"]
						region_type = seq_dict[in_file]["contigs"][contig]["alleles"][hit]["region_type"]
						cov = float(seq_dict[in_file]["contigs"][contig]["alleles"][hit]["cov"])
						cov = round(cov,2)
						cov = cov*100.0
						if new and edge_match:
							text += contig +"\t" + "pubMLST" + "\t" + "{}".format(region_type) + "\t" + start + "\t" + stop + "\t" + "." + "\t" + strand + "\t" + "0" + "\t" + "ID={}_{};gene={};allele_id={};inference={};locus_tag={}_{};flags={};product={};scheme={}\n".format(contig,count,allele_name,allele_id,"pubmlst",contig,count,flags,annotations,scheme_type)							
						elif new and not edge_match:
							text += contig +"\t" + "pubMLST" + "\t" + "{}".format(region_type) + "\t" + start + "\t" + stop + "\t" + "." + "\t" + strand + "\t" + "0" + "\t" + "ID={}_{};gene={};allele_id={};inference={};locus_tag={}_{};flags={};product={};scheme={}\n".format(contig,count,allele_name,allele_id,"pubmlst",contig,count,flags,annotations,scheme_type)
						elif edge_match and not new:
							text += contig +"\t" + "pubMLST" + "\t" + "{}".format(region_type) + "\t" + start + "\t" + stop + "\t" + "." + "\t" + strand + "\t" + "0" + "\t" + "ID={}_{};gene={};allele_id={}_edge_match_identity({})_coverage({}%);inference={};locus_tag={}_{};flags={};product={};scheme={}\n".format(contig,count,allele_name,allele_id,identity,cov,"pubmlst",contig,count,flags,annotations,scheme_type)
						elif cov != 100:
							text += contig +"\t" + "pubMLST" + "\t" + "{}".format(region_type) + "\t" + start + "\t" + stop + "\t" + "." + "\t" + strand + "\t" + "0" + "\t" + "ID={}_{};gene={};allele_id={}_partial_cov_({}%);inference={};locus_tag={}_{};flags={};product={};scheme={}\n".format(contig,count,allele_name,allele_id,cov,"pubmlst",contig,count,flags,annotations,scheme_type)
						else:
							text += contig +"\t" + "pubMLST" + "\t" + "{}".format(region_type) + "\t" + start + "\t" + stop + "\t" + "." + "\t" + strand + "\t" + "0" + "\t" + "ID={}_{};gene={};allele_id={};inference={};locus_tag={}_{};flags={};product={};scheme={}\n".format(contig,count,allele_name,allele_id,"pubmlst",contig,count,flags,annotations,scheme_type)
						count+=1
				header_info += "##sequence-region {} 1 {}\n".format(contig,str(seq_length))
				fasta_text += ">{}\n".format(contig)
				fasta_text += "{}\n".format(sequence)
			f.write("##gff-version 3\n")
			f.write(header_info)
			f.write(text)
			f.write("##FASTA\n")
			f.write(fasta_text)

### Parse blast results, filter for best allele
def analyze_blast(results,input_file,seq_dict):
	results_dict = {"species":seq_dict[input_file]["species"],"contigs":{}}
	final_dict = {}
	final_dict[input_file] = {"species":seq_dict[input_file]["species"],"contigs":{}}
	#final_dict[input_file]["species"] = seq_dict[input_file]["species"]
	for item in results:
		for line in item:
			if line.decode(encoding).strip() != "":
				pident = line.decode(encoding).split("\t")[2]
				ident = float(pident)
				#only keep hits with 90% or greater identity
				if ident < 90.0:
					continue				
				query_name = line.decode(encoding).split("\t")[0]
				contig = query_name
				if contig not in results_dict["contigs"]:
					results_dict["contigs"][contig] = {}
				subject = line.decode(encoding).split("\t")[1]
				allele_num = subject.split("_")[-1]
				subject_name = re.sub(r'(_{})$'.format(allele_num),r'',subject)
				if subject_name not in results_dict["contigs"][contig]:
					results_dict["contigs"][contig][subject_name] = {}
				if subject_name in gene_names:
					allele_name = gene_names[subject_name]
				else:
					allele_name = subject_name

				qcovs = line.decode(encoding).split("\t")[3]
				qstart = line.decode(encoding).split("\t")[6]
				qend = line.decode(encoding).split("\t")[7]
				score = float(line.decode(encoding).split("\t")[10])
				send = int(line.decode(encoding).split("\t")[9])
				sstart = int(line.decode(encoding).split("\t")[8])
				a_length = int(line.decode(encoding).split("\t")[11])
				s_len = int(line.decode(encoding).split("\t")[13])
				q_len = int(line.decode(encoding).split("\t")[12])
				qseq = line.decode(encoding).split("\t")[15]
				if int(qstart) > int(qend):
					qstart = line.decode(encoding).split("\t")[7]
					qend = line.decode(encoding).split("\t")[6]
				coordinates = qstart+"*"+qend
				if coordinates not in results_dict["contigs"][contig][subject_name]:
					results_dict["contigs"][contig][subject_name][coordinates] = {}	
				if allele_num not in results_dict["contigs"][contig][subject_name][coordinates]:
					results_dict["contigs"][contig][subject_name][coordinates][allele_num] = []					
				try:
					sframe = str(line.decode(encoding).split("\t")[14])
					if "-" in sframe:
						strand = "-"
					else:
						strand = "+"
				except:
					strand = "+"
				hit_info = {"allele_name":allele_name,"new":False,"length":a_length,"q_length":q_len,"identity":pident,"score":score,"qstart":qstart,"qend":qend,"contig":query_name,"strand":strand,"allele_id":allele_num,"s_length":s_len,"a_length":a_length,"qseq":qseq}
				results_dict["contigs"][contig][subject_name][coordinates][allele_num].append(hit_info)		
	complete_cov = {}
	for contig in results_dict["contigs"]:
		final_dict[input_file]["contigs"][contig] = {}
		false_hits = []
		contig_stop = int(seq_dict[input_file]["contigs"][contig]["length"])		
		for allele in sorted(results_dict["contigs"][contig]):
			final_dict[input_file]["contigs"][contig][allele] = []
			for coordinates in sorted(results_dict["contigs"][contig][allele]):
				current_score = 0
				current_identity = 0.0
				current_cov = 0.0
				current_factor = 0.0				
				for allele_id in sorted(results_dict["contigs"][contig][allele][coordinates]):
					for a_num in results_dict["contigs"][contig][allele][coordinates][allele_id]:
						score = int(a_num["score"])
						identity = float(a_num["identity"])
						align_length = float(a_num["a_length"])
						subject_length = float(a_num["s_length"])
						query_length = float(a_num["q_length"])
						allele_name = a_num["allele_name"]
						start = int(a_num["qstart"])
						end = int(a_num["qend"])
						start_stop = str(start)+"*"+str(end)											
						cov = float(align_length/subject_length)
						false_pos = False
						if cov > 1.0:
							cov = 1.0
						factor = float(cov * identity)
						ident_cutoff = 0
						cov_cutoff = .70
						# factor is custom scoring (cov * identity) - compares each hit's factor vs current best factor within each allele set that hit the same region				
						if factor >= current_factor:					
							current_score = score
							current_identity = identity
							current_cov = cov
							current_factor = factor
							if end == contig_stop and identity >= 95.0:
								edge_match = True						
							elif start == 1 and identity >= 95.0 and cov < 1.0:
								edge_match = True
							else:
								edge_match = False														
							if a_num["strand"] is "+":
								qseq = seq_dict[input_file]["contigs"][contig]["seq"][start-1:end]
							if a_num["strand"] is "-":												
								qseq = seq_dict[input_file]["contigs"][contig]["seq"][start-1:end].reverse_complement()						
							if (current_identity < ident_cutoff or current_cov < cov_cutoff) and not edge_match:
								#~ if allele in allele_exceptions:
									#~ if current_cov > .30:
										#~ false_pos = False
									#~ else:
										#~ false_pos = True
								#~ else:
									#~ false_pos = True
								if edge_match:
									if current_cov < 0.10:
										false_pos = True							
							else:
								false_pos = False
							if current_cov > cov_cutoff:
								if allele not in complete_cov:
									complete_cov[allele] = True							
							remove_these = []			
							best_hit = True				
							#Checks current best allele against previously seen best alleles for this locus, and decides whether or not it remains the best
							if final_dict[input_file]["contigs"][contig][allele]:
								for seen_hits in final_dict[input_file]["contigs"][contig][allele]:
									r_start = int(seen_hits["qstart"])
									r_stop = int(seen_hits["qend"])
									r_score = int(seen_hits["score"])
									r_factor = float(seen_hits["factor"])
									if start <= r_start and end <= r_stop and end >= r_start:
										if current_factor >= r_factor:
											remove_these.append(seen_hits)
										else:
											best_hit = False																							
									elif start >= r_start and start <= r_stop and end >= r_stop:
										if current_factor >= r_factor:
											remove_these.append(seen_hits)	
										else:
											best_hit = False																																					
									elif start >= r_start and start <= r_stop and end >= r_start and end <= r_stop:
										if current_factor >= r_factor:
											remove_these.append(seen_hits)
										else:
											best_hit = False																																						
									elif start <= r_start and start <= r_stop and end >= r_start and end >= r_stop:
										if current_factor >= r_factor:
											remove_these.append(seen_hits)		
										else:
											best_hit = False	
							if not best_hit:
								false_pos = True																										
							for hits in remove_these:
								final_dict[input_file]["contigs"][contig][allele].remove(hits)	
							if "Insertion_Element" in allele:
								region_type = "ISE"
							else:
								region_type = "CDS"
							#Stores metadata in final results dict					
							hit = {"region_type":region_type,"false_pos":false_pos,"allele":allele,"factor":factor,"allele_name":allele_name,"edge":edge_match,"new":False,"qseq":str(qseq),"cov":cov,"subject_length": subject_length, "align_length":align_length,"allele_id":a_num["allele_id"],"identity":current_identity,"score":current_score,"contig":a_num["contig"],"strand":a_num["strand"],"qstart":a_num["qstart"],"qend":a_num["qend"]}
							final_dict[input_file]["contigs"][contig][allele].append(hit)	

																	
	for allele in complete_cov:
		#~ if allele in allele_exceptions:
			#~ continue
		for contig in final_dict[input_file]["contigs"]:
			if contig == "species":
				continue			
			if allele in final_dict[input_file]["contigs"][contig]:
				for hit in final_dict[input_file]["contigs"][contig][allele]:
					if hit["cov"] < cov_cutoff:
						hit["false_pos"] = True		
	return final_dict


def sql_find_allele(edited_allele):
	c.execute("select * from Gene where name like'{}'".format(edited_allele))
	a_result = c.fetchone()
	return a_result
	
def sql_get_annotations(allele_db_id):
	c.execute("select * from Annotation where gene_id='{}'".format(allele_db_id))
	results = c.fetchall()
	return results
	
## Step 2 - finds best hit per coordinate range	and pulls pubmlst info for best hits
def analyze_results(results_dict,threads):
	print("Parsing BLAST data to identify alleles to extract from pubMLST")
	### Gets id-specific allele information such as sequence and flag information from pubMLST for each result ###
	alleles_to_grab = {}
	allele_info = {}
	allele_count = 0
	allele_annotations = {}	
	for in_file in sorted(results_dict):
		if "species" in results_dict[in_file]:
			species = results_dict[in_file]["species"]
		else:
			results_dict[in_file]["species"] = "hinfluenzae"
			species = results_dict[in_file]["species"]			
		if species not in alleles_to_grab:
			alleles_to_grab[species] = {}			
		false_pos = {}
		edge_cases = {}
		for contig in sorted(results_dict[in_file]["contigs"]):
			if contig == "species": #Legacy compatability
				continue
			edge_cases[contig] = {}
			seen_regions = {}				
			if contig not in false_pos:
				false_pos[contig] = {}
			for allele in sorted(results_dict[in_file]["contigs"][contig]):
				if allele not in false_pos[contig]:
					false_pos[contig][allele] = []
				if results_dict[in_file]["contigs"][contig][allele]:				
					for hit in results_dict[in_file]["contigs"][contig][allele]:					
						false_positive = hit["false_pos"]						
						if false_positive:
							false_pos[contig][allele].append(hit)
							continue				
						hits_to_remove = []		
						allele_name = hit["allele_name"]							
						allele_id = hit["allele_id"] 				
						start = int(hit["qstart"]) 					
						stop = int(hit["qend"])
						score = float(hit["factor"])
						edge_match = hit["edge"] 
						add = False
						found_overlap = False
						if seen_regions:
							for region in sorted(seen_regions):
								add_region = False
								remove_region = True
								replace = False
								current_hit_allele = seen_regions[region]["hit"]["allele_name"]
								r_start = int(region.split("*")[0])					
								r_stop = int(region.split("*")[1])
								ignore_flag = False
								if start <= r_start and stop <= r_stop and stop >= r_start:
									distance = int(abs(stop-r_start))
									r_distance = int(r_stop-r_start)
									overlap_frac = float(distance/r_distance)*100
									overlap_dist = max(overlap_frac,10)
									#print("1",start,stop,allele,contig,region,distance,score,seen_regions[region]["score"],overlap_dist)									
									if distance <= 10:
										add_region = True
										remove_region = False
									#~ if allele_name in longer_overlaps and current_hit_allele in longer_overlaps:
										#~ add_region = True
										#~ remove_region = False		
									if score > seen_regions[region]["score"] or add_region:
										add = True
										if remove_region:
											hits_to_remove.append(region)								
									else:
										false_positive = True
									found_overlap = True									
								if start >= r_start and start <= r_stop and stop >= r_stop:
									distance = int(abs(start-r_stop))
									r_distance = int(r_stop-r_start)
									overlap_dist = float(distance/r_distance)	
									overlap_dist = max(overlap_dist,10)																	
									if distance <= 10:
										add_region = True
										remove_region = False
									#~ if allele_name in longer_overlaps and current_hit_allele in longer_overlaps:
										#~ add_region = True
										#~ remove_region = False												
									if score > seen_regions[region]["score"] or add_region:
										add = True
										if remove_region:
											hits_to_remove.append(region)								
									else:
										false_positive = True
									found_overlap = True	
																									
								if start >= r_start and start <= r_stop and stop >= r_start and stop <= r_stop:
									#~ if allele_name in longer_overlaps and current_hit_allele in longer_overlaps:
										#~ add_region = True
										#~ remove_region = False											
									if score > seen_regions[region]["score"] or add_region:
										add = True
										if remove_region:
											hits_to_remove.append(region)								
									else:
										false_positive = True
									found_overlap = True	
																										
								if start <= r_start and start <= r_stop and stop >= r_start and stop >= r_stop:
									#~ if allele_name in longer_overlaps and current_hit_allele in longer_overlaps:
										#~ add_region = True
										#~ remove_region = False																		
									if score > seen_regions[region]["score"] or add_region:
										add = True
										if remove_region:
											hits_to_remove.append(region)								
									else:
										false_positive = True
									found_overlap = True	
																																
						else:
							new_region = str(start)+"*"+str(stop)								
							seen_regions[new_region] = {"score":score,"hit":hit}
						if not found_overlap:
							new_region = str(start)+"*"+str(stop)								
							seen_regions[new_region] = {"score":score,"hit":hit}
						if false_positive:
							false_pos[contig][allele].append(hit)
							continue
						if hit["region_type"] == "ISE":
							if float(hit["cov"]) < .05:
								if not hit["edge"]:
									false_pos[contig][allele].append(hit)
						if add:
							new_region = str(start)+"*"+str(stop)
							for to_remove in hits_to_remove:
								if to_remove in seen_regions:
									hit_to_remove = seen_regions[to_remove]["hit"]
									allele_to_remove = hit_to_remove["allele"]
									allele_id_to_remove = hit_to_remove["allele_id"]
									false_pos[contig][allele_to_remove].append(hit_to_remove)
									seen_regions.pop(to_remove,None)
							seen_regions[new_region] = {"score":score,"hit":hit}									
						if allele not in alleles_to_grab[species]:
							alleles_to_grab[species][allele] = []
							allele_annotations[allele] = ""			
						if allele_id not in alleles_to_grab[species][allele]:
							alleles_to_grab[species][allele].append(allele_id)
						else:
							continue
						allele_annotations[allele] = "None Found"
			for allele in false_pos[contig]:
				for hit in false_pos[contig][allele]:
					if hit in results_dict[in_file]["contigs"][contig][allele]:
						results_dict[in_file]["contigs"][contig][allele].remove(hit)
	for species in alleles_to_grab:
		for allele in alleles_to_grab[species]:
			for allele_id in alleles_to_grab[species][allele]:
				allele_count+=1									
	print("Grabbing allele information from pubMLST for {} alleles".format(allele_count))				
	q=0
	for species in alleles_to_grab:
		for allele in alleles_to_grab[species]:
			if allele not in allele_info:
				allele_info[allele] = {}		
			for allele_id in alleles_to_grab[species][allele]:
				url = "http://rest.pubmlst.org/db/pubmlst_{}_seqdef/loci/{}/alleles/{}".format(species,allele,allele_id)
				request = http.request('GET',url)
				request_data = json.loads(request.data.decode('utf-8'))					
				#pp.pprint(request_data)
				if allele_id not in allele_info[allele]:
					allele_info[allele][allele_id] = {}
				if "flags" in request_data:
					flags = request_data["flags"]
				else:
					flags = []
				if "comments" in request_data:
					comments = request_data["comments"]
				else:
					comments = []
				for field in request_data:
					if "mutation" in field:
						allele_info[allele][allele_id][field] = request_data[field]
				allele_info[allele][allele_id]["flags"] = flags
				allele_info[allele][allele_id]["comments"] = comments
				q+=1
				if (q%500) == 0:
						print("Completed {} so far".format(str(q)))
		
		
	print("Compiling results")
	for in_file in results_dict:	
		for contig in results_dict[in_file]["contigs"]:
			if contig == "species":
				continue			
			for allele in results_dict[in_file]["contigs"][contig]:
				if results_dict[in_file]["contigs"][contig][allele]:
					for hit in results_dict[in_file]["contigs"][contig][allele]:
						false_positive = hit["false_pos"]
						hit["annotations"] = allele_annotations[allele]
						if not false_positive:
							identity = float(hit["identity"])
							edge_match = hit["edge"]
							qseq = hit["qseq"]
							cov = hit["cov"]
							allele_id = hit["allele_id"]	
							sequence = qseq
							DNA_flag = True						
							if identity < 100.0:
								for letter in sequence:
									if letter not in DNA:
										DNA_flag = False
										break
								if DNA_flag:							
									seq_obj = Seq(sequence, IUPAC.unambiguous_dna)
								hit["flags"] = []					
								hit["new"] = True
								hit["allele_id"] = "new_allele_similar_to_{}_identity({})%_cov({})%".format(allele_id,identity,round(cov*100,2))
							elif edge_match:
								hit["flags"] = []
							else:
								if allele in allele_info:
									hit["flags"] = allele_info[allele][allele_id]["flags"]
									hit["comments"] = allele_info[allele][allele_id]["comments"]
									for field in allele_info[allele][allele_id]:
										if "mutation" in field:
											hit[field] = allele_info[allele][allele_id][field]			
								if sequence != "N/A":
									for letter in sequence:
										if letter not in DNA:
											DNA_flag = False
											break
								if DNA_flag:						
									seq_obj = Seq(sequence, IUPAC.unambiguous_dna)
									
							if not edge_match and DNA_flag and sequence != "N/A":
								protein_sequence = seq_obj.translate(table=11)
								protein_sequence_spliced = protein_sequence[:-1]
								internal_stops = int(protein_sequence.count("*"))
								internal_stops_spliced = int(protein_sequence_spliced.count("*"))
								if internal_stops > 1:
									if "internal stop codon" not in hit["flags"]:
										new_flag = hit["new"]
										hit["flags"].append("internal stop codon")	
										if "N/A" in hit["flags"]:
											hit["flags"].pop("N/A")
								elif internal_stops_spliced > 0:
									if "internal stop codon" not in hit["flags"]:
										new_flag = hit["new"]
										hit["flags"].append("internal stop codon")	
										if "N/A" in hit["flags"]:
											hit["flags"].pop("N/A")							
							 	
	return results_dict
			
			
def blast_command(blast_db,query_file,threads_to_use,nucl):
	if nucl:
		results = check_output(["blastn","-db",blast_db,"-outfmt","6 qseqid sseqid pident qcovus mismatch gapopen qstart qend sstart send score length qlen slen sframe qseq","-query",query_file,"-num_threads",threads_to_use,"-max_target_seqs","200"], shell=False)
	else:
		results = check_output(["blastx","-db",blast_db,"-outfmt","6 qseqid sseqid pident qcovus mismatch gapopen qstart qend sstart send score length qlen slen sframe qseq","-query",query_file,"-num_threads",threads_to_use,"-query_gencode","11","-max_target_seqs","200"], shell=False)
	return results


def run_blast(working_dir,in_file,file_name,threads_to_use,blast_dir,seq_dict):
	processes = int(threads_to_use)
	print("Blasting against Hinfluenzae capsule database with",threads_to_use,"workers")
	## Setup pool for multiprocessing ###
	pool = Pool(processes)	
	folders = {}
	i=0
	q=0
	### Setup 2d array to arrange jobs for multiprocessing ###
	folders[i] = {}
	folders = []
	for folder in os.listdir(blast_dir):
		folders.append(folder)
	### setup query file path and check if blast DB is protein or nucleotide ###
	query_file = os.path.join(working_dir,file_name)
	blast_results_list = []		
	nucl_dict = {}
	for folder in folders:
		nucl = True
		for item in os.listdir(os.path.join(blast_dir,folder)):
			if ".pin" in item:
				nucl = False
		nucl_dict[folder] = nucl
	### list comprehension to launch blast_command for each thread - sends all jobs to pool and pool updates jobs as others finish ###
	blast_time = [pool.apply_async(blast_command,args=(os.path.join(blast_dir,folder,folder),query_file,"1",nucl_dict[folder])) for folder in folders]
	output = [result.get() for result in blast_time]
	blast_final_results = [re.split(b"\n",item.rstrip()) for item in output]
	print("Completed BLAST for",in_file)	
	pool.terminate()
	#input_file = input_file.replace(".fasta","").replace(".fna","")	
	final_dict = analyze_blast(blast_final_results,in_file,seq_dict)
	return final_dict	



def generate_st_predictions(data):
	out_path = os.path.join(OUTPUT_DIR,"serotype")
	hi_set = {}
	st_results = {"Serotype":[]}
	for query in sorted(data):
		species = data[query]["species"]
		if species == "hinfluenzae":
			hi_set[query] = data[query]
			
	with open(os.path.join(out_path,"serotype_predictions_{}.tab".format(time.time())),"w") as f:
		f.write("Query\tST\tGenes_Present\tNotes\n")
		for query in sorted(hi_set):
			st_dict = {"sample_name":query,"predicted_st":"","baseST":"","genes":[]}                    
			found_result = False            
			partial_set = {}
			seen_genes = {}
			for contig in hi_set[query]["contigs"]:
				partial_set[contig] = {}
				for allele in sorted(hi_set[query]["contigs"][contig]):
					partial_set[contig][allele] = []
					if hi_set[query]["contigs"][contig][allele]:
						for hit in hi_set[query]["contigs"][contig][allele]:
							allele_name = hit["allele_name"]
							for st in serotypes:
								if st not in seen_genes:
									seen_genes[st] = {"essential":[],"nonessential":[]}
								for gene in serotypes[st]["essential"]:
									if allele_name == gene:
										found_result = True
										if hit["cov"] < 1 and not hit["edge"]:
											hit["partial"] = True
											partial_set[contig][allele].append(hit)
										else:
											hit["partial"] = False                                        
										seen_genes[st]["essential"].append(hit)
								for n_gene in serotypes[st]["nonessential"]:
									if allele_name == n_gene:    
										seen_genes[st]["nonessential"].append(hit)    
			to_remove = []
			for contig in partial_set:
				for allele in partial_set[contig]:
					for partial_hit in partial_set[contig][allele]:
						partial_hit_start = int(partial_hit["qstart"])
						partial_hit_stop = int(partial_hit["qend"])
						partial_hit_name = partial_hit["allele_name"]
						for allele in hi_set[query]["contigs"][contig]:
							for hit in hi_set[query]["contigs"][contig][allele]:
								hit_type = hit["region_type"]
								if hit_type == "ISE":
									hit_start = int(hit["qstart"])
									hit_name = hit["allele_name"]
									hit_id = str(hit["allele_id"])
									hit_stop = int(hit["qend"])
									start_dist = abs(partial_hit_start - hit_stop)
									stop_dist = abs(partial_hit_stop - hit_start)
									if start_dist < 20 or stop_dist < 20:
										for sg in seen_genes:
											for hit in seen_genes[sg]["essential"]:
												if hit["allele_name"] == partial_hit_name:
													if hit != partial_hit:
														pass
													else:
														hit["disrupted"] = hit_name +"_"+hit_id
			
			current_count = 0
			top_st = None
			matching_num_sg = 0
			unique_st = False
			current_st = ""
			first_st = True
			contamination = False
			#first_seen_unique_partial = False
			seen_sts = []
			if found_result:
				for sg in seen_genes:
					for hit in seen_genes[sg]["essential"]:
						allele_name = hit["allele_name"]
						partiality = hit["partial"]
						if allele_name in ST_unique:
							if first_st:
								current_st = ST_unique[allele_name]
								first_st = False
								if partiality:
									first_seen_unique_partial = True
								seen_sts.append(current_st)
							else:
								if ST_unique[allele_name] != current_st:
									if ST_unique[allele_name] not in seen_sts:
										seen_sts.append(ST_unique[allele_name])
									contamination = True
			if not contamination:
				if found_result:
					for st in seen_genes:
						gene_list = []
						for hit in seen_genes[st]["essential"]:
							allele_name = hit["allele_name"]
							partiality = hit["partial"]
							if allele_name not in gene_list:
								gene_list.append(allele_name)
						if len(gene_list) > current_count:
							top_st = st
							current_count = len(gene_list)
							current_gene_list = gene_list
							notes = []        
							seen_alleles = []
							final_st = top_st
							to_remove = []
							quality_check = {}
							flag_list = []
					for hit in seen_genes[top_st]["essential"]:
						allele_name = hit["allele_name"]
						partiality = hit["partial"]
						flags = hit["flags"]
						if allele_name not in quality_check:
							quality_check[allele_name] = {"full_length_match":False,"list":[],"flags":True}
						if not partiality:
							quality_check[allele_name]["full_length_match"] = True
						if len(flags) == 0:
							quality_check[allele_name]["flags"] = False
						quality_check[allele_name]["list"].append(hit)
					for allele in quality_check:
						if len(quality_check[allele]["list"]) > 1:
							full_length_match = quality_check[allele]["full_length_match"]
							if full_length_match:
								for hit in quality_check[allele]["list"]:
									allele_name = hit["allele_name"]
									if hit["partial"]:
										to_remove.append(hit)
							flag_hit = quality_check[allele]["flags"]
							if not flag_hit:
								for hit in quality_check[allele]["list"]:
									allele_name = hit["allele_name"]
									if len(hit["flags"]) > 0:
										to_remove.append(hit)
					for hit in to_remove:
						allele_name = hit["allele_name"] 
						if hit in seen_genes[top_st]["essential"]:
							seen_genes[top_st]["essential"].remove(hit)
					
					for hit in seen_genes[top_st]["essential"]:
						allele_name = hit["allele_name"]
						if "sequence" in hit:
							hit.pop("sequence")
						st_dict["genes"].append(hit)
						partiality = hit["partial"]
						flag_list = hit["flags"]
						if "disrupted" in hit:
							note_to_add = "{} disrupted by {}".format(allele_name,hit["disrupted"])
							if note_to_add not in notes:
								notes.append(note_to_add)
							final_st = "NT"
						elif partiality:
							if not hit["new"] or hit["cov"] < 0.95:
								final_st = "NT"
								hit_cov = round(hit["cov"]*100,2)
								notes.append("{} fragmented ({}% cov)".format(allele_name,hit_cov))
						
						if len(flag_list) > 0:
							flags = ",".join(flag_list)
							if "phase" in flags:
								final_st = "NT"
								notes.append("phase variable OFF in {}".format(allele_name))                
							elif "internal stop" in flags:
								final_st = "NT"
								notes.append("internal stop in {}".format(allele_name))
								
					for gene in current_gene_list:
						if gene in ST_unique:
							unique_st = True
			else:
				final_st = "Contaminated"
				serotypes_found = ",".join(seen_sts)
				current_gene_list = []
				final_notes = "Found genes for serotypes {}, possible contamination".format(serotypes_found)
			if not contamination:
				if top_st == None:
					notes = []
					current_gene_list = []
					final_st = "NT"
					top_st = None
					final_genes = None
					notes.append("No Capsule Genes Found")

				if unique_st == False and current_gene_list:
					final_st = "NT"
					top_st = "Unclear"
					notes.append("Capsule genes present shared across multiple STs")
				if unique_st != False:
					for st_gene in serotypes[top_st]["essential"]:
						if st_gene not in current_gene_list:
							notes.append("missing {}".format(st_gene))
							final_st = "NT"                
					
				if len(notes) == 0:
					notes.append("All essential capsule genes intact and present")
					
				if top_st != None:
					final_notes = "{} backbone: ".format(top_st)+",".join(sorted(notes))
				elif top_st == None:
					final_notes = "No Backbone: "+",".join(sorted(notes))
				else:
					final_notes = ",".join(sorted(notes))
			final_genes = ",".join(sorted(current_gene_list))
			result_line = "{}\t{}\t{}\t{}\n".format(query,final_st,final_genes,final_notes)
			st_dict["predicted_st"] = final_st
			st_dict["baseST"] = final_st
			st_results["Serotype"].append(st_dict)
			f.write(result_line)

	with open(os.path.join(OUTPUT_DIR,"serotype","serotype_results.json"),"w") as f:
		json.dump(st_results,f)
		
						
	
def main():
	parser = argparse.ArgumentParser(description="Script for predicting serotype of hinfluenzae genomes")
	parser.add_argument('-d','--indir',help="Input Dir", required=True)
	parser.add_argument('-o','--out',help="Output Directory", required=True)
	parser.add_argument('-t','--threads',help="Number of Threads to use (default=1)",default="1")
	args = vars(parser.parse_args())
	### Print Args ###
	print ("hinfluenzae capsule characterization")
	print("v1.0")
	print("#######################################################################################################")
	print ("Running with the following parameters:")
	for arg in args:
		print (arg,":",args[arg])
	##################

	### Get working directory, setup temp folder, set output ###
	working_dir = args["indir"]
	char_data = {}
	temp_dir = tempfile.mkdtemp()
	set_output(args["out"])
	main_blast_dir = PUBMLST_DB
	### Get pubMLST schemes and store in scheme_data dict ###	
	print("#######################################################################################################")	
	print("Step 1. BLASTing against hinfluenzae capsule DB")
	### Setup results variables ### 
	results = []
	results_dict = {}
	#final_results_dict = {}
	seq_dict = OrderedDict()
	#for existing_json in 		
	existing_results = []
	if not os.path.isdir(os.path.join(OUTPUT_DIR,"json")):
		os.system("mkdir {}".format(os.path.join(OUTPUT_DIR,"json")))
	else:
		for existing_json in os.listdir(os.path.join(OUTPUT_DIR,"json")):				
			if "raw" in existing_json:
				existing_json = existing_json.split("_raw")[0]
				existing_results.append(existing_json)
	for in_file in os.listdir(working_dir):
		include = False
		for ext in fasta_extensions:
			if in_file.endswith(ext):
				include = True
				file_ext = ext
				break
		if include:
			file_name = in_file
			in_file = in_file.replace(file_ext,"")
			seq_dict[in_file] = {"species":"","contigs":{},"file_name":file_name}			
			blast_species = "hinfluenzae"			
			blast_dir = os.path.join(ABS_PATH,main_blast_dir,blast_species)		
			with open(os.path.join(working_dir,file_name),"rU") as f:							
				for seq_record in SeqIO.parse(f,"fasta"):
					id = seq_record.id					
					if id not in seq_dict[in_file]["contigs"]:
						seq_dict[in_file]["contigs"][id] = {}
					seq = seq_record.seq
					length = len(seq)
					seq_dict[in_file]["contigs"][id]["seq"] = seq
					seq_dict[in_file]["contigs"][id]["file_name"] = in_file
					seq_dict[in_file]["contigs"][id]["length"] = length
					seq_dict[in_file]["species"] = blast_species
					seq_dict[in_file]["contigs"][id]["alleles"] = {}
			if in_file not in existing_results:
				final_dict = run_blast(working_dir,in_file,file_name,args["threads"],blast_dir,seq_dict)
				#with open(os.path.join(OUTPUT_DIR,"json","{}_raw_results.json".format(in_file)),"w") as f:
				with open(os.path.join(OUTPUT_DIR,"json","{}_raw_results.json".format(in_file)),"w") as f:
					json.dump(final_dict,f)					
			else:
				in_file_path = os.path.join(OUTPUT_DIR,"json",in_file+"_raw_results.json")
				with open(in_file_path) as f:
					final_dict = json.load(f)
				print("Found existing BLAST results for {}, using those results (use -fr flag to force overwriting of results)".format(in_file))
			results_dict.update(final_dict)				
	### Analyze results in results dict ###	
	
	print("#######################################################################################################")	
	print("Step 2. Parsing BLAST results")							
	final_results_dict=analyze_results(results_dict,args["threads"])
	print("#######################################################################################################")
	print("Step 3. Running analyses")		
	for in_file in final_results_dict:
		with open(os.path.join(OUTPUT_DIR,"json","{}_final_results.json".format(in_file)),"w") as f:		
			temp_dict = {}
			temp_dict[in_file] = final_results_dict[in_file]			
			json.dump(temp_dict,f)			
	if not os.path.isdir(os.path.join(OUTPUT_DIR,"gff")):
		os.system("mkdir {}".format(os.path.join(OUTPUT_DIR,"gff")))
	if not os.path.isdir(os.path.join(OUTPUT_DIR,"serotype")):
		os.system("mkdir {}".format(os.path.join(OUTPUT_DIR,"serotype")))
	generate_st_predictions(final_results_dict)
	create_gff(final_results_dict,seq_dict)
				



if __name__ == "__main__":
	main()
