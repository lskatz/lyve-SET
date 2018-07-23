#!/usr/local/bin/python
# Khushbu Patel | 05/15/2018
# Corrects the filter and Calculates number of bases falsely assigned incorrect filter when the bases pass the consensus and coverage. Prints out a new vcf with filters corrected. 
# Python 3.4
# Usage: ./correct_vcf_filter.py inputfile.vcf
# Output will be another vcf file, named correct_vcf.vcf; This script will not overwrite the original VCF

import sys
import os
import subprocess

infile= sys.argv[1]	# Takes the file name as a command line argument
f=open(infile,"r")

basename = os.path.basename(infile)
outfile_vcf = basename + "_corrected.vcf"
outfile_txt = basename + "_corrected.txt"
temp = []
count = 0
out = []
ID = ''
old_filter = 0


with open(outfile_vcf, 'a') as f1:
	with open(outfile_txt, 'a') as f2:
		f2.write( "CHROM\tPOS\tNEW_FILTER\tDEPTH\tFREQ\tOLD_FILTER\n")
		for line in f:
			line =line.rstrip()
			if line.startswith('#'):
				ID = line					# Printing the headers
				f1.write(ID)
				f1.write("\n")
		
			else:
				array = line.split()
				temp = array[9].split(':')
			
				if(array[6] != "PASS"):
					if(len(temp)> 6):
						old_filter = array[6]		# Stores old filter
						temp[6] = temp[6].replace('%','')
						if(int(temp[3]) >= 20 and float(temp[6]) < 5.0):		# if coverage and consensus both meet, change the filter to PASS
							#print(line, "--Incorrect Filter!")			# Sanity check!
							array[6] = "PASS"
							count += 1					# Count the number of sites that have been assigned wrong filter
							out = '\t'.join(array)
							f1.write(out)
							str = array[0]+'\t'+array[1]+'\t'+array[6]+'\t'+temp[3]+'\t'+temp[6]+'\t'+old_filter+'\n'
							f2.write(str)
						
						
						else:							# if allele frequency is greater than 5% or coverage does not meet; keep the filter
							out = '\t'.join(array)
							f1.write(out)
							f1.write("\n")
	
					else:									# else - print everything as is
						out = '\t'.join(array)
						f1.write(out)
						f1.write("\n")
				
				else:									# else - print everything as is
					out = '\t'.join(array)
					f1.write(out)
					f1.write("\n")
	
	
	
print('Number of sites that had been assigned wrong filters %d'%count)
