#!/usr/bin/python
# coding=utf-8

# This script takes the allele frequencies (RO and AO) estimated by Freebayes and converts them in the .SYNC format used 
# by Popoolation2 to estimate Fst
# It can handle multiple alleles as produce by Freebayes with the option -n-alleles=2 or 3 or 4. Even higher but it's no use, 
# because THIS SCRIPT HANDLE MAXIMUM 4 ALLELES per position (the reference, plus additional max 3. See below for explanation at line 80-90)


# Should be applied on vfc file filtered for coverage and missing genotypes.
# If missing genotype is found, a warning is issued. Remember then to save the log

# See below for details about Popoolation2 SYNC format. 
# Remove header (printed for reference) manually from output file before running Popoolation2

#Requires: 
#1: filein
# Typically a .vcf file filtered for Quality e depth of coverage

#2: fileout

# Example: 
# python trasform_vcf_2_PREsync.py input.vcf output.temp 1>1_run.log 2>1_run.err

import sys

filein=open(sys.argv[1],'r')
fileout=open(sys.argv[2],'w')

# ----> write header (for referece, but REMOVE before running PoPoolation2)
for line in filein:
	if line.startswith("#CHROM"):
		line=line.rstrip("\n")
		line=line.lstrip("#")	
		l=line.split("\t")
		p=len(l)-9
		fileout.write("\t".join(l[0:2])+"\t"+"REF"+"\t"+"\t".join(l[9:len(l)])+"\n")
		break	
print p


# ----> builds a dictionary to use in making changes to the list 'a'
# A:T:C:G:N:*
a=["0","0","0","0","0","0"]
z={'A':0,'T':1,'C':2,'G':3,'N':4,'*':5}


filein=open(sys.argv[1],'r')
for line in filein:
	if not line.startswith("#"):
		line=line.rstrip("\n")
		l=line.split("\t")
		
		
		a=["0","0","0","0","0","0"]
		REF=l[3].split(",")	# referece allele - can only have one genotype, but I still split it. 
		ALT=l[4].split(",") # alternative allele - can be as many as -n-alleles as sepcified in FreeBayes. THIS SCRIPT HANDLES MAXIMUM 3 ALTERNATIVE + 1 REFERENCE ALLELE (see here below) !!!
	
			
		# I consider only calls in which MAXIMUM 3 alleles have been called as ALTERNATIVE.
		if len(ALT)<4:
			fileout.write("\t".join(l[0:2])+"\t"+"N")
			
			
			# ----> this is to distinguish between invariant (ALT[0]==".") and variant sites.
			if ALT[0]==".":
				alleles=REF
			else:
				alleles=REF+ALT
			
	
			# ----->  when any allele has length > 1 (i.e. > SNP)
			if any(len(allele)>1 for allele in alleles):     # if an genotype of the list is > 1, i.e. a haplotype longer than a SNP
			#	print "COMPLEX HAPLO"
			
				# take the oservaton counts for the alleles
				for i in range(p):
					s=l[9+i]
					if s==".":
						fileout.write(("\t"+"."))
						print "missing genotype at"+"\t"+("\t".join(l[0:2]))
					else:
						m=s.split(":")	
						RO=m[2].split(",")	
						AO=m[4].split(",")
						OBS=RO+AO
				
						# I put in the list a in position n, the value of the observation counts in position n
						# THIS MEANS THAT I RANDOMLY ASSIGN A NUCLEOTIDE TO HAPLOTYPES > 1 SNP !!!!
						for n in range(0,len(alleles)):		
							a[n]=OBS[n]
						fileout.write("\t"+":".join(a))
					
			# ----->   when all genotypes of the list are SNPs
			if all(len(allele)==1 for allele in alleles):
				
				# take the oservaton counts for the alleles
				for i in range(p):
					s=l[9+i]
					if s==".":		
						fileout.write(("\t"+"."))
						print "missing genotype at"+"\t"+("\t".join(l[0:2]))
					else:
						m=s.split(":")	
						RO=m[2].split(",")	
						AO=m[4].split(",")
						OBS=RO+AO
					
						for n in range(0,len(alleles)):
							a[z[alleles[n]]]=OBS[n]
						fileout.write("\t"+":".join(a))
			fileout.write("\n")	
				
			
# ------------------------------------------------------------------------------
# A synchronized file, for example (from help of Popoolation2)

# Unknown_group_104	5943	N	0:0:10:0:0:0	0:0:10:0:0:0	0:0:10:0:0:0
# Unknown_group_104	5944	N	0:8:0:0:0:0	0:8:0:0:0:0	0:8:0:0:0:0
# Unknown_group_104	5945	N	0:0:0:8:0:0	0:0:0:8:0:0	0:0:0:8:0:0
# Unknown_group_104	5946	N	0:0:9:0:0:0	0:0:9:0:0:0	0:0:9:0:0:0
# Unknown_group_104	5947	N	0:7:0:0:0:0	0:7:0:0:0:0	0:7:0:0:0:0
 
# col1: reference contig (chromosome)
# col2: position in the reference contig
# col3: reference character
# col4: population 1
# col5: population 2
# coln: population n
 
# population data are in the form
# A:T:C:G:N:*
# A: count of character A
# T: count of character T
# C: count of character C
# G: count of character G
# N: count of character N
# *: deletion, count of deletion

