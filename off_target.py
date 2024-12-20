#!/usr/bin/env python3

import sys
import argparse
import pysam
import os

outfileslist=[]
		

def function(infile,outdir,bed):
	for file in infile:
		#Check if a BAM file is indexed, and index it if not.
		bai_file = file + '.bai'
		if not os.path.exists(bai_file):
			print(f"Index not found for {file}. Creating index...")
			pysam.index(file)
			print(f"Index created for {file}.")
			print(f"Continuing with analysis...")
		bamfile=pysam.AlignmentFile(file,"rb")
		#Foolprofing in case the script is not run in the BAM files directory:
		if "/" in file:
			outfilename=file.split("/")[-1][:-11]+"_off_target.tsv" #Default input format is samplename.sorted.bam
		else:
			outfilename=file[:-11]+"_off_target.tsv"
		outfileslist.append(outfilename)
		reads=0
		reads_ON=0
		with open(outdir+outfilename,'w') as g:
			g.write("Assay"+'\t'+"Reads"+'\t'+"Reads_percent"+'\n')	
			reads_total=0
			for bamline in bamfile.fetch(until_eof=True, multiple_iterators=True): # Total reads, not including seconday or supplementary alignments (multiple alignments), until_eof=True flag is needed to include unaligned reads.
				if bamline.is_secondary==True or bamline.is_supplementary==True:
					continue
				reads_total+=1
			
			
			with open(bed) as b:
				for bedline in b: # Count number of reads per annotated amplicon in the bed file
					reads=0
					bed_parts=bedline.split("\t")
					for bamline in bamfile.fetch(bed_parts[0],int(bed_parts[1]),int(bed_parts[2]),multiple_iterators=True): # Use bamline.fetch to extract each bamline matching the genomic coordinates of the amplicon
						if bamline.is_secondary==True or bamline.is_supplementary==True:
							continue
						bamline=str(bamline)
						bam_parts=bamline.split("\t")
						reads+=1
						reads_ON+=1
						
					if reads>150: # Write down the total number of aligned reads and percentage of total reads per amplicon:
						g.write(bed_parts[3].rstrip()+'\t'+str(reads)+'\t'+str(100*(reads/reads_total))+'\n')
						
			
			reads=0
			for bamline in bamfile.fetch(until_eof=True):
				if bamline.is_secondary==True or bamline.is_supplementary==True:
					continue
				bamline=str(bamline)
				bam_parts=bamline.split("\t")
				reads+=1

			g.write("Off-target"+'\t'+str(reads-reads_ON)+'\t'+str(100*(reads-reads_ON)/reads)+'\n')
		
		bamfile.close()
	return outfileslist

def main(infile,outdir,bedfile,merge=False):
	#Foolproofing so output directory is passed no matter of whether the diagonal bar is provided in the outdir argument or not: 
	if outdir is None:
		outdir=""
	else:
		if outdir[-1]!="/":
			outdir+="/"
	
	print("Starting analysis...")		
	outfiles=function(infile,outdir,bedfile)
	print("Done!")	
	
				
if __name__=='__main__':
	def parseArgs():
		parser=argparse.ArgumentParser(description="Calculate the percentage of off-target reads")
		parser.add_argument('-i','--infile',dest='infile',help='Path to input file(s)',nargs='+',required=True)
		parser.add_argument('-b','--bed',dest='bedfile',help='Path to bed file annotation',required=True)
		parser.add_argument('-o','--outdir',dest='outdir',help='Output directory', required=False)
		args=parser.parse_args(sys.argv[1:])
		return(args)
	args=parseArgs()	
	main(args.infile, args.outdir, args.bedfile)
