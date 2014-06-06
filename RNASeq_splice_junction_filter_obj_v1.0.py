# -*- coding: cp936 -*-
"""
The main purpose of this script is to filter the RNAseq variant artifact calls
due to misalignment caused by splice-junction or small deletion.

=============================
Usage: python RNASeq_splice_junction_filter_obj_v1.0.py

-h help

===vcf parameter setup=====
-i (vcf)files to processed                      *[No default value]

-a path for the input vcf file			[default value: current path]

-s suffix for vcf files to be processed         [default value "vcf"]

===bam parameter setup=====
-I bam input file paired with vcf file		*[No default value]

-A path for the paired bam input file		[default value: current path]

-S suffix for the bam files to be processed     [default value "bam"]

===filter parameter setup=====

-L bases to remove on both side of read         [default value 6]

-F exon splicing border length filter		[default value 10]

-m minimal coverage requirement			[default value 8]

-E mode of processing 				[default value 1]
	1 means one single sample-vcf vs one single sample bam
	2 means one compiled multi sample vcf vs multiple single sample bam

-X count all					[default value 0]
	0 means only count when they are unlikely to be false positive
	1 means count all that show up

-c count all files not just the ones with calls	[default value 0]
	0 means does not count the ones without the calls
	1 means does count the ones without the calls

Less important Parameters:
==========================

-O output path					[default value: current folder]

-l unique_id length for the infile		[default value 2]

-p prefix of output files                       [default value:"Remove_Splice_artifact(RSA)"]

-d sep_char within among columns                [default value '\t']

===================
input description:
input files:
1. filtered RNA editing VCF file or mutation call VCF files
2. bam files

======================
output files:
subset of the input file
the format of input file will be the same
the filtered one, filtered column will be marked as "Border_bases_Artifact"

============================
Module requirement:
No additional Python Module is required.

============================
Library file requirement:
Standalone version, no library file is required.
============================



"""

##Copyright
##By Liye Zhang
##Boston University
##Version: 1.0
##Compatible Python Version:2.4 or above
##Developed at: 2013-4-23
##Last Update: 2014-1-7

###Code Framework

if __name__ == "__main__":
	###Module Import	
	import sys, csv, getopt, re
	import os
	import math
	from itertools import ifilter
	import subprocess
	
	digit_list=["0","1","2","3","4","5","6","7","8","9"]
	
	##class definition
	class GeneralFile_class:
	
		def __init__(self,name):
			self.filename=name
			self.name_only=name
			self.SEP_CHAR='\t'
			self.SKIP_HEADER=0
			self.SAMPLE_ID_LEN=1 #this will determin how many section will be considered to be the unique ID
			self.SAMPLE_ID_POS=0
			self.UNIQUE_ID_COLUMN=0
			self.FILENAME_SPLIT_CHAR='_'
			self.RECORD=""
			self.AUTOSKIP_HEADER=True
			self.OUTPUT_PATH=os.getcwd()

			
		def generate_sample_id(self,POS=0):
			if '/' not in self.filename:
				filename_list=(self.name_only).split(self.FILENAME_SPLIT_CHAR)
			else:
				infile_path_list=self.filename.split('/')
				infile_name=infile_path_list[-1]
				filename_list=infile_name.split(self.FILENAME_SPLIT_CHAR)
			if self.SAMPLE_ID_LEN==1:
				sample_id=filename_list[POS]
			else:
				filename_list_len=len(filename_list)
				if self.SAMPLE_ID_LEN>filename_list_len:
					self.SAMPLE_ID_LEN=filename_list_len
					
				for i in range(self.SAMPLE_ID_LEN):
					if i == 0 :
						sample_id=filename_list[POS]
					else:
						sample_id+="_"+filename_list[POS+i]
			return sample_id
		
		def outputfilename_gen(self,name_fragment="std_out",suffix="txt",POS=0):
			##version2.0
			if '/' not in self.filename:
				
				#infile_name_list=(self.filename).split(self.FILENAME_SPLIT_CHAR)
				#sample_id=infile_name_list[POS]
				sample_id=self.generate_sample_id()
				output_filename=sample_id+"_"+name_fragment+'.'+suffix
				return output_filename
			else:
				infile_path_list=self.filename.split('/')
				infile_name=infile_path_list[-1]
				print "infile_name",infile_name
				#infile_name_list=(infile_name).split(self.FILENAME_SPLIT_CHAR)
				#sample_id=infile_name_list[POS]
				self.name_only=infile_name
				sample_id=self.generate_sample_id()
				output_filename=sample_id+"_"+name_fragment+'.'+suffix
				return output_filename
		
		def reader_gen(self,FILEPATH=os.getcwd()):
			if '/' in self.filename:
				complete_path=self.filename
			else:
				complete_path=FILEPATH + '/' + self.filename
			reader=csv.reader(open(complete_path,'rU'),delimiter=self.SEP_CHAR,quoting=csv.QUOTE_NONE)
			if self.AUTOSKIP_HEADER==True:
				## this will over-write provided default
				skip_number=0
				row=reader.next()
				while row[0][0]=="#":
					skip_number+=1
					row=reader.next()
				self.SKIP_HEADER=skip_number
			
			reader=csv.reader(open(complete_path,'rU'),delimiter=self.SEP_CHAR,quoting=csv.QUOTE_NONE)
			for i in range(self.SKIP_HEADER):
				reader.next()
			return reader
			
		def output_handle_gen(self,header_file=None,FILEPATH=os.getcwd(),sup_list=[],HEAD_LINE=1):
			##Version2.0
			##Updated 2012-10-31
			'''
			header_file is the file contains the header information
			FILEPATH is the path for the output file
			sup_list is the additional annotations added to the output file header
			HEAD_LINE is the number of header lines extracted from header file and writen into output file
			'''
			if '/' in self.filename:
				complete_path=self.filename
			else:
				complete_path=FILEPATH + '/' + self.filename
			
			self.handle=open(complete_path,'w')
			
			if self.RECORD=="":
				pass
			else:
				self.handle.write(self.RECORD)
			
			if header_file==None:
				pass
			else:
				output_header_file(header_file,self.handle,sup_list,eliminate=0)
	
	class VCF_File_class(GeneralFile_class):
		def __init__(self,name):
			GeneralFile_class.__init__(self,name)
			self.QUAL_COLUMN=5
			self.ALT_COLUMN=4
			self.REF_COLUMN=3
			self.ID_COLUMN=2
			self.COOR_COLUMN=1
			self.CHRO_COLUMN=0
			self.FILTER_COLUMN=6
			self.INFO_COLUMN=7
			self.FORMAT_COLUMN=8
			self.SKIP_HEADER=0
			self.SEP_INFO_START_COLUMN=9
			self.ALT_SEP_CHAR=','
		
		def sample_list_gen(self,FILEPATH=os.getcwd()):
			
			if '/' in self.filename:
				complete_path=self.filename
			else:
				complete_path=FILEPATH + '/' + self.filename
			reader=csv.reader(open(complete_path,'r'),delimiter=self.SEP_CHAR)
			for rows in reader:
				first_item=rows[0]
				if first_item[0]=='#' and first_item[1]!='#':
					sample_list=rows[self.FORMAT_COLUMN+1:]
					break
				else:
					pass
			return sample_list
					
		def reader_gen(self,FILEPATH=os.getcwd()):
			if '/' in self.filename:
				complete_path=self.filename
			else:
				complete_path=FILEPATH + '/' + self.filename
			
			reader=csv.reader(open(complete_path,'rU'),delimiter=self.SEP_CHAR)
			skip_count=0
			rows=reader.next()
			while rows[0][0]=="#":
				skip_count+=1
				rows=reader.next()
			
			reader=csv.reader(open(complete_path,'rU'),delimiter=self.SEP_CHAR)
			for i in range(skip_count):
				reader.next()
			return reader		 
		
		def output_sample_info(self,FILEPATH=os.getcwd()):
			if '/' in self.filename:
				complete_path=self.filename
			else:
				complete_path=FILEPATH + '/' + self.filename
			
			reader=csv.reader(open(complete_path,'rU'),delimiter=self.SEP_CHAR)
			
			
			for rows in reader:
				if rows[0][0:2]=="##":
					pass
				else:
					sample_info=rows[9:]
					break
			return sample_info
		
		def output_handle_gen(self,header_file=None,FILEPATH=os.getcwd(),sup_list=[],HEAD_LINE=1):
			##Version2.0
			##Updated 2012-10-31
			'''
			header_file is the file contains the header information
			FILEPATH is the path for the output file
			sup_list is the additional annotations added to the output file header
			HEAD_LINE is the number of header lines extracted from header file and writen into output file
			'''
			if '/' in self.filename:
				complete_path=self.filename
			else:
				complete_path=FILEPATH + '/' + self.filename
			
			self.handle=open(complete_path,'w')
			
			
			
			if header_file==None:
				self.handle.write("##fileformat=VCFv4\n")
				if self.RECORD=="":
					pass
				else:
					self.handle.write(self.RECORD)
				pass
			else:
				#output_header_file(header_file,self.handle,sup_list,eliminate=0)
				output_header_VCF_file(header_file,self.handle,self.RECORD,sup_list,eliminate=0)
				
		def sample_count(self,FILEPATH=os.getcwd()):
			if '/' in self.filename:
				complete_path=self.filename
			else:
				complete_path=FILEPATH + '/' + self.filename
			reader=csv.reader(open(complete_path,'rU'),delimiter=self.SEP_CHAR)
			sample_count=0
			for rows in reader:
				if rows[0][0]=="#":
					pass
				else:
					sample_info=rows[9:]
					for sample_data in sample_info:
						if sample_data.count(":")==2:
							sample_count+=1
						elif sample_data=="./.":
							sample_count+=1
						else:
							pass
					break
			return sample_count
	
	class SAM_File_class(GeneralFile_class):
		def __init__(self,name):
			GeneralFile_class.__init__(self,name)
			self.QNAME_COLUMN=0
			self.FLG_COLUMN=1
			self.CHRO_COLUMN=2
			self.COOR_COLUMN=3
			self.MAPQ_COLUMN=4
			self.CIGAR_COLUMN=5
			self.RNEXT_COLUMN=6
			self.PNEXT_COLUMN=7
			self.TLEN_COLUMN=8
			self.SEQ_COLUMN=9
			self.QUAL_COLUMN=10
			self.READGROUP_COLUMN=12
			self.SKIP_HEADER=0
		
		def reader_gen(self,FILEPATH=os.getcwd()):
			if '/' in self.filename:
				complete_path=self.filename
			else:
				complete_path=FILEPATH + '/' + self.filename
				
			reader=csv.reader(open(complete_path,'rU'),delimiter=self.SEP_CHAR,quoting=csv.QUOTE_NONE)
			if self.AUTOSKIP_HEADER==True:
				## this will over-write provided default
				skip_number=0
				row=reader.next()
				while row[0][0]=="@":
					skip_number+=1
					row=reader.next()
				self.SKIP_HEADER=skip_number
			
			reader=csv.reader(open(complete_path,'rU'),delimiter=self.SEP_CHAR,quoting=csv.QUOTE_NONE)
			for i in range(self.SKIP_HEADER):
				reader.next()
			return reader	

				
	##common function 
	def output_header_VCF_file(infile,output_handle,cmd_record,sup_list=[],eliminate=0):
		##Version1.0
		##write header into files
		reader=csv.reader(open(infile,"rU"),delimiter="\t")
		first_line=True
		for row in reader:
			
			#if first_line!=True and (row[0][0:17]=="##fileformat=VCFv"):
			#	continue
			
			if row[0][0:2]=="##":
				if first_line and row[0][0:18]!="##fileformat=VCFv4":
					output_handle.write("##fileformat=VCFv4.0\n")
					output_row(output_handle,row,eliminate)
					if cmd_record!="":
						output_handle.write(cmd_record)	
					
				elif first_line and row[0][0:18]=="##fileformat=VCFv4":
					output_row(output_handle,row,eliminate)
					if cmd_record!="":
						output_handle.write(cmd_record)
				
				else:
					#output_handle.write("##fileformat=VCFv4.0\n")
					output_row(output_handle,row,eliminate)	
						
					
			elif row[0][0]=="#" and row[0][1]!="#":
				if first_line==True:
					output_handle.write("##fileformat=VCFv4.0\n")
					if cmd_record!="":
						output_handle.write(cmd_record)
				combined_row=row+sup_list
				output_row(output_handle,combined_row,eliminate)
			else:
				if first_line==True:
					output_handle.write("##fileformat=VCFv4.0\n")
					if cmd_record!="":
						output_handle.write(cmd_record)
					print "quit early"
				break
			first_line=False
	
	def generate_paired_files_by_ID(infiles_combine,unique_ID_length,vcf_suffix="vcf",bam_suffix="bam",sep_char="_"):
		## this script will pair input files based on equal unique_ID
		## Version update, provide flexible pairing for the 2 file case
		
		paired_list=list()
		
		
		if len(infiles_combine)==2:
			##quick function
			vcf_infile=""
			bam_infile=""
			for infile in infiles_combine:
				
				if infile.endswith(vcf_suffix):
					vcf_infile=infile
				elif infile.endswith(bam_suffix):
					bam_infile=infile
				else:
					print "[ERROR] unsupported file suffix"
					print "Please check your input file"
					sys.exit(0)
					
			if vcf_infile!="" and bam_infile!="":
				result_list=[vcf_infile,bam_infile]
				paired_list.append(result_list)
				return paired_list
			else:
				
				print "[ERROR] not enough input, input files are missing"
				print "Please check your input file"
				sys.exit(0)	
			
		else:
			data_dict=dict()
			for infile in infiles_combine:
				if infile.count("/")>0:
					infile_temp=infile.split("/")
					infile_info=infile_temp[-1]
				else:
					infile_info=infile
				infile_list=infile_info.split(sep_char)
				for i in range(unique_ID_length):
					if i==0:
						unique_ID=infile_list[i]
					else:
						unique_ID=unique_ID+"_"+infile_list[i]
				data_dict[unique_ID]=[]
				
			for infile in infiles_combine:
				if infile.count("/")>0:
					infile_temp=infile.split("/")
					infile_info=infile_temp[-1]
				else:
					infile_info=infile
				infile_list=infile_info.split(sep_char)
				for i in range(unique_ID_length):
					if i==0:
						unique_ID=infile_list[i]
					else:
						unique_ID=unique_ID+"_"+infile_list[i]
				data_dict[unique_ID].append(infile)
			
			data_list=data_dict.keys()
			for data_pair in data_list:
				if len(data_dict[data_pair])!=2:
					print "incorrect data_pair id is", data_pair 
					print "Incorrect pairing, Please check your input"
					sys.exit(0)
				else:
					data1,data2=data_dict[data_pair]
					if data1.count(bam_suffix)>=1 and data2.endswith(vcf_suffix):
						ordered_list=[data2,data1]
					elif data1.endswith(vcf_suffix) and data2.count(bam_suffix)>=1:
						ordered_list=[data1,data2]
					else:
						print "incorrect pairing, Please check your input"
						sys.exit(0)
					paired_list.append(ordered_list)
		return paired_list
	
	def list2string(data_list):
		
		data_length=len(data_list)
		for index in range(data_length):
			data=str(data_list[index])
			if index==0:
				output=data
			else:
				output=output+"_"+data
		return output		
			
	def output_row(handle,row,eliminate=0):
		##write row into files
		len_row=len(row)-eliminate
		for i in xrange(len_row):
			if i==(len_row-1):
				handle.write(str(row[i])+'\n')
			else:
				handle.write(str(row[i])+'\t')
				
	def extract_files_with_suffix_from_path_to_full_path(CurrentFolder=os.getcwd(), suffix=""):
		"""
		getting all the files with certain suffix
		"""
		
		#debug
		#print "function exerted"
		tmp_result, result = [], []
		def checkSuffix(f):
			return f.endswith(suffix) and os.path.isfile(CurrentFolder+'/'+f)
	
		#CurrentFolder=os.getcwd()
		folder_level=0
		for dirpath, dirname, filename in os.walk(CurrentFolder):
			folder_level+=1
			#the os.walk will go into the all the folders and even daughter folders
			#and therefore make sure you are using it correctly
			#otherwise the correct value will be replaced by others
			#even if there are mutiple ones, the first one is the route number
			if folder_level==1:
				tmp_result = filename
			else:
				pass
				
		
		if suffix != "":
			for filename in tmp_result:
				if checkSuffix(filename)==True:
					result.append(CurrentFolder+'/'+filename)
			
			return result
		else:
			return tmp_result	

	def full_filter_retrieve_bases(read,cigar,read_start_coor,search_coor,read_quality,READEND_CUTOFF=6,SPLICEJUNCTION_CUTOFF=10):
		
		##update on the deletion criteria: 2013.11.11
		
		##[SectionI]
		##Variable Initialization	
		read_length=len(read)
		filter_out=False    ##filter_out indicates whether this read confidently support the SNV calls
		temp_base=""	    ##temp variable to store the cigar patterns	
		mapping_pattern_list=[]	  ##This list stores the separated cigar patterns for READ alignments
		border_list=[]	    ## This list stores the splice-junction coordinates	
		read_left_set=False	## this record whether left end read has been determined
		read_right_set=False    ## this record whether right end read has been determined
		search_coor_found=False	## this record whether SNV coordinate has been located
			
		##[SectionII]
		##Process the cigar string and separate it into Lists
		for character in cigar:
			if character in digit_list:   ## test whether the string is still a digit 0-9
				temp_base=temp_base+character     ##compile number 
			else:
				fragment_data=temp_base+character   ##other wise create cigar item
				temp_base=""
				mapping_pattern_list.append(fragment_data)
					
		##[Section III]
		##locate the coordinate bases
		current_read_index=-1			## record the read location
		current_coor_index=read_start_coor-1	## record the genome coordinate location
		SNV_base="N"
		
		for data in mapping_pattern_list:
			data_length=int(data[:-1])	## get the cigar length 
			data_type=data[-1]		## get the cigar type
			
			##update the coordinate and read location
			if data_type=="M":		## M means mapped
				current_coor_index+=data_length		##update genome coor
				current_read_index+=data_length		##update read coor
			elif data_type=="N":		## N means splicing (jumpped length)
				border_list.append(current_coor_index)
				current_coor_index+=data_length
				border_list.append(current_coor_index+1)
			elif data_type=="I":		## I means insertion
				current_read_index+=data_length
			elif data_type=="D":		## D means deletion
				current_coor_index+=data_length
			elif data_type=="S":
				current_read_index+=data_length
			else:
				print "incorrect data format,check input"
				sys.exit(0)
			
			## Setup left read border	
			if current_read_index>=(READEND_CUTOFF-1) and read_left_set==False:
				read_left_coor=current_coor_index-(current_read_index-READEND_CUTOFF+1)
				read_left_set=True
				
			## Setup right read border
			if current_read_index>=(read_length-READEND_CUTOFF) and read_right_set==False:
				read_right_coor=current_coor_index-(current_read_index-(read_length-READEND_CUTOFF))
				read_right_set=True
			
			if current_coor_index>=search_coor and search_coor_found==False:
				if data_type=="N":
					SNV_base="N"
					filtered_out=True
					False_Positive_Type="splicing_skip"
					
				elif data_type=="D":
					SNV_base="D"
					#filtered_out=True
					False_Positive_Type="deletion"
					
				else:
					SNV_read_index=current_read_index-(current_coor_index-search_coor)
					SNV_base=read[SNV_read_index]
					SNV_quality=read_quality[SNV_read_index]
					if (ord(SNV_quality)-33) < MINIMUM_QUAL_SCORE:
						filtered_out=True
						False_Positive_Type="low_quality_call"
				search_coor_found=True
		
		## check near splice junction	
		for border_coor in border_list:
			border_distance=abs(border_coor-search_coor)
			if border_distance<=(SPLICEJUNCTION_CUTOFF-1):
				False_Positive_Type="near_splice_junction"
				filter_out=True
				break
		
		## check near read border
		if search_coor<=read_left_coor or search_coor>=read_right_coor:
			False_Positive_Type="near_read_border"
			filter_out=True
			
		if filter_out==False and SNV_base!="D":
			False_Positive_Type="NA"
			
		return filter_out,SNV_base,False_Positive_Type
	
	
	### Specific Functions definiation
	def Filter_Alignment_Artifact_compiled(infiles,bam_infiles):
		
		for infile in infiles:
			
			print "processing compiled infile:", infile
			infile_obj=VCF_File_class(infile)   
			infile_obj.SAMPLE_ID_LEN=unique_id_length   ##unique ID length define
			infile_reader=infile_obj.reader_gen()      ## generate reader
			sample_list=infile_obj.sample_list_gen()  ##generate the sample list
			sample_COUNT=len(sample_list)              
			sample_dict=dict()
			
			#Associate vcf column with bam file
			for sample_ID in sample_list:
				sample_ID_info=sample_ID.split(".")   ##default generated by the vcf combined  
				sample_ID_simple=sample_ID_info[0]    
				for bam_file in bam_infiles:         
					if bam_file.count(sample_ID_simple)==1:
						sample_dict[sample_ID]=bam_file
						break
			
			##Qc check the Association between vcf columns and bam files
			pair_count=len(sample_dict.keys())
			if pair_count!=sample_COUNT:
				print "Error: Incorrect Sample Pairing"
				print "paired sample found is", pair_count
				print "total sample number is", sample_COUNT
				sys.exit(0)
			
			## output initilization
			final_prefix=prefix+"_R"+str(BORDER_FILTER_LENGTH)+"_S"+str()
			detailed_prefix=final_prefix+"_details"
			outfile_name=infile_obj.outputfilename_gen(detailed_prefix,"VCF")
			outfile_clean_name=infile_obj.outputfilename_gen(final_prefix,"vcf")
			outfile=OUTPUT_PATH+'/'+outfile_name
			outfile_clean=OUTPUT_PATH+'/'+outfile_clean_name
			
			##default is with full path
			##otherwise will use outfile_wo_path or outfile_name_only
			
			outfile_obj=VCF_File_class(outfile)
			outfile_clean_obj=VCF_File_class(outfile_clean)
			
			outfile_obj.output_handle_gen(infile,os.getcwd(),sample_list) ##include additional information
			outfile_clean_obj.output_handle_gen(infile)  ## keep original data only
			
			for rows in infile_reader:
				chro=rows[infile_obj.CHRO_COLUMN]
				coor=rows[infile_obj.COOR_COLUMN]
				ref=rows[infile_obj.REF_COLUMN]
				alt=rows[infile_obj.ALT_COLUMN]
				ref=ref.upper()
				alt=alt.upper()
				alt_filter=rows[infile_obj.FILTER_COLUMN]
				data_all=rows[infile_obj.SEP_INFO_START_COLUMN:]
				process_dict=dict()
				result_dict=dict()
				test_dict=dict()
				
				## function to figure out which one to process
				for index in range(sample_COUNT):
					sample_ID=sample_list[index]
					data=data_all[index]
					if data=="./." and count_all_samples == 0:
						process_dict[sample_ID]=0
						result_dict[sample_ID]=["NA"]
						test_dict[sample_ID]=["NA"]
					else:
						process_dict[sample_ID]=1
						result_dict[sample_ID]=[]
						test_dict[sample_ID]=[]
						
				##Processing the splice-junction removal
				
				for index in range(sample_COUNT):
					sample_ID=sample_list[index]
					if process_dict[sample_ID]==1:
						bam_infile=sample_dict[sample_ID]
						
						if alt.count(infile_obj.ALT_SEP_CHAR)>0:   ##check multiple variant per site
							alt_list=alt.split(infile_obj.ALT_SEP_CHAR) 
						else:
							alt_list=[alt]
							
						alt_final_list=[]
						alt_final_dict=dict()
						
						for base in "ADTCG":  ##handle deletion cases
							alt_final_dict[base]=0
							
						for alt_SNV in alt_list:
							if len(ref)==len(alt_SNV) and len(ref)==1:   ##only handle the point mutation(not indel)
								alt_final_list.append(alt_SNV)
								
						##generate the command to generate intermediate sam file
						coor_region=chro + ":" + coor + "-" + coor
						#outfile_sam=bam_infile + ".sam"
						outfile_sam=OUTPUT_PATH + "/" +outfile_name+".sam"
						sam_extract_cmd="samtools view "+bam_infile + " " + coor_region +"> "+  outfile_sam
						subprocess.call(sam_extract_cmd,shell=True)
						
						sam_obj=SAM_File_class(outfile_sam)
						total_correct_coverage=0
						try:	## mainly due to the fact that some samples do not have coverage at all
							samfile_reader=sam_obj.reader_gen()
							
							for alignment in samfile_reader:
								alignment_flag=int(alignment[sam_obj.FLG_COLUMN])
								if alignment_flag > 1024:   
									pass  ##do not consider the duplicates read
								else:
									alignment_coor=int(alignment[sam_obj.COOR_COLUMN])
									alignment_cigar=alignment[sam_obj.CIGAR_COLUMN]
									alignment_fasta=alignment[sam_obj.SEQ_COLUMN]
									alignment_quality=alignment[sam_obj.QUAL_COLUMN]
									filter_result,nt,false_positive_type=full_filter_retrieve_bases(alignment_fasta,alignment_cigar,alignment_coor,int(coor),alignment_quality,BORDER_FILTER_LENGTH)
									#print filter_result,nt
									if nt==ref:
										alt_final_dict[nt]+=1
										total_correct_coverage+=1
									elif nt=="D":
										#total_correct_coverage+=1
										alt_final_dict[nt]+=1
									elif nt=="N":
										pass
									else:
										total_correct_coverage+=1
										if filter_result==False and count_all==0:
											alt_final_dict[nt]+=1
										if count_all!=0:
											alt_final_dict[nt]+=1
						except:
							pass
										
						PASS_SNV=""
						freq_SNV=""
						test_outcome=""
						#print "Total coverage is", total_correct_coverage
						if total_correct_coverage<MINIMUM_COVERAGE:
							test_outcome="Not-enough-coverage"
						else:
							for alt in alt_final_list:
								alt_count=alt_final_dict[alt]
								freq=float(alt_count)*100/total_correct_coverage
								if freq > MINIMUM_ALLEL_FRQ:
									if PASS_SNV=="":
										PASS_SNV=alt
										freq_SNV=str(freq)+'%'
									else:
										PASS_SNV=PASS_SNV+","+alt
										freq_SNV=freq_SNV+","+str(freq)+'%'
									test_outcome="true-positive"
									
							##test deletion in the same coordinate
							deletion_count=alt_final_dict["D"]
							freq=int(float(deletion_count)*100/total_correct_coverage)
							if freq>DELETION_THRESHOLD:
								test_outcome="Candidate-Deletion-artifacts"
								PASS_SNV=""
						
						if test_outcome=="":
							test_outcome="Candidate-splice-artifacts"
						
						rows[infile_obj.FILTER_COLUMN]=test_outcome
						
						output_info=""
						for base in "ADTCG":
							base_count=str(alt_final_dict[base])
							base_info=base+"_"+base_count
							if base=="A":
								output_info=base_info
							else:
								output_info=output_info+":"+base_info
						
						test_dict[sample_ID].append(test_outcome)
						result_dict[sample_ID].append(output_info)
						
				
				
				
				test_list=[]
				count_list=[]
				for index in range(sample_COUNT):
					sample_ID=sample_list[index]
					test_data=test_dict[sample_ID][0]
					test_list.append(test_data)
					count_data=result_dict[sample_ID][0]
					count_list.append(count_data)
					
				filter_output=rows[infile_obj.FILTER_COLUMN]	
				test_output=list2string(test_list)
				
				##Process the data
				if test_list.count("true-positive")>=1:
					## output in the clean handle
					output_row(outfile_clean_obj.handle,rows)
					filter_output=filter_output+";"+test_output
				else:	
					filter_output="false-positive"+";"+test_output
				rows[infile_obj.FILTER_COLUMN]=filter_output
				output_row(outfile_obj.handle,rows+count_list)
				
				
			outfile_obj.handle.close()
			outfile_clean_obj.handle.close()
			
			rm_sam_cmd="rm "+ outfile_sam
			subprocess.call(rm_sam_cmd,shell=True)
			
	def Filter_Alignment_Artifact(infiles):
		
		##Section I: Generate the gene annotation dictionary
			
		for infile,bam_infile in infiles:
			print "Processing infile:", infile
			##Set up vcf infile object
			
			infile_obj=VCF_File_class(infile)
			infile_obj.SAMPLE_ID_LEN=unique_id_length
			infile_reader=infile_obj.reader_gen()
			
			final_prefix=prefix+"_"+str(BORDER_FILTER_LENGTH)
			
			detailed_prefix=final_prefix+"_details"
			outfile_name=infile_obj.outputfilename_gen(detailed_prefix,"VCF")
			outfile_clean_name=infile_obj.outputfilename_gen(final_prefix,"vcf")
			outfile=OUTPUT_PATH+'/'+outfile_name
			outfile_clean=OUTPUT_PATH+'/'+outfile_clean_name
			
			##default is with full path
			##otherwise will use outfile_wo_path or outfile_name_only
			
			outfile_obj=VCF_File_class(outfile)
			outfile_clean_obj=VCF_File_class(outfile_clean)
			
			outfile_obj.output_handle_gen(infile)
			outfile_clean_obj.output_handle_gen(infile)
			
			for rows in infile_reader:
				chro=rows[infile_obj.CHRO_COLUMN]
				coor=rows[infile_obj.COOR_COLUMN]
				ref=rows[infile_obj.REF_COLUMN]
				alt=rows[infile_obj.ALT_COLUMN]
				ref=ref.upper()
				alt=alt.upper()
				alt_filter=rows[infile_obj.FILTER_COLUMN]
				
				##Processing the splice-junction removal
				
				if alt.count(infile_obj.ALT_SEP_CHAR)>0:   ##check multiple variant per site
					alt_list=alt.split(infile_obj.ALT_SEP_CHAR) 
				else:
					alt_list=[alt]
					
				alt_final_list=[]
				alt_final_dict=dict()
				
				for base in "ADTCG":  ##handle deletion cases
					alt_final_dict[base]=0
					
				for alt_SNV in alt_list:
					if len(ref)==len(alt_SNV) and len(ref)==1:   ##only handle the point mutation(not indel)
						alt_final_list.append(alt_SNV)
						
				##generate the command to generate intermediate sam file
				coor_region=chro + ":" + coor + "-" + coor
				#outfile_sam=bam_infile + ".sam"
				outfile_sam=OUTPUT_PATH + "/" +outfile_name+".sam"
				sam_extract_cmd="samtools view "+bam_infile + " " + coor_region +"> "+  outfile_sam
				subprocess.call(sam_extract_cmd,shell=True)
				
				sam_obj=SAM_File_class(outfile_sam)
				total_correct_coverage=0
				try:
					
					samfile_reader=sam_obj.reader_gen()
					
					for alignment in samfile_reader:
						alignment_flag=int(alignment[sam_obj.FLG_COLUMN])
						if alignment_flag > 1024:   
							pass  ##do not consider the duplicates read
						else:
							alignment_coor=int(alignment[sam_obj.COOR_COLUMN])
							alignment_cigar=alignment[sam_obj.CIGAR_COLUMN]
							alignment_fasta=alignment[sam_obj.SEQ_COLUMN]
							alignment_quality=alignment[sam_obj.QUAL_COLUMN]
							filter_result,nt,false_positive_type=full_filter_retrieve_bases(alignment_fasta,alignment_cigar,alignment_coor,int(coor),alignment_quality,BORDER_FILTER_LENGTH)
							#print filter_result,nt
							if nt==ref:
								alt_final_dict[nt]+=1
								total_correct_coverage+=1
							elif nt=="D":
								#total_correct_coverage+=1
								alt_final_dict[nt]+=1
							elif nt=="N":
								pass
							else:
								total_correct_coverage+=1
								if filter_result==False and count_all==0:
									alt_final_dict[nt]+=1
								if count_all!=0:
									alt_final_dict[nt]+=1
				except:
					pass
								
				PASS_SNV=""
				freq_SNV=""
				test_outcome=""
				#print "Total coverage is", total_correct_coverage
				if total_correct_coverage<MINIMUM_COVERAGE:
					test_outcome="Not_enough_coverage"
					pass
				else:
					for alt in alt_final_list:
						alt_count=alt_final_dict[alt]
						freq=float(alt_count)*100/total_correct_coverage
						if freq > MINIMUM_ALLEL_FRQ:
							if PASS_SNV=="":
								PASS_SNV=alt
								freq_SNV=str(freq)+'%'
							else:
								PASS_SNV=PASS_SNV+","+alt
								freq_SNV=freq_SNV+","+str(freq)+'%'
							test_outcome="true-positive"
							
					##test deletion in the same coordinate
					deletion_count=alt_final_dict["D"]
					freq=int(float(deletion_count)*100/total_correct_coverage)
					if freq>DELETION_THRESHOLD:
						test_outcome="Candidate_Deletion_artifacts"
						PASS_SNV=""
				
				if test_outcome=="":
					test_outcome="Candidate_splice_artifacts"
				
				rows[infile_obj.FILTER_COLUMN]=test_outcome
				
				if test_outcome=="true-positive":
					output_row(outfile_clean_obj.handle,rows)
				
				output_info=""
				for base in "ADTCG":
					base_count=str(alt_final_dict[base])
					base_info=base+"_"+base_count
					if base=="A":
						output_info=base_info
					else:
						output_info=output_info+":"+base_info
				output_row(outfile_obj.handle,rows+[output_info])
				
				rm_sam_cmd="rm "+ outfile_sam
				subprocess.call(rm_sam_cmd,shell=True)
				
			outfile_obj.handle.close()
			outfile_clean_obj.handle.close()
			
				
        ###Start of the program       
	#exit if not enough arguments
	if len(sys.argv) < 3:
		print __doc__
		sys.exit(0)
	
	###set default value
	suffix="vcf"
	bam_suffix="bam"
	bam_infile=None
	infile=None
	sep_char='\t'
	unique_id_length=2
	BORDER_FILTER_LENGTH=6
	prefix="RAA"
	SPLICING_FILTER_LENGTH=10
	Filter_Disable=True
	VCF_PATH=os.getcwd()
	BAM_PATH=os.getcwd()
	OUTPUT_PATH=os.getcwd()
	MINIMUM_COVERAGE=8
	MINIMUM_ALLEL_FRQ=5.0
	DELETION_THRESHOLD=5.0
	MINIMUM_QUAL_SCORE=15
	mode_of_Processing=1
	count_all =0
	count_all_samples = 0
	###get arguments(parameters)
	optlist, cmd_list = getopt.getopt(sys.argv[1:], 'hi:s:a:A:O:c:r:D:d:L:X:E:I:j:Q:r:m:M:F:t:p:l:S:z')
	
	
        for opt in optlist:
        	if opt[0] == '-h':
        		print __doc__; sys.exit(0)
    		elif opt[0] == '-i': infile = opt[1]
    		elif opt[0] == '-I': bam_infile = opt[1]
    		elif opt[0] == '-s': suffix = opt[1]
        	elif opt[0] == '-S': bam_suffix = opt[1]
        	elif opt[0] == '-d': sep_char =opt[1]
        	elif opt[0] == '-L': BORDER_FILTER_LENGTH =int(opt[1])
		elif opt[0] == '-F': SPLICING_FILTER_LENGTH = int(opt[1])
        	elif opt[0] == '-p': prefix = opt[1]
        	elif opt[0] == '-l': unique_id_length = int(opt[1])
		elif opt[0] == '-a': VCF_PATH = opt[1]
		elif opt[0] == '-A': BAM_PATH = opt[1]
		elif opt[0] == '-O': OUTPUT_PATH = opt[1]
		elif opt[0] == '-m': MINIMUM_COVERAGE= int(opt[1])
		elif opt[0] == '-M': MINIMUM_ALLEL_FRQ= float(opt[1])
		elif opt[0] == '-Q': MINIMUM_QUAL_SCORE = int(opt[1])
		elif opt[0] == '-E': mode_of_Processing = int(opt[1])
		elif opt[0] == '-X': count_all = int(opt[1])
		elif opt[0] == '-c': count_all_samples = int(opt[1])
        	else:
        		pass
        
        if infile==None:
        	infiles=extract_files_with_suffix_from_path_to_full_path(VCF_PATH, suffix)
        else:
        	infiles=[infile]
		
	if bam_infile==None:
		bam_infiles=extract_files_with_suffix_from_path_to_full_path(BAM_PATH, bam_suffix)
	else:
		bam_infiles=[bam_infile]
	
	print "infile",infiles
	print "number of infiles is", len(infiles)
	print "bam_infiles",bam_infiles
	print "number of bam infiles is", len(bam_infiles)
    	
	if mode_of_Processing==1:
		combined_input=infiles + bam_infiles
		paired_filelist=generate_paired_files_by_ID(combined_input,unique_id_length,suffix,bam_suffix)
		Filter_Alignment_Artifact(paired_filelist)
		
	elif mode_of_Processing==2:
		Filter_Alignment_Artifact_compiled(infiles,bam_infiles)

	