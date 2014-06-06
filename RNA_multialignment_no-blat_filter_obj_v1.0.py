# -*- coding: cp936 -*-
"""
The main purpose of this script is to filter the RNA variant
call due to multiple alignment in the RNASeq data WITHOUT doing
realignment.

=============================
Usage: python RNA_multialignment_filter_obj_v0.0.py
-h help

-i vcf files to processed                       *[No default value]

-I bam input file paired with vcf file		*[No default value]

-a vcf input file path 				[default value:current folder]

-A bam input file path 				[default value:current folder]

-p prefix of output files                       [default value:"Remove_blat_artifact(RBA)"]

-s suffix for the files to be processed         [default value "vcf"]

-S suffix for the bam files to be processed     [default value "bam"]

-d sep_char within among columns                [default value '\t']

-l unique_id length for the infile		[default value 2]

-L minimum multiple aliged read frequency to discard  	[default value 5]

-P pair vcf and bam files under the current folder [default value "Yes"]

-E mode of processing				[default value 1]
	1 means one single sample-vcf vs one single sample bam
	2 means one compiled multi sample vcf vs multiple single sample bam

===================
input description:
input files:
1. filtered RNA editing VCF file or mutation call VCF files
2. bam file(s)

======================
output files:
subset of the input file
the format of input file will be the same
the filtered one, the filtered one will be discarded directly
the detailed information one: filtered column will be marked as "Multiple_alignment_Artifact": and info in the additional column

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
##Version: 1.0
##Compatible Python Version:2.4 or above
##Developed at: 2013-4-23
##Last Update: 2014-6-6
##Status: developing properly

###Code Framework

if __name__ == "__main__":
	###Module Import	
	import sys, csv, getopt, re
	import os
	import math
	from itertools import ifilter
	import subprocess
	

	##
	SAM_FLG_COLUMN=1
	SAM_READGROUP_COLUMN=12
	SAM_CHR1_COLUMN=2
	SAM_CHR2_COLUMN=6
	SAM_COOR_COLUMN=3
	SAM_CIGAR_COLUMN=5
	SAM_SEQ_COLUMN=9
	SAM_LEN_COLUMN=8
	SAM_QUALITY_COLUMN=10
	
	##VCF common constant
	VCF_CHRO_COLUMN=0
	VCF_COOR_COLUMN=1
	VCF_ALT_COLUMN=4
	VCF_REF_COLUMN=3
	VCF_INFO_COLUMN=7
	VCF_FILTER_COLUMN=6
	VCF_INFO_SEP=';'
	VCF_SAMPLE_COLUMN=9
	
	#from Constant_Library import *
	#from General_Library import *
	#from File_Class import *
	#from Sequencing_Library import *
	#from Specific_File_Class import MISO_class
	
	digit_list=["0","1","2","3","4","5","6","7","8","9"]
	full_human_chro_list=["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",\
		 "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM"]
	
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
			#self.count_column_number()
			
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
				
		
		def sampleID_gen(self):
			
			##version2.0
			if self.name_only.count("/")>0:
				tmp_list=self.name_only.split("/")
				self.name_only=tmp_list[-1]
			tmp_list=(self.name_only).split(self.FILENAME_SPLIT_CHAR)
			sampleID=tmp_list[self.SAMPLE_ID_POS]
			return sampleID
		
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
		
		def unique_ID_list_gen(self,reader,unique_ID_column):
			unique_ID_list=[]
			for row in reader:
				unique_ID=row[unique_ID_column]
				unique_ID_list.append(unique_ID)
			return unique_ID_list
		
		def unique_ID_list_gen_v2(self,reader,unique_ID_column):
			## under development
			unique_ID_list=[]
			for row in reader:
				unique_ID=row[unique_ID_column]
				unique_ID_list.append(unique_ID)
			return unique_ID_list
			
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
	
	
	OUTPUT_SEP_CHAR='\t'
	
	### Specific Functions definiation
	def test_100P(read_count_info):
		read_count_list=read_count_info.split(":")
		read_count_data=read_count_list[1]
		read_count_ref,read_count_alt=read_count_data.split(",")
		read_count_ref=int(read_count_ref)
		
		#read_count_alt=int(read_count_alt)
		total_count=int(read_count_list[2])
		freq=int(100.00*read_count_ref/total_count)
		if read_count_ref==0:
			return True
		elif read_count_ref==1 and total_count>=10:
			return True
		elif read_count_ref==1 and total_count<10:
			return False
		elif total_count<10 and read_count_ref>=2:
			return False
		elif total_count>=10 and freq>=5:
			return False
		else:
			return True
			
	def relaxed_P1_test(raw_vcf_data):
		## this is a less stringent criteria
		## this is based on not 100 percent artifacts
		result=False
		
		freq_list=raw_vcf_data[VCF_SAMPLE_COLUMN:]
		data_count=0
		data_100P_count=0
		for freq in freq_list:
			if freq==".":
				pass
			elif freq=="./.":
				pass
			elif freq.count(":")==2:
				##Process the data
				data_count+=1
				if test_100P(freq):
					data_100P_count+=1
			else:
				pass
		
		if data_count<=2 and data_100P_count==0:
			result=True
		elif data_count>=3 and data_100P_count<=1:
			result=True
		else:
			pass
		return result
	
	def Filter_multiple_alignment_Artifact_compiled(infiles,bam_infiles):
		##
		##
		print "process the 1 vs n type"
		for infile in infiles:
			
			print "Processing infile:", infile
			##Set up vcf infile object
			infile_obj=VCF_File_class(infile)
			infile_obj.SAMPLE_ID_LEN=unique_id_length
			infile_reader=infile_obj.reader_gen()
			
			sample_list=infile_obj.sample_list_gen()
			sample_COUNT=len(sample_list)
			sample_dict=dict()
			
			for sample_ID in sample_list:   ##link the bam file with vcf file
				sample_ID_info=sample_ID.split(".")  ##standard of combine vcf results
				sample_ID_simple=sample_ID_info[0]
				for bam_file in bam_infiles:
					if bam_file.count(sample_ID_simple)==1:
						sample_dict[sample_ID]=bam_file
						break
			
			pair_count=len(sample_dict.keys())  ##qc of the linked bam - vcf file
			if pair_count!=sample_COUNT:
				print "Error: Incorrect Sample Pairing"
				print "paired sample found is", pair_count
				print "total sample number is", sample_COUNT
				sys.exit(0)
			
			
			### output file standard output
			final_prefix=prefix+"_"+str(FILTER_MULTI_PERCENTAGE)
			detailed_prefix=final_prefix+"_details"
			outfile_name=infile_obj.outputfilename_gen(detailed_prefix,"VCF")
			outfile_clean_name=infile_obj.outputfilename_gen(final_prefix,"vcf")
			
			outfile_obj=VCF_File_class(outfile_name)
			outfile_clean_obj=VCF_File_class(outfile_clean_name)
			
			outfile_obj.output_handle_gen(infile,os.getcwd(),sample_list)
			outfile_clean_obj.output_handle_gen(infile)
			
			
			## data process
			for rows in infile_reader:
				chro=rows[infile_obj.CHRO_COLUMN]
				coor=rows[infile_obj.COOR_COLUMN]
				ref=rows[infile_obj.REF_COLUMN]
				alt=rows[infile_obj.ALT_COLUMN]
				#print chro,coor,ref,alt
				if len(alt)==1 and len(ref)==1:
					
					data_all=rows[infile_obj.SEP_INFO_START_COLUMN:]
					process_dict=dict()  ## indicator whether was covered in vcf or not
					result_list=[]    
					test_dict=dict()     ## 
					
					## function to figure out which one to process
					for index in range(sample_COUNT):
						sample_ID=sample_list[index]
						data=data_all[index]
						if data=="./.":
							process_dict[sample_ID]=0
							
							test_dict[sample_ID]=["NA"]
						else:
							process_dict[sample_ID]=1
							
							test_dict[sample_ID]=[]
					
					for index in range(sample_COUNT):
						sample_ID=sample_list[index]
						if process_dict[sample_ID]==1:
							bam_infile=sample_dict[sample_ID]
							
							##multiple chromosome pair indicator of multi
							## alignment issue
							human_chro_dict=dict()
							for chromosome in full_human_chro_list:
								human_chro_dict[chromosome]=0
							
							##generate the command to generate intermediate sam file
							coor_region=chro + ":" + coor + "-" + coor
							#outfile_sam=bam_infile + ".sam"
							outfile_sam=outfile_name+".sam"
							sam_extract_cmd="samtools view "+ bam_infile + " " + coor_region +"> "+ outfile_sam
							subprocess.call(sam_extract_cmd,shell=True)
						
							#sam_obj=GeneralFile_class(outfile_sam)
							samfile_reader=csv.reader(open(outfile_sam,'rU'),delimiter='\t',quoting=csv.QUOTE_NONE)
							
							total_alt_coverage=0
							total_multi_coverage=0
							diff_chr_count=0
							total_all_coverage=0
							same_chr_count=0
							
							for alignment in samfile_reader:
								alignment_flag=int(alignment[SAM_FLG_COLUMN])
								if alignment_flag > 1024:
									pass
								else:
									
									alignment_coor=int(alignment[SAM_COOR_COLUMN])
									alignment_cigar=alignment[SAM_CIGAR_COLUMN]
									alignment_fasta=alignment[SAM_SEQ_COLUMN]
									alignment_quality=alignment[SAM_QUALITY_COLUMN]
									filter_decision,coor_base,false_positive_type=full_filter_retrieve_bases(alignment_fasta,alignment_cigar,alignment_coor,int(coor),alignment_quality)
									
									if coor_base!="N":
										##ignore the splicing reads
										right_chr=alignment[SAM_CHR2_COLUMN]
										alignment_fragment_length=int(alignment[SAM_LEN_COLUMN])
										if right_chr=="=" and abs(alignment_fragment_length)>LONG_DISTANCE_CUTOFF:
											same_chr_count+=1
										if right_chr=="*" and abs(alignment_fragment_length)>LONG_DISTANCE_CUTOFF:
											same_chr_count+=1
										total_all_coverage+=1
										
									if coor_base==alt:	
										
										#left_chr=alignment[SAM_CHR1_COLUMN]
										
										if len(right_chr)>1:
											#diff_chr_count+=1
											human_chro_dict[right_chr]+=1

										alignment_info=alignment[SAM_READGROUP_COLUMN:]
										total_alt_coverage+=1
										for info in alignment_info:
											if info.count("NH:i:")==1:
												alignment_count=int(info[5:])
												break
										if alignment_count>1:
											total_multi_coverage+=1
										
										
							if total_multi_coverage<=2:
								total_multi_coverage=0
							chro_count_list=human_chro_dict.values()
							diff_chr_count=max(chro_count_list)
							#print "same_chr_count",same_chr_count
							#print total_
							
							if diff_chr_count<=2:
								diff_chr_count=0
							if same_chr_count<=2:
								same_chr_count=0
							##this eliminate the small sample random issue in here.
							
							#print diff_chr_count	
							if total_alt_coverage==0:
								result_list.append("NA")
							else:
								multi_alignment_artifact=False  ## consider it to be true positive until showed otherwise
								multi_percentage=int(100.00*total_multi_coverage/total_alt_coverage)
								diff_chr_percentage=int(100.00*diff_chr_count/total_all_coverage)
								same_chr_percentage=int(100.00*same_chr_count/total_all_coverage)
								if multi_percentage>=FILTER_MULTI_PERCENTAGE or diff_chr_percentage>FILTER_MULTI_PERCENTAGE or same_chr_percentage>FILTER_MULTI_PERCENTAGE:
									multi_alignment_artifact=True
								all_data=str(multi_percentage)+"_"+str(diff_chr_percentage)+"_"+str(same_chr_percentage)	
								#result_dict[bam_infile].append(all_data)
								if multi_alignment_artifact:
									result_list.append("false-positive")
	
								else:								
									result_list.append("true-positive")
						
					#result_list=result_dict.values()
					result_true_count=result_list.count("true-positive")
					result_false_count=result_list.count("true-positive")
					coveraged_count=result_true_count+result_false_count
					if result_false_count==0:
						filter_append="true-positive"
					elif result_false_count==1 and result_true_count>=1:
						filter_append="likely-false-positive"
					elif result_false_count==1 and result_true_count==0:
						filter_append="false-positive"
					elif result_false_count>=2:
						filter_append="false-positive"
					else:
						pass
					
					
					##if relaxed criteria will be applied
					if filter_append.count("false-positive")==0:
						output_row(outfile_clean_obj.handle,rows)
						
					elif relaxed_method==0:
						filter_result=rows[infile_obj.FILTER_COLUMN]
						rows[infile_obj.FILTER_COLUMN]=filter_result+";No-blat-test:"+filter_append
						
					elif relaxed_method==1:
						sample_data=rows[infile_obj.SEP_INFO_START_COLUMN:]
						candidate=relaxed_P1_test(sample_data)
						if candidate==True:
							output_row(outfile_clean_obj.handle,rows)
							filter_result=rows[infile_obj.FILTER_COLUMN]
							rows[infile_obj.FILTER_COLUMN]=filter_result+";No-blat-test:relaxed"
						else:
							filter_result=rows[infile_obj.FILTER_COLUMN]
							rows[infile_obj.FILTER_COLUMN]=filter_result+";No-blat-test:false-positive"
							
					output_row(outfile_obj.handle,rows)  ## output
						
					##Round 2: testing the multiple-alignment issue directly
					
			outfile_obj.handle.close()
			outfile_clean_obj.handle.close()
		
	def Filter_multiple_alignment_Artifact(paired_files):
		
		##Section I: Generate the gene annotation dictionary
		
			
		print "process the 1 vs 1 type"	
		for infile,bam_infile in paired_files:
			print "Processing infile:", infile
			##Set up vcf infile object
			infile_obj=VCF_File_class(infile)
			infile_obj.SAMPLE_ID_LEN=unique_id_length
			infile_reader=infile_obj.reader_gen()
			
			final_prefix=prefix+"_"+str(FILTER_MULTI_PERCENTAGE)
			outfile_name=infile_obj.outputfilename_gen(final_prefix,"VCF")
			outfile_clean_name=infile_obj.outputfilename_gen(final_prefix,"vcf")
			
			outfile_obj=VCF_File_class(outfile_name)
			outfile_clean_obj=VCF_File_class(outfile_clean_name)
			
			outfile_obj.output_handle_gen(infile)
			outfile_clean_obj.output_handle_gen(infile)
			
			for rows in infile_reader:
				chro=rows[infile_obj.CHRO_COLUMN]
				coor=rows[infile_obj.COOR_COLUMN]
				ref=rows[infile_obj.REF_COLUMN]
				alt=rows[infile_obj.ALT_COLUMN]
				#print chro,coor,ref,alt
				if len(alt)==1 and len(ref)==1:
					multi_alignment_artifact=False
					result="NA"
	
					human_chro_dict=dict()
					for chromosome in full_human_chro_list:
						human_chro_dict[chromosome]=0
					
					##generate the command to generate intermediate sam file
					coor_region=chro + ":" + coor + "-" + coor
					#outfile_sam=bam_infile + ".sam"
					outfile_sam=outfile_name+".sam"
					sam_extract_cmd="samtools view "+ bam_infile + " " + coor_region +"> "+ outfile_sam
					subprocess.call(sam_extract_cmd,shell=True)
				
					#sam_obj=GeneralFile_class(outfile_sam)
					samfile_reader=csv.reader(open(outfile_sam,'rU'),delimiter='\t',quoting=csv.QUOTE_NONE)
					total_alt_coverage=0
					total_multi_coverage=0
					diff_chr_count=0
					total_all_coverage=0
					same_chr_count=0
					for alignment in samfile_reader:
						
						alignment_flag=int(alignment[SAM_FLG_COLUMN])
						if alignment_flag > 1024:
							pass
						else:
							
							alignment_coor=int(alignment[SAM_COOR_COLUMN])
							alignment_cigar=alignment[SAM_CIGAR_COLUMN]
							alignment_fasta=alignment[SAM_SEQ_COLUMN]
							alignment_quality=alignment[SAM_QUALITY_COLUMN]
							filter_decision,coor_base,false_positive_type=full_filter_retrieve_bases(alignment_fasta,alignment_cigar,alignment_coor,int(coor),alignment_quality)
							
							if coor_base!="N":
								##ignore the splicing reads
								right_chr=alignment[SAM_CHR2_COLUMN]
								alignment_fragment_length=int(alignment[SAM_LEN_COLUMN])
								if right_chr=="=" and abs(alignment_fragment_length)>5000:
									same_chr_count+=1
								if right_chr=="*" and abs(alignment_fragment_length)>5000:
									same_chr_count+=1
								total_all_coverage+=1
								
							if coor_base==alt:	
								
								#left_chr=alignment[SAM_CHR1_COLUMN]
								
								if len(right_chr)>1:
									#diff_chr_count+=1
									human_chro_dict[right_chr]+=1

								alignment_info=alignment[SAM_READGROUP_COLUMN:]
								total_alt_coverage+=1
								for info in alignment_info:
									if info.count("NH:i:")==1:
										alignment_count=int(info[5:])
										break
								if alignment_count>1:
									total_multi_coverage+=1
									
									
						if total_multi_coverage<=2:
							total_multi_coverage=0
						chro_count_list=human_chro_dict.values()
						diff_chr_count=max(chro_count_list)
						#print "same_chr_count",same_chr_count
						#print total_
						
						if diff_chr_count<=2:
							diff_chr_count=0
						if same_chr_count<=2:
							same_chr_count=0
						##this eliminate the small sample random issue in here.
						
						#print diff_chr_count	
						if total_alt_coverage==0:
							result="NA"
						else:		
							multi_percentage=int(100.00*total_multi_coverage/total_alt_coverage)
							diff_chr_percentage=int(100.00*diff_chr_count/total_all_coverage)
							same_chr_percentage=int(100.00*same_chr_count/total_all_coverage)
							if multi_percentage>=FILTER_MULTI_PERCENTAGE or diff_chr_percentage>FILTER_MULTI_PERCENTAGE or same_chr_percentage>FILTER_MULTI_PERCENTAGE:
								multi_alignment_artifact=True
							all_data=str(multi_percentage)+"_"+str(diff_chr_percentage)+"_"+str(same_chr_percentage)
							
							if multi_alignment_artifact:
								result="false-positive"
							else:
								result="true-positive"
					
				
				
				if result=="true-positive":
					output_row(outfile_clean_obj.handle,rows)
				
				filter_result=rows[infile_obj.FILTER_COLUMN]
				rows[infile_obj.FILTER_COLUMN]=filter_result+";No-blat-test:"+result
	
				output_row(outfile_obj.handle,rows)
					
			
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
	reference="Linga"
	unique_id_length=2
	BORDER_FILTER_LENGTH=6
	prefix="RBA_with_retrieve"
	FILTER_MULTI_PERCENTAGE=5
	#FILTER_DIFF_CHRO_PERCENTAGE=10
	LONG_DISTANCE_CUTOFF=5000
	infile_path=os.getcwd()
	bam_path=os.getcwd()
	mode_of_Processing=1
	relaxed_method = 1
	MINIMUM_QUAL_SCORE=15
	###get arguments(parameters)
	optlist, cmd_list = getopt.getopt(sys.argv[1:], 'hi:a:A:s:r:d:L:I:j:c:b:B:r:F:t:p:l:S:E:R:z')

        for opt in optlist:
        	if opt[0] == '-h':
        		print __doc__; sys.exit(0)
    		elif opt[0] == '-i': infile = opt[1]
    		elif opt[0] == '-I': bam_infile = opt[1]
		elif opt[0] == '-a': infile_path = opt[1]
		elif opt[0] == '-A': bam_path = opt[1]
    		elif opt[0] == '-s': suffix = opt[1]
        	elif opt[0] == '-S': bam_suffix = opt[1]
        	elif opt[0] == '-d': sep_char =opt[1]
        	elif opt[0] == '-L': FILTER_MULTI_PERCENTAGE =int(opt[1])
        	elif opt[0] == '-r': reference = opt[1]
        	elif opt[0] == '-p': prefix = opt[1]
        	elif opt[0] == '-l': unique_id_length = int(opt[1])
		elif opt[0] == '-E': mode_of_Processing = int(opt[1])
		elif opt[0] == '-c' : LONG_DISTANCE_CUTOFF = int(opt[1])
		elif opt[0] == '-R' : relaxed_method = int(opt[1])
        	else:
        		pass
        
        if infile==None:
        	infiles=extract_files_with_suffix_from_path_to_full_path(os.getcwd(), suffix)
        else:
        	infiles=[infile]
	if bam_infile==None:
		bam_infiles=extract_files_with_suffix_from_path_to_full_path(bam_path, bam_suffix)
		print bam_infiles
	else:
		bam_infiles=[bam_infile]
	print "infile",infiles
	print "number of infiles is", len(infiles)
	print "bam_infiles",bam_infiles
	print "number of bam infiles is", len(bam_infiles)
    	
	#combined_input=infiles + bam_infiles
    	#paired_filelist=generate_paired_files_by_ID(combined_input,unique_id_length,suffix,bam_suffix)
	
	if mode_of_Processing==1:
		combined_input=infiles + bam_infiles
		paired_filelist=generate_paired_files_by_ID(combined_input,unique_id_length,suffix,bam_suffix)
		Filter_multiple_alignment_Artifact(paired_filelist)
		
	elif mode_of_Processing==2:
		Filter_multiple_alignment_Artifact_compiled(infiles,bam_infiles)
		

	