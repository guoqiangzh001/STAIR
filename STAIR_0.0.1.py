import glob
import os
import subprocess
import shutil
import re

__author__ = 'Zhang Guoqiang; email: zhanggq@big.ac.cn'


def CCSs_correcting(long_reads, short_reads, prefix):
	#utlizing proovread software to error-correct the FLNC CCSs with high confidence short reads (illumina sequening reads such as Hiseq2000. )
	cmd = 'proovread -l {long_reads} '.format(long_reads=long_reads) + '-s {short_reads} '.format(short_reads=short_reads) + '--pre {pre} '.format(pre=prefix)
	if subprocess.check_call(cmd, shell=True) != 0:
		raise SystemCommandError

def CCSs_classify(CCSs, sub_reads, FLNC):
	#obtaining the Full-Length Non-Chimeric (FLNC) CCSs by pbtranscript tools (Pacific Biosciences), only the non-artificial-concatemer CCS with all three signals (poly (A), signal 5' and 3' adaptors) was considered as FLNC.
	cmd='pbtranscript.py classify {ccs} '.format(ccs=CCSs) + ' {fasta} '.format(fasta=sub_reads) + '--min_seq_len {minLen} '.format(minLen=300) + '--cpus {cpus} '.format(cpus=12) + '--flnc {flnc} '.format(flnc=FLNC) + '--nfl {nfl} '.format(nfl='nfl')
	subprocess.check_call(cmd, shell=True)

def genome_gmapDatabase_building(genome_dir, genome_database, genome_seq):
	if not os.path.exists('genome_gmapdir'):
		os.mkdir('genome_gmapdir')
	cmd='gmap_build -D {genome_dir} '.format(genome_dir=genome_dir) + '-d {genome_database} '.format(genome_database=genome_database) + ' {fa} '.format(fa=genome_seq) + '> gmapGD_log'
	if subprocess.check_call(cmd, shell=True) != 0:
		raise SystemCommandError

def genome_StarDatabase_building(genome_dir, genome_seq):
	if not os.path.exists('genome_stardir'):
		os.mkdir('genome_stardir')
	cmd='STAR --runThreadN {T} '.format(T=12) + '--runMode {runMode} '.format(runMode='genomeGenerate') + '--genomeDir {genome_database} '.format(genome_database='genome_stardir') + '--genomeFastaFiles {genomeFasta} '.format(genomeFasta=genome_seq) + '> starGD_log'
	if subprocess.check_call(cmd, shell=True) != 0:
		raise SystemCommandError

def gmapping(genome_dir, genome_database, corr_flnc, sam_out, sam_log):
	cmd='gmap -D {genome_dir} '.format(genome_dir=genome_dir) + '-d {genome_database} '.format(genome_database=genome_database) + '-f {samse} '.format(samse='samse') + '-n {nPaths} '.format(nPaths=0) + '-t {thread} '.format(thread=12) + ' {fa} '.format(fa=corr_flnc) + '> {sam_out} '.format(sam_out=sam_out) + '2> {sam_log} '.format(sam_log=sam_log)
	if subprocess.check_call(cmd, shell=True) != 0:
		raise SystemCommandError

def SJ_decting_pe(genome_dir, short_reads1, short_reads2, outPre,Log):
	cmd='STAR --runThreadN {T} '.format(T=12) + '--genomeDir {genome_database} '.format(genome_database=genome_dir) + '--readFilesIn {short_reads1} {short_reads2} '.format(short_reads1=short_reads1, short_reads2=short_reads2) + '--outFileNamePrefix {outPre} '.format(outPre=outPre) + '--outSAMtype {bam} '.format(bam='BAM SortedByCoordinate') + '--outSAMstrandField {strandField} '.format(strandField='intronMotif') + '--outSAMattributes {SAMattributes} '.format(SAMattributes='NH HI AS nM MD') + '--outFilterMismatchNmax {MismatchNmax} '.format(MismatchNmax=10) + '--outFilterIntronMotifs {IntronMotifs} '.format(IntronMotifs='RemoveNoncanonical') + '--alignIntronMax 200000 {IntronMax} '.format(IntronMax=200000) + '--alignMatesGapMax {GapMax} '.format(GapMax=200000) + '--chimSegmentMin {chimSegmentMin} '.format(chimSegmentMin=15) + '--chimJunctionOverhangMin {chimJunctionOverhangMin} '.format(chimJunctionOverhangMin=15) + '> {log}'.format(log=Log)
	if subprocess.check_call(cmd, shell=True) != 0:
		raise SystemCommandError

def SJ_decting_se(genome_dir, short_reads, outPre,Log):
	cmd='STAR --runThreadN {T} '.format(T=12) + '--genomeDir {genome_database} '.format(genome_database=genome_dir) + '--readFilesIn {short_reads} '.format(short_reads=short_reads) + '--outFileNamePrefix {outPre} '.format(outPre=outPre) + '--outSAMtype {bam} '.format(bam='BAM SortedByCoordinate') + '--outSAMstrandField {strandField} '.format(strandField='intronMotif') + '--outSAMattributes {SAMattributes} '.format(SAMattributes='NH HI AS nM MD') + '--outFilterMismatchNmax {MismatchNmax} '.format(MismatchNmax='10') + '--outFilterIntronMotifs {IntronMotifs} '.format(IntronMotifs='RemoveNoncanonical') + '--alignIntronMax 200000 {IntronMax} '.format(IntronMax=200000) + '--alignMatesGapMax {GapMax} '.format(GapMax=200000) + '--chimSegmentMin {chimSegmentMin} '.format(chimSegmentMin=15) + '--chimJunctionOverhangMin {chimJunctionOverhangMin} '.format(chimJunctionOverhangMin=15) + '> {log}'.format(log=Log)
	if subprocess.check_call(cmd, shell=True) != 0:
		raise SystemCommandError

class SystemCommandError(Exception):
	"""System Command Error Class"""
	pass

def gff(gff_file):
	#Chr1	PacBio	transcript	1224015	1229923	.	-	.	gene_id "PB.29"; transcript_id "PB.29.2";
	#Chr1	PacBio	exon	1224015	1224326	.	-	.	gene_id "PB.29"; transcript_id "PB.29.2";
	#Chr1	PacBio	exon	1225054	1225269	.	-	.	gene_id "PB.29"; transcript_id "PB.29.2";
	#Chr1	PacBio	exon	1225578	1225662	.	-	.	gene_id "PB.29"; transcript_id "PB.29.2";
	introns={}
	strnd={}
	exons={}
	id_pre=''
	Gene={}
	transcripts={}
	for line in gff_file:
		line=line.rstrip().split('\t')
		if len(line) > 2:
			if line[2] == 'gene':
				id=line[-1].split(';')[1].split('=')[1]
				Gene[id]=[line[0], int(line[3]), int(line[4]), line[6]]
			elif line[2] == 'transcript':
				intron=[]
				exon=[]
				id=line[-1].split(';')[1].split(' ')[-1].replace('"','')
				transcripts[id]=[line[0], int(line[3]), int(line[4]), line[6]]
			elif line[2] == 'exon':
				id==line[-1].split('"')[1]
				exon.append([line[0], line[3], line[4]])
				if id_pre == id:
					intrn=(line[0], int(exn[2])+1, int(line[3])-1)
					#intrn=(line[0], int(exn[2])+1, int(line[3])-1, line[6])
					intron.append(intrn)
					introns[id]=intron
				id_pre=id
				exn=[line[0], line[3], line[4]]
				exons[id]=exon
				strnd[id]=line[6]
	return Gene, transcripts, exons, introns, strnd

def splicingJunction_validating(splicingJunction_file, dp):
	junCov={}
	strnd={'0':'NA', '1':'+', '2':'-'}
	for line in splicingJunction_file.readlines():
		line=line.rstrip().split()
		if int(line[6]) > dp:
			id=(line[0], int(line[1]), int(line[2]))
			#id=(line[0], int(line[1]), int(line[2]), strnd[line[3]])
			junCov[id]=line[6]
	return junCov

def isoforms_validating(juns,covs):
	Cov={}
	noAll_val={}
	for k,v in juns.items():
		cov=[]
		for i in range(len(v)):
			if v[i] in covs:
				cov.append((v[i],covs[v[i]]))
			else:
				cov.append((v[i], 'NA'))
				noAll_val[k]='noAll validated'
		Cov[k]=cov
	allVal={}
	for k1,v1 in Cov.items():
		if k1 not in noAll_val:
			allVal[k1]=juns[k1]
	return allVal

def files_merging(raw_files,mergedFile):
	merging=['cat {file} '.format(file=file) + '>> {mergedFile}'.format(mergedFile=mergedFile) for file in raw_files]
	for Merge in merging:
		os.system(Merge)

def fa_reader(seq_file):
	Seq={}
	Name={}
	for line in seq_file:
		#>PB.1.1|Chr1:16263-20316(+)|Chr1_c533/f10p6/2504
		#aagcaaaaggactctgtattttctgtgttgtctcctctagtcctctcctcccttctgtTC
		#>PB.3.1|Chr10:4639803-4644295(-)|m140625_210453_42213_c100631252550000001823104809121442_s1_p0/26872/31_1658_CCS
		#ATCCATCGAAACCCCACTTTATTAGCACCAACCCCGTCTCCTTCCTCCATTCCTCCGCTG
		#CTAGAACCTTCTAGAAGCTCTCCTCCCTCCCATGGCGGCGACGATCGTGTCCGTCAAGGC
		#GCGCCAGATCTTCGACAGCAGGGGCAACCCACCGTCGAGGTCGATGTGTGCTGCTCAGAT
		#Line=line.rstrip()
		if '>' in line:
			Id=line.replace('>','').split('|')[0]
			Name[Id]=line
		else:
			if Id in Seq:
				Seq[Id] +=line
			else:
				Seq[Id]=line
	return Name, Seq

def accur(Read):
	if 'NM:i:' in Read:
		Len=len(Read.split('\t')[9])
		MisMatch=int(re.match(".+NM:i:(\d+)", Read).group(1))
		Accur=MisMatch/Len
	return Accur

def sam_filter_sort(ori_sam, processed_sam):
	filtering_sam=open('filtering_sam.sam','w')
	for line in open(ori_sam, 'r'):
		if line[0]=='@':
			filtering_sam.write('{0}'.format(line))
		else:
			if int(line.split('\t')[4]) >= 30 and accur(line) < 0.15:
				filtering_sam.write('{0}'.format(line))
	filtering_sam.close()
	cmd='sort -k 3,3 -k 4,4n {0} > {1}'.format('filtering_sam.sam', processed_sam)
	if subprocess.check_call(cmd, shell=True) != 0:
		raise SystemCommandError
	return processed_sam

def fq2fa(fq_file,fa_file):
	fq = open(fq_file,'r')
	fa = open(fa_file,'w')
	i=1
	line_i = 0
	for line in fq:
		if line_i == 0:
			if line[0] !='@':
				print("Err: invalid LR fastq format")
			fa.write('>' + line[1:-1] + '\n')
			i+=1
			line_i = 1
		elif line_i == 1:
			fa.write(line)
			line_i = 2
		elif line_i == 2:
			line_i = 3
		elif line_i == 3:
			line_i = 0

def main(args):
	#filesList=os.listdir(os.getcwd())
	Cwd=os.getcwd()
	if args.temp:
		temp_dir=args.temp
	else:
		temp_dir=os.getcwd() + '/temp'
	
	if not os.path.exists(temp_dir):
		os.mkdir(temp_dir)
	os.chdir(temp_dir)

	if args.reads:
		if '/' in args.reads:
			Reads=args.reads
		else:
			Reads= Cwd + '/' + args.reads
	else:
		if args.read1:
			if '/' in args.read1:
				Reads1=args.read1
				Reads2=args.read2
			else:
				Reads1=Cwd + '/' + args.read1
				Reads2=Cwd + '/' + args.read2

	if '/' in args.ref:
		Ref=args.ref
	else:
		Ref= Cwd + '/' + args.ref

	if '/' in args.input:
		Input=args.input
	else:
		Input= Cwd + '/' + args.input

	if args.subreads:
		if '/' in args.subreads:
			Subreads=args.subreads
		else:
			Subreads= Cwd + '/' + args.subreads

	if args.gmapdir:
		if '/' in args.gmapdir:
			Gmapdir=args.gmapdir
		else:
			Gmapdir= Cwd + '/' + args.gmapdir
		if '/' in args.gmapdb:
			Gmapdb=args.gmapdb.split('/')[-1]
		else:
			Gmapdb=args.gmapdb

	if '/' in args.stardir:
		Stardir=args.stardir
	else:
		Stardir= Cwd + '/' + args.stardir

	#print(os.getcwd())
	if args.SJ:
		sj_file=args.SJ
	else:
		if not os.path.exists(args.o+'.SJ.out.tab'):
			Pre=args.o+'.'
			if not args.stardir:
				genome_StarDatabase_building('genome_stardir', Ref)
				if args.reads:
					SJ_decting_se('genome_stardir', Reads, Pre, 'RNASeq_starmap.log')
				else:
					SJ_decting_pe('genome_stardir', Reads1, Reads2, Pre, 'RNASeq_starmap.log')
			else:
				if args.reads:
					SJ_decting_se(Stardir, Reads, Pre, 'RNASeq_starmap.log')
				else:
					SJ_decting_pe(Stardir, Reads1, Reads2, Pre, 'RNASeq_starmap.log')
			sj_file=Pre+'SJ.out.tab'
		else:
			sj_file=args.o+'.SJ.out.tab'
	sjVal=splicingJunction_validating(open(sj_file),1)

	if not os.path.exists('flnc_ccs_corrected.untrimmed.fa'):
		'generating the FLNC CCSs and correcting them when the input is in fofn format which provides the pacbio raw reads files'
		if args.classify=='Y':
			#'classifying the FLNC CCSs into FLNC and NFL group'
			CCSs_classify(Input, Subreads, 'flnc_ccs.fasta')
			FLNC='flnc_ccs.fasta'
		else:
			FLNC=Input

		if args.corr=='Y':
			if os.path.exists('flnc_ccs_corrected'):
				os.system('rm -rf flnc_ccs_corrected')
			if args.read1 and args.read2:
				if not os.path.exists('mergedPE_reads.fasta'):
					pe_files=[Reads1, Reads2]
					SR='mergedPE_reads.fasta'
					mergedPE=files_merging(pe_files, SR)
				SR='mergedPE_reads.fasta'
			elif args.reads:
				SR=Reads
			Prefix='flnc_ccs_corrected'
			CCSs_correcting(Input, SR, Prefix)
			corr_flnc_fq=os.getcwd() + '/flnc_ccs_corrected/flnc_ccs_corrected.untrimmed.fq'
			corr_flnc_fa='flnc_ccs_corrected.untrimmed.fa'
			fq2fa(corr_flnc_fq,corr_flnc_fa)
		else:
			if args.input.split('.')[-1]=='fa':
				corr_flnc_fa=Input
			else:
				corr_flnc_fa='flnc_ccs.fa'
				fq2fa(Input,corr_flnc_fa)
	else:
		corr_flnc_fa='flnc_ccs_corrected.untrimmed.fa'

	if args.LR_sam:
		if "_processed.sam" in args.LR_sam:
			processed_sam=args.LR_sam
		else:
			processed_sam='flnc_processed.sam'
			sam_filter_sort(args.LR_sam, processed_sam)
	else:
		if not os.path.exists('flnc_processed.sam'):
			sam_out='flnc_ccs_corrected.sam'
			sam_log='flnc_ccs_corrected.sam.log'
			if args.gmapdir and args.gmapdb:
				gmapping(Gmapdir, Gmapdb, corr_flnc_fa, sam_out, sam_log)
			else:
				genome_gmapDatabase_building('genome_gmapdir', 'genome_gmapdb', Ref)
				gmapping('genome_gmapdir', 'genome_gmapdb', corr_flnc_fa, sam_out, sam_log)
			processed_sam='flnc_processed.sam'
			sam_filter_sort(sam_out, processed_sam)
		else:
			processed_sam='flnc_processed.sam'

	cmd='collapse_isoforms_by_sam.py --input {fa} '.format(fa=corr_flnc_fa) + '-s {sam} '.format(sam=processed_sam) + '-o {outPrefix} '.format(outPrefix=args.o + '.ori')
	if subprocess.check_call(cmd, shell=True) != 0:
		raise SystemCommandError
	else:
		gff_file=args.o + '.ori.collapsed.gff'
		fa_file=args.o + '.ori.collapsed.rep.fa'
		pb_isoforms=gff(open(gff_file))
		pb_isoforms_val=isoforms_validating(pb_isoforms[3],sjVal)
		val_gff=open(args.o + '.val.collapsed.gff','w')
		val_fa=open(args.o + '.val.collapsed.fa','w')
		pb_isoforms_fa=fa_reader(open(fa_file))
		singleExon_gff=open(args.o + '.singleExon.collapsed.gff','w')
		singleExon_fa=open(args.o + '.singleExon.collapsed.fa','w')
		Val={}
		Single={}
		for line in open(gff_file,'r'):
			#'Chr1	PacBio	transcript	16263	20316	.	+	.	gene_id "PB.1"; transcript_id "PB.1.1";
			#Chr1	PacBio	exon	16263	16976	.	+	.	gene_id "PB.1"; transcript_id "PB.1.1";
			T_id=line.rstrip(';\n').split('transcript_id ')[-1].replace('"','')
			if T_id in pb_isoforms_val:
				val_gff.write('{0}'.format(line))
				Val[T_id]='{0}{1}'.format(pb_isoforms_fa[0][T_id], pb_isoforms_fa[1][T_id])
			if T_id not in pb_isoforms[3]:
				singleExon_gff.write('{0}'.format(line))
				Single[T_id]='{0}{1}'.format(pb_isoforms_fa[0][T_id], pb_isoforms_fa[1][T_id])
		for k,v in Val.items():
			val_fa.write('{0}'.format(Val[k]))
		for k,v in Single.items():
			singleExon_fa.write('{0}'.format(Single[k]))
		for f in glob.glob('*singleExon.collapsed*'):
			shutil.move(f,Cwd)
		for f in glob.glob('*val.collapsed*'):
			shutil.move(f,Cwd)

if __name__ == "__main__":
	from argparse import ArgumentParser
	parser = ArgumentParser(description='collapsing the error-corrected FLNC CCSs into splicing isoforms and validated them by checking their chained splicing junctions with high confidence short reads')
	parser.add_argument('-input', required=True, help='a FA/FQ filename, a required fasta/fastq file contains the FLNC CCSs generated by RS_IsoSeq protocol (default is fasta).')
	parser.add_argument('-temp', help='Temp folder (default is "temp" under current directory).')
	parser.add_argument('-corr', help='"Y/N" refers to correct the FLNC CCSs with high confidence short reads or not.')
	parser.add_argument('-classify', help='"Y/N" refers to classify the CCSs into FLNC and NFL groups or not, the option "-subreads" (the sequences of subreads) must be provided in the case of "Y".')
	parser.add_argument('-subreads', help='the sequences of subreads (default is fasta) must be provided when you choose to classify the CCSs into FLNC and NFL groups.')
	parser.add_argument('-stardir', help='the directory where the genome indices are stored for star mapping.')
	parser.add_argument('-gmapdir', help='the directory where the genome indices are stored for gmap mapping.')
	parser.add_argument('-gmapdb', help='the database being used for gmap mapping.')
	parser.add_argument('-SJ', help='splicing junctions in the format are produced with STAR software.')
	parser.add_argument('-LR_sam', help='"sam filename", alignments of pacbio reads against reference genome with GMAP software in sam format (raw or processed alignments in suffix ".sam" or "_processed.sam", respectively).')
	parser.add_argument('-reads', help='"FA/FQ filename", High confidence short single-end reads (default is fasta).')
	parser.add_argument('-read1', help='"FA/FQ filename", High confidence pair-end reads1 (default is fasta).')
	parser.add_argument('-read2', help='"FA/FQ filename", High confidence pair-end reads2 (default is fasta).')
	parser.add_argument('-ref', required=True, help='"a FA filename, required", reference genome squences.')
	parser.add_argument('-o', required=True, help='an Output filename prefix, required.')
	args = parser.parse_args()
	main(args)
