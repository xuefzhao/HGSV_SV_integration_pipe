#############################################################################################################
##################read in vcf info from different algorithms, keep del and dup only for now##################
#############################################################################################################
#############################################################################################################

###################couple global variables###################
import os
import numpy
import time
import argparse

def path_modify(path):
	if not path[-1]=='/':
		path+='/'
	return path

def clear_parameters():
	parser = argparse.ArgumentParser(description='step2.integrate_ILL_svs.py')
	parser.add_argument('input_path', help='directory of input vcf')
	parser.add_argument('reference', help='reference genome used for bam alignment')
	parser.add_argument('blacklist', help='blacklist regions to be excluded in the integration')
	args = parser.parse_args()
	global ref,ppre,blacklist,size_cff,sample_names,Process_Steps,Process_INS
	ref = args.reference
	ppre = path_modify(args.input_path)
	blacklist = args.blacklist
	size_cff = 50
	sample_names = ['HG00512_HG00513_HG00514_HG00731_HG00732_HG00733_NA19238_NA19239_NA19240']
	Process_Steps = ['Illumina_readin','PacBio_readin','Origianl_Compare','Cluster_ILL','Cleanup_Merged_SVs','QCs','merge_INS_with_others']
	Process_INS = ['ILL_PB_Comparison']
	sv_list = ['DEL','DUP','INV','DEL;DUP','DEL;INV','DUP;INV','DEL;DUP;INV']

clear_parameters()

###################Step0a, readin SVs from each Illumina callset and integrate them into a united file###################
###################Step0a, output file='STEP0_ILL_Calls.bed'#############################################
#GATK???????
if Process_Steps[0]=='Illumina_readin':
	def algorithm_num_count(list):
		#eg of list=[['chr22', 50434319, 50434865, 'DEL', '0/1', 5], ['chr22', 50434319, 50434896, 'DEL', '0/1', 10]]
		out=[]
		for x in list:
			if not x[-1] in out:
				out.append(x[-1])
		return len(out)
	def algorithm_type_count(list):
		out=[]
		for x in list:
			if not x[-1] in out:
				out.append(x[-1])
		return out
	def bashir_pacbio_info_unify(bashir_pacbio_info):
		out={}
		for k1 in bashir_pacbio_info.keys():
			if not k1 in out.keys():	out[k1]={}
			for k2 in bashir_pacbio_info[k1].keys():
				if not k2 in out[k1].keys():	out[k1][k2]={}
				for k3 in bashir_pacbio_info[k1][k2]:
					if not k3[0] in out[k1][k2].keys():	out[k1][k2][k3[0]]={}
					if not k3[1] in out[k1][k2][k3[0]].keys():	out[k1][k2][k3[0]][k3[1]]={}
					if not k3[2] in out[k1][k2][k3[0]][k3[1]].keys():	out[k1][k2][k3[0]][k3[1]][k3[2]]=[]
					if not k3 in out[k1][k2][k3[0]][k3[1]][k3[2]]:	out[k1][k2][k3[0]][k3[1]][k3[2]].append(k3)
		out2={}
		for k1 in out.keys():
			out2[k1]={}
			for k2 in out[k1].keys():
				out2[k1][k2]=[]
				for k3 in sorted(out[k1][k2].keys()):
					for k4 in sorted(out[k1][k2][k3].keys()):
						for k5 in sorted(out[k1][k2][k3][k4].keys()):
							for k6 in sorted(out[k1][k2][k3][k4][k5]):
								if not k6 in out2[k1][k2]: out2[k1][k2].append(k6)
		return out2
	def bashir_pacbio_info_modify(bashir_pacbio_path):
		bashir_pacbio_filenames=vcf_name_readin(bashir_pacbio_path,'_',0)
		bashir_pacbio_info={}
		for k1 in bashir_pacbio_filenames.keys():
			bashir_pacbio_info[k1]={}
			for k2 in bashir_pacbio_filenames[k1]:
				test=vcf_readin_bashir(k2)
				for k3 in test.keys():
					if not k3 in bashir_pacbio_info[k1].keys():
						bashir_pacbio_info[k1][k3]=[]
					bashir_pacbio_info[k1][k3]+=test[k3]
		bashir_pacbio_info=bashir_pacbio_info_unify(bashir_pacbio_info)
		out={}
		for k1 in bashir_pacbio_info.keys():
			for k2 in bashir_pacbio_info[k1].keys():
				if not k2 in out.keys():
					out[k2]={}
				if not k1 in out[k2].keys():
					out[k2][k1]=bashir_pacbio_info[k1][k2]
		return out
	def eichler_pacbio_info_modify(eichler_pacbio_path):
		eichler_pacbio_filenames=vcf_name_readin(eichler_pacbio_path,'.',3)
		eichler_pacbio_info={}
		for k1 in eichler_pacbio_filenames.keys():
			eichler_pacbio_info[k1]=vcf_readin_eichler(eichler_pacbio_filenames[k1][0])
		out={}
		for k1 in eichler_pacbio_info.keys():
			for k2 in eichler_pacbio_info[k1].keys():
				if not k2 in out.keys():
					out[k2]={}
				if not k1 in out[k2].keys():
					out[k2][k1]=eichler_pacbio_info[k1][k2]
		return out
	def chr_start_end_extract(pin):
		out=[pin[0],int(pin[1])]
		sv_type=svtype_extract(pin)
		for x in pin[7].split(';'):
			if x.split('=')[0]=='END':
				out.append(int(x.split('=')[1]))
		if len(out)==3 and out[2]-out[1]<2:
			if sv_type in ['INS','ins','insertion']:
				sv_len=sv_len_extract(pin)
				out[2]=out[1]+sv_len
		elif len(out)==2:
			sv_len=sv_len_extract(pin)
			out.append(out[1]+sv_len)
		return [sv_type,out]
	def chr_start_end_extract_melt(pin):
		#eg of bp info inpin:ADJLEFT=0;ADJRIGHT=2995287
		start=int(pin[1])
		end='0'
		SVLEN=0
		for i in pin[7].split(';'):
			if i.split('=')[0]=='SVLEN': end=start+int(i.split('=')[1])
		if end=='0':
			for i in pin[7].split(';'):
				if i.split('=')[0]=='ADJRIGHT':			end=int(i.split('=')[1])
		return [pin[0],start,end]
	def chromos_readin(ref):
		fin=open(ref+'.fai')
		chromos=[]
		for line in fin:
				pin=line.strip().split()
				chromos.append(pin[0])
		fin.close()
		return chromos
	def chromos_readin_canonical(ref):
		fin=open(ref+'.fai')
		chromos=[]
		for line in fin:
				pin=line.strip().split()
				if not '-' in pin[0] and not '_' in pin[0]:
						chromos.append(pin[0])
		fin.close()
		return chromos
	def copynumber_extract_2(info_list):
		rec_pos=-1
		for x in info_list[0].split(':'):
			rec_pos+=1
			if x=='CN': break
		cn=int(info_list[1].split(':')[rec_pos])
		return cn
	def compare_between_bashir_and_eichler(eichler_pacbio_info,bashir_pacbio_info,match_sv_type):
		pacbio_vs_pacbio_stat={}
		for k1 in eichler_pacbio_info.keys():
			if k1 in bashir_pacbio_info.keys():
				pacbio_vs_pacbio_stat[k1]={}
				for k2 in eichler_pacbio_info[k1].keys():
					if k2 in bashir_pacbio_info[k1].keys():
						pacbio_vs_pacbio_stat[k1][k2]=sv_list_quick_compare(eichler_pacbio_info[k1][k2],bashir_pacbio_info[k1][k2],match_sv_type,overlap_rc)
		pacbio_vs_pacbio_stat_comparison_stat={}
		for k1 in pacbio_vs_pacbio_stat.keys():
			for k2 in pacbio_vs_pacbio_stat[k1].keys():
				if not k2 in pacbio_vs_pacbio_stat_comparison_stat.keys():
					pacbio_vs_pacbio_stat_comparison_stat[k2]=pacbio_vs_pacbio_stat[k1][k2]
				else:
					pacbio_vs_pacbio_stat_comparison_stat[k2][0]+=pacbio_vs_pacbio_stat[k1][k2][0]
					pacbio_vs_pacbio_stat_comparison_stat[k2][1]+=pacbio_vs_pacbio_stat[k1][k2][1]
					pacbio_vs_pacbio_stat_comparison_stat[k2][2]+=pacbio_vs_pacbio_stat[k1][k2][2]
		fo=open(illumina_path+'pacbio.vs.pacbio.stat.match_sv_type.'+match_sv_type,'w')
		for k1 in sorted(pacbio_vs_pacbio_stat_comparison_stat.keys()):
			print('\t'.join([str(i) for i in [k1]+pacbio_vs_pacbio_stat_comparison_stat[k1]]),file=fo)
		fo.close()
	def genotype_extract(pin):
		out=[0,0]
		rec_pos=-1
		for x in pin[8].split(':'):
			rec_pos+=1
			if x=='GT': break
		geno=pin[9].split(':')[rec_pos]
		if geno=='./.':	#take het for ./.
			out=[0,1]
		elif geno=='.':
			out=[0,1]
		else:
			if '/' in geno:
				out=[int(i) for i in geno.split('/')]
			elif '|' in geno:
				out=[int(i) for i in geno.split('|')]
		return out
	def genotype_extract_2(info_list):
		#eg of info_list=[ 'GT:CN:CNF:CNL:CNP:CNQ:GSPC', '.:3:3.0370:-1000.00,-1000.00,-1000.00,0.00,-186.39:-1000.00,-1000.00,-1000.00,0.00,-186.69:99.0:0']
		out='/'
		rec_pos=-1
		for x in info_list[0].split(':'):
			rec_pos+=1
			if x=='GT': break
		geno=info_list[1].split(':')[rec_pos]
		return geno
	def genotype_extract_3(info):
		#eg of info=pin[8:]=['GT:START:END:SVLEN:SPNUM', '0;NA;NA;NA;NA', '1;10530892;10530992;16;3', '1;10530802;10530902;65;2', '0;NA;NA;NA;NA', '0;NA;NA;NA;NA', '0;NA;NA;NA;NA', '0;NA;NA;NA;NA', '0;NA;NA;NA;NA', '1;10530829;10530929;35;3']
		geno_pos=info[0].split(':').index('GT')
		geno_info=[i.split(';')[geno_pos] for i in info[1:]]
		out=[]
		for x in geno_info:
			if x=='0': 
				out.append('1/1')
			elif x=='1': 
				out.append('0/1')
			elif x=='2': 
				out.append('0/0')
		return out
	def genotype_extract_3(pin):
		gt_pos=pin[8].split(':').index('GT')
		out=[i.split(':')[gt_pos] for i in pin[9:]]
		return out
	def genotype_extract_Manta(pin):
		geno_info_all=pin[9:]
		GT_pos=pin[8].split(':').index('GT')
		FT_pos=pin[8].split(':').index('FT')
		out=['0/0' for i in geno_info_all]
		rec=-1
		for i in geno_info_all:
			rec+=1
			if len(i.split(':'))==1 or i=='./.':	
				continue
			elif len(i.split(':'))>1 and i.split(':')[FT_pos]=='PASS':
				out[rec]=i.split(':')[GT_pos]
		for i in range(len(out)):
			if out[i] in ['./.','.']:	
				out[i]='0/0' 
		return out
	def geno_to_sv(geno):
		#eg of geno: 0/0,0/1,1/1
		#eg of output: if 0/0: 'no';
		temp=[0,0]
		if '/' in geno:
			temp=geno.split('/')
		elif '|' in geno:
			temp=geno.split('|')
		if geno=='.':
			temp=['.','.']
		if temp==['.','.']:
			temp=[0,1]	#for genotype = .  or genotype = ./., we report them as 0/1 for now. might change later ????????
		else:
			temp=[int(i) for i in temp]
		if sum(temp)==0:
			return 'no'
		else:
			return 'yes'
	def geno_to_sv_2(geno):
		#eg of geno: 0,1,2
		if geno in ['.',2]: 
			return 'no'
		elif geno in [0,1]: 
			return 'yes'
		else: print(geno)
	def geno_to_sv_3(geno):
		if geno_to_sv(geno)=='no':
			if geno.isdigit():
				if int(geno)>2:	
					return 'yes'
				else:	
					return 'no'
		else:
			return 'yes'
	def illumina_info_hash_to_list(illumina_info_hash,k1,ka):
		#eg of k1='chr22'; eg of ka='NA19238'
		out=[]
		for k2 in illumina_info_hash.keys():
			for k3 in illumina_info_hash[k2].keys():
				if k3==ka:
					num_index_appdix=algorithm_num_index[k2]
					for k4 in illumina_info_hash[k2][k3]:
						out.append(k4+[num_index_appdix])
		out_hash=list_reorder(out)
		return out_hash
	def illumina_info_list_write(illumina_info_list,file_out):
		if not os.path.isfile(file_out):
			fo=open(file_out,'w')
		else:
			fo=open(file_out,'a')
		for x in illumina_info_list:
			x[-1]=num_algorithm_index[x[-1]]
			print('\t'.join([str(i) for i in x]),file=fo)
		fo.close()
	def illumina_info_integrate(illumina_list):
		#integrate calls as long as they overlap with each other
		#eg of illumina_list=illumina_info_list
		start=0
		end=0
		out=[]
		for x in illumina_list:
			if  x[2]-x[1]>10**6: continue
			if out==[]:
				out.append([x])
				start=x[1]
				end=x[2]
			else:
				if x[1]>end:
					out.append([x])
					start=x[1]
					end=x[2]
				else:
					out[-1].append(x)
					start=min([start,x[1]])
					end=max([end,x[2]])
		return out
	def illumina_info_integrated_modify(illumina_list):
		#eg of illumina_list=illumina_left_list
		#integrate calls based no 50 % rc
		out1=illumina_info_integrate(illumina_list)
		out2=[]
		for k1 in out1:
			reciprocal_result=reciprocal_overlap_matrix(k1)
			for k2 in reciprocal_result:
				temp2=[k1[x] for x in k2]
				out2.append(temp2)
		return out2
	def illumina_info_hash_modify(illumina_info_hash):
		out={}
		for k1 in illumina_info_hash.keys():
			out[k1]={}
			for k2 in illumina_info_hash[k1].keys():
				out[k1][k2]=[]
				for k3 in illumina_info_hash[k1][k2]:
					k3[3]=k3[3].split(':')[0].upper()
					if k3[3] in ['DEL','DUP','INV']:
						out[k1][k2].append(k3)
		return out
	def illumina_reported_readin(chromo,ppre):
		file_hash={}
		illumina_info_hash={}
		size_cff=50
		file_hash['GenomeSTRiP_CNV']=ppre+'20161001_GenomeSTRiP_CNV/HGSVC_trios.illumina.hg38.GenomeSTRiP_cnv.20161001.genotypes.vcf'
		illumina_info_hash['GenomeSTRiP_CNV']=vcf_readin_GenomeSTRiP_cnv(file_hash['GenomeSTRiP_CNV'],chromo,size_cff)
		file_hash['dCGH_filt']=ppre+'20160930_dCGH_filt/ALL.dCGH_filt.illumina_high_coverage.20160930.cnv.genotypes.vcf'
		illumina_info_hash['dCGH_filt']=vcf_readin_dCGH_filt(file_hash['dCGH_filt'],chromo,size_cff)
		file_hash['Pindel']=ppre+'20170110_ALL_pindel_illumina_sites_and_genotypes/Pindel.trios.20170110.genotypes.vcf'
		illumina_info_hash['Pindel']=vcf_readin_Pindel(file_hash['Pindel'],chromo,size_cff)
		file_hash['VH']=ppre+'20161002_VH/ALL.VH.illumina_high_coverage.DEL.sorted.vcf'
		illumina_info_hash['VH']=vcf_readin_VH(file_hash['VH'],chromo,size_cff)
		file_hash['svelter']=ppre+'20170109_svelter_update/SVelter_UMich.alt_bwamem_GRCh38DH.CHS.high_coverage.vcf'
		illumina_info_hash['SVelter_UMich']=vcf_readin_SVelter_UMich(file_hash['svelter'],chromo,size_cff)
		file_hash['cloudSV']=ppre+'20161003_cloudSV/cloudSV_GRCh38_all_pb_ill.vcf'
		illumina_info_hash['cloudSV']=vcf_readin_cloudSV(file_hash['cloudSV'],chromo,size_cff)
		file_hash['retroCNV']=ppre+'20161003_retroCNV/20161003_retroCNV.vcf'
		illumina_info_hash['retroCNV']=vcf_readin_retroCNV(file_hash['retroCNV'],chromo,size_cff)
		file_hash['wham']=ppre+'20160930_wham/ALL.wham.20160930.genotypes.vcf'
		#file_hash['wham']=ppre+'20170125_wham_lumpy_revised_v2/ALL.wham.20170110.genotypes.vcf'
		illumina_info_hash['wham']=vcf_readin_wham(file_hash['wham'],chromo,size_cff)
		file_hash['lumpy']=ppre+'20160930_lumpy/ALL.lumpy.20160930.genotypes.vcf'
		#file_hash['lumpy']=ppre+'20170125_wham_lumpy_revised_v2/ALL.lumpy.20170110.genotypes.vcf'
		illumina_info_hash['lumpy']=vcf_readin_lumpy(file_hash['lumpy'],chromo,size_cff)
		file_hash['melt']='/scratch/remills_flux/xuefzhao/SV_discovery_index/download/different_callers/illumina/20160930_MELT/ALL.melt.MEI.illumina_high_coverage.20160930.genotypes.vcf'
		illumina_info_hash['melt']=melt_list_modfiy(vcf_readin_melt(file_hash['melt'],chromo,size_cff))
		delly_illumina_illumina_path='/scratch/remills_flux/xuefzhao/SV_discovery_index/download/different_callers/illumina/20161002_delly_illumina/'
		file_hash['delly_illumina']=[delly_illumina_illumina_path+i for i in ['hgsvc.delly.complex.GRCh38.20160931.high_coverage.vcf','hgsvc.delly.svs.GRCh38.20160931.high_coverage.vcf']]
		illumina_info_hash['delly_illumina']=vcf_readin_delly(file_hash['delly_illumina'],chromo,size_cff)
		#file_hash['sebat_TruSeqSLR_Moleculo']=ppre+'20161114_moleculo_pacbio_ucsd/ALL.wgs.UCSD_Moleculo.20161026.sv.Illumina_Moleculo.sites.vcf'
		#illumina_info_hash['sebat_TruSeqSLR_Moleculo']=vcf_readin_sebat_TruSeqSLR_Moleculo(file_hash['sebat_TruSeqSLR_Moleculo'],chromo,size_cff)
		UCSD_forest_SV_path='/scratch/remills_flux/xuefzhao/SV_discovery_index/download/different_callers/illumina/20161108_UCSD_ForestSV/'
		file_hash['UCSD_forestSV']=[UCSD_forest_SV_path+i for i in ['ALL.wgs.UCSD_ForestSV-gtCNV.20163110.sv.Illumina_high-coverage_PCR-free.genotypes.vcf']]
		illumina_info_hash['UCSD_forestSV']=vcf_readin_ForestSV(file_hash['UCSD_forestSV'],chromo,size_cff)
		UCSD_Manta_path='/scratch/remills_flux/xuefzhao/SV_discovery_index/download/different_callers/illumina/20161108_UCSD_Manta/'
		file_hash['UCSD_Manta']=[UCSD_Manta_path+'ALL.wgs.UCSD_Manta.20162710.sv.Illumina_high-coverage_PCR-free.genotypes.vcf']
		illumina_info_hash['UCSD_Manta']=vcf_readin_Manta(file_hash['UCSD_Manta'],chromo,size_cff)
		liWGS_illumina_path='/scratch/remills_flux/xuefzhao/SV_discovery_index/download/different_callers/illumina/20160501_liWGS_SV/'
		liWGS_illumina_filenames=liWGS_filenames_readin(liWGS_illumina_path)
		illumina_info_hash['liWGS_illumina']=liWGS_bed_readin(liWGS_illumina_filenames,chromo,size_cff)
		kchen_novobreak_illumina_path='/scratch/remills_flux/xuefzhao/SV_discovery_index/download/different_callers/illumina/20160930_kchen_novobreak/'
		file_hash['kchen_novobreak']=vcf_name_readin(kchen_novobreak_illumina_path,'_',1)
		illumina_info_hash['kchen_novobreak']=vcf_readin_kchen_novobreak(file_hash['kchen_novobreak'],chromo,size_cff)
		HySA_del_indel_path=ppre+'20170105_HySA_9samples_Pacbio_Illumina_INDEL/'
		file_hash['kchen_HySA_indel']=[HySA_del_indel_path+i for i in ['20170105_ALL.wgs.KchenLab-HySA.10122016.DEL.Illumina_Pacbio.deep_coverage.sites.vcf','20170105_ALL.wgs.KchenLab-HySA.10122016.INS.Illumina_Pacbio.deep_coverage.sites.vcf']]
		illumina_info_hash['kchen_HySA_indel']=vcf_readin_kchen_HySA(file_hash['kchen_HySA_indel'],chromo,size_cff)
		return illumina_info_hash
	def illumina_left_cluster_double_check(illumina_left_cluster):
		info_all=[]
		for x in illumina_left_cluster:
			temp=[]
			for y in x:
				temp+=y[1:3]
			info_all.append([min(temp),max(temp)])
		for x in range(len(info_all)-1):
			if info_all[x+1][0]<info_all[x][1]:
				print([info_all[x],info_all[x+1]])
	def liWGS_filenames_readin(liWGS_illumina_path):
		out=[]
		for k1 in os.listdir(liWGS_illumina_path):
			if 'hg38' in k1 and not 'lowquality' in k1:
				out.append(liWGS_illumina_path+k1)
		return out
	def liWGS_bed_readin(liWGS_illumina_filenames,chromo,size_cff):
		out={}
		for x in liWGS_illumina_filenames:
			if 'deletion' in x:
				fin=open(x)
				pin=fin.readline().strip().split()
				for line in fin:
					pin=line.strip().split()
					if not 'chr'+pin[0]==chromo: continue
					samples=[i.split('_')[0].replace('GM','NA') for i in pin[6].replace('[','').replace(']','').split(',')]
					pos_info=['chr'+pin[0],int(pin[1]),int(pin[2])]+['del','0/1']
					if not pos_info[2]-pos_info[1]>size_cff: continue
					for ka in samples:
						if not ka in out.keys():	out[ka]={}
						if not pos_info[0] in out[ka].keys():	out[ka][pos_info[0]]={}
						if not pos_info[1] in out[ka][pos_info[0]].keys():	out[ka][pos_info[0]][pos_info[1]]={}
						if not pos_info[2] in out[ka][pos_info[0]][pos_info[1]].keys():	out[ka][pos_info[0]][pos_info[1]][pos_info[2]]=[]
						if not pos_info in out[ka][pos_info[0]][pos_info[1]][pos_info[2]]:	out[ka][pos_info[0]][pos_info[1]][pos_info[2]].append(pos_info)
				fin.close()
			elif 'duplication' in x:
				fin=open(x)
				pin=fin.readline().strip().split()
				for line in fin:
					pin=line.strip().split()
					if not 'chr'+pin[0]==chromo: continue
					samples=[i.split('_')[0].replace('GM','NA') for i in pin[6].replace('[','').replace(']','').split(',')]
					pos_info=['chr'+pin[0],int(pin[1]),int(pin[2])]+['dup','0/1']
					if not pos_info[2]-pos_info[1]>size_cff: continue
					for ka in samples:
						if not ka in out.keys():	out[ka]={}
						if not pos_info[0] in out[ka].keys():	out[ka][pos_info[0]]={}
						if not pos_info[1] in out[ka][pos_info[0]].keys():	out[ka][pos_info[0]][pos_info[1]]={}
						if not pos_info[2] in out[ka][pos_info[0]][pos_info[1]].keys():	out[ka][pos_info[0]][pos_info[1]][pos_info[2]]=[]
						if not pos_info in out[ka][pos_info[0]][pos_info[1]][pos_info[2]]:	out[ka][pos_info[0]][pos_info[1]][pos_info[2]].append(pos_info)
				fin.close()
			elif 'inversion' in x:
				fin=open(x)
				pin=fin.readline().strip().split()
				for line in fin:
					pin=line.strip().split()
					if not 'chr'+pin[0]==chromo: continue
					samples=[i.split('_')[0].replace('GM','NA') for i in pin[6].replace('[','').replace(']','').split(',')]
					pos_info=['chr'+pin[0],int(pin[1]),int(pin[2])]+['inv','0/1']
					if not pos_info[2]-pos_info[1]>size_cff: continue
					for ka in samples:
						if not ka in out.keys():	out[ka]={}
						if not pos_info[0] in out[ka].keys():	out[ka][pos_info[0]]={}
						if not pos_info[1] in out[ka][pos_info[0]].keys():	out[ka][pos_info[0]][pos_info[1]]={}
						if not pos_info[2] in out[ka][pos_info[0]][pos_info[1]].keys():	out[ka][pos_info[0]][pos_info[1]][pos_info[2]]=[]
						if not pos_info in out[ka][pos_info[0]][pos_info[1]][pos_info[2]]:	out[ka][pos_info[0]][pos_info[1]][pos_info[2]].append(pos_info)
				fin.close()
			else:			print(x)
		return sv_hash_order_2(out)
	def left_cluster_to_stat(illumina_left_cluster):
		out={}
		for x in illumina_left_cluster:
			x_key=[]
			for y in x:
				if not y[-1] in x_key:
					x_key.append(y[-1])
			if not len(x_key) in out.keys():
				out[len(x_key)]=[]
			out[len(x_key)]+=x_key
		return out
	def list_reorder(sv_list):
		out_1={}
		for k1 in sv_list:
			if not k1[0] in out_1.keys():
				out_1[k1[0]]={}
			if not k1[1] in out_1[k1[0]].keys():
				out_1[k1[0]][k1[1]]={}
			if not k1[2] in out_1[k1[0]][k1[1]].keys():
				out_1[k1[0]][k1[1]][k1[2]]=[]
			if not k1 in out_1[k1[0]][k1[1]][k1[2]]:
				out_1[k1[0]][k1[1]][k1[2]].append(k1)
		out_2=[]
		for k1 in out_1.keys():
			for k2 in sorted(out_1[k1].keys()):
				for k3 in sorted(out_1[k1][k2].keys()):
					for k4 in sorted(out_1[k1][k2][k3]):
						out_2.append(k4)
		return out_2
	def melt_list_modfiy(melt_list):
		out={}
		for k1 in melt_list.keys():
			out[k1]=[]
			for k2 in melt_list[k1]:
				k2[-2]='DEL'
				out[k1].append(k2)
		return out
	def num_algorithm_index_build(illumina_info_hash):
		global num_algorithm_index
		num_algorithm_index={}
		global algorithm_num_index
		algorithm_num_index={}
		out=[]
		rec=0
		for kx in sorted(illumina_info_hash.keys()):
			rec+=1
			num_algorithm_index[rec]=kx
			algorithm_num_index[kx]=rec
	def path_mkdir(path):
		if not os.path.isdir(path):
				os.system(r'''mkdir %s'''%(path))
	def reciprocal_overlap_matrix(k1):
		out=[]
		for ka in range(len(k1)):
			for kb in range(ka,len(k1)):
				if ka==kb: continue
				else:
					test=round(overlap_calcu(k1[ka],k1[kb],'TRUE'),2)
					if test>0.5:	out.append([ka,kb])
		if out==[]: return [range(len(k1))]
		out_hash={}
		for ka in out:
			if not ka[0] in out_hash.keys():
				out_hash[ka[0]]=[]
			out_hash[ka[0]].append(ka[1])
			if not ka[1] in out_hash.keys():
				out_hash[ka[1]]=[]
			out_hash[ka[1]].append(ka[0])
		out2=[]
		while True:
			if out_hash=={}: break
			temp=[out_hash.keys()[0]]
			while True:
				temp_rec=[i for i in temp]
				for x in temp:
					if x in out_hash.keys():
						temp+=out_hash[x]
						del out_hash[x]
				temp=sorted(unify_list(temp))
				if temp==temp_rec: break
			out2.append(temp)
		out5=[]
		for i in range(len(k1)):
			flag=0
			for j in out2:
				if i in j: flag+=1
			if flag==0: out5.append([i])
		return out2+out5
	def rec_retrive(file_name):
		fin=open(file_name)
		rec=0
		for line in fin:
			pin=line.strip().split()
			if int(pin[-2]) >rec: rec=int(pin[-2])
		fin.close()
		return rec
	def stat_hash_to_hash(bashir_left_stat):
		out={}
		for k1 in bashir_left_stat.keys():
			for k2 in bashir_left_stat[k1].keys():
				if not k2 in out.keys():
					out[k2]={}
				for k3 in bashir_left_stat[k1][k2].keys():
					if not k3 in out[k2].keys():
						out[k2][k3]={}
					for k4 in bashir_left_stat[k1][k2][k3]:
						if not k4 in out[k2][k3].keys():
							out[k2][k3][k4]=1
						else:
							out[k2][k3][k4]+=1
		return out
	def stat_list_to_hash(bashir_stat):
		out={}
		for k1 in bashir_stat.keys():
			for k2 in bashir_stat[k1].keys():
				if not k2 in out.keys():
					out[k2]={}
				for k3 in bashir_stat[k1][k2]:
					if not k3 in out[k2].keys():
						out[k2][k3]=1
					else:
						out[k2][k3]+=1
		return out
	def stat_list_to_hash_2(bashir_stat_2):
		out={}
		for k1 in bashir_stat_2.keys():
			for k2 in bashir_stat_2[k1].keys():
				if not k2 in out.keys():
					out[k2]={}
				for k3 in bashir_stat_2[k1][k2].keys():
					if not k3 in out[k2].keys():
						out[k2][k3]={}
					for k4 in bashir_stat_2[k1][k2][k3]:
						if not k4 in out[k2][k3].keys():
							out[k2][k3][k4]=1
						else:
							out[k2][k3][k4]+=1
		return out
	def stat_list_write(bashir_comparison_stat,file_out):
		fo=open(file_out,'w')
		for k1 in sorted(bashir_comparison_stat.keys()):
			for k2 in sorted(bashir_comparison_stat[k1].keys()):
				print('\t'.join([str(i) for i in [k1,k2,bashir_comparison_stat[k1][k2]]]),file=fo)
		fo.close()
	def stat_2_list_write(bashir_comparison_stat_2,file_out_name):
		fo=open(file_out_name,'w')
		for k1 in bashir_comparison_stat_2.keys():
			for k2 in bashir_comparison_stat_2[k1].keys():
				for k3 in bashir_comparison_stat_2[k1][k2].keys():
					algorithm_name=num_algorithm_index[k3]
					print('\t'.join([str(i) for i in [k1,k2,k3,bashir_comparison_stat_2[k1][k2][k3],algorithm_name]]),file=fo)
		fo.close()
	def sv_list_quick_compare(list1,list2,match_sv_type,overlap_rc=0.5):
		out=[len(list1),len(list2),0]
		rec_a=0 #record pacbio_list
		rec_b=0 #record illumina_list
		for k1 in list1:
			for k2 in list2:
				if k1[0]==k2[0]:
					if k2[1]>k1[2]: break
					elif k2[2]<k1[1]: continue
					else:
						overlap=overlap_calcu(k1,k2,match_sv_type)
						if overlap>overlap_rc:
							out[2]+=1
		return out
	def svtype_extract(pin):
		svtype=''
		for x in pin[7].split(';'):
			if 'SVTYPE' in x:
				svtype=x.split('=')[1]
		if svtype=='':
			svtype=pin[4].replace('<','').replace('>','')
		return svtype
	def svtype_svpos_csv_extract(pin):
		out1=[]
		out2=[]
		for x in pin[7].split(';'):
			if x.split('=')[0]=='del':	
				out1.append('del')
				out2.append([x.split('=')[1].split(':')[0]]+[int(i) for i in x.split('=')[1].split(':')[1].split('-')])
			elif x.split('=')[0] in ['dup','disdup','dup_inv']:
				out1.append('dup')
				out2.append([x.split('=')[1].split(':')[0]]+[int(i) for i in x.split('=')[1].split(':')[1].split('-')])
			elif x.split('=')[0]=='inv':
				out1.append('inv')
				out2.append([x.split('=')[1].split(':')[0]]+[int(i) for i in x.split('=')[1].split(':')[1].split('-')])
		return [out1,out2]
	def sv_len_extract(pin):
		out=0
		for x in pin[7].split(';'):
			if 'SVLEN' in x:
				out=int(x.split('=')[1])
		return out
	def SV_num_stat(readin_hash):
		out={}
		stat_hash={}
		for k1 in readin_hash.keys():
			if not k1 in stat_hash.keys():			stat_hash[k1]={}
			for k2 in readin_hash[k1]:	
				if not k2[3] in stat_hash[k1].keys():		stat_hash[k1][k2[3]]=0
				stat_hash[k1][k2[3]]+=1
				if not k2[0] in out.keys():					out[k2[0]]={}
				if not k2[1] in out[k2[0]].keys():			out[k2[0]][k2[1]]={}
				if not k2[2] in out[k2[0]][k2[1]].keys():	out[k2[0]][k2[1]][k2[2]]=[]
				if not k2[3] in out[k2[0]][k2[1]][k2[2]]:	out[k2[0]][k2[1]][k2[2]].append(k2[3])
		stat_hash['all']={}
		for k1 in out.keys():
			for k2 in out[k1].keys():
				for k3 in out[k1][k2].keys():
					for k4 in out[k1][k2][k3]:
						if not k4 in stat_hash['all'].keys():	stat_hash['all'][k4]=0
						stat_hash['all'][k4]+=1
		return stat_hash
	def sv_hash_order(sv_hash):
		out={}
		for k1 in sv_hash.keys():
			out[k1]=[]
			for k2 in sorted(sv_hash[k1].keys()):
				for k3 in sorted(sv_hash[k1][k2].keys()):
					for k4 in sorted(sv_hash[k1][k2][k3]):
						out[k1].append(k4)
		return out
	def sv_hash_order_2(sv_hash):
		out={}
		for k1 in sv_hash.keys():
			out[k1]=[]
			for k2 in sorted(sv_hash[k1].keys()):
				for k3 in sorted(sv_hash[k1][k2].keys()):
					for k4 in sorted(sv_hash[k1][k2][k3].keys()):
						for k5 in sv_hash[k1][k2][k3][k4]:
							k5[3]=k5[3].upper()
							if not k5 in out[k1]:						out[k1].append(k5)
		return out
	def sv_info_comparison(pacbio_list,illumina_list,match_sv_type,file_out_name,overlap_rc=0.5):
		#eg of pacbio_list=bashir_pacbio_info[k1][ka]
		#eg of illumina_list=illumina_info_list
		#rec_a=0 #record pacbio_list
		#rec_b=0 #record illumina_list
		pacbio_stat_list={}
		illumina_kept=[]
		rec_a=-1
		for k1 in pacbio_list:
			rec_a+=1
			flag=0
			for k2 in illumina_list:
				overlap=overlap_calcu(k1,k2,match_sv_type)
				if overlap>overlap_rc:
					if not rec_a in pacbio_stat_list.keys():
						pacbio_stat_list[rec_a]=[]
					if not k2 in pacbio_stat_list[rec_a]:
						pacbio_stat_list[rec_a].append(k2)
					if not k2 in illumina_kept:	illumina_kept.append(k2)
		out_stat=[]
		out_stat_2={}
		file_initiate(file_out_name)
		fo=open(file_out_name,'a')
		for k1 in pacbio_stat_list.keys():
			temp_info=algorithm_type_count(pacbio_stat_list[k1])
			out_stat.append(algorithm_num_count(pacbio_stat_list[k1]))
			print('\t'.join([str(i) for i in pacbio_list[k1]+[out_stat[-1]]]),file=fo)
			if len(temp_info)>0: 
				if not len(temp_info) in out_stat_2.keys():
					out_stat_2[len(temp_info)]=[]
				out_stat_2[len(temp_info)]+=temp_info
		illumina_left=[]
		for x in illumina_list:
			if not x in illumina_kept:
				illumina_left.append(x)
		return [out_stat,illumina_left,out_stat_2]
	def sv_info_comparison_del(pacbio_list,illumina_list,match_sv_type,overlap_rc=0.5):
		#eg of pacbio_list=bashir_pacbio_info[k1][ka]
		#eg of illumina_list=illumina_info_list
		rec_a=0 #record pacbio_list
		rec_b=0 #record illumina_list
		pacbio_stat_list={}
		illumina_kept=[]
		while True:
			if rec_a==len(pacbio_list): break
			if rec_b==len(illumina_list):
				break
			else:
				pacbio_unit=pacbio_list[rec_a]
				if not rec_a in pacbio_stat_list.keys():
					pacbio_stat_list[rec_a]=[]
				illumina_unit=illumina_list[rec_b]
				if illumina_unit[2]<pacbio_unit[1]: rec_b+=1
				elif illumina_unit[1]>pacbio_unit[2]: rec_a+=1
				else:
					overlap=overlap_calcu_del(pacbio_unit,illumina_unit,match_sv_type)
					rec_b+=1
					if overlap>overlap_rc:
						pacbio_stat_list[rec_a].append(illumina_unit)
						illumina_kept.append(illumina_unit)
		out_stat=[]
		out_stat_2={}
		for k1 in pacbio_stat_list.keys():
			temp_info=algorithm_type_count(pacbio_stat_list[k1])
			out_stat.append(algorithm_num_count(pacbio_stat_list[k1]))
			if len(temp_info)>0: 
				if not len(temp_info) in out_stat_2.keys():
					out_stat_2[len(temp_info)]=[]
				out_stat_2[len(temp_info)]+=temp_info
		illumina_left=[]
		for x in illumina_list:
			if not x in illumina_kept:
				illumina_left.append(x)
		return [out_stat,illumina_left,out_stat_2]
	def overlap_calcu(list1,list2,match_sv_type):
		#eg of list1=['chr22', 10717135, 10717147, 'INS', 0, 1]
		#eg of list2=['chr22', 10527787, 10700173, 'CNV', '0', 7]
		if not list1[0]==list2[0]: return 0
		if match_sv_type=='TRUE':	
			if not list1[3].upper()==list2[3].upper(): return 0
		if list2[1]>list1[2]: return 0
		if list1[1]>list2[2]: return 0
		return float(sorted(list1[1:3]+list2[1:3])[2]-sorted(list1[1:3]+list2[1:3])[1])/float(max([list1[2]-list1[1],list2[2]-list2[1]]))
	def overlap_calcu_del(list1,list2,match_sv_type):
		#eg of list1=['chr22', 10717135, 10717147, 'INS', 0, 1]
		#eg of list2=['chr22', 10527787, 10700173, 'CNV', '0', 7]
		if not list1[0]==list2[0]: return 0
		if match_sv_type=='TRUE':	
			if not list1[3]==list2[3]: return 0
		if not list1[3]=='DEL': return 0
		if list2[1]>list1[2]: return 0
		if list1[1]>list2[2]: return 0
		return float(sorted(list1[1:3]+list2[1:3])[2]-sorted(list1[1:3]+list2[1:3])[1])/float(max([list1[2]-list1[1],list2[2]-list2[1]]))
	def pacbio_info_hash_modify(test):
		out={}
		for k1 in test.keys():
			out[k1]={}
			for k2 in test[k1].keys():
				out[k1][k2]=[]
				for k3 in test[k1][k2]:
					k3[3]=k3[3][:3].upper()
					if not k3 in out[k1][k2]:
						out[k1][k2].append(k3)
		return out
	def similarity_Lists(list1,list2):
		#decide if 2 lists share the same item
		out='FALSE'
		test=sorted(list1+list2)
		for x in test:
			if test.count(x)>list1.count(x) and test.count(x)>list2.count(x):
				out='TRUE'
		return out
	def vcf_name_readin(path,sep,key_pos):
		#eg of path='/scratch/remills_flux/xuefzhao/SV_discovery_index/download/different_callers/pacbio_pre_cshl/bashir_svs/'
		#eg of sep='_'
		#eg of key_pos=0
		out={}
		for k1 in os.listdir(path):
			if k1.split('.')[-1]=='vcf':
				key_info=k1.split(sep)[key_pos]
				if not key_info in out.keys():
					out[key_info]=[]
				if not path+k1 in out[key_info]:
					out[key_info].append(path+k1)
		for k1 in out.keys():
			out[k1].sort()
		return out
	def vcf_name_list_readin(path):
		out=[]
		for k1 in os.listdir(path):
			if k1.split('.')[-1]=='vcf':
				out.append(path+k1)
		return sorted(out)
	def vcf_readin_bashir(vcf_name):
		#eg of vcf_name='/scratch/remills_flux/xuefzhao/SV_discovery_index/download/different_callers/pacbio_pre_cshl/bashir_svs/HG00512_PacBio_structural_variants_20160502.vcf'
		fin=open(vcf_name)
		info_hash={}
		for line in fin:
			pin=line.strip().split()
			if not pin[0][0]=='#':
				[sv_type,sv_pos]=chr_start_end_extract(pin)
				if len(sv_pos)==3:
					sv_main_info=sv_pos+[sv_type]+genotype_extract(pin)
					if sum(sv_main_info[-2:])>0 and sv_main_info[2]-sv_main_info[1]>size_cff:
						if not sv_main_info[0] in info_hash.keys():
							info_hash[sv_main_info[0]]={}
						if not sv_main_info[1] in info_hash[sv_main_info[0]].keys():
							info_hash[sv_main_info[0]][sv_main_info[1]]={}
						if not sv_main_info[2] in info_hash[sv_main_info[0]][sv_main_info[1]].keys():
							info_hash[sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]=[]
						if not sv_main_info in info_hash[sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]:
							info_hash[sv_main_info[0]][sv_main_info[1]][sv_main_info[2]].append(sv_main_info)
		fin.close()
		return sv_hash_order(info_hash)
	def vcf_readin_eichler(vcf_name):
		#eg of vcf_name='/scratch/remills_flux/xuefzhao/SV_discovery_index/download/different_callers/pacbio_pre_cshl/eichler_svs/HG00512_SMRT_SV_structural_variants_20160502.vcf'
		fin=open(vcf_name)
		info_hash={}
		for line in fin:
			pin=line.strip().split()
			if not pin[0][0]=='#' and pin[6]=='PASS':  
				[sv_type,sv_pos]=chr_start_end_extract(pin)
				if len(sv_pos)==3 and sv_pos[2]-sv_pos[1]>size_cff:
					if len(pin)==10:
						sv_main_info=sv_pos+[sv_type]+genotype_extract(pin)
					elif len(pin)==8:
						sv_main_info=sv_pos+[sv_type]					
					if not sv_main_info[0] in info_hash.keys():
						info_hash[sv_main_info[0]]={}
					if not sv_main_info[1] in info_hash[sv_main_info[0]].keys():
						info_hash[sv_main_info[0]][sv_main_info[1]]={}
					if not sv_main_info[2] in info_hash[sv_main_info[0]][sv_main_info[1]].keys():
						info_hash[sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]=[]
					if not sv_main_info in info_hash[sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]:
						info_hash[sv_main_info[0]][sv_main_info[1]][sv_main_info[2]].append(sv_main_info)
		fin.close()
		return sv_hash_order(info_hash)

###################Step0b, readin SVs from integrated pacbio calls###################
###################Step0b, output file='STEP0_Pacbio_Calls.bed###################
if Process_Steps[1]=='PacBio_readin':
	def merge_hash(bed_info,vcf_info):
		out={}
		for k1 in bed_info.keys():
			out[k1]=bed_info[k1]
			if k1 in vcf_info.keys():
				for k2 in vcf_info[k1].keys():
					if not k2 in out[k1].keys():	out[k1][k2]={}
					for k3 in vcf_info[k1][k2].keys():
						if not k3 in out[k1][k2].keys():	out[k1][k2][k3]=[]
						out[k1][k2][k3]+=vcf_info[k1][k2][k3]
		for k1 in vcf_info.keys():
			if not k1 in out.keys():
				out[k1]=vcf_info[k1]
		return out
	def pacbio_bed_readin(bed_file,chromo):
		out={}
		fin=open(bed_file)
		for line in fin:
			pin=line.strip().split()
			if not pin[0][0]=='#':
				if not pin[0]==chromo: continue
				if not pin[0] in out.keys():	out[pin[0]]={}
				if not int(pin[1]) in out[pin[0]].keys():	out[pin[0]][int(pin[1])]={}
				if not int(pin[2]) in out[pin[0]][int(pin[1])].keys():	out[pin[0]][int(pin[1])][int(pin[2])]=[]
				if 'inversions' in bed_file:	
					if 'HOM' in pin[-1]:out[pin[0]][int(pin[1])][int(pin[2])].append([pin[0],int(pin[1]),int(pin[2]),'INV','1/1'])
					if 'HAP' in pin[-1]:out[pin[0]][int(pin[1])][int(pin[2])].append([pin[0],int(pin[1]),int(pin[2]),'INV','0/1'])
		fin.close()
		return out
	def pacbio_vcf_readin(vcf_file,chromo):
		out={}
		fin=open(vcf_file)
		for line in fin:
			pin=line.strip().split()
			if not pin[0][0]=='#':
				if not pin[0]==chromo: continue
				[sv_type,sv_pos]=chr_start_end_extract(pin)
				genotype=pin[-1]
				if not sv_pos[0] in out.keys():	out[sv_pos[0]]={}
				if not int(sv_pos[1]) in out[sv_pos[0]].keys():	out[sv_pos[0]][int(sv_pos[1])]={}
				if not int(sv_pos[2]) in out[sv_pos[0]][int(sv_pos[1])].keys():	out[sv_pos[0]][int(sv_pos[1])][int(sv_pos[2])]=[]
				out[sv_pos[0]][int(sv_pos[1])][int(sv_pos[2])].append(sv_pos+[sv_type,genotype])
		return out
	def pacbio_merged_readin(pacbio_path,chromo):
		bed_hash={}
		vcf_hash={}
		for k1 in os.listdir(pacbio_path):
			if k1.split('.')[-1]=='bed' and 'inver' in k1:
				bed_hash[k1.split('.')[0]]=pacbio_path+k1
			elif k1.split('.')[-1]=='vcf':
				vcf_hash[k1.split('.')[0]]=pacbio_path+k1
		bed_info={}
		vcf_info={}
		for k1 in bed_hash.keys():
			bed_info[k1]=pacbio_bed_readin(bed_hash[k1],chromo)
		for k1 in vcf_hash.keys():
			vcf_info[k1]=pacbio_vcf_readin(vcf_hash[k1],chromo)
		out={}
		for k1 in bed_info.keys():
			out[k1]=merge_hash(bed_info[k1],vcf_info[k1])
		return sv_hash_order_2(out)
	def main_0b(file_out,pacbio_path):
		fo=open(file_out,'w')
		print('\t'.join(['CHR','START','END','SVTYPE','GT','SAMPLE']), file=fo)
		for tmp in [i for i in range(23)[1:]]+['X']:
			chromo='chr'+str(tmp)
			pacbio_merged_info=pacbio_merged_readin(pacbio_path,chromo)
			for k1 in pacbio_merged_info.keys():
				for k2 in pacbio_merged_info[k1]:
					print ('\t'.join([str(i) for i in k2+[k1]]), file=fo)
		fo.close()


#Compare insertion point of ILL callers to PB calls:
#!python
if Process_INS[0]=='ILL_PB_Comparison':
	def PB_INS_readin(PB_file):
		out={}
		fin=open(PB_file)
		for line in fin:
			pin=line.strip().split()
			if not 'CHR' in pin:
				if pin[3]=='INS':
					if not pin[-1] in out.keys():	out[pin[-1]]={}
					if not pin[0] in out[pin[-1]].keys():	out[pin[-1]][pin[0]]={}
					if not int(pin[1]) in out[pin[-1]][pin[0]].keys():	out[pin[-1]][pin[0]][int(pin[1])]={}
					if not int(pin[2]) in out[pin[-1]][pin[0]][int(pin[1])].keys():	out[pin[-1]][pin[0]][int(pin[1])][int(pin[2])]=[]
					if not pin[4:] in out[pin[-1]][pin[0]][int(pin[1])][int(pin[2])]:	out[pin[-1]][pin[0]][int(pin[1])][int(pin[2])].append(pin[4:])
		fin.close()
		return out
	def ILL_INS_readin(ILL_file):
		out={}
		fin=open(ILL_file)
		for line in fin:
			pin=line.strip().split()
			if not 'CHR' in pin:
				if pin[3] in ['ALU','LINE1','SVA','HERVK','INS','DUP:RETROTRANSPOSITION', 'INS:ME:ALU','INS:ME:SVA','INS:ME:LINE1','INS:MEI:Alu','INS:MEI:L1','DUP:LINE:L1:L1MED','DUP:LINE:L1:L1MCA']:
					if not pin[-1] in out.keys():	out[pin[-1]]={}
					if not pin[0] in out[pin[-1]].keys():	out[pin[-1]][pin[0]]={}
					if not int(pin[1]) in out[pin[-1]][pin[0]].keys():	out[pin[-1]][pin[0]][int(pin[1])]={}
					if not int(pin[2]) in out[pin[-1]][pin[0]][int(pin[1])].keys():	out[pin[-1]][pin[0]][int(pin[1])][int(pin[2])]=[]
					if not pin[4:] in out[pin[-1]][pin[0]][int(pin[1])][int(pin[2])]:	out[pin[-1]][pin[0]][int(pin[1])][int(pin[2])].append(pin[4:])
		fin.close()
		return out
	def assign_ILL_to_PB(PB_INSs,ILL_INSs):
		out={}
		for k1 in ILL_INSs:
			temp=[]
			for k2 in PB_INSs:
				if temp==[]:
					temp.append(abs(k1-k2))
					temp.append(k2)
				else:
					if abs(k1-k2)<temp[0]:
						temp[0]=abs(k1-k2)
						temp[1]=k2
			if not temp[1] in out.keys():			out[temp[1]]=[]
			out[temp[1]].append(k1)
		return out
	def PB_ILL_INS_POS_Compare(PB_INS_list,ILL_INS_list):
		#assign ILL ins point to its cloest PB ins point
		out={}
		for k1 in PB_INS_list.keys():
			if k1 in ILL_INS_list.keys():
				out[k1]={}
				for k2 in PB_INS_list[k1].keys():
					if k2 in ILL_INS_list[k1].keys():
						out[k1][k2]={}
						PB_INSs=sorted(PB_INS_list[k1][k2].keys())
						ILL_INSs=sorted(ILL_INS_list[k1][k2].keys())
						ILL_PB_hash=assign_ILL_to_PB(PB_INSs,ILL_INSs)
						for ka in ILL_PB_hash.keys():
							out[k1][k2][ka]=[[],[]]	#PB_INS_Length,BPlist, ILL_INS_length
							for kb in PB_INS_list[k1][k2][ka].keys():	out[k1][k2][ka][0].append(kb-ka)
							for kb in ILL_PB_hash[ka]:
								for kc in ILL_INS_list[k1][k2][kb]:
									for kd in ILL_INS_list[k1][k2][kb][kc]:
										out[k1][k2][ka][1].append([kb,kc]+kd)
		return out
	def PB_ILL_INS_COMP_to_comp_list(PB_ILL_INS_COMP):
		bp_dis={}
		len_dis={}
		for k1 in PB_ILL_INS_COMP.keys():
			for k2 in PB_ILL_INS_COMP[k1].keys():
				for k3 in PB_ILL_INS_COMP[k1][k2].keys():
					for k4 in PB_ILL_INS_COMP[k1][k2][k3][1]:
						if not k4[3] in bp_dis.keys():	
							bp_dis[k4[3]]=[]
							len_dis[k4[3]]=[]
						bp_dis[k4[3]].append(k4[0]-k3)
						len_dis[k4[3]].append(k4[1]-k4[0]-PB_ILL_INS_COMP[k1][k2][k3][0][0])
		return [bp_dis,len_dis]
	def write_bp_len_dis(bp_dis,len_dis):
		for k1 in bp_dis.keys():
			fo=open(ppre+'Ins_dis.'+k1+'.rec','w')
			for k2 in range(len(bp_dis[k1])):
				print(' '.join([str(i) for i in [bp_dis[k1][k2],len_dis[k1][k2]]]),file=fo)
			fo.close()
	def split_ILL_file_to_INS_and_NONINS(ILL_file):
		fin=open(ILL_file)
		fo=open('.'.join(ILL_file.split('.')[:-1]+['INS',ILL_file.split('.')[-1]]),'w')
		fo2=open('.'.join(ILL_file.split('.')[:-1]+['Non_INS',ILL_file.split('.')[-1]]),'w')
		for line in fin:
			pin=line.strip().split()
			pin[3]=pin[3].upper()
			if not 'CHR' in pin:
				if pin[3] in ['INS','DUP:RETROTRANSPOSITION', 'INS:ME:ALU','INS:ME:SVA','INS:ME:LINE1','INS:MEI:Alu','INS:MEI:L1','DUP:LINE:L1:L1MED','DUP:LINE:L1:L1MCA','ALU','LINE1','SVA','HERVK']:
					if not pin[5] in ['wham','Wham']:
						print('\t'.join(pin),file=fo)
				else:
					print('\t'.join(pin),file=fo2)
		fin.close()
		fo.close()
		fo2.close()
	def split_ILL_Callers(ILL_file):
		fin=open(ILL_file.replace('_all.bed','_all.Non_INS.bed'))
		fo=open(ILL_file.replace('_all.bed','.Non_INS.bed'),'w')
		for line in fin:
			pin=line.strip().split()
			if not 'Moleculo' in pin and not 'tenX' in pin and not 'HySA' in pin and not 'Moleculo_PB' in pin:
				print('\t'.join(pin),file=fo)
		fin.close()
		fo.close()
	def split_ILL_INS_Callers(ILL_file):
		fin=open(ILL_file.replace('_all.bed','_all.INS.bed'))
		fo=open(ILL_file.replace('_all.bed','.INS.bed'),'w')
		for line in fin:
			pin=line.strip().split()
			if not 'Moleculo' in pin and not 'tenX' in pin and not 'HySA' in pin and not 'Moleculo_PB' in pin:
				print('\t'.join(pin),file=fo)
		fin.close()
		fo.close()
	def main_0c():
		PB_file=ppre+'STEP0_Pacbio_Calls.bed'
		ILL_file=ppre+'STEP0_ILL_Calls_all.bed'
		#os.system(r'''cp %s %s'''%(ppre+'STEP0_ILL_Calls.bed',ppre+'STEP0_ILL_Calls_all.bed'))
		PB_INS_list=PB_INS_readin(PB_file)
		ILL_INS_list=ILL_INS_readin(ILL_file)
		PB_ILL_INS_COMP=PB_ILL_INS_POS_Compare(PB_INS_list,ILL_INS_list)
		[bp_dis,len_dis]=PB_ILL_INS_COMP_to_comp_list(PB_ILL_INS_COMP)
		write_bp_len_dis(bp_dis,len_dis)
		split_ILL_file_to_INS_and_NONINS(ILL_file)
		split_ILL_file_to_INS_and_NONINS(PB_file)
		split_ILL_Callers(ILL_file)
		split_ILL_INS_Callers(ILL_file)

main_0c()

###################Step1. compare each illumina caller with emrged pacbio, first by 50% RC###################
if Process_Steps[2]=='Origianl_Compare':
	def bp_off_decide_both_end_separate_based_on_RC(illumina_list,pacbio_list,match_sv_type='FALSE',RC=0.5):
		out=[[],[]]
		rec=[]
		for k1 in illumina_list:
			for k2 in pacbio_list:
				RC_overlap=overlap_calcu(k1,k2,match_sv_type)
				if RC_overlap>RC:
					rec.append(RC_overlap)
					out[0].append(k2[1]-k1[1])
					out[1].append(k2[2]-k1[2])
		return [out,rec]
	def illumina_caller_stat(illumina_info_hash):
		out={}
		for k1 in illumina_info_hash.keys():
			out[k1]={}
			for k2 in illumina_info_hash[k1].keys():
				out[k1][k2]={}
				for k3 in illumina_info_hash[k1][k2]:
					if not k3[3] in out[k1][k2].keys():					out[k1][k2][k3[3]]=0
					out[k1][k2][k3[3]]+=1
		return out
	def overlap_calcu(list1,list2,match_sv_type):
		#eg of list1=['chr22', 10717135, 10717147, 'INS', 0, 1]
		#eg of list2=['chr22', 10527787, 10700173, 'CNV', '0', 7]
		if not list1[0]==list2[0]: return 0
		if match_sv_type=='TRUE':	
			if not list1[3].upper()==list2[3].upper(): return 0
		if list2[1]>list1[2]: return 0
		if list1[1]>list2[2]: return 0
		return float(sorted(list1[1:3]+list2[1:3])[2]-sorted(list1[1:3]+list2[1:3])[1])/float(max([list1[2]-list1[1],list2[2]-list2[1]]))
	def overlap_calcu_based_on_one_SV(list1,list2,match_sv_type):
		#eg of list1=['chr22', 10717135, 10717147, 'INS', 0, 1]
		#eg of list2=['chr22', 10527787, 10700173, 'CNV', '0', 7]
		#overlap percentage is reported based on list2
		if float(list2[2]-list2[1])==0:
			print(list2)
		else:
			if not list1[0]==list2[0]: return 0
			if match_sv_type=='TRUE':	
				if not list1[3].upper()==list2[3].upper(): return 0
			if list2[1]>list1[2]: return 0
			if list1[1]>list2[2]: return 0
			return float(sorted(list1[1:3]+list2[1:3])[2]-sorted(list1[1:3]+list2[1:3])[1])/float(list2[2]-list2[1])
	def overlap_based_on_pacbio_calculate(illumina_info_hash,pacbio_merged_info,file_out,chromo):
		file_initiate(file_out)
		fo=open(file_out,'a')
		out_list={}
		for k1 in illumina_info_hash.keys():
			out_list[k1]=[]
			for k2 in pacbio_merged_info.keys():
				if k2 in illumina_info_hash[k1].keys():
					for k3 in pacbio_merged_info[k2]:
						if not k3[2]-k3[1]>0:	
							print(k3)
							continue
						temp=[]
						for k4 in illumina_info_hash[k1][k2]:
							overlap_rc=overlap_calcu_based_on_one_SV(k4,k3,'FALSE')
							if overlap_rc>0:							temp.append(overlap_rc)
						if temp==[]:	temp=[0]
						out_list[k1]+=temp
		for k1 in out_list.keys():		print(' '.join([str(i) for i in [chromo,k1,k2,','.join([str(i) for i in out_list[k1]])]]),file=fo)
		fo.close()
	def overlap_hash_derive(illumina_info_hash,pacbio_merged_info):
		overlap_hash={}
		rec_all=[]
		for k1 in illumina_info_hash.keys():
			overlap_hash[k1]={}
			for k2 in illumina_info_hash[k1].keys():
				if k2 in pacbio_merged_info.keys():
					[overlap_hash[k1][k2],rec]=bp_off_decide_both_end_separate_based_on_RC(illumina_info_hash[k1][k2],pacbio_merged_info[k2])
					rec_all+=rec
		return [overlap_hash,rec_all]
	def write_rec(file_out,rec):
		file_initiate(file_out)
		fo=open(file_out,'a')
		for x in rec:
			print(x,file=fo)
		fo.close()
	def write_SV_report_number_stat(file_out,SV_report_number_stat):
		fo=open(file_out,'w')
		for k1 in SV_report_number_stat.keys():
			for k2 in SV_report_number_stat[k1].keys():
				for k3 in SV_report_number_stat[k1][k2].keys():
					for k4 in SV_report_number_stat[k1][k2][k3].keys():
						if 'ALU' in k4:
							print(' '.join([str(i) for i in [k1,k2,k3,'DEL',SV_report_number_stat[k1][k2][k3][k4]]]),file=fo)
						elif 'L1' in k4:
							print(' '.join([str(i) for i in [k1,k2,k3,'DEL',SV_report_number_stat[k1][k2][k3][k4]]]),file=fo)
						else:
							print(' '.join([str(i) for i in [k1,k2,k3,k4,SV_report_number_stat[k1][k2][k3][k4]]]),file=fo)
		fo.close()
	def write_total_stat_to_write_to_file(ppre,total_stat_to_write,file_out=ppre+'STEP1_bp_off_ILL_vs_PB_both_end_separate_RC50.txt'):
		fo=open(file_out,'w')
		for k1 in total_stat_to_write.keys():
			for k2 in total_stat_to_write[k1].keys():
				for k3 in range(len(total_stat_to_write[k1][k2][0])):
					print([k1,k2,k3])
					print(' '.join([str(i) for i in [k1,k2,total_stat_to_write[k1][k2][0][k3],total_stat_to_write[k1][k2][1][k3]]]),file=fo)
		fo.close()
	def ILL_rec_readin(chromo,ppre,ILL_file_name):
		fin=open(ppre+ILL_file_name)
		out={}
		for line in fin:
			pin=line.strip().split()
			if 'START' in pin: continue
			if not pin[0]==chromo: continue
			if not pin[5] in out.keys():	out[pin[5]]={}
			if not pin[6] in out[pin[5]].keys():	out[pin[5]][pin[6]]=[]
			out[pin[5]][pin[6]].append([pin[0]]+[int(pin[1]),int(pin[2])]+pin[3:5])
		fin.close()
		return out
	def PB_rec_readin(chromo,ppre,PB_file_name):
		fin=open(ppre+PB_file_name)
		out={}
		for line in fin:
			pin=line.strip().split()
			if 'START' in pin: continue
			if not pin[0]==chromo: continue
			if not pin[5] in out.keys():	out[pin[5]]=[]
			out[pin[5]].append([pin[0]]+[int(pin[1]),int(pin[2])]+pin[3:5])
		fin.close()
		return out
	def main_1(ppre):
		[total_off_bp_stat,SV_report_number_stat]=[{},{}]
		for k1 in ['chr'+str(i) for i in range(23)[1:]]+['chrX']:
			illumina_info_hash=ILL_rec_readin(k1,ppre,'STEP0_ILL_Calls.Non_INS.bed')
			pacbio_merged_info=PB_rec_readin(k1,ppre,'STEP0_Pacbio_Calls.Non_INS.bed')
			#illumina_info_hash=illumina_reported_readin(k1,ppre)
			SV_report_number_stat[k1]=illumina_caller_stat(illumina_info_hash)
			#pacbio_merged_info=pacbio_merged_readin(pacbio_path,k1)
			[total_off_bp_stat[k1],rec]=overlap_hash_derive(illumina_info_hash,pacbio_merged_info)
			overlap_based_on_pacbio_calculate(illumina_info_hash,pacbio_merged_info,ppre+'STEP1_consensus_bp_CI_original.RO_based_on_PB'+'.txt',k1)
			write_rec(ppre+'STEP1_consensus_bp_CI_original'+'.txt',rec)
		write_SV_report_number_stat(ppre+'STEP1_ILL_Caller_SVtype_NUM_report.txt',SV_report_number_stat)
		total_stat_to_write={}
		for k1 in total_off_bp_stat.keys():
			for k2 in total_off_bp_stat[k1].keys():
				if not k2 in total_stat_to_write.keys():	total_stat_to_write[k2]={}
				for k3 in total_off_bp_stat[k1][k2].keys():
					if not k3 in total_stat_to_write[k2].keys():	total_stat_to_write[k2][k3]=[[],[]]
					total_stat_to_write[k2][k3][0]+=total_off_bp_stat[k1][k2][k3][0]
					total_stat_to_write[k2][k3][1]+=total_off_bp_stat[k1][k2][k3][1]
		write_total_stat_to_write_to_file(ppre,total_stat_to_write)
		print (ppre+'STEP1_consensus_bp_CI_original.RO_based_on_PB'+'.txt')
		print (ppre+'STEP1_consensus_bp_CI_original'+'.txt')
		print (ppre+'STEP1_ILL_Caller_SVtype_NUM_report.txt')
	def file_initiate(file_out,header_list=[]):
		if not os.path.isfile(file_out):
			fo=open(file_out,'w')
			if not header_list==[]:		print >>fo, '\t'.join(header_list)
			fo.close()

main_1(ppre)


os.system(r''' Rscript scripts/R1.calcu_CI_INS.R %s'''%(ppre))
os.system(r''' Rscript scripts/R2.calcu_CI_SVs.R %s'''%(ppre))


###################Step2. cluster illumina events###################
if Process_Steps[3]=='Cluster_ILL':
	def bp_pair_to_cons_le_ri(bp_pairs,bp_cons_hash_le,bp_cons_hash_ri):
		out=[]
		for k1 in bp_pairs:
			temp=[]
			if not k1[0] in bp_cons_hash_le.keys() or not k1[1] in bp_cons_hash_ri.keys(): continue
			temp.append(bp_cons_hash_le[k1[0]])
			temp.append(bp_cons_hash_ri[k1[1]])
			if temp[1]-temp[0]>size_cff:
				out.append(temp)
		return unify_list(sorted(out))
	def bp_pair_to_consensus(bp_pairs,bp_reported_consensus_hash):
		out=[]
		for k1 in bp_pairs:
			temp=[]
			for k2 in k1:
				temp.append(bp_reported_consensus_hash[k2])
			out.append(sorted(temp))
		return unify_list(out)
	def bp_to_bprange_by_CI(bp,CI_list_single):
		#eg of CI_list=[-143.9, 203.8]
		out=[bp-int(max([abs(i) for i in CI_list_single])),bp+int(max([abs(i) for i in CI_list_single]))]
		return out
	def CI_list_readin(CI_list_file):
		fin=open(CI_list_file)
		out={}
		for line in fin:
			pin=line.strip().split()
			out[pin[0]]=[float(i) for i in pin[1:]]
		fin.close()
		return out
	def cluster_illumina_according_to_pacbio(chromo,CI_list,ILL_file_name):
		#eg of chromo='chr22'
		illumina_info_hash=ILL_rec_readin(chromo,ppre,ILL_file_name)
		write_stat_illumina_info_hash(illumina_info_hash,chromo,ppre+'different_algorithm_SVreport_summary.stat')
		pacbio_merged_info=PB_rec_readin(chromo,ppre,PB_file_name)
		out_match_all={}
		out_other_all={}
		for k1 in pacbio_merged_info.keys():
			out_match={}
			out_other={}
			for k2 in illumina_info_hash.keys():
				if k1 in illumina_info_hash[k2].keys():
					for k3 in pacbio_merged_info[k1]:
						if 'INS' in k3: continue
						pb_key='_'.join([str(i) for i in k3])
						if not pb_key in out_match.keys():	out_match[pb_key]=[]
						for k4 in illumina_info_hash[k2][k1]:
							if k4[1]>k3[2]: break
							if k4[2]<k3[1]: continue
							if not k1 in CI_list.keys():	continue
							num_match_bps=two_sv_list_CI_compare(k3,k4,CI_list[k2])
							if num_match_bps>0: out_match[pb_key].append(k4+[k2])
							else:
								if not k4[0] in out_other.keys():	out_other[k4[0]]={}
								if not k4[1] in out_other[k4[0]].keys():	out_other[k4[0]][k4[1]]={}
								if not k4[2] in out_other[k4[0]][k4[1]].keys():	out_other[k4[0]][k4[1]][k4[2]]=[]
								out_other[k4[0]][k4[1]][k4[2]].append(k4[3:]+[k2])
			out_match_all[k1]=out_match
			out_other_all[k1]=out_other
		return [out_match_all,out_other_all]
	def cluster_illumina_only_left_right_bp(pb_other,CI_list):
		ill_bp_left_rec={}
		ill_bp_right_rec={}
		for k1 in pb_other.keys():
			ill_bp_left_rec[k1]={}
			ill_bp_right_rec[k1]={}
			for k2 in pb_other[k1].keys():
				for k3 in pb_other[k1][k2].keys():
					for k4 in pb_other[k1][k2][k3]:
						if not k2 in ill_bp_left_rec[k1].keys():	ill_bp_left_rec[k1][k2]=[]
						if not k3 in ill_bp_right_rec[k1].keys():	ill_bp_right_rec[k1][k3]=[]
						ill_bp_left_rec[k1][k2].append(k4[-1])
						ill_bp_right_rec[k1][k3].append(k4[-1])
		out={}
		out['left']={}
		ill_bp_rec=ill_bp_left_rec
		for ka in ill_bp_rec.keys():
			consensus_bp_rec=[]
			consensus_range_rec=[]
			consensus_sv_finder_rec=[]
			for kb in sorted(ill_bp_rec[ka].keys()):
				for kc in unify_list(ill_bp_rec[ka][kb]):
					bp_range=bp_to_bprange_by_CI(kb,CI_list[kc])
					if consensus_range_rec==[]:
						consensus_range_rec.append(bp_range)
						consensus_bp_rec.append([kb])
						consensus_sv_finder_rec.append([kc])
					else:
						if bp_range[0]>consensus_range_rec[-1][1]:	
							consensus_range_rec.append(bp_range)
							consensus_bp_rec.append([kb])
							consensus_sv_finder_rec.append([kc])
						elif not bp_range[0]>consensus_range_rec[-1][1] and not bp_range[1]<consensus_range_rec[-1][0]:
							consensus_range_rec[-1]=merge_bp_range(bp_range,consensus_range_rec[-1])
							consensus_sv_finder_rec[-1].append(kc)
							consensus_bp_rec[-1].append(kb)
			consensus_bp=[]
			for i in range(len(consensus_bp_rec)):
				consensus_bp.append(consensus_bp_derive(consensus_bp_rec[i],consensus_range_rec[i]))
			out['left'][ka]=[consensus_bp,consensus_bp_rec,consensus_range_rec,consensus_sv_finder_rec]
		out['right']={}
		ill_bp_rec=ill_bp_right_rec
		for ka in ill_bp_rec.keys():
			consensus_bp_rec=[]
			consensus_range_rec=[]
			consensus_sv_finder_rec=[]
			for kb in sorted(ill_bp_rec[ka].keys()):
				for kc in unify_list(ill_bp_rec[ka][kb]):
					bp_range=bp_to_bprange_by_CI(kb,CI_list[kc])
					if consensus_range_rec==[]:
						consensus_range_rec.append(bp_range)
						consensus_bp_rec.append([kb])
						consensus_sv_finder_rec.append([kc])
					else:
						if bp_range[0]>consensus_range_rec[-1][1]:	
							consensus_range_rec.append(bp_range)
							consensus_bp_rec.append([kb])
							consensus_sv_finder_rec.append([kc])
						elif not bp_range[0]>consensus_range_rec[-1][1] and not bp_range[1]<consensus_range_rec[-1][0]:
							consensus_range_rec[-1]=merge_bp_range(bp_range,consensus_range_rec[-1])
							consensus_sv_finder_rec[-1].append(kc)
							consensus_bp_rec[-1].append(kb)
			consensus_bp=[]
			for i in range(len(consensus_bp_rec)):
				consensus_bp.append(consensus_bp_derive(consensus_bp_rec[i],consensus_range_rec[i]))
			out['right'][ka]=[consensus_bp,consensus_bp_rec,consensus_range_rec,consensus_sv_finder_rec]
		out_new=out['left']
		for k1 in out['right'].keys():
			out_new[k1][0]+=out['right'][k1][0]
			out_new[k1][1]+=out['right'][k1][1]
			out_new[k1][2]+=out['right'][k1][2]
			out_new[k1][3]+=out['right'][k1][3]
		return out_new
	def cluster_ill_only_bp_le_ri(pb_other,CI_list):
		ill_bp_le={}
		ill_bp_ri={}
		for k1 in pb_other.keys():
			ill_bp_le[k1]={}
			ill_bp_ri[k1]={}
			for k2 in pb_other[k1].keys():
				for k3 in pb_other[k1][k2].keys():
					for k4 in pb_other[k1][k2][k3]:
						if not k2 in ill_bp_le[k1].keys():	ill_bp_le[k1][k2]=[]
						if not k3 in ill_bp_ri[k1].keys():	ill_bp_ri[k1][k3]=[]
						if not k4[-1] in ill_bp_le[k1][k2]:	ill_bp_le[k1][k2].append(k4[-1])
						if not k4[-1] in ill_bp_ri[k1][k3]:	ill_bp_ri[k1][k3].append(k4[-1])
		out_left=ill_ori_bp_hash_to_cons_freq(ill_bp_le,CI_list)
		out_right=ill_ori_bp_hash_to_cons_freq(ill_bp_ri,CI_list)
		return [out_left,out_right]
	def cluster_illumina_only_event_bp(pb_other,CI_list):
		ill_bp_rec={}
		for k1 in pb_other.keys():
			ill_bp_rec[k1]={}
			for k2 in pb_other[k1].keys():
				for k3 in pb_other[k1][k2].keys():
					for k4 in pb_other[k1][k2][k3]:
						if not k2 in ill_bp_rec[k1].keys():	ill_bp_rec[k1][k2]=[]
						if not k3 in ill_bp_rec[k1].keys():	ill_bp_rec[k1][k3]=[]
						ill_bp_rec[k1][k2].append(k4[-1])
						ill_bp_rec[k1][k3].append(k4[-1])
		out={}
		for ka in ill_bp_rec.keys():
			consensus_bp_rec=[]
			consensus_range_rec=[]
			consensus_sv_finder_rec=[]
			for kb in sorted(ill_bp_rec[ka].keys()):
				for kc in unify_list(ill_bp_rec[ka][kb]):
					bp_range=bp_to_bprange_by_CI(kb,CI_list[kc])
					if consensus_range_rec==[]:
						consensus_range_rec.append(bp_range)
						consensus_bp_rec.append([kb])
						consensus_sv_finder_rec.append([kc])
					else:
						if bp_range[0]>consensus_range_rec[-1][1]:	
							consensus_range_rec.append(bp_range)
							consensus_bp_rec.append([kb])
							consensus_sv_finder_rec.append([kc])
						elif not bp_range[0]>consensus_range_rec[-1][1] and not bp_range[1]<consensus_range_rec[-1][0]:
							consensus_range_rec[-1]=merge_bp_range(bp_range,consensus_range_rec[-1])
							consensus_sv_finder_rec[-1].append(kc)
							consensus_bp_rec[-1].append(kb)
			consensus_bp=[]
			for i in range(len(consensus_bp_rec)):
				consensus_bp.append(consensus_bp_derive(consensus_bp_rec[i],consensus_range_rec[i]))
			out[ka]=[consensus_bp,consensus_bp_rec,consensus_range_rec,consensus_sv_finder_rec]
		return out
	def cluster_illumina_only_event_sv(pb_other,CI_list,file_out,file_out_CI):
		consensus_bp=cluster_illumina_only_event_bp(pb_other,CI_list)
		for chromo in consensus_bp.keys():
			if chromo in pb_other.keys():
				consensus_CI=unify_list(consensus_bp[chromo][2])
				cluster_Illumina_independently(pb_other,CI_list,file_out)
				write_consensus_CI(file_out_CI,consensus_CI,chromo)
	def cluster_Illumina_independently_le_ri(pb_other,CI_list,file_out):
		[cons_bp_le,cons_bp_ri]=cluster_ill_only_bp_le_ri(pb_other,CI_list)
		pb_reverse_other=pb_other_to_pb_reverse_other(pb_other)
		for chromo in cons_bp_le.keys():
			if chromo in pb_other.keys():
				bp_pairs=pb_other_hash_to_bp_clusters(pb_other,chromo)
				bp_cluster=pb_other_hash_to_ordered_list(pb_other,chromo)
				bp_cons_hash_le=consensus_bp_to_hash(cons_bp_le,chromo)
				bp_cons_hash_ri=consensus_bp_to_hash(cons_bp_ri,chromo)
				consensus_pairs=bp_pair_to_cons_le_ri(bp_pairs,bp_cons_hash_le,bp_cons_hash_ri)
				consensus_full_info=consensus_pairs_to_consensus_full_info_le_ri(consensus_pairs,cons_bp_le,cons_bp_ri,chromo,pb_other,pb_reverse_other)
				consensus_new_full= consensus_full_info_supp_1(consensus_full_info,pb_other,chromo)
				write_consensus_new_full(consensus_new_full,file_out,chromo)
	def cluster_Illumina_independently_INS(pb_other,CI_INS,file_out):
		[cons_bp_le,cons_bp_ri]=cluster_ill_only_bp_le_ri(pb_other,CI_INS)
		consensus_new_full=[]
		for chromo in cons_bp_le.keys():
			if chromo in pb_other.keys():
				for k1 in range(len(cons_bp_le[chromo][0])):
					temp=[cons_bp_le[chromo][0][k1],cons_bp_le[chromo][1][k1],cons_bp_le[chromo][2][k1],cons_bp_le[chromo][3][k1],[]]
					for k2 in pb_other[chromo][temp[0]].keys():
						for k3 in pb_other[chromo][temp[0]][k2]:
							temp[-1].append([temp[0],k2]+k3)
					temp_new=[[i,i] for i in temp]
					consensus_new_full.append(temp_new)
				write_consensus_new_full(consensus_new_full,file_out,chromo)
	def cluster_Illumina_independently(pb_other,CI_list,file_out):
		#eg of file_out=ppre+'consensus_BP_Calculated_by_CI.all_info.report'
		consensus_bp=cluster_illumina_only_event_bp(pb_other,CI_list)
		#consensus_bp=cluster_illumina_only_left_right_bp(pb_other,CI_list)
		for chromo in consensus_bp.keys():
			if chromo in pb_other.keys():
				bp_pairs=pb_other_hash_to_bp_clusters(pb_other,chromo)
				bp_cluster=pb_other_hash_to_ordered_list(pb_other,chromo)
				bp_reported_consensus_hash=consensus_bp_to_hash(consensus_bp,chromo)
				consensus_pairs=bp_pair_to_consensus(bp_pairs,bp_reported_consensus_hash)
				consensus_full_info=consensus_pairs_to_consensus_full_info(consensus_pairs,consensus_bp,chromo)
				consensus_new_full= consensus_full_info_supp_1(consensus_full_info,pb_other,chromo)
				write_consensus_new_full(consensus_new_full,file_out,chromo)
	def consensus_bp_derive(bp_list,range_list):
		if len(bp_list)==0:	return	bp_list[0]
		else:
			range_mean=numpy.median(range_list)
			dis_bp=[abs(i-range_mean) for i in bp_list]
			return bp_list[dis_bp.index(min(dis_bp))]
	def consensus_bp_to_mate_list(consensus_bp_list,pb_link_hash):
		consensus_mate_list=[]
		for ka in consensus_bp_list:
			kb=[]
			for i in ka:	kb+=pb_link_hash[i]
			consensus_mate_list.append(kb)
		return consensus_mate_list
	def consensus_pairs_to_consensus_full_info(consensus_pairs,consensus_bp,chromo):
		consensus_full_info=[]
		for k1 in sorted(consensus_pairs):
			temp=[k1,[],[],[]]
			for k2 in k1:
				pos_index=consensus_bp[chromo][0].index(k2)
				temp[1].append(consensus_bp[chromo][1][pos_index])
				temp[2].append(consensus_bp[chromo][2][pos_index])
				temp[3].append(consensus_bp[chromo][3][pos_index])
			consensus_full_info.append(temp)
		return consensus_full_info
	def consensus_pair_modify(temp,pb_other,pb_reverse_other,chromo):
		temp_both=[]
		for k1 in temp[1][0]:
			for k2 in pb_other[chromo][k1].keys():
				for k3 in pb_other[chromo][k1][k2]:
					temp_both.append([k1,k2]+k3)
		for k1 in temp[1][1]:
			for k2 in pb_reverse_other[chromo][k1].keys():
				for k3 in pb_reverse_other[chromo][k1][k2]:
					temp_both.append([k1,k2]+k3)
		out=[]
		for x in temp_both:
			if x[0] in temp[1][0] and x[1] in temp[1][1] and x[-1] in temp[3][0] and x[-1] in temp[3][1]:
				if not x in out:	out.append(x)
		out.sort()
		temp_new=[temp[0],[[i[0] for i in out],[i[1] for i in out]],temp[2],[[i[-1] for i in out],[i[-1] for i in out]]]
		return temp_new
	def pb_other_to_pb_reverse_other(pb_other):
		out={}
		for k1 in pb_other.keys():	
			out[k1]={}
			for k2 in pb_other[k1].keys():
				for k3 in pb_other[k1][k2].keys():
					if not k3 in out[k1].keys():		out[k1][k3]={}
					if not k2 in out[k1][k3].keys():	out[k1][k3][k2]=[]
					for k4 in pb_other[k1][k2][k3]:
						if not k4 in out[k1][k3][k2]:out[k1][k3][k2].append(k4)
		return out
	def consensus_pairs_to_consensus_full_info_le_ri(consensus_pairs,cons_bp_le,cons_bp_ri,chromo,pb_other,pb_reverse_other):
		consensus_full_info=[]
		for k1 in sorted(consensus_pairs):
			temp=[k1,[],[],[]]
			k2=k1[0]
			pos_index=cons_bp_le[chromo][0].index(k2)
			temp[1].append(cons_bp_le[chromo][1][pos_index])
			temp[2].append(cons_bp_le[chromo][2][pos_index])
			temp[3].append(cons_bp_le[chromo][3][pos_index])
			k2=k1[1]
			pos_index=cons_bp_ri[chromo][0].index(k2)
			temp[1].append(cons_bp_ri[chromo][1][pos_index])
			temp[2].append(cons_bp_ri[chromo][2][pos_index])
			temp[3].append(cons_bp_ri[chromo][3][pos_index])			
			temp_new=consensus_pair_modify(temp,pb_other,pb_reverse_other,chromo)
			consensus_full_info.append(temp_new)
		return consensus_full_info
	def consensus_full_info_supp_1(consensus_full_info,pb_other,chromo):
		pb_other_supp=pb_other_hash_supp(pb_other)
		out=[]
		for k1 in consensus_full_info:
			temp=[[],[]]
			temp[0]=query_bp_info(k1[1][0],pb_other,chromo,pb_other_supp)
			temp[1]=query_bp_info(k1[1][1],pb_other,chromo,pb_other_supp)
			out.append(k1+[temp])
		return out
	def consensus_bp_to_hash(consensus_bp,chromo):
		out={}
		for i in range(len(consensus_bp[chromo][1])):
			for k2 in unify_list(consensus_bp[chromo][1][i]):
				if not k2 in out.keys():	out[k2]=consensus_bp[chromo][0][i]
		return out
	def ill_ori_bp_hash_to_cons_freq(ill_bp_rec,CI_list):
		#CI_list=CI_list_readin(CI_list_file)
		out={}
		for ka in ill_bp_rec.keys():
			consensus_bp_rec=[]
			consensus_range_rec=[]
			consensus_sv_finder_rec=[]
			for kb in sorted(ill_bp_rec[ka].keys()):
				for kc in unify_list(ill_bp_rec[ka][kb]):
					if not kc in CI_list: continue
					bp_range=bp_to_bprange_by_CI(kb,CI_list[kc])
					if consensus_range_rec==[]:
						consensus_bp_rec.append([kb])
						consensus_range_rec.append(bp_range)
						consensus_sv_finder_rec.append([kc])
					else:
						if not bp_range[0]>max(consensus_range_rec[-1]):
							consensus_bp_rec[-1].append(kb)
							consensus_range_rec[-1]+=bp_range
							consensus_sv_finder_rec[-1].append(kc)
						else:
							consensus_bp_rec.append([kb])
							consensus_range_rec.append(bp_range)
							consensus_sv_finder_rec.append([kc])
			rec=-1
			out[ka]=[[],[],[],[]]
			for x in consensus_range_rec:
				rec+=1
				if len(x)==2:	
						out[ka][0].append(consensus_bp_rec[rec][0])
						out[ka][1].append(consensus_bp_rec[rec])
						out[ka][2].append(x)
						out[ka][3].append(consensus_sv_finder_rec[rec])
				else:
					merged_range=consensus_range_list_collaps(x)
					if merged_range[1]==x:	
						out[ka][0].append(consensus_bp_rec[rec][0])
						out[ka][1].append(consensus_bp_rec[rec])
						out[ka][2].append(merged_range[0])
						out[ka][3].append(consensus_sv_finder_rec[rec])
					else:
						[cons_CI,supp_CI]=merged_range
						x_new_rec=[[x[2*i],x[2*i+1]] for i in range(int(len(x)/2))]
						for y in range(len(cons_CI)):
							if len(supp_CI[y])>0:
								raw_range_rec=[[supp_CI[y][2*i],supp_CI[y][2*i+1]] for i in range(int(len(supp_CI[y])/2))]
								raw_bp_rec=[consensus_bp_rec[rec][x_new_rec.index(i)] for i in raw_range_rec]
								raw_caller_rec=[consensus_sv_finder_rec[rec][x_new_rec.index(i)] for i in raw_range_rec]
								cons_bp=consensus_bp_assign(raw_bp_rec,cons_CI[y])
								out[ka][0].append(cons_bp)
								out[ka][1].append(raw_bp_rec)
								out[ka][2].append(cons_CI[y])
								out[ka][3].append(raw_caller_rec)
		return out
	def consensus_bp_assign(raw_bp_rec,cons_CI_int):
		out=[]
		num=[]
		for x in raw_bp_rec:
			if not x in out:
				out.append(x)
				num.append(1)
			else:
				num[out.index(x)]+=1
		out_1=[out[i] for i in range(len(out)) if num[i]==max(num)]
		if len(out_1)==1:	return out_1[0]
		else:
			out_dis=[abs(x-numpy.median(cons_CI_int)) for x in out_1]
			out_2=[out_1[i] for i in range(len(out_1)) if out_dis[i]==min(out_dis)]
			return out_2[0]
	def CI_list_group(CI_list_temp):
		# eg of CI_list_temp=[[29940308, 29940334], [29940305, 29940337], [29940310, 29940332], [29940282, 29940368], [29940286, 29940366], [29940275, 29940377], [29940312, 29940340], [29940311, 29940341], [29941192, 29941218], [29941189, 29941221], [29941162, 29941248], [29941165, 29941245], [29941195, 29941217], [29941155, 29941257], [29941192, 29941220], [29941191, 29941221]]
		out=[]
		for i in CI_list_temp:
			if out==[]:	out.append(i)
			else:
				if i[0]<max(out[-1]):	out[-1]+=i
				else:	out.append(i)
		out_new=[[[i[2*j],i[2*j+1]] for j in range(int(len(i)/2))] for i in out]
		return out_new
	def consensus_CI_derive(CI_list_temp):
		# eg of CI_list_temp=[[29940308, 29940334], [29940305, 29940337], [29940310, 29940332], [29940282, 29940368], [29940286, 29940366], [29940275, 29940377], [29940312, 29940340], [29940311, 29940341], [29941192, 29941218], [29941189, 29941221], [29941162, 29941248], [29941165, 29941245], [29941195, 29941217], [29941155, 29941257], [29941192, 29941220], [29941191, 29941221]]
		le_s=[i[0] for i in CI_list_temp]
		ri_s=[i[1] for i in CI_list_temp]
		if numpy.max(le_s)<numpy.min(ri_s):		return [numpy.max(le_s),numpy.min(ri_s)]
		else:									return 'error'
	def consensus_CI_derive_2(consensus_range_list):
		#ef pg consensus_range_list=[49973679, 49973705, 49973681, 49973761, 49973695, 49973797, 49973727, 49973829, 49973763, 49973843, 49973769, 49973849, 49973788, 49973868, 49973826, 49973906, 49973815, 49973917, 49973827, 49973913]
		le_s=[consensus_range_list[2*i] for i in range(int(len(consensus_range_list)/2))]
		ri_s=[consensus_range_list[2*i+1] for i in range(int(len(consensus_range_list)/2))]
		if numpy.max(le_s)<numpy.min(ri_s):		return [numpy.max(le_s),numpy.min(ri_s)]
		else:									return 'error'
	def consensus_CI_peak_finding(cov_list,consensus_range_list):
		out=[]
		for i in range(len(cov_list))[1:-1]:
			if cov_list[i+1]-cov_list[i]<0 and cov_list[i]-cov_list[i-1]>0:
				out.append(i)
		if out==[]:out.append(cov_list[cov_list.index(max(cov_list))])
		supp_num=[cov_list[i] for i in out]
		cons_reg=[[sorted(consensus_range_list)[i],sorted(consensus_range_list)[i+1]] for i in out]
		return [cons_reg,supp_num]
	def consensus_CI_supp_CI_assign(consensus_range_list,cons_reg,supp_num):
		supp_hash={}
		for i in range(len(supp_num)):
			if not supp_num[i] in supp_hash.keys():	supp_hash[supp_num[i]]=[]
			if not cons_reg[i] in supp_hash[supp_num[i]]: supp_hash[supp_num[i]].append(cons_reg[i])
		cons_CI_out=[]
		supp_CI_list=[]
		supp_CI_sub_rec=[]
		for i in sorted(supp_hash.keys())[::-1]:
			cons_CI_list=supp_hash[i]
			for k in cons_CI_list:	
				supp_CI_sub=[]	
				for j in range(int(len(consensus_range_list)/2)):
					if consensus_range_list[2*j]<k[1] and consensus_range_list[2*j+1]>k[0]:
						if not consensus_range_list[2*j] in supp_CI_sub_rec and not consensus_range_list[2*j+1] in supp_CI_sub_rec:
							supp_CI_sub+=[consensus_range_list[2*j],consensus_range_list[2*j+1]]
							supp_CI_sub_rec+=[consensus_range_list[2*j],consensus_range_list[2*j+1]]
				supp_CI_list.append(supp_CI_sub)
				cons_CI_out.append(k)
		return [cons_CI_out,supp_CI_list]
	def consensus_CI_stack_up(consensus_range_list):
		le_s=[consensus_range_list[2*i] for i in range(int(len(consensus_range_list)/2))]
		ri_s=[consensus_range_list[2*i+1] for i in range(int(len(consensus_range_list)/2))]
		cov_list=[]
		cov_rec=0
		for i in sorted(consensus_range_list):
			if i in le_s and i in ri_s:	continue
			else:
				if i in le_s:	cov_rec+=1
				elif i in ri_s:	cov_rec-=1
			cov_list.append(cov_rec)
		return cov_list
	def consensus_range_list_collaps(consensus_range_list):
		#eg of consensus_range_list=[49973679, 49973705, 49973681, 49973761, 49973695, 49973797, 49973727, 49973829, 49973763, 49973843, 49973769, 49973849, 49973788, 49973868, 49973826, 49973906, 49973815, 49973917, 49973827, 49973913]
		test1=consensus_CI_derive_2(consensus_range_list)
		if not test1=='error':	return [test1,consensus_range_list] #ideal situation, a consensus CI can be derived by taking the overlap
		else:	#few long CIs
			cov_list=consensus_CI_stack_up(consensus_range_list)
			[cons_reg,supp_num]=consensus_CI_peak_finding(cov_list,consensus_range_list)
			[cons_CI,supp_CI]=consensus_CI_supp_CI_assign(consensus_range_list,cons_reg,supp_num)
			extra_CI=consensus_CI_extra_CI_extract(consensus_range_list,supp_CI)
			extra_CI_assignment=consensus_CI_extra_CI_assign(cons_CI,extra_CI)
			for ia in range(len(extra_CI)):	supp_CI[extra_CI_assignment[ia]]+=extra_CI[ia]
			return [cons_CI,supp_CI]
	def consensus_CI_extra_CI_assign(cons_CI,extra_CI):
		out=[]
		for ia in extra_CI:
			temp=[abs(numpy.median(ia)-numpy.median(j)) for j in cons_CI]
			out.append(temp.index(min(temp)))
		return out
	def consensus_CI_extra_CI_extract(consensus_range_list,supp_CI):
		temp=[]
		for ia in supp_CI:
			for ib in range(int(len(ia)/2)):
				temp.append([ia[2*ib],ia[2*ib+1]])
		all_temp=[]
		for ib in range(int(len(consensus_range_list)/2)):
				all_temp.append([consensus_range_list[2*ib],consensus_range_list[2*ib+1]])
		out=[]
		for i in all_temp:
			if not i in temp:	out.append(i)
		return out
	def ill_ori_bp_hash_to_cons_closest_median(ill_bp_rec,CI_list):
		out={}
		for ka in ill_bp_rec.keys():
			consensus_bp_rec=[]
			consensus_range_rec=[]
			consensus_sv_finder_rec=[]
			for kb in sorted(ill_bp_rec[ka].keys()):
				for kc in unify_list(ill_bp_rec[ka][kb]):
					bp_range=bp_to_bprange_by_CI(kb,CI_list[kc])
					if consensus_range_rec==[]:
						consensus_bp_rec.append([kb])
						consensus_range_rec.append(bp_range)
						consensus_sv_finder_rec.append([kc])
					else:
						if bp_range[0]>consensus_range_rec[-1][1]:	
							consensus_bp_rec.append([kb])
							consensus_range_rec.append(bp_range)
							consensus_sv_finder_rec.append([kc])
						elif not bp_range[0]>consensus_range_rec[-1][1] and not bp_range[1]<consensus_range_rec[-1][0]:
							consensus_range_rec[-1]=merge_bp_range(bp_range,consensus_range_rec[-1])
							consensus_sv_finder_rec[-1].append(kc)
							consensus_bp_rec[-1].append(kb)
			consensus_bp=[]
			for i in range(len(consensus_bp_rec)):
				consensus_bp.append(consensus_bp_derive(consensus_bp_rec[i],consensus_range_rec[i]))
			out[ka]=[consensus_bp,consensus_bp_rec,consensus_range_rec,consensus_sv_finder_rec]
		return out
	def illumina_info_hash_reorganize(illumina_info_hash):
		out={}
		child_names=['HG00512','HG00513','HG00514','HG00731','HG00732','HG00733','NA19238','NA19239','NA19240']
		for k1 in illumina_info_hash.keys():
			for k2 in illumina_info_hash[k1].keys():
				if k2 in child_names:
					if not k2 in out.keys():	out[k2]={}
					for k3 in illumina_info_hash[k1][k2]:
						if not k3[0] in out[k2]:	out[k2][k3[0]]={}
						if not k3[1] in out[k2][k3[0]].keys():	out[k2][k3[0]][k3[1]]={}
						if not k3[2] in out[k2][k3[0]][k3[1]].keys():	out[k2][k3[0]][k3[1]][k3[2]]=[]
						if not k3[3:]+[k1] in out[k2][k3[0]][k3[1]][k3[2]]:	out[k2][k3[0]][k3[1]][k3[2]].append(k3[3:]+[k1])
		return out
	def merge_bp_range(bp_range1,bp_range2):
		out=[max([bp_range1[0],bp_range2[0]]), min([bp_range1[1],bp_range2[1]])]
		return out
	def pb_link_index(pb_link):
		out={}
		for k1 in pb_link:
			for k2 in k1:
				if not k2 in out.keys():	out[k2]=[]
				out[k2]+=[i for i in k1 if not i==k2]
		return out
	def pb_link_stat(pb_link):
		out={}
		for k1 in pb_link:
			if not len(k1) in out.keys():
				out[len(k1)]=0
			out[len(k1)]+=1
	def pb_other_merge(sample_list,illumina_info_new_hash):
		out={}
		for i in sample_list:
			for j in illumina_info_new_hash[i].keys():
				if not j in out.keys():	out[j]={}
				for k in illumina_info_new_hash[i][j].keys():
					if not k in out[j].keys():	out[j][k]={}
					for m in illumina_info_new_hash[i][j][k].keys():
						if not m in out[j][k].keys():	out[j][k][m]=[]
						for n in illumina_info_new_hash[i][j][k][m]:
							if not n[:-1]+[i,n[-1]] in out[j][k][m]:out[j][k][m].append(n[:-1]+[i,n[-1]])
		return out
	def pb_other_hash_to_ordered_list(pb_other,k1):
		pb_link=[]
		for k2 in sorted(pb_other[k1].keys()):
			temp=[k2]
			for k3 in sorted(pb_other[k1][k2].keys()):
				temp.append(k3)
			pb_link.append(temp)
		return pb_link
	def pb_other_hash_to_bp_clusters(pb_other,k1):
		pb_cluster=[]
		for ka in pb_other[k1].keys():
			for kb in pb_other[k1][ka].keys():
				if not [ka,kb] in pb_cluster:
					pb_cluster.append([ka,kb])
		return pb_cluster
	def pb_other_hash_supp(pb_other):
		out={}
		for k1 in pb_other.keys():
			out[k1]={}
			for k2 in pb_other[k1].keys():
				for k3 in pb_other[k1][k2].keys():
					if not k3 in out[k1].keys():
						out[k1][k3]=[]
					out[k1][k3].append(k2)
		return out
	def query_bp_info(bp_list,pb_other,chromo,pb_other_supp):
		out=[]
		for k2 in bp_list:
			if k2 in pb_other[chromo].keys():
				for k3 in pb_other[chromo][k2].keys():
					for k4 in pb_other[chromo][k2][k3]:
						out.append([k2,k3]+k4)
			elif k2 in pb_other_supp[chromo].keys():
				for k3 in pb_other_supp[chromo][k2]:
					for k4 in pb_other[chromo][k3].keys():
						if not k4==k2:	continue
						for k5 in pb_other[chromo][k3][k4]:
							out.append([k3,k4]+k5)
		return out
	def stat_illumina_info_hash(illumina_info_hash):
		out_stat={}
		for k1 in illumina_info_hash.keys():	
			out_stat[k1]={}
			for k2 in illumina_info_hash[k1].keys():
				out_stat[k1][k2]={}
				for k3 in illumina_info_hash[k1][k2]:
					if not k3[3] in out_stat[k1][k2].keys():	out_stat[k1][k2][k3[3]]=0
					out_stat[k1][k2][k3[3]]+=1
		return out_stat
	def sv_info_retrive(temp2,pb_other_sub):
		#eg of temp2=[[50434305, 50434319, 50434319], [50434771, 50434865, 50434896, 50434865, 50434896]]
		out=[]
		for k1 in temp2[0]:
			if k1 in pb_other_sub.keys():
				for k2 in temp2[1]:
					if k2 in pb_other_sub[k1].keys():
						out+=pb_other_sub[k1][k2]
		for k1 in temp2[1]:
			if k1 in pb_other_sub.keys():
				for k2 in temp2[0]:
					if k2 in pb_other_sub[k1].keys():
						out+=pb_other_sub[k1][k2]
		return unify_list(out)
	def two_list_merge(list1,list2):
		if list1==list2:	return [list1]
		else:
			flag_list1=0
			for i in list1:
				if not i in list2:	flag_list1+=1
			if flag_list1==0:	return [list2]
			flag_list2=0
			for i in list2:
				if not i in list1:	flag_list2+=1
			if flag_list2==0:	 return [list1]
			return [list1,list2]
	def two_sv_list_CI_compare(list1,list2,CI_range,match_sv_type='TRUE'):
		#eg of list1=['chr22', 11946180, 11946241, 'DEL', '1|0']
		#eg of list2=['chr22', 10784610, 10785884, 'del', '1/1']
		#eg of CI_range=[-59.0, 32.0]
		if match_sv_type=='TRUE':
			if not list1[3].upper()==list2[3].upper(): return 0
		out_list=[0,0]
		if not overlap_calcu(list1,list2,match_sv_type)<0.5:
			if abs(list1[1]-list2[1])<max([abs(i) for i in CI_range]): out_list[0]=1
			if abs(list1[2]-list2[2])<max([abs(i) for i in CI_range]): out_list[1]=1
		return sum(out_list)
	def unify_list(list):
		out=[]
		for i in list:
			if not i in out:	out.append(i)
		return out
	def write_all_info_illumina_cluster_multi_sample(chromo,CI_list,file_out_prefix,ppre,sample_list):
		illumina_info_hash=ILL_rec_readin(chromo,ppre,ILL_file_name)
		illumina_info_new_hash=illumina_info_hash_reorganize(illumina_info_hash)
		pb_other=pb_other_merge(sample_list,illumina_info_new_hash)
		cluster_Illumina_independently(pb_other,CI_list,file_out_prefix+'.'+'_'.join(sample_list)+'.report')
	def write_all_info_illumina_cluster_multi_sample_le_ri(chromo,CI_list,file_out_prefix,ppre,sample_list):
		illumina_info_hash=ILL_rec_readin(chromo,ppre,ILL_file_name)
		illumina_info_new_hash=illumina_info_hash_reorganize(illumina_info_hash)
		pb_other=pb_other_merge(sample_list,illumina_info_new_hash)
		cluster_Illumina_independently_le_ri(pb_other,CI_list,file_out_prefix+'.'+'_'.join(sample_list)+'.report')
	def ILL_rec_readin(chromo,ppre,ILL_file_name):
		fin=open(ppre+ILL_file_name)
		out={}
		for line in fin:
			pin=line.strip().split()
			if 'START' in pin: continue
			if not pin[0]==chromo: continue
			if not pin[5] in out.keys():	out[pin[5]]={}
			if not pin[6] in out[pin[5]].keys():	out[pin[5]][pin[6]]=[]
			out[pin[5]][pin[6]].append([pin[0]]+[int(pin[1]),int(pin[2])]+pin[3:5])
		fin.close()
		return out
	def write_all_info_illumina_cluster_multi_sample_INS(chromo,CI_INS,file_out_prefix,ppre,sample_list):
		illumina_info_hash=ILL_rec_readin(chromo,ppre,ILL_INS_name)
		illumina_info_new_hash=illumina_info_hash_reorganize(illumina_info_hash)
		pb_other=pb_other_merge(sample_list,illumina_info_new_hash)
		cluster_Illumina_independently_INS(pb_other,CI_INS,file_out_prefix+'.'+'_'.join(sample_list)+'.report')
	def write_all_info_illumina_cluster(chromo,CI_list,file_out_prefix,ppre):
		illumina_info_hash=ILL_rec_readin(chromo,ppre,ILL_file_name)
		illumina_info_new_hash=illumina_info_hash_reorganize(illumina_info_hash)
		for k2 in illumina_info_new_hash.keys():
			pb_other=illumina_info_new_hash[k2]
			cluster_Illumina_independently(pb_other,CI_list,file_out_prefix+'.'+k2+'.report')
	def write_all_info_illumina_cluster_le_ri(chromo,CI_list,file_out_prefix,ppre):
		illumina_info_hash=ILL_rec_readin(chromo,ppre,ILL_file_name)
		illumina_info_new_hash=illumina_info_hash_reorganize(illumina_info_hash)
		for k2 in illumina_info_new_hash.keys():
			pb_other=illumina_info_new_hash[k2]
			cluster_Illumina_independently_le_ri(pb_other,CI_list,file_out_prefix+'.'+k2+'.report')
	def write_bp_match(pb_match,file_out):
		file_initiate(file_out)
		fo=open(file_out,'a')
		for k1 in pb_match.keys():
			key_new=':'.join(k1.split('_'))
			info_new=[':'.join([str(i) for i in j]) for j in pb_match[k1]]
			print('\t'.join([key_new]+info_new),file=fo)
		fo.close()
	def write_consensus_CI(file_out,consensus_CI,chromo):
		file_initiate(file_out)
		fo=open(file_out,'a')
		for k1 in consensus_CI:
			print('\t'.join([str(i) for i in [chromo]+k1]),file=fo)
		fo.close()
	def write_consensus_bp(file_out,consensus_bp_for_write,chromo):
		file_initiate(file_out)
		fo=open(file_out,'a')
		for ka in consensus_bp_for_write:
			print('\t'.join([str(i) for i in [chromo]+ka[:2]+[';'.join([':'.join([str(j) for j in i]) for i in ka[2:-2]])]+[':'.join([str(j) for j in i]) for i in ka[-2:]]]),file=fo)
		fo.close()
	def write_consensus_new_full(consensus_new_full,file_out,chromo):
		file_initiate(file_out)
		fo=open(file_out,'a')
		print(' '.join(['chr','start-consensus','end-consensus','star-CI','end-CI','Illumina_caller_reported_records','reported_bp','reported_algorithm']),file=fo)
		for k1 in consensus_new_full:
			print('\t'.join([str(x) for x in [chromo]+k1[0]+['-'.join([str(j) for j in k]) for k in k1[2]]+[';'.join([':'.join([str(i) for i in j]) for j in k]) for k in k1[-1]]+[';'.join([':'.join([str(i) for i in j]) for j in k1[1]])]+[';'.join([':'.join([str(i) for i in j]) for j in k1[3]])] ]),file=fo)
		fo.close()
	def write_stat_illumina_info_hash(illumina_info_hash,chromo,file_out):
		stat_hash=stat_illumina_info_hash(illumina_info_hash)
		file_initiate(file_out)
		fo=open(file_out,'a')
		for k1 in stat_hash.keys():
			for k2 in stat_hash[k1].keys():
				for k3 in stat_hash[k1][k2].keys():
					print(' '.join([str(i) for i in [chromo,k1,k2,k3,stat_hash[k1][k2][k3]]]),file=fo)
		fo.close()
	def main_2():
		chromos=chromos_readin_canonical(ref)
		for k1 in chromos[:23]:
			print(k1)
			write_all_info_illumina_cluster_multi_sample_le_ri(k1,CI_list,ppre+'consensus_BP_Calculated_by_CI.all_info',ppre,['HG00512','HG00513','HG00514','HG00731','HG00732','HG00733','NA19238','NA19239','NA19240'])
			[pb_match_all,pb_other_all]=cluster_illumina_according_to_pacbio(k1,CI_list,ILL_file_name)
			for k2 in pb_match_all.keys():
				pb_match=pb_match_all[k2]
				pb_other=pb_other_all[k2]
				file_out=ppre+'consensus_BP_Calculated_by_CI.Pacbio_Excluded.'+k2+'.bed'
				file_out_CI=ppre+'consensus_BP_Calculated_by_CI.Pacbio_Excluded.BP_CI.'+k2+'.bed'
				cluster_illumina_only_event_sv(pb_other,CI_list,file_out,file_out_CI)
				write_bp_match(pb_match,ppre+'consensus_BP_Calculated_by_CI.Pacbio_based.'+k2+'.bed')
	def main_2b():
		chromos=chromos_readin_canonical(ref)
		for k1 in chromos[:23]:
			print(k1)
			write_all_info_illumina_cluster_multi_sample_INS(k1,CI_INS,ppre+'consensus_BP_Calculated_by_CI.INS',ppre,['HG00512','HG00513','HG00514','HG00731','HG00732','HG00733','NA19238','NA19239','NA19240'])

global CI_list_file,CI_list,CI_INS_file,CI_INS,ILL_file_name
CI_list_file=ppre+'consensus_bp_CI_90Quantile.txt'
CI_list=CI_list_readin(CI_list_file)
CI_INS_file=ppre+'consensus_INS_CI_90Quantile.txt'
CI_INS=CI_list_readin(CI_INS_file)
ILL_file_name='STEP0_ILL_Calls.Non_INS.bed'
ILL_INS_name='STEP0_ILL_Calls.INS.bed'
PB_file_name='STEP0_Pacbio_Calls.Non_INS.bed'
main_2()
main_2b()



###################Steo2b. clean up the consensus_BP_Calculated_by_CI.all_info.HG00514.report file########################
if Process_Steps[4]=='Cleanup_Merged_SVs':
	def bp_CIs_to_consensus_CIs(left_CIs,range_cff=1000):
		#set 1000 as the default cutoff, if all CIs are below 1000, or all above 1000, it's fine;
		#but if some are below but some are above, larger ones will be removed
		CI_categorize=[[],[]]
		CI_categorize[0]=[i for i in left_CIs if i[1]-i[0]<range_cff]
		CI_categorize[1]=[i for i in left_CIs if not i[1]-i[0]<range_cff]
		if len(CI_categorize[0])>0 and len(CI_categorize[1])>0:
			new_CIs=CI_categorize[0]
		else:
			new_CIs=left_CIs
		out=[]
		for i in sorted(new_CIs):
			if out==[]:	out.append(i)
			else:
				if i[0]>out[-1][-1]:	out.append(i)
				else:	
					new_CI=[max([out[-1][0],i[0]]),min([out[-1][1],i[1]])]
					out[-1]=new_CI
		return out
	def bp_list_to_distance_rank(bp_list,median):
		out=[abs(i-median) for i in bp_list]
		out_hash={}
		for i in range(len(out)):
			if not out[i] in out_hash.keys():out_hash[out[i]]=[]
			out_hash[out[i]].append(bp_list[i])
		out_new=[]
		for i in sorted(out_hash.keys()):
			out_new+=out_hash[i]
		return out_new
	def bp_list_to_abundancy_rank(bp_list,bp_median):
		out=[]
		freq=[]
		for i in bp_list:
			if not i in out:	
				out.append(i)
				freq.append(1)
			else:
				freq[out.index(i)]+=1
		out_hash={}
		for i in range(len(freq)):
			if not freq[i] in out_hash.keys():			out_hash[freq[i]]=[]
			out_hash[freq[i]].append(out[i])
		out_new=[]
		for i in sorted(out_hash.keys())[::-1]:	
			out_new+=bp_list_to_distance_rank(out_hash[i],bp_median)
		return out_new
	def bp_list_to_consensus_bp(left_bps,left_CIs, deviate_dis=10):
		#sub func of sv_list_to_CIs
		left_consensus_CI=bp_CIs_to_consensus_CIs(left_CIs)
		left_bps_order1_frequency=bp_list_to_abundancy_rank(left_bps,numpy.median(left_consensus_CI))
		out=[]
		for i in left_consensus_CI:
			for j in left_bps_order1_frequency:
				if not j<i[0]-deviate_dis and not j>i[1]+deviate_dis:
					out.append([j,i])
					break
		return out
	def bp_to_BPCI(bp,sv,CI_list):
		return [bp-min([abs(i) for i in CI_list[sv]]),bp+min([abs(i) for i in CI_list[sv]])]
	def CI_list_to_CI(CI_input):
		return (max([i[0] for i in CI_input]),min([i[1] for i in CI_input]))
	def calcu_most_frequent_number_as_consensus_bp(support_hash,record_hash):
		support_modified_hash={}
		for k2 in sorted(support_hash.keys())[::-1]:
			for k3 in support_hash[k2]:
				if not int(k3.split('_')[0]) in support_modified_hash.keys():				
					support_modified_hash[int(k3.split('_')[0])]=record_hash[int(k3.split('_')[0])][:3]+[k2]
				if k3.split('_')[1]=='le':
					if not ';' in record_hash[int(k3.split('_')[0])][-1]: 
						freq_bp=most_frequent_number_pick(record_hash[int(k3.split('_')[0])][-2].split(';')[0].split(':'))
						support_modified_hash[int(k3.split('_')[0])][1]=freq_bp
				elif k3.split('_')[1]=='ri':
					if not ';' in record_hash[int(k3.split('_')[0])][-2]: 
						freq_bp=most_frequent_number_pick(record_hash[int(k3.split('_')[0])][-2].split(';')[1].split(':'))
						support_modified_hash[int(k3.split('_')[0])][2]=freq_bp
		return support_modified_hash
	def clean_record_step1_merge_overlap(record_new_hash,organize_hash):
		out=[]
		out_index=[]
		for k1 in organize_hash[1]:
			brief_rec=record_new_hash[k1][:3]+[record_to_name_of_supportive_algorithm(record_new_hash[k1])]
			if out==[]:	
				out.append(brief_rec)
				out_index.append([k1])
			else:
				if not brief_rec[0]==out[-1][0]:	
					out.append(brief_rec)
					out_index.append([k1])
				elif not brief_rec[-1]==out[-1][-1]:	
					out.append(brief_rec)
					out_index.append([k1])
				elif brief_rec[1]>out[-1][-2]:
					out.append(brief_rec)
					out_index.append([k1])
				else:
					out[-1]+=brief_rec
					out_index[-1]+=[k1]
		merged_list=clean_record_step1_sub1_merge_bp(out,out_index,record_new_hash)
		return merged_list
	def clean_record_step1_sub1_merge_bp(out,out_index,record_new_hash):
		merged_list=[]
		for i in range(len(out_index)):
			if len(out_index[i])==1:
				merged_list.append(record_new_hash[out_index[i][0]][:-1]+record_to_name_of_supportive_algorithm(record_new_hash[out_index[i][0]]))
			else:
				temp=[record_new_hash[j] for j in out_index[i]]
				merged_temp=clean_record_step1_sub2_merge_recs(temp)
				merged_list.append(merged_temp[:-1]+record_to_name_of_supportive_algorithm(merged_temp))
		return merged_list
	def clean_record_step1_sub2_merge_recs(temp):
		temp=unify_list(temp)
		bp_info_hash_left={}
		bp_info_hash_right={}
		for k1 in temp:
			if not k1[1] in bp_info_hash_left.keys():
				bp_info_hash_left[k1[1]]=[k1[0],k1[1],k1[3],k1[5],k1[7].split(';')[0],k1[8].split(';')[0]]
			if not k1[2] in bp_info_hash_right.keys():
				bp_info_hash_right[k1[2]]=[k1[0],k1[2],k1[4],k1[6],k1[7].split(';')[1],k1[8].split(';')[1]]
		temp2=[bp_info_hash_left[min(bp_info_hash_left.keys())],bp_info_hash_right[max(bp_info_hash_right.keys())]]
		out=[temp2[0][0],temp2[0][1],temp2[1][1],temp2[0][2],temp2[1][2],temp2[0][3],temp2[1][3],';'.join([temp2[0][4],temp2[1][4]]),';'.join([temp2[0][5],temp2[1][5]])]
		return out
	def clean_record_step1b_add_info_to_cleaned_hash_1(cleaned_hash_1,record_new_hash):
		cleaned_hash_1b={}
		for x in cleaned_hash_1.keys():
			if x>1:
				temp=[]
				for y in cleaned_hash_1[x]:				temp.append(record_new_hash[y])
				cleaned_hash_1b[x]=temp
			else:
				cleaned_hash_1b[x]=cleaned_hash_1[x]
		return cleaned_hash_1b
	def chromo_contig_prepare(ref_file):
		out=[]
		out.append('##reference='+ref_file)
		if os.path.isfile(ref_file+'.fai'):
			fin=open(ref_file+'.fai')
			for line in fin:
				pin=line.strip().split()
				if not '_' in pin[0] and not '-' in pin[0]:
					out.append('##contig=<ID='+pin[0]+',length='+pin[1]+'>')
			fin.close()
		return out
	def get_sv_len(vcf_info):
		for x in vcf_info[7].split(';'):
			if x.split('=')[0]=='SVLEN':
				return int(x.split('=')[1])
	def left_cons_bp_modify(left_cons_bp):
		#sub func of sv_list_to_CIs
		out=[]
		for x in left_cons_bp:
			if x[0]<x[1][0]:
				x[1][0]=x[0]
			elif x[0]>x[1][1]:
				x[1][1]=x[0]
			out.append([x[0],[int(i) for i in x[1]]])
		return out
	def most_frequent_number_pick(number_list):
		number_list=[int(i) for i in number_list]
		item_list=[]
		num_list=[]
		for i in number_list:
			if not i in item_list:
				item_list.append(i)
				num_list.append(1)
			else:
				num_list[item_list.index(i)]+=1
		out=[item_list[i] for i in range(len(item_list)) if num_list[i]==max(num_list)]
		if len(out)==1:	return out[0]
		else:
			mean_pos=numpy.median(number_list)
			mean_dis=[abs(i-mean_pos) for i in out]
			out_new=[out[i] for i in range(len(out)) if mean_dis[i]==min(mean_dis)]
			return out_new[0]
	def order_cleeaned_hash_by_pos(cleaned_hash_1b):
		#cleaned_hash_1b is a hash, with keys stands for number of supportive of each consensus bp.[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
		pos_hash={}
		for k1 in cleaned_hash_1b.keys():
			rec=-1
			for k2 in cleaned_hash_1b[k1]:
				rec+=1
				if not k2[0] in pos_hash.keys():	pos_hash[k2[0]]={}
				if not int(k2[1]) in pos_hash[k2[0]].keys():	pos_hash[k2[0]][int(k2[1])]={}
				if not int(k2[2]) in pos_hash[k2[0]][int(k2[1])].keys():	pos_hash[k2[0]][int(k2[1])][int(k2[2])]=[]
				if not [k1,rec] in pos_hash[k2[0]][int(k2[1])][int(k2[2])]: pos_hash[k2[0]][int(k2[1])][int(k2[2])].append([k1,rec])
		ref='/data/talkowski/xuefang/data/reference/GRCh38.1KGP/GRCh38_full_analysis_set_plus_decoy_hla.fa'
		chromos=chromos_readin(ref)
		pos_list=[]
		for k1 in chromos:
			if k1 in pos_hash.keys():
				for k2 in sorted(pos_hash[k1].keys()):
					for k3 in sorted(pos_hash[k1][k2].keys()):
						for k4 in pos_hash[k1][k2][k3]:
							if len(pos_list)>0:
								if not k4[:3]==pos_list[-1][:3]:			
									pos_list.append([k1,k2,k3,k4])
							else:											
								pos_list.append([k1,k2,k3,k4])
		return pos_list
	def pick_most_freq_bp(BP_list):
		bps=[]
		list=[]
		for i in BP_list:
			if not i in bps:
				bps.append(i)
				list.append(1)
			else:
				list[bps.index(i)]+=1
		out1=[bps[i] for i in range(len(bps)) if list[i]==max(list)]
		dis=[abs(i-numpy.median(BP_list)) for i in out1]
		out2=[out1[i] for i in range(len(out1)) if dis[i]==min(dis)]
		return out2[0]
	def record_hash_modify_replace_closest_with_most_freq_bp(record_hash,support_modified_hash):
		for k1 in support_modified_hash.keys():
			record_hash[k1][:3]=support_modified_hash[k1][:3]
		return record_hash
	def record_to_name_of_supportive_algorithm(record):
		supportive_algorithms=[i.split(':') for i in record[-1].split(';')]
		supportive_algorithms_both=unify_list(supportive_algorithms[0]+supportive_algorithms[1])
		supportive_algorithms_both.sort()
		return supportive_algorithms_both
	def record_new_hash_clean1_remove_overlap_called_by_single_algorithm(record_new_hash,individual_name,stat_out1=ppre+'consensus_bp_number_of_supportive_callers.stat'):
		organize_hash={}
		for k1 in sorted(record_new_hash.keys()):
			supportive_algorithms=[i.split(':') for i in record_new_hash[k1][-1].split(';')]
			supportive_algorithms_both=unify_list(supportive_algorithms[0]+supportive_algorithms[1])
			if not len(supportive_algorithms_both) in organize_hash.keys():
				organize_hash[len(supportive_algorithms_both)]=[]
			organize_hash[len(supportive_algorithms_both)].append(k1)
		file_initiate(stat_out1)
		fo=open(stat_out1,'a')
		print (' '.join(['individual','#supportive_algorithms','#consensus_events']),file=fo)
		for k1 in organize_hash.keys():
			print (' '.join([str(i) for i in [individual_name,k1,len(organize_hash[k1])]]),file=fo)
		fo.close()
		merged_list=clean_record_step1_merge_overlap(record_new_hash,organize_hash) #merge consensus SVs that overlap with each other and all contributed by one caller 
		organize_hash[1]=merged_list
		return organize_hash
	def ref_base_readin(ref,chr,pos):
		fin=os.popen(r'''samtools faidx %s %s:%s-%s'''%(ref,chr,pos,pos))
		pin=fin.readline().strip().split()
		pin=fin.readline().strip().split()
		fin.close()
		if len(pin)>0:	return pin[0]
		else:				return 'N'
	def report_file_readin(file_in):
		support_hash={}
		record_hash={}
		rec=-1
		fin=open(file_in)
		for line in fin:
			pin=line.strip().split()
			if not 'start-consensus' in pin:
				if len(pin)<9: continue
				rec+=1
				record_hash[rec]=pin
				supportive_algorithms=[i.split(':') for i in pin[-1].split(';')]
				if not len(unify_list(supportive_algorithms[0])) in support_hash.keys():	support_hash[len(unify_list(supportive_algorithms[0]))]=[]
				support_hash[len(unify_list(supportive_algorithms[0]))].append(str(rec)+'_le')
				if not len(unify_list(supportive_algorithms[1])) in support_hash.keys():	support_hash[len(unify_list(supportive_algorithms[1]))]=[]
				support_hash[len(unify_list(supportive_algorithms[1]))].append(str(rec)+'_ri')
		fin.close()
		return [support_hash,record_hash]
	def sv_type_collect(record):
		#ef of record=['chr1', 1, 136094, '-2492-2494', '133601-138587', '1:41000:DUP:18:dCGH_filt', '71326:136094:DUP:29:dCGH_filt', '1;136094', 'dCGH_filt']
		SV_type_list=[i.split(':')[2] for i in record[5].split(';')]+[i.split(':')[2] for i in record[6].split(';')]
		return '<'+';'.join(unify_list(SV_type_list))+'>'
	def sv_support_collect(record):
		out1=[','.join(unify_list(record[5].split(';'))),','.join(unify_list(record[6].split(';')))]
		if not out1[0]==out1[1]:
			out=['INFO_POS='+out1[0],'INFO_END='+out1[1]]
		else:
			out=['INFO_POS='+out1[0],'INFO_END='+'.']
		return out
	def sv_info_collect(sv_info,sv_type,kb):
		out=['MERGE_TYPE='+sv_type.replace('<','').replace('>','').replace(';',','),'END='+str(sv_info[2]),'SVLEN='+str(sv_info[2]-sv_info[1]),'CIPOS='+','.join([str(i) for i in [int(i)-sv_info[1] for i in sv_info[3].split('-')]]),'CIEND='+','.join([str(i) for i in [int(i)-sv_info[2] for i in sv_info[4].split('-')]])]
		if not ';' in sv_info[-1]:	supp_algorithm=sv_info[-1]
		else:						supp_algorithm=','.join(record_to_name_of_supportive_algorithm(sv_info))
		out.append('NUM_CALLER='+str(len(supp_algorithm.split(','))))
		out.append('CALLER='+supp_algorithm)
		out+=sv_support_collect(sv_info)
		return ';'.join(out)
	def sv_info_modify_1(sv_info):
		sv_info[5]=';'.join(sorted(unify_list(sv_info[5].split(';'))))
		sv_info[6]=';'.join(sorted(unify_list(sv_info[6].split(';'))))
		if sv_info[1]==sv_info[2] and not ';' in sv_info[5] and not ';' in sv_info[6] and sv_info[5]==sv_info[6]:
			sv_info[1]=int(sv_info[5].split(':')[0])
			sv_info[2]=int(sv_info[5].split(':')[1])
		return sv_info
	def sv_list_to_CIs(left_svs,CI_list):
		#sub func of sv_info_modify_3
		sv_method=[i.split(':')[-1] for i in left_svs]
		#left bp
		left_bps=[int(i.split(':')[0]) for i in left_svs]
		left_CIs=[ [left_bps[i]+j for j in CI_list[sv_method[i]]]   for i in range(len(sv_method))]
		left_cons_bp=bp_list_to_consensus_bp(left_bps,left_CIs, deviate_dis=10)
		#right bp
		right_bps=[int(i.split(':')[1]) for i in left_svs]
		right_CIs=[ [right_bps[i]+j for j in CI_list[sv_method[i]]]   for i in range(len(sv_method))]
		right_cons_bp=bp_list_to_consensus_bp(right_bps,right_CIs, deviate_dis=10)
		if [] in [left_cons_bp,right_cons_bp]:
			return []
		else:
			return [left_cons_bp_modify(left_cons_bp),left_cons_bp_modify(right_cons_bp)]
	def sv_app_info_modify(left_new,left_svs,CI_list,deviate_dis=10):
		#sub func of sv_info_modify_3
		out_combine=[]
		left_CIs_from_CI_list=[[int(i.split(':')[0])+j for j in CI_list[i.split(':')[-1]]]  for i in left_svs]
		right_CIs_from_CI_list=[[int(i.split(':')[1])+j for j in CI_list[i.split(':')[-1]]]  for i in left_svs]
		for i in left_new[0]:
			for j in left_new[1]:
				out_combine.append([i,j])
		out_support=[]
		for i in out_combine:
			rec=-1
			out_support.append([])
			for j in left_CIs_from_CI_list:
				rec+=1
				k=right_CIs_from_CI_list[rec]
				if not i[0][0]<j[0]-deviate_dis and not i[0][0]>j[1]+deviate_dis and not i[1][0]<k[0]-deviate_dis and not i[1][0]>k[1]+deviate_dis :
					out_support[-1].append(left_svs[rec])
		new_support=[';'.join(sorted(unify_list(i))) for i in out_support]
		out_total=[]
		for i in range(len(out_combine)):
			out_total.append([out_combine[i][0][0],out_combine[i][1][0],'-'.join([str(k) for k in out_combine[i][0][1]]),'-'.join([str(k) for k in out_combine[i][1][1]]) , new_support[i],new_support[i],
				 ';'.join([':'.join([j.split(':')[0] for j in out_support[i]]),':'.join([j.split(':')[1] for j in out_support[i]])]),
					 ';'.join([':'.join([j.split(':')[-1] for j in out_support[i]]),':'.join([j.split(':')[-1] for j in out_support[i]])])])
		return out_total
	def sv_info_modify_3(sv_info,CI_list,CI_fold_cff=10):
		left_svs=unify_list(sv_info[5].split(';'))
		right_svs=unify_list(sv_info[6].split(';'))
		left_CI=[int(i) for i in sv_info[3].split('-')]
		right_CI=[int(i) for i in sv_info[4].split('-')]
		out=[]
		if float(right_CI[1]-right_CI[0])>0 and float(left_CI[1]-left_CI[0])/float(right_CI[1]-right_CI[0])>CI_fold_cff and float(left_CI[1]-left_CI[0])>100:
			left_new=sv_list_to_CIs(left_svs,CI_list)
			right_new=sv_list_to_CIs(right_svs,CI_list)
			if left_new==right_new and not left_new==[]:
				out+=[[sv_info[0]]+i for i in sv_app_info_modify(left_new,left_svs,CI_list)]
			else:
				if not left_new==[]:			out+=[[sv_info[0]]+i for i in sv_app_info_modify(left_new,left_svs,CI_list)]
				if not right_new==[]:			out+=[[sv_info[0]]+i for i in sv_app_info_modify(right_new,right_svs,CI_list)]
		elif float(left_CI[1]-left_CI[0])>0 and float(right_CI[1]-right_CI[0])/float(left_CI[1]-left_CI[0])>CI_fold_cff and float(right_CI[1]-right_CI[0])>100:
			left_new=sv_list_to_CIs(left_svs,CI_list)
			right_new=sv_list_to_CIs(right_svs,CI_list)
			if left_new==right_new and not left_new==[]:
				out+=[[sv_info[0]]+i for i in sv_app_info_modify(left_new,left_svs,CI_list)]
			else:
				if not left_new==[]:			out+=[[sv_info[0]]+i for i in sv_app_info_modify(left_new,left_svs,CI_list)]
				if not right_new==[]:			out+=[[sv_info[0]]+i for i in sv_app_info_modify(right_new,right_svs,CI_list)]
		else:
			out.append(sv_info)
		return out
	def sv_info_modify_2(sv_info,CI_list):
		if sv_info[1]==sv_info[2]:
			#do arbitrary CI=100 for now
			left_bps=[int(i.split(':')[0]) for i in sv_info[5].split(';')]
			left_callers=[i.split(':')[-1] for i in sv_info[5].split(';')]
			left_CI=CI_list_to_CI([bp_to_BPCI(left_bps[i],left_callers[i],CI_list) for i in range(len(left_bps))])
			right_bps=[int(i.split(':')[1]) for i in sv_info[6].split(';')]
			right_callers=[i.split(':')[-1] for i in sv_info[6].split(';')]
			right_CI=CI_list_to_CI([bp_to_BPCI(right_bps[i],right_callers[i],CI_list) for i in range(len(right_bps))])
			if max(left_bps)-min(left_bps)<100 and max(right_bps)-min(right_bps)<100:
				sv_info[1]=pick_most_freq_bp(left_bps)
				sv_info[2]=pick_most_freq_bp(right_bps)
				sv_info[3]='-'.join([str(i) for i in [min(left_bps),max(left_bps)]])
				sv_info[4]='-'.join([str(i) for i in [min(right_bps),max(right_bps)]])
				return sv_info
			else:
				return 'FAIL'
		else:
			left_bps=[int(i.split(':')[0]) for i in sv_info[5].split(';')]
			if len(left_bps)>3:
				left_CI=[min(left_bps),max(left_bps)]
				if left_CI[0]>int(sv_info[3].split('-')[0]) and left_CI[1]<int(sv_info[3].split('-')[1]):
					sv_info[3]='-'.join(str(i) for i in left_CI)
			right_bps=[int(i.split(':')[1]) for i in sv_info[6].split(';')]
			if len(right_bps)>3:
				right_CI=[min(right_bps),max(right_bps)]
				if right_CI[0]>int(sv_info[4].split('-')[0]) and right_CI[1]<int(sv_info[4].split('-')[1]):
					sv_info[4]='-'.join(str(i) for i in right_CI)
			return sv_info
	def sv_info_modify_4(sv_info):
		#make sure for singletons, the bp are the caller's bp
		#sv_info=['chr1', 24602041, 24602590, '24601988-24602074', '24602577-24602591', '24602041:24602591:DEL:1/0:NA19238:SVelter_UMich', '24602041:24602591:DEL:1/0:NA19238:SVelter_UMich', '24602041;24602591', 'SVelter_UMich;SVelter_UMich']
		le_bp=unify_list([i.split(':')[0] for i in sv_info[5].split(';')])
		ri_bp=unify_list([i.split(':')[1] for i in sv_info[6].split(';')])	
		if len(le_bp)==1 and not sv_info[1]==int(le_bp[0]): sv_info[1]=int(le_bp[0])
		if len(ri_bp)==1 and not sv_info[2]==int(ri_bp[0]): sv_info[2]=int(ri_bp[0])
		sv_info[1]=int(sv_info[1])
		sv_info[2]=int(sv_info[2])
		return sv_info
	def sv_info_qc(sv_info):
		if sv_info[-3]=='':return 'FAIL'
		if sv_info[-4]=='':return 'FAIL'
		if sv_info[-2].split(';') in [[''],['', '']]: return 'FAIL'
		if sv_info[-1].split(';') in [[''],['', '']]: return 'FAIL'
		if len(sv_info[3].split('-'))>2 or len(sv_info[4].split('-'))>2: return 'FAIL'
		return 'PASS'
	def write_bed_most_frequent_number_as_consensus_bp(record_hash,support_modified_hash,file_out):
		fo=open(file_out,'w')
		print (' '.join(['chr','closest_to_median_bp_left','closest_to_median_bp_right','chr','most_frequently_reported_bp_left','most_frequently_reported_bp_right','number_of_supportive_algorthms']), file=fo)
		for k1 in support_modified_hash.keys():
			print (' '.join([str(i) for i in record_hash[k1][:3]+support_modified_hash[k1]]),file=fo)
		fo.close()
	def remove_redun(cleaned_hash_1b):
		out={}
		for k1 in cleaned_hash_1b.keys():
			out[k1]=unify_list(cleaned_hash_1b[k1])
		return out
	def write_cleeaned_hash_to_vcf(cleaned_hash_1b,vcf_prefix,individual_name):
		clean_hash_1c=remove_redun(cleaned_hash_1b)
		ordered_pos_list=order_cleeaned_hash_by_pos(clean_hash_1c)
		write_VCF_header(vcf_prefix+'.vcf',ref)
		fo1=open(vcf_prefix+'.vcf','a')
		fo2=open(vcf_prefix+'.report','w')
		print ('\t'.join(['chr','start-consensus','end-consensus','Illumina_caller_reported_records_left','Illumina_caller_reported_records_right','star-CI;end-CI','reported_bp:left;right','reported_algorithm:left;right']),file=fo2)
		for ka in ordered_pos_list:
			for kb in ka[3:]:
				sv_info=clean_hash_1c[kb[0]][kb[1]]
				if sv_info[3][0]=='-' or sv_info[4][0]=='-': continue
				if sv_info_qc(sv_info)=='FAIL': continue
				sv_info=sv_info_modify_1(sv_info)
				sv_info=sv_info_modify_2(sv_info,CI_list)
				if sv_info=='FAIL': continue
				sv_info_list=sv_info_modify_3(sv_info,CI_list)
				for sv_info in sv_info_list:
					if sv_info=='FAIL': continue
					if sv_info_qc(sv_info)=='FAIL': continue
					sv_info=sv_info_modify_4(sv_info)
					print ('\t'.join([str(i) for i in sv_info[:3]+sv_info[5:7]+[';'.join(sv_info[3:5])]+sv_info[7:]]),file=fo2)
					sv_type=sv_type_collect(sv_info)
					vcf_info=sv_info[:2]+['.',ref_base_readin(ref,sv_info[0],sv_info[1]),sv_type,'.','.',sv_info_collect(sv_info,sv_type,kb)]
					print ('\t' .join([str(i) for i in vcf_info]),file=fo1 )
		fo1.close()
		fo2.close()
	def write_VCF_header(output_file,ref):
		fo=open(output_file,'w')
		print('##fileformat=VCFv4.1',file=fo)
		print('##fileDate='+time.strftime("%Y%m%d"),file=fo)
		ref_info=chromo_contig_prepare(ref)
		for x in ref_info:
			print (x,file=fo)
	   	#print('##reference=hg19',file=fo)
		print('##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">',file=fo)
		print('##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">',file=fo)
		print('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">',file=fo)
		print('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the variant">',file=fo)
		print('##INFO=<ID=MERGE_TYPE,Number=.,Type=String,Description="Type of structural variant">',file=fo)
		print('##INFO=<ID=NUM_CALLER,Number=1,Type=Integer,Description="=Number of algorithms supporting this variant">',file=fo)
		print('##INFO=<ID=CALLER,Number=.,Type=String,Description="Name of algorithms supporting this variant">',file=fo)
		print('##INFO=<ID=INFO_POS,Number=.,Type=String,Description="Calls by individual algorithm supporting left breakpoint of this variant">',file=fo)
		print('##INFO=<ID=INFO_END,Number=.,Type=String,Description="Calls by individual algorithm supporting right breakpoint of this variant, if different from INFO_POS, otherwise .">',file=fo)
		print('##FILTER=<ID=LowQual,Description="Score of final structural - Theoretical Score <-50">',file=fo)
		print('##ALT=<ID=DEL,Description="Deletion">',file=fo)
		print('##ALT=<ID=DUP,Description="Duplication">',file=fo)
		print('##ALT=<ID=INV,Description="Inversion">',file=fo)
		print('##ALT=<ID=DEL;DUP,Description="merged calls with individual callers report both deletion and duplication">',file=fo)
		print('##ALT=<ID=DEL;INV,Description="merged calls with individual callers report both deletion and inversion">',file=fo)
		print('##ALT=<ID=DUP;INV,Description="merged calls with individual callers report both duplication and inversion">',file=fo)
		print('##ALT=<ID=DEL;DUP;INV,Description="merged calls with individual callers report deletion, duplication and inversion">',file=fo)
		print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',file=fo)
		print('##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype quality">',file=fo)
		print('##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">',file=fo)
		print('##FORMAT=<ID=CNQ,Number=1,Type=Float,Description="Copy number genotype quality for imprecise events">',file=fo)
		print('\t'.join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']),file=fo)
		fo.close()
	def unify_list(list):
		out=[]
		for i in list:
			if not i in out:	
				out.append(i)
		return out
	def main_2b(sample_list,prefix,prefix_new,ppre):
		for k1 in sample_list:
			individual_name=k1
			file_in=ppre+prefix+'.'+k1+'.report'
			[support_hash,record_hash]=report_file_readin(file_in)
			support_modified_hash=calcu_most_frequent_number_as_consensus_bp(support_hash,record_hash)
			#write_bed_most_frequent_number_as_consensus_bp(record_hash,support_modified_hash,ppre+prefix+'.'+k1+'.closest.bp.vs.most_freq.bp.report')
			record_new_hash=record_hash_modify_replace_closest_with_most_freq_bp(record_hash,support_modified_hash)
			cleaned_hash_1=record_new_hash_clean1_remove_overlap_called_by_single_algorithm(record_new_hash,individual_name)
			cleaned_hash_1b=clean_record_step1b_add_info_to_cleaned_hash_1(cleaned_hash_1,record_new_hash)
			write_cleeaned_hash_to_vcf(cleaned_hash_1b,ppre+prefix_new+'.'+individual_name,individual_name)


prefix='consensus_BP_Calculated_by_CI.all_info'
prefix_new='consensus_BP_Calculated_by_CI.new.all_info'
CI_list_oth=CI_list_readin(CI_list_file)
CI_list_INS=CI_list_readin(CI_INS_file)
CI_list=CI_list_oth
for i in CI_list_INS.keys():
	if not i in CI_list.keys():	
		CI_list[i]=CI_list_INS[i]

main_2b(['HG00512_HG00513_HG00514_HG00731_HG00732_HG00733_NA19238_NA19239_NA19240'],prefix,prefix_new,ppre)

prefix_INS='consensus_BP_Calculated_by_CI.INS'
prefix_INS_new='consensus_BP_Calculated_by_CI.new.INS'
CI_list_INS=CI_list_readin(CI_INS_file)
CI_list=CI_list_INS
main_2b(['HG00512_HG00513_HG00514_HG00731_HG00732_HG00733_NA19238_NA19239_NA19240'],prefix_INS,prefix_INS_new,ppre)


#######vcf to sorted.vcf#######
if Process_Steps[5]=='QCs':
	############################QC1############################
	#######to clean up redundancy from vcf files #######
	#!python
	def bp_cluster(bp_list,dis_merge=20):
		out=[]
		bp_list.sort()
		for x in bp_list:
			if out==[]:	out.append([x])
			else:
				if x-out[-1][-1]<dis_merge:	out[-1].append(x)
				else:	out.append([x])
		return out
	def check_cluster(bp_cluster):
		out=[]
		for x in bp_cluster:
			if len(x)>1:
				out.append(x[-1]-x[0])
		return out
	def choose_cluster(bp_cluster):
		out=[]
		for x in bp_cluster:
			if len(x)>1:
				out.append(x)
		return out
	def final_check_0_readin_vcf(vcf_file):
		fin=open(vcf_file)
		chromos=chromos_readin(ref)
		info_hash={}
		for line in fin:
			pin=line.strip().split()
			if not pin[0][0]=='#':
				[sv_type,sv_pos]=chr_start_end_extract(pin)
				if not pin[0] in info_hash.keys():	info_hash[sv_pos[0]]={}
				if not sv_pos[1] in info_hash[sv_pos[0]].keys():	info_hash[sv_pos[0]][sv_pos[1]]={}
				if not sv_pos[2] in info_hash[sv_pos[0]][sv_pos[1]].keys():	info_hash[sv_pos[0]][sv_pos[1]][sv_pos[2]]=[]
				if not pin in info_hash[sv_pos[0]][sv_pos[1]][sv_pos[2]]: info_hash[sv_pos[0]][sv_pos[1]][sv_pos[2]].append(pin)
		fin.close()
		return info_hash
	def final_check_1_vcf(final_hash_0):
		rec=0
		for k1 in final_hash_0.keys():
			for k2 in final_hash_0[k1].keys():
				for k3 in final_hash_0[k1][k2].keys():
					if len(final_hash_0[k1][k2][k3])>1:
						rec+=1
						temp=pick_one_from_exact_overlap(final_hash_0[k1][k2][k3])
						final_hash_0[k1][k2][k3]=[temp]
		return [rec,final_hash_0]
	def final_check_2_vcf(clean_1_hash,dis_merge=20):
		#merge down svs with both bp within 5bp by default
		for k1 in clean_1_hash.keys():
			left_bp=sorted(clean_1_hash[k1].keys())
			right_bp=[]
			ri_le_hash={}
			for k2 in left_bp:
				right_bp+=clean_1_hash[k1][k2].keys()
				for k3 in clean_1_hash[k1][k2].keys():
					if not k3 in ri_le_hash.keys():	ri_le_hash[k3]=[k2]
					else:							ri_le_hash[k3]+=[k2]
			right_bp.sort()
			vcf_rec_group=[]
			left_cluster=bp_cluster(left_bp,dis_merge)
			left_group=choose_cluster(left_cluster)
			for i in left_group:
				vcf_rec_group.append([])
				for j in i:
					for k in clean_1_hash[k1][j].keys():
						for l in clean_1_hash[k1][j][k]:
							vcf_rec_group[-1].append(l)
			right_cluster=bp_cluster(right_bp,dis_merge)
			right_group=choose_cluster(right_cluster)
			for i in right_group:
				vcf_rec_group.append([])
				for j in i:
					for ka in ri_le_hash[j]:
						for kb in clean_1_hash[k1][ka][j]:
							vcf_rec_group[-1].append(kb)
			merged_vcf_rec=vcf_rec_group_down(vcf_rec_group)
			for i in range(len(merged_vcf_rec[0])):
				for j in merged_vcf_rec[1][i]:
					if j[1] in clean_1_hash[j[0]].keys() and j[2] in clean_1_hash[j[0]][j[1]].keys():
						del  clean_1_hash[j[0]][j[1]][j[2]]
				j_new=chr_start_end_extract(merged_vcf_rec[2][i])[1]
				if not j_new[1] in clean_1_hash[k1].keys(): clean_1_hash[k1][j_new[1]]={}
				if not j_new[2] in clean_1_hash[k1][j_new[1]].keys():	clean_1_hash[k1][j_new[1]][j_new[2]]=[]
				if not merged_vcf_rec[2][i] in clean_1_hash[k1][j_new[1]][j_new[2]]:	clean_1_hash[k1][j_new[1]][j_new[2]].append(merged_vcf_rec[2][i])
		return clean_1_hash
	def final_check_3_vcf(clean_2_hash):
		#remove merges if INFO_POS and INFO_END do not have overlap
		out={}
		for k1 in clean_2_hash.keys():
			out[k1]={}
			for k2 in clean_2_hash[k1].keys():
				for k3 in clean_2_hash[k1][k2].keys():
					for k4 in clean_2_hash[k1][k2][k3]:
						info_supp_oc=INFO_SUPP_overlap_check(k4)
						if not info_supp_oc=='error':
							if not k2 in out[k1].keys():	out[k1][k2]={}
							if not k3 in out[k1][k2].keys():	out[k1][k2][k3]=[]
							if not k4 in out[k1][k2][k3]: out[k1][k2][k3].append(k4)
						else:
							num_caller=num_caller_extract(k4)
		return out
	def INFO_SUPP_overlap_check(pin):
		out=[[],[]]
		for x in pin[7].split(';'):
			if x.split('=')[0]=='INFO_POS':
				out[0]+=x.split('=')[1].split(',')
			if x.split('=')[0]=='INFO_END':
				out[1]+=x.split('=')[1].split(',')
				if x.split('=')[1]=='.':	out[1]=out[0]
		overlap=[i for i in out[0] if i in out[1]]
		if overlap==[]:	return 'error'
		else:
			return pin
	def merge_sv_num_caller_extract(pin):
		out=0
		for x in pin[7].split(';'):
			if x.split('=')[0]=='NUM_CALLER':
				out=int(x.split('=')[1])
		return out
	def merge_sv_CI_extract(pin):
		out=[[],[]]
		for x in pin[7].split(';'):
			if x.split('=')[0]=='CIPOS':
				out[0]=[int(i) for i in x.split('=')[1].split(',')]
			elif x.split('=')[0]=='CIEND':
				out[1]=[int(i) for i in x.split('=')[1].split(',')]
		return out
	def merge_sv_CI_modify(CI_list):
		#CI range should be [neg, pos], if not, correct to 0
		out=[]
		if CI_list[0]>0:	out.append(0)
		else:	out.append(CI_list[0])
		if CI_list[1]<0:	out.append(0)
		else:	out.append(CI_list[1])
		return out
	def num_caller_extract(pin):
		out=0
		for i in pin[7].split(';'):
			if i.split('=')[0]=='NUM_CALLER':
				out=int(i.split('=')[1])
		return out
	def pick_one_from_exact_overlap(info_list):
		#eg of info_list=[['chr8', '112066381', '.', 'T', '<DEL>', '.', 'PASS', 'MERGE_TYPE=DEL;END=112066471;SVLEN=90;CIPOS=-13,13;CIEND=-53,53;NUM_CALLER=9;CALLER=GenomeSTRiP_CNV,SVelter_UMich,UCSD_Manta,UCSD_forestSV,VH,dCGH_filt,delly_illumina,liWGS_illumina,lumpy;INFO_POS=112065516:112077905:DEL:0/1:NA19239:liWGS_illumina,112065516:112077905:DEL:0/1:NA19240:liWGS_illumina,112066054:112077146:DEL:0/1:NA19239:GenomeSTRiP_CNV,112066054:112077146:DEL:0/1:NA19240:GenomeSTRiP_CNV,112066157:112080377:DEL:0/1:NA19239:dCGH_filt,112066157:112080377:DEL:0/1:NA19240:dCGH_filt,112066201:112077600:DEL:0/1:NA19239:UCSD_forestSV,112066201:112077600:DEL:0/1:NA19240:UCSD_forestSV,112066381:112077101:DEL:./.:HG00512:UCSD_Manta,112066381:112077101:DEL:./.:HG00513:UCSD_Manta,112066381:112077101:DEL:./.:HG00514:UCSD_Manta,112066381:112077101:DEL:./.:HG00731:UCSD_Manta,112066381:112077101:DEL:./.:HG00732:UCSD_Manta,112066381:112077101:DEL:./.:HG00733:UCSD_Manta,112066381:112077101:DEL:0/1:NA19239:UCSD_Manta,112066381:112077101:DEL:0/1:NA19240:UCSD_Manta,112066381:112077102:DEL:0/1:NA19239:delly_illumina,112066381:112077102:DEL:0/1:NA19240:delly_illumina,112066384:112077101:DEL:0/1:NA19239:lumpy,112066384:112077101:DEL:0/1:NA19240:lumpy,112066386:112066471:DEL:0/1:NA19238:SVelter_UMich,112066386:112077102:DEL:0/1:NA19239:VH,112066386:112077102:DEL:0/1:NA19240:VH,112066386:112077102:DEL:1/0:NA19239:SVelter_UMich,112066386:112077102:DEL:1/0:NA19240:SVelter_UMich;INFO_END=112066386:112066471:DEL:0/1:NA19238:SVelter_UMich'],['chr8', '112066381', '.', 'T', '<DEL>', '.', 'PASS', 'MERGE_TYPE=DEL;END=112066471;SVLEN=90;CIPOS=-13,7;CIEND=-53,33;NUM_CALLER=1;CALLER=SVelter_UMich;INFO_POS=112066386:112066471:DEL:0/1:NA19238:SVelter_UMich;INFO_END=.']]
		#Selection Criteria: 1. pick the call with more caller supportive; 2. if same number of caller, pick the call with smaller CIs
		num_caller_list=[merge_sv_num_caller_extract(pin) for pin in info_list]
		out_index_1=[i for i in range(len(num_caller_list)) if num_caller_list[i]==max(num_caller_list)-3]
		if len(out_index_1)==1: return info_list[out_index_1[0]]
		else:
			CI_range_list=[merge_sv_CI_extract(pin) for pin in info_list]
			CI_modify_list=[[merge_sv_CI_modify(i) for i in j] for j in CI_range_list]
			CI_size=[numpy.median([j[0][1]-j[0][0],j[1][1]-j[1][0]]) for j in CI_modify_list]
			out_index_2=[i for i in range(len(CI_size)) if CI_size[i]==min(CI_size)]
			return info_list[out_index_2[0]]
	def sub_dis_cluster(sub_dis,dis_merge=20):
		#eg of sub_dis=[['chr1', 248795399, 248799154], ['chr1', 248798394, 248799154]]
		#this function clusters event if both bp are within dis_merge
		out=[]
		for i in sub_dis:
			if out==[]: out.append([i])
			else:	
				flag=0
				for j in out:
					for k in j:
						if abs(i[1]-k[1])<dis_merge and abs(i[2]-k[2])<dis_merge:
							flag+=1
							j.append(i)
						if flag>0: break
					if flag>0: break
				if flag==0:	out.append([i])
		return out
	def supp_num_extract(pin):
		out=0
		for x in pin[7].split(';'):
			if x.split('=')[0]=='NUM_CALLER':	out=int(x.split('=')[1])
		return out
	def supp_algs_extract(pin):
		out=[]
		for x in pin[7].split(';'):
			if x.split('=')[0]=='INFO_POS':	
				if not x.split('=')[1]=='.':			out+=[i.split(':')[-1] for i in x.split('=')[1].split(',')]
		for x in pin[7].split(';'):
			if x.split('=')[0]=='INFO_END':	
				if not x.split('=')[1]=='.':			out+=[i.split(':')[-1] for i in x.split('=')[1].split(',')]
		return unify_list(out)
	def vcf_rec_group_down(vcf_rec_group,dis_merge=20):
		out=[[],[],[]]#[old_pins,old_pos,new_pin]
		for i in vcf_rec_group:
			sub_dis=[]
			sub_rec=[]
			for j in unify_list(i):
				sub_rec.append(j)
				sub_dis.append(chr_start_end_extract(j)[1])
			sub_dis_clu=[i for i in sub_dis_cluster(sub_dis) if len(i)>1]
			for j in sub_dis_clu:
				sub_rec_clu=[sub_rec[sub_dis.index(k)] for k in j]
				j_new=pick_one_from_exact_overlap(sub_rec_clu)
				out[0].append(sub_rec_clu)
				out[1].append(j)
				out[2].append(j_new)
		return out
	def vcf_rec_reformat_1(pin):
		out=[]
		for x in pin[7].split(';'):
			if x.split('=')[0]=='INFO_POS':
				x_new='INFO_POS='+','.join(unify_list([i.replace('|','/').replace('1/0','0/1') for i in  x.split('=')[1].split(',')]))
				out.append(x_new)
			elif x.split('=')[0]=='INFO_END':
				x_new='INFO_END='+','.join(unify_list([i.replace('|','/').replace('1/0','0/1') for i in  x.split('=')[1].split(',')]))
				out.append(x_new)
			else:	out.append(x)
		pin[7]=';'.join(out)
		return pin
	def vcf_record_QC_2(pin):
		#this function does: 1. unify SV_type names; 2. make sure NUM_CALLER and CALLER matches INFO_POS and INFO_END
		sv_type=pin[4].replace('<','').replace('>','')
		sv_type_new=';'.join(sorted(sv_type.split(';')))
		pin[4]='<'+sv_type_new+'>'
		info_old=pin[7].split(';')
		info_new=[]
		supp_algs=supp_algs_extract(pin)
		for x in info_old:
			if x.split('=')[0]=='MERGE_TYPE':
				info_new.append('MERGE_TYPE='+sv_type)
			elif x.split('=')[0]=='NUM_CALLER':
				info_new.append('NUM_CALLER='+str(len(supp_algs)))
			elif x.split('=')[0]=='CALLER':
				info_new.append('CALLER='+','.join(supp_algs))
			else:
				info_new.append(x)
		pin[7]=';'.join(info_new)
		return pin
	def write_VCF_file(clean_2_hash,file_out,ref):
		write_VCF_header(file_out,ref)
		fo=open(file_out,'a')
		chromos=chromos_readin(ref)
		for k1 in chromos:
			if k1 in clean_2_hash.keys():
				for k2 in sorted(clean_2_hash[k1].keys()):
					for k3 in sorted(clean_2_hash[k1][k2].keys()):
						if len(clean_2_hash[k1][k2][k3])>1:	continue
						else:
							for k4 in clean_2_hash[k1][k2][k3]:	
								k4_new=vcf_record_QC_2(vcf_rec_reformat_1(k4))
								print('\t'.join([str(i) for i in k4_new]),file=fo )
		fo.close()
	def write_VCF_header(output_file,ref):
		fo=open(output_file,'w')
		print( '##fileformat=VCFv4.1',file=fo)
		print('##fileDate='+time.strftime("%Y%m%d"),file=fo)
		ref_info=chromo_contig_prepare(ref)
		for x in ref_info:
			print (x,file=fo)
	   	#print('##reference=hg19'
		print('##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">',file=fo)
		print('##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">',file=fo)
		print('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">',file=fo)
		print('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the variant">',file=fo)
		print('##INFO=<ID=MERGE_TYPE,Number=.,Type=String,Description="Type of structural variant">',file=fo)
		print('##INFO=<ID=NUM_CALLER,Number=1,Type=Integer,Description="=Number of algorithms supporting this variant">',file=fo)
		print('##INFO=<ID=CALLER,Number=.,Type=String,Description="Name of algorithms supporting this variant">',file=fo)
		print('##INFO=<ID=INFO_POS,Number=.,Type=String,Description="Calls by individual algorithm supporting left breakpoint of this variant">',file=fo)
		print('##INFO=<ID=INFO_END,Number=.,Type=String,Description="Calls by individual algorithm supporting right breakpoint of this variant, if different from INFO_POS, otherwise .">',file=fo)
		print('##FILTER=<ID=LowQual,Description="Score of final structural - Theoretical Score <-50">',file=fo)
		print('##ALT=<ID=DEL,Description="Deletion">',file=fo)
		print('##ALT=<ID=DUP,Description="Duplication">',file=fo)
		print('##ALT=<ID=INV,Description="Inversion">',file=fo)
		print('##ALT=<ID=DEL;DUP,Description="merged calls with individual callers report both deletion and duplication">',file=fo)
		print('##ALT=<ID=DEL;INV,Description="merged calls with individual callers report both deletion and inversion">',file=fo)
		print('##ALT=<ID=DUP;INV,Description="merged calls with individual callers report both duplication and inversion">',file=fo)
		print('##ALT=<ID=DEL;DUP;INV,Description="merged calls with individual callers report deletion, duplication and inversion">',file=fo)
		print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',file=fo)
		print('##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype quality">',file=fo)
		print('##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">',file=fo)
		print('##FORMAT=<ID=CNQ,Number=1,Type=Float,Description="Copy number genotype quality for imprecise events">',file=fo)
		print('\t'.join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']),file=fo)
		fo.close()

vcf_prefix_list=['consensus_BP_Calculated_by_CI.new.all_info','consensus_BP_Calculated_by_CI.new.INS']
for vcf_prefix in vcf_prefix_list:
	for s_name in sample_names:
		vcf_name=ppre+vcf_prefix+'.'+s_name+'.vcf'
		final_hash_0=final_check_0_readin_vcf(vcf_name)
		[rec1,clean_1_hash]=final_check_1_vcf(final_hash_0)
		clean_2_hash=final_check_2_vcf(clean_1_hash)
		clean_3_hash=final_check_3_vcf(clean_2_hash)
		write_VCF_file(clean_3_hash,ppre+vcf_prefix+'.'+s_name+'.sorted.vcf',ref)


############################QC3############################
#######check vcf, to see: 1. if CIPOS and CIEND are of similar range #######
#######check vcf, to see: 2. if one SV calls supports multiple merged call #######
#######sorted.vcf to sorted.QC1.vcf#######
############################QC3############################
if Process_Steps[5]=='QCs':
	def info_pos_extact(pin):
		out=[[],[]]
		for x in pin[7].split(';'):
			if x.split('=')[0]=='INFO_POS':
				out[0]+=x.split('=')[1].split(',')
			if x.split('=')[0]=='INFO_END':
				out[1]+=x.split('=')[1].split(',')
				if x.split('=')[1]=='.':	out[1]=out[0]
		return out
	def remove_item_from_INF(k1,pin):
		#eg of k1='23931957:23969515:DEL:0/1:HG00514:wham'
		#eg of pin=['chr22', '23931957', '.', 'A', '<DEL;DUP>', '.', 'PASS', 'MERGE_TYPE=DEL,DUP;END=23969109;SVLEN=37152;CIPOS=-13,9;CIEND=-11,11;NUM_CALLER=6;CALLER=UCSD_Manta,VH,dCGH_filt,delly_illumina,liWGS_illumina,lumpy;INFO_POS=23930885:23970175:DEL:0/1:HG00514:liWGS_illumina,23930885:23970175:DEL:0/1:HG00733:liWGS_illumina,23930885:23970175:DEL:0/1:NA19240:liWGS_illumina,23931245:23969408:DUP:0/1:HG00514:dCGH_filt,23931245:23969408:DUP:1/1:HG00733:dCGH_filt,23931245:23969408:DUP:1/1:NA19240:dCGH_filt,23931955:23969108:DEL:0/1:HG00514:UCSD_Manta,23931955:23969108:DEL:0/1:HG00733:UCSD_Manta,23931955:23969108:DEL:0/1:NA19240:UCSD_Manta,23931955:23969109:DEL:0/1:HG00733:delly_illumina,23931955:23969109:DEL:0/1:NA19240:delly_illumina,23931955:23969109:DEL:1/1:HG00514:delly_illumina,23931957:23969107:DEL:0/1:HG00733:VH,23931957:23969107:DEL:0/1:NA19240:VH,23931957:23969107:DEL:1/1:HG00514:VH,23931957:23969108:DEL:0/1:HG00733:lumpy,23931957:23969108:DEL:0/1:NA19240:lumpy,23931957:23969108:DEL:1/1:HG00514:lumpy,23931957:23969515:DEL:0/1:HG00514:wham;INFO_END=23930885:23970175:DEL:0/1:HG00514:liWGS_illumina,23930885:23970175:DEL:0/1:HG00733:liWGS_illumina,23930885:23970175:DEL:0/1:NA19240:liWGS_illumina,23931245:23969408:DUP:0/1:HG00514:dCGH_filt,23931245:23969408:DUP:1/1:HG00733:dCGH_filt,23931245:23969408:DUP:1/1:NA19240:dCGH_filt,23931541:23969109:DEL:0/1:HG00514:wham,23931955:23969108:DEL:0/1:HG00514:UCSD_Manta,23931955:23969108:DEL:0/1:HG00733:UCSD_Manta,23931955:23969108:DEL:0/1:NA19240:UCSD_Manta,23931955:23969109:DEL:0/1:HG00733:delly_illumina,23931955:23969109:DEL:0/1:NA19240:delly_illumina,23931955:23969109:DEL:1/1:HG00514:delly_illumina,23931957:23969107:DEL:0/1:HG00733:VH,23931957:23969107:DEL:0/1:NA19240:VH,23931957:23969107:DEL:1/1:HG00514:VH,23931957:23969108:DEL:0/1:HG00733:lumpy,23931957:23969108:DEL:0/1:NA19240:lumpy,23931957:23969108:DEL:1/1:HG00514:lumpy']
		info_pos=[]
		info_end=[]
		for i in pin[7].split(';'):
			if i.split('=')[0]=='INFO_POS':
				info_pos=i.split('=')[1].split(',')
			if i.split('=')[0]=='INFO_END':
				info_end=i.split('=')[1].split(',')
		info_pos_new=[i for i in info_pos if not i==k1]
		info_end_new=[i for i in info_end if not i==k1]
		pin_7=pin[7].split(';')
		pin_new_7=[]
		for i in pin_7:
			if i.split('=')[0]=='INFO_POS':	pin_new_7.append('INFO_POS='+','.join(info_pos_new))
			elif i.split('=')[0]=='INFO_END':	pin_new_7.append('INFO_END='+','.join(info_end_new))
			else:	pin_new_7.append(i)
		pin[7]=';'.join(pin_new_7)
		return pin
	def vcf_hash_modify_1(vcf_hash,remove_hash):
		for k1 in remove_hash.keys():
			for k2 in remove_hash[k1]:
				temp=[]
				for v1 in vcf_hash[k2[1]][k2[2]]:
					v1_new=remove_item_from_INF(k1,v1)
					temp.append(v1_new)
				vcf_hash[k2[1]][k2[2]]=temp
		return vcf_hash
	def vcf_rec_modify_2(pin):
		info_list=info_pos_extact(pin)
		sv_type=[]
		for i in info_list:
			if not i==['']:
				for j in i:
					sv_type.append(j.split(':')[2])
		sv_new=[i if not i=='DISDUP' else 'DUP' for i in sv_type]
		sv_type=sorted(unify_list(sv_new))
		return ','.join(sv_type)
	def vcf_hash_modify_2(vcf_hash):
		#correct for MERGE_TYPE
		out={}
		for k1 in vcf_hash.keys():
			if not k1 in out.keys():	out[k1]={}
			for k2 in vcf_hash[k1].keys():
				if not k2 in out[k1].keys():	out[k1][k2]=[]
				for k3 in vcf_hash[k1][k2]:
					sv_type_new=vcf_rec_modify_2(k3)
					if not sv_type_new==[]:
						k3_new=k3[:4]+['<'+sv_type_new.replace(',',';')+'>']+k3[5:7]
						k3_7_old=[i for i in k3[7].split(';') if '=' in i]
						k3_7_new=[]
						for i in k3_7_old:
							if i.split('=')[0]=='MERGE_TYPE':	k3_7_new.append('MERGE_TYPE='+sv_type_new)
							else:								k3_7_new.append(i)
						k3_new.append(';'.join(k3_7_new))
						if not k3_new in out[k1][k2]:	out[k1][k2].append(k3_new)
		return out
	def vcf_more_supp_merge_down(supp_list,vcf_hash,key_info):
		#eg of supp_list=[['chr22', 42705314, 42706355, 'ri'], ['chr22', 42705440, 42706177, 'le'], ['chr22', 42705458, 42706355, 'le'], ['chr22', 42705458, 42706355, 'ri']]
		#eg of key_info='44677206:44680250:DEL:0/1:HG00514:lumpy'
		temp1=[]
		temp_num1=[]
		for k1 in supp_list:
			if not k1[:3] in temp1:	
				temp1.append(k1[:3])
				temp_num1.append([k1[-1]])
			else:
				temp_num1[temp1.index(k1[:3])].append(k1[-1])
		temp_sta1=[len(i) for i in temp_num1]
		temp1_left=[i for i in range(len(temp_sta1)) if temp_sta1[i]==max(temp_sta1)]
		if len(temp1_left)==1:
			return [temp1[temp1_left[0]]+[i] for i in temp_num1[temp1_left[0]]]
		else:
			if max(temp_sta1)==1:
				mean_dis=[abs(i[1]-int(key_info.split(':')[0])) if i[-1]=='ri' else abs(i[2]-int(key_info.split(':')[1])) for i in supp_list]
				return [supp_list[mean_dis.index(min(mean_dis))]]
			elif max(temp_sta1)==2:
				mean_dis=[max([abs(i[1]-int(key_info.split(':')[0])),abs(i[2]-int(key_info.split(':')[1]))])  for i in temp1]
				pos=mean_dis.index(min(mean_dis))
				return [temp1[pos]+[i] for i in temp_num1[pos]]
			else:
				return 'error'
	def polish_test(test):
		out={}
		for k1 in test.keys():
			if len(test[k1])==2 and test[k1][0][:3]==test[k1][1][:3]: continue
			elif len(test[k1])<2:	continue
			else:	out[k1]=test[k1]
		return out
	def vcf_QC_3(vcf_name,chrom):
		fin=open(vcf_name)
		test={}
		vcf_hash={}
		for line in fin:
			pin=line.strip().split()
			if not pin[0][0]=='#' and pin[0]==chrom:
				info_pos=info_pos_extact(pin)
				pos=chr_start_end_extract(pin)[1]
				for x in info_pos[0]: 
					if not x in test.keys():	test[x]=[]
					if not pos+['le'] in test[x]: test[x].append(pos+['le'])
				for x in info_pos[1]: 
					if not x in test.keys():	test[x]=[]
					if not pos+['ri'] in test[x]: test[x].append(pos+['ri'])
				if not pos[1] in vcf_hash.keys():	vcf_hash[pos[1]]={}
				if not pos[2] in vcf_hash[pos[1]].keys():	vcf_hash[pos[1]][pos[2]]=[]
				if not pin in vcf_hash[pos[1]][pos[2]]:	vcf_hash[pos[1]][pos[2]].append(pin)
		fin.close()
		remove_hash={}
		test2=polish_test(test)
		for k1 in test2.keys():
				merge_1=vcf_more_supp_merge_down(test2[k1],vcf_hash,k1)
				if not merge_1=='error':
					remove_1=[i for i in test2[k1] if not i in merge_1]
					if not k1 in remove_hash.keys():	remove_hash[k1]=[]
					remove_hash[k1]+=remove_1
		vcf_hash_new_1=vcf_hash_modify_1(vcf_hash,remove_hash)
		vcf_hash_new_2=vcf_hash_modify_2(vcf_hash_new_1)
		return vcf_hash_new_2
	def CI_info_extract(pin):
		out=[0,0]
		for i in pin[7].split(';'):
			if i.split('=')[0]=='CIPOS':
				out[0]=i.split('=')[1]
			if i.split('=')[0]=='CIEND':
				out[1]=i.split('=')[1]
		return out
	def rep_bp_extract(pin):
		out=[[],[],[],[]]
		for i in pin[7].split(';'):
			if i.split('=')[0]=='INFO_POS':
				out[0]=[j.split(':')[0] for j in  i.split('=')[1].split(',')]
				out[1]=[j.split(':')[-1] for j in  i.split('=')[1].split(',')]
			if i.split('=')[0]=='INFO_END':
				out[2]=[j.split(':')[1] for j in  i.split('=')[1].split(',')]
				out[3]=[j.split(':')[-1] for j in  i.split('=')[1].split(',')]
		out_new=[':'.join([str(i) for i in unify_list(j)]) for j in out]
		return out_new
	def vcf_to_report(vcf_new):
		report_name='.'.join(vcf_new.split('.')[:-1]+['report'])
		fo=open(report_name,'w')
		print('\t'.join(['chr','start_cons','end_cons','ILL_report_le','ILL_report_ri','star-CI;end-CI','reported_bp:left;right','reported_algorithm:left;right']),file=fo)
		fin=open(vcf_new)
		for line in fin:
			pin=line.strip().split()
			if not pin[0][0]=='#':
				rep_new=[]
				pos=chr_start_end_extract(pin)[1]
				info_pos=info_pos_extact(pin)
				rep_bp=';'.join([':'.join(unify_list([i.split(':')[0] for i in info_pos[0]])),':'.join(unify_list([i.split(':')[1] for i in info_pos[1]]))])
				rep_al=';'.join([':'.join(unify_list([i.split(':')[-1] for i in info_pos[0]])),':'.join(unify_list([i.split(':')[-1] for i in info_pos[1]]))])
				print ('\t'.join([str(i) for i in pos+[';'.join(i) for i in info_pos]+[rep_bp,rep_al]]),file=fo)
		fo.close()

prefix='consensus_BP_Calculated_by_CI.new.all_info'
chromos=chromos_readin_canonical(ref)
for s_name in sample_names:
	appdix='sorted.vcf'
	appdis_new='sorted.QC1.vcf'
	vcf_name=ppre+'.'.join([prefix,s_name,appdix])
	vcf_new=ppre+'.'.join([prefix,s_name,appdis_new])
	write_VCF_header(vcf_new,ref)
	fo=open(vcf_new,'a')
	for chrom in chromos[:23]:
		test=vcf_QC_3(vcf_name,chrom)
		for k1 in sorted(test.keys()):
			for k2 in sorted(test[k1].keys()):
				for k3 in sorted(test[k1][k2]):
					if '<DUP;DUP>' in k3:	print(k3)
					print ('\t'.join(k3),file=fo)
	fo.close()
	#vcf_to_report(vcf_new)

if Process_Steps[6]=='merge_INS_with_others':
	for s_name in sample_names:
		os.system(r'''vcf-concat %s %s > %s'''%(ppre+'consensus_BP_Calculated_by_CI.new.all_info.'+s_name+'.sorted.QC1.vcf',ppre+'consensus_BP_Calculated_by_CI.new.INS.'+s_name+'.sorted.vcf',	ppre+'consensus_BP_Calculated_by_CI.new.all_info.'+s_name+'.QC2.vcf'))
		os.system(r'''vcf-sort %s > %s'''%(ppre+'consensus_BP_Calculated_by_CI.new.all_info.'+s_name+'.QC2.vcf', ppre+'consensus_BP_Calculated_by_CI.new.all_info.'+s_name+'.QC2.sorted.vcf'))
		os.system(r'''bgzip %s'''%(ppre+'consensus_BP_Calculated_by_CI.new.all_info.'+s_name+'.QC2.sorted.vcf'))
		os.system(r'''tabix %s'''%(ppre+'consensus_BP_Calculated_by_CI.new.all_info.'+s_name+'.QC2.sorted.vcf.gz'))


if 'size_control_and_blacklist'=='size_control_and_blacklist':
	def pin_cha_extract(pin,chr):
		out=''
		for i in pin[7].split(';'):
			if i.split('=')[0]==chr:
				out=i.split('=')[1]
		return out
	def vcf_readin_size_and_pos_control(file_name):
		fin=os.popen(r'''zcat %s'''%(file_name))
		blaclist_hash=blacklist_readin(blacklist)
		header=[]
		info=[]
		for line in fin:
			pin=line.strip().split()
			if pin[0][:2]=='##': 	header.append(pin)
			else:
				if pin[0][0]=='#': 	info.append(pin)
				else:
					if pin_cha_extract(pin, 'INFO_POS')=='': continue
					svpos=[pin[0],int(pin[1]),int(pin_cha_extract(pin,'END')),int(pin_cha_extract(pin,'SVLEN')), pin_cha_extract(pin,'CALLER')]
					if svpos[2]-svpos[1]<2:	svpos[2]=svpos[1]+svpos[3]
					if svpos[3]==100 and svpos[4] in ['svelter','SVelter']:		
						pin[6]='LowQual'	#remove SVelter only calls that are of size 100 , those are SVelter artifacts
					else:
						if svpos[3]>1000000: 									
							pin[6]='LowQual'		#remove SVs over 1Mb
						else:
							if svpos[0] in blaclist_hash.keys():
								blackref=blaclist_hash[svpos[0]]
								for j in blackref:
									if (j[1]-svpos[1])*(j[2]-svpos[1])>0 and (j[1]-svpos[2])*(j[2]-svpos[2])>0:
										continue
									else:									
										pin[6]='LowQual'	#remove SVs with breakpoints within blacklist
					if not pin[6]=='LowQual':
						pin[6]='PASS'
					info.append(pin)
		fin.close()
		return [header, info]
	def blacklist_readin(blacklist):
		out={}
		fin=open(blacklist)
		for line in fin:
			pin=line.strip().split()
			if not pin[0] in out.keys():	out[pin[0]]=[]
			if len(pin)>2:
				out[pin[0]].append([pin[0],int(pin[1]),int(pin[2])])
		return out
	def vcf_write(file_out, header, info):
		fo=open(file_out,'w')
		for i in header: 
			print(' '.join(i), file=fo)
		for i in info:
			print('\t'.join(i), file=fo)
		fo.close()
	s_name=sample_names[0]
	file_name=ppre+'consensus_BP_Calculated_by_CI.new.all_info.'+s_name+'.QC2.sorted.vcf.gz'
	file_out=ppre+'consensus_BP_Calculated_by_CI.new.all_info.'+s_name+'.QC3.sorted.vcf'
	[header,info]=vcf_readin_size_and_pos_control(file_name)
	vcf_write(file_out, header, info)


if 'redundancy_control'=='redundancy_control':
	def unify_list(list):
		out=[]
		for i in list:
				if not i in out:
						out.append(i)
		return out
	def file_readin(filename):
		fin=open(filename)
		header=[]
		info=[]
		for line in fin:
				pin=line.strip().split()
				if pin[0]=='CHR': 
						header.append(pin)
				else:   
						info.append([pin[0],int(pin[1]),int(pin[2])]+pin[3:])
		fin.close()
		return [header,info]
	def file_write(header, cluster_new, fileout):
			fo=open(fileout,'w')
			for i in header:
					print('\t'.join(i), file=fo)
			for i in cluster_new:
					print('\t'.join([str(j) for j in i]), file=fo)
			fo.close()
	def sv_cha_exractr(pin, character):
			out=''
			for i in pin[7].split(';'):
					if i.split('=')[0]==character:
							out=i.split('=')[1]
			return out
	def vcf_readin(filename):
			fin=open(filename)
			header=[]
			info=[]
			info_new=[]
			for line in fin:
					pin=line.strip().split()
					if pin[0][:2]=='##': 
							header.append(pin)
					else:
							if pin[0][0]=='#':	
									info.append(pin)
							else:
									if pin[6]=='PASS':
											svpos=[pin[0],int(pin[1]),int(sv_cha_exractr(pin,'END')),int(sv_cha_exractr(pin,'SVLEN')),sv_cha_exractr(pin,'NUM_CALLER'),sv_cha_exractr(pin,'MERGE_TYPE'),pin]
											info_new.append(svpos)
									else:
											info.append(pin)
			return [header, info, info_new]
	def write_vcf(header, info, info_new, cluster_new, fileout):
			fo=open(fileout,'w')
			for i in header:		print(' '.join(i), file=fo)
			for i in info:		  print('\t'.join(i), file=fo)
			for i in cluster_new:   print('\t'.join(i[-1]), file=fo)
			for i in info_new:
				if not i in cluster_new:
					i[-1][6]='LowQual'
					print('\t'.join(i[-1]), file=fo)
			fo.close()	
	def main(ppre):
		s_name=sample_names[0]
		filename=ppre+'consensus_BP_Calculated_by_CI.new.all_info.'+s_name+'.QC3.sorted.vcf'
		fileout=ppre+'consensus_BP_Calculated_by_CI.new.all_info.'+s_name+'.QC4.vcf'
		[header, info, info_new]=vcf_readin(filename)
		cluster=[[info_new[0]]]
		max_rec=info_new[0][2]
		chr_rec=info_new[0][0]
		for i in info_new[1:]:
				if i[0]==chr_rec:
						if i[1]>max_rec:	
								cluster.append([i])
								max_rec=i[2]
						else:
								cluster[-1].append(i)
								max_rec=max(max_rec, i[2])
				else:
						chr_rec=i[0]
						max_rec=i[2]
						cluster.append([i])
		cluster_new=[]
		for i in cluster:
				if len(unify_list(i))==1:		cluster_new.append(i[0])
				else:
						tmp_new=[j for j in i if not j[4]=='1']
						if len(tmp_new)==1:	 cluster_new+=tmp_new
						elif len(tmp_new)>1:
								for j in i:	 
										print(j)
								print('')
								caller_num=[int(j[4]) for j in tmp_new]
								tmp_out=tmp_new[caller_num.index(max(caller_num))]
								cluster_new.append(tmp_out)
						elif len(tmp_new)==0:   
								continue
		write_vcf(header, info, info_new, cluster_new, fileout)

main(ppre)


