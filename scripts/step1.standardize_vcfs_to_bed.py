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
	#if '/' in geno:
	#	out=[int(i) for i in geno.split('/')]
	#elif '|' in geno:
	#	out=[int(i) for i in geno.split('|')]
	return geno

def genotype_extract_3(info):
	#eg of info=pin[8:]=['GT:START:END:SVLEN:SPNUM', '0;NA;NA;NA;NA', '1;10530892;10530992;16;3', '1;10530802;10530902;65;2', '0;NA;NA;NA;NA', '0;NA;NA;NA;NA', '0;NA;NA;NA;NA', '0;NA;NA;NA;NA', '0;NA;NA;NA;NA', '1;10530829;10530929;35;3']
	geno_pos=info[0].split(':').index('GT')
	geno_info=[i.split(';')[geno_pos] for i in info[1:]]
	out=[]
	for x in geno_info:
		if x=='0': out.append('1/1')
		elif x=='1': out.append('0/1')
		elif x=='2': out.append('0/0')
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
		if len(i.split(':'))==1 or i=='./.':	continue
		elif len(i.split(':'))>1 and i.split(':')[FT_pos]=='PASS':
			out[rec]=i.split(':')[GT_pos]
	for i in range(len(out)):
		if out[i] in ['./.','.']:	out[i]='0/0' 
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
		print ('\t'.join([str(i) for i in x]), file=fo)
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
			print ([info_all[x],info_all[x+1]])

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

def path_modify(path):
    if not path[-1]=='/':
        path+='/'
    return path

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
			print ('\t'.join([str(i) for i in [k1,k2,bashir_comparison_stat[k1][k2]]]), file=fo)
	fo.close()

def stat_2_list_write(bashir_comparison_stat_2,file_out_name):
	fo=open(file_out_name,'w')
	for k1 in bashir_comparison_stat_2.keys():
		for k2 in bashir_comparison_stat_2[k1].keys():
			for k3 in bashir_comparison_stat_2[k1][k2].keys():
				algorithm_name=num_algorithm_index[k3]
				print('\t'.join([str(i) for i in [k1,k2,k3,bashir_comparison_stat_2[k1][k2][k3],algorithm_name]]), file=fo)
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
		print ('\t'.join([str(i) for i in pacbio_list[k1]+[out_stat[-1]]]), file=fo)
		if len(temp_info)>0: 
			if not len(temp_info) in out_stat_2.keys():
				out_stat_2[len(temp_info)]=[]
			out_stat_2[len(temp_info)]+=temp_info
		#if algorithm_num_count(pacbio_stat_list[k1])==0: print pacbio_stat_list[k1]
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
		#if algorithm_num_count(pacbio_stat_list[k1])==0: print pacbio_stat_list[k1]
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

def unify_list(list):
	out=[]
	for x in list:
		if not x in out:
			out.append(x)
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

def vcf_readin_GenomeSTRiP_cnv(vcf_name,chromo,size_cff,upper_size_cff=249250621):
	#eg of vcf_name='/scratch/remills_flux/xuefzhao/SV_discovery_index/download/different_callers/illumina/20161001_GenomeSTRiP_CNV/HGSVC_trios.illumina.hg38.GenomeSTRiP_cnv.20161001.genotypes.vcf'
	#eg of chromo='chr1'
	fin=open(vcf_name)
	info_hash={}
	test=[]
	for line in fin:
		pin=line.strip().split()
		if not pin[0][:2]=='##': 
			if pin[0][0]=='#':
				sample_names=pin[9:]
				for x in sample_names:
					info_hash[x]={}
			elif pin[0]==chromo and pin[6]=='PASS': 
				[sv_type,sv_pos]=chr_start_end_extract(pin)
				if not sv_type in test:	test.append(sv_type)
				#if sv_type in ['ins','INS','insertion']:	print pin
				if len(sv_pos)==3 and sv_pos[2]-sv_pos[1]>size_cff and sv_pos[2]-sv_pos[1]<upper_size_cff:
					if sv_type=='CNV':  
						sv_main_info=sv_pos+[sv_type]+[copynumber_extract_2([pin[8],pin[i]]) for i in range(9,len(pin))]
						temp_sv_info=[]
						for x in range(len(sv_main_info[4:])):
							if sv_main_info[x+4]>2:
								temp_sv_info.append('dup')
							elif sv_main_info[x+4]<2:
								temp_sv_info.append('del')
							else:
								temp_sv_info.append('NA')
						rec=-1
						for x in sv_main_info[4:]:
							rec+=1
							if not x==2:
								if not sv_main_info[0] in info_hash[sample_names[rec]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]]={}
								if not sv_main_info[1] in info_hash[sample_names[rec]][sv_main_info[0]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]]={}
								if not sv_main_info[2] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]=[]
								if not sv_main_info[:4]+[sv_main_info[rec+4]] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]:
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]].append(sv_main_info[:3]+[temp_sv_info[rec]]+[sv_main_info[rec+4]])
					else: 
						sv_main_info=sv_pos+[sv_type]+[genotype_extract_2([pin[8],pin[i]]) for i in range(9,len(pin))]
						rec=-1
						for x in sv_main_info[4:]:
							rec+=1
							if geno_to_sv(x)=='yes':
								if not sv_main_info[0] in info_hash[sample_names[rec]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]]={}
								if not sv_main_info[1] in info_hash[sample_names[rec]][sv_main_info[0]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]]={}
								if not sv_main_info[2] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]=[]
								if not sv_main_info[:4]+[sv_main_info[rec+4]] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]:
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]].append(sv_main_info[:4]+[sv_main_info[rec+4]])				
	fin.close()
	out_info=sv_hash_order_2(info_hash)
	print('\t'.join(['GenomeStrip']+test))
	return out_info

def vcf_readin_dCGH_filt(vcf_name,chromo,size_cff,upper_size_cff=249250621):
	#eg of vcf_name='/scratch/remills_flux/xuefzhao/SV_discovery_index/download/different_callers/illumina/20160930_dCGH_filt/ALL.dCGH_filt.illumina_high_coverage.20160930.cnv.genotypes.vcf'
	#eg of chromo='chr1'
	fin=open(vcf_name)
	info_hash={}
	test=[]
	for line in fin:
		pin=line.strip().split()
		if not pin[0][:2]=='##': 
			if pin[0][0]=='#':
				sample_names=pin[9:]
				for x in sample_names:
					info_hash[x]={}
			elif pin[0]==chromo and pin[6]=='PASS': 
				[sv_type,sv_pos]=chr_start_end_extract(pin)
				if not sv_type in test:	test.append(sv_type)
				if len(sv_pos)==3 and sv_pos[2]-sv_pos[1]>size_cff and sv_pos[2]-sv_pos[1]<upper_size_cff:
					if sv_type=='CNV': 
						sv_main_info=sv_pos+[sv_type]+[copynumber_extract_2([pin[8],pin[i]]) for i in range(9,len(pin))]
						temp_sv_info=[]
						for x in range(len(sv_main_info[4:])):
							if sv_main_info[x+4]>2:
								temp_sv_info.append('dup')
							elif sv_main_info[x+4]<2:
								temp_sv_info.append('del')
							else:
								temp_sv_info.append('NA')
						rec=-1
						for x in sv_main_info[4:]:
							rec+=1
							if not x==2:
								if not sv_main_info[0] in info_hash[sample_names[rec]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]]={}
								if not sv_main_info[1] in info_hash[sample_names[rec]][sv_main_info[0]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]]={}
								if not sv_main_info[2] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]=[]
								if not sv_main_info[:4]+[sv_main_info[rec+4]] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]:
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]].append(sv_main_info[:3]+[temp_sv_info[rec]]+[sv_main_info[rec+4]])
					else:
						sv_main_info=sv_pos+[sv_type]+[genotype_extract_2([pin[8],pin[i]]) for i in range(9,len(pin))]
						rec=-1
						for x in sv_main_info[4:]:
							rec+=1
							if geno_to_sv(x)=='yes':
								if not sv_main_info[0] in info_hash[sample_names[rec]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]]={}
								if not sv_main_info[1] in info_hash[sample_names[rec]][sv_main_info[0]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]]={}
								if not sv_main_info[2] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]=[]
								if not sv_main_info[:4]+[sv_main_info[rec+4]] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]:
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]].append(sv_main_info[:4]+[sv_main_info[rec+4]])				
	fin.close()
	print('\t'.join(['dCGH']+test))
	return sv_hash_order_2(info_hash)

def vcf_readin_Pindel(vcf_name,chromo,size_cff,upper_size_cff=249250621):
	#eg of vcf_name='/scratch/remills_flux/xuefzhao/SV_discovery_index/download/different_callers/illumina/20160930_Pindel/Pindel.trios.20160930.genotypes.vcf'
	#eg of chromo='chr1'
	fin=open(vcf_name)
	info_hash={}
	test=[]
	for line in fin:
		pin=line.strip().split()
		if not pin[0][:2]=='##': 
			if pin[0][0]=='#':
				sample_names=pin[9:]
				for x in sample_names:
					info_hash[x]={}
			elif pin[0]==chromo and pin[6]=='PASS': 
				[sv_type,sv_pos]=chr_start_end_extract(pin)
				if not sv_type in test:	test.append(sv_type)
				if len(sv_pos)==3 and sv_pos[2]-sv_pos[1]>size_cff and sv_pos[2]-sv_pos[1]<upper_size_cff:
					if sv_type=='CNV': 
						sv_main_info=sv_pos+[sv_type]+[copynumber_extract_2([pin[8],pin[i]]) for i in range(9,len(pin))]
						rec=-1
						for x in sv_main_info[4:]:
							rec+=1
							if not x==2:
								if not sv_main_info[0] in info_hash[sample_names[rec]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]]={}
								if not sv_main_info[1] in info_hash[sample_names[rec]][sv_main_info[0]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]]={}
								if not sv_main_info[2] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]=[]
								if not sv_main_info[:4]+[sv_main_info[rec+4]] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]:
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]].append(sv_main_info[:4]+[sv_main_info[rec+4]])
					else:
						sv_main_info=sv_pos+[sv_type]+[genotype_extract_2([pin[8],pin[i]]) for i in range(9,len(pin))]
						rec=-1
						for x in sv_main_info[4:]:
							rec+=1
							if geno_to_sv(x)=='yes':
								if not sv_main_info[0] in info_hash[sample_names[rec]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]]={}
								if not sv_main_info[1] in info_hash[sample_names[rec]][sv_main_info[0]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]]={}
								if not sv_main_info[2] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]=[]
								if not sv_main_info[:4]+[sv_main_info[rec+4]] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]:
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]].append(sv_main_info[:4]+[sv_main_info[rec+4]])				
	fin.close()
	print('\t'.join(['Pindel']+test))
	return sv_hash_order_2(info_hash)

def vcf_readin_VH(vcf_name,chromo,size_cff,upper_size_cff=249250621):
	#eg of vcf_name='/scratch/remills_flux/xuefzhao/SV_discovery_index/download/different_callers/illumina/20161002_VH/ALL.VH.illumina_high_coverage.DEL.sorted.vcf'
	#eg of chromo='chr1'
	fin=open(vcf_name)
	info_hash={}
	test=[]
	for line in fin:
		pin=line.strip().split()
		if not pin[0][:2]=='##':  
			if pin[0][0]=='#':
				sample_names=pin[9:]
				for x in sample_names:
					info_hash[x]={}
			elif pin[0]==chromo and pin[6]=='PASS':
				[sv_type,sv_pos]=chr_start_end_extract(pin)
				if not sv_type in test:	test.append(sv_type)
				if len(sv_pos)==3 and sv_pos[2]-sv_pos[1]>size_cff and sv_pos[2]-sv_pos[1]<upper_size_cff:
					if sv_type=='CNV':
						sv_main_info=sv_pos+[sv_type]+[copynumber_extract_2([pin[8],pin[i]]) for i in range(9,len(pin))]
						rec=-1
						for x in sv_main_info[4:]:
							rec+=1
							if not x==2:
								if not sv_main_info[0] in info_hash[sample_names[rec]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]]={}
								if not sv_main_info[1] in info_hash[sample_names[rec]][sv_main_info[0]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]]={}
								if not sv_main_info[2] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]=[]
								if not sv_main_info[:4]+[sv_main_info[rec+4]] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]:
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]].append(sv_main_info[:4]+[sv_main_info[rec+4]])
					else:
						sv_main_info=sv_pos+[sv_type]+[genotype_extract_2([pin[8],pin[i]]) for i in range(9,len(pin))]
						rec=-1
						for x in sv_main_info[4:]:
							rec+=1
							if geno_to_sv(x)=='yes':
								if not sv_main_info[0] in info_hash[sample_names[rec]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]]={}
								if not sv_main_info[1] in info_hash[sample_names[rec]][sv_main_info[0]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]]={}
								if not sv_main_info[2] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]=[]
								if not sv_main_info[:4]+[sv_main_info[rec+4]] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]:
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]].append(sv_main_info[:4]+[sv_main_info[rec+4]])				
	fin.close()
	print('\t'.join(['VH']+test))
	return sv_hash_order_2(info_hash)

def vcf_readin_SVelter_UMich_with_Complex(vcf_name,chromo,size_cff,upper_size_cff=249250621):
	#eg of vcf_name='/scratch/remills_flux/xuefzhao/SV_discovery_index/download/different_callers/illumina/20161003_SVelter_UMich/SVelter.all.individual.genotypes.vcf'
	#eg of chromo='chr1'
	fin=open(vcf_name)
	info_hash={}
	test=[]
	for line in fin:
		pin=line.strip().split()
		if not pin[0][:2]=='##': 
			if pin[0][0]=='#':
				sample_names=pin[9:]
				for x in sample_names:
					info_hash[x]={}
			elif pin[0]==chromo and pin[6]=='PASS':	#readin simple event only 
				if pin[4] in ['<CANNOT_CLASSIFY_FOR_NOW>']:continue
				elif pin[4] in ['<TANDUP>', '<DEL>', '<DISDUP>','<DUP_INV>','<INV>']:
					[sv_type,sv_pos]=chr_start_end_extract(pin)
					if not sv_type in test:	test.append(sv_type)
					if len(sv_pos)==3  and sv_pos[2]-sv_pos[1]<upper_size_cff:
						if sv_pos[2]-sv_pos[1]>size_cff:
								sv_main_info=sv_pos+[sv_type]+[genotype_extract_2([pin[8],pin[i]]) for i in range(9,len(pin))]
								rec=-1
								for x in sv_main_info[4:]:
									rec+=1
									if geno_to_sv(x)=='yes':
										if not sv_main_info[0] in info_hash[sample_names[rec]].keys():
											info_hash[sample_names[rec]][sv_main_info[0]]={}
										if not sv_main_info[1] in info_hash[sample_names[rec]][sv_main_info[0]].keys():
											info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]]={}
										if not sv_main_info[2] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]].keys():
											info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]=[]
										if not sv_main_info[:4]+[sv_main_info[rec+4]] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]:
											info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]].append(sv_main_info[:4]+[sv_main_info[rec+4]])				
				else:
					[sv_type_list,sv_pos_list]=svtype_svpos_csv_extract(pin)
					for i in range(len(sv_type_list)):
						sv_type=sv_type_list[i]
						sv_pos=sv_pos_list[i]
						if len(sv_pos)==3  and sv_pos[2]-sv_pos[1]<upper_size_cff:
							if sv_pos[2]-sv_pos[1]>size_cff:
									sv_main_info=sv_pos+[sv_type]+[genotype_extract_2([pin[8],pin[i]]) for i in range(9,len(pin))]
									rec=-1
									for x in sv_main_info[4:]:
										rec+=1
										if geno_to_sv(x)=='yes':
											if not sv_main_info[0] in info_hash[sample_names[rec]].keys():
												info_hash[sample_names[rec]][sv_main_info[0]]={}
											if not sv_main_info[1] in info_hash[sample_names[rec]][sv_main_info[0]].keys():
												info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]]={}
											if not sv_main_info[2] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]].keys():
												info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]=[]
											if not sv_main_info[:4]+[sv_main_info[rec+4]] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]:
												info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]].append(sv_main_info[:4]+[sv_main_info[rec+4]])				
	fin.close()
	print('\t'.join(['SVelter']+test))
	return sv_hash_order_2(info_hash)

def vcf_readin_SVelter_UMich(vcf_name,chromo,size_cff,upper_size_cff=249250621):
	#eg of vcf_name='/scratch/remills_flux/xuefzhao/SV_discovery_index/download/different_callers/illumina/20161003_SVelter_UMich/SVelter.all.individual.genotypes.vcf'
	#eg of chromo='chr1'
	fin=open(vcf_name)
	info_hash={}
	test=[]
	for line in fin:
		pin=line.strip().split()
		if not pin[0][:2]=='##': 
			if pin[0][0]=='#':
				sample_names=pin[9:]
				for x in sample_names:
					info_hash[x]={}
			elif pin[0]==chromo and pin[6]=='PASS':	#readin simple event only 
				if pin[4] in ['<CANNOT_CLASSIFY_FOR_NOW>']:continue
				elif pin[4] in ['<TANDUP>', '<DEL>', '<DISDUP>','<INV>']:
					[sv_type,sv_pos]=chr_start_end_extract(pin)
					if not sv_type in test:	test.append(sv_type)
					if len(sv_pos)==3  and sv_pos[2]-sv_pos[1]<upper_size_cff:
						if sv_pos[2]-sv_pos[1]>size_cff:
							sv_main_info=sv_pos+[sv_type]+[genotype_extract_2([pin[8],pin[i]]) for i in range(9,len(pin))]
							rec=-1
							for x in sv_main_info[4:]:
								rec+=1
								if geno_to_sv_3(x)=='yes':
									if not sv_main_info[0] in info_hash[sample_names[rec]].keys():
										info_hash[sample_names[rec]][sv_main_info[0]]={}
									if not sv_main_info[1] in info_hash[sample_names[rec]][sv_main_info[0]].keys():
										info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]]={}
									if not sv_main_info[2] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]].keys():
										info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]=[]
									if not sv_main_info[:4]+[sv_main_info[rec+4]] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]:
										info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]].append(sv_main_info[:4]+[sv_main_info[rec+4]])				
	fin.close()
	print('\t'.join(['SVelter']+test))
	return sv_hash_order_2(info_hash)

def vcf_readin_melt(vcf_name,chromo,size_cff,upper_size_cff=249250621):
	###only deletions from melt are read in
	#eg of vcf_name='/scratch/remills_flux/xuefzhao/SV_discovery_index/download/different_callers/illumina/20160930_MELT/ALL.melt.MEI.illumina_high_coverage.20160930.genotypes.vcf'
	fin=open(vcf_name)
	info_hash={}
	test=[]
	for line in fin:
		pin=line.strip().split()
		if not pin[0][:2]=='##': 
			if pin[0][0]=='#':
				sample_names=pin[9:]
				for x in sample_names:
					info_hash[x]={}
			elif pin[0]==chromo and pin[6] in ['PASS','.']:	#readin simple event only 
				if 'CN' in pin[4]:
					sv_type='DEL'
					sv_pos=chr_start_end_extract_melt(pin)
					if not sv_type in test:	test.append(sv_type)
					if len(sv_pos)==3 and sv_pos[2]-sv_pos[1]>size_cff and sv_pos[2]-sv_pos[1]<upper_size_cff and not '0' in sv_pos:
						sv_main_info=sv_pos+[sv_type]+[genotype_extract_2([pin[8],pin[i]]) for i in range(-len(sample_names),0)]
						rec=-1
						for x in sv_main_info[4:]:
							rec+=1
							if geno_to_sv(x)=='yes':
								if not sv_main_info[0] in info_hash[sample_names[rec]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]]={}
								if not sv_main_info[1] in info_hash[sample_names[rec]][sv_main_info[0]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]]={}
								if not sv_main_info[2] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]=[]
								if not sv_main_info[:4]+[sv_main_info[rec+4]] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]:
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]].append(sv_main_info[:4]+[sv_main_info[rec+4]])				
				else:
					[sv_type,sv_pos]=chr_start_end_extract(pin)
					sv_type=pin[4].replace('<','').replace('>','')
					if not sv_type in test:	test.append(sv_type)
					if len(sv_pos)==3  and sv_pos[2]-sv_pos[1]>size_cff and sv_pos[2]-sv_pos[1]>size_cff:
						sv_main_info=sv_pos+[sv_type]+[genotype_extract_2([pin[8],pin[i]]) for i in range(-len(sample_names),0)]
						rec=-1
						for x in sv_main_info[4:]:
							rec+=1
							if geno_to_sv(x)=='yes':
								if not sv_main_info[0] in info_hash[sample_names[rec]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]]={}
								if not sv_main_info[1] in info_hash[sample_names[rec]][sv_main_info[0]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]]={}
								if not sv_main_info[2] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]=[]
								if not sv_main_info[:4]+[sv_main_info[rec+4]] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]:
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]].append(sv_main_info[:4]+[sv_main_info[rec+4]])				
	fin.close()
	print('\t'.join(['MELT']+test))
	return sv_hash_order_2(info_hash)

def tardis_INS_info_extract(pin):
	out=[pin[0],int(pin[1])]
	sv_type=pin[4].replace('<','').replace('>','')
	sample_names=[i for i in ['HG00512', 'HG00513', 'HG00514', 'HG00731', 'HG00732', 'HG00733', 'NA19238', 'NA19239', 'NA19240'] if i in pin[7]]
	return [sv_type,out,sample_names]

def vcf_readin_tardis(vcf_name,chromo,size_cff,upper_size_cff=249250621):
	fin=open(vcf_name)
	info_hash={}
	for line in fin:
		pin=line.strip().split()
		if not pin[0][:2]=='##': 
			if pin[0]==chromo and pin[6] in ['PASS','.']:
				[sv_type,INS_pos,sample_names]=tardis_INS_info_extract(pin)
				if not INS_pos[0] in info_hash.keys():	info_hash[INS_pos[0]]={}
				if not int(INS_pos[1]) in info_hash[INS_pos[0]].keys():info_hash[INS_pos[0]][int(INS_pos[1])]=[]
				info_hash[INS_pos[0]][int(INS_pos[1])].append([sv_type]+sample_names)
	fin.close()
	out={}
	for k1 in sorted(info_hash.keys()):
		for k2 in sorted(info_hash[k1].keys()):
			for k3 in info_hash[k1][k2]:
				for k4 in k3[1:]:
					if not k4 in out.keys():	out[k4]=[]
					if not [k1,k2,k3[0]] in out[k4]: out[k4].append([k1,k2,k2,k3[0],'./.'])
	return out

def vcf_readin_cloudSV(vcf_name,chromo,size_cff,upper_size_cff=249250621):
	#eg of vcf_name='/scratch/remills_flux/xuefzhao/SV_discovery_index/download/different_callers/illumina/20161003_cloudSV/cloudSV_GRCh38_all_pb_ill.vcf'
	#eg of chromo='chr1'
	fin=open(vcf_name)
	info_hash={}
	for line in fin:
		pin=line.strip().split()
		if not pin[0][:2]=='##': 
			if pin[0][0]=='#':
				sample_names=pin[9:]
				for x in sample_names:
					info_hash[x]={}
			elif pin[0]==chromo:
				[sv_type,sv_pos]=chr_start_end_extract(pin)
				if len(sv_pos)==3 and sv_pos[2]-sv_pos[1]>size_cff and sv_pos[2]-sv_pos[1]<upper_size_cff:
					if sv_type=='CNV':
						sv_main_info=sv_pos+[sv_type]+[copynumber_extract_2([pin[8],pin[i]]) for i in range(9,len(pin))]
						temp_sv_info=[]
						for x in range(len(sv_main_info[4:])):
							if sv_main_info[x+4]>2:
								temp_sv_info.append('dup')
							elif sv_main_info[x+4]<2:
								temp_sv_info.append('del')
						rec=-1
						for x in sv_main_info[4:]:
							rec+=1
							if not x==2:
								if not sv_main_info[0] in info_hash[sample_names[rec]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]]={}
								if not sv_main_info[1] in info_hash[sample_names[rec]][sv_main_info[0]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]]={}
								if not sv_main_info[2] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]=[]
								if not sv_main_info[:4]+[sv_main_info[rec+4]] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]:
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]].append(sv_main_info[:3]+[temp_sv_info[rec]]+[sv_main_info[rec+4]])
					else:
						sv_main_info=sv_pos+[sv_type]+[genotype_extract_2([pin[8],pin[i]]) for i in range(9,len(pin))]
						rec=-1
						for x in sv_main_info[4:]:
							rec+=1
							if geno_to_sv(x)=='yes':
								if not sv_main_info[0] in info_hash[sample_names[rec]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]]={}
								if not sv_main_info[1] in info_hash[sample_names[rec]][sv_main_info[0]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]]={}
								if not sv_main_info[2] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]=[]
								if not sv_main_info[:4]+[sv_main_info[rec+4]] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]:
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]].append(sv_main_info[:4]+[sv_main_info[rec+4]])				
	fin.close()
	return sv_hash_order_2(info_hash)

def retroCNV_INS_info_extract(pin):
	out=[pin[0],int(pin[1])]
	INS_region=[0,0,0]
	for i in pin[7].split(';'):
		if i.split('=')[0]=='CHRPAR':	INS_region[0]=i.split('=')[1]
		if i.split('=')[0]=='STARTPAR':	INS_region[1]=int(i.split('=')[1])
		if i.split('=')[0]=='ENDPAR':	INS_region[2]=int(i.split('=')[1])
	sv_type=pin[4].replace('<','').replace('>','')
	return [sv_type,out,INS_region]

def vcf_delly_ins_len_readin(pin):
	out=0
	for k1 in pin[7].split(';'):
		if k1.split('=')[0]=='INSLEN':	out=int(k1.split('=')[1])
	return out

def vcf_readin_retroCNV(vcf_name,chromo,size_cff,upper_size_cff=249250621):
	#eg of vcf_name='/scratch/remills_flux/xuefzhao/SV_discovery_index/download/different_callers/illumina/20161003_retroCNV/20161003_retroCNV.vcf'
	#eg of chromo='chr1'
	fin=open(vcf_name)
	info_hash={}
	for line in fin:
		pin=line.strip().split()
		if not pin[0][:2]=='##': 
			if pin[0][0]=='#':
				sample_names=pin[9:]
				for x in sample_names:
					info_hash[x]={}
			elif pin[0]==chromo and pin[6]=='PASS': 
				[sv_type,INS_pos,INS_length]=retroCNV_INS_info_extract(pin)
				if not INS_pos[0] in info_hash.keys():	info_hash[INS_pos[0]]={}
				if not int(INS_pos[1]) in info_hash[INS_pos[0]].keys():	info_hash[INS_pos[0]][int(INS_pos[1])]=[sv_type,'_'.join([str(i) for i in INS_length])]
	fin.close()
	out={}
	for k0 in ['HG00512','HG00513','HG00514','HG00731','HG00732','HG00733','NA19238','NA19239','NA19240']:
		out[k0]=[]
		for k1 in info_hash.keys():
			for k2 in sorted(info_hash[k1].keys()):
				out[k0].append([k1,k2,k2]+info_hash[k1][k2])
	return out

def vcf_readin_wham(vcf_name,chromo,size_cff,upper_size_cff=249250621):
	#eg of vcf_name='/scratch/remills_flux/xuefzhao/SV_discovery_index/download/different_callers/illumina/20160930_wham/ALL.wham.20160930.genotypes.vcf'
	#eg of chromo='chr1'
	fin=open(vcf_name)
	info_hash={}
	test=[]
	for line in fin:
		pin=line.strip().split()
		if not pin[0][:2]=='##': 
			if pin[0][0]=='#':
				sample_names=pin[9:]
				for x in sample_names:
					info_hash[x]={}
			elif pin[0]==chromo:
				[sv_type,sv_pos]=chr_start_end_extract(pin)
				if not sv_type in test:	test.append(sv_type)
				if len(sv_pos)==3 and sv_pos[2]-sv_pos[1]>size_cff and sv_pos[2]-sv_pos[1]<upper_size_cff:
					if sv_type=='CNV':
						sv_main_info=sv_pos+[sv_type]+[copynumber_extract_2([pin[8],pin[i]]) for i in range(9,len(pin))]
						temp_sv_info=[]
						for x in range(len(sv_main_info[4:])):
							if sv_main_info[x+4]>2:
								temp_sv_info.append('dup')
							elif sv_main_info[x+4]<2:
								temp_sv_info.append('del')
						rec=-1
						for x in sv_main_info[4:]:
							rec+=1
							if not x==2:
								if not sv_main_info[0] in info_hash[sample_names[rec]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]]={}
								if not sv_main_info[1] in info_hash[sample_names[rec]][sv_main_info[0]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]]={}
								if not sv_main_info[2] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]=[]
								if not sv_main_info[:4]+[sv_main_info[rec+4]] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]:
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]].append(sv_main_info[:3]+[temp_sv_info[rec]]+[sv_main_info[rec+4]])
					else:
						sv_main_info=sv_pos+[sv_type]+[genotype_extract_2([pin[8],pin[i]]) for i in range(9,len(pin))]
						rec=-1
						for x in sv_main_info[4:]:
							rec+=1
							if geno_to_sv(x)=='yes':
								if not sv_main_info[0] in info_hash[sample_names[rec]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]]={}
								if not sv_main_info[1] in info_hash[sample_names[rec]][sv_main_info[0]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]]={}
								if not sv_main_info[2] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]=[]
								if not sv_main_info[:4]+[sv_main_info[rec+4]] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]:
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]].append(sv_main_info[:4]+[sv_main_info[rec+4]])				
	fin.close()
	print('\t'.join(['Wham']+test))
	return sv_hash_order_2(info_hash)

def vcf_readin_wham_ins(vcf_name,chromo,size_cff,upper_size_cff=249250621):
	#eg of vcf_name='/scratch/remills_flux/xuefzhao/SV_discovery_index/download/different_callers/illumina/20160930_wham/ALL.wham.insertions.20160930.sites.vcf'
	#eg of chromo='chr1'
	fin=open(vcf_name)
	info_hash={}
	test=[]
	for line in fin:
		pin=line.strip().split()
		if not pin[0][:2]=='##': 
			if pin[0][0]=='#':
				sample_names=pin[9:]
				for x in sample_names:
					info_hash[x]={}
			elif pin[0]==chromo:
				[sv_type,sv_pos]=chr_start_end_extract(pin)
				if not sv_type in test:	test.append(sv_type)
				if len(sv_pos)==3 and sv_pos[2]-sv_pos[1]>size_cff and sv_pos[2]-sv_pos[1]<upper_size_cff:
					if sv_type=='CNV':
						sv_main_info=sv_pos+[sv_type]+[copynumber_extract_2([pin[8],pin[i]]) for i in range(9,len(pin))]
						temp_sv_info=[]
						for x in range(len(sv_main_info[4:])):
							if sv_main_info[x+4]>2:
								temp_sv_info.append('dup')
							elif sv_main_info[x+4]<2:
								temp_sv_info.append('del')
						rec=-1
						for x in sv_main_info[4:]:
							rec+=1
							if not x==2:
								if not sv_main_info[0] in info_hash[sample_names[rec]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]]={}
								if not sv_main_info[1] in info_hash[sample_names[rec]][sv_main_info[0]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]]={}
								if not sv_main_info[2] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]=[]
								if not sv_main_info[:4]+[sv_main_info[rec+4]] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]:
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]].append(sv_main_info[:3]+[temp_sv_info[rec]]+[sv_main_info[rec+4]])
					else:
						sv_main_info=sv_pos+[sv_type]+[genotype_extract_2([pin[8],pin[i]]) for i in range(9,len(pin))]
						rec=-1
						for x in sv_main_info[4:]:
							rec+=1
							if geno_to_sv(x)=='yes':
								if not sv_main_info[0] in info_hash[sample_names[rec]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]]={}
								if not sv_main_info[1] in info_hash[sample_names[rec]][sv_main_info[0]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]]={}
								if not sv_main_info[2] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]=[]
								if not sv_main_info[:4]+[sv_main_info[rec+4]] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]:
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]].append(sv_main_info[:4]+[sv_main_info[rec+4]])				
	fin.close()
	print('\t'.join(['Wham']+test))
	return sv_hash_order_2(info_hash)

def vcf_readin_lumpy(vcf_name,chromo,size_cff,upper_size_cff=249250621):
	#eg of vcf_name='/scratch/remills_flux/xuefzhao/SV_discovery_index/download/different_callers/illumina/20160930_lumpy/ALL.lumpy.20160930.genotypes.vcf'
	#eg of chromo='chr1'
	fin=open(vcf_name)
	info_hash={}
	test=[]
	for line in fin:
		pin=line.strip().split()
		if not pin[0][:2]=='##': 
			if pin[0][0]=='#':
				sample_names=pin[9:]
				for x in sample_names:
					info_hash[x]={}
			elif pin[0]==chromo:
				[sv_type,sv_pos]=chr_start_end_extract(pin)
				if not sv_type in test:	test.append(sv_type)
				if len(sv_pos)==3 and sv_pos[2]-sv_pos[1]>size_cff and sv_pos[2]-sv_pos[1]<upper_size_cff:
					if sv_type=='CNV':
						sv_main_info=sv_pos+[sv_type]+[copynumber_extract_2([pin[8],pin[i]]) for i in range(9,len(pin))]
						temp_sv_info=[]
						for x in range(len(sv_main_info[4:])):
							if sv_main_info[x+4]>2:
								temp_sv_info.append('dup')
							elif sv_main_info[x+4]<2:
								temp_sv_info.append('del')
						rec=-1
						for x in sv_main_info[4:]:
							rec+=1
							if not x==2:
								if not sv_main_info[0] in info_hash[sample_names[rec]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]]={}
								if not sv_main_info[1] in info_hash[sample_names[rec]][sv_main_info[0]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]]={}
								if not sv_main_info[2] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]=[]
								if not sv_main_info[:4]+[sv_main_info[rec+4]] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]:
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]].append(sv_main_info[:3]+[temp_sv_info[rec]]+[sv_main_info[rec+4]])
					else:
						sv_main_info=sv_pos+[sv_type]+[genotype_extract_2([pin[8],pin[i]]) for i in range(9,len(pin))]
						rec=-1
						for x in sv_main_info[4:]:
							rec+=1
							if geno_to_sv(x)=='yes':
								if not sv_main_info[0] in info_hash[sample_names[rec]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]]={}
								if not sv_main_info[1] in info_hash[sample_names[rec]][sv_main_info[0]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]]={}
								if not sv_main_info[2] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]].keys():
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]=[]
								if not sv_main_info[:4]+[sv_main_info[rec+4]] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]:
									info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]].append(sv_main_info[:4]+[sv_main_info[rec+4]])				
	fin.close()
	print('\t'.join(['Lumpy']+test))
	return sv_hash_order_2(info_hash)

def vcf_readin_10x(vcf_hash,chromo,size_cff,upper_size_cff=249250621):
	#eg of vcf_hash=file_hash['10x_genomics_svs']
	info_hash={}
	test=[]
	for k1 in vcf_hash.keys():
		info_hash[k1]={}
		for k2 in vcf_hash[k1]:
			fin=open(k2)
			for line in fin:
				pin=line.strip().split()
				if not pin[0][:2]=='##': 
					if pin[0][0]=='#':	continue	
					elif pin[0]==chromo and 'PASS' in pin[6]:  
						genotype_info=genotype_extract_2(pin[8:10])
						[sv_type,sv_pos]=chr_start_end_extract(pin)
						if not sv_type in test:	test.append(sv_type)
						if geno_to_sv(genotype_info)=='yes' and len(sv_pos)==3 and sv_pos[2]-sv_pos[1]>size_cff and sv_pos[2]-sv_pos[1]<upper_size_cff:
							sv_main_info=sv_pos+[sv_type]+[genotype_info]
							if not sv_main_info[0] in info_hash[k1].keys():
								info_hash[k1][sv_main_info[0]]={}
							if not sv_main_info[1] in info_hash[k1][sv_main_info[0]].keys():
								info_hash[k1][sv_main_info[0]][sv_main_info[1]]={}
							if not sv_main_info[2] in info_hash[k1][sv_main_info[0]][sv_main_info[1]].keys():
								info_hash[k1][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]=[]
							if not sv_main_info in info_hash[k1][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]:
								info_hash[k1][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]].append(sv_main_info)				
			fin.close()
	print ('\t'.join(['10x']+test))
	return sv_hash_order_2(info_hash)

def vcf_readin_kchen_novobreak(vcf_hash,chromo,size_cff,upper_size_cff=249250621):
	#eg of vcf_hash=file_hash['kchen_novobreak']
	info_hash={}
	test=[]
	for k1 in vcf_hash.keys():
		info_hash[k1.split('.')[0]]={}
		for k2 in vcf_hash[k1]:
			fin=open(k2)
			for line in fin:
				pin=line.strip().split()
				if not pin[0][:2]=='##': 
					if pin[0][0]=='#':	continue	
					elif pin[0]==chromo and 'PASS' in pin[6]:
						genotype_info=genotype_extract_2(pin[8:10])
						[sv_type,sv_pos]=chr_start_end_extract(pin)
						if not sv_type in test:	test.append(sv_type)
						if geno_to_sv(genotype_info)=='yes' and len(sv_pos)==3 and sv_pos[2]-sv_pos[1]>size_cff and sv_pos[2]-sv_pos[1]<upper_size_cff:
							sv_main_info=sv_pos+[sv_type]+[genotype_info]
							if not sv_main_info[0] in info_hash[k1.split('.')[0]].keys():
								info_hash[k1.split('.')[0]][sv_main_info[0]]={}
							if not sv_main_info[1] in info_hash[k1.split('.')[0]][sv_main_info[0]].keys():
								info_hash[k1.split('.')[0]][sv_main_info[0]][sv_main_info[1]]={}
							if not sv_main_info[2] in info_hash[k1.split('.')[0]][sv_main_info[0]][sv_main_info[1]].keys():
								info_hash[k1.split('.')[0]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]=[]
							if not sv_main_info in info_hash[k1.split('.')[0]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]:
								info_hash[k1.split('.')[0]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]].append(sv_main_info)				
			fin.close()
	print('\t'.join(['novoBreak']+test))
	return sv_hash_order_2(info_hash)

def vcf_readin_kchen_HySA(vcf_name,chromo,size_cff,upper_size_cff=249250621):
	#eg of vcf_name=['/scratch/remills_flux/xuefzhao/SV_discovery_index/download/different_callers/illumina/20170126_kchen_calls_hysa_regenotype/ALL.wgs.KchenLab-HySA.01262017.DEL.Illumina_Pacbio.deep_coverage.sites.vcf', '/scratch/remills_flux/xuefzhao/SV_discovery_index/download/different_callers/illumina/20170126_kchen_calls_hysa_regenotype/ALL.wgs.KchenLab-HySA.01262017.INS.Illumina_Pacbio.deep_coverage.sites.vcf']
	info_hash={}
	for k2 in vcf_name:
		fin=open(k2)
		for line in fin:
			pin=line.strip().split()
			if not pin[0][:2]=='##': 
				if pin[0][0]=='#':	
					sample_names=pin[9:]
					for x in sample_names:
						if not x in info_hash.keys():	info_hash[x]={}
				elif pin[0]==chromo and 'PASS' in pin[6]: 
					pos_info=[pin[9+i].split(';')+[sample_names[i]] for i in range(len(sample_names))]
					sv_type=svtype_extract(pin)
					for i_2 in pos_info:
						i=[]
						for j in i_2:	i+=j.split(':')
						if 'NA' in i:	continue
						else:
							sv_pos=[pin[0],int(i[1]),int(i[1])+int(i[3])]
							if sv_pos[2]-sv_pos[1]>size_cff and sv_pos[2]-sv_pos[1]<upper_size_cff:
								if not i[0]=='0/0':
									if not sv_pos[0] in info_hash[i[-1]].keys():	info_hash[i[-1]][sv_pos[0]]={}
									if not sv_pos[1] in info_hash[i[-1]][sv_pos[0]].keys():	info_hash[i[-1]][sv_pos[0]][sv_pos[1]]={}
									if not sv_pos[2] in info_hash[i[-1]][sv_pos[0]][sv_pos[1]].keys():	info_hash[i[-1]][sv_pos[0]][sv_pos[1]][sv_pos[2]]=[]
									if not sv_pos+[sv_type,i[0]] in info_hash[i[-1]][sv_pos[0]][sv_pos[1]][sv_pos[2]]: info_hash[i[-1]][sv_pos[0]][sv_pos[1]][sv_pos[2]].append(sv_pos+[sv_type,i[0]])
		fin.close()
	return sv_hash_order_2(info_hash)

def vcf_readin_delly(vcf_list,chromo,size_cff,upper_size_cff=249250621):
	#vcf_list=file_hash['delly_illumina']
	#eg of chromo='chr1'
	info_hash={}
	test=[]
	for vcf_name in vcf_list:
		fin=open(vcf_name)
		for line in fin:
			pin=line.strip().split()
			if not pin[0][:2]=='##': 
				if pin[0][0]=='#':
					sample_names=pin[9:]
					for x in sample_names:
						if not x in info_hash.keys():
							info_hash[x]={}
				elif pin[0]==chromo and pin[6]=='PASS':
					[sv_type,sv_pos]=chr_start_end_extract(pin)
					if not sv_type in test:	test.append(sv_type)
					if sv_type=='INS':	
						sv_len=vcf_delly_ins_len_readin(pin)
						sv_pos[2]=sv_pos[1]+sv_len
					if len(sv_pos)==3:
						if sv_pos[2]-sv_pos[1]>size_cff and sv_pos[2]-sv_pos[1]<upper_size_cff:
							if sv_type=='CNV':
								sv_main_info=sv_pos+[sv_type]+[copynumber_extract_2([pin[8],pin[i]]) for i in range(9,len(pin))]
								temp_sv_info=[]
								for x in range(len(sv_main_info[4:])):
									if sv_main_info[x+4]>2:
										temp_sv_info.append('dup')
									elif sv_main_info[x+4]<2:
										temp_sv_info.append('del')
								rec=-1
								for x in sv_main_info[4:]:
									rec+=1
									if not x==2:
										if not sv_main_info[0] in info_hash[sample_names[rec]].keys():
											info_hash[sample_names[rec]][sv_main_info[0]]={}
										if not sv_main_info[1] in info_hash[sample_names[rec]][sv_main_info[0]].keys():
											info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]]={}
										if not sv_main_info[2] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]].keys():
											info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]=[]
										if not sv_main_info[:4]+[sv_main_info[rec+4]] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]:
											info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]].append(sv_main_info[:3]+[temp_sv_info[rec]]+[sv_main_info[rec+4]])
							else:
								sv_main_info=sv_pos+[sv_type]+[genotype_extract_2([pin[8],pin[i]]) for i in range(9,len(pin))]
								rec=-1
								for x in sv_main_info[4:]:
									rec+=1
									if geno_to_sv(x)=='yes':
										if not sv_main_info[0] in info_hash[sample_names[rec]].keys():
											info_hash[sample_names[rec]][sv_main_info[0]]={}
										if not sv_main_info[1] in info_hash[sample_names[rec]][sv_main_info[0]].keys():
											info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]]={}
										if not sv_main_info[2] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]].keys():
											info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]=[]
										if not sv_main_info[:4]+[sv_main_info[rec+4]] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]:
											info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]].append(sv_main_info[:4]+[sv_main_info[rec+4]])				
		fin.close()
	print('\t'.join(['delly']+test)) 
	return sv_hash_order_2(info_hash)

def vcf_readin_svs_PCR_free(vcf_list,chromo,size_cff,upper_size_cff=249250621):
	#vcf_list=file_hash['sebat_forest_SV']
	#eg of chromo='chr1'
	info_hash={}
	for vcf_name in vcf_list:
		fin=open(vcf_name)
		for line in fin:
			pin=line.strip().split()
			if not pin[0][:2]=='##': 
				if pin[0][0]=='#':
					sample_names=pin[9:]
					for x in sample_names:
						if not x in info_hash.keys():
							info_hash[x]={}
				elif pin[0]==chromo and pin[6]=='PASS': 
					[sv_type,sv_pos]=chr_start_end_extract(pin)
					if len(sv_pos)==3:
						if sv_pos[2]-sv_pos[1]>size_cff and sv_pos[2]-sv_pos[1]<upper_size_cff:
							if sv_type=='CNV':
								sv_main_info=sv_pos+[sv_type]+[copynumber_extract_2([pin[8],pin[i]]) for i in range(9,len(pin))]
								temp_sv_info=[]
								for x in range(len(sv_main_info[4:])):
									if sv_main_info[x+4]>2:
										temp_sv_info.append('dup')
									elif sv_main_info[x+4]<2:
										temp_sv_info.append('del')
								rec=-1
								for x in sv_main_info[4:]:
									rec+=1
									if not x==2:
										if not sv_main_info[0] in info_hash[sample_names[rec]].keys():
											info_hash[sample_names[rec]][sv_main_info[0]]={}
										if not sv_main_info[1] in info_hash[sample_names[rec]][sv_main_info[0]].keys():
											info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]]={}
										if not sv_main_info[2] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]].keys():
											info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]=[]
										if not sv_main_info[:4]+[sv_main_info[rec+4]] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]:
											info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]].append(sv_main_info[:3]+[temp_sv_info[rec]]+[sv_main_info[rec+4]])
							else:
								sv_main_info=sv_pos+[sv_type]+[genotype_extract_2([pin[8],pin[i]]) for i in range(9,len(pin))]
								rec=-1
								for x in sv_main_info[4:]:
									rec+=1
									if geno_to_sv(x)=='yes':
										if not sv_main_info[0] in info_hash[sample_names[rec]].keys():
											info_hash[sample_names[rec]][sv_main_info[0]]={}
										if not sv_main_info[1] in info_hash[sample_names[rec]][sv_main_info[0]].keys():
											info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]]={}
										if not sv_main_info[2] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]].keys():
											info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]=[]
										if not sv_main_info[:4]+[sv_main_info[rec+4]] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]:
											info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]].append(sv_main_info[:4]+[sv_main_info[rec+4]])				
		fin.close()
	return sv_hash_order_2(info_hash)

def vcf_readin_Manta(vcf_name,chromo,size_cff,upper_size_cff=249250621):
	#vcf_list=['/scratch/remills_flux/xuefzhao/SV_discovery_index/download/different_callers/illumina/20161108_UCSD_Manta/ALL.wgs.UCSD_Manta.20162710.sv.Illumina_high-coverage_PCR-free.genotypes.vcf']
	#eg of chromo='chr1'
	info_hash={}
	test=[]
	fin=open(vcf_name)
	for line in fin:
		pin=line.strip().split()
		if not pin[0][:2]=='##': 
			if pin[0][0]=='#':
				sample_names=pin[9:]
				for x in sample_names:
					if not x in info_hash.keys():
						info_hash[x]={}
			elif pin[0]==chromo and pin[6]=='PASS': 
				[sv_type,sv_pos]=chr_start_end_extract(pin)
				if not sv_type in test:	test.append(sv_type)
				if len(sv_pos)==3:
					if sv_pos[2]-sv_pos[1]>size_cff and sv_pos[2]-sv_pos[1]<upper_size_cff:
							sv_main_info=sv_pos+[sv_type]+genotype_extract_Manta(pin)
							rec=-1
							for x in sv_main_info[4:]:
								rec+=1
								if geno_to_sv(x)=='yes':
									if not sv_main_info[0] in info_hash[sample_names[rec]].keys():
										info_hash[sample_names[rec]][sv_main_info[0]]={}
									if not sv_main_info[1] in info_hash[sample_names[rec]][sv_main_info[0]].keys():
										info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]]={}
									if not sv_main_info[2] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]].keys():
										info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]=[]
									if not sv_main_info[:4]+[sv_main_info[rec+4]] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]:
										info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]].append(sv_main_info[:4]+[sv_main_info[rec+4]])				
	fin.close()
	print('\t'.join(['Manta']+test))
	return sv_hash_order_2(info_hash)

def vcf_readin_ForestSV(vcf_name,chromo,size_cff,upper_size_cff=249250621):
	#eg of vcf_list=['/scratch/remills_flux/xuefzhao/SV_discovery_index/download/different_callers/illumina/20161108_UCSD_ForestSV/ALL.wgs.UCSD_ForestSV-gtCNV.20163110.sv.Illumina_high-coverage_PCR-free.genotypes.vcf']
	#eg of chromo='chr1'
	test=[]
	info_hash={}
	fin=open(vcf_name)
	for line in fin:
		pin=line.strip().split()
		if not pin[0][:2]=='##': 
			if pin[0][0]=='#':
				sample_names=pin[9:]
				for x in sample_names:
					if not x in info_hash.keys():
						info_hash[x]={}
			elif pin[0]==chromo and pin[6]=='PASS': 
				[sv_type,sv_pos]=chr_start_end_extract(pin)
				if not sv_type in test:	test.append(sv_type)
				if len(sv_pos)==3:
					if sv_pos[2]-sv_pos[1]>size_cff and sv_pos[2]-sv_pos[1]<upper_size_cff:
							sv_main_info=sv_pos+[sv_type]+genotype_extract_Manta(pin)
							rec=-1
							for x in sv_main_info[4:]:
								rec+=1
								if geno_to_sv(x)=='yes':
									if not sv_main_info[0] in info_hash[sample_names[rec]].keys():
										info_hash[sample_names[rec]][sv_main_info[0]]={}
									if not sv_main_info[1] in info_hash[sample_names[rec]][sv_main_info[0]].keys():
										info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]]={}
									if not sv_main_info[2] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]].keys():
										info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]=[]
									if not sv_main_info[:4]+[sv_main_info[rec+4]] in info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]:
										info_hash[sample_names[rec]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]].append(sv_main_info[:4]+[sv_main_info[rec+4]])				
	fin.close()
	print('\t'.join(['ForestSV']+test))
	return sv_hash_order_2(info_hash)

def vcf_readin_Moleculo(Moleculo_ucsd,chromo,size_cff,upper_size_cff=249250621):
	fin=open(Moleculo_ucsd)
	info_hash={}
	test=[]
	for line in fin:
		pin=line.strip().split()
		if not pin[0][:2]=='##': 
			if pin[0][0]=='#':	
				sample_names=pin[9:]
				for x in sample_names:
					if not x in info_hash.keys():	info_hash[x]={}
			else: 
				if pin[0]==chromo:
					[sv_type,sv_pos]=chr_start_end_extract(pin)
					if not sv_type in test:test.append(sv_type)
					if sv_pos[2]-sv_pos[1]>size_cff and sv_pos[2]-sv_pos[1]<upper_size_cff:
						sv_main_info=sv_pos+[sv_type]+genotype_extract_3(pin)
						for i in range(len(sample_names)):
							if not sv_main_info[4+i]=='0/0':
								if not sv_main_info[0] in info_hash[sample_names[i]].keys():
									info_hash[sample_names[i]][sv_main_info[0]]={}
								if not sv_main_info[1] in info_hash[sample_names[i]][sv_main_info[0]].keys():
									info_hash[sample_names[i]][sv_main_info[0]][sv_main_info[1]]={}
								if not sv_main_info[2] in info_hash[sample_names[i]][sv_main_info[0]][sv_main_info[1]].keys():
									info_hash[sample_names[i]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]=[]
								if not sv_main_info[:4]+[sv_main_info[4+i]] in info_hash[sample_names[i]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]]:
									info_hash[sample_names[i]][sv_main_info[0]][sv_main_info[1]][sv_main_info[2]].append(sv_main_info[:4]+[sv_main_info[4+i]])
	fin.close()
	print('\t'.join(['Moleculo']+test))
	return sv_hash_order_2(info_hash)

def write_illumina_left_cluster(illumina_left_cluster,file_out,ka):
	modified_hash={}
	for k1 in illumina_left_cluster:
		algorithm_List=unify_list([i[-1] for i in k1])
		if not len(algorithm_List) in modified_hash.keys():
			modified_hash[len(algorithm_List)]=[]
		modified_hash[len(algorithm_List)].append(k1)
	for k1 in modified_hash.keys():
		file_initiate(file_out+'.num_callsers.'+str(k1))
		fo=open(file_out+'.num_callsers.'+str(k1),'a')
		rec=rec_retrive(file_out+'.num_callsers.'+str(k1))
		for k2 in modified_hash[k1]:
			rec+=1
			for k3 in k2:
				print('\t'.join([str(i) for i in k3[:-1]]+[ka,str(rec)]+[num_algorithm_index[k3[-1]]]), file=fo)
		fo.close()

def write_illumina_algorithm_index(illumina_path):
	fo=open(illumina_path+'illumina_algorithm_index','w')
	for k1 in sorted(num_algorithm_index.keys()):
		print ('\t'.join([str(i) for i in [k1,num_algorithm_index[k1]]]), file=fo)
	fo.close()

def chr_list_readin(contig):
	fin=open(contig)
	out=[]
	for line in fin:
		pin=line.strip().split()
		out+=pin
	fin.close()
	return out

def file_initiate(file_out,header=[]):
	if not os.path.isfile(file_out):
		fo=open(file_out,'w')
		if not header==[]:		print('\t'.join(header), file=fo)
		fo.close()

def vcf_file_readin(path, appdix='vcf'):
	path=path_modify(path)
	out=[]
	for i in os.listdir(path):
		if i.split('.')[-1]==appdix:
			out.append(path+i)
	return out

def vcf_name_readin_2(path,sep,key_pos):
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

def vcf_num_stat(chromo,ppre, caller_name,illumina_info_hash,illumina_stat_hash ,size_cff=50):
	#illumina_info_hash={}
	#illumina_stat_hash={}
	#########
	if caller_name=='dCGH':
		dCGH=vcf_file_readin(ppre+'dCGH/')[0]
		#dCGH=ppre+'20160930_dCGH_filt/ALL.dCGH_filt.illumina_high_coverage.20160930.cnv.genotypes.vcf'
		illumina_info_hash['dCGH']=vcf_readin_dCGH_filt(dCGH,chromo,size_cff)	
		illumina_stat_hash['dCGH']=SV_num_stat(illumina_info_hash['dCGH'])
	#########	
	if caller_name in ['GenomeStrip','GenomeSTRiP','genomestrip','gstrip']:
		GenomeSTRiP=vcf_file_readin(ppre+caller_name+'/')[0]
		#GenomeStrip=ppre+'20161001_GenomeSTRiP_CNV/HGSVC_trios.illumina.hg38.GenomeSTRiP_cnv.20161001.genotypes.vcf'
		illumina_info_hash['GenomeSTRiP']=vcf_readin_GenomeSTRiP_cnv(GenomeSTRiP,chromo,size_cff)
		illumina_stat_hash['GenomeSTRiP']=SV_num_stat(illumina_info_hash['GenomeSTRiP'])
	#########
	if caller_name in ['delly','Delly','DELLY']:
		delly=vcf_file_readin(ppre+caller_name+'/')
		#Delly=[ppre+'20161002_delly_illumina/'+i for i in ['hgsvc.delly.complex.GRCh38.20160931.high_coverage.vcf','hgsvc.delly.svs.GRCh38.20160931.high_coverage.vcf']]
		illumina_info_hash['delly']=vcf_readin_delly(delly,chromo,size_cff)
		illumina_stat_hash['delly']=SV_num_stat(illumina_info_hash['delly'])
	#########
	if caller_name in ['Manta','manta','MANTA']:
		manta=vcf_file_readin(ppre+caller_name+'/')[0]
		#Manta=ppre+'20161108_UCSD_Manta/ALL.wgs.UCSD_Manta.20162710.sv.Illumina_high-coverage_PCR-free.genotypes.vcf'
		illumina_info_hash['manta']=vcf_readin_Manta(manta,chromo,size_cff)
		illumina_stat_hash['manta']=SV_num_stat(illumina_info_hash['manta'])
	#########
	if caller_name in ['Lumpy','lumpy','LUMPY']:
		lumpy=vcf_file_readin(ppre+caller_name+'/')[0]
		#Lumpy=ppre+'20160930_lumpy/ALL.lumpy.20160930.genotypes.vcf'
		illumina_info_hash['lumpy']=vcf_readin_lumpy(lumpy,chromo,size_cff)
		illumina_stat_hash['lumpy']=SV_num_stat(illumina_info_hash['lumpy'])
	#########
	if caller_name in ['pindel','Pindel','PINDEL']:
		pindel=vcf_file_readin(ppre+caller_name+'/')[0]
		#Pindel=ppre+'20170110_ALL_pindel_illumina_sites_and_genotypes/Pindel.trios.20170110.genotypes.vcf'
		illumina_info_hash['pindel']=vcf_readin_Pindel(pindel,chromo,size_cff)
		illumina_stat_hash['pindel']=SV_num_stat(illumina_info_hash['pindel'])	
	#########
	if caller_name in ['SVelter','svelter']:
		svelter=vcf_file_readin(ppre+caller_name+'/')[0]
		#SVelter=ppre+'20170109_svelter_update/SVelter_UMich.alt_bwamem_GRCh38DH.CHS.high_coverage.vcf'
		illumina_info_hash['svelter']=vcf_readin_SVelter_UMich(svelter,chromo,size_cff)
		illumina_stat_hash['svelter']=SV_num_stat(illumina_info_hash['svelter'])
	#########
	if caller_name in ['VH','VariantHunter']:
		VH=vcf_file_readin(ppre+caller_name+'/')[0]
		#VH=ppre+'20161002_VH/ALL.VH.illumina_high_coverage.DEL.sorted.vcf'
		illumina_info_hash['VH']=vcf_readin_VH(VH,chromo,size_cff)
		illumina_stat_hash['VH']=SV_num_stat(illumina_info_hash['VH'])
	#########
	if caller_name in ['Wham','wham','WHAM']:
		Wham=vcf_file_readin(ppre+caller_name+'/')[0]
		#Wham=ppre+'20160930_wham/ALL.wham.20160930.genotypes.vcf'
		illumina_info_hash['wham']=vcf_readin_wham(Wham,chromo,size_cff)
		illumina_stat_hash['wham']=SV_num_stat(illumina_info_hash['wham'])
	#########
	if caller_name in ['ForestSV','forestsv']:
		ForestSV=vcf_file_readin(ppre+caller_name+'/')[0]
		#ForestSV=ppre+'20161108_UCSD_ForestSV/ALL.wgs.UCSD_ForestSV-gtCNV.20163110.sv.Illumina_high-coverage_PCR-free.genotypes.vcf'
		illumina_info_hash['ForestSV']=vcf_readin_ForestSV(ForestSV,chromo,size_cff)
		illumina_stat_hash['ForestSV']=SV_num_stat(illumina_info_hash['ForestSV'])
	#########
	if caller_name in ['MELT','melt']:
		MELT=vcf_file_readin(ppre+caller_name+'/')[0]
		#MELT=ppre+'20160930_MELT/ALL.melt.MEI.illumina_high_coverage.20160930.genotypes.vcf'
		illumina_info_hash['MELT']=vcf_readin_melt(MELT,chromo,size_cff)
		illumina_stat_hash['MELT']=SV_num_stat(illumina_info_hash['MELT'])
	#########
	if caller_name in ['Tardis','tardis']:
		Tardis=vcf_file_readin(ppre+caller_name+'/')[0]
		#Tardis=ppre+'20161011_Tardis_MEI_Calls/ALL.MEI.Tardis.Sorted.vcf'
		illumina_info_hash['Tardis']=vcf_readin_tardis(Tardis,chromo,size_cff)
		illumina_stat_hash['Tardis']=SV_num_stat(illumina_info_hash['Tardis'])
	#########
	if caller_name in ['novoBreak','novobreak']:
		novoBreak=vcf_name_readin_2(ppre+caller_name+'/','_',1)
		illumina_info_hash['novoBreak']=vcf_readin_kchen_novobreak(novoBreak,chromo,size_cff)
		illumina_stat_hash['novoBreak']=SV_num_stat(illumina_info_hash['novoBreak'])
	#########
	if caller_name in ['liWGS']:
		liWGS=vcf_file_readin(ppre+caller_name+'/','bed')
		#liWGS=[ppre+'20160501_liWGS_SV/'+i for i in [ 'ALL.wgs.hg38_liftOver.liWGS_SV.20160506.complexSV.3500bp_JumpingLibrary.sites.bedpe', 'ALL.wgs.hg38_liftOver.liWGS_SV.20160506.deletion.3500bp_JumpingLibrary.sites.bed', 'ALL.wgs.hg38_liftOver.liWGS_SV.20160506.duplication.3500bp_JumpingLibrary.sites.bed', 'ALL.wgs.hg38_liftOver.liWGS_SV.20160506.insertion.3500bp_JumpingLibrary.sites.bedpe', 'ALL.wgs.hg38_liftOver.liWGS_SV.20160506.inversion.3500bp_JumpingLibrary.sites.bed']]
		illumina_info_hash['liWGS']=liWGS_bed_readin(liWGS,chromo,size_cff)
		illumina_stat_hash['liWGS']=SV_num_stat(illumina_info_hash['liWGS'])
	#########
	if caller_name in ['retroCNV']:
		retroCNV=vcf_file_readin(ppre+caller_name+'/')[0]
		#retroCNV=ppre+'20161003_retroCNV/20161003_retroCNV.vcf'
		illumina_info_hash['retroCNV']=vcf_readin_retroCNV(retroCNV,chromo,size_cff)
		illumina_stat_hash['retroCNV']=SV_num_stat(illumina_info_hash['retroCNV'])
	#########
		#Moleculo_ucsd=ppre+'20161114_moleculo_pacbio_ucsd/ALL.wgs.UCSD_Moleculo.20161026.sv.Illumina_Moleculo.sites.vcf'
		#illumina_info_hash['Moleculo']=vcf_readin_Moleculo(Moleculo_ucsd,chromo,size_cff)
		#illumina_stat_hash['Moleculo']=SV_num_stat(illumina_info_hash['Moleculo'])
	#########
		#Moleculo_pb_supp=ppre+'20161114_moleculo_pacbio_ucsd/ALL.wgs.UCSD_Moleculo_PacBio-support.20161026.sv.Illumina_Moleculo.sites.vcf'
		#illumina_info_hash['Moleculo_PB']=vcf_readin_Moleculo(Moleculo_pb_supp,chromo,size_cff)
		#illumina_stat_hash['Moleculo_PB']=SV_num_stat(illumina_info_hash['Moleculo_PB'])
	#########
		#tenX_genomics=vcf_name_readin(ppre+'20160915_10x_genomics_svs/','_',0)
		#illumina_info_hash['tenX']=vcf_readin_10x(tenX_genomics,chromo,size_cff)
		#illumina_stat_hash['tenX']=SV_num_stat(illumina_info_hash['tenX'])
	#########
		#HySA=vcf_name_list_readin(ppre+'20170126_kchen_calls_hysa_regenotype/')
		#illumina_info_hash['HySA']=vcf_readin_kchen_HySA(HySA,chromo,size_cff)
		#illumina_stat_hash['HySA']=SV_num_stat(illumina_info_hash['HySA'])
	#########
	return [illumina_info_hash,illumina_stat_hash]

def write_Illumina_Info(fileout,statout, chromo,illumina_info_hash,illumina_stat_hash):
	fo=open(fileout,'a')
	fo2=open(statout,'a')
	for ka in illumina_info_hash.keys():
		for kb in illumina_info_hash[ka].keys():
			for kc in illumina_info_hash[ka][kb]:
				print('\t'.join([str(i) for i in kc+[ka,kb]]), file=fo)
	for ka in illumina_stat_hash.keys():
		for kb in illumina_stat_hash[ka].keys():
			for kc in illumina_stat_hash[ka][kb].keys():
				print('\t'.join([str(i) for i in [ka,kb,kc,illumina_stat_hash[ka][kb][kc],chromo]]), file=fo2)
	fo.close()
	fo2.close()

def main():
	parser = argparse.ArgumentParser(description='step1.standardize_vcfs_to_bed.py')
	parser.add_argument('input_path', help='directory of input vcf')
	parser.add_argument('output_file', help='name of output bed file')
	parser.add_argument('reference', help='reference genome used for bam alignment')
	parser.add_argument('caller_names', help='file containing names of algorithms to be standardized')
	parser.add_argument('contig_names', help='file containing names of contigs to be standardized')
	args = parser.parse_args()
	fileout=args.output_file
	statout=fileout.replace('.bed','.num_svs.stat')
	#eg: callers='dCGH,GenomeSTRiP,Manta,Pindel,VH,delly,liWGS,MELT,retroCNV,wham,ForestSV,lumpy,novobreak,Tardis'
	chrom_list=chr_list_readin(args.contig_names)
	callers_list=chr_list_readin(args.caller_names)
	file_initiate(fileout,['CHR','POS','END','SVTYPE','GT','CALLER','SAMPLE'])
	file_initiate(statout)
	[illumina_info_hash,illumina_stat_hash]=[{},{}]
	for chromo in chrom_list:
		for caller_name in callers_list:
			[illumina_info_hash,illumina_stat_hash]=vcf_num_stat(chromo,ppre, caller_name,illumina_info_hash,illumina_stat_hash)
		write_Illumina_Info(fileout,statout, chromo,illumina_info_hash,illumina_stat_hash)

import os
import argparse

if __name__ == '__main__':
	main()


