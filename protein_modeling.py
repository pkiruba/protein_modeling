import os
import glob
import pandas as pd
import csv
from collections import defaultdict
import sys
import pyrosetta
pyrosetta.init()

EMBOSS = '/usr/local/emboss/bin/needle'
template_sequence_library = glob.glob(os.path.join(os.getcwd(), 'pdb_sequeuce_fasta/*.fasta'))
aligned_hits = 'template_hits.csv'
template_pdb_files = os.path.join(os.getcwd(), 'pdb_structure_files')

def sequence_alignment(target_sequence, template_sequence):
	'''Protein sequence alignment using EMBOSS program''' 
	for template_seq in template_sequence:
		target_seq_id = os.path.basename(target_sequence).split('.')[0]
		template_seq_id = os.path.basename(template_seq).split('.')[0]
		os.system('{0} -sid1 {1} -asequence {2} -sid2 {3} -bsequence {4} -gapopen 10.0 -gapextend 0.5 -aformat3 markx3 -outfile {1}_{3}.needle'.format(EMBOSS, target_seq_id, target_sequence, template_seq_id, template_seq))
		
def top_hits_from_sequence_alignment(alignend_sequence):
	'''Selecting the template hits based on the ranking of score'''
	result_table = pd.DataFrame(columns = ('query', 'template', 'length', 'identity', 'similarity', 'gaps', 'score'))
	emboss_ind_title = 1
	for alignment_file in alignend_sequence:
		with open(alignment_file, 'r') as readFILE:
			for line in readFILE:
				line = line.strip()
				if '# 1:' in line:
					query = line.split()[-1]
				elif '# 2:' in line:
					template = line.split()[-1]
				elif '# Length:' in line:
					length = line.split()[-1]
				elif '# Identity:' in line:
					identity = line.split()[-1].replace('(', '').replace(')', '').replace('%', '')
				elif '# Similarity:' in line:
					similarity = line.split()[-1].replace('(', '').replace(')', '').replace('%', '')
				elif '# Gaps:' in line:
					gaps = line[24:].replace('(', '').replace(')', '').replace('%', '')
				elif '# Score:' in line:
					score = line.split()[-1].replace('(', '').replace(')', '').replace('%', '')
			result_table.loc[emboss_ind_title] = query, template, float(length), float(identity), float(similarity), float(gaps), float(score)
			emboss_ind_title += 1
	result_table['score_rank'] = result_table['score'].rank(method = 'dense', ascending=0)
	top_hit_table = result_table.sort_values('score_rank', ascending=True)
	top_hit_table[:10].to_csv(aligned_hits, sep=',', index=False)
	
def modeling(template_hits, template_pdb_path, target_seq_path):
	'''Modeling of protein structure using PyRosetta pose manipulation application'''
	top_hit_template_file_path = []
	with open(template_hits, newline='') as csvFile:
		reader = csv.DictReader(csvFile)
		for row in reader:
			target_seq = os.path.basename(target_seq_path).split('.')[0]
			alignment_file_name = '{}_{}.needle'.format(target_seq, row['template'])
			alignment_file_path_with_name = os.path.join(os.getcwd(), alignment_file_name)
			top_hit_template_file_path.append(alignment_file_path_with_name)
	aligned_seq = defaultdict(list)
	for path in top_hit_template_file_path:
		target_template_file_name = os.path.splitext(os.path.basename(path))[0]
		target_name_fasta_format = '>{} ..'.format('_'.join(target_template_file_name.split('_')[0:2]))
		template_name_fasta_format = '>{} ..'.format('_'.join(target_template_file_name.split('_')[2:]))
		target_aligned_seq = ''
		template_aligned_seq = ''
		with open (path, 'r') as readFile:
			parse = False
			parse2 = False
			for line in readFile:
				line = line.strip()
				if not parse:
					if line.startswith(target_name_fasta_format):
						parse = True
				elif line.startswith(template_name_fasta_format):
					parse = False
				else:
					target_aligned_seq+=line

				if not parse2:
					if line.startswith(template_name_fasta_format):
						parse2 = True
				elif line.startswith('#'):
					parse2 = False
				else:
					template_aligned_seq += line
		aligned_seq[target_template_file_name].append(target_aligned_seq)
		aligned_seq[target_template_file_name].append(template_aligned_seq)
	target_seq_for_modeling = {}
	for name, alignment_file in aligned_seq.items():
		# top_hits_alignment = '{}\n{}\n{}\n\n'.format(name, alignment_file[0], alignment_file[1])
		# with open('top_hits_alignment.txt', 'a') as writeFile:
		# 	writeFile.write(top_hits_alignment)
		target_seq_based_on_temp_pdb = ''
		for i in range(len(alignment_file[0])):
			if not alignment_file[1][i] == '-':
				target_seq_based_on_temp_pdb += alignment_file[0][i]
		target_seq_for_modeling[name]=target_seq_based_on_temp_pdb
	final_target_template_for_modeling = {}
	for target_template, target_final_seq in target_seq_for_modeling.items():
		template_name = '_'.join(target_template.split('_')[2:])
		temp_list_dir = os.listdir(template_pdb_path)
		for template_hit in temp_list_dir:
			if template_name in template_hit:
				final_target_template_for_modeling[template_hit] = target_final_seq
	for template_pdb, target_seq in final_target_template_for_modeling.items():
		output_model_name = '{}_{}.pdb'.format(target_seq_path.split('.')[0], '_'.join(template_pdb.split('_')[0:2]))
		join_apo_dir_path = os.path.join(template_pdb_path, template_pdb)
		pose = pyrosetta.pose_from_file(join_apo_dir_path)
		assert(pose.size() == len(target_seq))
		scorefxn = pyrosetta.get_fa_scorefxn()
		for i in range(len(target_seq)):
			seqpos = i + 1
			name1 = target_seq[i]
			if name1 == "-":
				continue
			pyrosetta.rosetta.protocols.toolbox.pose_manipulation.repack_this_residue(seqpos, pose, scorefxn, True, name1)
		pose.dump_pdb(output_model_name)

if __name__ == "__main__":
	sequence_alignment(sys.argv[1], template_sequence_library)
	top_hits_from_sequence_alignment(glob.glob('*.needle'))
	modeling(aligned_hits, template_pdb_files, sys.argv[1])
