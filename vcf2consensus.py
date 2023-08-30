#!/usr/bin/python3

import os
import sys
import argparse
import random

def int_to_chr(int, chrom_pos):
    '''
    Convert integer coordinate to chromosome and position
    '''
    for chrom in chrom_pos:
        if int >= chrom_pos[chrom][0] and int <= chrom_pos[chrom][1]:
            break
    return chrom, int - chrom_pos[chrom][0] + 1

def chr_to_int(chrom, pos, chrom_pos):
    '''
    Convert chromosome and position to integer coordinate
    '''
    return chrom_pos[chrom][0] + pos - 1

def load_vcf(vcf, fasta, diff):
    '''
    Parse vcf-file and return:
    vcf_data - dictionary [integer coordinate] -> [[alleles] [sample1] [sample2] ...]
    chrom_pos - dictionary with chromosomes coordinates in integer
    genome_len - genome length
    samples - list with samples names
    '''
    vcf_file = open(vcf, 'r')
    vcf_data = {}
    
    chrom_pos = {}
    genome_len = 0
    # Get chromosomes and genome length from genome file
    if diff:
        chrom_pos, genome_len = load_genome(fasta)
    
    for line in vcf_file:
        if line.startswith('#'):
            # Get chromosomes and genome length from vcf-header
            if line.startswith('##contig') and not diff:
                contig_info = {}
                for item in line[line.find('<') + 1 : line.rfind('>')].split(','):
                    item = item.split('=')
                    contig_info[item[0]] = item[1]
                chrom_pos[contig_info['ID']] = [genome_len + 1, genome_len + int(contig_info['length'])]
                genome_len += int(contig_info['length'])
            # Get samples name
            if line.startswith('#CHROM'):
                samples = line.split()[9:]
            continue
        
        line = line.split()
        # Skip if no alt allele
        if line[4] == '.':
            continue
        # Skip if no allele depth info
        try:
            pos_AD = line[8].split(':').index('AD')
        except:
            continue
        
        int_position = chr_to_int(line[0], int(line[1]), chrom_pos)
        vcf_data[int_position] = [ [line[3], *line[4].split(',')] ]
        for sample_data in line[9:]:
            var_count = list(map(int, sample_data.split(':')[pos_AD].split(',')))
            # Return list with percent. If snp isn't covered, return 0 for all alleles
            try:
                var_perc = list(map(lambda x:x/sum(var_count), var_count))
            except ZeroDivisionError:
                var_perc = list(map(lambda x:0, var_count))
            vcf_data[int_position].append(var_perc)
    
    vcf_data[genome_len + 1] = [['N'], [0]]
    vcf_file.close()
    
    return vcf_data, chrom_pos, genome_len, samples

def get_vcf_snp(vcf_keys, cons_len):
    '''
    Create list of consensuses starts that contain alter alleles
    '''
    snps = []
    for int_position in vcf_keys:
        snps.extend(range(int_position-cons_len,int_position+1))
    return list(set(snps))

def load_genome(fasta):
    '''
    Read genome and return chrom_pos and genome_len
    '''
    chrom_pos = {}
    genome_len = 0
    chrom = ''
    fasta_file = open(fasta, 'r')
    
    for line in fasta_file:
        if line.startswith('>'):
            if chrom != '':
                chrom_pos[chrom][1] = genome_len
            chrom = line.split()[0].strip('>')
            chrom_pos[chrom] = [genome_len + 1, 0]
            continue
        genome_len += len(line.strip())
    
    chrom_pos[chrom][1] = genome_len
    fasta_file.close()
    
    return chrom_pos, genome_len

def filter_freq(perc, cutoff_freq):
    '''
    Returning allele's indices with requested frequency
    '''
    index_pass = []
    for i, value in enumerate(perc):
        if i != 0 and value > cutoff_freq:
            index_pass.append(i)
    return index_pass

if __name__ == "__main__":
    # Read and check arguments
    parser = argparse.ArgumentParser(description='Create consensus sequences from reference genome and VCF.')
    parser.add_argument('-v', dest='vcf', type=str, help='VCF file [required]', required=True)
    parser.add_argument('-g', dest='genome', type=str, help='Genome in fasta format [required]', required=True)
    parser.add_argument('-l', dest='length', type=int, help='Consensus length [required, >=10]', required=True)
    parser.add_argument('-c', dest='count', type=int, help='Consensus count [required, >=1]', required=True)
    parser.add_argument('-f', dest='freq', type=float, help='Min allele frequence [required, 0<x<=1]', required=True)
    parser.add_argument('-a', dest='alt', action='store_false', help='''Consensus must contain at least 1 snp from vcf (can be reference 
                              for some samples) [default: on]''')
    parser.add_argument('-d', dest='diff', action='store_true', help='''Use if vcf-header doesn\'t contain list of 
                              chromosomes and their length, or if chromosomes order is different from fasta file [default: off]''')
    parser.add_argument('-n', dest='n', type=int, help='Split consensuses into rows with length N [default: off, >=60]')
    parser.add_argument('-o', dest='output', type=str, help='Output name [default: output.fa]', default='output.fa')
    args = parser.parse_args()
    
    if not os.path.exists(args.vcf):
        sys.exit('VCF file not found.\nExit.')
    if not os.path.exists(args.genome):
        sys.exit('Genome file not found.\nExit.')
    if args.length < 10:
        sys.exit('Consensus length must be equal ot greater than 10.\nExit.')
    if args.count < 1:
        sys.exit('Consensus count must be equal ot greater than 1.\nExit.')
    if args.freq <= 0 or args.freq > 1:
        sys.exit('Frequency must be between 0 and 1.\nExit.')
    if args.n and args.n < 60:
        sys.exit('Rows length must be at least 60.\nExit.')
    
    # Read vcf
    vcf_data, chrom_pos, genome_len, samples = load_vcf(args.vcf, args.genome, args.diff)
    vcf_keys = sorted(vcf_data.keys())
    chrom_keys = chrom_pos.keys()
    
    # Create snps list
    if args.alt:
        snps_list = get_vcf_snp(vcf_keys, args.length)
        snps = sorted(random.sample(snps_list, min(args.count, len(snps_list))))
    else:
        snps = sorted(random.sample(range(1, genome_len + 1), min(args.count, genome_len)))
    
    fasta = open(args.genome, 'r')
    output = open(args.output, 'w')
    ref = [0, '']
    prev_snp = 0
    cur_vcf_pos = 0
    cur_fasta_pos = [0, 0]
    
    for snp in snps:
        prev_ref = ref
        ref = [snp, '']
        chrom, pos = int_to_chr(snp, chrom_pos)
        # Choose random sample for consensus
        sample_id = random.sample(range(0, len(samples)), 1)[0]
        sample_name = samples[sample_id]
        header = ['>consensus_', chrom, '_', str(pos), ' sample:', sample_name]
        
        # Search first snp after consensuses start position
        for vcf_i in range(cur_vcf_pos, len(vcf_keys)):
            if snp <= vcf_keys[vcf_i]:
                break
        cur_vcf_pos = vcf_i
        
        # Search fasta-line with consensuses start position
        if snp > cur_fasta_pos[1]:
            for line in fasta:
                line = line.strip()
                if line.startswith('>'):
                    chrom = line.split()[0].strip('>')
                    cur_fasta_chr = [chrom, chrom_pos[chrom][1]]
                    continue
                cur_fasta_pos = [cur_fasta_pos[1]+1, cur_fasta_pos[1]+len(line)]
                if cur_fasta_pos[0] <= snp <= cur_fasta_pos[1]:
                    break
        
        alt_pos = 0
        ref_pos = snp
        alt_seq = []
        ref_seq = []
        print_flag = 1
        # Create consensus
        while alt_pos <= args.length:
            # Ignore consensus if chromosome ended
            if cur_fasta_pos[1] == cur_fasta_chr[1]:
                print_flag = 0
                break
            # Get next snps from vcf
            if vcf_keys[vcf_i] < ref_pos:
                for vcf_i in range(vcf_i, len(vcf_keys)):
                    if ref_pos <= vcf_keys[vcf_i]:
                       break
            # Get next line from fasta
            if ref_pos >= cur_fasta_pos[1] + 1:
                line = fasta.readline().strip()
                cur_fasta_pos = [cur_fasta_pos[1] + 1, cur_fasta_pos[1] + len(line)]
            
            # If start was before and genomes line doesn't contain it, get sequence from previous snp reference
            if ref_pos < cur_fasta_pos[0]:
                seq = prev_ref[1]
                start = prev_ref[0]
                end = start + len(seq)
            else:
                seq = line
                start = cur_fasta_pos[0]
                end = cur_fasta_pos[1] + 1
            
            # Alternative allele or reference
            if vcf_keys[vcf_i] == ref_pos:
                ref_add = vcf_data[vcf_keys[vcf_i]][0][0]
                allele_ind = filter_freq(vcf_data[vcf_keys[vcf_i]][sample_id + 1], args.freq)
                # Check alleles with the requested frequency
                if allele_ind:
                    alt_add = vcf_data[vcf_keys[vcf_i]][0][random.sample(allele_ind,1)[0]]
                    chrom_i, pos_i = int_to_chr(vcf_keys[vcf_i], chrom_pos)
                    header.append(' ' + str(pos_i) + ':' + ref_add + '>' + alt_add)
                else:
                    alt_add = ref_add
            else:
                ref_add = seq[ref_pos-start:min(end,vcf_keys[vcf_i])-start]
                alt_add = ref_add
            
            ref_seq.append(ref_add)
            alt_seq.append(alt_add)
            ref_pos += len(ref_seq[-1])
            alt_pos += len(alt_seq[-1])
        
        # If current reference (consensus with insertions) ended before the previous, then keep the previous
        ref[1] = ''.join(ref_seq)
        if ref[0] + len(ref[1]) < prev_ref[0] + len(prev_ref[1]):
            ref = prev_ref
        
        if print_flag:
            alt_seq_join = ''.join(alt_seq)[0:args.length]
            if args.n:
                alt_seq_final = '\n'.join([alt_seq_join[i:i+args.n] for i in range(0, len(alt_seq_join), args.n)])
            else:
                alt_seq_final = alt_seq_join
            output.write(''.join([*header, '\n', alt_seq_final, '\n']))
    fasta.close()
    output.close()

