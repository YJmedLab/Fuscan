import argparse
import os
import pysam
import pickle
import gzip
import sys
from multiprocessing import Pool
from operator import itemgetter
from Fuscan_pre import *
from loguru import logger
import time

__author__ = "Zhaoying Liu"
__version__ = "1.0.0"
__date__ = "2024-12-01"
__email__ = "liuzhaoying361@126.com"

class breakpoints:
    def __init__(self):
        self.dict = {}
         
    def bulid_breakpoints(self, split_reads_dict):
        for read_name, read_list in split_reads_dict.items():
            for read_info in read_list:
                targeted_start = read_info[0]
                targeted_end = read_info[1]
                targeted_cigar = read_info[2]
                targeted_cigar_type = list(targeted_cigar.keys())
                ref_chr = read_info[3]
                ref_start = read_info[4]
                ref_end = read_info[5]
                is_reverse = read_info[6]
                if targeted_cigar_type == [0,4]:
                    targeted_breakpoint = targeted_start
                    if is_reverse:
                        ref_breakpoint = (ref_chr, ref_end)
                    else:
                        ref_breakpoint = (ref_chr, ref_start)
                elif targeted_cigar_type == [4,0]:
                    targeted_breakpoint = targeted_end
                    if is_reverse:
                        ref_breakpoint = (ref_chr, ref_start)
                    else:
                        ref_breakpoint = (ref_chr, ref_end)          
                breakpoint_tuple = (targeted_breakpoint, ref_breakpoint)
                if breakpoint_tuple not in self.dict:
                    self.dict[breakpoint_tuple] = [1, {breakpoint_tuple:1}]
                else:
                    self.dict[breakpoint_tuple][0] += 1
                    self.dict[breakpoint_tuple][1][breakpoint_tuple] += 1
    
    def show_breakpoints(self, count_threshold):
        print("targeted_breakpoint_pos\tref_chrom\tref_breakpoint_pos\tsplit_reads_count\tdiscordant_reads_count")
        for key, info in self.dict.items():
            try:
                all_count = info[2]['count'] + info[0]
            except:
                all_count = info[0]
            if all_count < count_threshold:
                continue
            else:
                try:
                    print(f"{key[0]}\t{key[1][0]}\t{key[1][1]}\t{info[0]}\t{info[2]['count']}")
                except:
                    print(f"{key[0]}\t{key[1][0]}\t{key[1][1]}\t{info[0]}\t0")
                    
    def write_results(self, targeted_info, out_file, count_threshold):
        targeted_region = targeted_info.split('::')[1]
        targeted_name = targeted_info.split('::')[0]
        if targeted_name == '':
            targeted_name = targeted_region
        targeted_chr = targeted_region.split(':')[0]
        targeted_start = int(targeted_region.split(':')[1].split('-')[0])
        with open(out_file, "w") as f:
            f.write("targeted_name\ttargeted_chr\ttargeted_breakpoint_pos\tref_chr\tref_breakpoint_pos\tsplit_reads_count\tdiscordant_reads_count\n")
            for key, info in self.dict.items():
                try:
                    all_count = info[2]['count'] + info[0]
                except:
                    all_count = info[0]
                if all_count < count_threshold:
                    continue
                else:
                    targeted_pos = targeted_start + key[0]
                    try:
                        f.write(f"{targeted_name}\t{targeted_chr}\t{targeted_pos}\t{key[1][0]}\t{key[1][1]}\t{info[0]}\t{info[2]['count']}\n")
                    except:
                        f.write(f"{targeted_name}\t{targeted_chr}\t{targeted_pos}\t{key[1][0]}\t{key[1][1]}\t{info[0]}\t0\n")

    def merge_breakpoints(self, merge_len):
        def is_merge(breakpoint_tuple1, breakpoint_tuple2, merge_len):
            target1_pos = breakpoint_tuple1[0]
            target2_pos = breakpoint_tuple2[0]
            ref1_chr = breakpoint_tuple1[1][0]
            ref2_chr = breakpoint_tuple2[1][0]
            ref1_pos = breakpoint_tuple1[1][1]
            ref2_pos = breakpoint_tuple2[1][1]
            if ref1_chr != ref2_chr:
                return False
            elif abs(target1_pos - target2_pos) > merge_len:
                return False
            elif abs(ref1_pos - ref2_pos) > merge_len:
                return False
            else:
                return True
                
        def merge(tmp_key_list):
            count = 0
            current_count = 0
            all_dict = {}
            for key in tmp_key_list:
                count += self.dict[key][0]
                all_dict[key] = self.dict[key][0]
                if self.dict[key][0] > current_count:
                    current_key = key
                    current_count = self.dict[key][0]
            return current_key, [count, all_dict]
        
        def final_merge(key_list):
            merged_dict = {}
            tmp_key_list = []
            for i in range(len(key_list)):
                if i == 0 or is_merge(key_list[i-1], key_list[i], merge_len):
                    tmp_key_list.append(key_list[i])
                else:
                    merged_key, merged_info = merge(tmp_key_list)
                    merged_dict[merged_key] = merged_info
                    tmp_key_list = [key_list[i]]
            if len(tmp_key_list) > 0:
                merged_key, merged_info = merge(tmp_key_list)
                merged_dict[merged_key] = merged_info
            self.dict = merged_dict

        final_merge(sorted(self.dict.keys(), key=itemgetter(0,1)))
        final_merge(sorted(self.dict.keys(), key=itemgetter(1,0)))
    
    def filter_breakpoints(self, discordant_reads_dict, homo_list):
        def is_support_by_discordant(key, discordant_reads_dict, threshold):
            is_support = False
            discordant_dict = {"count":0}
            targeted_breakpoint, ref_breakpoint = key
            for targeted, ref in discordant_reads_dict.items():
                targeted_start, targeted_end = targeted
                targeted_start = targeted_start - threshold
                targeted_end = targeted_end + threshold
                count = ref[0]
                ref_chr = ref[1]
                ref_start = ref[2] - threshold
                ref_end = ref[3] + threshold
                if ref_chr != ref_breakpoint[0]:
                    continue
                if (targeted_start <= targeted_breakpoint <= targeted_end) and (ref_start <= ref_breakpoint[1] <= ref_end):
                    is_support = True
                    discordant_dict["count"] += count
                    discordant_dict[targeted] = ref
            return is_support, discordant_dict

        filtered_count = 0
        del_key_list = []
        for key, info in self.dict.items():
            targeted_breakpoint, ref_breakpoint = key
            ref_chr, ref_breakpoint = ref_breakpoint
            if homo_list.is_in_model(targeted_breakpoint - 1, targeted_breakpoint, ref_chr, ref_breakpoint - 1 , ref_breakpoint):
                filtered_count += 1
                del_key_list.append(key)
            is_support, discordant_dict = is_support_by_discordant(key, discordant_reads_dict, 100)
            if is_support:
                self.dict[key].append(discordant_dict)
            else:
                filtered_count += 1
                del_key_list.append(key)
        del_key_list = list(set(del_key_list))
        for key in del_key_list:
            del self.dict[key]
        logger.info(f"[Pid: {os.getpid()}] Filtered {filtered_count} breakpoints")

def cigar_extract(seq, cigar):
    if cigar == []:
        return seq, {}
    index = 0
    out_seq = ""
    out_cigar = {}
    for op, num in cigar:
        if op == 0:
            out_seq += seq[index:index+num]
            index += num
            if 0 in out_cigar:
                out_cigar[0] += num
            else:
                out_cigar[0] = num
        elif op == 4:
            out_seq += seq[index:index+num]
            index += num
            if 4 in out_cigar:
                if num > out_cigar[4]:
                    out_cigar[0] += out_cigar[4]
                    del out_cigar[4]
                    out_cigar[4] = num
                else:
                    out_cigar[0] += num
            else:
                out_cigar[4] = num
        elif op == 1:
            out_seq += seq[index:index+num]
            index += num
            if 0 in out_cigar:
                out_cigar[0] += num
            else:
                out_cigar[0] = num
        elif op == 2:
            out_seq += "N" * num
            if 0 in out_cigar:
                out_cigar[0] += num
            else:
                out_cigar[0] = num
        else:
            if op in out_cigar:
                out_cigar[op] += num
            else:
                out_cigar[op] = num
    return out_seq, out_cigar

def cigar_str2dict(cigar_str):
    cigar_dict = {}
    i = 0
    while i < len(cigar_str):
        num = ""
        while cigar_str[i].isdigit():
            num += cigar_str[i]
            i += 1
        cigar_dict[{'M':0,'I':1,'D':2,'N':3,'S':4,'H':5,'P':6}[cigar_str[i]]] = int(num)
        i += 1
    return cigar_dict

def cigar_count_len(cigar_str):
    cigar_count = {'M':0,'I':0,'D':0,'N':0,'S':0,'H':0,'P':0}
    cigar_len = {'M':0,'I':0,'D':0,'N':0,'S':0,'H':0,'P':0}
    S_list = []
    i = 0
    while i < len(cigar_str):
        num = ""
        while cigar_str[i].isdigit():
            num += cigar_str[i]
            i += 1
        cigar_count[cigar_str[i]] += 1
        if cigar_str[i] == 'S' and cigar_count['S'] == 2:
            S_list.append(cigar_len[cigar_str[i]])
            S_list.append(int(num))
        cigar_len[cigar_str[i]] += int(num)
        i += 1
    return cigar_count, cigar_len, S_list

def split_reads_process(split_reads_reference_bam, targeted_ref_chr, targeted_ref_start, targeted_ref_end):
    split_reads_dict = {}
    logger.info(f"[Pid: {os.getpid()}] Processing: Extracting split reads...")
    for read in pysam.AlignmentFile(split_reads_reference_bam, "rb"):
        if read.is_unmapped:
            continue
        if read.mapq < 30:
            continue
        is_reverse = read.is_reverse
        targeted_read_name = read.qname.split(';')[0]
        targeted_start = int(read.qname.split(';')[2])
        targeted_end = int(read.qname.split(';')[3])
        targeted_cigar = read.qname.split(';')[1]
        targeted_cigar = cigar_str2dict(targeted_cigar)
        ref_chr = read.reference_name
        ref_start = read.reference_start
        ref_end = read.reference_end
        if (ref_chr == targeted_ref_chr) and is_intersect(ref_start, ref_end, targeted_ref_start, targeted_ref_end):
            continue
        ref_cigar = cigar_str2dict(read.cigarstring)
        if ref_cigar[0]/len(read.seq) < 0.5:
            continue
        threshold = 0.5
        if ref_cigar[0] > (targeted_cigar[4] * (1 - threshold)) and ref_cigar[0] < (targeted_cigar[4] * (1 + threshold)):
            if targeted_read_name not in split_reads_dict:
                split_reads_dict[targeted_read_name] = [[targeted_start, targeted_end, targeted_cigar, ref_chr, ref_start, ref_end, is_reverse]]
            else:
                split_reads_dict[targeted_read_name].append([targeted_start, targeted_end, targeted_cigar, ref_chr, ref_start, ref_end, is_reverse])
    return split_reads_dict

def is_intersect(start1, end1, start2, end2):
    return max(start1, start2) <= min(end1, end2)

def discordant_reads_process(discordant_reads_reference_bam, targeted_ref_chr, targeted_ref_start, targeted_ref_end):
    def is_merge(key1, list1, key2, list2, threshold):
        targeted_start1 = key1[0] - threshold
        targeted_end1 = key1[1] + threshold
        targeted_start2 = key2[0]
        targeted_end2 = key2[1]
        ref_chr1 = list1[1]
        ref_chr2 = list2[1]
        ref_start1 = list1[2] - threshold
        ref_end1 = list1[3] + threshold
        ref_start2 = list2[2]
        ref_end2 = list2[3]
        if ref_chr1 != ref_chr2:
            return False
        elif is_intersect(targeted_start1, targeted_end1, targeted_start2, targeted_end2) and is_intersect(ref_start1, ref_end1, ref_start2, ref_end2):
            return True
        else:
            return False

    def merge(tmp_key_list, reads_dict):
        count = 0
        targeted_start = []
        targeted_end = []
        ref_start = []
        ref_end = []
        for key in tmp_key_list:
            count += reads_dict[key][0]
            targeted_start.append(key[0])
            targeted_end.append(key[1])
            ref_start.append(reads_dict[key][2])
            ref_end.append(reads_dict[key][3])
        targeted_start = min(targeted_start)
        targeted_end = max(targeted_end)
        ref_start = min(ref_start)
        ref_end = max(ref_end)
        ref_chr = reads_dict[tmp_key_list[0]][1]
        return (targeted_start, targeted_end), [count, ref_chr, ref_start, ref_end]
        
    def final_merge(key_list, reads_dict):
        merged_dict = {}
        tmp_key_list = []
        for i in range(len(key_list)):
            if i == 0 or is_merge(key_list[i-1], reads_dict[key_list[i-1]], key_list[i], reads_dict[key_list[i]], 50):
                tmp_key_list.append(key_list[i])
            else:
                merged_key, merged_info = merge(tmp_key_list, reads_dict)
                merged_dict[merged_key] = merged_info
                tmp_key_list = [key_list[i]]
        if len(tmp_key_list) > 0:
            merged_key, merged_info = merge(tmp_key_list, reads_dict)
            merged_dict[merged_key] = merged_info
        return merged_dict
 
    discordant_reads_dict = {}
    logger.info(f"[Pid: {os.getpid()}] Processing: Extracting discordant reads...")
    for read in pysam.AlignmentFile(discordant_reads_reference_bam, "rb"):
        if read.is_unmapped:
            continue
        if read.mapq < 30: 
            continue
        targeted_start = int(read.qname.split(';')[1])
        targeted_end = int(read.qname.split(';')[2])
        ref_chr = read.reference_name
        ref_start = int(read.reference_start)
        ref_end = int(read.reference_end)
        if (ref_chr == targeted_ref_chr) and is_intersect(ref_start, ref_end, targeted_ref_start, targeted_ref_end):
            continue
        if (targeted_start, targeted_end) not in discordant_reads_dict:
            discordant_reads_dict[(targeted_start, targeted_end)] = [1, ref_chr, ref_start, ref_end]
        else:
            if ref_chr == discordant_reads_dict[(targeted_start, targeted_end)][1]:
                if is_intersect(ref_start, ref_end, discordant_reads_dict[(targeted_start, targeted_end)][2], discordant_reads_dict[(targeted_start, targeted_end)][3]):
                    discordant_reads_dict[(targeted_start, targeted_end)][0] += 1
                    discordant_reads_dict[(targeted_start, targeted_end)][2] = min(discordant_reads_dict[(targeted_start, targeted_end)][2], ref_start)
                    discordant_reads_dict[(targeted_start, targeted_end)][3] = max(discordant_reads_dict[(targeted_start, targeted_end)][3], ref_end)
                else:
                    discordant_reads_dict[(targeted_start, targeted_end)] = [1, ref_chr, ref_start, ref_end]
    discordant_reads_dict = final_merge(sorted(discordant_reads_dict.keys()), discordant_reads_dict)
    by_ref = sorted(discordant_reads_dict, key = lambda x:(discordant_reads_dict[x][1], discordant_reads_dict[x][2], discordant_reads_dict[x][3]))
    discordant_reads_dict = final_merge(by_ref, discordant_reads_dict)
    return discordant_reads_dict    
            
def filter_targeted_bam(out_dir, targeted_name, threads, fasta, homo_list ,targeted_bam):
    split_reads_reference_bam = os.path.join(out_dir, f'{targeted_name}.split_reads.reference.bam')
    discordant_reads_reference_bam = os.path.join(out_dir, f'{targeted_name}.discordant_reads.reference.bam')
    filtered_targeted_bam = os.path.join(out_dir, f'{targeted_name}.filtered.targeted.bam')
    tmp_fa1 = os.path.join(out_dir, f"{targeted_name}_split_tmp.fa")
    tmp_fa2 = os.path.join(out_dir, f"{targeted_name}_discordant_tmp.fa")
    targeted_ref_chr = targeted_name.split("_")[-2]
    targeted_ref_start = int(targeted_name.split("_")[-1].split("-")[0]) - 1000
    targeted_ref_end = int(targeted_name.split("_")[-1].split("-")[1]) + 1000
    mapped_reads = {}
    with open(tmp_fa1, "w") as f1, open(tmp_fa2, "w") as f2:
        for read in pysam.AlignmentFile(targeted_bam, "rb"):
            if not read.is_unmapped:
                cigar_count, cigar_len, S_list = cigar_count_len(read.cigarstring)
                if cigar_count['I'] > 2 or cigar_count['D'] > 2:
                    continue
                M_ratio = cigar_len['M']/len(read.seq)
                if cigar_count['S'] > 1:
                    if M_ratio < 0.3:
                        continue
                    if min(S_list) > 20:
                        if M_ratio < 0.5:
                            continue
                else:
                    if M_ratio < 0.15:
                        continue
                if read.mapq < 40:
                    continue
            cigar = read.cigar
            seq = read.seq
            out_seq, out_cigar = cigar_extract(seq, cigar)
            if read.qname not in mapped_reads:
                mapped_reads[read.qname] = [read.is_mapped, read.reference_start, read.reference_end, out_seq]
            else:
                paired_reads_type = [mapped_reads[read.qname][0], read.is_mapped]
                if paired_reads_type == [True, False]:
                    f2.write(f'>{read.qname};{mapped_reads[read.qname][1]};{mapped_reads[read.qname][2]}\n{out_seq}\n')
                elif paired_reads_type == [False, True]:
                    f2.write(f'>{read.qname};{read.reference_start};{read.reference_end}\n{mapped_reads[read.qname][3]}\n')
            if len(out_cigar) != 2:
                continue
            if list(out_cigar.keys()) not in [[4,0],[0,4]]:
                continue
            if out_cigar[4] < 20:
                continue
            if list(out_cigar.keys()) == [4,0]:
                split_seq = out_seq[0:out_cigar[4]]
                split_start = read.pos
                split_end = split_start + out_cigar[4]
            elif list(out_cigar.keys()) == [0,4]:
                split_seq = out_seq[out_cigar[0]:]
                split_start = read.pos + out_cigar[0]
                split_end = split_start + out_cigar[4]
            cigar_str = "".join([str(num) + {0:"M",4:"S"}[op] for op, num in out_cigar.items()])
            f1.write(f'>{read.qname};{cigar_str};{split_start};{split_end};{read.cigarstring};{read.mapq}\n{split_seq}\n')
    os.system(f"bwa mem -t {threads} {fasta} {tmp_fa1} 2>{os.devnull} | samtools view -buS -F 260 | samtools sort -o {split_reads_reference_bam}")
    os.system(f"bwa mem -t {threads} {fasta} {tmp_fa2} 2>{os.devnull} | samtools view -buS -F 260 | samtools sort -o {discordant_reads_reference_bam}")
    os.system(f"samtools index {split_reads_reference_bam}")
    os.system(f"samtools index {discordant_reads_reference_bam}")
    os.remove(tmp_fa1)
    os.remove(tmp_fa2)

    split_reads_dict = split_reads_process(split_reads_reference_bam, targeted_ref_chr, targeted_ref_start, targeted_ref_end)
    discordant_reads_dict = discordant_reads_process(discordant_reads_reference_bam, targeted_ref_chr, targeted_ref_start, targeted_ref_end)

    # Generate filtered targeted.bam
    with pysam.AlignmentFile(filtered_targeted_bam, "wb", header=pysam.AlignmentFile(targeted_bam, "rb").header) as out_bam:
        for read in pysam.AlignmentFile(targeted_bam, "rb"):
            if read.qname not in split_reads_dict:
                continue
            else:
                out_bam.write(read)
    os.system(f"samtools index {filtered_targeted_bam}")
    return split_reads_dict, discordant_reads_dict

def targeted_fusion(key):
    targeted_info = key[0]
    targeted_name = key[0].replace("::", "_").replace(":", "_")
    seq = key[1]
    Homo_list = Homo_dict[key]
    
    # Define output file names
    targeted_fasta = os.path.join(args.OUTDIR, f'{targeted_name}.targeted_seq.fasta')
    targeted_bam = os.path.join(args.OUTDIR, f'{targeted_name}.targeted.bam')

    # Generate targeted sequence fasta file
    with open(targeted_fasta, "w") as f:
        f.write(f">{targeted_info}\n{seq}\n")
    logger.info(f"[Pid: {os.getpid()}] Processing: BWA indexing for {targeted_info}...")     
    os.system(f"bwa index {targeted_fasta} > {os.devnull} 2>&1")
    
    # Index reference sequence
    if os.path.exists(args.FASTA + ".bwt") and os.path.exists(args.FASTA + ".pac") and os.path.exists(args.FASTA + ".amb") and os.path.exists(args.FASTA + ".ann") and os.path.exists(args.FASTA + ".sa"):
        pass
    else:
        logger.warning(f"[Pid: {os.getpid()}] BWA index of reference sequence not found, generating index...")
        os.system(f"bwa index {args.FASTA} > {os.devnull} 2>&1")
    logger.info(f"[Pid: {os.getpid()}] Processing: Mapping reads to targeted regions...")
    
    # Align reads to targeted regions
    os.system(f"bwa mem -t {each_thread} -B 13 -O [18,18] {targeted_fasta} {tmp_R1} {tmp_R2} 2>{os.devnull} | samtools view -buS -F 256 | samtools sort -@ {each_thread} -o {targeted_bam}")
    os.system(f"samtools index {targeted_bam}")

    # Filtering for targeted.bam
    split_reads_dict, discordant_reads_dict = filter_targeted_bam(args.OUTDIR, targeted_name, each_thread, args.FASTA, Homo_list, targeted_bam)
    
    # Run breakpoints detection
    raw_out_file = os.path.join(args.OUTDIR, f'{targeted_name}.breakpoints.raw_results.txt')
    Breakpoints = breakpoints()
    Breakpoints.bulid_breakpoints(split_reads_dict)
    Breakpoints.merge_breakpoints(10)
    Breakpoints.write_results(targeted_info, raw_out_file, 0)
    Breakpoints.filter_breakpoints(discordant_reads_dict, Homo_list)
    out_file = os.path.join(args.OUTDIR, f'{targeted_name}.breakpoints.results.txt')
    Breakpoints.write_results(targeted_info, out_file, 3)
    Breakpoints.show_breakpoints(3)

def write_reads_name(read_name):
    with open(os.path.join(args.OUTDIR,"tmp_reads_name.txt"), "a") as f:
        f.write(read_name)

def filter_reads(R):
    each_thread = int(args.THREADS) // 2
    reads_name = os.path.join(args.OUTDIR, "tmp_reads_name.txt")
    if R == "R1":
        os.system(f'seqkit grep -f {reads_name} -j {each_thread} {args.R1} > {tmp_R1}')
    elif R == "R2":
        os.system(f'seqkit grep -f {reads_name} -j {each_thread} {args.R2} > {tmp_R2}')

def annotate_fusion(key):
    def getRank(rankList):
        if rankList == "NA":
            return "NA"
        try:
            if "Intron" in rankList[0]:
                rankInfo = rankList[0] + rankList[1]
            elif "Exon" in rankList[0]:
                rankInfo = rankList[0] + rankList[1]
            elif "Transcript" in rankList[0]:
                rankInfo = "NA"
            else:
                rankInfo = rankList[0]
        except:
            rankInfo = "NA"
        return rankInfo

    def getGene(geneList):
        if "Gene" in geneList:
            geneName = geneList[1]
        else:
            geneName = "NA"
        return geneName

    def getInfo(closest):
        elements = closest.split('|')[1:]
        deeply_nested_elements = []
        for element in elements:
            sub_elements = element.split(',')
            deeply_split_elements = [sub_elem.split(':') for sub_elem in sub_elements]
            deeply_nested_elements.append(deeply_split_elements)
        def filter_anno_list(anno_list):
            protein_coding = [item for item in anno_list if len(item) > 2 and len(item[2]) > 0 and item[2][-1] == 'protein_coding']
            pseudogene = [item for item in anno_list if len(item) > 2 and len(item[2]) > 0 and item[2][-1] == 'pseudogene']
            combined_list = protein_coding + pseudogene
            intron = [item for item in protein_coding if len(item) > 0 and item[0][0] == 'Intron']
            exon = [item for item in protein_coding if len(item) > 0 and item[0][0] == 'Exon']
            retained_intron = [item for item in intron if len(item[0]) > 3 and 'RETAINED-RETAINED' in item[0]]
            if retained_intron:
                return retained_intron[0]
            elif intron:
                return intron[0]
            elif exon:
                return exon[0]
            elif protein_coding:
                return protein_coding[0]
            elif pseudogene:
                return pseudogene[0]
            elif combined_list:
                return combined_list[0]
            elif len(anno_list) > 0:
                return anno_list[0]
            else:
                return ['NA','NA']
        return filter_anno_list(deeply_nested_elements)
        
    targeted_name = key[0].replace("::", "_").replace(":", "_")
    input_file = os.path.join(args.OUTDIR, f'{targeted_name}.breakpoints.results.txt')
    output_file = os.path.join(args.OUTDIR, f'{targeted_name}.breakpoints.annotated.results.txt')
    tmp_vcf = os.path.join(args.OUTDIR, f'{targeted_name}.tmp.breakpoints.vcf')
    tmp_anno_vcf = os.path.join(args.OUTDIR, f'{targeted_name}.tmp.annotated.breakpoints.vcf')
    with open(input_file, "r") as in_f:
        if len(in_f.readlines()) == 1:
            return
    with open(input_file, "r") as in_f, open(tmp_vcf, "w") as out_f:
        out_f.write("##fileformat=VCFv4.1\n##INFO=<ID=EDC,Number=.,Type=String,Description=\"Evidence\">\n##INFO=<ID=FREQ,Number=.,Type=String,Description=\"Frequency\">\n##INFO=<ID=NOTE,Number=.,Type=String,Description=\"Hottag\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for line in in_f.readlines():
            if line.startswith('targeted_name'):
                continue
            if line == '':
                return
            targeted_chr = line.split('\t')[1]
            targeted_breakpoint_pos = int(line.split('\t')[2])
            ref_chr = line.split('\t')[3]
            ref_breakpoint_pos = int(line.split('\t')[4])
            split_reads_count = int(line.split('\t')[5])
            discordant_reads_count = int(line.split('\t')[6])
            note = f"{split_reads_count},{discordant_reads_count}"
            out_f.write(f"{targeted_chr}\t{targeted_breakpoint_pos}\t.\t.\t<INV>\t.\t.\tEDC=NA;FREQ=100%;NOTE={note}\n")
            out_f.write(f"{ref_chr}\t{ref_breakpoint_pos}\t.\t.\t<INV>\t.\t.\tEDC=NA;FREQ=100%;NOTE={note}\n")
    snpeff_dir = os.path.join(os.path.dirname(sys.argv[0]),"../snpeff")
    snpeff_java = os.path.join(snpeff_dir,"snpEff.jar")
    snpeff_config = os.path.join(snpeff_dir,"snpEff.config")
    os.system(f"java -jar {snpeff_java} closest -canon -c {snpeff_config} hg19 {tmp_vcf} > {tmp_anno_vcf} 2>{os.devnull}")
    
    i = 0
    with open(output_file, "w") as f:
        f.write("targeted_name\ttargeted_chr\ttargeted_breakpoint_pos\ttargeted_gene\ttargeted_rank\tref_chr\tref_breakpoint_pos\tref_gene\tref_rank\tsplit_reads_count\tdiscordant_reads_count\n")
        for vcf in pysam.VariantFile(tmp_anno_vcf, "r"):
            i += 1
            if vcf.info.get("CLOSEST") == None:
                continue
            closest = ",".join(vcf.info.get("CLOSEST"))
            annoList = getInfo(closest)
            if i % 2 == 1:
                gene1 = getGene(annoList[-1])
                rank1 = getRank(annoList[0])
                targeted_chr = vcf.chrom
                targeted_breakpoint_pos = vcf.pos
            else:
                gene2 = getGene(annoList[-1])
                rank2 = getRank(annoList[0])
                ref_chr = vcf.chrom
                ref_breakpoint_pos = vcf.pos
                if gene1 == gene2:
                    continue
                split_reads_count, discordant_reads_count = vcf.info.get("NOTE")
                f.write(f"{targeted_name}\t{targeted_chr}\t{targeted_breakpoint_pos}\t{gene1}\t{rank1}\t{ref_chr}\t{ref_breakpoint_pos}\t{gene2}\t{rank2}\t{split_reads_count}\t{discordant_reads_count}\n")
    os.remove(tmp_vcf)
    os.remove(tmp_anno_vcf)

def depth_count(bam_file, pos):
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        ref = bam.references[0]
        if pos == 0:
            pos += 1
        acgt_depth = bam.count_coverage(ref, pos-1, pos, quality_threshold=0)
    return sum([i[0] for i in acgt_depth])

def results_summary(out_dir, background_dict, fusion_pair_dict):
    split_interest = int(args.THRESHOLD.split(",")[0])
    discordant_interest = int(args.THRESHOLD.split(",")[1])
    split_intron = int(args.THRESHOLD.split(",")[2])
    discordant_intron = int(args.THRESHOLD.split(",")[3])
    split_exon = int(args.THRESHOLD.split(",")[4])
    discordant_exon = int(args.THRESHOLD.split(",")[5])
    split_genenic = int(args.THRESHOLD.split(",")[6])
    discordant_genenic = int(args.THRESHOLD.split(",")[7])

    with open(os.path.join(out_dir, "results_summary.txt"), "w") as f:
        f.write("Targeted_name\tTargeted_chr\tTargeted_breakpoint_pos\tTargeted_gene\tTargeted_rank\tRef_chr\tRef_breakpoint_pos\tRef_gene\tRef_rank\tSplit_reads_count\tDiscordant_reads_count\t\tImproper_ratio\tBackground\n")
        for line in os.popen(f"cat {out_dir}/*.breakpoints.annotated.results.txt | sort -k6,6nr").read().split("\n"):
            if line.startswith('targeted_name'):
                continue
            if line == '':
                continue
            targeted_name = line.split('\t')[0]
            targeted_chr = line.split('\t')[1]
            targeted_breakpoint_pos	= int(line.split('\t')[2])
            targeted_gene = line.split('\t')[3]
            targeted_rank = line.split('\t')[4]
            ref_chr = line.split('\t')[5]
            ref_breakpoint_pos = int(line.split('\t')[6])
            ref_gene = line.split('\t')[7]
            ref_rank = line.split('\t')[8]
            split_reads_count = int(line.split('\t')[9])
            discordant_reads_count = int(line.split('\t')[10])
            if fusion_pair_dict == None:
                if "Intron" in ref_rank:
                    if split_reads_count < split_intron or discordant_reads_count < discordant_intron:
                        continue
                elif "Exon" in ref_rank:
                    if split_reads_count < split_exon or discordant_reads_count < discordant_exon:
                        continue
                else:
                    if split_reads_count < split_genenic or discordant_reads_count < discordant_genenic:
                        continue
            else:
                if targeted_gene in fusion_pair_dict and ref_gene in fusion_pair_dict[targeted_gene]:
                    if split_reads_count < split_interest or discordant_reads_count < discordant_interest:
                        continue
                elif ref_gene in fusion_pair_dict and targeted_gene in fusion_pair_dict[ref_gene]:
                    if split_reads_count < split_interest or discordant_reads_count < discordant_interest:
                        continue
                else:
                    if "Intron" in ref_rank:
                        if split_reads_count < split_intron or discordant_reads_count < discordant_intron:
                            continue
                    elif "Exon" in ref_rank:
                        if split_reads_count < split_exon or discordant_reads_count < discordant_exon:
                            continue
                    else:
                        if split_reads_count < split_genenic or discordant_reads_count < discordant_genenic:
                            continue
            targeted_bam = os.path.join(out_dir, f"{targeted_name}.targeted.bam")
            targeted_pos = targeted_breakpoint_pos - int(targeted_name.split("_")[-1].split("-")[0])
            try:
                improper_ratio = round(split_reads_count / depth_count(targeted_bam, targeted_pos),4)
            except:
                improper_ratio = 0
            line += f"\t{improper_ratio}"
            if improper_ratio < 0.05:
                continue
            if background_dict == None:
                f.write(line + "\tPASS\n")
                continue
            if (ref_chr, ref_breakpoint_pos) in background_dict['2R']:
                f.write(line + "\t2R\n")
            elif (ref_chr, ref_breakpoint_pos) in background_dict:
                if (targeted_chr, targeted_breakpoint_pos) in background_dict[(ref_chr, ref_breakpoint_pos)][1]:
                    f.write(line + "\tBIB\n")
                else:
                    f.write(line + "\tIB\n")
            else:
                f.write(line + "\tPASS\n")
            
# Run main
if __name__ == '__main__':
    start_time = time.time()
    def usage(name=None):                                                            
        return 'Fuscan -f <FASTA> -b <BED> -R1 <R1> -R2 <R2> -o <OUTDIR>'
    # Parser command line arguments
    parser = argparse.ArgumentParser(description='An ultra-sensitive and noise-free fusion detector for genomic breakpoints discovery in DNA sequencing data.', usage=usage())
    parser.add_argument('-f','--FASTA', type=str, metavar="", required=True, default='', help='FASTA file of reference sequence')
    parser.add_argument('-b','--BED', type=str, metavar="", required=False, default='', help='BED file of targeted regions of interested genes')
    parser.add_argument('-p','--PKL', type=str, metavar="", required=False, default='', help='PKL file generated from Fuscan_pre')
    parser.add_argument('-R1', type=str, metavar="", required=True, default='', help='FASTQ or FASTQ.gz file of R1 reads')
    parser.add_argument('-R2', type=str, metavar="", required=True, default='', help='FASTQ or FASTQ.gz file of R2 reads')
    parser.add_argument('-o','--OUTDIR', metavar="", type=str, default=os.getcwd(), help='Output directory for results [default: current directory]')
    parser.add_argument('-bg','--BG', metavar="", type=str, default=None, help='Background PKL file generated from Fuscan_bg [default: None]')
    parser.add_argument('-fp','--FUSION_PAIR', metavar="", type=str, default=None, help='Tab-separated file of interested fusion pairs [default: None]')
    parser.add_argument('-ts','--THRESHOLD', metavar="", type=str, default='1,1,10,10,15,15,20,20', help='Threshhold of split reads and discordant reads count for interested fusion pairs, intron, exon and intergenic region [default: 1,1,10,10,15,15,20,20]')
    parser.add_argument('-t','--THREADS', metavar="", type=str, default='1', help='Number of threads to use [default: 1]')
    parser.add_argument('-v','--version', action='version', version="Fuscan version: " + __version__)
    args = parser.parse_args()
    logger.add(os.path.join(args.OUTDIR, "processing.log"), format='{time:YYYY-MM-DD HH:mm:ss} [{level}] {message}', enqueue=True)
    logger.add(sys.stderr, format='{time:YYYY-MM-DD HH:mm:ss} [{level}] {message}', enqueue=True)
    logger.remove(0)
    # Generate output directory
    if not os.path.exists(args.OUTDIR):
        os.makedirs(args.OUTDIR)
    
    # Load homologous region list
    if args.PKL == '':
        if args.BED == '':
            logger.error("Please provide a BED file and a FASTA file for reference sequence.")
            exit()
        else:
            Homo_dict = targeted_prepare(args.FASTA, args.BED, args.OUTDIR)
    else:
        with open(args.PKL, "rb") as f:
            Homo_dict = pickle.load(f)
    
    # Filter reference sequence
    logger.info("Processing: Mapping reads to reference sequence...")
    tmp_ref_bam = os.path.join(args.OUTDIR, "tmp.ref.bam")
    os.system(f"bwa mem -t {args.THREADS} -L [60,60] {args.FASTA} {args.R1} {args.R2} 2>>{os.path.join(args.OUTDIR, "processing.log")} | samtools sort - -@ {args.THREADS} -n -o {tmp_ref_bam}")

    logger.info("Processing: Mapping reads finished.")
    logger.info("Processing: Selecting improper reads...")
    save_set = set()
    with pysam.AlignmentFile(tmp_ref_bam, "rb") as f:
        last_read_chr = ""
        last_read_name = ""
        for read in f:
            if read.qname in save_set:
                last_read_chr = read.reference_name
                last_read_name = read.qname
                continue 
            if abs(read.isize) > 1000 or read.mapq < 30 or (read.qname == last_read_name and read.reference_name != last_read_chr and last_read_name != ""):
                save_set.add(read.qname + "\n")
            last_read_chr = read.reference_name
            last_read_name = read.qname
    os.remove(tmp_ref_bam)
    logger.info("Processing: Writing improper reads to file...")
    tmp_R1 = os.path.join(args.OUTDIR, "tmp.R1.fq.gz")
    tmp_R2 = os.path.join(args.OUTDIR, "tmp.R2.fq.gz")

    # Run targeted fusion detection
    each_thread = int(args.THREADS) // len(Homo_dict.keys())
    pool = Pool(int(args.THREADS))
    pool.map(write_reads_name, save_set)
    pool.map(filter_reads, ["R1", "R2"])
    logger.info("Processing: Running targeted fusion detection...")
    pool.map(targeted_fusion, Homo_dict.keys())
    logger.info("Processing: Annotating fusion results...")
    pool.map(annotate_fusion, Homo_dict.keys())
    pool.close()
    pool.join()
    os.remove(os.path.join(args.OUTDIR,"tmp_reads_name.txt"))
    os.remove(tmp_R1)
    os.remove(tmp_R2)
    
    # Load background_dict
    if args.BG == None:
        background_dict = None
    else:
        with open(args.BG, "rb") as f:
            background_dict = pickle.load(f)
    
    # Load fusion_pair_dict
    if args.FUSION_PAIR == None:
        fusion_pair_dict = None
    else:
        fusion_pair_dict = {}
        with open(args.FUSION_PAIR, "r") as f:
            for line in f:
                line = line.strip().split("\t")
                pair_set = set()
                for pair in line[1].split(","):
                    pair_set.add(pair)
                fusion_pair_dict[line[0]] = pair_set
        
    # Generate results summary
    logger.info("Processing: Generating results summary...")
    results_summary(args.OUTDIR, background_dict, fusion_pair_dict)
    end_time = time.time()
    round_time = round(end_time - start_time, 2)
    logger.info(f"Elapsed time: {round_time} seconds.")
