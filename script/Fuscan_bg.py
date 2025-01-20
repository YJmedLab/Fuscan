import os
import argparse
import pickle

__author__ = "Zhaoying Liu"
__version__ = "1.0.0"
__date__ = "2024-12-01"
__email__ = "liuzhaoying361@126.com"

def Fuscan_bg(fusion_pair_dict):
    background_dict = {"2R":set()}
    lines = os.popen(f'cat {os.path.join(args.INPUTDIR,"*/*.breakpoints.annotated.results.txt")}').read().split("\n")
    for line in lines:
        if line.startswith('targeted_name'):
            continue
        else:
            if line == '':
                continue
            targeted_chr = line.split('\t')[1]
            targeted_breakpoint_pos	= int(line.split('\t')[2])
            targeted_gene = line.split('\t')[3]
            ref_chr = line.split('\t')[5]
            ref_breakpoint_pos = int(line.split('\t')[6])
            ref_gene = line.split('\t')[7]
            split_reads_count = int(line.split('\t')[9])
            discordant_reads_count = int(line.split('\t')[10])
            try:
                if targeted_gene in fusion_pair_dict and ref_gene in fusion_pair_dict[targeted_gene]:
                    continue
                elif ref_gene in fusion_pair_dict and targeted_gene in fusion_pair_dict[ref_gene]:
                    continue
            except:
                pass
            if (ref_chr, ref_breakpoint_pos) not in background_dict:
                background_dict[(ref_chr, ref_breakpoint_pos)] = [1, set(), [split_reads_count], [discordant_reads_count], [split_reads_count + discordant_reads_count]]
                background_dict[(ref_chr, ref_breakpoint_pos)][1].add((targeted_chr, targeted_breakpoint_pos))
            else:
                background_dict[(ref_chr, ref_breakpoint_pos)][0] += 1
                background_dict[(ref_chr, ref_breakpoint_pos)][1].add((targeted_chr, targeted_breakpoint_pos))
                background_dict[(ref_chr, ref_breakpoint_pos)][2].append(split_reads_count)
                background_dict[(ref_chr, ref_breakpoint_pos)][3].append(discordant_reads_count)
                background_dict[(ref_chr, ref_breakpoint_pos)][4].append(split_reads_count + discordant_reads_count)
                if background_dict[(ref_chr, ref_breakpoint_pos)][0] == 2:
                    background_dict["2R"].add((ref_chr, ref_breakpoint_pos))
    return background_dict

if __name__ == '__main__':
    def usage(name=None):                                                            
        return 'Fuscan_bg -i <INPUTDIR> -fp <FUSION_PAIR> -o <OUTDIR>'
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Preparation of background PKL file for Fuscan', usage=usage())
    parser.add_argument('-fp','--FUSION_PAIR', metavar="", type=str, default=None, help='Interested fusion pairs to be excluded in background [default: None]')
    parser.add_argument('-i','--INPUTDIR', type=str, metavar="", required=True, default='', help='Input directory containing results folders from Fuscan')
    parser.add_argument('-o','--OUTDIR', metavar="", type=str, default=os.getcwd(), help='Output directory for background_dict.pkl [default: current directory]')
    args = parser.parse_args()
    
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
    # Run Fuscan_bg function
    background_dict = Fuscan_bg(fusion_pair_dict)

    # Save background_dict as a pickle file
    with open(os.path.join(args.OUTDIR, 'background_dict.pkl'), 'wb') as f:
        pickle.dump(background_dict, f)
        