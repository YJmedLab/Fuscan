import argparse
import os
import pickle

__author__ = "Zhaoying Liu"
__version__ = "1.0.0"
__date__ = "2024-12-01"
__email__ = "liuzhaoying361@126.com"

class homo_list:
    def __init__(self):
        self.list = []
    
    def build_list(self, outdir, targeted_fasta, fasta):
        tmp_blat = os.path.join(outdir, "blat_tmp.psl")
        os.system(f"blat -stepSize=5 -minIdentity=90 -minScore=20 -noHead {fasta} {targeted_fasta} {tmp_blat}")
        with open(tmp_blat, "r") as f:
            for line in f.readlines():
                array = line.strip().split('\t')
                targeted_start = int(array[11])
                targeted_end = int(array[12])
                ref_chr = array[13]
                ref_start = int(array[15])
                ref_end = int(array[16])
                self.list.append([targeted_start, targeted_end, ref_chr, ref_start, ref_end])
        f.close()
        os.remove(tmp_blat)
    
    def is_in_model(self, targeted_start, targeted_end, ref_chr, ref_start, ref_end):
        is_in = False
        for region in self.list:
            threshold = 10
            homo_targeted_start = region[0] - threshold
            homo_targeted_end = region[1] + threshold
            homo_ref_chr = region[2]
            homo_ref_start = region[3] - threshold
            homo_ref_end = region[4] + threshold            
            if (targeted_end - homo_targeted_start) < 0 or (homo_targeted_end - targeted_start) < 0:
                continue
            elif ref_chr != homo_ref_chr:
                continue
            elif (ref_end - homo_ref_start) < 0 or (homo_ref_end - ref_start) < 0:
                continue
            else:
                is_in = True
        return is_in

def get_targeted_fasta(bed, fasta):
    ref_dict = {}
    lines = os.popen(f"bedtools getfasta -fi {fasta} -bed {bed} -name").read().splitlines()
    i = 0
    while i < len(lines) -1:
        if lines[i].startswith(">"):
            key = lines[i].strip(">")               
            seq = lines[i+1]
            if seq.startswith(">"):
                exit("Error: BED file contains invalid region.")
            ref_dict[key] = seq
            i += 2
        else:
            exit("Error: BED file contains invalid region.")
    return ref_dict

def targeted_prepare(fasta, bed, outdir):
    # Get targeted sequence from BED file
    print_log("Getting targeted sequence from BED file...")
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    ref_dict = get_targeted_fasta(bed, fasta)

    # Build homologous regions list for each targeted sequence
    print_log("Building homologous regions list for each targeted sequence...")
    Homo_dict = {}
    for targeted_tuple in ref_dict.items():
        key, seq = targeted_tuple
        name = key.replace("::", "_").replace(":", "_")
        targeted_fasta = os.path.join(outdir, f'{name}.targeted_seq.fasta')
        with open(targeted_fasta, "w") as f:
            f.write(f">{key}\n{seq}\n")
        print_log(f"Processing {len(seq)/1000} kb sequence...")
        Homo_list = homo_list()
        Homo_list.build_list(outdir, targeted_fasta, fasta)
        Homo_dict[targeted_tuple] = Homo_list
        os.system(f"rm {targeted_fasta}")

    # Generate pickle file for Homo_dict
    print_log("Generating pickle file for Homo_dict...")
    prefix = os.path.basename(bed).split(".")[0]
    Homo_dict_file = os.path.join(outdir,f"{prefix}.homo_dict.pkl")
    with open(Homo_dict_file, "wb") as f:
        pickle.dump(Homo_dict, f)
    print_log(f"Homo_dict file path: {Homo_dict_file}.")    
    return Homo_dict

# Run main
if __name__ == '__main__':
    def usage(name=None):                                                            
        return 'Fuscan_pre -f <FASTA> -b <BED> -o <OUTDIR>'
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Preparation of homologous-regions PKL file for Fuscan', usage=usage())
    parser.add_argument('-f','--FASTA', type=str, metavar="", required=True, default='', help='FASTA file of reference sequence')
    parser.add_argument('-b','--BED', type=str, metavar="", required=True, default='', help='BED file of targeted regions')
    parser.add_argument('-o','--OUTDIR', metavar="", type=str, default=os.getcwd(), help='Output directory for homo_dict.pkl file [default: current directory]')
    args = parser.parse_args()
    
    # Run targeted_prepare
    Homo_dict = targeted_prepare(args.FASTA, args.BED, args.OUTDIR)
