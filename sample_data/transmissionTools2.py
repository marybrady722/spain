import subprocess
from subprocess import Popen, PIPE
import ntpath
import os
import shutil
import glob
import datetime
import errno
from pathlib import Path
import random


import vcf
import matplotlib.pyplot as plt
import numpy as np
import toytree
import toyplot
import toyplot.svg
import Bio 
from Bio import AlignIO


def filter_empty(vcf_path):
    '''
    goes through folder with vcfs and moves empty vcfs to
    new folder within current folder
    '''
    vcf_folder = Path(vcf_path)
    directory = "empty_vcfs"
    path = os.path.join(vcf_folder, directory)
    os.mkdir(path)
    empty_vcf_list = []
    for filename in os.listdir(vcf_folder):
        file_path = os.path.join(vcf_folder, filename)
        if filename != 'empty_vcfs':
            data = vcf.Reader(open(file_path, 'r'))
            counter = 0
            while counter == 0:
                for dummy in data:
                    counter += 1
                counter += 1
            if counter == 1:
                empty_vcf_list.append(file_path)
    return path, empty_vcf_list

def empty_move(path, empty_vcf_list):
    for filename in empty_vcf_list:
        shutil.move(filename, path)


def choose_file(input_folder):
    invalid = ['empty_vcfs', 'standard_bar_plots', 'weighted_bar_plots', 'msv', 'positions', 'bottleneck']
    vcfs = os.listdir(input_folder)
    notSelected = True
    while notSelected == True:
        randomNum1 = random.randrange(len(vcfs))
        randomNum2 = random.randrange(len(vcfs))
        if randomNum1 != randomNum2 and vcfs[randomNum1] not in invalid and vcfs[randomNum2] not in invalid:
            file1 = vcfs[randomNum1]
            file2 = vcfs[randomNum2]
            notSelected = False
    file1 = os.path.join(input_folder, file1)
    file2 = os.path.join(input_folder, file2)
    return file1, file2

def _is_valid_lfv(min_read_depth, max_AF, var_reads, total_reads):
    """
    Boolean function called as a helper for determining if a variant fits the
    user-defined parameters for parsing the data.
    """
    # Calculate allele frequency
    freq = var_reads / total_reads

    #If the allele passes restrictions, return True
    if var_reads >= min_read_depth and freq <= max_AF:
        return True

    #Otherwise, return False
    return False

def mask_parse(masks):
    """
    Parses a txt file of masked positions
    """
    mask_positions = []
    with open(masks, "r") as mask_file:
            for line in mask_file:
                comma_split = line.split(",")
                for item in comma_split:
                    nums = item.split("-")
                    if len(nums) == 1:
                        mask_positions.append(int(nums[0]))
                    else:
                        for num in range(int(nums[0]), int(nums[1])):
                            mask_positions.append(int(num))
    return mask_positions




def extract_lfv(filename, input_folder=None, min_read_depth=10, min_freq=0, max_AF=1, parse_type="biallelic", store_reference=True, masks=None, mask_type='hide'):
    """
    Extracts variant data from VCF and creates a dictionary storing data
    in the form: {position: {variant: [frequency, depth]}}.
    """

    #### Handle Errors #####
    PARSE_TYPES = {"biallelic", "multiallelic"}
    MASK_TYPES = {"hide", "highlight"}
    if min_read_depth < 0 or int(max_AF) > 1 or parse_type not in PARSE_TYPES or str(mask_type) not in MASK_TYPES:
        raise ValueError("Invalid input.")

    #########################

    if masks != None:
        mask_positions = mask_parse(masks)

    #lfv_data = Biallelic() if parse_type == "biallelic" else Multiallelic()
    lfv_data = {}
    if input_folder != None:
        filename2 = os.path.join(input_folder, filename)
        data = vcf.Reader(open(filename2, 'r'))
    else:
        data = vcf.Reader(open(filename, 'r'))
    ref_data = {}

    for record in data:

        # Calculate amount of reads support variant
        var_depth = record.INFO['RS']

        # Calculate amount of reads at particular location
        raw_depth = record.INFO['DP']

        # Find data in VCF file
        pos, var, freq = record.POS, str(record.ALT[0]), float(float(var_depth) / float(raw_depth))
        
 
        #doesn't include masked positions based on user settings
        if masks != None and pos in mask_positions and mask_type == 'hide':
            continue
            
        

        # If variant passes restrictions, store data
        if _is_valid_lfv(min_read_depth, max_AF, var_depth, raw_depth) and freq > min_freq:
            #lfv_data.store(pos, var, freq, var_depth)
            if pos not in lfv_data:
                lfv_data[pos] = {var: [freq, var_depth]}
            else:
                lfv_data[pos].update({var: [freq, var_depth]})

        if store_reference and pos not in ref_data:
            # Calculate how many reads support reference
            ref_depth = float(record.INFO['DP']) - float(record.INFO['RS'])

            ref_freq = float(ref_depth)/float(raw_depth)
            
            #Find reference allele
            ref = str(record.REF[0])
            
            # If ref allele passes restrictions, store the data
            if _is_valid_lfv(min_read_depth, max_AF, ref_depth, raw_depth) and ref_freq > min_freq:
                ref_data[pos] = {ref: [(ref_depth / raw_depth), ref_depth]}

    # After parsing is complete, make object into a dictionary
    #lfv_data = dict(lfv_data)

    # If we collected reference data, update lfv_data
    if store_reference:
        for pos in ref_data:
            if pos not in lfv_data:
                lfv_data[pos] = ref_data[pos]
            else:
                lfv_data[pos].update(ref_data[pos])

    return lfv_data

#RENAME to donor_recip
def bb_input_data(donor, recip, input_folder, min_read_depth=0, max_AF=1, parse_type="multiallelic", store_reference=True, masks=None, mask_type="hide"):
    """
    Stores info from parsing VCF files to dictionary.
    """
    #donor_file = os.path.join(input_folder, donor)
    #recip_file = os.path.join(input_folder, recip)
    donor_file = donor
    recip_file = recip
    
    donor_data = extract_lfv(donor_file, input_folder, min_read_depth=min_read_depth, max_AF=0.5, min_freq=0.02, parse_type=parse_type, store_reference=store_reference, masks=masks, mask_type=mask_type)
    recip_data = extract_lfv(recip_file, input_folder, min_read_depth=min_read_depth, max_AF=1, min_freq=0, parse_type=parse_type, store_reference=store_reference, masks=masks, mask_type=mask_type)
    
    shared_count = 0

    # Stored as {pos: {var: [donor freq., recip. freq]}} bc Maria had two bb input files 
    # and one required all this info, might change later tho
    bb_data = {} 

    if masks != None:
        mask_list = mask_parse(masks)

    # Iterate through each variant at each position
    for pos in donor_data: 

        if masks != None and pos in mask_list and mask_type == "hide":
            continue

        bb_data[pos] = {}

        for var in donor_data[pos]:

            # Store donor data
            donor_freq = donor_data[pos][var][0]
            donor_depth = donor_data[pos][var][1]
            if pos not in bb_data:
                bb_data[pos] = {var: [donor_freq, 0.0, donor_depth, 0.0]}
            else:
                bb_data[pos].update({var: [donor_freq, 0.0, donor_depth, 0.0]})

            # If recipient has same variant at same location, store it
            if pos in recip_data and var in recip_data[pos] and donor_data[pos][var][0] > 0.02 and donor_data[pos][var][0] <= 0.5:
                recip_freq = recip_data[pos][var][0]
                bb_data[pos][var][1] = recip_freq
                shared_count += 1
                recip_depth = recip_data[pos][var][1]
                bb_data[pos][var][3] = recip_depth


    return (bb_data, shared_count)

def all_pairs_parse(vcf_path, masks=None, mask_type='hide', min_read_depth=10, max_AF=1, parse_type='biallelic'):
    all_pairs = []
    invalid = ['empty_vcfs', 'standard_bar_plots', 'weighted_bar_plots', 'msv', 'positions', 'bottleneck']
    vcf_folder = Path(vcf_path)
    for donor_filename in os.listdir(vcf_folder):
        for recipient_filename in os.listdir(vcf_folder):
            if donor_filename != recipient_filename and donor_filename not in invalid and recipient_filename not in invalid:
                all_pairs.append((bb_input_data(donor_filename, recipient_filename, vcf_path, min_read_depth=min_read_depth, max_AF=max_AF, parse_type=parse_type, masks=masks, mask_type=mask_type), donor_filename, recipient_filename))
    return all_pairs

def build_majority_consensus(
    vcf_file, 
    reference, 
    masks=None, 
    mask_status='hide'
    ):
    """
    With a reference file specified, builds consensus sequence.
    """
    
    consensus = list(get_seq(reference))
    variants = extract_lfv(
        vcf_file, 
        min_freq=0, 
        max_AF=1, 
        parse_type='biallelic', 
        store_reference=True, 
        masks=None, 
        mask_type='hide'
    )

    highest_variants = {}

    for pos in variants:
        highest_freq = 0
        kept_var = None
        for var in variants[pos]:
            if variants[pos][var][0] > highest_freq:
                highest_freq = variants[pos][var][0]
                kept_var = var
        if kept_var:
            highest_variants[pos] = str(kept_var)

    for pos in highest_variants:
        consensus[pos - 1] = highest_variants[pos]

    return ''.join(consensus)

def build_minor_consensus(
    vcf_file, 
    reference, 
    min_AF=0, 
    max_AF=1,  
    masks=None, 
    mask_status='hide'
    ):
    """
    Builds consensus with variants that have second highest frequency.
    """
    consensus = list(get_seq(reference))
    data = extract_lfv(
        vcf_file, 
        min_freq=min_AF, 
        max_AF=max_AF, 
        parse_type='multiallelic', 
        store_reference=True, 
        masks=None, 
        mask_type='hide'
    )

    variants = {}

    for pos in data:
        highest_freq = 0
        max_freq = max([data[pos][var][0] for var in data[pos]])
        kept_var = None
        for var in data[pos]:
            freq = data[pos][var][0]
            if freq > highest_freq and freq != max_freq:
                highest_freq = data[pos][var][0]
                kept_var = var
        if kept_var:
            variants[pos] = str(kept_var)

    for pos in variants:
        consensus[pos - 1] = variants[pos]

    return ''.join(consensus)


def bb_file_writer(donor, recipient, input_folder, parse_type="biallelic", min_read_depth=0, max_AF=1, masks='default_mask.txt'):
    """
    Writes input file to BB_bottleneck software.
    """
    # Parse data and store outputs
    vcf_data, shared_count = bb_input_data(donor, recipient, input_folder, parse_type=parse_type, min_read_depth=min_read_depth, max_AF=max_AF, masks=masks)

    # Create shortened filenames
    fname1 = donor.split("_")[0] 
    fname2 = recipient.split("_")[0]

    # Write data to file in BB_Bottleneck input format
    filename = f"new_{fname1}_{fname2}_thred{shared_count}_complete_nofilter_bbn.txt"
    with open(filename, "w") as f:
        for pos in vcf_data:
            for var in vcf_data[pos]:
                f.write(str(vcf_data[pos][var][0])+'\t'+str(vcf_data[pos][var][1])+'\n')
                
    return filename

def all_pairs(vcf_path, input_folder, min_read_depth=0, max_AF=1, min_shared=1, output_folder=None):
    """
    Inputs:
        vcf_folder: path to folder with vcf files
    """
    invalid = ['empty_vcfs', 'standard_bar_plots', 'weighted_bar_plots', 'msv', 'positions', 'bottleneck']
    bottleneck_sizes, left_bound_intervals, right_bound_intervals, pair_nums = [], [], [], []
    pair_num = 1
    vcf_folder = Path(vcf_path)
    all_pairs_input = []
    for donor_filename in os.listdir(vcf_folder):
        for recipient_filename in os.listdir(vcf_folder):
            if donor_filename != recipient_filename and donor_filename not in invalid and recipient_filename not in invalid:
                shared_count = bb_input_data(donor_filename, recipient_filename, input_folder, min_read_depth=min_read_depth, max_AF=max_AF, masks='default_mask.txt')[1]
                if shared_count < min_shared:
                    continue
                else:
                    bb_file_writer(donor_filename, recipient_filename, input_folder, min_read_depth=min_read_depth, max_AF=max_AF)
                    donor = donor_filename.split("_")[0] 
                    recipient = recipient_filename.split("_")[0] 
                    filename = f"new_{donor}_{recipient}_thred{shared_count}_complete_nofilter_bbn.txt"
                    bb_command = str('Rscript Bottleneck_size_estimation_approx.r --file ' + filename + ' --plot_bool TRUE --var_calling_threshold 0.0001 --Nb_min 1 --Nb_max 200 --confidence_level .95')
                    pipe = subprocess.Popen(bb_command, shell=True, stdout=PIPE).stdout
                    output = str(pipe.read())
                    #print(output)
                    output_lst = output.split("\\r\\n[1]")
                    if output_lst[0] == 'b\'[1] "Bottleneck size"':
                        bottleneck_size = int(output_lst[1])
                        left_bound_interval = int(output_lst[3])
                        right = output_lst[5].split('\\r\\n')
                        right_bound_interval = int(right[0])
                        all_pairs_input.append([bottleneck_size, left_bound_interval, right_bound_interval])
                    pair_num += 1
    
    all_pairs_input.sort(key=lambda x:int(x[0]))
    pair_num = 1
    for read in all_pairs_input:
        bottleneck_sizes.append(read[0])
        left_bound_intervals.append(read[1])
        right_bound_intervals.append(read[2])
        pair_nums.append(pair_num)
        pair_num += 1
        
        
    plt.figure(figsize=(10,8))
    plt.errorbar(pair_nums,bottleneck_sizes, np.asarray([left_bound_intervals,right_bound_intervals]), ls = "None", color = "gray")
    plt.scatter(pair_nums,bottleneck_sizes,c='blue', s=32,marker="o")
    plt.title('Complete Genome (=4 shared SNVs)', fontweight='bold',fontsize=30)
    plt.xlabel('Pair ID', fontweight='bold',fontsize=26)
    plt.ylabel('Bottleneck Size', fontweight='bold',fontsize=26)
    plt.tight_layout()
    #plt.savefig('%s%s_bottleneck.png'%(vcf_folder,1), dpi=300,bbox_inches='tight')
    plt.yticks(fontsize = 24)
    plt.xticks(fontsize = 24)
    plt.show()
    if output_folder == None:
        path = os.path.join(input_folder, 'bottleneck')
        os.mkdir(path)
        filename = 'bottleneck_chart'
        save = os.path.join(path, filename)
        plt.savefig(save, dpi=300, bbox_inches='tight')
        plt.show()
    else:
        filename = '%s_%s_allele_freq.png'%(donor_filename, recipient_filename)
        save = os.path.join(output_folder, filename)
        plt.savefig(save, dpi=300, bbox_inches='tight')

def sv_counter(all_pairs_list, mask_pairs_list=None, high_shared=6):
    """
    Counts the number shared variants in all possible pairs
    Inputs: 
    - all_pairs_list: parsed list of all possible pairs
    - mask_pairs_list: parsed list of all possible pairs with masking
    - high_shared: threshold that shows pairs with this number or more vairants shared
    Outputs:
    - sv_count: a dictionary in the form {number of shared variants: number of pairs with that number of shared variants}
    - important_pairs: a list of pairs that have the threshold amount of shared variants or higher
    """
    sv_count = {}
    important_pairs = []

    #Unmasked list
    for pair in all_pairs_list:
        dict_pair = pair[0][0]
        pair_sv = 0
        for pos in dict_pair:
            for var in dict_pair[pos]:
                
                if dict_pair[pos][var][1] != 0:
                    pair_sv += 1
        if mask_pairs_list == None:
            if pair_sv in sv_count:
                sv_count[pair_sv] += 1
            else:
                sv_count[pair_sv] = 1
        else: 
            if pair_sv in sv_count:
                sv_count[pair_sv][0] += 1
            else:
                sv_count[pair_sv] = [1, 0]
        if pair_sv > high_shared:
            important_pairs.append((pair[1], pair[2]))

    #Masked list
    if mask_pairs_list != None:
        for pair in mask_pairs_list:
            dict_pair = pair[0][0]
            mask_pair_sv = 0
            for pos in dict_pair:
                for var in dict_pair[pos]:
                    if dict_pair[pos][var][1] != 0:
                        mask_pair_sv += 1
            if mask_pair_sv in sv_count:
                if len(sv_count[mask_pair_sv]) == 1:
                    sv_count[mask_pair_sv][1] = 1
                else:
                    sv_count[mask_pair_sv][1] += 1
            if mask_pair_sv > high_shared:
                important_pairs.append((pair[1], pair[2]))

    return sv_count, important_pairs

def position_counter(all_pairs_list):
    """
    Counts number of pairs that have variants at each position
    Inputs: 
    - all_pairs_list: a parsed list of all possible pairs
    Outputs: 
    - position_count: a dictionary in the form {position: number of pairs with a variant at position}
    """
    position_count = {}
    for pair in all_pairs_list:
        dict_pair = pair[0][0]
        for pos in dict_pair:
            for var in dict_pair[pos]:
                if dict_pair[pos][var][1] != 0:
                    if pos not in position_count:
                        position_count[pos] = 1
                    else:
                        position_count[pos] += 1
    return position_count



def bar_plots(vcf_path, plot_type, masks='default_mask.txt', mask_type='hide', min_read_depth=10, max_AF=1, parse_type='biallelic', min_shared=4, output_path=None):
    """
    Generates all possible bar plots of a data set and places into folder
    """
    invalid = ['empty_vcfs', 'standard_bar_plots', 'weighted_bar_plots', 'msv', 'positions', 'bottleneck']
    if output_path == None:
        folder = str(plot_type + "_bar_plots")
        path = os.path.join(vcf_path, folder)
        os.mkdir(path)
    else:
        folder = str(plot_type + "_bar_plots")
        path = os.path.join(output_path, folder)
        os.mkdir(path)
    all_pairs = []
    vcf_folder = Path(vcf_path)
    shown = 0
    for donor_filename in os.listdir(vcf_folder):
        for recipient_filename in os.listdir(vcf_folder):
            if donor_filename != recipient_filename and donor_filename not in invalid and recipient_filename not in invalid:
                if bb_input_data(donor_filename, recipient_filename, vcf_path)[1] > min_shared:
                    make_standard_bar_plot(donor_filename, recipient_filename, vcf_path, min_read_depth=min_read_depth, max_AF=1, output_folder=path, parse_type=parse_type, store_reference=True, plot_type=plot_type, masks=masks, mask_type='hide')
                    if shown <= 5:
                        plt.show()
                        shown += 1

                

def make_standard_bar_plot(donor_filename, recipient_filename, input_folder, min_read_depth=0, max_AF=1, output_folder=None, parse_type="biallelic", store_reference=True, plot_type='standard', masks='default_mask.txt', mask_type='hide'):
    """
    Saves a barplot figure showing allele frequencies vs. position.
    Inputs:
    - donor_filename: a string representation of the name of the donor's vcf file
    - recipient_filename: a string representation of the name of the recipient's vcf file
    - min_AF:  the minimum frequency at which we consider the allele to be a low frequency 
        variant; if not specified, it is defaulted to be at 0.03
    - max_AF: the maximum frequency at which we consider the allele to be a low frequency 
        variant; if not specified, it is defaulted to be at 0.5
    """
    PARSE_TYPES = {"biallelic", "multiallelic"}
    PLOT_TYPES = {"standard", "weighted"}
    MASK_TYPES = {"hide", "highlight"}
    if parse_type not in PARSE_TYPES or plot_type not in PLOT_TYPES or mask_type not in MASK_TYPES:
        raise ValueError("Invalid input.")


    vcf_data = bb_input_data(donor_filename, recipient_filename, input_folder, min_read_depth=min_read_depth, max_AF=max_AF, parse_type=parse_type, store_reference=store_reference, masks=masks, mask_type=mask_type)[0]
    positions, donor_freqs, recipient_freqs, donor_depths, recipient_depths = [], [], [], [], []

    donor_filename = donor_filename.split("_")
    recipient_filename = recipient_filename.split("_")
    donor_label = str('Donor: ' + donor_filename[0])
    recipient_label = str('Recipient: ' + recipient_filename[0])

    max_read = 0
    min_read = 1000000000000000000000000
    
    for pos in sorted(vcf_data):
        for nuc in vcf_data[pos]:
            positions.append((pos, nuc))
            donor_freqs.append(round(vcf_data[pos][nuc][0], 4))
            recipient_freqs.append(round(vcf_data[pos][nuc][1], 4))
            donor_depths.append((vcf_data[pos][nuc][2]))
            recipient_depths.append((vcf_data[pos][nuc][3]))
    
    for num in donor_depths:
        if num > max_read:
            max_read = num
        if num < min_read and num != 0:
            min_read = num
    for num in recipient_depths:
        if num > max_read:
            max_read = num
        if num < min_read and num != 0:
            min_read = num
            
    width_divider = max_read*3
    
    for num in range(0, len(donor_depths)):
        donor_depths[num] = donor_depths[num]/width_divider
    for num in range(0, len(recipient_depths)):
        recipient_depths[num] = recipient_depths[num]/width_divider
                    
    textstr = str("Allele Frequencies      Max Depth: " + str(int(max_read)) + " Min Depth: " + str(int(min_read)))

    x = np.arange(len(positions))  # the label locations
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots()

    if plot_type == 'standard':
        rects1 = ax.bar(x - width/2, donor_freqs,  color='#7f6d5f', width=width, \
            edgecolor='white', label=donor_label)
        rects2 = ax.bar(x + width/2, recipient_freqs, color='#557f2d', width=width, \
            edgecolor='white', label=recipient_label)

    if plot_type == 'weighted':
        rects1 = ax.bar(x - width/2, donor_freqs,  color='#7f6d5f', width=donor_depths, \
            edgecolor='white', label=donor_label)
        rects2 = ax.bar(x + width/2, recipient_freqs, color='#557f2d', width=recipient_depths, \
            edgecolor='white', label=recipient_label) 
    
    # Add some text for labels, title and custom x-axis tick labels, etc.
    if plot_type == 'weighted':
        ax.set_ylabel(textstr)
    else:
        ax.set_ylabel("Allele Frequencies")
    ax.set_xticks(x)
    ax.set_xticklabels(positions)
    ax.tick_params(axis = 'x', rotation = 50)
    ax.legend()

    
    maxn = len(positions)
    fig.tight_layout() 
    plt.legend(loc='lower left', bbox_to_anchor= (0.0, 1.01), ncol=2, borderaxespad=0, frameon=False, fontsize=(8))
    plt.xlim(-0.5,maxn-0.5)
    plt.tight_layout()
    if output_folder == None:
        folder = str(plot_type + "_bar_plots_samples")
        if input_folder != None:
            path = os.path.join(input_folder, folder)
        else: 
            path = folder
        os.mkdir(path)
        filename = '%s_%s_allele_freq.png'%(donor_filename, recipient_filename)
        save = os.path.join(path, filename)
        plt.savefig(save, dpi=300, bbox_inches='tight')
        plt.show()
    else:
        filename = '%s_%s_allele_freq.png'%(donor_filename, recipient_filename)
        save = os.path.join(output_folder, filename)
        plt.savefig(save, dpi=300, bbox_inches='tight')


def masked_shared_variants(sv_count, important_pairs, input_folder, output_folder=None):
    """
    Creates chart showing number of sharede variants all pairs have
    """
    
    
    shared_variants = sorted(sv_count.items())
    sv, masked_pairs, complete_pairs, cell_text, columns, masked, complete= [], [], [], [], [], [], []
    for num in shared_variants:
        sv.append(num[0])
        complete_pairs.append(num[1][0])
        complete.append(num[1][0])
        masked_pairs.append(num[1][1])
        masked.append(num[1][1])
    
    cell_text.append(masked)
    cell_text.append(complete)

    columns = tuple(sv)
    rows = ['Masked Genome', 'Complete Genome']
    colors = ['red', 'blue']


    masked_label = "Masked Genome"
    complete_label = "Complete Genome"

    x = np.arange(len(sv))  # the label locations
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots()


    rects1 = ax.bar(x - width/2, masked_pairs,  color='red', width=width, \
        edgecolor='white', label=masked_label)
    rects2 = ax.bar(x + width/2, complete_pairs, color='blue', width=width, \
        edgecolor='white', label=complete_label)

    plt.table(cellText=cell_text, rowLabels=rows, rowColours=colors, colLabels=columns)
    
    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Number of Pairs')
    ax.set_xlabel('Number of Shared Variants')
    #ax.set_xticks(x)
    #ax.set_xticklabels(sv)
    #ax.tick_params(axis = 'x')
    ax.legend()

    
    maxn = len(sv)
    fig.tight_layout()
    #plt.figure(figsize=(10,8)) 
    plt.legend(loc='lower left', bbox_to_anchor= (0.0, 1.01), ncol=2, borderaxespad=0, frameon=False, fontsize=(12))
    plt.xlim(-0.5,maxn-0.5)
    plt.xticks([])
    plt.xlabel("Shared Variants", labelpad=50)
    #plt.tight_layout()
    #plt.savefig('%s_%s_allele_freq.png'%(donor_filename, recipient_filename), dpi=300, bbox_inches='tight')
    print("Important Pairs: ", important_pairs)
    if output_folder == None:
        path = os.path.join(input_folder, 'msv')
        os.mkdir(path)
        filename = 'msv_chart'
        save = os.path.join(path, filename)
        plt.savefig(save, dpi=300, bbox_inches='tight')
        plt.show()
    else:
        filename = 'msv_chart'
        save = os.path.join(output_folder, filename)
        plt.savefig(save, dpi=300, bbox_inches='tight')


def shared_positions(position_count, mask_file, input_folder, output_folder=None, shown_variants=20):
    """
    Creates chart showing top 20 positions that have the most pairs with a variant at that position
    """
    mask_positions = mask_parse(mask_file)

    sorted_pos = sorted(position_count.items(), key=lambda x: x[1], reverse=True)
    position_data, frequency_data, colors = [], [], []
    sorted_pos = sorted_pos[0:shown_variants]
    #print(sorted_pos)
    for pos in sorted_pos:
        position_data.append(pos[0])
        if pos[0] in mask_positions:
            colors.append("red")
        else:
            colors.append("blue")
        frequency_data.append(pos[1])

    x = np.arange(len(position_data))  # the label locations
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots()

  
    rects1 = ax.bar(x - width/2, frequency_data,  color=colors, width=width, \
        edgecolor='white')
    
    masked_label = "Masked Genome"
    complete_label = "Complete Genome"
    
    labels = (complete_label, masked_label)


    
    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Number of Pairs')
    ax.set_xlabel("iSNV Position")
    ax.set_xticks(x)
    ax.set_xticklabels(position_data)
    ax.tick_params(axis = 'x', rotation = 50)
    line, = ax.plot([0, 0], label=complete_label, color='blue')
    line2, = ax.plot([0, 0], label=masked_label, color='red')
    ax.legend()

    
    maxn = len(position_data)
    fig.tight_layout() 
    plt.legend(loc='lower left', bbox_to_anchor= (0.0, 1.01), ncol=2, borderaxespad=0, frameon=False, fontsize=(12))
    plt.xlim(-0.5,maxn-0.5)
    plt.tight_layout()    
    if output_folder == None:
        path = os.path.join(input_folder, 'positions')
        os.mkdir(path)
        filename = 'positions_chart'
        save = os.path.join(path, filename)
        plt.savefig(save, dpi=300, bbox_inches='tight')
        plt.show()
    else:
        filename = 'positions_chart'
        save = os.path.join(output_folder, filename)
        plt.savefig(save, dpi=300, bbox_inches='tight')

THREADS = 2

class FastaAligner:
    """
    Class for aligning FASTA files.

    This class is meant for handling FASTA files with one sequence per file, and will raise
    an error if given a multi-FASTA file.
    """
    def __init__(self, fasta_dir):
        if os.path.isdir(fasta_dir):
            self.dir = fasta_dir
        else:
            raise ValueError("The specified path should be a directory of fasta files.")

    def align(self, reference, output_dir='', threads=THREADS):
        """
        Given a reference genome file, aligns all fasta files in directory using Parsnp.
        """
        if not os.path.exists(reference):
            raise FileNotFoundError(f"{reference}")
        if output_dir and not os.path.exists(output_dir):
            os.mkdir(output_dir)
        
        output = os.path.join(os.getcwd(), output_dir)
        cmnd = f"parsnp -d {self.dir} -r {reference} -o {output} -p {threads}"
        subprocess.call(cmnd.split())

        path2xmfa = os.path.join(output, 'parsnp.xmfa')
        cmnd2 = f"harvesttools -x {path2xmfa} -M " + os.path.join(output_dir,'parsnp.mfa')
        subprocess.call(cmnd2.split())

class FastaRecord:
    """
    Simple class for storing and accessing data from a FASTA record.
    """
    def __init__(self, name, sequence):
        if not isinstance(name, str):
            raise TypeError('Name should be a string.')
        if not isinstance(sequence, str):
            raise TypeError('Sequence should be a string.')

        self.name = name
        self.seq = sequence
    
    def set_id(self, identifier):
        if not isinstance(identifier, (str, int, float)):
            raise TypeError('ID should be a string, integer, or float.')
        self.id = identifier

class MultiFastaParser:
    """
    Makes it easier to access data in a multi-FASTA file.
    """
    def __init__(self, multifasta):
        if os.path.isfile(multifasta):
            self.fasta = multifasta
        else:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), multifasta)
        
        self.records = []

        # Populates self.records with FastaRecord objects
        with open(self.fasta, 'r') as f:
            lines = [line.strip() for line in f.readlines()]
            length = len(lines)
            for i in range(0, length - 1):
                if lines[i][0] == '>':
                    j = 1
                    while i+j < length:
                        if lines[i+j][0] == '>':
                            break
                        j += 1
                    consensus = ''.join(line for line in lines[i+1:i+j])
                    name = lines[i][1:]
                    self.records.append(FastaRecord(name, consensus))

    def get_groups(self):
        """
        Method that parses multi-FASTA file and groups records by sequence in dictionary.
        """
        groups = dict()
        for record in self.records:
            if record.seq in groups:
                groups[record.seq].add(record.name)
            else:
                groups[record.seq] = {record.name}
        return list(groups.values())

    def infer_phylogeny(self, output_dir='', label='tree', threads=THREADS, custom_cmd=''):
        """
        Runs RAxML on the multi-fasta file and stores output in output_dir.
        """
        # If directory already exists, delete all files with same label
        if output_dir:
            if os.path.exists(output_dir):
                for fname in glob.glob(os.path.join(output_dir, f"*.{label}")):
                    os.remove(os.path.join(output_dir, fname))
            else:
                os.mkdir(output_dir)

        # If user specifies command, then we run that in command line
        if custom_cmd:
            cmnd = custom_cmd

        # Otherwise, run RaxML with default arguments
        else:
            cwd = os.getcwd()
            path = os.path.join(cwd, self.fasta)
            cmnd = f"raxmlHPC -s {path} -w {os.path.join(cwd, output_dir)} \
                -n {label} -m GTRGAMMA -p {threads}"
        subprocess.call(cmnd.split())

        # Print out when finished
        msg = f'Ran RAxML on {self.fasta} and stored files in directory: {output_dir}' + '.'
        print(msg)

LINE_WRAP = 80 # Max. length of each line in fasta file 

def vcf2fasta(
    vcf_path, 
    reference, 
    output_dir="", 
    line_length=LINE_WRAP, 
    min_AF=0,
    max_AF=1,
    parse_type='biallelic',
    masks=None,
    mask_status='hide',
    consensus_type='majority'
    ):
    """
    Writes a FASTA file with a VCF file and reference FASTA file as input.
    """
    # If user specifies output directory
    if output_dir:

        # Check if directory already exists
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

    name = getpathleaf(vcf_path).split('.')[0]
    path = os.path.join(output_dir, name + '.fna')
    if consensus_type == 'majority':
        seq = build_majority_consensus(
            vcf_path, 
            reference, 
            masks=None, 
            mask_status='hide'
        )
    elif consensus_type == 'minor':
        seq = build_minor_consensus(
            vcf_path, 
            reference, 
            min_AF=0, 
            max_AF=1, 
            masks=None, 
            mask_status='hide'
        )
    else:
        raise ValueError(f'Unexpected consensus_type: {consensus_type}.')

    with open(path, 'w') as f:
        f.write(f'>{name}\n')
        for i in range(0, len(seq), line_length):
            f.write(seq[i: i+line_length] + '\n')
 
def convert_file(input_file, input_type, output_file, output_type):
    """
    Converts files using BioPython AlignIO module.
    """
    with open(input_file, 'r') as f1, open(output_file, 'w') as f2:
        seq = AlignIO.parse(f1, input_type)
        AlignIO.write(seq, f2, output_type)

def get_seq(sequence):
    """
    Returns string representation of sequence genome given a FASTA file.
    Assumes only one sequence in the file
    """    
    consensus = str()

    with open(sequence, 'r') as f:
        lines = [line.strip() for line in f.readlines()]
        for line in lines[1:]:
            consensus += line
            if line and line[0] == '>':
                raise ValueError('File should only include one sequence.')

    return consensus
    
def getpathleaf(path):
    '''
    Returns the leaf of a given path.
    For example, if inputted, home/user/.../file, it
    will return file
    '''
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)


COLORS = [
    '#DC050C', '#E8601C', '#F1932D', '#F6C141', '#F7F056', '#CAE0AB',
    '#90C987', '#4EB265', '#7BAFDE', '#5289C7', '#1965B0', '#882E72'
]

class PhyloTree:
    """
    Class for visualizing phylogenies using toytree
    """
    def __init__(self, newick, root=None, colors=COLORS):
        tree = toytree.tree(newick, tree_format=1)
        self.tree = tree
        self.colors = colors
        if root:
            self.remove_node(root)
        default_style = {
            'layout': 'r',
            'edge_type': 'p',
            'tip_labels_align': True,
            'node_labels': False,
            'node_sizes': [0 if i else 8 for i in self.tree.get_node_values(None, 1, 0)],
            'node_colors': '#FFFAFA',
            'node_style': {
                'stroke': 'black'
            },
        }
        self.style = default_style
        
    
    def update(self, **kwargs):
        static_attributes = {'layout', 'node_sizes'}
        for key in kwargs.keys():
            if key in static_attributes:
                raise ValueError(f'Cannot change attribute: {key}.')
        self.style.update(kwargs)
    
    def remove_node(self, node):
        self.tree = self.tree.drop_tips(wildcard=node)

    def color_groups(self, multifasta):
        # Parse data and group records with same sequence
        data = MultiFastaParser(multifasta)
        assign_color = dict() #{name: color}
        ref = data.records[0].name
        max_idx = len(self.colors) - 1
        idx = 0
        for group in data.get_groups():
            if ref in group:
                group.remove(ref)
            if len(group) >= 2:
                for name in group:
                    assign_color[name.split('.')[0]] = self.colors[idx%max_idx]
                idx += 1
        self.group_colors = assign_color
        
        # Traverse tree coloring nodes with same sequence
        for node in self.tree.treenode.traverse('postorder'):
            node.add_feature('color', '#FFFAFA') # default color white
            node_name = node.name.split('.')[0]
            if node_name in assign_color:
                node.add_feature('color', assign_color[node_name])
        color_data = self.tree.get_node_values('color', show_root=1, show_tips=1)
        self.update(node_colors=color_data)

    def add_heatmap(
        self, 
        vcf_dir,  
        position_range=None, 
        width=1000, 
        height=450, 
        min_AF=0, 
        masks=None, 
        mask_status='hide',
        filter_columns=True,
        store_ref=False,
        variant_type='major'
    ):
        if position_range and not isinstance(position_range, tuple):
            raise TypeError("position_range parameter must be a tuple.")

        # Create heatmap matrix from low frequency variants
        variant_pos = set() #keeps track of var. positions in all files
        matrix_data = dict() # maps filename to LFV data

        for fname in os.listdir(vcf_dir):
            invalid = ['empty_vcfs', 'standard_bar_plots', 'weighted_bar_plots', 'msv', 'positions', 'bottleneck']
            if fname not in invalid:
                path = os.path.join(vcf_dir, fname)
                data = extract_lfv(
                    path, 
                    min_freq=min_AF, 
                    max_AF=1,
                    parse_type='multiallelic',
                    store_reference=store_ref, 
                    masks=masks, 
                    mask_type=mask_status
                )
            for pos in data:
                variant_pos.add(pos)
            matrix_data[fname.split('.')[0]] = data
        variant_pos = sorted(variant_pos)

        # Create matrix and keep track of which row carries which file's data
        rows, cols = len(matrix_data), len(variant_pos)
        matrix = [[0 for j in range(cols)] for i in range(rows)]
        row_map = {idx: name.split('.')[0] for idx, name in enumerate(self.tree.get_tip_labels()[::-1])}
        column2position = {idx: pos for idx, pos in enumerate(variant_pos)}

        # Populate matrix with data
        for i, data in enumerate(matrix_data):
            name = row_map[i]
            if name[0] == "'":
                name = name[1:]
            for j, pos in column2position.items():
                if pos in matrix_data[name]:
                    for var in matrix_data[name][pos]:
                        matrix[i][j] = matrix_data[name][pos][var][0]
                else:
                    matrix[i][j] = 0.0
        
        # Keep track of columns that should be removed
        bad_columns = set()
        for j in range(cols):
            
            # filters out columns with less than 2 nonzero frequencies
            if filter_columns:
                count = 0 
                for i in range(rows):
                    if matrix[i][j] != 0:
                        count += 1
                if count < 2:
                    bad_columns.add(j)
                    
            # filters out columns not within range
            if position_range:
                if not position_range[0] <= column2position[j] <= position_range[1]:
                    bad_columns.add(j)
        
        # Label positions on heatmap
        position_labels = []
        for j in range(cols):
            if not j in bad_columns:
                position_labels.append(column2position[j])
        
        # Modify matrix to not include filtered columns
        for i in range(rows):
            new_row = []
            for j in range(cols):
                if not j in bad_columns:
                    new_row.append(matrix[i][j])

            matrix[i] = new_row
    
        rows = len(matrix)
        cols = len(matrix[0])

        if cols == 0:
            raise IOError('There are no variants to make a heatmap in given VCF files.')

        # create a canvas
        canvas = toyplot.Canvas(width=width, height=height)

        # add tree 
        tree_bounds = ('1%', '20%', '15%', '75%')
        axes = canvas.cartesian(bounds=tree_bounds)
        self.tree.draw(axes=axes, tip_labels=False, **self.style)
        
        colormap = toyplot.color.brewer.map("BlueRed", domain_min=0, domain_max=1)

        # add matrix
        matrix_bounds = ('21%', '88%', '6%', '84%')
        tlocator = toyplot.locator.Explicit(range(cols), position_labels)
        rlocator = toyplot.locator.Explicit(range(rows), self.tree.get_tip_labels()[::-1])
        canvas.matrix(
            (matrix, colormap),
            tshow=True,
            tlabel='Variant Positions',
            lshow=False,
            bounds=matrix_bounds,
            rshow=True,
            rlocator=rlocator,
            tlocator=tlocator
        )
        # add the color scale, maybe should make this optional
        xmax_range = .95 * width
        ymax_range = .8 * height
        ymin_range = .1 * height
        canvas.color_scale(
                colormap=colormap,
                x1=xmax_range,
                y1=ymax_range,
                x2=xmax_range,
                y2=ymin_range,
                width=10,
                padding=10,
                show=True,
                label="Variant Frequency",
                scale="linear",
            )

        # hide axes coordinates
        axes.show = False

        return canvas

    def draw(self):
        return self.tree.draw(**self.style)
    
    def save(self, canvas, output, render_type='html'):
        if render_type=='svg':
            import toyplot.svg
            toyplot.svg.render(canvas, output)
        elif render_type=='pdf':
            import toyplot.pdf
            toyplot.pdf.render(canvas, output)
        elif render_type=='html':
            import toyplot
            toyplot.html.render(canvas, output)
        else:
            raise ValueError('Invalid render type.')

def visualize(
    vcfdir, 
    ref, 
    colors=COLORS, 
    output_dir='TransmissionViz', 
    render_type='html', 
    threads=2,
    min_AF=0,
    masks=None,
    mask_status='hide',
    position_range=None
    ):
    '''
    Function that generates figures on the full tree and subtree phylogenies given
    a directory of VCF files and a reference genome.

    It should be noted that running this provides less customizablity than using the
    PhyloTree class methods themselves.
    '''

    # Check if path exists
    if not os.path.exists(vcfdir):
        raise FileNotFoundError(f"Directory does not exist: {vcfdir}")
    if not os.path.exists(ref):
        raise FileNotFoundError(f"File does not exist: {ref}")

    # Make new directory
    os.mkdir(output_dir)

    # Generate fasta file and put those files in tmp folder
    tmp_path = os.path.join(output_dir, "tmp")
    os.mkdir(tmp_path)
    for fn in os.listdir(vcfdir):
        fp = os.path.join(vcfdir, fn)
        vcf2fasta(
            fp, 
            ref, 
            output_dir=tmp_path,
            min_AF=min_AF,
            max_AF=1,
            parse_type='biallelic',
            masks=masks,
            mask_status=mask_status,
        )

    # Align fasta files using parsnp and put them it in parsnp folder
    parsnp = os.path.join(output_dir, 'parsnp')
    tmpdata = FastaAligner(tmp_path)
    tmpdata.align(ref, output_dir=parsnp)
    newick = os.path.join(parsnp, 'parsnp.tree')

    # Parse multifasta
    multifasta = os.path.join(parsnp, 'parsnp.mfa')
    seqs = MultiFastaParser(multifasta)
    refname = seqs.records[0].name #assumes ref seq is first record in multifasta (is the case w/ parsnp)

    # Run normal tree stuff or whatever
    tree = PhyloTree(newick, root=refname)
    tree.color_groups(multifasta)
    fig = tree.draw()[0]
    tree.save(fig,os.path.join(output_dir, 'full_tree.' + render_type))
    tree_fig = tree.add_heatmap(vcfdir, height=400, width=1000, position_range=position_range, min_AF=min_AF)
    tree.save(tree_fig, os.path.join(output_dir, 'full_tree_heatmap.' + render_type))

    # store colors assigned
    assigned_colors = tree.group_colors

    # Group vcf files into list of groups
    old_groups = list(seqs.get_groups())
    groups = list()
    for group in old_groups:
        if refname in group:
            group.remove(refname)
        if len(group) > 2:
            groups.append(group)
            
    dir_map = {name.split('.')[0]: name for name in os.listdir(vcfdir)} #maps nodes back to vcf files

    # For each subgroup, generate subtree heatmap figure
    i = 1
    for group in groups:
        group_dir = os.path.join(output_dir, f'group_{i}')
        group_tmp = os.path.join(group_dir, 'tmp')
        tmp_vcf = os.path.join(group_dir, 'vcf')
        os.mkdir(group_dir)
        os.mkdir(group_tmp)
        os.mkdir(tmp_vcf)
        for name in group:
            splitname = name.split('.')[0]
            color = assigned_colors[splitname]
            vcf = dir_map[splitname]
            vcfpath = os.path.join(vcfdir, vcf)
            shutil.copy(vcfpath, tmp_vcf)
            vcf2fasta(
                vcfpath, 
                ref, 
                output_dir=group_tmp, 
                min_AF=0, 
                max_AF=1, 
                parse_type='biallelic',
                masks=masks, 
                mask_status=mask_status,
                consensus_type='minor'
            )

        parsnp = os.path.join(group_dir, 'parsnp')
        tmpdata = FastaAligner(group_tmp)
        tmpdata.align(ref, output_dir=parsnp)

        newick = os.path.join(parsnp, 'parsnp.tree')
        multifasta = os.path.join(parsnp, 'parsnp.mfa')
        seqs = MultiFastaParser(multifasta)
        refname = seqs.records[0].name #assumes ref seq is first record in multifasta (is the case w/ parsnp)

        tree = PhyloTree(newick, root=refname)
        tree.update(node_colors=color)
        tree.draw()

        try:
            heatmap_fig = tree.add_heatmap(
                tmp_vcf, 
                height=450, 
                width=1250, 
                position_range=position_range, 
                filter_columns=True,
                store_ref=False,
                variant_type='minor'
            )
        except IOError:
            print('Not enough variants to make heatmap... continue!')
            i += 1
            continue

        tree.save(heatmap_fig, os.path.join(group_dir, f'subtree_heatmap{i}.' + render_type))
        shutil.rmtree(group_tmp, ignore_errors=True)
        i += 1

    # Get rid of tmp directory
    shutil.rmtree(tmp_path, ignore_errors=True)
