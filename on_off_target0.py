#!/usr/bin/python
# -*- coding: utf-8 -*-

#import modules
from collections import namedtuple, defaultdict
import pysam
import re
import sys
import matplotlib.pyplot as plt
import numpy as np
import glob


#Load files
REFERENCE = "/media/partition/hg19_broadinstitute/ucsc.hg19.fasta"
#BAM_FILE = "/media/partition/Haloplex_Test_1_Late_January/Velona/15038519_S3.bam"
INTERVALS_BED = "/media/partition/00100-1407755742_Regions.bed"
#BAM_FILE = "/media/partition/TST15_Test_1_Early_February/Base_Space/BRAF_15018040/Libraries/15018040_S1.bam"
#INTERVALS_BED = "/media/partition/TST15_Test_1_Early_February/TST_15-A-manifest.bed"
#BAM_FILE = "/media/partition/15001181_S1.bam"
#BAM_FILE = "/media/partition/Haloplex/Haloplex_Test_2_Mid_February/Velona/15028422_S6.bam"
directory = "/media/partition/hpx_csc_velona/"
#BAM_FILE = "/media/partition/Haloplex_Test_2_Mid_February/Velona/15051669_S1.bam"
#BAM_FILE = "/media/partition/Haloplex_Test_2_Mid_February/Velona/15020056_S2.bam"
#BAM_FILE = "/media/partition/Haloplex_Test_2_Mid_February/Velona/15010800_S3.bam"
#BAM_FILE = "/media/partition/Haloplex_Test_2_Mid_February/Velona/15039121_S7.bam"

IntervalColumns = namedtuple('bed', ['chr', 'start', 'end'])
intervals_list = []

if INTERVALS_BED:
    with open(INTERVALS_BED, "r") as fin:
        for line in fin.readlines():
            line = line.rstrip('\n')
            line = re.split(r'\t+', line.rstrip('\t'))
            prc = ()

            #For Agilent Haloplex BED file
            if len(line) == 4:
                line[1] = int(line[1])
                line[2] = int(line[2])
                bed_line = IntervalColumns(*(line[0], int(line[1]), int(line[2])))
                intervals_list.append(bed_line)

            #For Illumina Trusight Tumor 15 BED files
            if len(line) == 12:
                line[1] = int(line[1])
                line[2] = int(line[2])
                bed_line = IntervalColumns(*(line[0], int(line[1]), int(line[2])))
                intervals_list.append(bed_line)
else:
    print("ERROR: Provide an interval list (bed format)")

fig,ax1 = plt.subplots()
plt.hold = True

bam_file_list = [f for f in glob.iglob(directory+"/*.bam")]

boxes1 = []
boxes2 = []
collected = defaultdict(list)

for file in bam_file_list:

    samfile = pysam.AlignmentFile(file, "rb")
    chr_lengths = samfile.lengths

    chr_list = ['chrM','chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

    total_mapped_list = []
    on_target_list = []
    off_target_ratios = []
    on_target_ratios = []

    total_mapped_chr_list = []
    total_mapped_chr_counter = 0

    for i, chr in enumerate(chr_list):

        on_target_chr_counter = 0
        total_mapped_chr = samfile.count(chr_list[i], 1, chr_lengths[i])
        total_mapped_chr_counter += total_mapped_chr
        total_mapped_chr_list.append(total_mapped_chr)

        for interval in intervals_list:
            if chr == interval.chr:
                on_target_chr = samfile.count(interval.chr, interval.start, interval.end)
                on_target_chr_counter += on_target_chr / total_mapped_chr

    boxes1.append(total_mapped_chr_list)
    boxes2.append(on_target_chr_counter)

plt.boxplot(boxes1,sym='')
plt.boxplot(boxes2,sym='')
#xtickNames = plt.setp(ax1, xticklabels = collected.keys())
#plt.setp(xtickNames, rotation=90, fontsize=7)
plt.ylabel('Coverage (x)')
plt.title('Normalized Coverage Distribution per Sample')

plt.show()
'''
fig, ax = plt.subplots()
ind = np.arange(len(on_target_list))
width = 0.35

if sys.argv[1] == 'counts':
    total_map = ax.bar(ind, total_mapped_list, color = 'r')
    chr_map = ax.bar(ind + width, on_target_list, color = 'b')
    ax.set_ylabel('Reads')
    ax.set_title('Comparison of total mapped reads \n per chromosome and mapped reads in regions')
    ax.set_xticks(ind + width)
    ax.set_xticklabels(chr_list)
elif sys.argv[1] == 'frequencies':
    N = len(total_mapped_list)
    x = range(N)
    p1 = plt.bar(x, off_target_ratios, color = 'y')
    p2 = plt.bar(x, on_target_ratios, color = 'r', bottom = off_target_ratios)
    ax.set_ylabel('Frequency')
    ax.set_title('Comparison of on-target and off-target \n read frequencies per chromosome')
    ax.set_xticks(ind + width)
    ax.set_xticklabels(chr_list)
    plt.legend((p1[0], p2[0]), ('On Target', 'Off Target'))
'''
plt.show()
