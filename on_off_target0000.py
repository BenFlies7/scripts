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
import seaborn as sns
import pandas as pd

#Load files
REFERENCE = "/media/partition/hg19_broadinstitute/ucsc.hg19.fasta"
BED_HPX = "/media/partition/00100-1407755742_Regions.bed"
BED_TST_A = "/media/partition/TST_15-A-manifest.bed"
BED_TST_B = "/media/partition/TST_15-B-manifest.bed"

hpx_surecall = '/media/partition/collected/hpx_csc_surecall'
hpx_csc_velona = '/media/partition/collected/hpx_csc_velona'
tst15_app_A = '/media/partition/collected/tst15_app/mixA'
tst15_app_B = '/media/partition/collected/tst15_app/mixB'
tst15_velona_A = '/media/partition/collected/tst15_velona/mixA'
tst15_velona_B = '/media/partition/collected/tst15_velona/mixB'

#Prepare BED file
IntervalColumns_hpx = namedtuple('bed', ['chr', 'start', 'end'])
intervals_list_hpx = []

IntervalColumns_tstA = namedtuple('bed', ['chr', 'start', 'end'])
intervals_list_tstA = []

IntervalColumns_tstB = namedtuple('bed', ['chr', 'start', 'end'])
intervals_list_tstB = []

if BED_HPX:
    with open(BED_HPX, "r") as fin:
        for line in fin.readlines():
            line = line.rstrip('\n')
            line = re.split(r'\t+', line.rstrip('\t'))

            if len(line) == 4:
                line[1] = int(line[1])
                line[2] = int(line[2])
                bed_line = IntervalColumns_hpx(*(line[0], int(line[1]), int(line[2])))
                intervals_list_hpx.append(bed_line)

if BED_TST_A:
    with open(BED_TST_A, "r") as fin:
        for line in fin.readlines():
            line = line.rstrip('\n')
            line = re.split(r'\t+', line.rstrip('\t'))

            #For Illumina Trusight Tumor 15 BED files
            if len(line) == 12:
                line[1] = int(line[1])
                line[2] = int(line[2])
                bed_line = IntervalColumns_tstA(*(line[0], int(line[1]), int(line[2])))
                intervals_list_tstA.append(bed_line)

if BED_TST_A:
    with open(BED_TST_A, "r") as fin:
        for line in fin.readlines():
            line = line.rstrip('\n')
            line = re.split(r'\t+', line.rstrip('\t'))

            #For Illumina Trusight Tumor 15 BED files
            if len(line) == 12:
                line[1] = int(line[1])
                line[2] = int(line[2])
                bed_line = IntervalColumns_tstB(*(line[0], int(line[1]), int(line[2])))
                intervals_list_tstB.append(bed_line)

fig = plt.figure()
ax = fig.add_subplot(111)

list_hpx_surecall = [f for f in glob.iglob(hpx_surecall+"/*.bam")]
list_hpx_csc_velona = [f for f in glob.iglob(hpx_csc_velona+"/*.bam")]
list_tst15_app_A = [f for f in glob.iglob(tst15_app_A+"/*.bam")]
list_tst15_app_B = [f for f in glob.iglob(tst15_app_B+"/*.bam")]
list_tst15_velona_A = [f for f in glob.iglob(tst15_velona_A+"/*.bam")]
list_tst15_velona_B = [f for f in glob.iglob(tst15_velona_B+"/*.bam")]

collected = {
'HPX_Surecall' : [],
'HPX_Velona' : [],
'TST15_App' : [],
'TST15_Velona' : []
}


for file in list_hpx_surecall:

    print('currently in file %s') %file

    total_mapped_counter = 0
    on_target_chr_counter = 0

    samfile = pysam.AlignmentFile(file, "rb")
    chr_lengths = samfile.lengths

    chr_list = ['chrM','chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

    for i, chr in enumerate(chr_list):

        on_target_chr_list = []

        total_mapped_chr = samfile.count(chr_list[i], 1, chr_lengths[i])
        total_mapped_counter += total_mapped_chr

        for interval in intervals_list_hpx:
            if chr == interval.chr:
                on_target_chr = samfile.count(interval.chr, interval.start, interval.end)
                on_target_chr_counter += on_target_chr

    if total_mapped_chr > 0:
        ratio = on_target_chr_counter / float(total_mapped_counter)
        print ratio
        collected['HPX_Surecall'].append(ratio)

print collected

'''
boxes = []
collected = defaultdict(list)
ratio_list = []

file_list = []

collected = {
'Sample' : [],
'Ratio' : []
}

for file in bam_file_list:

    print('currently in file %s') %file

    basename = re.sub('.bam$','',file)
    basename = re.sub(directory,'',basename)

    file_list.append(basename)

    total_mapped_counter = 0
    on_target_chr_counter = 0

    samfile = pysam.AlignmentFile(file, "rb")
    chr_lengths = samfile.lengths

    chr_list = ['chrM','chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

    for i, chr in enumerate(chr_list):

        on_target_chr_list = []

        total_mapped_chr = samfile.count(chr_list[i], 1, chr_lengths[i])
        total_mapped_counter += total_mapped_chr
        #total_mapped_chr_list.append(total_mapped_chr)

        for interval in intervals_list:
            if chr == interval.chr:
                on_target_chr = samfile.count(interval.chr, interval.start, interval.end)
                on_target_chr_counter += on_target_chr
                #on_target_chr_list.append(float(on_target_chr / total_mapped_chr))

    if total_mapped_chr > 0:
        ratio = on_target_chr_counter / float(total_mapped_counter)
        collected['Sample'].append(basename)
        collected['Ratio'].append(ratio)

bla = pd.DataFrame(collected)

print bla

g = sns.barplot('Sample','Ratio',data=bla,color='b')
g, labels = plt.xticks()
plt.setp(labels, rotation=45)
plt.show()

g = sns.boxplot(y='Ratio',data=bla)
g, labels = plt.xticks()
plt.setp(labels)
plt.show()
'''
