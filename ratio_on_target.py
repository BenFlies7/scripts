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
hpx_velona = '/media/partition/collected/hpx_csc_velona'
tst15_app_A = '/media/partition/collected/tst15_app/mixA'
tst15_app_B = '/media/partition/collected/tst15_app/mixB'
tst15_velona_A = '/media/partition/collected/tst15_velona/mixA'
tst15_velona_B = '/media/partition/collected/tst15_velona/mixB'

directory_list = []
directory_list.extend(hpx_surecall, hpx_velona, tst15_velona_A, tst15_velona_B, tst15_app_A, tst15_app_B)

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

if BED_TST_B:
    with open(BED_TST_B, "r") as fin:
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
list_hpx_velona = [f for f in glob.iglob(hpx_csc_velona+"/*.bam")]
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

fig,ax1 = plt.subplots()
plt.hold = True
boxes = []

chr_list = ['chrM','chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

for file in list_hpx_surecall:

    print('currently in file %s') %file

    total_mapped_counter = 0
    on_target_chr_counter = 0

    samfile = pysam.AlignmentFile(file, "rb")
    chr_lengths = samfile.lengths

    for i, chr in enumerate(chr_list):

        on_target_chr_list = []

        total_mapped_chr = samfile.count(chr_list[i], 1, chr_lengths[i])
        total_mapped_counter += total_mapped_chr

        for interval in intervals_list_hpx:
            if chr == interval.chr:
                on_target_chr = samfile.count(interval.chr, interval.start, interval.end)
                on_target_chr_counter += on_target_chr

    if total_mapped_counter > 0:
        ratio = on_target_chr_counter / float(total_mapped_counter)
        print ratio
        collected['HPX_Surecall'].append(ratio)

print collected

for file in list_hpx_velona:

    print('currently in file %s') %file

    total_mapped_counter = 0
    on_target_chr_counter = 0

    samfile = pysam.AlignmentFile(file, "rb")
    chr_lengths = samfile.lengths

    for i, chr in enumerate(chr_list):

        on_target_chr_list = []

        total_mapped_chr = samfile.count(chr_list[i], 1, chr_lengths[i])
        total_mapped_counter += total_mapped_chr

        for interval in intervals_list_hpx:
            if chr == interval.chr:
                on_target_chr = samfile.count(interval.chr, interval.start, interval.end)
                on_target_chr_counter += on_target_chr

    if total_mapped_counter > 0:
        ratio = on_target_chr_counter / float(total_mapped_counter)
        print ratio
        collected['HPX_Velona'].append(ratio)

print collected

for file in list_tst15_velona_A:

    print('currently in file %s') %file

    total_mapped_counter = 0
    on_target_chr_counter = 0

    samfile = pysam.AlignmentFile(file, "rb")
    chr_lengths = samfile.lengths

    for i, chr in enumerate(chr_list):

        on_target_chr_list = []

        total_mapped_chr = samfile.count(chr_list[i], 1, chr_lengths[i])
        total_mapped_counter += total_mapped_chr

        for interval in intervals_list_tstA:
            if chr == interval.chr:
                on_target_chr = samfile.count(interval.chr, interval.start, interval.end)
                on_target_chr_counter += on_target_chr

    if total_mapped_counter > 0:
        ratio = on_target_chr_counter / float(total_mapped_counter)
        print ratio
        collected['TST15_Velona'].append(ratio)

print collected

for file in list_tst15_velona_B:

    print('currently in file %s') %file

    total_mapped_counter = 0
    on_target_chr_counter = 0

    samfile = pysam.AlignmentFile(file, "rb")
    chr_lengths = samfile.lengths

    for i, chr in enumerate(chr_list):

        on_target_chr_list = []

        total_mapped_chr = samfile.count(chr_list[i], 1, chr_lengths[i])
        total_mapped_counter += total_mapped_chr

        for interval in intervals_list_tstB:
            if chr == interval.chr:
                on_target_chr = samfile.count(interval.chr, interval.start, interval.end)
                on_target_chr_counter += on_target_chr

    if total_mapped_counter > 0:
        ratio = on_target_chr_counter / float(total_mapped_counter)
        print ratio
        collected['TST15_Velona'].append(ratio)

print collected

for file in list_tst15_app_A:

    print('currently in file %s') %file

    total_mapped_counter = 0
    on_target_chr_counter = 0

    samfile = pysam.AlignmentFile(file, "rb")
    chr_lengths = samfile.lengths

    for i, chr in enumerate(chr_list):

        on_target_chr_list = []

        total_mapped_chr = samfile.count(chr_list[i], 1, chr_lengths[i])
        total_mapped_counter += total_mapped_chr

        for interval in intervals_list_tstA:
            if chr == interval.chr:
                on_target_chr = samfile.count(interval.chr, interval.start, interval.end)
                on_target_chr_counter += on_target_chr

    if total_mapped_counter > 0:
        ratio = on_target_chr_counter / float(total_mapped_counter)
        print ratio
        collected['TST15_App'].append(ratio)

print collected

for file in list_tst15_app_B:

    print('currently in file %s') %file

    total_mapped_counter = 0
    on_target_chr_counter = 0

    samfile = pysam.AlignmentFile(file, "rb")
    chr_lengths = samfile.lengths

    for i, chr in enumerate(chr_list):

        on_target_chr_list = []

        total_mapped_chr = samfile.count(chr_list[i], 1, chr_lengths[i])
        total_mapped_counter += total_mapped_chr

        for interval in intervals_list_tstB:
            if chr == interval.chr:
                on_target_chr = samfile.count(interval.chr, interval.start, interval.end)
                on_target_chr_counter += on_target_chr

    if total_mapped_counter > 0:
        ratio = on_target_chr_counter / float(total_mapped_counter)
        print ratio
        collected['TST15_App'].append(ratio)

print collected

for key, value in collected.items():
    boxes.append(collected[key])

bp = plt.boxplot(boxes,patch_artist=True,sym='')

for box in bp['boxes']:
    box.set(linewidth=0.1)

for whisker in bp['whiskers']:
    whisker.set(linewidth=1)

xtickNames = plt.setp(ax1, xticklabels = collected.keys())

plt.setp(xtickNames, rotation=90, fontsize=7)
plt.ylabel('Coverage (x)')
plt.xlabel('Target ID')
#plt.title('Comparison of Amplicon Depths Across Samples')
plt.show()
