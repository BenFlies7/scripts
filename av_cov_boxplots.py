#!/usr/bin/python
# -*- coding: utf-8 -*-

#import modules
import pybedtools
from collections import namedtuple, defaultdict
import re
import matplotlib.pyplot as plt
import glob

IntervalColumns = namedtuple('bed', ['chr', 'start', 'end', 'gene'])
intervals_list = []

INTERVALS_BED = "/media/partition/Haloplex/00100-1407755742_Regions.bed"

if INTERVALS_BED:
    with open(INTERVALS_BED, "r") as fin:
        for line in fin.readlines():
            line = line.rstrip('\n')
            line = re.split(r'\t+', line.rstrip('\t'))
            if len(line) == 4:
                line[1] = int(line[1])
                line[2] = int(line[2])
                bed_line = IntervalColumns(*(line[0], int(line[1]), int(line[2]), line[3]))
                intervals_list.append(bed_line)
            if len(line) == 12:
                line[1] = int(line[1])
                line[2] = int(line[2])
                bed_line = IntervalColumns(*(line[0], int(line[1]), int(line[2]), line[3]))
                intervals_list.append(bed_line)
else:
    print("ERROR: Provide an interval list (bed format)")

directory = '/media/partition/Haloplex/Haloplex_Test_2_Mid_February/Velona'

fig,ax1 = plt.subplots()
plt.hold = True

bam_file_list = [f for f in glob.iglob(directory+"/*.bam")]

collected = defaultdict(list)

boxes = []
collected = defaultdict(list)

for file in bam_file_list:

    coverage_list = []
    almnt = pybedtools.BedTool(file)

    coverage_result = almnt.coverage(intervals_list).sort()

    basename = re.sub('.bam$','',file)
    basename = re.sub('/media/partition/Haloplex/Haloplex_Test_2_Mid_February/Velona/','',basename)

    collected[basename] = {
    'Coverages' : []
    }

    for interval in coverage_result:
        collected[basename]['Coverages'].append(float(interval[4]))

for key, value in collected.items():
    boxes.append(collected[key]['Coverages'])

plt.boxplot(boxes)
xtickNames = plt.setp(ax1, xticklabels = collected.keys())
plt.setp(xtickNames, rotation=90, fontsize=7)
plt.ylabel('Coverage (x)')
plt.title('Comparison of Amplicon Depths per Sample')

plt.show()
