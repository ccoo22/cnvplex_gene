#!/home/genesky/software/python/3.9.4/bin/python3
import os
import re

import xlsxwriter
from utils import (design_blue, design_green, design_yellow,
                   judge_green_to_blue, match_probe_by_region,
                   match_probe_by_region2, output_blue, output_green,
                   output_readme, output_yellow, parse_probe_info,
                   read_mrna_gene_info, read_probe, workbook_format)

genome = '/home/genesky/database_new/ucsc/fasta/hg19/samtools_index/hg19.fa'

# 读入转录本的结构数据
mrna_gene_info = read_mrna_gene_info('./mrna_gene_region.txt')

# 读取转录本预设的外显子信息
mrna_prepare = parse_probe_info('info.blue.txt', 'info.yellow.txt', 'info.green.fixed.txt')

# 判定转录本是否转换为blue
judge_green_to_blue(mrna_gene_info, mrna_prepare)

# 开始设计
designed_region = {}

# 黄色区域设计
design_yellow(designed_region, mrna_gene_info, mrna_prepare['yellow'])
design_green(designed_region, mrna_gene_info, mrna_prepare['green'])
design_blue(designed_region, mrna_gene_info, mrna_prepare['blue'], genome)


# 读取探针
probe_list = read_probe('./probe')

# 探针匹配
## 1. 对蓝色、黄色区域，直接在区域内匹配探针
match_probe_by_region(designed_region, probe_list)
## 2. 对绿色区域，匹配探针
match_probe_by_region2(designed_region, probe_list)

# 结果输出
excel = './probe/probe_whole_gene.xlsx'
workbook = xlsxwriter.Workbook(excel)
wb_format = workbook_format(workbook)
output_blue(designed_region, workbook, wb_format)
output_yellow(designed_region, workbook, wb_format)
output_green(designed_region, workbook, wb_format)

output_readme(workbook, wb_format)

workbook.close()

