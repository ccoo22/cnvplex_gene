#!/home/genesky/software/python/3.9.4/bin/python3
import re
import sys

from utils import (calculate_intron, calculate_UTR3, calculate_UTR5, parse_gtf,
                   parse_gtf2, parse_mrna_gene)

# 获取所有转录本
mrna_gene = parse_mrna_gene(['info.blue.txt', 'info.green.txt', 'info.yellow.txt'])

# 获取所有mrna结构
mrna_info = parse_gtf2(['/home/pub/output2/research_and_customized_project/cnvplex_gene/GCF_000001405.25_GRCh37.p13_genomic.gtf', '/home/pub/output2/research_and_customized_project/cnvplex_113/GCF_000001405.25_GRCh37.p13_genomic.gtf'])
 

# 提取需要的转录本信息
# 染色体名称转换
chrom_ncbi_to_ucsc = {}
with open('/home/pub/output2/research_and_customized_project/cnvplex_100/GRCh37_latest_assembly_report.txt', 'r') as fh:
    for line in fh:
        if re.match('#', line):
            continue
        values = line.strip().split('\t')
        chrom_ncbi_to_ucsc[values[6]] = values[9]


features = ['exon', 'CDS', 'start_codon', 'stop_codon']
mrna_gene_region = "./mrna_gene_region.txt"

with open(mrna_gene_region, 'w') as fh:
    fh.write("\t".join(['gene', 'mrna', 'gene_id', 'mrna_id', 'if_mrna_match', 'chrom', 'start', 'end', 'width', 'strand', 'feature', 'feature_number']) + "\n")
    # gene/mrna  是苏州提供的基因、转录本名称
    # gene_id/mrna_id  是最终从gtf中提取的基因、转录本名称
    for mrna in mrna_gene.keys():
        gene = mrna_gene[mrna]
        mrna_clean = mrna.split('.')[0]  # 剔除版本号
        
        # 1. 转录本不存在
        if mrna_clean not in mrna_info:
            fh.write("\t".join([gene, mrna, '.', '.', '.', '.', '.', '.', '.', '.', '.', '.']) + "\n")
            continue
        else:
            # 2. 转录本存在
            ## 判断转录本原始版本号是否存在
            mrna_choose = mrna
            if_mrna_match = True
            # 如果版本号对不上，则从中选取最接近的一个
            if mrna_choose not in mrna_info[mrna_clean]:
                for mrna_tmp in sorted(mrna_info[mrna_clean].keys()):
                    if re.match(mrna_clean + "\.\d+$", mrna_tmp):
                        mrna_choose = mrna_tmp
                        if_mrna_match = False
                        break
                if if_mrna_match:
                    print(f"[Error] {mrna} 没有找到 ...\d+ 格式的转录本编号")
                    sys.exit()
            
            for feature in features:
                # gtf_info[mrna][feature].append({'chrom': chrom, 'start': start, 'end': end, 'feature': feature, 'feature_number': feature_number})
                if feature in mrna_info[mrna_clean][mrna_choose]:
                    for region in mrna_info[mrna_clean][mrna_choose][feature]:
                        fh.write("\t".join([gene, mrna, region['gene_id'], region['mrna_id'], str(if_mrna_match), chrom_ncbi_to_ucsc[region['chrom']], str(region['start']), str(region['end']), str(region['end'] - region['start'] + 1), region['strand'] ,region['feature'], region['feature_number']]) + "\n")
            # 计算出intron区域
            introns = calculate_intron(mrna_info[mrna_clean][mrna_choose]['exon'])
            for region in introns:
                fh.write("\t".join([gene, mrna, region['gene_id'], region['mrna_id'], str(if_mrna_match),  chrom_ncbi_to_ucsc[region['chrom']], str(region['start']), str(region['end']), str(region['end'] - region['start'] + 1), region['strand'] ,region['feature'], region['feature_number']]) + "\n")
            # 计算出UTR区域
            utr5s = calculate_UTR5(mrna_info[mrna_clean][mrna_choose]['exon'], mrna_info[mrna_clean][mrna_choose]['start_codon'][0])
            for region in utr5s:
                    fh.write("\t".join([gene, mrna, region['gene_id'], region['mrna_id'], str(if_mrna_match), chrom_ncbi_to_ucsc[region['chrom']], str(region['start']), str(region['end']), str(region['end'] - region['start'] + 1), region['strand'] ,region['feature'], region['feature_number']]) + "\n")
            utr3s = calculate_UTR3(mrna_info[mrna_clean][mrna_choose]['exon'], mrna_info[mrna_clean][mrna_choose]['stop_codon'][0])
            for region in utr3s:
                    fh.write("\t".join([gene, mrna, region['gene_id'], region['mrna_id'], str(if_mrna_match), chrom_ncbi_to_ucsc[region['chrom']], str(region['start']), str(region['end']), str(region['end'] - region['start'] + 1), region['strand'] ,region['feature'], region['feature_number']]) + "\n")
   