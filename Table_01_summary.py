#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
示例：从 SC3 文件夹的 Aligned、Deduplicated、Methylation 子文件夹中
自动解析 Bismark 的三个报告文件，提取关键信息并输出。
"""

import os
import re
import argparse

def parse_args():
    parser = argparse.ArgumentParser(
        description="从 Bismark 的对齐、去重、甲基化拆分报告文件中提取统计信息并生成表格。"
    )
    parser.add_argument(
        "-i", "--input_folder",
        required=True,
        help="包含 Aligned、Deduplicated、Methylation 三个子目录的根目录。例如：SC3"
    )
    parser.add_argument(
        "-o", "--output_file",
        help="输出结果的文件路径；若不指定，则直接打印到屏幕"
    )
    return parser.parse_args()

def main():
    args = parse_args()

    base_dir = args.input_folder
    aligned_dir = os.path.join(base_dir, "Aligned")
    dedup_dir = os.path.join(base_dir, "Deduplicated")
    meth_dir = os.path.join(base_dir, "Methylation")

    # 准备输出
    output_lines = []
    # 先输出表头
    header = (
        "ID\tSample\tRaw Seq.\tSeq. For Mapping1\tSeq. Mapped\t% Mapped\t"
        "% Duplication\tNb. CpGs\t% of total CpGs\t% mCpG\t"
        "% mCpG chr.MT\t% mCHH/mCHG\tNb. CHH/G"
    )
    output_lines.append(header)

    # 扫描 Deduplicated 文件夹，找出所有 '.deduplicated.bam' 文件
    # 这类文件名形如: XM-ZX-HiME-4_L3_1_val_1_bismark_bt2_pe.deduplicated.bam
    dedup_files = [
        f for f in os.listdir(dedup_dir) if f.endswith(".deduplicated.bam")
    ]

    # 对每个样本 bam 文件进行处理
    for bam_file in dedup_files:
        # 根据文件名提取样本 ID
        #   形如: "XM-ZX-HiME-4_L3_1_val_1_bismark_bt2_pe"
        #   只要匹配前缀 (.+?)_bismark_bt2_pe 即可
        sample_match = re.match(r"(.+?)_bismark_bt2_pe", bam_file)
        if not sample_match:
            continue
        sample_id = sample_match.group(1)

        # 组装对应报告文件的完整路径：
        #   1) Aligned/XM-ZX-HiME-4_L3_1_val_1_bismark_bt2_PE_report.txt
        #   2) Deduplicated/XM-ZX-HiME-4_L3_1_val_1_bismark_bt2_pe.deduplication_report.txt
        #   3) Methylation/XM-ZX-HiME-4_L3_1_val_1_bismark_bt2_pe.deduplicated_splitting_report.txt
        pe_report = os.path.join(aligned_dir, f"{sample_id}_bismark_bt2_PE_report.txt")
        dedup_report = os.path.join(dedup_dir, f"{sample_id}_bismark_bt2_pe.deduplication_report.txt")
        splitting_report = os.path.join(meth_dir, f"{sample_id}_bismark_bt2_pe.deduplicated_splitting_report.txt")

        # 初始化提取信息的变量
        sequence_pairs_total = 0      # 读取自 PE_report
        sequence_pairs_mapped = 0     # 读取自 PE_report
        pct_mapped = 0.0              # 读取自 PE_report
        pct_duplication = 0.0         # 读取自 Deduplication report
        m_cpg = 0                     # 读取自 splitting_report
        u_cpg = 0                     # 读取自 splitting_report
        pct_mCpG = 0.0                # 读取自 splitting_report
        pct_chg = 0.0
        pct_chh = 0.0
        m_chg = 0
        m_chh = 0

        # ============ 1) 解析 PE_report 文件 ============
        if os.path.exists(pe_report):
            with open(pe_report, "r") as f:
                for line in f:
                    line = line.strip()
       
                    if line.startswith("Sequence pairs analysed in total:"):
                        sequence_pairs_total = int(line.split(":")[1].strip())

                    elif line.startswith("Number of paired-end alignments with a unique best hit:"):
                        sequence_pairs_mapped = int(line.split(":")[1].strip())
            
                    elif line.startswith("Mapping efficiency:"):
                        pct_mapped = float(line.split(":")[1].strip().replace("%", ""))

        # ============ 2) 解析 Deduplication_report 文件 ============
        if os.path.exists(dedup_report):
            with open(dedup_report, "r") as f:
                for line in f:
                    line = line.strip()
                    # Total number duplicated alignments removed:	2746545 (44.47%)
                    if line.startswith("Total number duplicated alignments removed:"):
                        match = re.search(r"\(([\d\.]+)%\)", line)
                        if match:
                            pct_duplication = float(match.group(1))

        # ============ 3) 解析 splitting_report 文件 ============
        if os.path.exists(splitting_report):
            with open(splitting_report, "r") as f:
                for line in f:
                    line = line.strip()
           
                    if line.startswith("Total methylated C's in CpG context:"):
                        m_cpg = int(line.split(":")[1].strip())
          
                    elif line.startswith("Total C to T conversions in CpG context:"):
                        u_cpg = int(line.split(":")[1].strip())
             
                    elif line.startswith("C methylated in CpG context:"):
                        pct_mCpG = float(line.split(":")[1].strip().replace("%", ""))
           
                    elif line.startswith("C methylated in CHG context:"):
                        pct_chg = float(line.split(":")[1].strip().replace("%", ""))
 
                    elif line.startswith("C methylated in CHH context:"):
                        pct_chh = float(line.split(":")[1].strip().replace("%", ""))

                    elif line.startswith("Total methylated C's in CHG context:"):
                        m_chg = int(line.split(":")[1].strip())
                    # Total methylated C's in CHH context:	3340348
                    elif line.startswith("Total methylated C's in CHH context:"):
                        m_chh = int(line.split(":")[1].strip())

        # ============ 衍生值计算 ============
        # 假设 Raw Seq = "sequence_pairs_total × 2" (因为是 PE reads)
        raw_seq = sequence_pairs_total * 2

        # 假设 "Nb. CpGs" = m_cpg + u_cpg
        nb_cpgs = m_cpg + u_cpg

        # % mCHH/mCHG, 这里示例做个简单平均
        pct_mCH = (pct_chg + pct_chh) / 2 if (pct_chg and pct_chh) else 0.0

        # Nb. CHH/G, 这里示例为 (m_chg + m_chh)，
        # 如果需要把未甲基化CHH/CHG也算进去，需要再结合 splitting_report 里的 conversions。
        nb_CHHG = m_chg + m_chh

        # 组合输出
        # 注意："% of total CpGs"、"% mCpG chr.MT" 在本例中用 "-" 占位
        result_line = (
            f"{sample_id}\t{sample_id}_oocyte\t"
            f"{raw_seq}\t{sequence_pairs_total}\t{sequence_pairs_mapped}\t"
            f"{pct_mapped:.1f}\t{pct_duplication:.1f}\t"
            f"{nb_cpgs}\t-\t{pct_mCpG:.1f}\t-\t{pct_mCH:.1f}\t{nb_CHHG}"
        )
        output_lines.append(result_line)

    # ============ 输出到文件或屏幕 ============
    if args.output_file:
        with open(args.output_file, "w") as fo:
            fo.write("\n".join(output_lines))
            fo.write("\n")
    else:
        print("\n".join(output_lines))

if __name__ == "__main__":
    main()
