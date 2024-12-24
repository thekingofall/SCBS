
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import sys

def main():
    """
    用法:
      python extract_methyl_stat.py /path/to/SC3 > result.tsv

    功能:
      1. 在 SC3/Deduplicated 下寻找 *.deduplicated.bam 文件(自动识别样本名)。
      2. 从 Aligned, Deduplicated, Methylation 中找到对应的 3 个 Bismark 报告文件，
         提取所需字段，并按照指定表格格式输出。

    注意:
      - Raw Seq. 和 Seq. For Mapping1 相同(= sequence pairs analysed in total)。
      - 双端测序不再乘以 2。
      - % mCHH/mCHG = (甲基化 CHG) / (甲基化 CHH)。
      - Nb. CHH/G = 甲基化CHH + 未甲基化CHH + 甲基化CHG + 未甲基化CHG。
    """

    if len(sys.argv) < 2:
        print(f"用法: {sys.argv[0]} <SC3目录>", file=sys.stderr)
        sys.exit(1)

    base_dir = sys.argv[1]
    aligned_dir = os.path.join(base_dir, "Aligned")
    dedup_dir   = os.path.join(base_dir, "Deduplicated")
    meth_dir    = os.path.join(base_dir, "Methylation")

    # 输出表头
    print(
        "ID\tSample\tRaw Seq.\tSeq. For Mapping1\tSeq. Mapped\t% Mapped\t% Duplication2\tNb. CpGs\t% of total CpGs\t% mCpG\t% mCpG chr.MT\t% mCHH/mCHG2\tNb. CHH/G"
    )

    if not os.path.isdir(dedup_dir):
        print(f"错误: {dedup_dir} 目录不存在或不可访问", file=sys.stderr)
        sys.exit(1)

    # 找到所有 *.deduplicated.bam 文件
    bam_files = [
        f for f in os.listdir(dedup_dir)
        if f.endswith(".deduplicated.bam")
    ]

    for bam_file in bam_files:
        # 假设文件名形如: "XM-ZX-HiME-4_L3_1_val_1_bismark_bt2_pe.deduplicated.bam"
        # 去掉 ".deduplicated.bam" => 样本 ID
        sample_id = bam_file.replace(".deduplicated.bam", "")

        # 构造对应的 3 个报告文件路径
        # 1) PE_report
        aligned_id = sample_id.replace("_pe", "_PE")
        pe_report = os.path.join(aligned_dir, f"{aligned_id}_report.txt")

        # 2) Dedup report
        dedup_report = os.path.join(dedup_dir, f"{sample_id}.deduplication_report.txt")

        # 3) Splitting report
        # 假设叫: "<sample_id>.deduplicated_splitting_report.txt"
        splitting_report = os.path.join(meth_dir, f"{sample_id}.deduplicated_splitting_report.txt")

        # 初始化变量
        seq_pairs_total = 0     # sequence pairs analysed in total
        seq_mapped = 0
        pct_mapped = 0.0
        pct_dup = 0.0

        total_c = 0             # total number of C's analysed
        meth_cpg = 0
        unmeth_cpg = 0
        pct_mCpG = 0.0

        # CHG
        meth_chg = 0
        unmeth_chg = 0  # 需要解析 "Total C to T conversions in CHG context"
        # CHH
        meth_chh = 0
        unmeth_chh = 0  # 需要解析 "Total C to T conversions in CHH context"

        pct_chg = 0.0
        pct_chh = 0.0   # 可能不再需要，但保留也行

        # ------------------ 读取 PE_report --------------------
        if os.path.exists(pe_report):
            with open(pe_report, "r") as f:
                for line in f:
                    line = line.strip()
                    if line.startswith("Sequence pairs analysed in total:"):
                        seq_pairs_total = int(line.split(":")[1].strip())
                    elif line.startswith("Number of paired-end alignments with a unique best hit:"):
                        seq_mapped = int(line.split(":")[1].strip())
                    elif line.startswith("Mapping efficiency:"):
                        val = line.split(":")[1].strip().replace("%","")
                        pct_mapped = float(val)

        # ----------------- 读取 Dedup report ------------------
        if os.path.exists(dedup_report):
            with open(dedup_report, "r") as f:
                for line in f:
                    line = line.strip()
                    if line.startswith("Total number duplicated alignments removed:"):
                        # 例如: "Total number duplicated alignments removed: 2746545 (44.47%)"
                        match = re.search(r"\(([\d\.]+)%\)", line)
                        if match:
                            pct_dup = float(match.group(1))

        # ---------------- 读取 Splitting report --------------
        if os.path.exists(splitting_report):
            with open(splitting_report, "r") as f:
                for line in f:
                    line = line.strip()
                    if line.startswith("Total number of C's analysed:"):
                        total_c = int(line.split(":")[1].strip())
                    elif line.startswith("Total methylated C's in CpG context:"):
                        meth_cpg = int(line.split(":")[1].strip())
                    elif line.startswith("Total C to T conversions in CpG context:"):
                        unmeth_cpg = int(line.split(":")[1].strip())

                    elif line.startswith("Total methylated C's in CHG context:"):
                        meth_chg = int(line.split(":")[1].strip())
                    elif line.startswith("Total methylated C's in CHH context:"):
                        meth_chh = int(line.split(":")[1].strip())
                    elif line.startswith("Total C to T conversions in CHG context:"):
                        unmeth_chg = int(line.split(":")[1].strip())
                    elif line.startswith("Total C to T conversions in CHH context:"):
                        unmeth_chh = int(line.split(":")[1].strip())

                    elif line.startswith("C methylated in CpG context:"):
                        # 用不到也行, 如果只需最后 % mCpG, 这里就可取
                        pct_mCpG = float(line.split(":")[1].strip().replace("%",""))
                    elif line.startswith("C methylated in CHG context:"):
                        pct_chg = float(line.split(":")[1].strip().replace("%",""))
                    elif line.startswith("C methylated in CHH context:"):
                        pct_chh = float(line.split(":")[1].strip().replace("%",""))

        # --------------- 计算表格各列 ---------------
        # 1) Raw Seq. & Seq. For Mapping = seq_pairs_total
        raw_seq = seq_pairs_total
        seq_for_mapping = seq_pairs_total

        # 2) Seq. Mapped = seq_mapped
        # 3) % Mapped = pct_mapped

        # 4) % Duplication2 = pct_dup

        # 5) Nb. CpGs = meth_cpg + unmeth_cpg
        nb_cpg = meth_cpg + unmeth_cpg

        # 6) % of total CpGs = (Nb. CpGs / total_c)*100
        pct_of_total_cpgs = "-"
        if total_c > 0:
            pct_of_total_cpgs = f"{(nb_cpg / 28000000) * 100:.2f}"

        # 7) % mCpG = pct_mCpG  (之前已从 splitting_report 中获取)
        # 8) % mCpG chr.MT => '-'
        pct_mCpG_chrMT = "-"

        # 9) % mCHH/mCHG2 = (甲基化CHG) / (甲基化CHH)，若 meth_chh=0 则给 "-"
        pct_mch_ratio = "-"
        if meth_chh > 0:
            ratio = meth_chh/meth_chg  
            pct_mch_ratio = f"{ratio:.2f}"

        # 10) Nb. CHH/G = (meth_chh + unmeth_chh + meth_chg + unmeth_chg)
        nb_chhg = meth_chh + unmeth_chh + meth_chg + unmeth_chg

        # 自定义你的 ID 与 Sample
        ID = sample_id
        SAMPLE = "K562"  # 你可根据需要硬编码样本名，或从 sample_id 里解析

        # 输出一行
        print(
            f"{ID}\t"
            f"{SAMPLE}\t"
            f"{raw_seq}\t"
            f"{seq_for_mapping}\t"
            f"{seq_mapped}\t"
            f"{pct_mapped:.1f}\t"
            f"{pct_dup:.1f}\t"
            f"{nb_cpg}\t"
            f"{pct_of_total_cpgs}\t"
            f"{pct_mCpG:.1f}\t"
            f"{pct_mCpG_chrMT}\t"
            f"{pct_mch_ratio}\t"
            f"{nb_chhg}"
        )

if __name__ == "__main__":
    main()

