
import pandas as pd
from os.path import join

# Parse cisbp database for human chipseq experiment
input_dir = "/Users/suryachhetri/Dropbox/for_genemodels/revision_2"
input_file = join(input_dir, "TF_Information_all_motifs.txt")
df = pd.read_csv(input_file, sep="\t")
chip_homo = df[(df.loc[:,"MSource_Type"]=="ChIP-seq") & (df.loc[:,"TF_Species"]=="Homo_sapiens")].iloc[:,[6,10,13,14,15,16]]

# chip_homo = df[(df.loc[:,"MSource_Type"]=="ChIP-seq") & (df.loc[:,"TF_Species"]=="Homo_sapiens")].iloc[:,[6,10,13,14,15,16,17,18]]
all_chip = chip_homo.copy()
all_chip["TF_Name"].nunique() # 131 TF
chip_reset = all_chip.reset_index(drop=True)
chip_reset.to_csv(join(input_dir, "CisBP_chipseq_human.txt"), sep="\t", header=True, index=True)

hepG2_chip = chip_reset.loc[chip_reset["DBID.1"].str.contains("HepG2")].reset_index(drop=True)
hepG2_chip.to_csv(join(input_dir, "CisBP_HepG2_chipseq_human.txt"), sep="\t", header=True, index=True)
hepG2_chip["TF_Name"].nunique() # 45 TF
