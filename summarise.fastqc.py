import os
import pandas as pd

trimmed_rv = []
trimmed_fwd = []
reads_rv = []
reads_fwd = []

db = pd.DataFrame()
for subdir, dirs, files in os.walk("data/fastqc/"):
    for file in files:
        filepath = subdir + os.sep + file
        if filepath.endswith("data.txt"):
            stats_file = open (filepath)
            for line in stats_file:
                if "Total Bases" in line:
                    if "trimmed" in filepath:
                        if "_2.trimmed_fastqc" in filepath:
                            trimmed_rv = trimmed_rv + [filepath.replace("data/fastqc/trimmed/","").replace("_2.trimmed_fastqc/fastqc_data.txt","") + "," + line.replace("Total Bases\t","").replace("\n","")]
                        if "_1.trimmed_fastqc" in filepath:
                            trimmed_fwd = trimmed_fwd + [filepath.replace("data/fastqc/trimmed/","").replace("_1.trimmed_fastqc/fastqc_data.txt","") + "," + line.replace("Total Bases\t","").replace("\n","")]
                    if "reads" in filepath:
                        if "_2_fastqc" in filepath:
                            reads_rv = reads_rv + [filepath.replace("data/fastqc/reads/","").replace("_2_fastqc/fastqc_data.txt","") + "," + line.replace("Total Bases\t","").replace("\n","")]
                        if "_1_fastqc" in filepath:
                            reads_fwd = reads_fwd + [filepath.replace("data/fastqc/reads/","").replace("_1_fastqc/fastqc_data.txt","") + "," + line.replace("Total Bases\t","").replace("\n","")]

trimmed_rv_1 = pd.DataFrame(trimmed_rv)
trimmed_fwd_1 = pd.DataFrame(trimmed_fwd)
reads_fwd_1 = pd.DataFrame(reads_fwd)
reads_rv_1 = pd.DataFrame(reads_rv)
trimmed_rv_1 = trimmed_rv_1[0].str.split(',',expand=True)
trimmed_rv_1 = trimmed_rv_1.rename(columns={0: "isolate", 1: "trim_rv"})
trimmed_fwd_1 = trimmed_fwd_1[0].str.split(',',expand=True)
trimmed_fwd_1 = trimmed_fwd_1.rename(columns={0: "isolate", 1: "trimmed_fwd"})
reads_fwd_1 = reads_fwd_1[0].str.split(',',expand=True)
reads_fwd_1 = reads_fwd_1.rename(columns={0: "isolate", 1: "reads_fwd"})
reads_rv_1 = reads_rv_1[0].str.split(',',expand=True)
reads_rv_1 = reads_rv_1.rename(columns={0: "isolate", 1: "reads_rv"})
trimmed_df = pd.merge(trimmed_rv_1, trimmed_fwd_1, on="isolate", how="inner")
reads_df = pd.merge(reads_rv_1, reads_fwd_1, on="isolate", how="inner")
merged_df = pd.merge(trimmed_df, reads_df, on="isolate", how="inner")
merged_df.to_csv("data/stats/fastqc_stats.csv", sep=',', index=False)
