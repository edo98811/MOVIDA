import glob
import re
from pathlib import Path

# --- Detect samples ---
folder = Path("data")
fastq_files = glob.glob(Path("folder", "*_R0[12]*.fastq.gz"))

def sample_name(f):
    m = re.search(r"(.+)_R0[12]", Path(f).stem)
    return m.group(1) if m else None

samples = sorted(set(sample_name(f) for f in fastq_files if sample_name(f)))


rule all:
    input: "samples.csv"


rule make_sample_csv:
    output: "samples.csv"
    run:
        with open(output[0], "w") as f:
            f.write("sample_id,fastq1,fastq2\n")
            for s in samples:
                r1 = sorted(glob.glob(f"data/{s}_R01*.fastq.gz"))
                r2 = sorted(glob.glob(f"data/{s}_R02*.fastq.gz"))
                if r1 and r2:
                    f.write(f"{s},{r1[0]},{r2[0]}\n")
