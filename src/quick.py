import os
from ftplib import FTP
from tqdm import tqdm

FTP_HOST = "ftp-trace.ncbi.nlm.nih.gov"
OUTPUT_DIR = "out"

SEQUENCING_DATA = {
    "Illumina": {
        "ftp_dir": "/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/reads",
        "extensions": ["_R1_", "_R2_"],
        "max_pairs": 0,  # set to 0 to skip
        "mode": "paired"
    },
    "PacBio": {
        "ftp_dir": "/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/",
        "extensions": [".Q20.fastq"],
        "max_files": 2,
        "mode": "single"
    },
    "Nanopore": {
        "ftp_dir": "/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/UCSC_Ultralong_OxfordNanopore_Promethion/",
        "extensions": [".fastq.gz"],
        "max_files": 0,  # set to 0 to skip
        "mode": "single"
    }
}

# Track downloaded files
download_summary = []

def format_size(size_bytes):
    if size_bytes >= 1024 ** 3:
        return f"{size_bytes / (1024 ** 3):.2f} GB"
    elif size_bytes >= 1024 ** 2:
        return f"{size_bytes / (1024 ** 2):.2f} MB"
    elif size_bytes >= 1024:
        return f"{size_bytes / 1024:.2f} KB"
    else:
        return f"{size_bytes} B"


def download_file(ftp, remote_file, local_file):
    if os.path.exists(local_file):
        print(f"✅ {local_file} already exists, skipping download.\n")
        download_summary.append({
            'filename': local_file,
            'size': os.path.getsize(local_file)
        })
        return

    try:
        print(f"⬇️  Downloading {remote_file} to {local_file}...")

        with open(local_file, 'wb') as f:
            pbar = tqdm(unit='B', unit_scale=True, desc=os.path.basename(local_file))

            def callback(data):
                f.write(data)
                pbar.update(len(data))

            ftp.retrbinary(f"RETR {remote_file}", callback)
            pbar.close()
            print(f"✅ {local_file} downloaded successfully.\n")

        download_summary.append({
            'filename': local_file,
            'size': os.path.getsize(local_file)
        })

    except Exception as e:
        print(f"❌ Error downloading {remote_file}: {e}")


def list_read_pairs(files, max_pairs):
    r1_files = sorted([f for f in files if "_R1_" in f])
    r2_files = sorted([f for f in files if "_R2_" in f])
    pairs = []

    for r1 in r1_files:
        r2 = r1.replace("_R1_", "_R2_")
        if r2 in r2_files:
            pairs.append((r1, r2))
            if len(pairs) >= max_pairs:
                break
    return pairs


def list_single_files(files, extensions, max_files):
    matching = sorted([f for f in files if any(f.endswith(ext) for ext in extensions)])
    return matching[:max_files]


def main():
    with FTP(FTP_HOST) as ftp:
        ftp.login()

        for platform, config in SEQUENCING_DATA.items():
            print(f"\n📁 Processing {platform} data...")

            if config["mode"] == "paired" and config.get("max_pairs", 0) == 0:
                print("⚠️  Skipping — max_pairs = 0")
                continue
            if config["mode"] == "single" and config.get("max_files", 0) == 0:
                print("⚠️  Skipping — max_files = 0")
                continue

            ftp.cwd(config["ftp_dir"])
            files = ftp.nlst()

            if config["mode"] == "paired":
                read_pairs = list_read_pairs(files, config["max_pairs"])
                for r1, r2 in read_pairs:
                    suffix = r1.split("_R1_")[1].replace(".fastq.gz", "")
                    sample_name = r1.split("_R1_")[0] + "_" + suffix
                    sample_dir = os.path.join(OUTPUT_DIR, sample_name)
                    os.makedirs(sample_dir, exist_ok=True)

                    local_r1 = os.path.join(sample_dir, f"{sample_name}_1.fastq.gz")
                    local_r2 = os.path.join(sample_dir, f"{sample_name}_2.fastq.gz")

                    download_file(ftp, r1, local_r1)
                    download_file(ftp, r2, local_r2)

            elif config["mode"] == "single":
                matching_files = list_single_files(files, config["extensions"], config["max_files"])
                for filename in matching_files:
                    sample_name = os.path.splitext(filename)[0]
                    sample_dir = os.path.join(OUTPUT_DIR, sample_name)
                    os.makedirs(sample_dir, exist_ok=True)

                    local_path = os.path.join(sample_dir, filename)
                    download_file(ftp, filename, local_path)

    # Final Report
    print("\n📄 Final Download Report")
    print("=" * 50)

    total_size = sum(f['size'] for f in download_summary)
    total_files = len(download_summary)

    for file in download_summary:
        print(f"📦 {file['filename']} ({format_size(file['size'])})")

    print("=" * 50)
    print(f"✅ Total files downloaded: {total_files}")
    print(f"📁 Total size: {format_size(total_size)}")
    print("\n🎉 All selected files downloaded successfully!")


if __name__ == "__main__":
    main()
