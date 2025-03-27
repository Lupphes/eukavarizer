import os
from ftplib import FTP
from tqdm import tqdm

# FTP details
FTP_HOST = "ftp.ncbi.nlm.nih.gov"
FTP_DIR = "/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/reads"
FILES_TO_DOWNLOAD = [
    ("D1_S1_L001_R1_001.fastq.gz", "D1_S1_L001_R2_001.fastq.gz"),
    ("D1_S1_L001_R1_002.fastq.gz", "D1_S1_L001_R2_002.fastq.gz"),
    ("D1_S1_L001_R1_003.fastq.gz", "D1_S1_L001_R2_003.fastq.gz")
]

# Output directory structure
OUTPUT_DIR = "./data/9606/sequences"

def download_file(ftp, remote_file, local_file):
    remote_size = ftp.size(remote_file)

    if os.path.exists(local_file):
        local_size = os.path.getsize(local_file)

        if local_size == remote_size:
            print(f"{local_file} already exists and is complete, skipping...")
            return
        else:
            print(f"{local_file} is incomplete. Resuming from {local_size / (1024 * 1024):.2f} MB...")

    print(f"Downloading {remote_file} to {local_file} ({remote_size / (1024 * 1024):.2f} MB)...")

    mode = 'ab' if os.path.exists(local_file) else 'wb'  # Append if file exists
    with open(local_file, mode) as f:
        with tqdm(total=remote_size, initial=os.path.getsize(local_file) if os.path.exists(local_file) else 0,
                    unit='B', unit_scale=True, desc=os.path.basename(local_file)) as pbar:

            def callback(data):
                f.write(data)
                pbar.update(len(data))

            if os.path.exists(local_file):
                # Resume download if file is incomplete
                ftp.sendcmd(f"REST {os.path.getsize(local_file)}")

            ftp.retrbinary(f"RETR {remote_file}", callback)

    print(f"{local_file} downloaded successfully!\n")

def main():
    with FTP(FTP_HOST) as ftp:
        ftp.login()
        ftp.cwd(FTP_DIR)

        for file1, file2 in FILES_TO_DOWNLOAD:
            # Extract base sample name and suffix (like _001, _002)
            base_sample = file1.split("_R1_")[0]
            suffix = file1.split("_R1_")[1].replace(".fastq.gz", "")

            # Create sample name using suffix
            sample_name = f"{base_sample}_{suffix}"

            # Create output directory for each sample
            sample_dir = os.path.join(OUTPUT_DIR, sample_name)
            os.makedirs(sample_dir, exist_ok=True)

            # Define output filenames inside the sample directory
            out_file1 = os.path.join(sample_dir, f"{sample_name}_1.fastq.gz")
            out_file2 = os.path.join(sample_dir, f"{sample_name}_2.fastq.gz")

            # Download files with progress and resume support
            download_file(ftp, file1, out_file1)
            download_file(ftp, file2, out_file2)

    print("\n✅ All files downloaded and organized successfully!")

if __name__ == "__main__":
    main()
