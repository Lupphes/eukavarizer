import os
import gzip

from .refseq_retriver import RefSeqRetriver
from .ena_searcher import ENASearcher

class Pipeline:
    """
    A main pipeline class to run all steps:
    1. Fetch reference genome from RefSeq.
    2. Search ENA for sequencing runs.
    3. Download FASTQ files.
    """

    def __init__(self, taxonomy_id: int, outdir: str, mode: str, genome_file: str = None,
                 sequences_dir: str = None, genome_size_ungapped: int = None, **kwargs):
        self.taxonomy_id = taxonomy_id
        self.outdir = os.path.abspath(outdir)
        self.mode = mode
        self.genome_file = genome_file
        self.sequences_dir = sequences_dir
        self.genome_size_ungapped = genome_size_ungapped

        os.makedirs(self.outdir, exist_ok=True)

        self.refseq_retriver = RefSeqRetriver(taxonomy_id=self.taxonomy_id, outdir=self.outdir)

        self.ena_searcher = ENASearcher(
            taxonomy_id=self.taxonomy_id,
            library_strategy=kwargs.get("library_strategy"),
            instrument_platform=kwargs.get("instrument_platform"),
            max_results=kwargs.get("max_results", 10),
            min_coverage=kwargs.get("minimum_coverage"),
            max_coverage=kwargs.get("maximum_coverage"),
            assembly_quality=kwargs.get("assembly_quality"),
            sort=True
        )


    def run(self):
        """Runs the appropriate pipeline based on the selected mode."""
        if self.mode == "refseq":
            return self.run_refseq()
        elif self.mode == "ena":
            return self.run_ena()
        elif self.mode == "both":
            refseq_results = self.run_refseq()
            self.genome_size = refseq_results.get("genome_size", self.genome_size)
            self.genome_size_ungapped = refseq_results.get("genome_size_ungapped", self.genome_size_ungapped)
            self.genome_file = refseq_results.get("genome_file", self.genome_file)
            self.run_ena()
        else:
            raise ValueError(f"Unsupported mode: {self.mode}")

    def run_refseq(self) -> dict[str, str]:
        """Fetches the reference genome from RefSeq."""
        if self.genome_file:
            print(f"Using local reference genome: {self.genome_file}")
            self.genome_size, self.genome_size_ungapped = self.calculate_genome_size_from_file(self.genome_file)
        else:
            print("Fetching RefSeq genome data...")
            self.genome_file, self.genome_size, self.genome_size_ungapped = self.refseq_retriver.get_refseq_genomes(self.taxonomy_id)

        print(f"Genome size: {self.genome_size} bp")
        print(f"Ungapped genome size: {self.genome_size_ungapped} bp")
        print(f"RefSeq genome saved to: {self.genome_file}")

        return {"genome_file": self.genome_file, "genome_size": self.genome_size, "genome_size_ungapped": self.genome_size_ungapped}

    def run_ena(self):
        """Executes the ENA search and FASTQ file download pipeline."""
        if self.sequences_dir:
            print(f"Using local sequences from directory: {self.sequences_dir}")
            files = self.list_sequence_files(self.sequences_dir)
            return {"sequences_dir": self.sequences_dir, "sequence_files": files}
        else:
            self.sequences_dir = os.path.join(self.outdir, str(self.taxonomy_id), "sequences")

        print("Searching for sequence data in ENA...")
        if not self.genome_size_ungapped and self.genome_file:
            self.genome_size, self.genome_size_ungapped = self.calculate_genome_size_from_file(self.genome_file)
        elif self.genome_size_ungapped:
            print(f"Using provided ungapped genome size: {self.genome_size_ungapped}")
        else:
            print("Genome size not provided. Please provide the ungapped genome size for ENA search for better results.")

        sequence_data = self.ena_searcher.search_sequence_data(genome_size_ungapped=self.genome_size_ungapped)

        if not sequence_data:
            print("No sequence data found.")
            return {"sequences_dir": self.sequences_dir, "sequence_data": []}

        print("Downloading FASTQ files...")

        self.ena_searcher.fetch_fastq_files(sequence_data, self.sequences_dir)

        files = self.list_sequence_files(self.sequences_dir)
        return {"sequences_dir": self.sequences_dir, "sequence_files": files}

    def list_sequence_files(self, directory: str, print_files: bool = True) -> list[str]:
        """Lists all sequence files in the specified directory and returns them as a list."""
        if os.path.exists(directory):
            files = os.listdir(directory)
            if files:
                if print_files:
                    print("\nSequence files in the directory:")
                    print("\n".join(f"- {file}" for file in files))
                return files
            else:
                print("No sequence files found in the directory.")
                return []
        else:
            print(f"Directory {directory} does not exist.")
            return []

    def calculate_genome_size_from_file(self, genome_file: str) -> tuple[int, int]:
        """Calculates genome size by summing sequence lengths from a genome file (handles both .gz and plain text)."""
        genome_size = 0
        genome_size_ungapped = 0

        # Determine if file is gzipped
        open_func = gzip.open if genome_file.endswith(".gz") else open

        with open_func(genome_file, 'rt', encoding="utf-8") as f:  # 'rt' ensures reading text mode
            for line in f:
                if not line.startswith('>'):
                    sequence = line.strip()
                    genome_size += len(sequence)
                    genome_size_ungapped += len(sequence.replace('N', ''))

        return genome_size, genome_size_ungapped
