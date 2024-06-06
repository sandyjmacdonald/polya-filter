# polya-filter
A simple and fast Python script to filter a fastq file and keep/discard reads with polyA/T stretches, designed for long read (e.g. ONT/Nanopore) transcriptomic reads where polyA/T tails are expected.

# Examples

Filter an input file called `reads.fastq`, filter to **keep** reads with a polyA/T stretch of 20 or more, with a mismatch rate of 0.2, and adaptor length of up to 30, and then output filtered fastq file to `reads.filtered.fastq`:

```bash
polya-filter.py -i reads.fastq -o reads.filtered.fastq -p 20 -m 0.2 -a 30 --keep
```

Filter an input file called `reads.fastq`, filter to **discard** reads with a polyA/T stretch of 12 or more, with a mismatch rate of 0.1, and adaptor length of up to 75, and then output to stdout (good for piping to e.g. seqkit):

```bash
polya-filter.py -i reads.fastq -p 12 -m 0.1 -a 75 --discard
```