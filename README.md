# polya-filter
A simple and fast Python script to filter a fastq file and keep/discard reads with polyA/T stretches, designed for long read (e.g. ONT/Nanopore) transcriptomic reads where polyA/T tails are expected.

Runs *really* quickly! Takes just over 5 seconds to filter a set of 1 million ONT reads in a ~500MB fastq file. This is because of the efficient way that reads are parsed, four lines at a time, and using only the python standard library components.

## Arguments

```
usage: polya-filter.py [-h] -i INPUT [-o OUTPUT] [-p POLYA_LENGTH] [-m MISMATCH] [-a ADAPTOR_LENGTH] [--keep] [--discard]

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input fastq filename
  -o OUTPUT, --output OUTPUT
                        output fastq filename, omit to write to stdout
  -p POLYA_LENGTH, --polya_length POLYA_LENGTH
                        minimum lenth of polyA/T to filter by
  -m MISMATCH, --mismatch MISMATCH
                        proportion of base mismatches within polyA/T (0.0-1.0)
  -a ADAPTOR_LENGTH, --adaptor_length ADAPTOR_LENGTH
                        length of adaptor to pad with before polyA/T
  --keep                Keep reads with polyA/T match
  --discard             Discard reads with polyA/T match
  ```

## Examples

Read an input file called `reads.fastq`, filter to **keep** reads with a polyA/T stretch of 20 or more, with a mismatch rate of 0.2, and adaptor length of up to 30, and then output filtered fastq file to `reads.filtered.fastq`:

```bash
./polya-filter.py -i reads.fastq -o reads.filtered.fastq -p 20 -m 0.2 -a 30 --keep
```

Read an input file called `reads.fastq`, filter to **discard** reads with a polyA/T stretch of 12 or more, with a mismatch rate of 0.1, and adaptor length of up to 75, and then output to stdout (good for piping to e.g. seqkit):

```bash
./polya-filter.py -i reads.fastq -p 12 -m 0.1 -a 75 --discard
```