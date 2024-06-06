# polya-filter
A simple and fast Python script to filter a fastq file and keep/discard reads with polyA/T stretches, designed for long read (e.g. ONT/Nanopore) transcriptomic reads where polyA/T tails are expected.

**Note:** requires a recent version of Python 3.

Runs *really* quickly! Takes less than 20 seconds to filter a set of 1 million ONT reads in a >2GB fastq file. This is because of the efficient way that reads are parsed, four lines at a time, and using only the Python standard library components.

![](https://media1.tenor.com/m/GNPs4yC-wYgAAAAd/44.gif)

## Arguments

```
usage: polya-filter.py [-h] -i INPUT [-o OUTPUT] [-p POLYA_LENGTH] [-m MISMATCH] [-a ADAPTOR_LENGTH] [--keep] [--discard] [--stats] [--histo]

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
  --keep                Keep reads with polyA/T match (default)
  --discard             Discard reads with polyA/T match
  --stats               Write stats to text file
  --histo               Write histogram of polyA/T lengths to csv file
  ```

## Examples

Read an input file called `reads.fastq`, filter to **keep** reads with a polyA/T stretch of 20 or more, with a mismatch rate of 0.2, and adaptor length of up to 30, and then output the filtered fastq file to `reads.filtered.fastq`, and write stats to `reads.filtered.stats`:

```bash
./polya-filter.py -i reads.fastq -o reads.filtered.fastq -p 20 -m 0.2 -a 30 --keep --stats
```

Read an input file called `reads.fastq`, filter to **discard** reads with a polyA/T stretch of 12 or more, with a mismatch rate of 0.1, and adaptor length of up to 75, and then output to stdout (good for piping to e.g. seqkit):

```bash
./polya-filter.py -i reads.fastq -p 12 -m 0.1 -a 75 --discard
```

The program can also estimate polyA/T lengths at the ends (i.e. 5' and 3') of reads and output a csv-formatted histogram of the lengths and counts, which can then be plotted in your favourite plotting program. Simply add `--histo` to do this, as follows:

```bash
./polya-filter.py -i reads.fastq -o reads.filtered.fastq -p 20 -m 0.2 -a 30 --keep --histo
```