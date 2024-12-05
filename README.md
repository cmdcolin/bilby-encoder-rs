## bilby-encoder-rs

Adaptation of https://github.com/ianholmeslab/bilby_encoder to rust

## Usage

```
## Pass BED file of regions
bilby-encoder-rs --bed regions.bed --bam rnaseq.bam > out

## Pass chrom,start,end
bilby-encoder-rs --chrom chr1 --start 1 --end 1000 --bam rnaseq.bam > out

## Pass whole bam to stdin, processes whole file
cat rnaseq.bam | bilby-encoder-rs > out
```

Reports output to stdout, looks like this currently, <-- annotations added for clarity

```
refName	pos	n_match	n_skip	n_match	n_match_to_n_skip
chr1	153357876	8001	0	0 <-- exon
chr1	153357877	8002	0	0
chr1	153357878	8003	0	0
chr1	153357879	8004	0	0
chr1	153357880	8005	0	-8005 <-- exon intron boundary
chr1	153357881	0	8005	0
chr1	153357882	1	8005	0
chr1	153357883	1	8005	0
chr1	153357884	1	8005	0 <-- intron
chr1	153357885	1	8005	0
chr1	153357886	1	8005	0
chr1	153357887	1	8005	0
chr1	153357888	1	8005	0
chr1	153357889	2	8005	0
```


## Changes compared to python bilby_encoder

- Use Rust and rust-htslib for speed
- Use "pileup" command instead of manual CIGAR tracking per read. This allows easy checking of ranges, no sophisticated custom CIGAR tracking. Transitions between two states are identified by counting up the difference between two columns of the pilup e.g. The count of a transition column from M->N can therefore be reporting for example as TOTAL_COVERAGE - N_SKIP at each position (and indeed, N->M can be reported the same way, but instead report a negative value. This allows usage of a single row, reducing number of channels needed)
- No read strandedness reporting right now

## Motivation

- Try to optimize the encoder as much as possible

## Todos

- Perhaps make compatible with current bilby_encoder test suite (would need to restore strand)
- Add more strandedness options that reflect biological strand rather than just read strand  https://github.com/ianholmeslab/bilby_encoder/issues/19
- System for processing multiple BAM files, I anticipate running this program in parallel per BAM file

## Alternatives and other options

- I considered using https://github.com/sstadick/perbase/ which is multithreaded but dos not accept BAM from stdin, and does not report strandedness, and is a bit advanced

- I also considered parsing the text output of `samtools mpileup` directly

- I also considered porting the bilby code as-is to rust

- I also considered using the pileup technique in python instead of rewrite

- I also considered no rewrite, just parallelizing the python script (e.g. running with GNU parallel...don't need to necessarily write any parallel code in python)

- Do we even need the transitions? would have to see at machine learning stage.
  Maybe that transition spike is very helpful, but also it is 'derivable' as basically an edge detection from coverage to n_skip

- Interfacing python and rust could also be considered via pyO3 if it is desirable to have very fast data handoff, but it just writes results to disk here


## Benchmark

```
## Python version on 1000 regions, ENCODE RNA-seq
$ time python src/encoder.py -b subset.bed -f hg19.mod.fa -a ENCFF712DCV.bam
533s

## Rust version on 1000 regions, ENCODE RNA-seq
$ time bilby-encoder-rs --bed subset.bed --bam ENCFF712DCV.bam
34s


## Parallelize by region with ENCODE RNA-seq
time parallel --colsep '\t' bilby-encoder-rs --bam ENCFF712DCV.bam --chrom {1} --start {2} --end {3}  :::: subset.bed
7s
```

Note: It's not an apples to apples comparison since we use pilup instead of CIGAR tracking, and don't write direct numpy arrays, and a couple other things, but it could be a possible avenue for speedup

