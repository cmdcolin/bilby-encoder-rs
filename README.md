## bilby_encoder

Adaptation of https://github.com/ianholmeslab/bilby_encoder


## Changes compared to python bilby_encoder

- Use Rust and rust-htslib for speed
- Use "pileup" command instead of manual CIGAR tracking per read. This allows easy checking of ranges, no sophisticated custom CIGAR tracking. Transitions between two states are identified by counting up the difference between e.g. The count of a transition column from M->N can therefore be reporting for example as TOTAL_COVERAGE - N_SKIP at each position (and indeed, N->M can be reported the same way, but instead report a negative value. This allows usage of a single row, reducing number of channels needed)
- No read strandedness reporting right now

## Motivation

- Try to optimize the encoder as much as possible

## Todos

- Perhaps make compatible with current bilby_encoder test suite (would need to restore strand)
- Add more strandedness options that reflect biological strand rather than just read strand  https://github.com/ianholmeslab/bilby_encoder/issues/19

## Alternatives considered

- I considered using https://github.com/sstadick/perbase/ (which although multithreaded, does not accept BAM from stdin, and does not report strandedness)

- I also considered parsing the text output of `samtools mpileup` directly

- I also considered porting the bilby code as is to rust, or using the pileup
technique in python, but i'm not actually that good at python and slightly
better a rust
