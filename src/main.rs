#![feature(iter_map_windows)]
use std::fs::File;
use std::io::{prelude::*, BufReader};

use clap::Parser;
use rust_htslib::bam::{Header, HeaderView, IndexedReader, Read, Reader};

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    #[arg(short, long)]
    bam: Option<String>,

    #[arg(short, long)]
    start: Option<u32>,

    #[arg(short, long)]
    end: Option<u32>,

    #[arg(short, long)]
    chrom: Option<String>,

    #[arg(long)]
    bed: Option<String>,
}

fn main() -> Result<(), &'static str> {
    let args = Args::parse();
    if let Some(bam) = args.bam {
        if let (Some(s), Some(e), Some(c)) = (args.start, args.end, args.chrom) {
            pileup_region(&bam, &c, s, e);
        } else if let Some(bed) = args.bed {
            let file = File::open(bed).unwrap();
            let reader = BufReader::new(file);

            for line in reader.lines() {
                let x = line.unwrap();
                let ll = (x.as_str()).split("\t").collect::<Vec<&str>>();
                let c = ll[0];
                let s = ll[1];
                let e = ll[2];
                pileup_region(
                    &bam,
                    &c,
                    s.parse::<u32>().unwrap(),
                    e.parse::<u32>().unwrap(),
                );
            }
        }
    } else {
        pileup();
    }

    Ok(())
}
struct BasePileupPosition {
    tid: u32,
    pos: u32,
    depth: u32,
    nskips: u32,
}

struct PileupPosition {
    chrom: String,
    pos: u32,
    coverage: u32,
    nskips: u32,
    delta_skip: i64,
}

fn pileup() {
    let mut bam = Reader::from_stdin().unwrap();
    let header = Header::from_template(&bam.header());
    let head_view = HeaderView::from_header(&header);
    let _ = bam
        .pileup()
        .flat_map(|p| {
            let pileup = p.unwrap();
            let mut nskips: u32 = 0;
            for a in pileup.alignments() {
                if a.is_refskip() {
                    nskips += 1;
                }
            }
            Some(BasePileupPosition {
                tid: pileup.tid(),
                pos: pileup.pos(),
                depth: pileup.depth(),
                nskips,
            })
        })
        .map_windows(|[x, y]| PileupPosition {
            chrom: std::str::from_utf8(head_view.tid2name(x.tid))
                .unwrap()
                .to_string(),
            pos: x.pos,
            delta_skip: i64::from(x.nskips) - i64::from(y.nskips),
            coverage: x.depth - x.nskips,
            nskips: x.nskips,
        })
        .for_each(|r| {
            println!(
                "{}\t{}\t{}\t{}\t{}",
                r.chrom, r.pos, r.coverage, r.nskips, r.delta_skip
            );
        });
}

fn pileup_region(bam: &str, chrom: &str, start: u32, end: u32) {
    let mut bam = IndexedReader::from_path(bam).unwrap();
    let _ = bam.fetch((chrom, start, end));
    let header = Header::from_template(&bam.header());
    let head_view = HeaderView::from_header(&header);
    let _ = bam
        .pileup()
        .flat_map(|p| {
            let pileup = p.unwrap();
            if pileup.pos() >= start && pileup.pos() < end {
                let mut nskips: u32 = 0;
                for a in pileup.alignments() {
                    if a.is_refskip() {
                        nskips += 1;
                    }
                }
                Some(BasePileupPosition {
                    tid: pileup.tid(),
                    pos: pileup.pos(),
                    depth: pileup.depth(),
                    nskips,
                })
            } else {
                None
            }
        })
        .map_windows(|[x, y]| PileupPosition {
            chrom: std::str::from_utf8(head_view.tid2name(x.tid))
                .unwrap()
                .to_string(),
            pos: x.pos,
            delta_skip: i64::from(x.nskips) - i64::from(y.nskips),
            coverage: x.depth - x.nskips,
            nskips: x.nskips,
        })
        .for_each(|r| {
            println!(
                "{}\t{}\t{}\t{}\t{}",
                r.chrom, r.pos, r.coverage, r.nskips, r.delta_skip
            );
        });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn simple_bam() {
        // let (lineno, err) = pileup_region("test/volvox.fa", "test/volvox.filtered.vcf");
        // assert_eq!(lineno, 73);
        // assert_eq!(err, 0);
    }
}
