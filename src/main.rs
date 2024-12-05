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
}

fn main() -> Result<(), &'static str> {
    let args = Args::parse();
    let result;
    if let (Some(s), Some(e), Some(c), Some(b)) = (args.start, args.end, args.chrom, args.bam) {
        result = pileup_region(&b, &c, s, e);
    } else {
        result = pileup();
    }
    for entry in result {
        println!(
            "{}\t{}\t{}\t{}\t{}",
            entry.chrom, entry.pos, entry.coverage, entry.nskips, entry.delta_skip
        )
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

fn pileup() -> Vec<PileupPosition> {
    let mut bam = Reader::from_stdin().unwrap();
    let header = Header::from_template(&bam.header());
    let head_view = HeaderView::from_header(&header);
    let rres: Vec<BasePileupPosition> = bam
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
        .collect();

    let result = rres
        .iter()
        .zip(rres.iter().skip(1))
        .map(|(a, b)| PileupPosition {
            chrom: std::str::from_utf8(head_view.tid2name(a.tid))
                .unwrap()
                .to_string(),
            pos: a.pos,
            delta_skip: i64::from(a.nskips) - i64::from(b.nskips),
            coverage: a.depth - a.nskips,
            nskips: a.nskips,
        })
        .collect::<Vec<_>>();

    result
}

fn pileup_region(bam: &str, chrom: &str, start: u32, end: u32) -> Vec<PileupPosition> {
    let mut bam = IndexedReader::from_path(bam).unwrap();
    let _ = bam.fetch((chrom, start, end));
    let rres: Vec<BasePileupPosition> = bam
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
        .collect();

    let result = rres
        .iter()
        .zip(rres.iter().skip(1))
        .map(|(a, b)| PileupPosition {
            chrom: chrom.to_string(),
            pos: a.pos,
            delta_skip: i64::from(a.nskips) - i64::from(b.nskips),
            coverage: a.depth - a.nskips,
            nskips: a.nskips,
        })
        .collect::<Vec<_>>();

    result
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
