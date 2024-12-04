use clap::Parser;
use rust_htslib::bam::{IndexedReader, Read};

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    #[arg(short, long)]
    bam: String,

    #[arg(short, long)]
    start: u32,

    #[arg(short, long)]
    end: u32,

    #[arg(short, long)]
    chrom: String,
}

fn main() -> Result<(), &'static str> {
    let args = Args::parse();
    let result = verify(&args.bam, &args.chrom, args.start, args.end);
    for entry in result {
        println!(
            "{}\t{}\t{}\t{}\t{}",
            entry.chrom, entry.pos, entry.coverage, entry.nskips, entry.delta_skip
        )
    }
    Ok(())
}
struct BasePileupPosition {
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

fn verify(bam: &str, chrom: &str, start: u32, end: u32) -> Vec<PileupPosition> {
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
        // let (lineno, err) = verify("test/volvox.fa", "test/volvox.filtered.vcf");
        // assert_eq!(lineno, 73);
        // assert_eq!(err, 0);
    }
}
