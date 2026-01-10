use anyhow::{Context, Result};
use clap::Parser;
use memchr::memmem;
use needletail::parse_fastx_file;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(author, version, about = "Split stitched Hi-C reads by enzyme site", long_about = None)]
struct Args {
    /// Input stitched FASTQ file
    #[arg(short, long)]
    input: PathBuf,

    /// Enzyme cut site sequence (e.g., GATCGATC)
    #[arg(short, long)]
    cutsite: String,

    /// Output file for CUT reads (Interleaved R1/R2)
    #[arg(long)]
    out_cut: PathBuf,

    /// Output file for UNCUT reads
    #[arg(long)]
    out_uncut: PathBuf,
}

fn process_fastq(
    input: &PathBuf,
    cutsite: &str,
    out_cut: &PathBuf,
    out_uncut: &PathBuf,
) -> Result<()> {
    let cutsite_bytes = cutsite.as_bytes();
    let cut_finder = memmem::Finder::new(cutsite_bytes);

    // Open buffers for writing
    let mut cut_writer = BufWriter::new(File::create(out_cut)
        .with_context(|| format!("Could not create output file {:?}", out_cut))?);
    
    let mut uncut_writer = BufWriter::new(File::create(out_uncut)
        .with_context(|| format!("Could not create output file {:?}", out_uncut))?);

    // Zero-copy parser
    let mut reader = parse_fastx_file(input)
        .with_context(|| format!("Could not open input file {:?}", input))?;

    while let Some(record) = reader.next() {
        let seqrec = record.context("Invalid FASTQ record")?;
        let seq = seqrec.seq();

        // Search for the enzyme site
        if let Some(pos) = cut_finder.find(&seq) {
            // --- CUT FOUND: Split into R1 and R2 ---
            
            // Fragment 1 (Left of cutsite)
            let r1_seq = &seq[..pos];
            let r1_qual = seqrec.qual().map(|q| &q[..pos]);

            // Fragment 2 (Right of cutsite)
            let start_r2 = pos + cutsite_bytes.len();
            let r2_seq = &seq[start_r2..];
            let r2_qual = seqrec.qual().map(|q| &q[start_r2..]);

            // Check if fragments are not empty (optional filtering could happen here)
            if !r1_seq.is_empty() && !r2_seq.is_empty() {
                // Write R1 (Fragment 1)
                cut_writer.write_all(b"@")?;
                cut_writer.write_all(seqrec.id())?;
                cut_writer.write_all(b"_1\n")?; // Append _1 suffix
                cut_writer.write_all(r1_seq)?;
                cut_writer.write_all(b"\n+\n")?;
                if let Some(q) = r1_qual { cut_writer.write_all(q)?; }
                cut_writer.write_all(b"\n")?;

                // Write R2 (Fragment 2)
                cut_writer.write_all(b"@")?;
                cut_writer.write_all(seqrec.id())?;
                cut_writer.write_all(b"_2\n")?; // Append _2 suffix
                cut_writer.write_all(r2_seq)?;
                cut_writer.write_all(b"\n+\n")?;
                if let Some(q) = r2_qual { cut_writer.write_all(q)?; }
                cut_writer.write_all(b"\n")?;
            }
        } else {
            // --- NO CUT: Write original record to uncut file ---
            seqrec.write(&mut uncut_writer, None)?;
        }
    }

    Ok(())
}

fn main() -> Result<()> {
    let args = Args::parse();
    process_fastq(&args.input, &args.cutsite, &args.out_cut, &args.out_uncut)
}
#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::io::Write;
    use tempfile::NamedTempFile;

    /// Helper struct to manage temporary files and run assertions cleanly
    struct TestFixture {
        input: NamedTempFile,
        out_cut: NamedTempFile,
        out_uncut: NamedTempFile,
    }

    impl TestFixture {
        fn new(content: &str) -> Self {
            let mut input = NamedTempFile::new().unwrap();
            write!(input, "{}", content).unwrap();
            Self {
                input,
                out_cut: NamedTempFile::new().unwrap(),
                out_uncut: NamedTempFile::new().unwrap(),
            }
        }

        fn run(&self, cutsite: &str) {
            process_fastq(
                &self.input.path().to_path_buf(),
                cutsite,
                &self.out_cut.path().to_path_buf(),
                &self.out_uncut.path().to_path_buf(),
            )
            .unwrap();
        }

        /// Read the "Cut" output file and compare exactly
        fn assert_cut_eq(&self, expected: &str) {
            let content = fs::read_to_string(self.out_cut.path()).expect("Failed to read cut file");
            assert_eq!(content, expected, "Cut output mismatch");
        }

        /// Read the "Uncut" output file and compare exactly
        fn assert_uncut_eq(&self, expected: &str) {
            let content = fs::read_to_string(self.out_uncut.path()).expect("Failed to read uncut file");
            assert_eq!(content, expected, "Uncut output mismatch");
        }
    }

    #[test]
    fn test_standard_cut() {
        // Case: Standard read with site in the middle
        // Read: AAAA (4) + GATC (4) + TTTT (4) = 12 bp
        let fastq = "@read1\nAAAAGATCTTTT\n+\nIIIIIIIIIIII\n";
        let fixture = TestFixture::new(fastq);
        
        fixture.run("GATC");

        // Expect: R1 (AAAA), R2 (TTTT)
        let expected_cut = "\
@read1_1
AAAA
+
IIII
@read1_2
TTTT
+
IIII
";
        fixture.assert_cut_eq(expected_cut);
        fixture.assert_uncut_eq("");
    }

    #[test]
    fn test_no_cut_found() {
        // Case: Read does not contain the site
        let fastq = "@read_uncut\nAAAAGGGGTTTT\n+\nIIIIIIIIIIII\n";
        let fixture = TestFixture::new(fastq);
        
        fixture.run("GATC");

        fixture.assert_cut_eq("");
        fixture.assert_uncut_eq(fastq); // Should match input exactly
    }

    #[test]
    fn test_multiple_cutsites() {
        // Case: Read has TWO cut sites.
        // Logic: Should split at the FIRST occurrence only.
        // Seq: AAAA GATC TTTT GATC CCCC
        let fastq = "@multi\nAAAAGATCTTTTGATCCCCC\n+\n11112222333344445555\n";
        let fixture = TestFixture::new(fastq);

        fixture.run("GATC");

        // Expect: 
        // R1: AAAA (before 1st GATC)
        // R2: TTTTGATCCCCC (everything after 1st GATC)
        let expected_cut = "\
@multi_1
AAAA
+
1111
@multi_2
TTTTGATCCCCC
+
333344445555
";
        fixture.assert_cut_eq(expected_cut);
    }

    #[test]
    fn test_edge_case_cut_at_start() {
        // Case: Cut site is at the very beginning.
        // Seq: GATC TTTT
        // R1 would be empty. Current logic requires !r1.is_empty(), so this should be DROPPED.
        let fastq = "@start_cut\nGATCTTTT\n+\nIIIIIIII\n";
        let fixture = TestFixture::new(fastq);

        fixture.run("GATC");

        // Should be empty in both because it failed the "non-empty fragment" check
        // but was found by finder, so didn't go to uncut.
        fixture.assert_cut_eq(""); 
        fixture.assert_uncut_eq(""); 
    }

    #[test]
    fn test_edge_case_cut_at_end() {
        // Case: Cut site is at the very end.
        // Seq: AAAA GATC
        // R2 would be empty. Should be DROPPED.
        let fastq = "@end_cut\nAAAAGATC\n+\nIIIIIIII\n";
        let fixture = TestFixture::new(fastq);

        fixture.run("GATC");

        fixture.assert_cut_eq("");
        fixture.assert_uncut_eq("");
    }

    #[test]
    fn test_short_reads() {
        // Case: Read is shorter than cut site
        let fastq = "@short\nGAT\n+\nIII\n";
        let fixture = TestFixture::new(fastq);

        fixture.run("GATC");

        // Should go to uncut
        fixture.assert_cut_eq("");
        fixture.assert_uncut_eq(fastq);
    }

    #[test]
    fn test_mixed_batch() {
        // Case: Mixed content in one file
        let fastq = "\
@valid
AAAAGATCTTTT
+
IIIIIIIIIIII
@uncut
AAAAAAAAAAAA
+
IIIIIIIIIIII
@dropped
GATCTTTT
+
IIIIIIII
";
        let fixture = TestFixture::new(fastq);
        fixture.run("GATC");

        let expected_cut = "\
@valid_1
AAAA
+
IIII
@valid_2
TTTT
+
IIII
";
        let expected_uncut = "\
@uncut
AAAAAAAAAAAA
+
IIIIIIIIIIII
";
        fixture.assert_cut_eq(expected_cut);
        fixture.assert_uncut_eq(expected_uncut);
    }
}