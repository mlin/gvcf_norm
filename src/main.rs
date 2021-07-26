use bio::io::fasta;
use clap::Clap;
use std::collections::HashMap;
use std::io::BufRead;
use std::{fs, io};

#[derive(Clap)]
struct Opts {
    // Reference genome uncompressed FASTA
    #[clap(short, long)]
    pub ref_fasta: String,

    // Uncompressed gVCF file/pipe [omit or - for standard input]
    #[clap(default_value = "-")]
    pub gvcf: String,
}

fn main() {
    let opts = Opts::parse();

    if opts.gvcf == "-" && atty::is(atty::Stream::Stdin) {
        panic!("pipe in data or supply input filename")
    }

    let ref_genome = read_ref_fasta(&opts.ref_fasta).unwrap();

    let (_, _, _, mut last_records) = fold_tsv(
        process_gvcf_line,
        (&ref_genome, vec![], String::from(""), vec![]),
        &opts.gvcf,
    )
    .unwrap();
    emit(&mut last_records)
}

fn read_ref_fasta(filename: &str) -> Result<HashMap<String, bio::io::fasta::Record>, io::Error> {
    let mut ans = HashMap::new();
    let mut reader = fasta::Reader::new(fs::File::open(filename)?).records();
    while let Some(Ok(record)) = reader.next() {
        ans.insert(String::from(record.id()), record);
    }
    return Ok(ans);
}

/// Fold over tab-separated lines of the file or standard input.
fn fold_tsv<F, X>(mut f: F, x0: X, filename: &str) -> Result<X, io::Error>
where
    F: FnMut(usize, X, &Vec<&str>) -> X,
{
    // https://stackoverflow.com/a/49964042/13393076
    let reader: Box<dyn io::BufRead> = if filename.is_empty() || filename == "-" {
        Box::new(io::BufReader::new(io::stdin()))
    } else {
        Box::new(io::BufReader::new(fs::File::open(filename)?))
    };

    let mut x = x0;
    let mut line_num = 0;
    for readline in reader.lines() {
        let line = readline?;
        line_num += 1;
        x = f(line_num, x, &line.split('\t').collect())
    }
    Ok(x)
}

fn process_gvcf_line<'a>(
    line_num: usize,
    state: (
        &'a HashMap<String, bio::io::fasta::Record>,
        Vec<String>,
        String,
        Vec<(u64, String)>,
    ),
    fields: &Vec<&str>,
) -> (
    &'a HashMap<String, bio::io::fasta::Record>,
    Vec<String>,
    String,
    Vec<(u64, String)>,
) {
    let (ref_genome, mut header, mut chrom, mut chrom_records) = state;
    if !fields.is_empty() && fields[0].chars().nth(0) == Some('#') {
        if chrom.len() > 0 {
            panic!("gvcf_norm: out-of-place header line {}", line_num)
        }
        header.push(fields.join("\t"));
        return (ref_genome, header, chrom, chrom_records);
    } else if header.is_empty() {
        panic!("gvcf_norm: no header found; check input is uncompressed")
    }
    if fields.len() < 10 || fields[0].len() == 0 {
        panic!("gvcf_norm: malformed input line")
    }
    if fields[0] != chrom {
        if chrom.len() == 0 {
            // emit header
            for header_line in header.iter() {
                println!("{}", header_line)
            }
        }
        emit(&mut chrom_records);
        chrom_records.clear();
        chrom = String::from(fields[0])
    }

    chrom_records.push(normalize_gvcf_record(ref_genome, line_num, fields));
    return (ref_genome, header, chrom, chrom_records);
}

fn emit(records: &mut Vec<(u64, String)>) {
    // sort records by pos, then write them to standard output
    records.sort();
    for (_, line) in records.iter() {
        println!("{}", line)
    }
}

fn normalize_gvcf_record(
    ref_genome: &HashMap<String, bio::io::fasta::Record>,
    line_num: usize,
    fields: &Vec<&str>,
) -> (u64, String) {
    let original_pos = fields[1].parse::<u64>().unwrap();

    // normalize variant records; add INFO field to indicate gvcf_norm_originalPOS=
    // ignore reference bands

    return (original_pos, fields.join("\t"));
}
