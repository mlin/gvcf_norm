use clap::Clap;
use std::io::BufRead;
use std::{fs, io};

#[derive(Clap)]
struct Opts {
    // Uncompressed [g]VCF file/pipe [omit or - for standard input]
    #[clap(default_value = "-")]
    pub input_filename: String,
}

// read gvcf into vector chromosome-by-chromosome
// normalize variant records; add INFO field to indicate gvcf_norm_originalPOS=
// ignore reference bands
// sort chromosome and emit
fn main() {
    let opts = Opts::parse();

    if opts.input_filename == "-" && atty::is(atty::Stream::Stdin) {
        panic!("pipe in data or supply input filename")
    }

    let (_, _, mut last_records) = fold_tsv(
        process_line,
        (vec![], String::from(""), vec![]),
        &opts.input_filename,
    )
    .unwrap();
    emit(&mut last_records)
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

fn process_line(
    line_num: usize,
    state: (Vec<String>, String, Vec<(u64, String)>),
    fields: &Vec<&str>,
) -> (Vec<String>, String, Vec<(u64, String)>) {
    let (mut header, mut chrom, mut chrom_records) = state;
    if !fields.is_empty() && fields[0].chars().nth(0) == Some('#') {
        if chrom.len() > 0 {
            panic!("gvcf_norm: out-of-place header line {}", line_num)
        }
        header.push(fields.join("\t"));
        return (header, chrom, chrom_records);
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
    let pos = fields[1].parse::<u64>().unwrap();

    // TODO actually realign
    chrom_records.push((pos, fields.join("\t")));
    return (header, chrom, chrom_records);
}

fn emit(records: &mut Vec<(u64, String)>) {
    // sort records by pos, then write them to standard output
    records.sort();
    for (_, line) in records.iter() {
        println!("{}", line)
    }
}
