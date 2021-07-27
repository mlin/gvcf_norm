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
            for i in 0..header.len() {
                if i == header.len() - 1 {
                    println!("##INFO=<ID=gvcf_norm_originalPOS,Number=1,Type=Integer,Description=\"POS before gvcf_norm left-aligned the variant\">");
                }
                println!("{}", header[i]);
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
    // parse gVCF fields
    let chrom = fields[0];
    let chrom_seq = ref_genome
        .get(chrom)
        .expect(&format!("line {} unknown CHROM: {}", line_num, chrom))
        .seq();
    let original_pos = fields[1].parse::<usize>().unwrap() - 1;
    let mut alleles = vec![Vec::from(fields[3].as_bytes())];

    if alleles[0].is_empty() {
        panic!("line {} invalid REF", line_num);
    }
    if original_pos + alleles[0].len() > chrom_seq.len() {
        panic!("line {} POS & REF beyond end of CHROM {}", line_num, chrom);
    }
    let ref_seq = &chrom_seq[original_pos..(original_pos + alleles[0].len())];
    if alleles[0] != ref_seq {
        panic!(
            "line {} {}:{} REF {} inconsistent with reference {}",
            line_num,
            chrom,
            original_pos + 1,
            std::str::from_utf8(&alleles[0]).unwrap(),
            std::str::from_utf8(ref_seq).unwrap()
        );
    }

    let mut real_alt = false;
    for alt in fields[4].split(",") {
        if alt.is_empty() {
            panic!("line {} invalid ALT", line_num);
        }
        let altb = Vec::from(alt.as_bytes())[0] as char;
        if altb != '<' && altb != '*' && altb != '.' {
            real_alt = true;
        }
        alleles.push(Vec::from(alt.as_bytes()));
    }

    if !real_alt {
        // pass through reference band (only symbolic ALT alleles)
        return (original_pos as u64, fields.join("\t"));
    }

    // perform normalization
    let mut pos = original_pos as u64;
    if !normalize_gvcf_alleles(chrom_seq, &mut pos, &mut alleles) {
        return (original_pos as u64, fields.join("\t"));
    }

    // generate normalized gVCF record
    let mut fields2 = vec![chrom]; // CHROM

    let postxt = format!("{}", pos + 1); //POS
    fields2.push(&postxt);

    fields2.push(fields[2]); // ID

    fields2.push(std::str::from_utf8(&alleles[0]).unwrap()); // REF

    // ALT
    let mut alts = Vec::new();
    for i in 1..alleles.len() {
        alts.push(std::str::from_utf8(&alleles[i]).unwrap());
    }
    let altstxt = alts.join(",");
    fields2.push(&altstxt);

    // QUAL
    fields2.push(fields[5]);
    // FILTER
    fields2.push(fields[6]);

    // INFO
    let mut infotxt = format!("gvcf_norm_originalPOS={}", original_pos + 1);
    if fields[7] != "." {
        infotxt += ";";
        infotxt += fields[7];
    }
    fields2.push(&infotxt);

    // remaining
    for i in 8..fields.len() {
        fields2.push(fields[i]);
    }

    return (pos as u64, fields2.join("\t"));
}

fn normalize_gvcf_alleles(chrom_seq: &[u8], pos: &mut u64, alleles: &mut Vec<Vec<u8>>) -> bool {
    // Multiallelic normalization algorithm, ref:
    // https://genome.sph.umich.edu/wiki/Variant_Normalization
    let mut changed = false;
    loop {
        // if alleles end with the same nucleotide...
        let end_nuc = alleles[0][alleles[0].len() - 1];
        let mut all_end_nuc = true;
        for al in alleles.iter() {
            if al[0] != '<' as u8 && al[0] != '*' as u8 && al[al.len() - 1] != end_nuc {
                all_end_nuc = false;
            }
        }
        if all_end_nuc {
            // ...then truncate rightmost nucleotide of each allele
            let mut any_empty = false;
            for i in 0..alleles.len() {
                if alleles[i][0] != '<' as u8 && alleles[i][0] != '*' as u8 {
                    alleles[i] = Vec::from(&alleles[i][0..alleles[i].len() - 1]);
                    if alleles[i].is_empty() {
                        any_empty = true;
                    }
                }
            }
            // if there exists an empty allele...
            if any_empty {
                // ...then extend alleles 1 nucleotide to the left
                assert!(*pos > 0);
                *pos -= 1;
                let left_nuc = chrom_seq[*pos as usize];
                if match left_nuc as char {
                    'A' | 'G' | 'C' | 'T' => false,
                    _ => true,
                } {
                    // (abort if not a nucleotide)
                    return false;
                }
                for i in 0..alleles.len() {
                    if alleles[i].is_empty()
                        || (alleles[i][0] != '<' as u8 && alleles[i][0] != '*' as u8)
                    {
                        alleles[i].insert(0, left_nuc);
                    }
                }
            }
            changed = true;
        }
        // repeat until convergence
        if !all_end_nuc {
            break;
        }
    }

    // while leftmost nucleotide of each allele are the same and all alleles have length 2 or more
    loop {
        let left_nuc = if alleles[0].len() >= 2 {
            alleles[0][0]
        } else {
            break;
        };
        let mut left_padded = true;
        for al in alleles.iter() {
            if al[0] != '<' as u8 && al[0] != '*' as u8 && (al.len() < 2 || al[0] != left_nuc) {
                left_padded = false;
                break;
            }
        }
        if !left_padded {
            break;
        }
        // truncate leftmost nucleotide of each allele
        for i in 0..alleles.len() {
            if alleles[i][0] != '<' as u8 && alleles[i][0] != '*' as u8 {
                alleles[i].remove(0);
                assert_eq!(alleles[i].is_empty(), false);
            }
        }
        *pos += 1;
        changed = true;
    }

    return changed;
}
