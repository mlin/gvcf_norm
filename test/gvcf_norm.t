#!/bin/bash

set -o pipefail

REPO="$(dirname "$0")/.."
cd "$REPO"
export BASH_TAP_ROOT=test/bash-tap
source test/bash-tap/bash-tap-bootstrap

plan tests 10

cargo build
is "$?" "0" "cargo build --release"
gvcf_norm="cargo run -q --release --"

if [[ -z $TMPDIR ]]; then
TMPDIR=/tmp
fi
TMPDIR=$(mktemp -d "${TMPDIR}/gvcf_norm_test_XXXXXX")
export TMPDIR=$(realpath "$TMPDIR")

bgzip -dc test/GRCh38.chr21.fa.gz > "${TMPDIR}/GRCh38.chr21.fa"
samtools faidx "${TMPDIR}/GRCh38.chr21.fa"

cat << EOF > "${TMPDIR}/ex1.g.vcf"
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=RefCall,Description="Genotyping model thinks this site is reference.">
##FILTER=<ID=LowQual,Description="Confidence in this variant being real is below calling threshold.">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position (for use with symbolic alleles)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Conditional genotype quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block.">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele">
##FORMAT=<ID=VAF,Number=A,Type=Float,Description="Variant allele fractions.">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer">
##contig=<ID=chr21,length=46709983>
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA12878
chr21	29848774	.	T	<*>	0	.	END=29848778	GT:GQ:MIN_DP:PL	0/0:50:30:0,99,989
chr21	29848779	.	AT	ATATATTT,T,<*>	40.3	PASS	.	GT:GQ:DP:AD:VAF:PL	1/2:16:32:1,11,12,0:0.34375,0.375,0:40,19,47,19,0,57,990,990,990,990
chr21	29848781	.	T	<*>	0	.	END=29848791	GT:GQ:MIN_DP:PL	0/0:50:19:0,57,569
EOF

$gvcf_norm -r "${TMPDIR}/GRCh38.chr21.fa" "${TMPDIR}/ex1.g.vcf" > "${TMPDIR}/ex1.norm.g.vcf"
is "$?" "0" "ex1 status"

cat "${TMPDIR}/ex1.norm.g.vcf" | tr $'\t' '|' | grep -F 'chr21|29848778|.|TA|TATATATT,T,<*>|40.3|PASS|gvcf_norm_originalPOS=29848779'
is "$?" "0" "ex1 correct"

cat << EOF > "${TMPDIR}/ex2.g.vcf"
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=RefCall,Description="Genotyping model thinks this site is reference.">
##FILTER=<ID=LowQual,Description="Confidence in this variant being real is below calling threshold.">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position (for use with symbolic alleles)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Conditional genotype quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block.">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele">
##FORMAT=<ID=VAF,Number=A,Type=Float,Description="Variant allele fractions.">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer">
##contig=<ID=chr21,length=46709983>
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA12878
chr21	26193733	.	T	<*>	0	.	END=26193733	GT:GQ:MIN_DP:PL	0/0:39:33:0,39,869
chr21	26193734	.	G	T,<*>	29	PASS	.	GT:GQ:DP:AD:VAF:PL	0/1:29:33:21,12,0:0.363636,0:29,0,48,990,990,990
chr21	26193735	.	T	<*>	0	.	END=26193740	GT:GQ:MIN_DP:PL	0/0:50:32:0,69,929
chr21	26193741	.	TTTTTT	T,<*>	29.1	PASS	.	GT:GQ:DP:AD:VAF:PL	0/1:29:32:22,8,0:0.25,0:29,0,55,990,990,990
chr21	26193747	.	T	<*>	0	.	END=26193751	GT:GQ:MIN_DP:PL	0/0:50:22:0,66,659
EOF

bgzip "${TMPDIR}/ex2.g.vcf"
bgzip -dc "${TMPDIR}/ex2.g.vcf.gz" | $gvcf_norm -r "${TMPDIR}/GRCh38.chr21.fa" - > "${TMPDIR}/ex2.norm.g.vcf"
is "$?" "0" "ex2 status"

cat "${TMPDIR}/ex2.norm.g.vcf" | tr $'\t' '|' | grep -F 'chr21|26193734|.|GTTTTT|G,<*>|29.1|PASS|gvcf_norm_originalPOS=26193741|'
is "$?" "0" "ex2 correct"

aria2c -c -d /tmp -s 4 -x 4 --retry-wait 2 \
    ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
is "$?" "0" "download GRCh38"
bgzip -dc /tmp/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz > /tmp/GRCh38.fa
samtools faidx /tmp/GRCh38.fa

gsutil -m cp -n gs://brain-genomics-public/research/cohort/1KGP/dv_vcf/v1/NA12878.dv0.8.0.g.vcf.gz /tmp/
is "$?" "0" "download NA12878.dv0.8.0.g.vcf.gz"

bgzip -dc /tmp/NA12878.dv0.8.0.g.vcf.gz | time $gvcf_norm -r /tmp/GRCh38.fa - > "${TMPDIR}/NA12878.norm.g.vcf"
is "$?" "0" "normalize NA12878 status"

is "$(grep originalPOS "${TMPDIR}/NA12878.norm.g.vcf" | wc -l)" "1443" "normalize NA12878 correct"
grep -v \# "${TMPDIR}/NA12878.norm.g.vcf" | cut -f 1,2 | sort -V -k 1,2 --check
is "$?" "0" "normalize NA12878 sorted"

rm -rf "$TMPDIR"
