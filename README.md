# gvcf_norm

**Command-line tool for [left-aligning and normalizing](https://genome.sph.umich.edu/wiki/Variant_Normalization#Algorithm_for_Normalization) gVCF variants**

Comparable to `vt normalize` and `bcftools norm -f`, but (i) tolerates symbolic alleles (`<NON_REF>`, `<*>`, `*`) in variant records and (ii) passes through reference bands unchanged. Because it can reposition variant records *but not* reference bands, a repositioned variant record will probably end up in the middle of an overlapping reference band, with a small coverage gap between the reference bands. But at least the variant record will be normalized.

Build: `cargo build --release` to build `target/release/gvcf_norm`

Usage: `bgzip -dc my.g.vcf.gz | ./gvcf_norm -r ref.fa - | bgzip -c > my.norm.g.vcf.gz`
