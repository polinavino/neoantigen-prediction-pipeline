#!/usr/bin/env python3
"""Convert a GDC MAF file to VCF format for pVACseq input."""
import sys
import csv

def maf_to_vcf(maf_file, vcf_file, sample_name):
    with open(maf_file, 'r') as maf_fh:
        # Skip comment lines
        lines = [l for l in maf_fh if not l.startswith('#')]
    
    reader = csv.DictReader(lines, delimiter='\t')
    
    with open(vcf_file, 'w') as vcf_fh:
        # Write VCF header
        vcf_fh.write('##fileformat=VCFv4.2\n')
        vcf_fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        vcf_fh.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n')
        vcf_fh.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n')
        vcf_fh.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\t{}\n'.format(sample_name))
        
        for row in reader:
            # Only process SNPs and indels
            if row.get('Variant_Classification') in [
                'Missense_Mutation', 'Nonsense_Mutation', 'Silent',
                'Splice_Site', 'Frame_Shift_Del', 'Frame_Shift_Ins',
                'In_Frame_Del', 'In_Frame_Ins', 'Nonstop_Mutation'
            ]:
                chrom = row['Chromosome']
                pos = row['Start_Position']
                ref = row['Reference_Allele']
                alt = row['Tumor_Seq_Allele2']
                
                # Get depth info if available
                t_depth = row.get('t_depth', '100')
                t_ref = row.get('t_ref_count', '50')
                t_alt = row.get('t_alt_count', '50')
                n_depth = row.get('n_depth', '100')
                n_ref = row.get('n_ref_count', '100')
                
                if not t_depth or t_depth == '.': t_depth = '100'
                if not t_ref or t_ref == '.': t_ref = '50'
                if not t_alt or t_alt == '.': t_alt = '50'
                if not n_depth or n_depth == '.': n_depth = '100'
                if not n_ref or n_ref == '.': n_ref = '100'
                
                vcf_fh.write('{}\t{}\t.\t{}\t{}\tPASS\tPASS\t.\tGT:AD:DP\t0/0:{},{}:{}\t0/1:{},{}:{}\n'.format(
                    chrom, pos, ref, alt,
                    n_ref, '0', n_depth,
                    t_ref, t_alt, t_depth
                ))

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: maf_to_vcf.py input.maf output.vcf sample_name")
        sys.exit(1)
    maf_to_vcf(sys.argv[1], sys.argv[2], sys.argv[3])
    print("Done. Written to {}".format(sys.argv[2]))
