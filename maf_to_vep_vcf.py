#!/usr/bin/env python3
"""Convert GDC MAF to VEP-annotated VCF for pVACseq."""
import sys
import csv

def maf_to_vep_vcf(maf_file, vcf_file, tumor_sample, normal_sample=None):
    with open(maf_file, 'r') as fh:
        lines = [l for l in fh if not l.startswith('#')]
    
    reader = list(csv.DictReader(lines, delimiter='\t'))
    
    # Get normal sample name from first row if not provided
    if not normal_sample and reader:
        normal_sample = reader[0].get('Matched_Norm_Sample_Barcode', 'NORMAL')

    with open(vcf_file, 'w') as out:
        # VCF headers
        out.write('##fileformat=VCFv4.2\n')
        out.write('##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|PICK|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|UNIPROT_ISOFORM|SOURCE|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|FrameshiftSequence|WildtypeProtein">\n')
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n')
        out.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n')
        out.write(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{normal_sample}\t{tumor_sample}\n')

        kept = 0
        skipped = 0
        for row in reader:
            vc = row.get('Variant_Classification', '')
            vt = row.get('Variant_Type', '')
            
            if vc not in ['Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation',
                          'Splice_Site', 'Frame_Shift_Del', 'Frame_Shift_Ins',
                          'In_Frame_Del', 'In_Frame_Ins']:
                skipped += 1
                continue

            chrom = row['Chromosome']
            pos = row['Start_Position']
            ref = row['Reference_Allele']
            alt = row['Tumor_Seq_Allele2']
            
            if ref == '-':  # insertion
                ref = 'N'
                alt = 'N' + alt
            elif alt == '-':  # deletion
                alt = 'N'
                ref = 'N' + ref

            # Build CSQ field from MAF columns
            allele = alt
            consequence = row.get('One_Consequence', row.get('Consequence', ''))
            impact = row.get('IMPACT', '')
            symbol = row.get('SYMBOL', row.get('Hugo_Symbol', ''))
            gene = row.get('Gene', row.get('Entrez_Gene_Id', ''))
            feature_type = 'Transcript'
            feature = row.get('Feature', row.get('Transcript_ID', ''))
            biotype = row.get('BIOTYPE', 'protein_coding')
            exon = row.get('EXON', '')
            intron = row.get('INTRON', '')
            hgvsc = row.get('HGVSc', '')
            hgvsp = row.get('HGVSp', '')
            cdna_pos = row.get('cDNA_position', '')
            cds_pos = row.get('CDS_position', '')
            protein_pos = row.get('Protein_position', '')
            amino_acids = row.get('Amino_acids', '')
            codons = row.get('Codons', '')
            existing = row.get('Existing_variation', '')
            distance = row.get('DISTANCE', '')
            strand = row.get('TRANSCRIPT_STRAND', row.get('STRAND', ''))
            canonical = row.get('CANONICAL', '')
            ensp = row.get('ENSP', '')
            swissprot = row.get('SWISSPROT', '')
            wildtype_protein = row.get('WildtypeProtein', '')
            frameshift_seq = row.get('FrameshiftSequence', '')

            csq = '|'.join([
                allele, consequence, impact, symbol, gene,
                feature_type, feature, biotype, exon, intron,
                hgvsc, hgvsp, cdna_pos, cds_pos, protein_pos,
                amino_acids, codons, existing, distance, strand,
                '', '1', '', '', '', canonical, '', '', '', '', '', ensp,
                swissprot, '', '', '', '', '', '', '', '', '', '',
                '', '', '', '', '', '', '', '', '', '', '', '', '',
                '', '', '', '', '', '', '', '', '', '', '', '',
                frameshift_seq, wildtype_protein
            ])

            info = f'CSQ={csq}'

            # Depth info
            t_depth = row.get('t_depth', '100') or '100'
            t_ref = row.get('t_ref_count', '50') or '50'
            t_alt_count = row.get('t_alt_count', '50') or '50'
            n_depth = row.get('n_depth', '100') or '100'
            n_ref = row.get('n_ref_count', '100') or '100'

            normal_gt = f'0/0:{n_ref},0:{n_depth}'
            tumor_gt = f'0/1:{t_ref},{t_alt_count}:{t_depth}'

            out.write(f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info}\tGT:AD:DP\t{normal_gt}\t{tumor_gt}\n')
            kept += 1

        print(f"Done. Kept {kept} variants, skipped {skipped}. Written to {vcf_file}")

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print("Usage: maf_to_vep_vcf.py input.maf output.vcf tumor_sample_name [normal_sample_name]")
        sys.exit(1)
    normal = sys.argv[4] if len(sys.argv) > 4 else None
    maf_to_vep_vcf(sys.argv[1], sys.argv[2], sys.argv[3], normal)
