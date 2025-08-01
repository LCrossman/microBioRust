use microBioRust::embl::{gbk_write, gff_write, RangeValue, Record};
use std::collections::BTreeMap;


 /// Test to create a new record
     /// We require a source, features, sequence features and a sequence
     /// The source is top level, a single genbank file has one source, multi-genbank has one per contig
     /// The SourceAttributes construct has a name (counter), start, stop, organism, moltype, strain, type material and db_xref
     /// The FeatureAttributes construct has a locus tag (counter), gene, product, start, stop, codon start, strand
     /// SourceAttribute start and stop are the coordinates of the source feature or per contig, FeatureAttributes start and stop are per coding sequence (CDS)
     /// The SequenceAttributes construct has a locus tag (counter), start, stop, sequence_ffn, sequence_faa, codon start, and strand
     /// SequenceAttribute start and stop, codon start and strand are duplicates of those in the FeatureAttributes
     /// To add an entry requires using the set_ values such as set_start, set_stop, set_counter, set_strand
     /// To write in GFF format requires gff_write(seq_region, record_vec, filename and true/false
     /// The seq_region is the region of interest with name and DNA coordinates such as ``` "source_1".to_string(), (1,897) ```
     /// record_vec is a list of the records.  If there is only one record ``` vec![record] ``` will suffice
     /// filename is the required filename string, true/false is whether the DNA sequence should be included in the GFF3 file
     /// Some GFF3 files have the DNA sequence, whilst others do not.  Some tools require the DNA sequence included.

#[test]
fn create_new_record() -> Result<(), anyhow::Error> {
    //let filename = format!("new_record.gff");
    let mut record = Record::new();
    let mut seq_region: BTreeMap<String, (u32, u32)> = BTreeMap::new();
    seq_region.insert("source_1".to_string(), (1, 910));
    record
        .source_map
        .set_counter("source_1".to_string())
        .set_start(RangeValue::Exact(1))
        .set_stop(RangeValue::Exact(910))
        .set_organism("Escherichia coli".to_string())
        .set_mol_type("DNA".to_string())
        .set_strain("K-12 substr. MG1655".to_string())
        // culture_collection.clone()
        .set_type_material("type strain of Escherichia coli K12".to_string())
        .set_db_xref("PRJNA57779".to_string());
    record
        .cds
        .set_counter("b3304".to_string())
        .set_start(RangeValue::Exact(1))
        .set_stop(RangeValue::Exact(354))
        .set_gene("rplR".to_string())
        .set_product("50S ribosomal subunit protein L18".to_string())
        .set_codon_start(1)
        .set_strand(-1);
    record
        .cds
        .set_counter("b3305".to_string())
        .set_start(RangeValue::Exact(364))
        .set_stop(RangeValue::Exact(897))
        .set_gene("rplF".to_string())
        .set_product("50S ribosomal subunit protein L6".to_string())
        .set_codon_start(1)
        .set_strand(-1);
    record
        .seq_features
        .set_counter("b3304".to_string())
        .set_start(RangeValue::Exact(1))
        .set_stop(RangeValue::Exact(354))
        .set_sequence_ffn(
            "ATGGATAAGAAATCTGCTCGTATCCGTCGTGCGACCCGCGCACGCCGCAAGCTCCAGGAG
CTGGGCGCAACTCGCCTGGTGGTACATCGTACCCCGCGTCACATTTACGCACAGGTAATT
GCACCGAACGGTTCTGAAGTTCTGGTAGCTGCTTCTACTGTAGAAAAAGCTATCGCTGAA
CAACTGAAGTACACCGGTAACAAAGACGCGGCTGCAGCTGTGGGTAAAGCTGTCGCTGAA
CGCGCTCTGGAAAAAGGCATCAAAGATGTATCCTTTGACCGTTCCGGGTTCCAATATCAT
GGTCGTGTCCAGGCACTGGCAGATGCTGCCCGTGAAGCTGGCCTTCAGTTCTAA"
                .to_string(),
        )
        .set_sequence_faa(
            "MDKKSARIRRATRARRKLQELGATRLVVHRTPRHIYAQVIAPNGSEVLVAASTVEKAIAE
QLKYTGNKDAAAAVGKAVAERALEKGIKDVSFDRSGFQYHGRVQALADAAREAGLQF"
                .to_string(),
        )
        .set_codon_start(1)
        .set_strand(-1);
    record
        .seq_features
        .set_counter("b3305".to_string())
        .set_start(RangeValue::Exact(364))
        .set_stop(RangeValue::Exact(897))
        .set_sequence_ffn(
            "ATGTCTCGTGTTGCTAAAGCACCGGTCGTTGTTCCTGCCGGCGTTGACGTAAAAATCAAC
GGTCAGGTTATTACGATCAAAGGTAAAAACGGCGAGCTGACTCGTACTCTCAACGATGCT
GTTGAAGTTAAACATGCAGATAATACCCTGACCTTCGGTCCGCGTGATGGTTACGCAGAC
GGTTGGGCACAGGCTGGTACCGCGCGTGCCCTGCTGAACTCAATGGTTATCGGTGTTACC
GAAGGCTTCACTAAGAAGCTGCAGCTGGTTGGTGTAGGTTACCGTGCAGCGGTTAAAGGC
AATGTGATTAACCTGTCTCTGGGTTTCTCTCATCCTGTTGACCATCAGCTGCCTGCGGGT
ATCACTGCTGAATGTCCGACTCAGACTGAAATCGTGCTGAAAGGCGCTGATAAGCAGGTG
ATCGGCCAGGTTGCAGCGGATCTGCGCGCCTACCGTCGTCCTGAGCCTTATAAAGGCAAG
GGTGTTCGTTACGCCGACGAAGTCGTGCGTACCAAAGAGGCTAAGAAGAAGTAA"
                .to_string(),
        )
        .set_sequence_faa(
            "MSRVAKAPVVVPAGVDVKINGQVITIKGKNGELTRTLNDAVEVKHADNTLTFGPRDGYAD
GWAQAGTARALLNSMVIGVTEGFTKKLQLVGVGYRAAVKGNVINLSLGFSHPVDHQLPAG
ITAECPTQTEIVLKGADKQVIGQVAADLRAYRRPEPYKGKGVRYADEVVRTKEAKKK"
                .to_string(),
        )
        .set_codon_start(1)
        .set_strand(-1);
    record.sequence = "TTAGAACTGAAGGCCAGCTTCACGGGCAGCATCTGCCAGTGCCTGGACACGACCATGATA
TTGGAACCCGGAACGGTCAAAGGATACATCTTTGATGCCTTTTTCCAGAGCGCGTTCAGC
GACAGCTTTACCCACAGCTGCAGCCGCGTCTTTGTTACCGGTGTACTTCAGTTGTTCAGC
GATAGCTTTTTCTACAGTAGAAGCAGCTACCAGAACTTCAGAACCGTTCGGTGCAATTAC
CTGTGCGTAAATGTGACGCGGGGTACGATGTACCACCAGGCGAGTTGCGCCCAGCTCCTG
GAGCTTGCGGCGTGCGCGGGTCGCACGACGGATACGAGCAGATTTCTTATCCATAGTGTT
ACCTTACTTCTTCTTAGCCTCTTTGGTACGCACGACTTCGTCGGCGTAACGAACACCCTT
GCCTTTATAAGGCTCAGGACGACGGTAGGCGCGCAGATCCGCTGCAACCTGGCCGATCAC
CTGCTTATCAGCGCCTTTCAGCACGATTTCAGTCTGAGTCGGACATTCAGCAGTGATACC
CGCAGGCAGCTGATGGTCAACAGGATGAGAGAAACCCAGAGACAGGTTAATCACATTGCC
TTTAACCGCTGCACGGTAACCTACACCAACCAGCTGCAGCTTCTTAGTGAAGCCTTCGGT
AACACCGATAACCATTGAGTTCAGCAGGGCACGCGCGGTACCAGCCTGTGCCCAACCGTC
TGCGTAACCATCACGCGGACCGAAGGTCAGGGTATTATCTGCATGTTTAACTTCAACAGC
ATCGTTGAGAGTACGAGTCAGCTCGCCGTTTTTACCTTTGATCGTAATAACCTGACCGTT
GATTTTTACGTCAACGCCGGCAGGAACAACGACCGGTGCTTTAGCAACACGAGACA"
        .to_string();
    gff_write(
        seq_region.clone(),
        vec![record.clone()],
        "new_output_embl.gff",
        true,
    )?;
    gbk_write(seq_region, vec![record], "new_output_embl.gbk")?;
    return Ok(());
}
