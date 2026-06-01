import json
import os
import textwrap

import microbiorust
import pytest


# FIXTURES: Generating Mock Data ---
@pytest.fixture
def mock_gbk(tmp_path):
    """Creates a minimal valid GenBank file."""
    path = tmp_path / "test.gbk"
    content = textwrap.dedent("""
        LOCUS       source_1                 910 bp    DNA     linear   CON 01-NOV-2024
        DEFINITION  Escherichia coli K-12 substr. MG1655.
        ACCESSION   source_1
        VERSION     source_1
        KEYWORDS    .
        SOURCE      Escherichia coli K-12 substr. MG1655
          ORGANISM  Escherichia coli K-12 substr. MG1655
        FEATURES             Location/Qualifiers
             source          1..910
                             /organism="K-12 substr. MG1655"
                             /mol_type="DNA"
             gene            complement(1..354)
                             /locus_tag="b3304"
             CDS             complement(1..354)
                             /locus_tag="b3304"
                             /codon_start=1
                             /gene="rplR"
                             /translation="MDKKSARIRRATRARRKLQELGATRLVVHRTPRHIYAQVIAPNGSLVAASTVEKAIAEQLKYTGNKDAAAAVGKAVAERALEKGIKDVSFDRSGFQYHGRVQALDAAREAGLQ"
                             /product="50S ribosomal subunit protein L18"
             gene            complement(364..897)
                             /locus_tag="b3305"
             CDS             complement(364..897)
                             /locus_tag="b3305"
                             /codon_start=1
                             /gene="rplF"
                             /translation="MSRVAKAPVVVPAGVDVKINGQVITIKGKNGELTRTLNDAVEVKHNTLTFGPRDGYADGWAQAGTARALLNSMVIGVTEGFTKKLQLVGVGYRAAVKGNVINLSGFSHPVDHQLPAGITAECPTQTEIVLKGADKQVIGQVAADLRAYRRPEPYKGKGVRYADVVRTKEAKK"
                             /product="50S ribosomal subunit protein L6"
        ORIGIN
                1 TTAGAACTGA AGGCCAGCTT CACGGGCAGC ATCTGCCAGT GCCTGGACAC GACCATGATA
               61 TTGGAACCCG GAACGGTCAA AGGATACATC TTTGATGCCT TTTTCCAGAG CGCGTTCAGC
              121 GACAGCTTTA CCCACAGCTG CAGCCGCGTC TTTGTTACCG GTGTACTTCA GTTGTTCAGC
              181 GATAGCTTTT TCTACAGTAG AAGCAGCTAC CAGAACTTCA GAACCGTTCG GTGCAATTAC
              241 CTGTGCGTAA ATGTGACGCG GGGTACGATG TACCACCAGG CGAGTTGCGC CCAGCTCCTG
              301 GAGCTTGCGG CGTGCGCGGG TCGCACGACG GATACGAGCA GATTTCTTAT CCATAGTGTT
              361 ACCTTACTTC TTCTTAGCCT CTTTGGTACG CACGACTTCG TCGGCGTAAC GAACACCCTT
              421 GCCTTTATAA GGCTCAGGAC GACGGTAGGC GCGCAGATCC GCTGCAACCT GGCCGATCAC
              481 CTGCTTATCA GCGCCTTTCA GCACGATTTC AGTCTGAGTC GGACATTCAG CAGTGATACC
              541 CGCAGGCAGC TGATGGTCAA CAGGATGAGA GAAACCCAGA GACAGGTTAA TCACATTGCC
              601 TTTAACCGCT GCACGGTAAC CTACACCAAC CAGCTGCAGC TTCTTAGTGA AGCCTTCGGT
              661 AACACCGATA ACCATTGAGT TCAGCAGGGC ACGCGCGGTA CCAGCCTGTG CCCAACCGTC
              721 TGCGTAACCA TCACGCGGAC CGAAGGTCAG GGTATTATCT GCATGTTTAA CTTCAACAGC
              781 ATCGTTGAGA GTACGAGTCA GCTCGCCGTT TTTACCTTTG ATCGTAATAA CCTGACCGTT
              841 GATTTTTACG TCAACGCCGG CAGGAACAAC GACCGGTGCT TTAGCAACAC GAGACA
        //
    """)
    path.write_text(content)
    return str(path)


@pytest.fixture
def mock_msa(tmp_path):
    """Creates a mock FASTA alignment file."""
    path = tmp_path / "align.fasta"
    content = ">Seq1\nATGC--AT\n>Seq2\nATGC--TT\n>Seq3\nATGCGGTT\n"
    path.write_text(content)
    return str(path)


@pytest.fixture
def mock_blast_tab(tmp_path):
    """Creates a mock BLAST tabular file."""
    path = tmp_path / "results.tab"
    # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    content = "seqA\tseqB\t99.0\t100\t1\t0\t1\t100\t1\t100\t1e-10\t200.0\n"
    path.write_text(content)
    return str(path)


def test_gbk_loader_aliases(mock_gbk):
    """Verify all specialized GBK loaders return a valid collection"""
    # Test all three aliases to ensure they point to the correct InternalRecord::Gbk logic
    loaders = [
        microbiorust.gbk_to_faa,
        microbiorust.gbk_to_fna,
        microbiorust.gbk_to_ffn,
    ]

    for load_func in loaders:
        collection = load_func(mock_gbk)
        keys = list(collection.keys())

        # This checks if "source_1" is either the key itself (FNA)
        # or the prefix of the key (FAA/FFN)
        assert any(k.startswith("source_1") for k in keys), (
            f"Failed on {load_func.__name__}"
        )
        assert len(collection.keys()) > 0


def test_gbk_count_function(mock_gbk):
    """Test the standalone count function (verifies re-parsing works if needed)."""
    # This calls the genbank! macro directly again
    count = microbiorust.gbk_to_faa_count(mock_gbk)
    assert isinstance(count, int)
    assert count >= 2  # Based on b3304 and b3305 in your mock


def test_record_access_and_metadata(mock_gbk):
    """Test retrieving a record and its lazy-loaded metadata (PyFeatureInfo)."""
    collection = microbiorust.parse_gbk(mock_gbk)
    first_key = list(collection.keys())[0]
    record = collection[first_key]

    # create the features on the first record
    feature_collection = record.features()
    feat = feature_collection["b3304"]

    # Test Metadata Proxy (PyFeatureInfo)
    assert feat.gene == "rplR"
    assert feat.product == "50S ribosomal subunit protein L18"

    # Verify numeric extractions from our InternalFeatureAttributes match expected logic
    assert feat.strand == -1
    assert feat.start == 1
    assert feat.stop == 354
    assert feat.codon_start == 1


def test_sequence_getitem(mock_gbk):
    collection = microbiorust.parse_gbk(mock_gbk)
    first_key = list(collection.keys())[0]
    rec = collection[first_key]
    assert "source_1" in rec.id()
    # this triggers rec.__getitem__ to return a PySequenceInfo
    sequence_collection = rec.sequences()

    # check __getitem__ works directly
    first_tag = list(sequence_collection.keys())[0]
    first_seq = sequence_collection[first_tag]
    assert first_seq is not None

    # use of values()
    all_proteins = [info.faa for info in sequence_collection.values() if info.faa]
    assert len(all_proteins) > 0
    print("len ", len(all_proteins))

    # check items() gives consistent key-value pairs
    for tag, info in sequence_collection.items():
        assert sequence_collection[tag].faa == info.faa


def test_embl_loader_integration():
    """Verify EMBL files flow through the same Record"""
    collection = microbiorust.parse_embl("example.embl")
    first_key = list(collection.keys())[0]
    record = collection[first_key]

    # access the id and features
    feature_collection = record.features()
    assert record.id() == "AM236082"
    feat = feature_collection["pRL80002"]
    assert feat is not None
    # Test Metadata Proxy (PyFeatureInfo)
    assert feat.gene == "repBp8"
    assert feat.product == "replication protein RepB"

    # Verify numeric extractions from our InternalFeatureAttributes match expected logic
    assert feat.strand == 1
    assert feat.start == 1321
    assert feat.stop == 2280
    assert feat.codon_start == 1


def test_missing_keys_raise_errors(mock_gbk):
    """Verify we get clean Python KeyErrors instead of Rust panics."""
    collection = microbiorust.parse_gbk(mock_gbk)

    first_key = list(collection.keys())[0]
    record = collection[first_key]
    with pytest.raises(KeyError):
        _ = record.features()["non_existent_gene"]


def test_gbk_to_gff(mock_gbk):
    # this function writes to {filename}.gff and reads again to check
    microbiorust.gbk_to_gff(mock_gbk, dna=True)
    gff_path = f"{mock_gbk}.gff"
    assert os.path.exists(gff_path)
    with open(gff_path, "r") as f:
        assert "source_source_1_1" in f.read()


def test_faa_no_ffn_attribute(mock_gbk):
    """PyFaaInfo should not have ffn, PyFfnInfo should not have faa - clean types."""
    faa_collection = microbiorust.gbk_to_faa(mock_gbk)
    ffn_collection = microbiorust.gbk_to_ffn(mock_gbk)

    faa_info = list(faa_collection.values())[0]
    ffn_info = list(ffn_collection.values())[0]

    assert hasattr(faa_info, "faa")
    assert not hasattr(faa_info, "ffn")

    assert hasattr(ffn_info, "ffn")
    assert not hasattr(ffn_info, "faa")


def test_collection_len_and_contains(mock_gbk):
    """Test __len__ and __contains__ on all collection types."""
    collection = microbiorust.parse_gbk(mock_gbk)
    assert len(collection) == 1  # one record in mock

    first_key = list(collection.keys())[0]
    record = collection[first_key]

    feat_col = record.features()
    assert len(feat_col) == 2  # b3304 and b3305
    assert "b3304" in feat_col
    assert "b3305" in feat_col
    assert "fake" not in feat_col

    faa_col = record.faa()
    assert len(faa_col) == 2
    assert any("b3304" in k for k in faa_col.keys())


def test_collection_to_json(mock_gbk):
    """Test that Collections successfully serialize to valid JSON strings."""
    collection = microbiorust.parse_gbk(mock_gbk)
    # test RecordCollection serialization
    record_json_str = collection.to_json()
    assert isinstance(record_json_str, str)
    # verify it parses into a native Python dictionary
    record_dict = json.loads(record_json_str)
    assert len(record_dict) > 0
    # test a sub-collection (FaaCollection) serialization
    first_key = list(collection.keys())[0]
    faa_collection = collection[first_key].faa()
    faa_json_str = faa_collection.to_json()
    assert isinstance(faa_json_str, str)
    faa_dict = json.loads(faa_json_str)
    assert len(faa_dict) > 0


def test_collection_keyerror_message(mock_gbk):
    """Verify KeyError message is helpful and contains the requested tag."""
    collection = microbiorust.parse_gbk(mock_gbk)
    first_key = list(collection.keys())[0]
    record = collection[first_key]
    faa_collection = record.faa()

    bad_tag = "totally_fake_tag"
    with pytest.raises(KeyError) as exc_info:
        _ = faa_collection[bad_tag]
    # error message should contain the tag name and a hint
    assert bad_tag in str(exc_info.value)
    assert "keys()" in str(exc_info.value)


def test_gbk_to_ffn_returns_ffn_collection(mock_gbk):
    """Test gbk_to_ffn returns FfnCollection not SequenceCollection."""
    ffn_collection = microbiorust.gbk_to_ffn(mock_gbk)

    assert len(ffn_collection) > 0
    values = list(ffn_collection.values())
    assert all(hasattr(v, "ffn") for v in values)
    assert all(v.ffn is not None for v in values)
    # ffn should not have faa attribute
    assert not hasattr(values[0], "faa")


def test_gbk_to_faa_returns_faa_collection(mock_gbk):
    """Test gbk_to_faa returns FaaCollection not SequenceCollection."""
    faa_collection = microbiorust.gbk_to_faa(mock_gbk)

    assert len(faa_collection) > 0

    # verify values() items() keys() all work
    values = list(faa_collection.values())
    assert len(values) > 0
    assert all(hasattr(v, "faa") for v in values)
    assert all(hasattr(v, "locus_tag") for v in values)
    # faa is non-Optional - should never be None
    assert all(v.faa is not None for v in values)

    # items() consistency
    for tag, info in faa_collection.items():
        assert faa_collection[tag].faa == info.faa

    # __contains__
    first_key = list(faa_collection.keys())[0]
    assert first_key in faa_collection
    assert "definitely_not_a_real_tag" not in faa_collection


def test_record_faa_method(mock_gbk):
    """Test record.faa() returns FaaCollection of PyFaaInfo with non-Optional faa."""
    collection = microbiorust.parse_gbk(mock_gbk)
    first_key = list(collection.keys())[0]
    record = collection[first_key]

    faa_collection = record.faa()
    assert len(faa_collection) > 0

    # values are PyFaaInfo - faa is String not Optional
    for info in faa_collection.values():
        assert info.locus_tag is not None
        assert info.faa is not None
        assert len(info.faa) > 0


# tests for multiple sequence alignment
def test_subset_msa(mock_msa):
    # subset mock alignment: Rows 0-2 (Seq1 & Seq2), Cols 0-4 (ATGC)
    subset = microbiorust.subset_msa_alignment(mock_msa, (0, 2), (0, 4))
    print("subset", subset)
    assert len(subset) == 2  # 2 headers + 2 sequences
    assert ">Seq1" in subset[0]
    assert "ATGC" in subset[1]


def test_purge_gaps(mock_msa, tmp_path):
    out_path = str(tmp_path / "purged.fasta")
    # threshold 0.5 should remove the '--' columns in Seq1 and Seq2 and write to file
    microbiorust.purge_gaps(mock_msa, out_path, 0.5)
    assert os.path.exists(out_path)


def test_get_consensus(mock_msa):
    # given 'ATGC' is constant in the mock, it should be in consensus
    consensus = microbiorust.get_consensus(mock_msa)
    assert consensus.startswith("ATGC")


# tests for Sequence Metrics
def test_hydrophobicity():
    seq = "MALWMRLLPLLALLALWGPDPAAAFVN"
    scores = microbiorust.hydrophobicity(seq, window_size=3)
    assert len(scores) > 0
    assert all(isinstance(s, float) for s in scores)


def test_amino_counts():
    seq = "MATAG"
    counts = microbiorust.amino_counts(seq)
    assert counts["M"] == 1
    assert counts["A"] == 2


# test for Async Tabular Parser


def test_parse_tabular(mock_blast_tab):
    results = microbiorust.parse_tabular(mock_blast_tab)
    assert len(results) == 1
    assert results[0]["qseqid"] == "seqA"
    assert results[0]["bitscore"] == 200.0
