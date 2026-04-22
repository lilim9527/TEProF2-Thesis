"""Tests for GeneAnnotator using real dict structure from force_build_dic.py.

Dict structure: {chrom: {(tx_start, tx_end, attr): {(ex_start, ex_end): "10-col-string"}}}
10-col string: chrom,source,feature,start,end,score,strand,frame,attributes,exon_num
"""
from teprof3.annotation.gene_annotator import GeneAnnotator

ATTR = 'gene_id "GENE1"; transcript_id "ENST001"; gene_type "protein_coding"; gene_name "GENE1";'

PLUS_DIC = {
    "chr1": {
        (1000, 2000, ATTR): {
            (1000, 1100): f"chr1,.,exon,1000,1100,.,+,.,{ATTR} exon_number 1;,1",
            (1100, 1200): f"chr1,.,intron,1100,1200,.,+,.,{ATTR} exon_number 1;,1",
            (1200, 2000): f"chr1,.,exon,1200,2000,.,+,.,{ATTR} exon_number 2;,2",
        }
    }
}

MINUS_DIC = {
    "chr2": {
        (3000, 5000, ATTR): {
            (3000, 4000): f"chr2,.,exon,3000,4000,.,-,.,{ATTR} exon_number 2;,2",
            (4000, 5000): f"chr2,.,exon,4000,5000,.,-,.,{ATTR} exon_number 1;,1",
        }
    }
}


def test_annotate_tss_exon():
    ann = GeneAnnotator(PLUS_DIC, MINUS_DIC)
    result = ann.annotate_tss("chr1", 1050, "+")
    assert result["tss_gene_name"] == "GENE1"
    assert result["tss_element"] == "exon"
    assert result["tss_element_num"] == "1"


def test_annotate_tss_intergenic():
    ann = GeneAnnotator(PLUS_DIC, MINUS_DIC)
    result = ann.annotate_tss("chr1", 5000, "+")
    assert result["tss_gene_name"] == "None"
    assert result["tss_element"] == "intergenic"


def test_annotate_tss_wrong_chrom_is_intergenic():
    """Chromosome mismatch must not produce false annotations."""
    ann = GeneAnnotator(PLUS_DIC, MINUS_DIC)
    result = ann.annotate_tss("chr99", 1050, "+")
    assert result["tss_element"] == "intergenic"


def test_annotate_tend_intron():
    ann = GeneAnnotator(PLUS_DIC, MINUS_DIC)
    result = ann.annotate_tend("chr1", 1150, "+")
    assert result["tend_element"] == "intron"


def test_annotate_tss_prefers_exon_over_intron():
    """Position 1100 overlaps both exon [1000,1100] and intron [1100,1200]; exon wins."""
    ann = GeneAnnotator(PLUS_DIC, MINUS_DIC)
    result = ann.annotate_tss("chr1", 1100, "+")
    assert result["tss_element"] == "exon"


def test_annotate_tss_minus_strand():
    """Minus strand: prefer element with max e_end (most 3' genomic = most 5' transcript)."""
    ann = GeneAnnotator(PLUS_DIC, MINUS_DIC)
    result = ann.annotate_tss("chr2", 4500, "-")
    assert result["tss_element"] == "exon"
    assert result["tss_element_num"] == "1"  # exon 1 is [4000,5000] on minus strand


def test_all_tss_keys_present():
    ann = GeneAnnotator(PLUS_DIC, MINUS_DIC)
    result = ann.annotate_tss("chr1", 1050, "+")
    expected_keys = [
        "tss_gene_type", "tss_gene_name", "tss_chr", "tss_start", "tss_end",
        "tss_element", "tss_element_num", "tss_start_codon",
        "tss_transcript_start", "tss_transcript_end",
        "tss_element_list", "tss_ref_transcript_id",
    ]
    for k in expected_keys:
        assert k in result, f"Missing key: {k}"


def test_all_tend_keys_present():
    ann = GeneAnnotator(PLUS_DIC, MINUS_DIC)
    result = ann.annotate_tend("chr1", 5000, "+")
    expected_keys = [
        "tend_gene_type", "tend_gene_name", "tend_chr", "tend_start", "tend_end",
        "tend_element", "tend_element_num", "tend_start_codon",
        "tend_transcript_start", "tend_transcript_end",
        "tend_element_list", "tend_ref_transcript_id",
    ]
    for k in expected_keys:
        assert k in result, f"Missing key: {k}"
