"""
Gene Annotator - TSS and tend annotation against gencode dictionaries.

Dict structure (from force_build_dic.py):
  {chrom: {(tx_start, tx_end, attr): {(ex_start, ex_end): "10-col-string"}}}

10-col string format (tab-separated GTF cols 1-8 + attr + exon_num):
  col0=chrom, col1=source, col2=feature, col3=start, col4=end,
  col5=score, col6=strand, col7=frame, col8=attributes, col9=exon_num
"""
from __future__ import annotations

_NONE_ANNOTATION = {
    "gene_type": "None", "gene_name": "None", "chr": "None",
    "start": "None", "end": "None", "element": "intergenic",
    "element_num": "None", "start_codon": "None",
    "transcript_start": "None", "transcript_end": "None",
    "element_list": "None", "ref_transcript_id": "None",
}


def _parse_attr(attr_str: str) -> dict:
    """Extract gene_type, gene_name, transcript_id from GTF attribute string."""
    result = {"gene_type": "None", "gene_name": "None", "ref_transcript_id": "None"}
    for field in attr_str.split(";"):
        field = field.strip()
        if not field:
            continue
        if "gene_type" in field:
            val = field.split('"')[1] if '"' in field else field.split()[-1]
            result["gene_type"] = val
        elif "gene_name" in field:
            val = field.split('"')[1] if '"' in field else field.split()[-1]
            result["gene_name"] = val
        elif "transcript_id" in field:
            val = field.split('"')[1] if '"' in field else field.split()[-1]
            result["ref_transcript_id"] = val
    return result


def _parse_record(record: str, tx_start: int, tx_end: int) -> dict | None:
    """Parse 10-col comma-joined record string into annotation dict."""
    parts = record.split(",", 9)
    if len(parts) < 10:
        return None
    # parts: chrom,source,feature,start,end,score,strand,frame,attributes,exon_num
    attr_info = _parse_attr(parts[8])
    return {
        "chr": parts[0],
        "start": parts[3],
        "end": parts[4],
        "element": parts[2],          # "exon" or "intron"
        "element_num": parts[9].strip(),
        "start_codon": "None",
        "transcript_start": str(tx_start),
        "transcript_end": str(tx_end),
        "element_list": "None",
        **attr_info,
    }


class GeneAnnotator:
    """Annotate transcript TSS and tend positions against gencode dictionaries."""

    def __init__(self, plus_dic: dict, minus_dic: dict) -> None:
        self._plus = plus_dic
        self._minus = minus_dic

    def annotate_tss(self, chrom: str, tss: int, strand: str) -> dict:
        dic = self._plus if strand == "+" else self._minus
        candidates = self._find_candidates(chrom, tss, dic)
        best = self._select_best(candidates, strand)
        return self._prefix(best, "tss_")

    def annotate_tend(self, chrom: str, tend: int, strand: str) -> dict:
        dic = self._plus if strand == "+" else self._minus
        candidates = self._find_candidates(chrom, tend, dic)
        best = self._select_best(candidates, strand)
        return self._prefix(best, "tend_")

    def _find_candidates(self, chrom: str, pos: int, dic: dict) -> list[dict]:
        chrom_dic = dic.get(chrom, {})
        results = []
        for (t_start, t_end, *_), elements in chrom_dic.items():
            if not (t_start <= pos <= t_end):
                continue
            for (e_start, e_end), val in elements.items():
                if isinstance(val, str) and e_start <= pos <= e_end:
                    parsed = _parse_record(val, t_start, t_end)
                    if parsed:
                        parsed["_e_start"] = e_start
                        parsed["_e_end"] = e_end
                        results.append(parsed)
        return results

    def _select_best(self, candidates: list[dict], strand: str) -> dict:
        if not candidates:
            return dict(_NONE_ANNOTATION)
        exons = [c for c in candidates if c["element"] == "exon"]
        pool = exons if exons else candidates
        if strand == "+":
            return min(pool, key=lambda c: c["_e_start"])
        else:
            return max(pool, key=lambda c: c["_e_end"])

    def _prefix(self, ann: dict, prefix: str) -> dict:
        return {
            f"{prefix}gene_type": ann.get("gene_type", "None"),
            f"{prefix}gene_name": ann.get("gene_name", "None"),
            f"{prefix}chr": ann.get("chr", "None"),
            f"{prefix}start": ann.get("start", "None"),
            f"{prefix}end": ann.get("end", "None"),
            f"{prefix}element": ann.get("element", "intergenic"),
            f"{prefix}element_num": ann.get("element_num", "None"),
            f"{prefix}start_codon": ann.get("start_codon", "None"),
            f"{prefix}transcript_start": ann.get("transcript_start", "None"),
            f"{prefix}transcript_end": ann.get("transcript_end", "None"),
            f"{prefix}element_list": ann.get("element_list", "None"),
            f"{prefix}ref_transcript_id": ann.get("ref_transcript_id", "None"),
        }
