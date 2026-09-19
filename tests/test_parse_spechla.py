"""SpecHLA's wide table must reach the consensus.

SpecHLA is the only caller with no parse step: modules/spechla.nf copied its native
layout straight through.

    # version: IPD-IMGT/HLA 3.38.0
    Sample   HLA_A_1         HLA_A_2         HLA_B_1   ...
    HG00096  A*01:01:01:01   A*29:02:01:01   B*08:01:01:01 ...

modules/aggregation.nf reads only long format -- field 1 the gene, fields 2 and 3
the alleles -- so it saw the literal gene "Sample" on one line and a sample id on
the next, matched neither against mv_genes, and dropped both. SpecHLA ran, called
correctly, and contributed nothing to the vote on every run.

The fixtures here are the real outputs from the Roihu WGS and WES runs.
"""
from __future__ import annotations

import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "bin" / "parse_spechla_results.py"

HEADER = ("Sample\tHLA_A_1\tHLA_A_2\tHLA_B_1\tHLA_B_2\tHLA_C_1\tHLA_C_2\t"
          "HLA_DPA1_1\tHLA_DPA1_2\tHLA_DPB1_1\tHLA_DPB1_2\tHLA_DQA1_1\tHLA_DQA1_2\t"
          "HLA_DQB1_1\tHLA_DQB1_2\tHLA_DRB1_1\tHLA_DRB1_2")

# From the WGS run.
WGS = ("# version: IPD-IMGT/HLA 3.38.0\n" + HEADER + "\n"
       "HG00096\tA*01:01:01:01\tA*29:02:01:01\tB*08:01:01:01\tB*44:03:01:13\t"
       "C*16:01:01:01\tC*07:01:01:01\tDPA1*02:01:02:02\tDPA1*01:03:01:16\t"
       "DPB1*01:01:01:01\tDPB1*02:01:02:29\tDQA1*05:03:02\tDQA1*02:01:01:01\t"
       "DQB1*02:01:01\tDQB1*02:01:01\tDRB1*03:01:01:01\tDRB1*07:01:01:01\n")

# From the WES run, where SpecHLA called B homozygous.
WES = ("# version: IPD-IMGT/HLA 3.38.0\n" + HEADER + "\n"
       "HG00096\tA*29:02:04\tA*01:01:46\tB*08:01:01:01\tB*08:01:01:01\t"
       "C*07:01:49\tC*07:01:01:01\tDPA1*01:03:01:01\tDPA1*02:01:02:01\t"
       "DPB1*01:01:01:01\tDPB1*02:01:02:01\tDQA1*05:01:01:01\tDQA1*02:01:01:01\t"
       "DQB1*02:01:01\tDQB1*02:02:01\tDRB1*03:01:01:01\tDRB1*07:01:01:01\n")


def parse(text, sample="HG00096"):
    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp)
        src = d / "hla.result.txt"
        src.write_text(text, encoding="utf-8")
        out = d / "out.txt"
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--input", str(src),
             "--sample", sample, "--output", str(out)],
            capture_output=True, text=True)
        return result, (out.read_text(encoding="utf-8") if out.exists() else "")


def rows(text):
    parsed = {}
    for line in text.splitlines():
        if line.startswith("#") or line.startswith("Gene"):
            continue
        fields = line.split("\t")
        parsed[fields[0]] = (fields[1], fields[2])
    return parsed


class WideToLongTests(unittest.TestCase):
    def test_the_wgs_calls_come_through(self):
        result, text = parse(WGS)
        self.assertEqual(result.returncode, 0, result.stderr)
        got = rows(text)
        self.assertEqual(got["A"], ("A*01:01:01:01", "A*29:02:01:01"))
        self.assertEqual(got["B"], ("B*08:01:01:01", "B*44:03:01:13"))
        self.assertEqual(got["C"], ("C*16:01:01:01", "C*07:01:01:01"))

    def test_a_homozygous_call_keeps_both_copies(self):
        """SpecHLA reported B twice on WES; both fields must survive."""
        _, text = parse(WES)
        self.assertEqual(rows(text)["B"], ("B*08:01:01:01", "B*08:01:01:01"))

    def test_the_version_comment_is_not_taken_as_a_header(self):
        _, text = parse(WGS)
        self.assertNotIn("version", rows(text))

    def test_neither_sample_nor_the_sample_id_becomes_a_gene(self):
        """The exact failure: 'Sample' and 'HG00096' were being read as genes."""
        _, text = parse(WGS)
        got = rows(text)
        self.assertNotIn("Sample", got)
        self.assertNotIn("SAMPLE", got)
        self.assertNotIn("HG00096", got)

    def test_every_vote_gene_is_present(self):
        """params.mv_genes defaults to A,B,C,DQA1,DQB1,DRB1."""
        _, text = parse(WGS)
        got = rows(text)
        for gene in ("A", "B", "C", "DQA1", "DQB1", "DRB1"):
            self.assertIn(gene, got, gene)

    def test_a_header_only_file_fails_rather_than_writing_nothing(self):
        result, _ = parse("# version: IPD-IMGT/HLA 3.38.0\n" + HEADER + "\n")
        self.assertNotEqual(result.returncode, 0)

    def test_an_unrecognisable_file_fails_loudly(self):
        result, _ = parse("Gene\tAllele1\tAllele2\nA\tA*01:01\tA*29:02\n")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("column pairs", result.stderr)


class SurvivesAggregationTests(unittest.TestCase):
    """The parsed rows must pass modules/aggregation.nf's own filter."""

    VOTE_GENES = {"A", "B", "C", "DQA1", "DQB1", "DRB1"}
    MISSING = {"", "-", "NA", "None", "none", "."}

    def test_aggregation_would_keep_these_rows(self):
        _, text = parse(WGS)
        kept = 0
        for line in text.splitlines():
            # This mirrors the loop in modules/aggregation.nf.
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            gene = parts[0].strip().replace("HLA-", "").upper()
            if gene.lower() == "gene" or gene not in self.VOTE_GENES:
                continue
            a1, a2 = parts[1].strip(), parts[2].strip()
            self.assertNotIn(a1, self.MISSING)
            self.assertNotIn(a2, self.MISSING)
            kept += 1
        self.assertEqual(kept, 6, "aggregation would keep %d of 6 vote genes" % kept)


class AgreesWithTheAnalysisPathTests(unittest.TestCase):
    """bin/hla_benchmark.py reads SpecHLA too. The two must not drift apart."""

    def test_same_gene_column_pairing(self):
        spec = importlib.util.spec_from_file_location("hb", REPO / "bin" / "hla_benchmark.py")
        hb = importlib.util.module_from_spec(spec)
        try:
            spec.loader.exec_module(hb)
        except Exception as exc:                      # pragma: no cover
            self.skipTest("hla_benchmark.py not importable here: %s" % exc)

        sys.path.insert(0, str(REPO / "bin"))
        try:
            import parse_spechla_results as ps
        finally:
            sys.path.pop(0)

        header = HEADER.split("\t")
        theirs = set(hb.extract_wide_gene_columns(header))
        ours = set(ps.gene_columns(header))
        self.assertEqual(ours, theirs)


if __name__ == "__main__":
    unittest.main()
