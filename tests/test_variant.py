"""Regression tests for the variant classifier (scientific core).

Self-contained: builds a synthetic primary, derives candidates with controlled
differences, and checks each lands in the right category. Skipped if minimap2
is not on PATH (the classifier shells out to it).
"""

import random
import shutil

import pytest

from easy45 import variant

pytestmark = pytest.mark.skipif(shutil.which("minimap2") is None,
                                reason="minimap2 not on PATH")

_NEXT = {"A": "C", "C": "G", "G": "T", "T": "A"}


def _primary():
    random.seed(1)
    seq = [random.choice("ACGT") for _ in range(3000)]
    seq[100:108] = list("AAAAAAAA")          # a homopolymer run to slip in
    return "".join(seq)


def _mut_subs(seq, n):
    s = list(seq)
    step = len(s) // (n + 1)
    for i in range(1, n + 1):
        s[i * step] = _NEXT[s[i * step]]
    return "".join(s)


def _verdict(tmp_path, primary, candidate, min_sites=2):
    pf = tmp_path / "p.fa"; pf.write_text(">p\n" + primary + "\n")
    cf = tmp_path / "c.fa"; cf.write_text(">c\n" + candidate + "\n")
    return variant.confirm_variant(pf, cf, tmp_path / "cmp.paf", min_sites)


def test_real_variant(tmp_path):
    p = _primary()
    v = _verdict(tmp_path, p, _mut_subs(p, 15))
    assert v.category == "variant"
    assert v.substitutions >= 10 and v.identity >= 0.99


def test_single_sub_is_noise(tmp_path):
    p = _primary()
    v = _verdict(tmp_path, p, _mut_subs(p, 1))
    assert v.category == "noise"


def test_homopolymer_indel_is_noise(tmp_path):
    p = _primary()
    cand = p[:100] + p[100] * 2 + p[100:]     # expand the homopolymer run
    v = _verdict(tmp_path, p, cand)
    assert v.category == "noise"
    assert v.substitutions == 0


def test_low_identity_is_foreign(tmp_path):
    p = _primary()
    v = _verdict(tmp_path, p, _mut_subs(p, 600))   # ~20% divergence
    assert v.category == "foreign"
