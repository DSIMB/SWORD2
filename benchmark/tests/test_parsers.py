import json
from pathlib import Path

import pytest

from benchmark.datasets import CathEntry, parse_cath_domain_string, read_merizo_csv
from benchmark.numbering import (
    ResidueKey,
    StructureNumbering,
    map_author_chopping,
    map_one_based_chopping,
)
from benchmark.runners.chainsaw import chainsaw_to_common_chopping, parse_chainsaw_tsv
from benchmark.runners.merizo import parse_merizo_tsv
from benchmark.runners.sword2_rust import parse_summary_partitions


def test_cath_domain_parser_removes_superfamily_suffix_and_masks():
    domains = parse_cath_domain_string("100-129_168-292:40_1078|1-99_130-167:999_999")

    assert domains == [["100-129", "168-292"], ["1-99", "130-167"]]


def test_read_merizo_cath663_csv_row_shape():
    entries = read_merizo_csv(
        Path("/home/chili/cretin/PROJECTS/Merizo/datasets/merizo_domains/CATH-663.csv"),
        dataset="cath663",
        limit=1,
    )

    assert entries == [
        CathEntry(
            pdb_id="19hc",
            chain_id="A",
            entry_id="19hcA",
            n_domains=2,
            n_residues=292,
            chopping="100-129_168-292:40_1078|1-99_130-167:40_1078",
            dataset="cath663",
        )
    ]


def test_numbering_maps_author_and_sequential_choppings_to_zero_based_space():
    numbering = StructureNumbering(
        residues=[
            ResidueKey("A", "10", ""),
            ResidueKey("A", "11", ""),
            ResidueKey("A", "12", ""),
            ResidueKey("A", "20", ""),
        ]
    )

    assert map_author_chopping("10-11|20", numbering) == "0-1,3"
    assert map_one_based_chopping("1-2|4") == "0-1,3"


def test_numbering_maps_negative_author_ranges():
    numbering = StructureNumbering(
        residues=[
            ResidueKey("A", "-2", ""),
            ResidueKey("A", "-1", ""),
            ResidueKey("A", "0", ""),
            ResidueKey("A", "1", ""),
            ResidueKey("A", "2", ""),
        ]
    )

    assert map_author_chopping("-2-1|2", numbering, chain_id="A") == "0-3,4"


def test_author_chopping_trims_missing_terminal_range_endpoint():
    numbering = StructureNumbering(
        residues=[
            ResidueKey("A", "1", ""),
            ResidueKey("A", "2", ""),
            ResidueKey("A", "3", ""),
            ResidueKey("A", "4", ""),
        ]
    )

    assert map_author_chopping("1-5", numbering, chain_id="A") == "0-3"


def test_chainsaw_chopping_falls_back_to_author_numbering_when_sequential_is_impossible():
    numbering = StructureNumbering(
        residues=[
            ResidueKey("B", "700", ""),
            ResidueKey("B", "1", ""),
            ResidueKey("B", "2", ""),
            ResidueKey("B", "3", ""),
            ResidueKey("B", "4", ""),
            ResidueKey("B", "5", ""),
        ]
    )

    assert chainsaw_to_common_chopping("700-2,3-5", numbering, chain_id="B") == "0-2,3-5"


def test_sword2_summary_parser_extracts_optimal_and_alternatives():
    summary = {
        "Optimal partition": {
            "Nb. domains": 2,
            "Domains": {
                "Domain 1": {"PUs": {"1-3": {}, "8-9": {}}},
                "Domain 2": {"PUs": {"4-7": {}}},
            },
        },
        "Alternative partition 1": {
            "Nb. domains": 1,
            "Domains": {"Domain 1": {"PUs": {"1-9": {}}}},
        },
    }

    partitions = parse_summary_partitions(summary, numbering=None)

    assert [(p.name, p.variant, p.chopping) for p in partitions] == [
        ("Optimal partition", "optimal", "0-2_7-8,3-6"),
        ("Alternative partition 1", "alternative", "0-8"),
    ]


def test_merizo_and_chainsaw_tsv_parsers():
    merizo_tsv = (
        "chain_id\tsequence_md5\tnres\tndom\tchopping\tconfidence\ttime_sec\n"
        "x\tmd5\t20\t2\t7-10_15-20,1-6\t0.9\t1.2\n"
    )
    chainsaw_tsv = "chain_id\tchopping\nA\t0-9|10-19\n"

    assert parse_merizo_tsv(merizo_tsv).chopping == "7-10_15-20,1-6"
    assert parse_chainsaw_tsv(chainsaw_tsv).chopping == "0-9|10-19"


def test_sword2_summary_parser_rejects_missing_optimal_partition():
    with pytest.raises(KeyError):
        parse_summary_partitions(json.loads("{}"), numbering=None)
