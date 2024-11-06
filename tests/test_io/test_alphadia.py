"""Tests for psm_utils.io.alphadia."""

from psm_utils.io.alphadia import AlphaDIAReader
from psm_utils.psm import PSM

test_psm = PSM(
    peptidoform="LMIEM[Oxidation]DGTENK/2",
    spectrum_id="79426",
    run="LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_03",
    collection=None,
    spectrum="79426",
    is_decoy=False,
    score=170.287918,
    qvalue=0.000000,
    pep=0.000005,
    precursor_mz=648.794128,
    retention_time=3111.141602,
    ion_mobility=0.000001,
    protein_list=["P06733"],
    rank=1,
    source="AlphaDIA",
    metadata={},
    rescoring_features={
        "rt_observed": 3111.141602,
        "mobility_observed": 0.000001,
        "mz_observed": 648.794128,
        "charge": 2,
        "delta_rt": -39.528809,
    },
)


class TestAlphaDIAReader:
    def test_iter(self):
        with AlphaDIAReader("./tests/test_data/test_alphadia.tsv") as reader:
            for psm in reader:
                psm.provenance_data = {}
                assert psm == test_psm

    def test__parse_peptidoform(self):
        test_cases = [
            {
                "sequence": "DNTTSGCGSDLQSATGTAR",
                "mods": "Carbamidomethyl@C",
                "mod_sites": "7",
                "charge": 2,
                "expected": "DNTTSGC[Carbamidomethyl]GSDLQSATGTAR/2",
            },
            {
                "sequence": "STCTEGEIACSADGK",
                "mods": "Carbamidomethyl@C;Carbamidomethyl@C",
                "mod_sites": "3;10",
                "charge": 2,
                "expected": "STC[Carbamidomethyl]TEGEIAC[Carbamidomethyl]SADGK/2",
            },
            {
                "sequence": "MLGETCADCGTILLQDK",
                "mods": "Oxidation@M;Carbamidomethyl@C;Carbamidomethyl@C",
                "mod_sites": "1;6;9",
                "charge": 2,
                "expected": "M[Oxidation]LGETC[Carbamidomethyl]ADC[Carbamidomethyl]GTILLQDK/2",
            },
            {
                "sequence": "VGLIGSCTNSSYEDMSR",
                "mods": "Oxidation@M;Carbamidomethyl@C",
                "mod_sites": "15;7",
                "charge": 2,
                "expected": "VGLIGSC[Carbamidomethyl]TNSSYEDM[Oxidation]SR/2",
            },
            {
                "sequence": "STATTTVTTSDQASHPTK",
                "mods": "Acetyl@Protein_N-term",
                "mod_sites": "0",
                "charge": 2,
                "expected": "[Acetyl]-STATTTVTTSDQASHPTK/2",
            },
            {
                "sequence": "MEPGPDGPAASGPAAIR",
                "mods": "Acetyl@Protein_N-term;Oxidation@M",
                "mod_sites": "0;1",
                "charge": 2,
                "expected": "[Acetyl]-M[Oxidation]EPGPDGPAASGPAAIR/2",
            },
            {
                "sequence": "AEPQPPSGGLTDEAALSCCSDADPSTK",
                "mods": "Acetyl@Protein_N-term;Carbamidomethyl@C;Carbamidomethyl@C",
                "mod_sites": "0;18;19",
                "charge": 3,
                "expected": "[Acetyl]-AEPQPPSGGLTDEAALSC[Carbamidomethyl]C[Carbamidomethyl]SDADPSTK/3",
            },
            {
                "sequence": "EPLISAPYLTTTKMSAPATLDAACIFCK",
                "mods": "Acetyl@Protein_N-term;Oxidation@M;Carbamidomethyl@C;Carbamidomethyl@C",
                "mod_sites": "0;14;24;27",
                "charge": 4,
                "expected": "[Acetyl]-EPLISAPYLTTTKM[Oxidation]SAPATLDAAC[Carbamidomethyl]IFC[Carbamidomethyl]K/4",
            },
            {
                "sequence": "GDIDANAFQHK",
                "mods": "Amidated@Any_C-term",
                "mod_sites": "-1",
                "charge": 2,
                "expected": "GDIDANAFQHK-[Amidated]/2",
            },
            {
                "sequence": "MNNPAMTIKGEQAK",
                "mods": "Acetyl@Protein_N-term;Oxidation@M;Oxidation@M;Amidated@Any_C-term",
                "mod_sites": "0;1;6;-1",
                "charge": 4,
                "expected": "[Acetyl]-M[Oxidation]NNPAM[Oxidation]TIKGEQAK-[Amidated]/4",
            },
        ]

        for test_case in test_cases:
            assert AlphaDIAReader._parse_peptidoform(
                test_case["sequence"], test_case["mods"], test_case["mod_sites"], test_case["charge"]
            ) == test_case["expected"]
