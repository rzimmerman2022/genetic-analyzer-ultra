"""Common test fixtures for the genetic analyzer tests."""

import pytest

@pytest.fixture
def toy_vcf_content():
    # A minimal VCF content string including a header and one variant (rs429358)
    # This variant is used in test_direction_conflict
    return (
        "##fileformat=VCFv4.2\n"
        "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n"
        "1\t1000\trs429358\tC\tT\t.\tPASS\t.\tGT\t0/1\n" # Example genotype for rs429358
        "1\t2000\trs7412\tC\tT\t.\tPASS\t.\tGT\t0/0\n"  # Another APOE variant
        "1\t3000\trs123\tA\tG\t.\tPASS\t.\tGT\t1/1\n"   # A generic variant for PRS
        "1\t4000\trs456\tC\tT\t.\tPASS\t.\tGT\t0/1\n"   # Another generic variant for PRS
    )

@pytest.fixture
def toy_vcf(tmp_path, toy_vcf_content):
    vcf_file = tmp_path / "toy.vcf"
    vcf_file.write_text(toy_vcf_content)
    return str(vcf_file) # Return path as string, as expected by AdvancedGeneticAnalyzer

@pytest.fixture
def mini_prs_model():
    """
    Provides a small Polygenic Risk Score (PRS) model.
    Keys are RSIDs, values are dicts with 'weight' and optionally 'se_weight'.
    """
    return {
       "rs123": {"weight": 0.25, "se_weight": 0.05, "risk_allele": "G"},
       "rs456": {"weight": -0.30, "se_weight": 0.07, "risk_allele": "T"},
       # Add more variants as needed for comprehensive PRS testing
    }

@pytest.fixture
def mock_variant_db():
    """
    Provides a mock variant database similar to self.known_variants.
    This can be used to monkeypatch the analyzer's internal database for specific tests.
    """
    return {
        'rs429358': {
            'gene': 'APOE',
            'trait': 'Alzheimer\'s disease risk',
            'risk_allele': 'T', # Assuming T is risk for this test setup
            'effect_size': 2.0, 
            'ci_95': (1.5, 2.5),
            'pmid': '12345678',
            'mechanism': 'Test mechanism'
        },
        # Add other variants as needed for tests
    }
