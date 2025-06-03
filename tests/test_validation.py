import pytest
from genetic_analyzer_ultra import AdvancedGeneticAnalyzer
# Assuming your validation rules and logic are accessible for monkeypatching
# or that the test can modify the analyzer's state appropriately.

def test_direction_conflict(toy_vcf, monkeypatch, mock_variant_db):
    """
    Tests if the validation harness flags a "direction conflict" when a known
    variant's Odds Ratio (OR) is flipped from risk (>1) to protective (<1)
    or vice-versa, against its literature-consistent direction.
    """
    analyzer = AdvancedGeneticAnalyzer(toy_vcf)
    analyzer.load_data()

    # Original OR for rs429358 in mock_variant_db is 2.0 (risk)
    # We will use this mock_variant_db for the test.
    # The actual known_variants in the analyzer might be different,
    # so we patch it for this specific test.
    
    # Patch the analyzer's internal variant database for this test
    analyzer.known_variants = mock_variant_db
    
    # Simulate the scenario: Flip the OR of rs429358 to be protective (0.2)
    # This variant is expected to be a risk variant (OR > 1) based on mock_variant_db
    # The key is that the validation logic should compare the *calculated/observed* OR
    # (which would be 0.2 after this flip) against an *expected direction* (risk).
    
    # To make this test work, the validation.py's `validate` function
    # needs to be designed to detect such conflicts.
    # For this example, let's assume `analyze_disease_risk` populates results
    # and then `validation.validate` checks these results.
    
    # Modify the OR directly in the analyzer's representation of this variant
    # This simulates if the input data or an intermediate step led to an inverted OR
    # for a variant that the literature (mock_variant_db) says should be risk.
    
    # If the analyzer directly uses its self.known_variants for risk calculation AND validation expectation:
    # We need a way for the validation to know the "original" expected direction.
    # Let's assume the validation logic in `validation.py` has access to a
    # "ground truth" or can infer expected direction.
    
    # For the sake of this test, we'll directly modify the 'effect_size' (OR)
    # that would be used in `_calculate_variant_risk`.
    # The `validate` function in `validation.py` would then compare this against
    # an internal or external benchmark.
    
    # Simulate that the variant's data in the analyzer now shows a protective OR
    # This is what `_calculate_variant_risk` would use.
    if 'rs429358' in analyzer.known_variants:
        analyzer.known_variants['rs429358']['effect_size'] = 0.2 # Flipped OR
        # If CI is also used by validation for consistency, flip it too (optional for this test)
        # analyzer.known_variants['rs429358']['ci_95'] = (0.1, 0.3)


    # Run the disease risk analysis, which should use the flipped OR for rs429358
    analyzer.analyze_disease_risk() # This populates self.results['disease_risk'] and self.results['validation_summary_report']
    
    validation_summary = analyzer.results.get('validation_summary_report', [])
    
    found_conflict = False
    conflict_details = ""
    for report_item in validation_summary:
        if report_item.get('rule_name') == 'APOE_Alz_Direction' and report_item.get('status') == 'DIRECTION_CONFLICT': # Hypothetical rule name and status
            found_conflict = True
            conflict_details = report_item.get('details', '')
            break
        # A more generic check if the exact rule name/status isn't predefined for this test
        elif "direction conflict" in str(report_item.get('details', '')).lower() and report_item.get('status', '').upper() in ['FAIL', 'CONCERN', 'CONFLICT', 'DIRECTION_CONFLICT']:
             found_conflict = True
             conflict_details = report_item.get('details', '')
             break


    assert found_conflict, f"Validation harness did not flag 'direction conflict'. Details: {conflict_details}. Validation Report: {validation_summary}"
    print(f"Direction conflict correctly flagged for rs429358. Details: {conflict_details}")

# Note: This test assumes that `validation.py` contains logic to compare
# the calculated/observed effect direction of a variant against an expected
# direction (e.g., from literature or a benchmark database).
# The `VALIDATION_RULES` in `validation.py` would need a rule like 'APOE_Alz_Direction'
# that checks if rs429358 (APOE) is showing a protective effect when it's expected to be risk,
# or vice-versa, and sets status to 'DIRECTION_CONFLICT'.
# The current `validation.py` provided earlier focuses on value deviation, not direction conflict.
# It would need to be extended for this test to pass as written.
