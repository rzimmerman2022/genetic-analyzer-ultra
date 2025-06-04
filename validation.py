import effect_utils # Added import

VALIDATION_RULES = {
    'APOE_Alz_Value_Check': { # Renamed for clarity, this is the value check
        'rsids': ['rs429358', 'rs7412'],
        'expected_risk_allele_counts': { # Genotype to risk allele count mapping
            # Assuming rs429358 'C' is risk, rs7412 'T' is risk (for APOE e4/e2)
            # This needs to be specific to how APOE genotype is determined
            # For e3/e4 (one C at rs429358, no T at rs7412 for e2)
            # For e4/e4 (two C at rs429358, no T at rs7412 for e2)
            # This is a simplified example; actual APOE genotyping is more complex
        },
        'metric_path': ['disease_risk', 'neurological'], # Path to find the relevant APOE result
        'target_rsid_for_metric': 'rs429358', # Which rsid's result to check if multiple are involved
        'expected_values': { # Expected relative_risk based on a simplified APOE e4 count
            1: 3.0,  # Approx for one e4 allele (e.g., e3/e4)
            2: 12.0  # Approx for two e4 alleles (e.g., e4/e4)
        },
        'value_key': 'relative_risk', # Key in the finding dict that holds the value to check
        'check_type': 'value_deviation' # Specify check type
    },
    'APOE_Alz_Direction': { # New rule for direction conflict
        'rsids': ['rs429358'], # Focus on one variant for direction
        'metric_path': ['disease_risk', 'neurological'],
        'target_rsid_for_metric': 'rs429358',
        'value_key': 'relative_risk', # The OR value to check direction of
        'expected_direction': 'risk', # Literature expects this to be a risk variant (OR > 1)
        'check_type': 'direction_conflict' # Specify check type
    },
    'C282Y_HFE_Penetrance': {
        'rsids': ['rs28940279'],
        'metric_path': ['rare_variants'], # Path to find HFE C282Y result
        'target_rsid_for_metric': 'rs28940279',
        'genotype_to_check': 'AA', # Assuming 'A' is the risk allele for C282Y
        'expected_values': { # Expected penetrance %
            'AA': 0.28 # 28% penetrance for AA homozygotes
        },
        'value_key': 'penetrance_estimate' # Hypothetical key, would need to be added to HFE results
    }
}

def get_nested_value(data_dict, path_list, target_rsid, value_key_in_finding):
    """
    Retrieves a value from a nested dictionary structure,
    specifically looking for a finding related to target_rsid.
    """
    current_level = data_dict
    for key in path_list:
        if isinstance(current_level, dict) and key in current_level:
            current_level = current_level[key]
        elif isinstance(current_level, list): # Handle lists of findings
            found_item = None
            for item in current_level:
                if isinstance(item, dict) and item.get('rsid') == target_rsid:
                    found_item = item
                    break
            if found_item:
                current_level = found_item
            else:
                return None # Target RSID not found in list
        else:
            return None # Path broken

    # If current_level is now the specific finding dict
    if isinstance(current_level, dict) and current_level.get('rsid') == target_rsid:
        return current_level.get(value_key_in_finding)
    # If current_level is a list (e.g. disease_risk categories are lists of findings)
    elif isinstance(current_level, list):
         for item in current_level:
            if isinstance(item, dict) and item.get('rsid') == target_rsid:
                return item.get(value_key_in_finding)
    return None


def validate(results: dict) -> list:
    """
    Validates analysis results against predefined rules.
    Returns a list of validation outcomes.
    """
    report = []
    for rule_name, config in VALIDATION_RULES.items():
        rule_outcome = {
            'rule_name': rule_name,
            'status': 'NOT_APPLICABLE',
            'details': 'Rule conditions not met or data not found in results.'
        }

        check_type = config.get('check_type')

        if check_type == 'value_deviation':
            # Logic for APOE_Alz_Value_Check (formerly APOE_Alz)
            if rule_name == 'APOE_Alz_Value_Check':
                apoe_e4_alleles = 0 # Placeholder
                rs429358_genotype = None
                if results.get('disease_risk') and results['disease_risk'].get('neurological'):
                    for finding in results['disease_risk']['neurological']:
                        if finding.get('rsid') == 'rs429358': # This rule specifically checks rs429358 based on e4 count
                            rs429358_genotype = finding.get('genotype')
                            if rs429358_genotype:
                                # This is a simplified e4 allele count based on rs429358's 'C' allele.
                                # Real APOE genotyping is more complex (involves rs7412).
                                # For this rule, we assume 'C' at rs429358 contributes to e4.
                                apoe_e4_alleles = rs429358_genotype.count('C') 
                            break
                
                if apoe_e4_alleles in config['expected_values']:
                    expected_value = config['expected_values'][apoe_e4_alleles]
                    calculated_value = get_nested_value(results, config['metric_path'], config['target_rsid_for_metric'], config['value_key'])

                    if calculated_value is not None:
                        deviation = abs(calculated_value - expected_value) / expected_value if expected_value != 0 else float('inf')
                        if deviation <= 0.10: # 10% tolerance
                            rule_outcome['status'] = 'PASS'
                            rule_outcome['details'] = f"Calculated OR: {calculated_value:.2f}, Expected OR: {expected_value:.2f} (Deviation: {deviation:.2%})"
                        else:
                            rule_outcome['status'] = 'CONCERN'
                            rule_outcome['details'] = f"Calculated OR: {calculated_value:.2f}, Expected OR: {expected_value:.2f} (Deviation: {deviation:.2%})"
                    else:
                        rule_outcome['details'] = f"Calculated value for {config['target_rsid_for_metric']} not found."
                else:
                    rule_outcome['details'] = f"APOE e4 allele count ({apoe_e4_alleles}) derived from rs429358 genotype ({rs429358_genotype}) not in expected values for rule."

            # Logic for C282Y_HFE_Penetrance
            elif rule_name == 'C282Y_HFE_Penetrance':
                hfe_genotype = None
                if results.get('rare_variants'):
                    for finding in results['rare_variants']:
                        if finding.get('rsid') == config['target_rsid_for_metric']:
                            hfe_genotype = finding.get('genotype')
                            break
                
                if hfe_genotype == config['genotype_to_check']:
                    expected_value = config['expected_values'][hfe_genotype]
                    calculated_value = get_nested_value(results, config['metric_path'], config['target_rsid_for_metric'], config['value_key'])

                    if calculated_value is not None:
                        deviation = abs(calculated_value - expected_value) / expected_value if expected_value != 0 else float('inf')
                        if deviation <= 0.10:
                            rule_outcome['status'] = 'PASS'
                            rule_outcome['details'] = f"Calculated: {calculated_value:.2f}, Expected: {expected_value:.2f} (Deviation: {deviation:.2%})"
                        else:
                            rule_outcome['status'] = 'CONCERN'
                            rule_outcome['details'] = f"Calculated: {calculated_value:.2f}, Expected: {expected_value:.2f} (Deviation: {deviation:.2%})"
                    else:
                        rule_outcome['details'] = f"Calculated value for {config['target_rsid_for_metric']} not found."
                else:
                    rule_outcome['details'] = f"Genotype {hfe_genotype} for {config['target_rsid_for_metric']} does not match rule's target genotype."

        elif check_type == 'direction_conflict':
            if rule_name == 'APOE_Alz_Direction':
                calculated_or = get_nested_value(results, config['metric_path'], config['target_rsid_for_metric'], config['value_key'])
                
                if calculated_or is not None:
                    observed_effect_category = effect_utils.categorize_or(calculated_or)
                    expected_direction = config.get('expected_direction') # 'risk' or 'protective'

                    observed_is_risk = 'risk' in observed_effect_category.lower()
                    observed_is_protective = 'protective' in observed_effect_category.lower()

                    conflict = False
                    if expected_direction == 'risk' and observed_is_protective:
                        conflict = True
                    elif expected_direction == 'protective' and observed_is_risk:
                        conflict = True
                    
                    if conflict:
                        rule_outcome['status'] = 'DIRECTION_CONFLICT'
                        rule_outcome['details'] = (f"Observed OR {calculated_or:.2f} ({observed_effect_category}) "
                                                   f"conflicts with expected direction '{expected_direction}'.")
                    else:
                        rule_outcome['status'] = 'PASS'
                        rule_outcome['details'] = (f"Observed OR {calculated_or:.2f} ({observed_effect_category}) "
                                                   f"is consistent with expected direction '{expected_direction}'.")
                else:
                    rule_outcome['details'] = f"Calculated OR for {config['target_rsid_for_metric']} not found for direction check."
        
        # Fallback for other rule types or if rule_name doesn't match specific handlers above
        # This part handles the original APOE_Alz and C282Y_HFE_Penetrance if they were not refactored into check_types
        # However, with the refactor, this 'else' block might become less relevant if all rules have a 'check_type'.
        # For now, keeping the original APOE_Alz logic as a fallback if rule_name is 'APOE_Alz' and no check_type is defined.
        elif rule_name == 'APOE_Alz': # Original rule name as fallback
            apoe_e4_alleles = 0 
            # This is a simplified APOE check. Real APOE genotyping is complex.
            # This block is now mostly redundant if APOE_Alz_Value_Check is used.
            # Kept for safety if old rule name 'APOE_Alz' is still somehow used without 'check_type'.
            if results.get('disease_risk') and results['disease_risk'].get('neurological'):
                 for finding in results['disease_risk']['neurological']:
                    if finding.get('rsid') == 'rs429358':
                        rs429358_genotype = finding.get('genotype')
                        if rs429358_genotype:
                            apoe_e4_alleles = rs429358_genotype.count('C')
                        break
            if apoe_e4_alleles in config.get('expected_values', {}): # Check if expected_values exists
                expected_value = config['expected_values'][apoe_e4_alleles]
                calculated_value = get_nested_value(results, config['metric_path'], config['target_rsid_for_metric'], config['value_key'])
                if calculated_value is not None:
                    deviation = abs(calculated_value - expected_value) / expected_value if expected_value != 0 else float('inf')
                    if deviation <= 0.10:
                        rule_outcome['status'] = 'PASS'
                        rule_outcome['details'] = f"Fallback APOE_Alz: Calculated: {calculated_value:.2f}, Expected: {expected_value:.2f}"
                    else:
                        rule_outcome['status'] = 'CONCERN'
                        rule_outcome['details'] = f"Fallback APOE_Alz: Calculated: {calculated_value:.2f}, Expected: {expected_value:.2f}"
            # else: rule_outcome['details'] remains 'NOT_APPLICABLE' or specific error from above

        report.append(rule_outcome)
    return report
