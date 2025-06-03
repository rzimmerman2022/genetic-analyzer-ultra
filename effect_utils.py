EFFECT_BINS = [
    (0.9, 1.1, 'negligible'),
    (0.8, 0.9, 'small'), 
    (1.1, 1.25, 'small'),
    (0.67, 0.8, 'moderate'), 
    (1.25, 1.5, 'moderate'),
    (0.5, 0.67, 'large'), 
    (1.5, 2.0, 'large')
]

def categorize_or(or_val: float) -> str:
    """
    Categorizes an odds ratio (OR) into predefined effect size bins.
    """
    if or_val is None:
        return 'unknown'
        
    # Handle protective effects by inverting if OR < 1, then categorizing magnitude
    # This ensures that OR=0.5 (strong protective) is treated same magnitude as OR=2.0 (strong risk)
    effective_or = or_val
    if 0 < or_val < 1:
        effective_or = 1 / or_val
    elif or_val == 1: # Exactly 1 is negligible
        return 'negligible'
    elif or_val <= 0: # Non-positive ORs are problematic or indicate strong protection
        return 'very_large_protective_effect_or_error'


    for low, high, label in EFFECT_BINS:
        # Check against the effective_or for magnitude
        if low <= effective_or <= high:
            # If original OR was protective, append that info
            if 0 < or_val < 1 and label != 'negligible':
                return f"{label}_protective"
            return label
            
    # If effective_or is outside all defined bins, it's very large
    if 0 < or_val < 1 and label != 'negligible': # Check original or_val for direction
        return 'very_large_protective'
    return 'very_large_risk'
