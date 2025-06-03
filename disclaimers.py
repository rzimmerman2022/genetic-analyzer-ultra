BASE = "Research & educational use only – not medical advice. This analysis is based on current scientific understanding, which is constantly evolving. Genetic predispositions do not determine health outcomes, as environmental and lifestyle factors play significant roles."

PSYCH = """Psychological trait predictions (including cognitive abilities, personality, etc.) are based on statistical associations found in large population studies. These genetic scores typically explain a small fraction of the variance in these complex traits (often <15%). 
Environmental factors (e.g., upbringing, education, life experiences, personal choices, socioeconomic status) and gene-environment interactions have a much larger impact. 
These predictions should NOT be used to assess an individual's intelligence, character, potential, or worth. They are not deterministic and have very limited predictive power for any single individual."""

ANCESTRY = """Accuracy of genetic predictions (for diseases, traits, and drug responses) may be significantly reduced for individuals whose ancestry is not well-represented in the underlying research studies. 
Most large-scale genetic studies to date have predominantly involved participants of European ancestry. Therefore, findings may not generalize well to individuals of African, Asian, Hispanic/Latinx, Native American, Pacific Islander, or mixed ancestries. 
This is a known limitation in current genomics research, and efforts are underway to improve diversity in genetic studies."""

FUNCTIONAL_PREDICTION = """Computational predictions of a variant's functional impact (e.g., on protein function) are based on algorithms and models. While these tools are valuable for research, they are not definitive. 
Experimental validation is required to confirm the actual biological effect of a variant. These predictions should not be used for clinical decision-making without further confirmatory evidence."""

RARE_DISEASE = """Findings related to rare diseases or carrier status should be confirmed by clinical-grade genetic testing in a certified laboratory. 
Consult with a healthcare provider or genetic counselor to understand the implications of any such findings for yourself and your family members."""

PHARMACOGENOMICS = """Pharmacogenomic information can help predict how you might respond to certain medications. However, these are predictions, and actual drug response can be influenced by many other factors including other medications, diet, age, and overall health. 
Always discuss medication decisions with your healthcare provider. Do not change or stop any medication based solely on this genetic report."""

def build_disclaimer(analysis_type: str = None, ancestry_flag: str = 'EU', has_functional_predictions: bool = False, has_rare_disease_findings: bool = False, has_pharmacogenomics: bool = False) -> str:
    """
    Builds a contextual disclaimer string.
    analysis_type can be 'psychological', 'disease_risk', etc.
    ancestry_flag can be 'EU', 'AFR', 'ASN', 'AMR', 'MIX', 'OTH' etc.
    """
    parts = [BASE]

    if analysis_type == 'psychological_traits': # Match key used in AdvancedGeneticAnalyzer
        parts.append(PSYCH)

    if ancestry_flag and ancestry_flag.upper() != 'EU':
        parts.append(ANCESTRY)
    
    if has_functional_predictions:
        parts.append(FUNCTIONAL_PREDICTION)
    
    if has_rare_disease_findings:
        parts.append(RARE_DISEASE)

    if has_pharmacogenomics:
        parts.append(PHARMACOGENOMICS)
        
    return "\n\n".join(parts)
