"""Validation utilities for :mod:`genetic_analyzer_ultra`.

This module provides a very small validation harness used by the unit tests.
It checks a few expected effect directions for key variants and reports
potential conflicts.  The implementation here is intentionally lightweight and
only covers what is necessary for the test-suite.
"""

from __future__ import annotations

from typing import Any, Dict, List


# Minimal lookup of expected effect direction.  In a real application this would
# likely come from a curated database.  For the tests we only care about APOE
# rs429358 where a risk effect (>1) is expected.
EXPECTED_DIRECTIONS = {
    "rs429358": "risk",
}


def validate(results: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Validate analysis results.

    Parameters
    ----------
    results:
        The results dictionary produced by ``AdvancedGeneticAnalyzer``.

    Returns
    -------
    List[Dict[str, Any]]
        A list of validation findings.  Each finding is a dictionary with at
        least ``rule_name``, ``status`` and ``details`` keys.
    """

    report: List[Dict[str, Any]] = []

    disease_risk = results.get("disease_risk", {})
    for findings in disease_risk.values():
        for item in findings:
            rsid = item.get("rsid")
            expected = EXPECTED_DIRECTIONS.get(rsid)
            if expected is None:
                continue

            rr = item.get("relative_risk")
            if rr is None:
                continue

            if expected == "risk" and rr < 1:
                report.append(
                    {
                        "rule_name": "APOE_Alz_Direction",
                        "status": "DIRECTION_CONFLICT",
                        "details": f"{rsid} shows protective effect (OR {rr}) but expected risk.",
                    }
                )
            elif expected == "protective" and rr > 1:
                report.append(
                    {
                        "rule_name": "APOE_Alz_Direction",
                        "status": "DIRECTION_CONFLICT",
                        "details": f"{rsid} shows risk effect (OR {rr}) but expected protective.",
                    }
                )

    return report
