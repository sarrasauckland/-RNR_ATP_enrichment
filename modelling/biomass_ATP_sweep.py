"""
biomass_ATP_sweep.py

Run a biomass ATP-cost sweep for the Buchnera iLG240 model.
For each extra ATP load (ΔATP), we modify the BIOMASS reaction
to require (20 + ΔATP) ATP per unit biomass and record the
resulting growth rate.

This script expects iLG240.xml to be in the same folder where you run it.
"""

import cobra
import numpy as np
import pandas as pd


def load_buchnera_model(path: str = "iLG240.xml") -> cobra.Model:
    """Load the Buchnera iLG240 model from an SBML file."""
    model = cobra.io.read_sbml_model(path)
    return model


def get_original_atp_coeffs(biomass_rxn: cobra.Reaction, model: cobra.Model) -> dict:
    """
    Return the original ATP/ADP/Pi/H2O/H+ stoichiometric coefficients
    from the biomass reaction so we can always reset to baseline.
    """
    mets = model.metabolites
    return {
        "atp_c": biomass_rxn.metabolites[mets.atp_c],
        "h2o_c": biomass_rxn.metabolites[mets.h2o_c],
        "adp_c": biomass_rxn.metabolites[mets.adp_c],
        "pi_c": biomass_rxn.metabolites[mets.pi_c],
        "h_c": biomass_rxn.metabolites[mets.h_c],
    }


def set_atp_load(model: cobra.Model, orig: dict, delta: float) -> None:
    """
    Reset the BIOMASS reaction to its original ATP-related stoichiometry,
    then apply an extra ATP cost of 'delta' (ΔATP).

    This guarantees that ATP changes are NOT cumulative between runs.
    """
    b = model.reactions.get_by_id("BIOMASS")
    mets = model.metabolites

    # 1. Reset to original coefficients
    b.add_metabolites(
        {
            mets.atp_c: -(b.metabolites[mets.atp_c] - orig["atp_c"]),
            mets.h2o_c: -(b.metabolites[mets.h2o_c] - orig["h2o_c"]),
            mets.adp_c: -(b.metabolites[mets.adp_c] - orig["adp_c"]),
            mets.pi_c: -(b.metabolites[mets.pi_c] - orig["pi_c"]),
            mets.h_c: -(b.metabolites[mets.h_c] - orig["h_c"]),
        },
        combine=True,
    )

    # 2. Apply ΔATP, if any
    if delta != 0:
        b.add_metabolites(
            {
                mets.atp_c: -delta,
                mets.h2o_c: -delta,
                mets.adp_c: delta,
                mets.pi_c: delta,
                mets.h_c: delta,
            },
            combine=True,
        )


def main():
    # Load model and identify BIOMASS
    model = load_buchnera_model("iLG240.xml")
    biomass = model.reactions.get_by_id("BIOMASS")

    # Store original ATP-related stoichiometry
    orig = get_original_atp_coeffs(biomass, model)

    # ΔATP values to test
    loads = [0, 5, 10, 20, 50, 80, 100, 160, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000]

    results = []

    for d in loads:
        # Work on a fresh copy each time for safety
        m_copy = model.copy()
        set_atp_load(m_copy, orig, d)
        sol = m_copy.optimize()
        growth = sol.objective_value if sol.status == "optimal" else 0.0
        results.append({"delta_ATP": d, "growth": growth})

    df = pd.DataFrame(results)
    df.to_csv("buchnera_biomass_ATP_sweep.csv", index=False)
    print(df)


if __name__ == "__main__":
    main()
