"""
atpm_forced_ATP_drain.py

Run an ATPM (non-growth–associated ATP maintenance) sweep for the
Buchnera iLG240 model.

This forces the model to hydrolyse a minimum amount of ATP through
the ATP maintenance reaction (ATPM), without altering any reaction
stoichiometry in the BIOMASS reaction.

For each ATPM lower bound value (5…1000), the model is re-optimised
and the resulting biomass flux is recorded.
"""

import cobra
import numpy as np
import pandas as pd


def load_buchnera_model(path: str = "iLG240.xml") -> cobra.Model:
    """Load the Buchnera iLG240 model from an SBML file."""
    return cobra.io.read_sbml_model(path)


def main():

    # Load model
    model = load_buchnera_model("iLG240.xml")

    # Identify ATPM reaction
    atpm = None
    for r in model.reactions:
        if r.id.lower() == "atpm":
            atpm = r
            break

    if atpm is None:
        raise ValueError("ATPM reaction not found in the model!")

    # Store baseline ATPM lower bound
    baseline_lb = atpm.lower_bound

    # Forced ATP drain values
    loads = [0, 5, 10, 20, 40, 80, 160, 300, 600, 1000]

    results = []

    for lb in loads:
        m = model.copy()
        r = m.reactions.get_by_id(atpm.id)

        # Reset ATPM to baseline
        r.lower_bound = baseline_lb

        # Apply new forced ATP drain
        r.lower_bound = lb

        sol = m.optimize()
        growth = sol.objective_value if sol.status == "optimal" else 0.0

        results.append({"ATPM_forced": lb, "growth": growth})

    df = pd.DataFrame(results)
    df.to_csv("buchnera_ATPM_sweep.csv", index=False)
    print(df)


if __name__ == "__main__":
    main()
