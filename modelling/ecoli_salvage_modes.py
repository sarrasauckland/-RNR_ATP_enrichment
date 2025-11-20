#!/usr/bin/env python3
"""
E. coli ATP Stress Simulations with Purine Salvage
--------------------------------------------------

This script extends the WT minimal-medium ATP stress analysis
by enabling purine salvage import at different capacities.

Simulations performed:
 - ATP Δ (extra ATP cost per unit biomass)
 - ATPM (minimum ATP hydrolysis for maintenance)
 - Under two salvage modes: LOW salvage and HIGH salvage

Output:
 - results_low_salvage_atpdelta
 - results_high_salvage_atpdelta
 - results_low_salvage_atpm
 - results_high_salvage_atpm
"""

import cobra

# ------------------------------------------------------------
# Load and configure WT E. coli iML1515
# ------------------------------------------------------------

model_ec = cobra.io.read_sbml_model("iML1515.xml")

# Set biomass objective
biomass = model_ec.reactions.get_by_id("BIOMASS_Ec_iML1515_core_75p37M")
model_ec.objective = biomass

# Close all exchange reactions
for ex in model_ec.exchanges:
    ex.lower_bound = 0.0

# Minimal medium components
minimal_media = [
    "EX_glc__D_e", "EX_nh4_e", "EX_pi_e", "EX_so4_e",
    "EX_h2o_e", "EX_h_e", "EX_k_e", "EX_na1_e", "EX_cl_e",
    "EX_mg2_e", "EX_ca2_e", "EX_fe2_e", "EX_mn2_e", "EX_zn2_e",
    "EX_cobalt2_e", "EX_mobd_e", "EX_cu2_e", "EX_ni2_e", "EX_o2_e"
]

# Open minimal components
for ex_id in minimal_media:
    if ex_id in model_ec.reactions:
        model_ec.reactions.get_by_id(ex_id).lower_bound = -1000

# Limit glucose uptake
model_ec.reactions.EX_glc__D_e.lower_bound = -10


# ------------------------------------------------------------
# Utility functions
# ------------------------------------------------------------

def set_atp_delta(model, delta_atp):
    """
    Injects an increased ATP cost into the biomass equation.
    delta_atp = number of ATP required per unit biomass.
    """
    biomass_rxn = model.reactions.get_by_id("BIOMASS_Ec_iML1515_core_75p37M")
    biomass_rxn.add_metabolites({
        model.metabolites.atp_c: -delta_atp,
        model.metabolites.adp_c: +delta_atp,
        model.metabolites.pi_c: +delta_atp,
        model.metabolites.h_c: +delta_atp
    }, combine=True)


def set_atpm(model, value):
    """
    Apply an ATP maintenance (ATPM) lower bound.
    """
    if "ATPM" in model.reactions:
        model.reactions.ATPM.lower_bound = value


def set_salvage(model, flux):
    """
    Open purine salvage import pathways at defined flux.
    """
    salvage_imports = ["EX_ade_e", "EX_hxan_e", "EX_gua_e", "EX_xan_e"]
    for ex_id in salvage_imports:
        if ex_id in model.reactions:
            model.reactions.get_by_id(ex_id).lower_bound = -flux


# ------------------------------------------------------------
# Salvage Simulation Parameters
# ------------------------------------------------------------

LOW_SALVAGE = 1        # extremely limited salvage import
HIGH_SALVAGE = 1000    # unconstrained salvage

# ATPΔ stress levels
ATP_DELTA_VALUES = [0, 5, 10, 20, 50, 80, 100, 160, 200, 250,
                    300, 400, 500, 600, 700, 800, 900, 1000]

# ATPM stress levels
ATPM_VALUES = [0, 10, 20, 30, 50, 80, 100, 150, 200, 250, 300]


# ------------------------------------------------------------
# Run ATPΔ under LOW salvage
# ------------------------------------------------------------

results_low_salvage_atpdelta = {}

for d in ATP_DELTA_VALUES:
    m = model_ec.copy()
    set_salvage(m, LOW_SALVAGE)
    set_atp_delta(m, d)
    sol = m.optimize()
    results_low_salvage_atpdelta[d] = sol.objective_value


# ------------------------------------------------------------
# Run ATPΔ under HIGH salvage
# ------------------------------------------------------------

results_high_salvage_atpdelta = {}

for d in ATP_DELTA_VALUES:
    m = model_ec.copy()
    set_salvage(m, HIGH_SALVAGE)
    set_atp_delta(m, d)
    sol = m.optimize()
    results_high_salvage_atpdelta[d] = sol.objective_value


# ------------------------------------------------------------
# Run ATPM under LOW salvage
# ------------------------------------------------------------

results_low_salvage_atpm = {}

for a in ATPM_VALUES:
    m = model_ec.copy()
    set_salvage(m, LOW_SALVAGE)
    set_atpm(m, a)
    sol = m.optimize()
    results_low_salvage_atpm[a] = sol.objective_value


# ------------------------------------------------------------
# Run ATPM under HIGH salvage
# ------------------------------------------------------------

results_high_salvage_atpm = {}

for a in ATPM_VALUES:
    m = model_ec.copy()
    set_salvage(m, HIGH_SALVAGE)
    set_atpm(m, a)
    sol = m.optimize()
    results_high_salvage_atpm[a] = sol.objective_value


# ------------------------------------------------------------
# Print summary to console (or save to file)
# ------------------------------------------------------------

print("\n=== ATP Δ under LOW salvage ===")
print(results_low_salvage_atpdelta)

print("\n=== ATP Δ under HIGH salvage ===")
print(results_high_salvage_atpdelta)

print("\n=== ATPM under LOW salvage ===")
print(results_low_salvage_atpm)

print("\n=== ATPM under HIGH salvage ===")
print(results_high_salvage_atpm)
