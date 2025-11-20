import cobra

# ---------------------------------------------------
# 1. Load E. coli model
# ---------------------------------------------------
model = cobra.io.read_sbml_model("iML1515.xml")

# Use the core biomass reaction
biomass = model.reactions.get_by_id("BIOMASS_Ec_iML1515_core_75p37M")
model.objective = biomass

# ---------------------------------------------------
# 2. Set minimal medium (model-specific)
# ---------------------------------------------------
# Close all imports first
for ex in model.exchanges:
    ex.lower_bound = 0.0

minimal = [
    "EX_glc__D_e",
    "EX_nh4_e",
    "EX_pi_e",
    "EX_so4_e",
    "EX_h2o_e",
    "EX_h_e",
    "EX_k_e",
    "EX_na1_e",
    "EX_cl_e",
    "EX_mg2_e",
    "EX_ca2_e",
    "EX_fe2_e",
    "EX_mn2_e",
    "EX_zn2_e",
    "EX_cobalt2_e",
    "EX_mobd_e",
    "EX_cu2_e",
    "EX_ni2_e",
    "EX_o2_e"
]

# Apply minimal medium
for ex_id in minimal:
    if ex_id in model.reactions:
        model.reactions.get_by_id(ex_id).lower_bound = -1000

# Limit glucose
model.reactions.EX_glc__D_e.lower_bound = -10

# Baseline growth
print("Baseline growth:", model.optimize().objective_value)

# ---------------------------------------------------
# 3. BIOMASS ATP HYDROLYSIS function
# ---------------------------------------------------
BASE_ATP  = 75.55223
BASE_H2O  = 70.028756
BASE_ADP  = 75.37723
BASE_PI   = 75.37323
BASE_H    = 75.37723

def set_atp_load(model, delta_atp):
    m = model.copy()
    b = m.reactions.get_by_id("BIOMASS_Ec_iML1515_core_75p37M")

    # Remove old ATP terms
    for met_id in ["atp_c", "adp_c", "pi_c", "h2o_c", "h_c"]:
        met = m.metabolites.get_by_id(met_id)
        b.add_metabolites({met: 0}, combine=False)

    # Add new ATP cost
    b.add_metabolites({
        m.metabolites.atp_c: -(BASE_ATP + delta_atp),
        m.metabolites.h2o_c: -(BASE_H2O + delta_atp),
        m.metabolites.adp_c: +(BASE_ADP + delta_atp),
        m.metabolites.pi_c:  +(BASE_PI + delta_atp),
        m.metabolites.h_c:   +(BASE_H  + delta_atp),
    }, combine=True)

    return m.optimize().objective_value

# ---------------------------------------------------
# 4. Run ΔATP series
# ---------------------------------------------------
deltas = [0,5,10,20,50,80,100,160,200,250,300,400,500,600,700,800,900,1000]
results = {}

for d in deltas:
    results[d] = set_atp_load(model, d)

print("\nBIOMASS ATP HYDROLYSIS results (E. coli):")
print(results)
