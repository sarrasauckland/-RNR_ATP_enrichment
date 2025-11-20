import cobra

# ---------------------------------------------------
# 1. Load E. coli model
# ---------------------------------------------------
model = cobra.io.read_sbml_model("iML1515.xml")

# Core biomass
biomass = model.reactions.get_by_id("BIOMASS_Ec_iML1515_core_75p37M")
model.objective = biomass

# ---------------------------------------------------
# 2. Set minimal medium (same as other script)
# ---------------------------------------------------
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

for ex_id in minimal:
    if ex_id in model.reactions:
        model.reactions.get_by_id(ex_id).lower_bound = -1000

model.reactions.EX_glc__D_e.lower_bound = -10

print("Baseline growth:", model.optimize().objective_value)

# ---------------------------------------------------
# 3. ATPM drain analysis
# ---------------------------------------------------
demands = [0,5,10,20,40,80,160,300,600,1000]
results = {}

for d in demands:
    m = model.copy()
    atpm = m.reactions.get_by_id("ATPM")
    atpm.lower_bound = d
    sol = m.optimize()
    results[d] = (sol.status, sol.objective_value)

print("\nATPM minimal medium results (E. coli):")
print(results)
