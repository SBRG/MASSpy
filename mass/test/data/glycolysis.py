import mass
print("Version: " + mass.__version__)
from mass import MassMetabolite, MassReaction, MassModel

# Glycolysis Participants
glc__D_c = MassMetabolite(
    "glc__D_c",
    name="D-Glucose",
    formula="C6H12O6",
    charge=None,
    compartment="c")
g6p_c = MassMetabolite(
    "g6p_c",
    name="D-Glucose 6-phosphate",
    formula="C6H11O9P",
    charge=-2,
    compartment="c")
f6p_c = MassMetabolite(
    "f6p_c",
    name="D-Fructose 6-phosphate",
    formula="C6H11O9P",
    charge=-2,
    compartment="c")
fdp_c = MassMetabolite(
    "fdp_c",
    name="D-Fructose 1,6-bisphosphate",
    formula="C6H10O12P2",
    charge=-4,
    compartment="c")
g3p_c = MassMetabolite(
    "g3p_c",
    name="Glyceraldehyde 3-phosphate",
    formula="C3H5O6P",
    charge=-2,
    compartment="c")
dhap_c = MassMetabolite(
    "dhap_c",
    name="Dihydroxyacetone phosphate",
    formula="C3H5O6P",
    charge=-2,
    compartment="c")
_13dpg_c = MassMetabolite(
    "13dpg_c",
    name="3-Phospho-D-glyceroyl phosphate",
    formula="C3H4O10P2",
    charge=-4,
    compartment="c")
_3pg_c = MassMetabolite(
    "3pg_c",
    name="3-Phospho-D-glycerate",
    formula="C3H4O7P",
    charge=-3,
    compartment="c")
_2pg_c = MassMetabolite(
    "2pg_c",
    name="D-Glycerate 2-phosphate",
    formula="C3H4O7P",
    charge=-3,
    compartment="c")
pep_c = MassMetabolite(
    "pep_c",
    name="Phosphoenolpyruvate",
    formula="C3H2O6P",
    charge=-3,
    compartment="c")
pyr_c = MassMetabolite(
    "pyr_c",
    name="Pyruvate",
    formula="C3H3O3",
    charge=-1,
    compartment="c")
lac__L_c = MassMetabolite(
    "lac__L_c",
    name="L-Lactate",
    formula="C3H5O3",
    charge=-1,
    compartment="c")

# Co-Participants
atp_c = MassMetabolite(
    "atp_c",
    name="ATP",
    formula="C10H12N5O13P3",
    charge=-4,
    compartment="c")
adp_c = MassMetabolite(
    "adp_c",
    name="ADP",
    formula="C10H12N5O10P2",
    charge=-4,
    compartment="c")
amp_c = MassMetabolite(
    "amp_c",
    name="AMP",
    formula="C10H12N5O7P",
    charge=-2,
    compartment="c")
h_c = MassMetabolite(
    "h_c",
    name="H+",
    formula="H",
    charge=1,
    compartment="c")
nad_c = MassMetabolite(
    "nad_c",
    name="Nicotinamide adenine dinucleotide",
    formula="C21H26N7O14P2",
    charge=-1,
    compartment="c")
nadh_c = MassMetabolite(
    "nadh_c",
    name="Nicotinamide adenine dinucleotide - reduced",
    formula="C21H27N7O14P2",
    charge=-2,
    compartment="c")
pi_c = MassMetabolite(
    "pi_c",
    name="Phosphate",
    formula="HPO4",
    charge=-2,
    compartment="c")
h2o_c = MassMetabolite(
    "h20_c",
    name="H2O",
    formula="H2O",
    charge=0,
    compartment="c")

### D-Glucose to D-Glucose 6-phosphate
HEX1 = MassReaction(
    "HEX1",
    name="Hexokinase (D-glucose:ATP)",
    subsystem="Glycolysis")

HEX1.add_metabolites({
    glc__D_c: -1,
    atp_c: -1,
    adp_c: 1,
    g6p_c: 1,
    h_c: 1
})

### D-Glucose 6-phosphate to D-Fructose 6-phosphate
PGI = MassReaction(
    "PGI",
    name="Glucose-6-phosphate isomerase",
    subsystem="Glycolysis")

PGI.add_metabolites({
    g6p_c: -1,
    f6p_c: 1
})

### D-Fructose 6-phosphate to D-Fructose 1,6-bisphosphate
PFK = MassReaction(
    "PFK",
    name="Phosphofructokinase",
    subsystem="Glycolysis")

PFK.add_metabolites({
    f6p_c: -1,
    atp_c: -1,
    fdp_c: 1,
    adp_c: 1,
    h_c: 1
})

### D-Fructose 1,6-bisphosphate to
### Dihydroxyacetone phosphate + Glyceraldehyde 3-phosphate
FBA = MassReaction(
    "FBA",
    name="Fructose-bisphosphate aldolase",
    subsystem="Glyolysis")

FBA.add_metabolites({
    fdp_c: -1,
    dhap_c: 1,
    g3p_c: 1
})

### Dihydroxyacetone phosphate to Glyceraldehyde 3-phosphate
TPI = MassReaction(
    "TPI",
    name="Triose-phosphate isomerase",
    subsystem="Glycolysis")

TPI.add_metabolites({
    dhap_c: -1,
    g3p_c: 1
})

### Glyceraldehyde 3-phosphate to 3-Phospho-D-glyceroyl phosphate
GAPD = MassReaction(
    "GAPD",
    name="Glyceraldehyde-3-phosphate dehydrogenase",
    subsystem="Glycolysis")

GAPD.add_metabolites({
    g3p_c: -1,
    nad_c: -1,
    pi_c: -1,
    _13dpg_c: 1,
    h_c: 1,
    nadh_c: 1
})

### 3-Phospho-D-glycerate to 3-Phospho-D-glyceroyl phosphate
PGK = MassReaction(
    "PGK",
    name="Phosphoglycerate kinase",
    subsystem="Glycolysis")

PGK.add_metabolites({
    _13dpg_c: 1,
    adp_c: 1,
    _3pg_c: -1,
    atp_c: -1
})

### D-Glycerate 2-phosphate to 3-Phospho-D-glycerate
PGM = MassReaction(
    "PGM",
    name="Phosphoglycerate mutase",
    subsystem="Glycolysis")

PGM.add_metabolites({
    _3pg_c: 1,
    _2pg_c: -1
})

### D-Glycerate 2-phosphate to Phophoenolpyruvate
ENO = MassReaction(
    "ENO",
    name="Enolase",
    subsystem="Glycolysis")

ENO.add_metabolites({
    _2pg_c: -1,
    h2o_c: 1,
    pep_c: 1
})

### Phophoenolpyruvate to Pyruvate
PYK = MassReaction(
    "PYK",
    name="Pyruvate kinase",
    subsystem="Glycolysis")

PYK.add_metabolites({
    pep_c: -1,
    h_c: -1,
    adp_c: -1,
    atp_c: 1,
    pyr_c: 1
})

### L-lactate to Pyruvate
LDH_L = MassReaction(
    "LDH_L",
    name="L-lactate dehydrogenase",
    subsystem="Misc.")

LDH_L.add_metabolites({
    lac__L_c: 1,
    nad_c: 1,
    h_c: -1,
    nadh_c: -1,
    pyr_c: -1
})

### AMP + ATP <-> 2 ADP
ADK1 = MassReaction(
    "ADK1",
    name="Adenylate kinase",
    subsystem="Misc.")

ADK1.add_metabolites({
    amp_c: -1,
    atp_c: -1,
    adp_c: 2
})

### ATP + H2O <-> ADP + H + Pi
### ### Note: This is a pseudoreaction
ATPM = MassReaction(
    "ATPM",
    name="ATP maintenance requirement",
    subsystem="Misc.")

ATPM.add_metabolites({
    atp_c: -1,
    h2o_c: -1,
    adp_c: 1,
    h_c: 1,
    pi_c: 1
})

### NADH <-> NAD + H
### ### Note: BiGG reaction is different and more complex
NADH = MassReaction(
    "NADH",
    name="NADH dehydrogenase",
    subsystem="Misc.")

NADH.add_metabolites({
    nadh_c: -1,
    nad_c: 1,
    h_c: 1
})

# Exchange Reactions below

### (null) <-> glc__D_c
### ### Note: BiGG reaction is different
GLUIN = MassReaction(
    "GLUIN",
    name="D-Glucose exchange",
    subsystem="Transport/Exchange")

GLUIN.add_metabolites({
    glc__D_c: 1
})

### (null) <-> amp_c
### ### Source for AMP
### ### Note: BiGG reaction is different
AMPIN = MassReaction(
    "AMPIN",
    name="AMP exchange",
    subsystem="Transport/Exchange")

AMPIN.add_metabolites({
    amp_c: 1
})

### amp_c <-> (null)
### ### Sink for AMP
### ### Note: BiGG reaction does not exist
AMP = MassReaction(
    "AMP",
    name="AMP sink",
    subsystem="Transport/Exchange")

AMP.add_metabolites({
    amp_c: -1
})

### lac__L_c <-> (null)
### ### Note: BiGG reaction is different
LAC = MassReaction(
    "LAC",
    name="L-Lactate exchange",
    subsystem="Transport/Exchange")

LAC.add_metabolites({
    lac__L_c: -1
})

### pyr_c <-> (null)
### ### Note: BiGG reaction is different
PYR = MassReaction(
    "PYR",
    name="Pyruvate exchange",
    subsystem="Transport/Exchange")

PYR.add_metabolites({
    pyr_c: -1
})

### h_c <-> (null)
### ### Note: BiGG reaction is different
H = MassReaction(
    "H",
    name="H+ exchange",
    subsystem="Transport/Exchange")

H.add_metabolites({
    h_c: -1
})

### h2o_c <-> (null)
### ### Note: BiGG reaction is different
H2O = MassReaction(
    "H2O",
    name="H2O exchange",
    subsystem="Transport/Exchange")

H2O.add_metabolites({
    h2o_c: -1
})

### Generate Model object
glycolysis = MassModel("Glycolysis")

### Add reactions to model
### Note: This will add all associated metabolites as well
rxns_list = [HEX1, PGI, PFK, FBA, TPI, GAPD, PGK,
             PGM, ENO, PYK, LDH_L, ADK1, ATPM, NADH,
            GLUIN, AMPIN, AMP, LAC, PYR, H, H2O]

glycolysis.add_reactions(rxns_list)


### Set initial conditions
### All values in units of millimoles/Liter
glc__D_c.initial_condition = 1.0
g6p_c.initial_condition = 0.0486
f6p_c.initial_condition = 0.0198

fdp_c.initial_condition = 0.0146
g3p_c.initial_condition = 0.00728
dhap_c.initial_condition = 0.16

_13dpg_c.initial_condition = 0.000243
_3pg_c.initial_condition = 0.0773
_2pg_c.initial_condition = 0.0113

pep_c.initial_condition = 0.017
pyr_c.initial_condition = 0.060301
lac__L_c.initial_condition = 1.36

atp_c.initial_condition = 1.6
adp_c.initial_condition = 0.29
amp_c.initial_condition = 0.0867281

h_c.initial_condition = 0.0000899757
nad_c.initial_condition = 0.0589
nadh_c.initial_condition = 0.0301

pi_c.initial_condition = 2.5
h2o_c.initial_condition = 1.0

ic_dict = {}
for metab in glycolysis.metabolites:
    ic_dict[metab] = metab.ic

glycolysis.update_initial_conditions(ic_dict)

### Add gibbs_formation
### All values in units of ""
