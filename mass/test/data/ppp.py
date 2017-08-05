import mass
print("Version: " + mass.__version__)
from mass import MassMetabolite, MassReaction, MassModel

# PPP Participants
g6p_c = MassMetabolite(
    "g6p_c",
    name="D-Glucose 6-phosphate",
    formula="C6H11O9P",
    charge=-2,
    compartment="c")
_6pgl_c = MassMetabolite(
    "6pgl_c",
    name="6-Phospho-D-gluco-1,5-lactone",
    formula="C6H9O9P",
    charge=-2,
    compartment="c")
_6pgc_c = MassMetabolite(
    "6pgc_c",
    name="6-Phospho-D-gluconate",
    formula="C6H10O10P",
    charge=-3,
    compartment="c")
ru5p__D_c = MassMetabolite(
    "ru5p__D_c",
    name="D-Ribulose 5-phosphate",
    formula="C5H9O8P",
    charge=-2,
    compartment="c")
r5p_c = MassMetabolite(
    "r5p_c",
    name="Alpha-D-Ribose 5-phosphate",
    formula="C5H9O8P",
    charge=-2,
    compartment="c")
xu5p__D_c = MassMetabolite(
    "xu5p__D_c",
    name="D-Xylulose 5-phosphate",
    formula="C5H9O8P",
    charge=-2,
    compartment="c")
g3p_c = MassMetabolite(
    "g3p_c",
    name="Glyceraldehyde 3-phosphate",
    formula="C3H5O6P",
    charge=-2,
    compartment="c")
s7p_c = MassMetabolite(
    "s7p_c",
    name="Sedoheptulose 7-phosphate",
    formula="C7H13O10P",
    charge=-2,
    compartment="c")
f6p_c = MassMetabolite(
    "f6p_c",
    name="D-Fructose 6-phosphate",
    formula="C6H11O9P",
    charge=-2,
    compartment="c")
e4p_c = MassMetabolite(
    "e4p",
    name="D-Erythrose 4-phosphate",
    formula="C4H7O7P",
    charge=-2,
    compartment="c")

# Co-Participants
h_c = MassMetabolite(
    "h_c",
    name="H+",
    formula="H",
    charge=1,
    compartment="c")
nadp_c = MassMetabolite(
    "nadp_c",
    name="Nicotinamide adenine dinucleotide phosphate",
    formula="C21H25N7O17P3",
    charge=-3,
    compartment="c")
nadph_c = MassMetabolite(
    "nadph_c",
    name="Nicotinamide adenine dinucleotide phosphate - reduced",
    formula="C21H26N7O17P3",
    charge=-4,
    compartment="c")
h2o_c = MassMetabolite(
    "h20_c",
    name="H2O",
    formula="H2O",
    charge=0,
    compartment="c")
co2_c = MassMetabolite(
    "co2_c",
    name="CO2",
    formula="CO2",
    charge=0,
    compartment="c")
gthox_c = MassMetabolite(
    "gthox_c",
    name="Oxidized glutathione",
    formula="C20H30N6O12S2",
    charge=-2,
    compartment="c")
gthrd_c = MassMetabolite(
    "gthrd_c",
    name="Reduced glutathione",
    formula="C10H15N3O6S,C10H16N3O6S",
    charge=-1,
    compartment="c")

# Oxidative Phase

### D-Glucose 6-phosphate to 6-Phospho-D-gluco-1,5-lactone
G6PDH2r = MassReaction(
    "G6PDH2r",
    name="Glucose 6-phosphate dehydrogenase",
    subsystem="Pentose Phosphate Pathway")

G6PDH2r.add_metabolites({
    g6p_c: -1,
    nadp_c: -1,
    _6pgl_c: 1,
    nadph_c: 1,
    h_c: 1
})

### 6-Phospho-D-gluco-1,5-lactone to 6-Phospho-D-gluconate
PGL = MassReaction(
    "PGI",
    name="6-phosphogluconolactonase",
    subsystem="Pentose Phosphate Pathway")

PGL.add_metabolites({
    _6pgl_c: -1,
    h2o_c: -1,
    _6pgc_c: 1,
    h_c: 1
})

### 6-Phospho-D-gluconate to D-Ribulose 5-phosphate
GND = MassReaction(
    "GND",
    name="Phosphogluconate dehydrogenase",
    subsystem="Pentose Phosphate Pathway")

GND.add_metabolites({
    _6pgc_c: -1,
    nadp_c: -1,
    nadph_c: 1,
    co2_c: 1,
    ru5p__D_c: 1
})


# Non-oxidative Phase

### Alpha-D-Ribose 5-phosphate to D-Ribulose 5-phosphate
RPI = MassReaction(
    "RPI",
    name="Ribulose 5-Phosphate Isomerase",
    subsystem="Glyolysis")

RPI.add_metabolites({
    r5p_c: -1,
    ru5p__D_c: 1
})

### D-Ribulose 5-phosphate to D-Xylulose 5-phosphate
RPE = MassReaction(
    "RPE",
    name="Ribulose 5-phosphate 3-epimerase",
    subsystem="Pentose Phosphate Pathway")

RPE.add_metabolites({
    ru5p__D_c: -1,
    xu5p__D_c: 1
})

### Alpha-D-Ribose 5-phosphate + D-Xylulose 5-phosphate to
### Glyceraldehyde 3-phosphate + Sedoheptulose 7-phosphate
TKT1 = MassReaction(
    "TKT1",
    name="Transketolase",
    subsystem="Pentose Phosphate Pathway")

TKT1.add_metabolites({
    r5p_c: -1,
    xu5p__D_c: -1,
    g3p_c: 1,
    s7p_c: 1
})

### Glyceraldehyde 3-phosphate + Sedoheptulose 7-phosphate to
### D-Erythrose 4-phosphate + D-Fructose 6-phosphate
TALA = MassReaction(
    "TALA",
    name="Transaldolase",
    subsystem="Pentose Phosphate Pathway")

TALA.add_metabolites({
    g3p_c: -1,
    s7p_c: -1,
    e4p_c: 1,
    f6p_c: 1
})

### D-Erythrose 4-phosphate + D-Xylulose 5-phosphate to
### D-Fructose 6-phosphate + Glyceraldehyde 3-phosphate
TKT2 = MassReaction(
    "TKT2",
    name="Transketolase",
    subsystem="Pentose Phosphate Pathway")

TKT2.add_metabolites({
    e4p_c: -1,
    xu5p__D_c: -1,
    f6p_c: 1,
    g3p_c: 1
})

# Misc.

### GSSH (aka gthox) to 2 GSH (aka gthrd)
GTHOr = MassReaction(
    "GTHOr",
    name="Glutathione oxidoreductase",
    subsystem="Misc.")

GTHOr.add_metabolites({
    gthox_c: -1,
    h_c: -1,
    nadph_c: -1,
    gthrd_c: 2,
    nadp_c: 1
})

### GSSG + 2H <-> 2GSH
### ### Note: BiGG reaction is different and more complex
GSHR = MassReaction(
    "GSHR",
    name="Glutathione-disulfide reductase",
    subsystem="Misc.")

GSHR.add_metabolites({
    gthox_c: -1,
    h_c: -2,
    gthrd_c: 2
})

# Exchange Reactions

### co2_c <-> (null)
### ### Note: BiGG reaction is different
CO2 = MassReaction(
    "CO2",
    name="CO2 exchange",
    subsystem="Transport/Exchange")

CO2.add_metabolites({
    co2_c: -1
})

### g6p_c <-> (null)
### ### Note: BiGG reaction is different
EX_G6P_C = MassReaction(
    "EX_G6P_C",
    name="D-Glucose 6-phosphate exchange",
    subsystem="Tranpsort/Exchange")

EX_G6P_C.add_metabolites({
    g6p_c: -1
})

### f6p_c <-> (null)
### ### Sink for D-Fructose 6-phosphate
### ### Note: BiGG reaction is different
EX_F6P_C = MassReaction(
    "EX_F6P_C",
    name="D-fructose 6-phosphate exchange",
    subsystem="Transport/Exchange")

EX_F6P_C.add_metabolites({
    f6p_c: -1
})

### r5p_c <-> (null)
### ### Sink for Alpha-D-Ribose 5-phosphate
### ### Note: BiGG reaction is different
EX_R5P_C = MassReaction(
    "EX_R5P_C",
    name="Alpha-D-Ribose 5-phosphate exchange",
    subsystem="Transport/Exchange")

EX_R5P_C.add_metabolites({
    r5p_c: -1
})

### h_c <-> (null)
### ### Note: BiGG reaction is different
EX_H_C = MassReaction(
    "EX_H_C",
    name="H+ exchange",
    subsystem="Transport/Exchange")

EX_H_C.add_metabolites({
    h_c: -1
})

### h2o_c <-> (null)
### ### Note: BiGG reaction is different
EX_H2O_C = MassReaction(
    "EX_H2O_C",
    name="H2O exchange",
    subsystem="Transport/Exchange")

EX_H2O_C.add_metabolites({
    h2o_c: -1
})

### g3p_c <-> (null)
### ### Sink for Glyceraldehyde 3-phosphate
### ### Note: BiGG reaction not found
EX_G3P_C = MassReaction(
    "EX_G3P_C",
    name="Glyceraldehyde 3-phosphate exchange",
    subsystem="Transport/Exchange")

EX_G3P_C.add_metabolites({
    g3p_c: -1
})

### Generate Model object
ppp = MassModel("Pentose_Phosphate_Pathway")

### Add reactions to model
### Note: This will add all associated metabolites as well
# FIX WHEN SOURCES AND SINKS ADDED
rxns_list = [G6PDH2r, PGL, GND, RPI, RPE,
             TKT1, TALA, TKT2, GTHOr, GSHR,
            CO2, EX_G6P_C, EX_F6P_C, EX_R5P_C,
            EX_H_C, EX_H2O_C, EX_G3P_C]

ppp.add_reactions(rxns_list)


### Set initial conditions
### All values in units of millimoles/Liter
g6p_c.initial_condition = 0.0486
f6p_c.initial_condition = 0.0198
g3p_c.initial_condition = 0.00728

_6pgl_c.initial_condition = 0.00175424
_6pgc_c.initial_condition = 0.0374753
ru5p__D_c.initial_condition = 0.00493679

xu5p__D_c.initial_condition = 0.0147842
r5p_c.initial_condition = 0.0126689
s7p_c.initial_condition = 0.023988

e4p_c.initial_condition = 0.00507507
nadp_c.initial_condition = 0.0002
nadph_c.initial_condition = 0.0658

gthrd_c.initial_condition = 3.2
gthox_c.initial_condition = 0.12
co2_c.initial_condition = 1.

h_c.initial_condition = 0.0000714957
h2o_c.initial_condition = 0.999998

ic_dict = {}
for metab in ppp.metabolites:
    ic_dict[metab] = metab.ic

ppp.update_initial_conditions(ic_dict)

### Add gibbs_formation
### All values in units of ""
