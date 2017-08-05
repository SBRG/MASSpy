import mass
print("Version: " + mass.__version__)
from mass import MassMetabolite, MassReaction, MassModel

# AMPS Participants
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
adn_c = MassMetabolite(
    "adn_c",
    name="Adenosine",
    formula="C10H13N5O4",
    charge=0,
    compartment="c")
ade_c = MassMetabolite(
    "ade_c",
    name="Adenine",
    formula="C5H5N5",
    charge=0,
    compartment="c")
imp_c = MassMetabolite(
    "imp_c",
    name="IMP",
    formula="C10H11N4O8P",
    charge=-2,
    compartment="c")
ins_c = MassMetabolite(
    "ins_c",
    name="Inosince",
    formula="C10H12N4O5",
    charge=0,
    compartment="c")
hxan_c = MassMetabolite(
    "hxan_c",
    name="Hypoxanthine",
    formula="C5H4N4O",
    charge=0,
    compartment="c")
r5p_c = MassMetabolite(
    "r5p_c",
    name="Alpha-D-Ribose 5-phosphate",
    formula="C5H9O8P",
    charge=-2,
    compartment="c")
r1p_c = MassMetabolite(
    "r1p_c",
    name="Alpha-D-Ribose 1-phosphate",
    formula="C5H9O8P",
    charge=-2,
    compartment="c")
prpp_c = MassMetabolite(
    "prpp_c",
    name="5-Phospho-alpha-D-ribose 1-diphosphate",
    formula="C5H8O14P3",
    charge=-5,
    compartment="c")

# Co-Participants
h_c = MassMetabolite(
    "h_c",
    name="H+",
    formula="H",
    charge=1,
    compartment="c")
h2o_c = MassMetabolite(
    "h20_c",
    name="H2O",
    formula="H2O",
    charge=0,
    compartment="c")
pi_c = MassMetabolite(
    "pi_c",
    name="Phosphate",
    formula="HPO4",
    charge=-2,
    compartment="c")
nh3_c = MassMetabolite(
    "nh3_c",
    name="Ammonia",
    formula="H3N",
    charge=0,
    compartment="c")

### Adenosine + ATP to ADP + AMP
### ### Note: BiGG reaction is different and more complex
AK = MassReaction(
    "AK",
    name="Adenosine kinase",
    subsystem="AMP Salvage Network")

AK.add_metabolites({
    adn_c: -1,
    atp_c: -1,
    adp_c: 1,
    amp_c: 1
})

### AMP to Adenosine
NTD7 = MassReaction(
    "NTD7",
    name="5'-nucleotidase (AMP)",
    subsystem="AMP Salvage Network")

NTD7.add_metabolites({
    amp_c: -1,
    h2o_c: -1,
    adn_c: 1,
    pi_c: 1
})

### AMP to IMP
### ### Note: BiGG reaction is different and more complex
AMPDA = MassReaction(
    "AMPDA",
    name="Adenosine monophosphate deaminase",
    subsystem="AMP Salvage Network")

AMPDA.add_metabolites({
    amp_c: -1,
    h2o_c: -1,
    imp_c: 1,
    nh3_c: 1
})

### IMP to Inosine
NTD11 = MassReaction(
    "NTD11",
    name="5'-nucleotidase (IMP)",
    subsystem="AMP Salvage Network")

NTD11.add_metabolites({
    imp_c: -1,
    h2o_c: -1,
    ins_c: 1,
    pi_c: 1
})

### Adenosine to Inosine
### ### Note: BiGG reaction is different and more complex
ADA = MassReaction(
    "ADA",
    name="Adenosine deaminase",
    subsystem="AMP Salvage Network")

ADA.add_metabolites({
    adn_c: -1,
    h2o_c: -1,
    ins_c: 1,
    nh3_c: 1
})

### Inosine to Hypoxanthine + Alpha-D-Ribose 1-phosphate
PUNP5 = MassReaction(
    "PUNP5",
    name="Purine-nucleoside phosphorylase (Inosine)",
    subsystem="AMP Salvage Network")

PUNP5.add_metabolites({
    ins_c: -1,
    pi_c: -1,
    hxan_c: 1,
    r1p_c: 1
})

### Alpha-D-Ribose 1-phosphate to Alpha-D-Ribose 5-phosphate
PPM = MassReaction(
    "PPM",
    name="Phosphopentomutase",
    subsystem="AMP Salvage Network")

PPM.add_metabolites({
    r1p_c: -1,
    r5p_c: 1
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

### 2ATP + Alpha-D-Ribose 5-phosphate to
### 2ADP + H + 5-Phospho-alpha-D-ribose 1-diphosphate
### ### Note: BiGG reaction is different and more complex
PRPPSYN = MassReaction(
    "PRPPSYN",
    name="Phosphoribosylpyrophosphate synthetase",
    subsystem="AMP Salvage Network")

PRPPSYN.add_metabolites({
    atp_c: -2,
    r5p_c: -1,
    adp_c: 2,
    h_c: 1,
    prpp_c: 1
})

### Adenosine + H2O + 5-Phospho-alpha-D-ribose 1-diphosphate to
### AMP + H + 2Phosphate
### ### Note: BiGG reaction is different and more complex
ADPRT = MassReaction(
    "ADPRT",
    name="Adenine phosphoribosyltransferase",
    subsystem="AMP Salvage Network")

ADPRT.add_metabolites({
    ade_c: -1,
    h2o_c: -1,
    prpp_c: -1,
    amp_c: 1,
    h_c: 1,
    pi_c: 2
})

# Exchange Reactions

### amp_c <-> (null)
### ### Note: BiGG reaction is different
AMP = MassReaction(
    "AMP",
    name="AMP exchange",
    subsystem="Transport/Exchange")

AMP.add_metabolites({
    amp_c: -1
})

### adn_c <-> (null)
### ### Note: BiGG reaction is different
ADO = MassReaction(
    "ADO",
    name="Adenosine exchange",
    subsystem="Transport/Exchange")

ADO.add_metabolites({
    adn_c: -1
})

### ade_c <-> (null)
### ### Note: BiGG reaction is different
ADE = MassReaction(
    "ADE",
    name="Adenine exchange",
    subsystem="Tranpsort/Exchange")

ADE.add_metabolites({
    ade_c: -1
})

### ins_c <-> (null)
### ### Note: BiGG reaction is different
INO = MassReaction(
    "INO",
    name="Inosine exchange",
    subsystem="Tranpsort/Exchange")

INO.add_metabolites({
    ins_c: -1
})

### hxan_c <-> (null)
### ### Note: BiGG reaction is different
HYP = MassReaction(
    "HYP",
    name="Hypoxanthine exchange",
    subsystem="Tranpsort/Exchange")

HYP.add_metabolites({
    hxan_c: -1
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

### pi_c <-> (null)
### ### Note: BiGG reaction is different
PHOS = MassReaction(
    "PHOS",
    name="Phosphate exchange",
    subsystem="Transport/Exchange")

PHOS.add_metabolites({
    pi_c: -1
})

### nh3_c <-> (null)
### ### Note: BiGG reaction does not exist
NH3 = MassReaction(
    "NH3",
    name="NH3 exchange",
    subsystem="Transport/Exchange")

NH3.add_metabolites({
    nh3_c: -1
})

### Generate Model object
amps = MassModel("AMPSalvageNetwork")

### Add reactions to model
### Note: This will add all associated metabolites as well
# FIX WHEN SOURCES AND SINKS ADDED
rxns_list = [AK, NTD7, AMPDA, NTD11, ADA, 
             PUNP5, PPM, ATPM, PRPPSYN, ADPRT, 
             ADO, ADE, INO, HYP, AMP, H, H2O, PHOS, NH3]

amps.add_reactions(rxns_list)


### Set initial conditions
### All values in units of millimoles/Liter

atp_c.ic = 1.6
adp_c.ic = 0.29
amp_c.ic = 0.0867281

adn_c.ic = 0.0012
ade_c.ic = 0.001

imp_c.ic = 0.01
ins_c.ic = 0.001

r5p_c.ic = 0.00494
r1p_c.ic = 0.06

hxan_c.ic = 0.002
prpp_c.ic = 0.005

h_c.ic = 0.0000630957
h2o_c.ic = 1.

pi_c.ic = 2.5
nh3_c.ic = 0.091

ic_dict = {}
for metab in amps.metabolites:
    ic_dict[metab] = metab.ic

amps.update_initial_conditions(ic_dict)

### Add gibbs_formation
### All values in units of ""