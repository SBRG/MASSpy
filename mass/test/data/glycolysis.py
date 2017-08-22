import mass
from mass import MassMetabolite, MassReaction, MassModel

### Generate Metabolites ###

# Glycolysis Participants
glc__D = MassMetabolite(
    id="glc__D",
    name="D-Glucose",
    formula="C6H12O6",
    charge=0,
    compartment="c"
)
g6p = MassMetabolite(
    id="g6p",
    name="D-Glucose 6-phosphate",
    formula="C6H11O9P",
    charge=-2,
    compartment="c"
)
f6p = MassMetabolite(
    id="f6p",
    name="D-Fructose 6-phosphate",
    formula="C6H11O9P",
    charge=-2,
    compartment="c"
)
fdp = MassMetabolite(
    id="fdp",
    name="D-Fructose 1,6-bisphosphate",
    formula="C6H10O12P2",
    charge=-4,
    compartment="c"
)
g3p = MassMetabolite(
    id="g3p",
    name="Glyceraldehyde 3-phosphate",
    formula="C3H5O6P",
    charge=-2,
    compartment="c"
)
dhap = MassMetabolite(
    id="dhap",
    name="Dihydroxyacetone phosphate",
    formula="C3H5O6P",
    charge=-2,
    compartment="c"
)
_13dpg = MassMetabolite(
    id="13dpg",
    name="3-Phospho-D-glyceroyl phosphate",
    formula="C3H4O10P2",
    charge=-4,
    compartment="c"
)
_3pg = MassMetabolite(
    id="3pg",
    name="3-Phospho-D-glycerate",
    formula="C3H4O7P",
    charge=-3,
    compartment="c"
)
_2pg = MassMetabolite(
    id="2pg",
    name="D-Glycerate 2-phosphate",
    formula="C3H4O7P",
    charge=-3,
    compartment="c"
)
pep = MassMetabolite(
    id="pep",
    name="Phosphoenolpyruvate",
    formula="C3H2O6P",
    charge=-3,
    compartment="c"
)
pyr = MassMetabolite(
    id="pyr",
    name="Pyruvate",
    formula="C3H3O3",
    charge=-1,
    compartment="c"
)
lac__L = MassMetabolite(
    id="lac__L",
    name="L-Lactate",
    formula="C3H5O3",
    charge=-1,
    compartment="c"
)

# Co-Participants
atp = MassMetabolite(
    id="atp",
    name="ATP",
    formula="C10H12N5O13P3",
    charge=-4,
    compartment="c"
)
adp = MassMetabolite(
    "adp",
    name="ADP",
    formula="C10H12N5O10P2",
    charge=-3,
    compartment="c"
)
amp = MassMetabolite(
    id="amp",
    name="AMP",
    formula="C10H12N5O7P",
    charge=-2,
    compartment="c"
)
h = MassMetabolite(
    id="h",
    name="H+",
    formula="H",
    charge=1,
    compartment="c"
)
nad = MassMetabolite(
    id="nad",
    name="Nicotinamide adenine dinucleotide",
    formula="C21H26N7O14P2",
    charge=-1,
    compartment="c"
)
nadh = MassMetabolite(
    id="nadh",
    name="Nicotinamide adenine dinucleotide - reduced",
    formula="C21H27N7O14P2",
    charge=-2,
    compartment="c"
)
pi = MassMetabolite(
    id="pi",
    name="Phosphate",
    formula="HPO4",
    charge=-2,
    compartment="c"
)
h2o = MassMetabolite(
    id="h2o",
    name="H2O",
    formula="H2O",
    charge=0,
    compartment="c"
)

### Generate Reactions ###

### D-Glucose to D-Glucose 6-phosphate
HEX1 = MassReaction(
    id="HEX1",
    name="Hexokinase (D-glucose:ATP)",
    subsystem="Glycolysis"
)
HEX1.add_metabolites({
    glc__D: -1,
    atp: -1,
    adp: 1,
    g6p: 1,
    h: 1
})

### D-Glucose 6-phosphate to D-Fructose 6-phosphate
PGI = MassReaction(
    id="PGI",
    name="Glucose-6-phosphate isomerase",
    subsystem="Glycolysis"
)
PGI.add_metabolites({
    g6p: -1,
    f6p: 1
})

### D-Fructose 6-phosphate to D-Fructose 1,6-bisphosphate
PFK = MassReaction(
    id="PFK",
    name="Phosphofructokinase",
    subsystem="Glycolysis"
)
PFK.add_metabolites({
    f6p: -1,
    atp: -1,
    fdp: 1,
    adp: 1,
    h: 1
})

### D-Fructose 1,6-bisphosphate to 
### Dihydroxyacetone phosphate + Glyceraldehyde 3-phosphate
FBA = MassReaction(
    id="FBA",
    name="Fructose-bisphosphate aldolase",
    subsystem="Glyolysis"
)
FBA.add_metabolites({
    fdp: -1,
    dhap: 1,
    g3p: 1
})

### Dihydroxyacetone phosphate to Glyceraldehyde 3-phosphate
TPI = MassReaction(
    id="TPI",
    name="Triose-phosphate isomerase",
    subsystem="Glycolysis"
)
TPI.add_metabolites({
    dhap: -1,
    g3p: 1
})

### Glyceraldehyde 3-phosphate to 3-Phospho-D-glyceroyl phosphate
GAPD = MassReaction(
    id="GAPD",
    name="Glyceraldehyde-3-phosphate dehydrogenase",
    subsystem="Glycolysis"
)
GAPD.add_metabolites({
    g3p: -1,
    nad: -1,
    pi: -1,
    _13dpg: 1,
    h: 1,
    nadh: 1
})

### 3-Phospho-D-glycerate to 3-Phospho-D-glyceroyl phosphate
PGK = MassReaction(
    id="PGK",
    name="Phosphoglycerate kinase",
    subsystem="Glycolysis"
)
PGK.add_metabolites({
    _13dpg: 1,
    adp: 1,
    _3pg: -1,
    atp: -1
})

### D-Glycerate 2-phosphate to 3-Phospho-D-glycerate
PGM = MassReaction(
    id="PGM",
    name="Phosphoglycerate mutase",
    subsystem="Glycolysis"
)
PGM.add_metabolites({
    _3pg: 1,
    _2pg: -1
})

### D-Glycerate 2-phosphate to Phophoenolpyruvate
ENO = MassReaction(
    id="ENO",
    name="Enolase",
    subsystem="Glycolysis"
)
ENO.add_metabolites({
    _2pg: -1,
    h2o: 1,
    pep: 1
})

### Phophoenolpyruvate to Pyruvate
PYK = MassReaction(
    id="PYK",
    name="Pyruvate kinase",
    subsystem="Glycolysis"
)
PYK.add_metabolites({
    pep: -1,
    h: -1,
    adp: -1,
    atp: 1,
    pyr: 1
})

### L-lactate to Pyruvate
LDH_L = MassReaction(
    id="LDH_L",
    name="L-lactate dehydrogenase",
    subsystem="Glycolysis"
)
LDH_L.add_metabolites({
    lac__L: 1,
    nad: 1,
    h: -1,
    nadh: -1,
    pyr: -1
})

### AMP + ATP <-> 2 ADP
ADK1 = MassReaction(
    id="ADK1",
    name="Adenylate kinase",
    subsystem="Misc."
)
ADK1.add_metabolites({
    amp: -1,
    atp: -1,
    adp: 2
})

### ATP + H2O <-> ADP + H + Pi
### ### Note: This is a pseudoreaction
ATPM = MassReaction(
    id="ATPM",
    name="ATP maintenance requirement",
    subsystem="Pseudo",
    reversible=False
)
ATPM.add_metabolites({
    atp: -1,
    h2o: -1,
    adp: 1,
    h: 1,
    pi: 1
})

### NADH <-> NAD + H
### ### Note: This is a pseudoreaction
DM_nadh = MassReaction(
    id="DM_nadh",
    name="Demand NADH",
    subsystem="Pseudo",
    reversible=False
)
DM_nadh.add_metabolites({
    nadh: -1,
    nad: 1,
    h: 1
})

# Exchange Reactions below

### (null) <-> D-Glucose
### ### Source for D-Glucose
EX_glc__D = MassReaction(
    id="EX_glc__D",
    name="D-Glucose exchange",
    subsystem="Transport/Exchange",
    reversible=False
)
EX_glc__D.add_metabolites({
    glc__D: 1
})

### (null) <-> AMP
### ### Source for AMP 
### ### Note: BiGG reaction does not exist
EX_AMPIN = MassReaction(
    id="EX_AMPIN",
    name="AMP source",
    subsystem="Transport/Exchange",
    reversible=False
)
EX_AMPIN.add_metabolites({
    amp: 1
})

### AMP <-> (null)
### ### Sink for AMP
### ### Note: BiGG reaction does not exist
EX_AMPOUT = MassReaction(
    id="EX_AMPOUT",
    name="AMP sink",
    subsystem="Transport/Exchange",
    reversible=False
)
EX_AMPOUT.add_metabolites({
    amp: -1
})

### L-Lactate <-> (null)
EX_lac__L = MassReaction(
    id="EX_lac__L",
    name="L-Lactate exchange",
    subsystem="Transport/Exchange"
)
EX_lac__L.add_metabolites({
    lac__L: -1
})

### Pyruvate <-> (null)
### ### Note: This is a pseudoreaction
EX_pyr = MassReaction(
    id="EX_pyr",
    name="Pyruvate exchange",
    subsystem="Transport/Exchange"
)
EX_pyr.add_metabolites({
    pyr: -1
})

### H <-> (null)
### ### Note: This is a pseudoreaction
EX_h = MassReaction(
    id="EX_h",
    name="H+ exchange",
    subsystem="Transport/Exchange"
)
EX_h.add_metabolites({
    h: -1
})

### H2O <-> (null)
### ### Note: BiGG reaction is different
### ### Note: This is a pseudoreaction
EX_h2o = MassReaction(
    id="EX_h2o",
    name="H2O exchange",
    subsystem="Transport/Exchange"
)
EX_h2o.add_metabolites({
    h2o: -1
})

### Generate Model ###

### Generate Model object
glycolysis = MassModel("glycolysis")

### Add reactions to model
### Note: This will add all associated metabolites as well
rxns_list = [HEX1, PGI, PFK, FBA, TPI, GAPD, PGK,
             PGM, ENO, PYK, LDH_L, ADK1, ATPM, DM_nadh,
             EX_glc__D, EX_AMPIN, EX_AMPOUT, EX_lac__L,
             EX_pyr, EX_h, EX_h2o]

glycolysis.add_reactions(rxns_list)

### Populate Model ###

metab_list = [glc__D, g6p, f6p, fdp, g3p, dhap, 
              _13dpg, _3pg, _2pg, pep, pyr, lac__L, 
              atp, adp, amp, h, nad, nadh, pi, h2o]


### Set initial conditions
### All values in units of millimoles/Liter
glc__D.initial_condition = 1.0
g6p.initial_condition = 0.0486
f6p.initial_condition = 0.0198

fdp.initial_condition = 0.0146
g3p.initial_condition = 0.00728
dhap.initial_condition = 0.16

_13dpg.initial_condition = 0.000243
_3pg.initial_condition = 0.0773
_2pg.initial_condition = 0.0113

pep.initial_condition = 0.017
pyr.initial_condition = 0.060301
lac__L.initial_condition = 1.36

atp.initial_condition = 1.6
adp.initial_condition = 0.29
amp.initial_condition = 0.0867281

h.initial_condition = 0.0000899757
nad.initial_condition = 0.0589
nadh.initial_condition = 0.0301

pi.initial_condition = 2.5
h2o.initial_condition = 1.0


### Set gibbs_formations
### All values in units of ""
# Add gibbs energy here #

### Add to model
glycolysis.set_initial_conditions(metab_list)
# Add gibbs energy to model here #

### Set rate/equilibrium constants for reactions in model
HEX1.equilibrium_constant = 850
HEX1.forward_rate_constant = 0.700007 # Liter * Hour^-1 * Millimole^-1
HEX1.reverse_rate_constant = HEX1.kf/HEX1.Keq # Liter * Hour^-1 * Millimole^-1

PGI.equilibrium_constant = 0.41
PGI.forward_rate_constant = 3644.44 # Hour^-1
PGI.reverse_rate_constant = PGI.kf/PGI.Keq # Hour^-1

PFK.equilibrium_constant = 310
PFK.forward_rate_constant = 35.3688 # Liter * Hour^-1 * Millimole^-1
PFK.reverse_rate_constant = PFK.kf/PFK.Keq # Liter * Hour^-1 * Millimole^-1

FBA.equilibrium_constant = 0.082 # Millimole * Liter^-1
FBA.forward_rate_constant = 2834.57 # Hour^-1
FBA.reverse_rate_constant = FBA.kf/FBA.Keq # Liter * Hour^-1 * Millimole^-1

TPI.equilibrium_constant = 0.0571429
TPI.forward_rate_constant = 34.3558 # Hour^-1
TPI.reverse_rate_constant = TPI.kf/TPI.Keq # Hour^-1

GAPD.equilibrium_constant = 0.0179 # Liter * Millimole^-1
GAPD.forward_rate_constant = 3376.75 # Liter^2 * Hour^-1 * Millimole*-2
GAPD.reverse_rate_constant = GAPD.kf/GAPD.Keq # Liter * Hour^-1 * Millimole^-1

### BiGG Reaction is reverse of Mathematica
### kf_BiGG = kr_mathetmatica, kr_BiGG = kf_mathetmatica
### Keq_BiGG = 1/Keq_mathetmatica
PGK.equilibrium_constant = 1/1800
PGK.reverse_rate_constant = 1.27353*10**6 # Liter * Hour^-1 * Millimole^-1
PGK.forward_rate_constant = PGK.Keq*PGK.kr # Liter * Hour^-1 * Millimole^-1

### BiGG Reaction is reverse of Mathematica
### kf_BiGG = kr_mathetmatica, kr_BiGG = kf_mathetmatica
### Keq_BiGG = 1/Keq_mathetmatica
PGM.equilibrium_constant = 1/0.147059
PGM.reverse_rate_constant = 4869.57 # Hour^-1
PGM.forward_rate_constant = PGM.Keq*PGM.kr # Hour^-1

ENO.equilibrium_constant = 1.69492
ENO.forward_rate_constant = 1763.78 # Hour^-1
ENO.reverse_rate_constant = ENO.kf/ENO.Keq # Hour^-1

PYK.equilibrium_constant = 363000
PYK.forward_rate_constant = 454.386 # Liter * Hour^-1 * Millimole^-1
PYK.reverse_rate_constant = PYK.kf/PYK.Keq # Liter * Hour^-1 * Millimole^-1

### BiGG Reaction is reverse of Mathematica
### kf_BiGG = kr_mathetmatica, kr_BiGG = kf_mathetmatica
### Keq_BiGG = 1/Keq_mathetmatica
LDH_L.equilibrium_constant = 1/26300
LDH_L.reverse_rate_constant = 1112.57 # Liter * Hour^-1 * Millimole^-1
LDH_L.forward_rate_constant = LDH_L.Keq*LDH_L.kr # Liter * Hour^-1 * Millimole^-1

### BiGG Reaction is reverse of Mathematica
### kf_BiGG = kr_mathetmatica, kr_BiGG = kf_mathetmatica
### Keq_BiGG = 1/Keq_mathetmatica
ADK1.equilibrium_constant = 1/1.65
ADK1.reverse_rate_constant = 100000 # Liter * Hour^-1 * Millimole^-1
ADK1.forward_rate_constant = ADK1.Keq*ADK1.kr # Liter * Hour^-1 * Millimole^-1

### Irreversible reaction (Keq has units of  Millimolde * Liter^-1)
ATPM.forward_rate_constant = 1.4 # Hour^-1

### Irreversible reactions (Keq is dimensionless)
DM_nadh.forward_rate_constant = 7.44186 # Hour^-1
EX_glc__D.forward_rate_constant = 1.12 # Hour^-1
EX_AMPIN.forward_rate_constant = 0.014 # Hour^-1
EX_AMPOUT.forward_rate_constant = 0.161424 # Hour^-1

EX_lac__L.equilibrium_constant = 1
EX_lac__L.forward_rate_constant = 5.6 # Hour^-1
EX_lac__L.reverse_rate_constant = EX_lac__L.kf # Hour^-1

EX_pyr.equilibrium_constant = 1
EX_pyr.forward_rate_constant = 744.186 # Hour^-1
EX_pyr.reverse_rate_constant = EX_pyr.kf # Hour^-1

EX_h.equilibrium_constant = 1
EX_h.forward_rate_constant = 100000 # Hour^-1
EX_h.reverse_rate_constant = EX_h.kf # Hour^-1

EX_h2o.equilibrium_constant = 1
EX_h2o.forward_rate_constant = 100000 # Hour^-1
EX_h2o.reverse_rate_constant = EX_h2o.kf # Hour^-1

HEX1.ssflux = 1.12 # Millimole * Hour^-1 * Liter^1
PGI.ssflux = 1.12 # Millimole * Hour^-1 * Liter^1
PFK.ssflux = 1.12 # Millimole * Hour^-1 * Liter^1
FBA.ssflux = 1.12 # Millimole * Hour^-1 * Liter^1
TPI.ssflux = 1.12 # Millimole * Hour^-1 * Liter^1

GAPD.ssflux = 2.24 # Millimole * Hour^-1 * Liter^1
PGK.ssflux = 2.24 # Millimole * Hour^-1 * Liter^1
PGM.ssflux = 2.24 # Millimole * Hour^-1 * Liter^1
ENO.ssflux = 2.24 # Millimole * Hour^-1 * Liter^1
PYK.ssflux = 2.24 # Millimole * Hour^-1 * Liter^1

LDH_L.ssflux = 2.016 # Millimole * Hour^-1 * Liter^1

ADK1.ssflux = 0 # Millimole * Hour^-1 * Liter^1
ATPM.ssflux = 2.24 # Millimole * Hour^-1 * Liter^1

DM_nadh.ssflux = 0.224 # Millimole * Hour^-1 * Liter^1
EX_glc__D.ssflux = 1.12 # Millimole * Hour^-1 * Liter^1

EX_AMPIN.ssflux = 0.014 # Millimole * Hour^-1 * Liter^1
EX_AMPOUT.ssflux = 0.014 # Millimole * Hour^-1 * Liter^1

EX_lac__L.ssflux = 2.016 # Millimole * Hour^-1 * Liter^1
EX_pyr.ssflux = 0.224 # Millimole * Hour^-1 * Liter^1

EX_h.ssflux = 2.688 # Millimole * Hour^-1 * Liter^1
EX_h2o.ssflux = 0 # Millimole * Hour^-1 * Liter^1