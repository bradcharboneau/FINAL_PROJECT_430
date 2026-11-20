VACUUM = 0
WATER = 1
FUEL = 2
#fuel properties, im just using pre determined 2 group values right now for the solver verification, ill add the enrichment function later
# used values from L. A. Hageman & J. B. Yasinsky Nucl. Sci. Eng. 38, 8–32, 1969, Fig. 3, apparently this is the common benchmark example for 2 group diffusuion, found it on ans website
seed_props = {
    "D": [1.4, 0.4],
    "sig_a": [0.01, 0.15],
    "sig_s12": 0.01,
    "nu_sig_f": [0.00765703, 0.218772]
}
materials = {FUEL: seed_props}


