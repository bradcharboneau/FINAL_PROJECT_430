VACUUM = 0
WATER = 1
FUEL = 2
#fuel properties, im just using pre determined 2 group values right now for the solver verification, ill add the enrichment function later
# used values from Nucl. Sci. Eng. 38, 8–32, 1969, Fig. 3
seed_props = {
    "D": [1.4, 0.4],
    "sig_a": [0.01, 0.15],
    "sig_s12": 0.01,
    "nu_sig_f": [0.00765703, 0.218772]
}
materials = {FUEL: seed_props}


