import numpy as np
"""
Library of editable parameters for each event type that a user may wish to assess
"""

# dimensions of data in real space
dim = 6

# ordered list of event types to consider from the various population options listed below
events = ["explosions", "earthquakes"]

# priors for each population
prior = {"explosions": 0.9999, "earthquakes": 0.0001}

# gamma values from spherical location model
gamma = {
    "explosions": np.array([
        [ 0.36740822,  0.28912206,  0.07658822,  0.36036246,  0.28245966,  0.75226945],
        [ 0.41243448,  0.37831385,  0.09172119,  0.57011912, -0.13295873, -0.57935293],
        [ 0.83208951, -0.28548006, -0.11076200, -0.45832301, -0.02211478, -0.05754047],
        [-0.00453744, -0.43562773, -0.77987658,  0.44715794,  0.04086863,  0.01949217],
        [ 0.03888890, -0.70802624,  0.60379566,  0.36350937, -0.00777659,  0.0204391 ],
        [-0.03167187, -0.02674923,  0.02601096, -0.05434901,  0.94885127, -0.30713557]
        ]),
    "earthquakes": np.array([
        [-0.42034001,  0.13441548,  0.80750027,  0.29492754,  0.25714967,  0.00904755],
        [ 0.88547429,  0.22260801,  0.30730861,  0.18300082,  0.16013458, -0.11318128],
        [-0.04782012, -0.01290958, -0.36344511,  0.92298289,  0.00720227,  0.11621081],
        [ 0.10057037,  0.04350852,  0.02754056, -0.10749762,  0.14407665,  0.97720024],
        [-0.06296248, -0.00478161, -0.29012272, -0.10784493,  0.93907639, -0.13544993],
        [-0.15127501,  0.96451835, -0.19100448, -0.06666910, -0.07454231, -0.01833563]
        ]),
}

# optimal PDF parametes for the SvMF distribution
kappa = {"explosions": 54.65367804, "earthquakes": 2.60879588}
a = {
    "explosions": np.array([2.0725349917445777, 2.0527225979792427, 1.2762063394892806, 0.47857163522325175, 0.3848575516754583]),
    "earthquakes": np.array([1.5258071357847298, 1.0744906950919515, 0.9805314488231607, 0.8035416739833638, 0.7741547209131909]),
}

# optimal PDF parametes for the vMF distribution (turn on to use a vMF instead of SvMF)
#kappa = {"explosions": 34.49525335, "earthquakes": 2.54531175}
#a = {"explosions": np.array([1, 1, 1, 1, 1]), "earthquakes": np.array([1, 1, 1, 1, 1]),}
