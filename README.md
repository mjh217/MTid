# MTid
`MTid` is a Python3 routine for assessing the likelihood that a seismic moment tensor belongs to any number of predefined event populations. For example, it can determine the probability that an event was an explosion as opposed to an earthquake. It does so by reading in the six independent components of the seismic moment tensor, normalising them, projecting them onto a hypersphere, and evaluating any number of probability distribution functions that have previously been defined on that hypersphere. It follows the methodology outlined in [*Hoggard et al. (2024)*](https://doi.org/10.1093/gji/ggae011), which builds on work by [*Ford et al. (2020)*](https://doi.org/10.1093/gji/ggz578). 

## Simple usage example
The file `example_MTs.csv` contains moment tensor solutions for nine seismic events that were published in [*Alvizuri & Tape (2018)*](https://doi.org/10.1785/0220180158). The first six are for nuclear tests carried out by the Democratic People's Republic of Korea from 2006 to 2017, the seventh is believed to be a *"cavity collapse"* that immediately followed the 2017 test, and the final two are earthquakes that ocurred in South Korea in 2016 and 2017.

The six components of the moment tensor must be comma-separated, can be in any unit (provided that it is consistent for all six components), and must be in a specific order: *M<sub>xx</sub>* is in the first column, followed by *M<sub>xy</sub>*, *M<sub>xz</sub>*, *M<sub>yy</sub>*, *M<sub>yz</sub>* and *M<sub>zz</sub>* (*x* represents east, *y* is north, and *z* is the vertical/radial component). The program can be run on the command line using

`python MTid.py example_MTs.csv`

which will state that only the first six events belong to the explosion population and will also output two files. The first is `prob.csv` and contains a comma-separated list of the probabilities that each moment tensor in the input file belongs to each of the possible event populations. The second is `lune_coords.csv` and contains a comma-separated list of the coordinates (longitude, latitude) of the moment tensors when projected onto the fundamental lune of [*Tape & Tape (2012)*](https://doi.org/10.1111/j.1365-246X.2012.05491.x), which can be useful for plotting purposes. A help message can be accessed using

`python MTid.py -h`

## Aspects that a user may wish to alter
The most likely modification that a user may wish to make is to alter the prior probability that an event was an explosion, collapse or earthquake from the default values (see Section 3.2 of [*Hoggard et al. 2024*](https://doi.org/10.1093/gji/ggae011)). This modification can be acheived by locating the library `populations.py` and updating the `prior =` variables, staying faithful to the fact that the sum of all priors must be equal to 1. It is also within this same library that a user may wish to update parameters that control the shape of the probability distribution functions or add new population descriptions for other types of seismic event.

## References
- *Hoggard, M.J., J.L. Scealy & B.G. Delbridge (2024). Seismic moment tensor classification using elliptical distribution functions on the hypersphere. Geophysical Journal International, 237, 1–13.*
- *Ford, S.R., G.D. Kraft & G.A. Ichinose (2020). Seismic moment tensor event screening. Geophysical Journal International, 221, 77–88.*
- *Alvizuri, C. & C. Tape (2018). Full moment tensor analysis of nuclear explosions in North Korea. Seismological Research Letters, 89, 2139–2151.*
- *Tape, W. & C. Tape (2012). A geometric setting for moment tensors. Geophysical Journal International, 190, 476–498.*
