snapper
=======

PyMC implementation of Smith et al. (2012) linked occupancy-abundance model

Created by M. Aaron MacNeil on 20/07/12.

"Smith, Adam Nicholas Howard, Marti Jane Anderson, and Russell Brian Millar. In press. Incorporating the intraspecific occupancy-abundance relationship into zero-inflated models. Ecology. http://dx.doi.org/10.1890/12-0460.1
Incorporating the intraspecific occupancy-abundance relationship into zero-inflated models

Adam Nicholas Howard Smith1,*, Marti Jane Anderson2, and Russell Brian Millar3
1Massey University, Albany Campus, New Zealand Institute for Advanced Study

2Massey University, New Zealand Institute for Advanced Study

3University of Auckland, Department of Statistics

Zero-inflated versions of standard distributions for count data are often required in order to account for excess zeros when modeling the abundance of organisms. Such distributions typically have as parameters λ, the mean of the count distribution, and π, the probability of an excess zero. Implementations of zero-inflated models in ecology typically model λ using a set of predictor variables, and π is fit either as a constant or with its own separate model. Neither of these approaches makes use of any relationship that might exist between π and λ. However, for many species, the rate of occupancy is closely and positively related to its average abundance. Here, this relationship was incorporated into the model by functionally linking π to λ, and demonstrated in a study of snapper (Pagrus auratus) in and around a marine reserve. This approach has several potential practical advantages, including better computational performance and more straightforward model interpretation. It is concluded that, where appropriate, directly linking π to λ can produce more ecologically accurate and parsimonious statistical models of species abundance data."