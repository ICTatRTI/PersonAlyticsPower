# this file will contain the the *Data functions
#
#  polyData
#  hyperData
#  invPolyData
#  expData
#  negExpData
#  logisticData
#
#  that will be used by
#
#  polyICT$makeData()
#  hyperICT$makeData()
#  invPolyICT$makeData()
#  expICT$makeData()
#  negExpICT$makeData()
#  logisticICT$makeData()
#
# These functions do not simulate data, rather they construct
# latent growth curve model data for one group within one phase
# from inputs of random effects and errrors.
#
# Mutligroup and/or multiphase data are achieved by concatenating data
# from multiple calls to *Data functions within *makeData methods.

# level 2 covariates -
# -- categorical: acheived via multiple groups
# -- continuous: implemented directly in lgmSim

# time varying covariates -
# -- continuous: implemented directly in lgmSim
# -- categorical: achieve by categarizing continuous covariates

polyData <- function()
{

}