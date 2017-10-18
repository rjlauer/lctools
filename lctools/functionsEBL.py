from hawc import hawcnest, data_structures
from hawc.hawcnest import HAWCUnits as U
from hawc.data_structures import *
from hawc import grmodel_services as grmodel_services
from HAWCNest import HAWCNest

from numpy import *

#powerlaw with cutoff
def plco(x,n,si,co,p):
    return n*pow(x/p,-si)*exp(-x/co)

class EBL:
    
    def __init__(self, redshift, model="Gilmore12FiducialEBLModel"):
        self.setModel(model)
        self.redshift = redshift
    def __del__(self):
        pass
        #print "destructed EBL"

    def setModel(self,model="Gilmore12FiducialEBLModel"):
        self.nest = HAWCNest()
        hawcnest.SetLoggingLevel(5, False)
        self.nest.Service(model, model)
        self.nest.Configure()
        self.eblModel = grmodel_services.GetEBLAbsorptionService(model)
        
    def attenuation(self,energy):
        return self.eblModel.GetAttenuation(energy * U.TeV, self.redshift)
    
    def tau(self,energy):
        return self.eblModel.GetOpticalDepth(energy * U.TeV, self.redshift)

def integratedFluxEBL(n,s,c=1000000.,e_th=1.,z=0., pivot=1.):
    upperEn = 200.
    nbins = 10000
    eTeV = linspace(e_th,upperEn,nbins)
    ebl = EBL(z)
    #integral:
    if isinstance(n, ndarray):
        intflux = zeros(len(n))
    else:
        intflux = 0.
    for i,e in enumerate(eTeV):
        eblatten = ebl.attenuation(e)
        fluxes  = plco(e,n,s,c,pivot)
        fluxes *= eblatten
        intflux += fluxes*(upperEn-e_th)/nbins
    ebl = 0
    return intflux





