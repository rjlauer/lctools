import sys
from numpy import *
from functionsLC import *

class Flux(object):

    def __init__(self):
        self._diffflux = 1.e-11
        return

    @property
    def diffflux(self):
        return self._diffflux

    @property
    def intflux(self):
        return integratedflux(self._diffflux,2.63,1000000,1.0)

class HAWC:
    
    def __init__(self,lcfile, 
                      mjd_begin = -1, 
                      mjd_end   = -1, 
                      mintransits = 0.5, 
                      crabflux  = 1.889e-11, 
                      intthresh = 1.,
                      redshift  = 0.,
                      nooverlap = True, 
                      logging = False):
         
        self.logging = logging
        
        """
        self._fitLL = array([], dtype=float)
        self._bgLL  = array([], dtype=float)
        self._fixLL = array([], dtype=float)
        self._flux  = array([], dtype=float)
        self._flux_err = array([], dtype=float)
        self._mjd   = array([], dtype=float) 
        self._mjd_tlerr = array([], dtype=float)
        self._mjd_tuerr  = array([], dtype=float)
        self._transits = array([], dtype=float)
        self._conversion = 0.
        self._fixNormIntegrated = 0.
        self._sdperiod = 0.
        """
        
        self.loadLC(lcfile, mjd_begin, mjd_end, mintransits, 
                    crabflux, intthresh, redshift, nooverlap)
    
    @property
    def fitLL(self):
        return self._fitLL
    
    @property
    def bgLL(self):
        return self._bgLL
    
    @property
    def fixLL(self):
        return self._fixLL
    
    @property
    def bgTS(self):
        return self._bgTS        
    
    @property
    def fixTS(self):
        return self._fixTS
    
    @property
    def sumBgTS(self):
        return self.bgTS.sum()
    
    @property
    def sumFixTS(self):
        return self.fixTS.sum()
    
    @property
    def flux(self):
        return self._flux
    
    @property
    def flux_err(self):
        return self._flux_err
    
    @property
    def norm1TeV(self):
        return self._norm1TeV
    
    @property
    def norm1TeV_err(self):
        return self._norm1TeV_err
    
    @property
    def mjd(self):
        return self._mjd
    
    @property
    def mjd_tlerr(self):
        return self._mjd_tlerr
    
    @property
    def mjd_tuerr(self):
        return self._mjd_tuerr
    
    @property
    def transits(self):
        return self._transits
    
    @property
    def index(self):
        return self._index
    
    @property
    def cutoff(self):
        return self._cutoff
    
    @property
    def conversion(self):
        return self._conversion
    
    @property
    def fixNormIntegrated(self):
        return self._fixNormIntegrated
    
    @property
    def sdperiod(self):
        return self._sdperiod
    
    @property
    def redshift(self):
        return self.redshift
    
    
    def loadLC(self, lcfile, mjd_begin, mjd_end, mintransits, 
               crabflux, intthresh, redshift, nooverlap=True):
        
        self.lcfile = lcfile
        self.mjd_begin = mjd_begin 
        self.mjd_end   = mjd_end
        self.mintransits = mintransits 
        self.crabflux  = crabflux 
        self.intthresh = intthresh
        self.redshift = redshift
        
        file = open(lcfile,"r")
        
        firstpeak = []
        lastpeak = []
        startmjd = []
        stopmjd = []
        sdperiod = []
        transits = []
        slope = []
        cutoff = []
        norm = []
        norm_err = []
        fitLL = []
        bgLL = []
        fixNorm = []
        fixLL = []
        lerr = []
        uerr = []
        
        one_slope = 0.
        same_slope = True
        one_cutoff = 0.
        same_cutoff = True
        n = 0
        for row in file:
            try:
                if (row[0]=="#"):
                    continue
                field = row.split()
                #ignore overlapping periods
                if (n>0):
                    if ((nooverlap) & (float(field[2])<=stopmjd[n-1])):
                        continue
                firstpeak.append(field[0])
                lastpeak.append(field[1])
                startmjd.append(float(field[2]))
                stopmjd.append(float(field[3]))
                sdperiod.append(int(field[4]))
                transits.append(float(field[5]))
                if (n==0):
                    one_slope = float(field[6])
                elif (float(field[6])!=one_slope):
                    #print "  [loadLC: %s] WARNING: different slopes in light curve, do not plot 2 axes!"%(lcfile)
                    same_slope = False
                    #sys.exit()
                slope.append(float(field[6]))
                if (n==0):
                    one_cutoff = float(field[7])
                elif (float(field[7])!=one_cutoff):
                    #print "  [loadLC: %s] ERROR: different slopes in light curve, do not plot 2 axes."%(lcfile)
                    same_cutoff = False
                    #sys.exit()
                cutoff.append(float(field[7]))
                norm.append(float(field[8]))
                norm_err.append(float(field[9]))
                fitLL.append(float(field[10]))
                bgLL.append(float(field[11]))
                fixNorm.append(float(field[12]))
                fixLL.append(float(field[13]))
                if (len(field)==16):
                    lerr.append(float(field[14]))
                    uerr.append(float(field[15]))
                    #norm_err[n] = max([lerr[n],uerr[n]]) 
                    if (norm[n]<0): norm_err[n] = uerr[n]
                n += 1
            except:
                raise
                sys.exit()
        
        file.close()
        
        if self.logging:
            print "  [loadLC: %s] Found %d light curve data points in file." %(lcfile,len(norm))
        
        if  ((same_slope==False) | (same_cutoff==False)):  
            print "  [loadLC: %s] WARNING: different slopes or cut-offs in light curve, do not plot 2 axes!"%(lcfile)
        #this is only a universal conversion factor if same_slope and same_cutoff are True
        intflux = integratedflux(1.,one_slope,one_cutoff,self.intthresh,self.redshift)
       
        if self.logging:
            print "  [loadLC: %s] Crab flux, integrated above 1 TeV: %.3e ph cm^-2 s^-1 "%(lcfile,self.crabflux)
        #calculate crab flux above intthreshold:
        self._conversion = float(self.crabflux/integratedflux(1.,2.63,1000000,1.)*integratedflux(1.,2.63,1000000.,self.intthresh))
        
        # from file:
        firstmjd = startmjd[0]
        lastmjd = stopmjd[-1]
        
        #if not given as input: start/stop at first/last date stored
        if (mjd_begin == -1): self.mjd_begin = firstmjd
        if (mjd_end == -1):   self.mjd_end   = lastmjd
        
        #reset first/last to the min. range of input dates and available data
        if (firstmjd < self.mjd_begin):
            firstmjd = self.mjd_begin
        if (lastmjd > self.mjd_end):
            lastmjd = self.mjd_end
        
        self.mjd_begin = floor(firstmjd)
        self.mjd_end = ceil(lastmjd)
        if self.logging:
            print "  [loadLC: %s] extracting MJDs from %.2f to %.2f"%(
                    lcfile,self.mjd_begin, self.mjd_end)
        
        flux = zeros(len(norm))
        flux_err = zeros(len(norm))
        TS = zeros(len(norm))
        sigma = zeros(len(norm))
        
        good = []
        
        minsigma = 1000.
        maxsigma = -1000.
        colors = []
        ignoredTransits = 0
        totaltransits = 0.
        for d in xrange(0,len(norm)):
            currentmjd = startmjd[d]
            if (currentmjd<self.mjd_begin):
                continue
            if (currentmjd>self.mjd_end):
                break
            if (transits[d]<self.mintransits):
                ignoredTransits += 1
                continue
            if (norm[d]!=norm[d]):
                print "  [loadLC: %s] Skipping MJDs %f - %f due to nan: norm = %.2f (fixTS = %.2f)" %(lcfile,startmjd[d],stopmjd[d],norm[d],2.*(fitLL[d]-fixLL[d]))
                continue
            if (fitLL[d]>=0.):
                print "  [loadLC: %s] Skipping MJDs %f - %f due to pos. fitLL value, i.e. failed fit." %(lcfile,startmjd[d],stopmjd[d])
                continue
            if ((fitLL[d]-bgLL[d]<-0.01) | (fitLL[d]-fixLL[d]<-0.01)):
                print "  [loadLC: %s] Skipping MJDs %f - %f due to fitLL (%.1f) more than 0.01 below bgLL (%.1f) or fixLL (%.1f)." %(lcfile,startmjd[d],stopmjd[d],fitLL[d],bgLL[d],fixLL[d])
                continue
            #experimental:
            #if (norm_err[d]<(0.01*fabs(norm[d]))):
            #    print "  [loadLC: %s] Skipping MJDs %f - %f due to norm-error (%.2e) less than 1/100 of norm (%.2e)." %(lcfile,startmjd[d],stopmjd[d],norm_err[d],norm[d])
            #    continue
                  
            #integrate above intthresh
            if  ((same_slope==False) | (same_cutoff==False)):  
                #individual integration, otherwise already calculated above
                intflux = integratedflux(1.,slope[d],cutoff[d],self.intthresh,self.redshift)
            
            flux[d] = norm[d]*intflux
            flux_err[d] = norm_err[d]*intflux
            if (len(lerr)>0):
                lerr[d]*=intflux
                uerr[d]*=intflux
            
            #TS values:
            bgTS = 2*(fitLL[d]-bgLL[d])
            sigma[d] = sqrt(bgTS)
            
            mjd = (startmjd[d]+stopmjd[d])/2.
            if ((mjd>=self.mjd_begin) & (mjd<=self.mjd_end)): 
                good.append(d)
                totaltransits += transits[d]
                if (sigma[d]<minsigma): minsigma = sigma[d]
                if (sigma[d]>maxsigma): maxsigma = sigma[d]
                
        if self.logging:
            print "  [loadLC: %s] Ignored %d transits with signal coverage fraction below %.2f"%(lcfile,ignoredTransits,self.mintransits)
            print "  [loadLC: %s] Total transit signal fraction: %f"%(lcfile,totaltransits)
            print "  [loadLC: %s] Significances: min = %5.2f , max = %5.2f"%(lcfile,minsigma,maxsigma)
        
        self._fitLL = array(fitLL, dtype=float)[good]
        self._bgLL = array(bgLL, dtype=float)[good]
        self._fixLL = array(fixLL, dtype=float)[good]
        self._norm1TeV = array(norm, dtype=float)[good]
        self._norm1TeV_err = array(norm_err, dtype=float)[good]
        self._flux = array(flux, dtype=float)[good]
        self._flux_err = array(flux_err, dtype=float)[good]
        self._mjd = (array(startmjd, dtype=float) + array(stopmjd, dtype=float))[good]/2.
        self._mjd_tlerr = self.mjd - array(startmjd, dtype=float)[good]
        self._mjd_tuerr = array(stopmjd, dtype=float)[good] - self.mjd
        self._transits = array(transits, dtype=float)[good]
        self._index = array(slope, dtype=float)[good]
        self._cutoff = array(cutoff, dtype=float)[good]
       
        self._bgTS = 2*(self._fitLL-self._bgLL)
        for d,t in enumerate(self._bgTS):
            if (t<0):
                print "  [loadLC: %s] bgTS for MJD %.1f is negative (%.1f), setting to zero."%(lcfile,self.mjd[d],t)
                self._bgTS[d] = 0.

        self._fixTS = 2*(self._fitLL-self._fixLL)
        for d,t in enumerate(self._fixTS):
            if (t<0):
                print "  [loadLC: %s] fixTS for MJD %.1f is negative (%.1f), setting to zero."%(lcfile,self.mjd[d],t)
                self._fixTS[d] = 0.
        self._fixNormIntegrated = fixNorm[0]*integratedflux(1.,one_slope,one_cutoff,self.intthresh,self.redshift)

        if ((same_slope==False) | (same_cutoff==False)):  
            print "  [loadLC: %s] Different slopes or cut-offs used in light curve, don't use fixTS, fixLL, or fixNormIntegrated."%(lcfile)
            self._fixLL = zeros(len(self._fixLL ))
            self._fixTS = zeros(len(self._fixTS ))
            self._fixNormIntegrated = 0.

        self._sdperiod = sdperiod[-1]
        
        return


class HAWCint:
    
    def __init__(self,lcfile, 
                      mjd_begin = -1, 
                      mjd_end   = -1, 
                      logging = False):
        
        self.logging = logging
        self._intThreshold = 1.
        self._minTransits = 0.
        #default, >1 TeV:
        self._crabflux = 1.889e-11
        
        self.loadLC(lcfile, mjd_begin, mjd_end)
    
    @property
    def mjd(self):
        return self._mjd
    
    @property
    def mjd_tlerr(self):
        return self._mjd_tlerr
    
    @property
    def mjd_tuerr(self):
        return self._mjd_tuerr
    
    @property
    def flux(self):
        return self._flux
    
    @property
    def flux_err(self):
        return self._flux_err

    @property
    def fluxCU(self):
        return self._fluxCU
    
    @property
    def fluxCU_err(self):
        return self._fluxCU_err

    @property
    def transits(self):
        return self._transits
    
    @property
    def TS(self):
        return self._TS
    
    @property
    def intThreshold(self):
        return self._intThreshold
    
    @property
    def minTransits(self):
        return self._minTransits
    
    @property
    def crabflux(self):
        return self._crabflux
    
    def loadLC(self,lcfile,mjd_begin,mjd_end):
        
        #if input dates given, use those as limits for scanning file
        if (mjd_begin>=0):
            self.mjd_begin = mjd_begin 
            if self.logging:
                print "  [loadLC: %s] extracting MJDs from %.2f"%(
                        lcfile,self.mjd_begin)
        if (mjd_end>0):
            self.mjd_end   = mjd_end 
            if self.logging:
                print "  [loadLC: %s] extracting MJDs to %.2f"%(
                        lcfile, self.mjd_end)
        else:
            self.mjd_end   = 100000. 
        
        Mstart = []
        Mstop = []
        flux = []
        flux_err = []
        fluxCU = []
        fluxCU_err = []
        transits = []
        TS = []
        
        ffile = open(lcfile,"r")
        for row in ffile:
            try:
                if (row[0]=="#"):
                    cf = row.split()
                    if (len(cf)==4):
                        if ((cf[0] =="#integration") & (cf[1]=="threshold:")):
                            self._intThreshold = float(cf[2])
                    if (len(cf)==2):
                        if (cf[0] =="#minimum-transit-coverage:"):
                            self._minTransits = float(cf[1])
                    if (len(cf)==6):
                        if ((cf[0] =="#Crab") & (cf[1]=="flux")):
                            self._crabflux = float(cf[5])
                    continue
                elif (row==""):
                    continue
                field = row.split()
                if (len(field)<4):
                    continue
                if ((float(field[1])<self.mjd_begin) | 
                    (float(field[0])>self.mjd_end)):
                    continue
                Mstart.append(float(field[0]))
                Mstop.append(float(field[1]))
                flux.append(float(field[2]))
                flux_err.append(float(field[3]))
                fluxCU.append(float(field[4]))
                fluxCU_err.append(float(field[5]))
                transits.append(float(field[6]))
                TS.append(float(field[9]))
            except:
                raise
                sys.exit()
        ffile.close()
        if self.logging:
            print "  [loadLC: %s] Loaded %d data points from light curve"%(lcfile,len(fflux))
        
        self._flux = array(flux, dtype=float)
        self._flux_err = array(flux_err, dtype=float)
        self._fluxCU = array(fluxCU, dtype=float)
        self._fluxCU_err = array(fluxCU_err, dtype=float)
        self._transits = array(transits, dtype=float)
        self._TS = array(TS, dtype=float)
         
        # from file:
        firstmjd = Mstart[0]
        lastmjd = Mstop[-1]
        
        #if not given as input: start/stop at first/last date stored
        if (mjd_begin == -1): self.mjd_begin = firstmjd
        if (mjd_end == -1):   self.mjd_end   = lastmjd
        
        #reset first/last to the min. range of input dates and available data
        if (firstmjd < self.mjd_begin):
            firstmjd = self.mjd_begin
        if (lastmjd > self.mjd_end):
            lastmjd = self.mjd_end
        
        self.mjd_begin = floor(firstmjd)
        self.mjd_end = ceil(lastmjd)
        if self.logging:
            print "  [loadLC: %s] extracted MJDs from %.2f to %.2f"%(lcfile,self.mjd_begin, self.mjd_end)
        
        
        fmjd = []
        fmjderr = []
        
        for f in xrange(0,len(Mstart)):
            fmjderr.append((Mstop[f]-Mstart[f])/2.)
            fmjd.append(Mstart[f]+fmjderr[f])
        
        self._mjd       = array(fmjd, dtype=float)
        self._mjd_tlerr = array(fmjderr, dtype=float)
        self._mjd_tuerr = array(fmjderr, dtype=float)
         
        return 
    
class FACT:
    
    def __init__(self,lcfile, 
                      mjd_begin = -1, 
                      mjd_end   = -1, 
                      logging = False):
        
        self.logging = logging
        
        self.loadLC(lcfile, mjd_begin, mjd_end)
    
    @property
    def mjd(self):
        return self._mjd
    
    @property
    def mjd_tlerr(self):
        return self._mjd_tlerr
    
    @property
    def mjd_tuerr(self):
        return self._mjd_tuerr
    
    @property
    def excess(self):
        return self._excess
    
    @property
    def excess_err(self):
        return self._excess_err

    @property
    def flux(self):
        return self._flux
    
    @property
    def flux_err(self):
        return self._flux_err

    @property
    def significance(self):
        return self._sig

    def loadLC(self,lcfile,mjd_begin,mjd_end):
       
        if (mjd_begin>0):
            self.mjd_begin = mjd_begin 
        else:
            self.mjd_begin = 0.
        if (mjd_end>0):
            self.mjd_end   = mjd_end 
        else:
            self.mjd_end   = 100000. 
        
        Mtime = []
        Mstart = []
        Mstop = []
        fexcess = []
        fexcess_err = []
        fflux = []
        fflux_err = []
        fsig = []
        
        ffile = open(lcfile,"r")
        for row in ffile:
            try:
                if (row[0]=="#"):
                   continue
                elif (row==""):
                    continue
                field = row.split()
                if (len(field)<4):
                    continue
                if ((float(field[1])<self.mjd_begin) | (float(field[0])>self.mjd_end)):
                    continue
                Mtime.append(float(field[0]))
                Mstart.append(float(field[1]))
                Mstop.append(float(field[2]))
                fexcess.append(float(field[3]))
                fexcess_err.append(float(field[7]))
                fflux.append(float(field[4])*1.e-11)
                fflux_err.append(float(field[8])*1.e-11)
                fsig.append(float(field[9]))
            except:
                raise
                sys.exit()
        ffile.close()
        if self.logging:
            print "  [loadLC: %s] Loaded %d data points from light curve"%(lcfile,len(fflux))
        
        self._excess = array(fexcess, dtype=float)
        self._excess_err = array(fexcess_err, dtype=float)
        self._flux = array(fflux, dtype=float)
        self._flux_err = array(fflux_err, dtype=float)
        
        fmjd = []
        fmjd_tlerr = []
        fmjd_tuerr = []
        
        for f in xrange(0,len(Mstart)):
            fmjd.append(Mtime[f])
            fmjd_tlerr.append(Mtime[f]-Mstart[f])
            fmjd_tuerr.append(Mstop[f]-Mtime[f])
        
        self._mjd       = array(fmjd, dtype=float)
        self._mjd_tlerr = array(fmjd_tlerr, dtype=float)
        self._mjd_tuerr = array(fmjd_tuerr, dtype=float)
         
        self._sig       = array(fsig, dtype=float)
        
        return 

