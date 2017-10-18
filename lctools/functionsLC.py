#Some functions for dealing with dates and fluxes

from numpy import *
from scipy.stats import chi2
from scipy.stats import norm as statsnorm
from scipy.optimize import curve_fit
from mpmath import gammainc
import datetime as dt
from astropy.time import *
#from lctools.functionsEBL import *

#powerlaw with cutoff
def plco(x,n,si,co,p):
    return n*pow(x/p,-si)*exp(-x/co)

def cuflux(n,s,c=1000000.,e_th=1.0,norm_crab=3.45e-11,si_crab=2.63):
    f = -1.
    if ((c==0) | (c>=100000)):
        f = n/norm_crab*(1.-si_crab)/(1.-s)*e_th**(si_crab-s)
    else:
        f = - n/norm_crab*c**(1.-s)*(1.-si_crab)*e_th**(-1.+si_crab)*gammainc(1.-s,e_th/c)
    return f

def integratedflux(n,s,c=1000000.,e_th=1.0, z=0.):
    if (z>0):
        f = integratedFluxEBL(n ,s ,c ,e_th, z)
    else:
        if ((c==0) | (c>=100000)):
            f = - n/(1.-s)*e_th**(1.-s)
        else:
            f = n*c**(1.-s)*gammainc(1.-s,float(e_th)/c)
    return f

def sdname2date(sdname):
    try:
        d = (sdname.strip()).split("T")
        f = d[0].split("-")    
        try:
            t = d[1].split(":")    
            date = dt.datetime(int(f[0]),int(f[1]),int(f[2]),int(t[0]),int(t[1]))
        except:
            date = dt.datetime(int(f[0]),int(f[1]),int(f[2]))
    except:
        raise
    return date


def date2sdname(date,ra):
    try:
        lst = UTC2LST(date)
        lstmidnight = date - dt.timedelta(days=lst/solar2sidereal)
    except:
        raise
    if ((ra>45) & (ra<=315)):
        z = LST2UTC(0,dt.datetime(lstmidnight.year,lstmidnight.month,lstmidnight.day))
        return "a_%d-%02d-%02d_%02dh%02d"%(z.year,z.month,z.day,z.hour,z.minute)
    else:
        z = LST2UTC(12,dt.datetime(lstmidnight.year,lstmidnight.month,lstmidnight.day))
        return "b_%d-%02d-%02d_%02dh%02d"%(z.year,z.month,z.day,z.hour,z.minute)

        

def dt2mjd(date):
    #HAWC:
    #utc = UTCDateTime(date.year,date.month,date.day,date.hour,date.minute,date.second)
    #mjd = ModifiedJulianDate(utc)
    #astropy:
    utc = Time(date)
    mjd = utc.mjd
    return mjd

def mjd2dt(mjd):
    #HAWC:
    #m = ModifiedJulianDate(mjd*U.day)
    #d = dt.datetime(m.datetime.year,m.datetime.month,m.datetime.day,m.datetime.hour,m.datetime.minute,m.datetime.second)
    m = Time(mjd, format='mjd')
    d = m.datetime
    return d


##################################################################
# The following functions allow conversion between UTC and
# sidereal time at the HAWC location

longHAWC = -(97.+(18.+27./60.)/60.)/360.*24.
solar2sidereal = 24.06570982441908

def DaysSince2000Jan1(d):
	d2000 = dt.datetime(2000,1,1,12)
	t = d-d2000
	r = t.days+t.seconds/3600./24.
	return r

def LST2UTC(hour,d):
	#hour is floating point LST hour
	#d (datetime obj) is UTC date (year,month,day) 
	gmst = (hour-longHAWC)%24
	offset = (18.697374558 + solar2sidereal*DaysSince2000Jan1(d))%24
	days = (gmst-offset)/solar2sidereal
	if (days<0):
		rd = d+dt.timedelta(days+24/solar2sidereal)
	else:
		rd = d+dt.timedelta(days)
	return rd

def UTC2LST(d):
	#from Wikipedia, good to 1 sec per century and
	#assuming UTC = UT1 (good to <1 sec accuracy)
	D =  DaysSince2000Jan1(d)
	GMST = (18.697374558 + solar2sidereal*D)
	LST = (GMST+longHAWC)%24.
	return LST
