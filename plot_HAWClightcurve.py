#from pylab import *
import matplotlib.pyplot as plt
import sys
from numpy import *
from scipy.stats import chi2
from scipy.stats import norm as statsnorm
from scipy.optimize import curve_fit
from scipy.special import erf
from scipy.stats import kstest
from mpmath import gammainc
import datetime as dt
import argparse
from hawc import hawcnest, data_structures
from hawc.data_structures import *
from hawc.hawcnest import HAWCUnits as U
from lctools import *

p = argparse.ArgumentParser(description="Plot light curve")
p.add_argument("lcfile", help="light curve ascii file")
p.add_argument("-m", "--minflux", dest="minflux", type=float, default = -1.5, 
                 help="minimum flux value [Crab Units] (default = -2)")
p.add_argument("-M", "--maxflux", dest="maxflux", type=float, default = 6.5, 
                 help="maximum flux value [Crab Units] (default = 6.5)")
p.add_argument("-i", "--phThreshold", dest="phThreshold", type=float, 
                 default = -1., 
                 help="integration threshold for photon flux values "
                 +"(default: same as in input file)")
p.add_argument("-c", "--cuThreshold", dest="cuThreshold", type=float, 
                 default = -1., 
                 help="integration threshold for CU flux values "
                 +"(default: same as in input file)")
p.add_argument("--zero-off", dest="zero", action="store_false", default=True,
                 help="Show line at zero flux.")
p.add_argument("-s","--nocrab", dest="showcrab", action="store_false", 
                 default=True,
                 help="Turn off right axis with Crab units and grid lines "
                 +"(default: show)")
p.add_argument("--no-errorbars", dest="errorbars", action="store_false", 
                 default=True,
                 help="Do not show errorbars on data points (default: show)")
p.add_argument("-e", "--events", dest="events", type=str, default= "", 
                 help="file with event MJD and label")
p.add_argument("--showavg", dest="showavg", action="store_true", default=False,
                 help="Show line at average flux (default: false)")
p.add_argument("--showfixednorm",dest="showfixednorm", action="store_true", 
                 default=False,
                 help="Show const flux best fit")
p.add_argument("--showCharFlux",dest="showCharFlux", action="store_true", 
                 default=False,
                 help="Show const flux best fit (default: false)")
p.add_argument("-o", "--output", dest="output", type=str, 
                 help="output image file name")
p.add_argument("-p", "--preliminary", dest="preliminary", action="store_true",
                 default=False,
                 help="Show PRELIMINARY on plot (default = false)")
p.add_argument("-P","--pass", dest="passx",  type=int, default = 4, 
                 help="Production Pass #")
p.add_argument("--noPass", dest="nopass", default = False, action="store_true", 
                 help="Do not show Pass # under PRELIMINARY")
p.add_argument("--showspectrum", dest="spectrum", action="store_true", 
                 default=False,
                 help="Don't show spectrum parameters on plot (default: show)")
p.add_argument("-d","--startdate",dest="startdate", type=str, 
                 default = "2000-01-01",
                 help="Start date 'YYYY-MM-DD' (default: first date in light "
                 +"curve file)")
p.add_argument("-D","--enddate",dest="enddate", type=str, 
                 default = "2050-01-01",
                 help="End date 'YYYY-MM-DD' (default: last date in light "
                 +"curve file)")
p.add_argument("--fluxlines", dest="fluxlines", type=str, default= "", 
                 help="file with flux lines (in TeV^-1 cm^-2 s^-1)")
p.add_argument("--culines", dest="culines", type=str, default= "", 
                 help="file with flux lines (in Crab units) to be added")
p.add_argument("--showdates", dest="showdates", type=str, 
                 default= "2015-01-01,2015-04-01,2015-07-01,2015-10-01,"
                         +"2016-01-01,2016-04-01,2016-07-01,2016-10-01," 
                         +"2017-01-01,2017-04-01,2017-07-01,2017-10-01",
                 help="comma-separated list of dates (YYYY-MM-DD) to be shown "
                 +"as vertical lines (default: every 3 months, beginning "
                 +"2015-01-01')")
p.add_argument("--figwidth", dest="figwidth", type=int, default= 18, 
                 help="Width of figure (default = 18)")
p.add_argument("--ncp", dest="ncpprior", type=float, default = 6.0, 
                 help="ncp-prior for Bayesian Blocks (default = 6.0)")
p.add_argument("--noBB", dest="noBB", action="store_true", default=False,
               help="Do not calculate/show Baysian Blocks (default = show)")
args = p.parse_args()

showcrab = args.showcrab

minCU = args.minflux
maxCU = args.maxflux
phThreshold = args.phThreshold
cuThreshold = args.cuThreshold
errorbars = args.errorbars
lcfile = args.lcfile

#information from the lc file name:

fname = lcfile.split("/")[-1]
ftitle= fname[:-4]
title = ftitle.split("_")
object = ""
for t in xrange(2,len(title)):
    if (title[t][:3] == "201"): break 
    object += "%s "%(title[t])


#argument start/end:
try:
    startd = args.startdate.split("-")    
    endd = args.enddate.split("-")    
except:
    print "wrong --startdate/--enddate format, use 'YYYY-MM-DD'"
    sys.exit()

startdate = dt.datetime(int(startd[0]),int(startd[1]),int(startd[2]))
mjd_begin = dt2mjd(startdate)
enddate = dt.datetime(int(endd[0]),int(endd[1]),int(endd[2]))
mjd_end = dt2mjd(enddate+dt.timedelta(days=1))

#loading the light curve:
hlc = LightCurve.HAWCint(lcfile, mjd_begin, mjd_end, logging = False)

redshift = hlc.redshift

mjd_begin = hlc.mjd_begin
mjd_end   = hlc.mjd_end
mjdrange  = mjd_end - mjd_begin

crabflux = hlc.crabflux
minFlux = minCU*crabflux
maxFlux = maxCU*crabflux

fluxes = hlc.flux
fluxes_err = hlc.flux_err
fluxesCU = hlc.fluxCU
fluxesCU_err = hlc.fluxCU_err
significances = sqrt(hlc.TS)

intThresh = hlc.intThreshold
#for photon flux, check if different threshold requested:
if ((phThreshold>-1.) & (phThreshold!=intThresh)):
    #convert to new threshold
    scale = integratedflux(ones(len(fluxes)),hlc.index,hlc.cutoff,phThreshold,redshift)
    scale /= integratedflux(ones(len(fluxes)),hlc.index,hlc.cutoff,intThresh,redshift)
    fluxes *= scale
    fluxes_err *= scale
    if isinstance(fluxes, ndarray):
        if not ( all(in1d(hlc.index,hlc.index[0])) & all(in1d(hlc.cutoff,hlc.cutoff[0])) ):
            print "Different spectra at different times prohibit simple scaling "
            +"between photon and CU fluxes. Will not show CU axis"
            showcrab = False
        else:
            minFlux*=scale[0]
            maxFlux*=scale[0]
    else:
        minFlux*=scale
        maxFlux*=scale
else:
    phThreshold = intThresh
#for photon flux, check if different threshold requested:
if ((cuThreshold>-1.) & (cuThreshold!=intThresh)):
    #recalculate based on simple PL for Crab
    scale = integratedflux(1,2.63,1000000,cuThreshold)
    scale /= integratedflux(1,2.63,1000000,intThresh)
    #crabflux *= scale
    fluxesCU *= scale
    fluxesCU_err *= scale
else:
    cuThreshold = intThresh

mjd = hlc.mjd
mjd_tlerr = hlc.mjd_tlerr
mjd_tuerr = hlc.mjd_tuerr

#consider max value of transits the sd=period:
sdperiod = max(hlc.transits)

#basic summary:
ntbins = fluxes.size
#ph
weights = 1./fluxes_err**2
avgflux = (weights*fluxes).sum()/weights.sum()
avgerr = average(fluxes_err)
#CU:
weightsCU = 1./fluxesCU_err**2
avgfluxCU = (weightsCU*fluxesCU).sum()/weightsCU.sum()
avgerrCU = average(fluxesCU_err)
print "Average flux:          %5.3e  ( %5.3f CU)"%(
        avgflux,avgfluxCU)
print "Average 1-sigma error: %5.3e  ( %5.3f CU)"%(
        avgerr,avgerrCU)

###########################################################
# the plot:
###########################################################

fig, ax_f = plt.subplots(figsize=(args.figwidth,6), tight_layout=True)
plt.xticks(fontsize=15) #, rotation=45)
plt.yticks(fontsize=15)

ax_f.set_ylim(minFlux,maxFlux)
ax_f.set_xlim(mjd_begin, mjd_end)

ax_f.set_xlabel('MJD [days]', fontsize=15)
ax_f.set_ylabel(r'Flux >%.0f TeV [ph cm$^{-2}$s$^{-1}$]'%(phThreshold), 
        fontsize=15)

plt.yticks(fontsize=15)


if (args.zero):
    ax_f.plot([mjd_begin-1,mjd_end+1],[0.,0.],color="grey",ls="-",alpha=0.4,
            linewidth=1.)

if (args.showavg):
    ax_f.plot([mjd_begin-1,mjd_end+1],[avgflux,avgflux],color="red",ls="--",
            alpha=0.5,linewidth=1.5,label='average flux')

if (args.showfixednorm):
    ax_f.plot([mjd_begin-1,mjd_end+1],[fixNormIntegrated,fixNormIntegrated],
            color="red",ls="--",alpha=0.8,linewidth=1.5,
            label='constant flux best fit')

if (args.spectrum):
    if (args.showdates!=""):
        sheight = maxFlux - (maxFlux-minFlux)*0.05
    else:
        sheight = minFlux+(maxFlux-minFlux)/80.
    if ((one_cutoff>100000) | (one_cutoff<=0)): 
        ax_f.text(mjd_begin+mjdrange*0.01,sheight,
              r'dF/dE = F$_0\cdot$(E/E$_{0}$)$^{-\gamma}$, $\gamma$=%.2f'%(
                  one_slope),color="black",fontsize=14)
    else:
        ax_f.text(mjd_begin+mjdrange*0.005,sheight,
              r'dF/dE = F$_0\cdot$(E/E$_{0}$)$^{-\gamma}\cdot '
              +'e^{-\mathrm{E/E_c}}$, $\gamma$=%.2f, E$_c$=%.1f TeV'%(
                  one_slope,one_cutoff),color="black",fontsize=14)


# external, additional plot elements:
# -----------------------------------

# events:
if (args.events!=""):
    evfile = open(args.events,"r")

    evmjd = []
    evwidth = []
    evtext = []
    evcolor = []
    evstyle = []
    evheight = []
    evtextmjd = []

    for row in evfile:
        try:
            field = row.split(";")
            thismjd = float(field[0])
            if ((thismjd<(mjd_begin-sdperiod[-1]/2.)) | (thismjd>mjd_end)):
                continue
            evmjd.append(thismjd)
            evw = float(field[1])
            evwidth.append(evw)
            evtext.append((field[2].strip()).split('\\\\'))
            
            if (len(field)>3):
                evcolor.append(field[3].strip())
            else:
                evcolor.append("black")
            if (len(field)>4):
                evheight.append(float(field[4].strip()))
            else:
                evheight.append(0.95)
            if (len(field)>5):
                evtextmjd.append(thismjd+float(field[5]))
            else:
                evtextmjd.append(thismjd+evw+1.)
        except:
                continue
    evfile.close()

    for i in xrange(0,len(evmjd)):
        ax_f.bar(evmjd[i],2.,evwidth[i],bottom=-1.,color=evcolor[i],alpha=0.4,
                linewidth=0,zorder=-1)
        for l in xrange(0,len(evtext[i])):
            ax_f.text(evtextmjd[i], 
                    minFlux+(maxFlux-minFlux)*(evheight[i]-l*0.045), 
                    '%s'%(evtext[i][l]),color=evcolor[i])

#flux lines
if (args.fluxlines!=""):
    flfile = open(args.fluxlines,"r")

    flval = []
    flup = []
    fllow = []
    fltext = []
    flcolor = []
    fltext_hor = []

    for row in flfile:
        try:
            field = row.split(";")
            f = float(field[0])
            ferr = float(field[1])
            fup = f  + ferr
            flow = f - ferr
            s = float(field[2])
            c = float(field[3])
            f = integratedflux(f,s,c,phThresh)
            fup = integratedflux(fup,s,c,phThresh)
            flow = integratedflux(flow,s,c,phThresh)
            if ((f<(minFlux)) | (f>maxFlux)):
                continue
            flval.append(float(f))
            flup.append(float(fup))
            fllow.append(float(flow))
            fltext.append((field[4].strip()).split('\\\\'))
            if (len(field)>5):
                flcolor.append(field[5].strip())
            else:
                flcolor.append("gray")
            if (len(field)>6):
                fltext_hor.append(float(field[6].strip()))
            else:
                fltext_hor.append(0.8)
        except:
                continue
    flfile.close()

    for i in xrange(0,len(flval)):
        ax_f.fill_between(array([mjd_begin-1,mjd_end+1]),
                array([fllow[i],fllow[i]]),array([flup[i],flup[i]]),
                color=flcolor[i],alpha=0.5)
        for l in xrange(0,len(fltext[i])):
            ax_f.text(mjd_begin+fltext_hor[i]*mjdrange, 
                    fllow[i]-(maxFlux-minFlux)*(1.5+l)*0.045, 
                    r'%s'%(fltext[i][l]),color=flcolor[i])

#cu lines
if (args.culines!=""):
    cufile = open(args.culines,"r")

    cuval = []
    cuup = []
    culow = []
    cutext = []
    cucolor = []
    cutext_hor = []

    for row in cufile:
        try:
            field = row.split(";")
            f = float(field[0])
            ferr = float(field[1])
            fup = f  + ferr
            flow = f - ferr
            if ((f<minCU) | (f>maxCU)):
                continue
            cuval.append(float(f))
            cuup.append(float(fup))
            culow.append(float(flow))
            cutext.append((field[2].strip()).split('\\\\'))
            if (len(field)>3):
                cucolor.append(field[3].strip())
            else:
                cucolor.append("gray")
            if (len(field)>4):
                cutext_hor.append(float(field[4].strip()))
            else:
                cutext_hor.append(0.8)
        except:
                continue
    cufile.close()

    for i in xrange(0,len(cuval)):
        ax_c.fill_between(array([mjd_begin-1,mjdrange]),
                array([culow[i],culow[i]]),
                array([cuup[i],cuup[i]]),color=cucolor[i],
                alpha=0.5)
        for l in xrange(0,len(cutext[i])):
            ax_f.text(mjd_begin+cutext_hor[i]*mjdrange, 
                    (culow[i]-(maxFlux-minFlux)*(1+l)*0.045), 
                    r'%s'%(cutext[i][l]),color=cucolor[i])

dates = args.showdates.split(",")
for d in dates:
    try:
        da = dt.datetime(int(d[:4]),int(d[5:7]),int(d[8:]))
    except:
        print "Wrong date format in --showdates: %s . Use YYYY/MM/DD"%(d)
        raise
    m = dt2mjd(da)
    nd = da.strftime("%-d %b. %Y") #d[:4]+"/"+d[5:7]+"/"+d[8:]
    ax_f.plot([m,m],[-1,1],color="black",alpha=0.5,lw=0.7)
    ax_f.text(m-mjdrange*0.011,minFlux+(maxFlux-minFlux)*0.18,nd,
            color="black",alpha=0.5,rotation="vertical",fontsize=12)

if (errorbars) :
    edgewidth=0.7
    ax_f.errorbar(mjd,fluxes,xerr=[mjd_tlerr,mjd_tuerr],yerr=fluxes_err,
            color="black",ls="",label="%d-transit HAWC flux"%(sdperiod),
            mew=edgewidth,lw=1,capsize=edgewidth,elinewidth=edgewidth)
else :
    msize = average(array(mjd_tlerr))*1.3
    medge = msize/2.3/1.3
    mewcorrection = 1.-0.03*medge
    if (msize<0.5):
        medge=0
    if (msize<3.): 
        msize=3.
    ax_f.errorbar(mjd,fluxes,
            xerr=[mjd_tlerr*mewcorrection,mjd_tuerr*mewcorrection],
            yerr=0,color="black",mec="black",ls="",
            label="HAWC %d-transit"%(sdperiod),marker="o",ms=msize,
            mew=medge,lw=1.5)
    ax_f.errorbar(lcg_d,lcg_v,
            xerr=[lcg_tlerr*mewcorrection,lcg_tuerr*mewcorrection],
            yerr=0,lolims=lcg_l,color="grey",mec="grey",ls="",
            marker="o",ms=msize,mew=medge,lw=1.5)




######################################################################
#Baysian Blocks
######################################################################

if not (args.noBB):
     
    ncp_prior = args.ncpprior

    (bbpoints, bbtimes, bbamps, bbampserr, cp, amplitudes, 
      amplitudes_err) = makeblocks(mjd,fluxes,fluxes_err,ncp_prior)
    ncp = len(bbpoints)
    N = len(fluxes)

    print "\n"
    print "Baysian Blocks:"
    print "---------------"
    print "Found %d blocks, at change points: "%(ncp)
    print bbpoints
    print "with durations (in days):"
    print bbtimes
    print "and amplitudes:"
    print bbamps
    print "with errors:"
    print bbampserr
    print "     MJD               duration     amplitude      error" 
    for c in xrange(0,len(bbpoints)):
        cp1 = bbpoints[c]
        if (c != (len(bbpoints)-1)):
            cp2 = bbpoints[c+1]-1
        else:
            cp2 = N-1
        print " %.3f - %.3f     %.3f          %5.3e         %5.3e " %(
                mjd[cp1]-0.125,mjd[cp2]+0.125,
                bbtimes[c],
                bbamps[c],
                bbampserr[c])                   

    ax_f.plot(mjd[cp],amplitudes,color="blue",lw=1,label = "Bayesian blocks")
    #plot errorbands
    ax_f.fill_between(mjd[cp],amplitudes-amplitudes_err,
            amplitudes+amplitudes_err,color="blue",alpha=0.4,linewidth=0)


  

#LC plotting:

if showcrab:
    ax_c = ax_f.twinx()
    ax_c.set_axisbelow(True)
    ax_c.set_ylim(minCU,maxCU)
    ax_c.set_ylabel(r'Flux >%.1f TeV [Crab Units]'%(cuThreshold), fontsize=15, labelpad=15)
    lclabels = ax_f.get_legend_handles_labels()
    ax_c.grid(axis="y")
    legend = ax_c.legend(lclabels[0],lclabels[1],loc='upper right',fontsize=15,numpoints=1)
else:
    legend = ax_f.legend(loc='upper right',fontsize=15,numpoints=1)

if (args.preliminary):
    ax_f.text(mjd_begin+mjdrange*0.35, minFlux+(maxFlux-minFlux)*0.92, 
            'PRELIMINARY',color="0.7",fontsize=22,fontweight="bold")
    passx = args.passx
    if not args.nopass:
        ax_f.text(mjd_begin+mjdrange*0.41, minFlux+(maxFlux-minFlux)*0.87, 
                'Pass %d'%(passx),color="0.7",fontsize=15,fontweight="bold")


legend.set_title(object)
legend.get_title().set_size(15)
legend.get_title().set_weight("bold")
plt.setp(legend.get_title(), x=-50)



#plotting or storing
if args.output:
     fig.savefig("LC_"+args.output, tight_layout=True, dpi=300)
     figFlux.savefig("histFlux_"+args.output, tight_layout=True)
     if not (args.noBB):
         figBB.savefig("histBB_"+args.output, tight_layout=True)
else:
     plt.show()


