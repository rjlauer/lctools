from numpy import *

######################################################################
#Baysian Blocks
######################################################################


def makeblocks(smjd,srate,srate_err,ncp_prior):

    N = len(srate)

    best = []
    last = []

    for R in xrange(1,N+1):
        #print "R = ",R
        #print "  srate[0:R] = ", srate[0:R]
        #print "  srate_err[0:R] = ", srate_err[0:R]
        sum_a = 1/2.*array([ sum (1./srate_err[k:R]**2.) for k in xrange(0,R)])
        sum_b = -array([ sum(srate[k:R]/srate_err[k:R]**2.) for k in xrange(0,R)])
        #print "  sum_a = ", sum_a
        #print "  sum_b = ", sum_b
        fit_vec = sum_b**2./(4.*sum_a)
        #print "  fit_vec = ", fit_vec
        last.append( argmax( concatenate([[0.], best]) + fit_vec - ncp_prior ) )
        #print "  last = ",last
        best.append( (concatenate([[0.], best]) + fit_vec -ncp_prior )[ last[R-1] ] )
        #print "  best = ",best

    s_index = last[N-1]
    s_ncp = 0
    #one entry per change point:
    s_bbpoints = array([],dtype=int)
    s_bbamps =  array([],dtype=float)
    s_bbampserr = array([],dtype=float)
    s_bbtimes = array([],dtype=float)
    # two entries per change point for start and stop:
    s_cp = array([N-1],dtype=int)
    s_amplitudes = array([],dtype=float)
    s_amplitudes_err = array([],dtype=float)
    while (s_index > 0):
        s_bbpoints = concatenate( [[s_index],s_bbpoints] ) 
        s_cp = concatenate( [ [s_index, s_index] , s_cp ] )
        #amplitudes
        sum_a = 1/2. * sum(1./srate_err[s_index:s_cp[-(1+2*s_ncp)]]**2.)
        sum_b = -sum(srate[s_index:s_cp[-(1+2*s_ncp)]]/srate_err[s_index:s_cp[-(1+2*s_ncp)]]**2.)
        amp = -sum_b/(2.*sum_a)
        s_bbamps = concatenate( [[amp], s_bbamps ] )
        s_amplitudes = concatenate( [ [ amp, amp ] , s_amplitudes ] )
        # errors from error on average: sqrt( sum (1/sigma**2) )
        aerr = 1./sqrt(2.*sum_a)
        s_bbampserr = concatenate( [[aerr], s_bbampserr ] )
        s_amplitudes_err = concatenate( [ [ aerr, aerr ] , s_amplitudes_err ] )

        if (s_ncp==0):
            s_bbtimes = concatenate([[smjd[N-1]-smjd[s_index]], s_bbtimes])
        else:
            s_bbtimes = concatenate([[smjd[s_bbpoints[1]]-smjd[s_index]], s_bbtimes])
        s_ncp += 1
        s_index = last[s_index -1 ]
    s_bbpoints = concatenate([[0],s_bbpoints])
    s_cp = concatenate( [ [0] , s_cp] ) 
    sum_a = 1/2. * sum(1./srate_err[0:s_cp[1]]**2.)
    sum_b = -sum(srate[0:s_cp[1]]/srate_err[0:s_cp[1]]**2.)
    firstamp = -sum_b/(2.*sum_a)
    s_amplitudes = concatenate( [ [ firstamp,firstamp ] , s_amplitudes ] )
    s_bbamps = concatenate( [[firstamp],s_bbamps ] )
    firstamp_err = 1./sqrt(2.*sum_a)
    s_amplitudes_err = concatenate( [ [ firstamp_err ,firstamp_err ] , s_amplitudes_err ] )
    s_bbampserr = concatenate( [[firstamp_err],s_bbampserr ] )
    s_bbtimes = concatenate([[smjd[s_cp[1]]-smjd[0]],s_bbtimes])
    #RMS
    #s_bbRMS.append(sqrt(1/(len(srate[s_index:s_cp[-(1+2*s_ncp)]])-1)*sum(srate[s_index:s_cp[-(1+2*s_ncp)]]**2.)))

    s_bbamps = array(s_bbamps)
    s_bbampserr = array(s_bbampserr)
    s_bbtimes = array(s_bbtimes)
    s_cp = array(s_cp)
    return (s_bbpoints,s_bbtimes,s_bbamps,s_bbampserr,s_cp,s_amplitudes,s_amplitudes_err)

