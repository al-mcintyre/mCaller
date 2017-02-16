import cython
cimport numpy as np
import numpy as np
import scipy.stats

cdef inline int int_max(int a, int b): return a if a >= b else b
cdef inline int int_min(int a, int b): return a if a <= b else b

def cython_viterbiish(np.ndarray[double,ndim=1] event_means, float event_var, str sequence, np.ndarray[double,ndim=1] transition_probs, dict kmer_model, int k, str m_sequence):
    def emission_prob(kmer_model,obs_mean,var):
        """returns probability density for an event given a Gaussian model for a kmer"""
        P = scipy.stats.norm(kmer_model[0],kmer_model[1]*var).pdf(obs_mean)
        return P

    cdef int nStates = len(sequence)-k+1
    cdef int nSamples = len(event_means)
    cdef int i, j
    cdef object states
    states = [sequence[i:i+k] for i in range(nStates)]

    if m_sequence:
        m_states = [m_sequence[i:i+k] for i in range(nStates)] 
        final_states = m_states
    else:
        final_states = states
    
    cdef double[:,:] vit = np.zeros((nSamples+1,nStates)) # initialize viterbi table
    cdef int[:,:] bpt = np.zeros((nSamples+1,nStates),dtype=np.intc) # initialize the best path table
    cdef int[:] best_path = np.zeros(nSamples+1,dtype=np.intc) # initialize output vector

    #initialization 
    vit[0,0] = 1 

    cdef double[:] results_across_i2 = np.zeros(nSamples+1)
    
    #first phase
    for j in xrange(1,nSamples):
        for i in xrange(0,nStates):
            start_ind = int_max(0,i-3)
            end_ind = int_min(nStates,i+1)
            for ii in xrange(0,start_ind):
                results_across_i2[ii] = 0
            for ii in xrange(start_ind,end_ind):
                results_across_i2[ii] = vit[j-1,ii]*transition_probs[i-ii]*emission_prob(kmer_model[states[i]],event_means[j],event_var)
            bpt[j,i] = np.argmax(results_across_i2)
            vit[j,i] = results_across_i2[bpt[j,i]]
    
    #termination step
    vit[-1,-1] = 1 
    bpt[-1,-1] = nStates-1 
    
    best_path[-1] = bpt[-1,-1]
    for j in range(nSamples-1,-1,-1):
        best_path[j] = bpt[j,best_path[j+1]]
    
    return [(x,final_states[x],event_means[i]) for i,x in enumerate(best_path[1:])]
