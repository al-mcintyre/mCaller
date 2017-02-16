import numpy as np
import scipy.stats

class Sequence:
    def __init__(self,ntseq):
        self.seq = ntseq.upper()
        self.length = len(self.seq)
        
    def methylate(self,ind,replacement='M'):
        self.methseq = self.seq[:ind-1]+replacement+self.seq[ind:]
        return self.methseq 
    
    def rev_comp(self):
        comp = {'A':'T','C':'G','T':'A','G':'C','N':'N','M':'M'}
        comp_seq = [comp[nt] for nt in list(self.seq)]
        self.revcomp = ''.join([comp[nt] for nt in list(self.seq)][::-1])
        try: 
            self.revmethseq = ''.join([comp[nt] for nt in list(self.methseq)][::-1])
        except AttributeError:
            pass
        return self.revcomp

class Signal:
    #R7 example: (56.40406723484848, 8851.159030544488, 0.7926077052243929, 0.01095617529880478, 'GAAAAA', 57.16981709401989, 1, 0.9999896287918091, 0.31872811913490295, 'GAAAAA', 0.31872811913490295, 0.94491046667099, 0.012281710281968117, 0.04044119268655777, 0.002366444794461131)
    def __init__(self,events,shift,drift,scale,var,starttime):
        self.means = [x[0] for x in events]
        #self.stds = [x[2] for x in events]
        self.event_times = [x[1] for x in events]
        self.kmers = [x[4] for x in events]
        self.advances = [x[6] for x in events]
        self.read_start = starttime
        self.drift = drift
        self.shift = shift
        self.scale = scale
        self.var = var
        self.transformed = False
        
    def transform(self):
        """double time = event.start_time - events[si][0].start_time;
        event.mean -= (time * pore_model[si].drift); """   
        if not self.transformed:
            self.means = np.array(self.means)/self.scale - (np.array(self.event_times)-self.read_start)*self.drift - self.shift
            self.transformed = True

def check_alignment(indices,seq1,seq2):
    """verifies whether index of alignment contains indel or mismatch (read = seq1, ref = seq2)"""
    if all([index != None for index in indices]) and (seq1[indices[0]] == seq2[indices[1]]): 
        return True
    else:
        return False

def find_perfect_match(bam_alignment,k,seq1,seq2):
    """finds next closest perfect k-mer in alignment, avoiding any indels or mismatches (read = seq1, ref = seq2)"""
    match = []
    i = -1
    last_ind = 0
    """ 
    while len(match) < k and i<len(bam_alignment)-k:
        align_ind1,align_ind2 = bam_alignment[i]
        if all([base1 == base2 for base1,base2 in zip(seq1[align_ind1:align_ind1+k],seq2[align_ind2:align_ind2+k])]):
            match = bam_alignment[i:i+k]
    """
    while len(match) < k and i < len(bam_alignment)-1:
        i += 1
        if check_alignment(bam_alignment[i],seq1,seq2):
            if len(match) > 0 and i == last_ind+1:
                match.append(bam_alignment[i]) #,seq1[bam_alignment[i][0]],seq2[bam_alignment[i][1]]))
            else:
                match = [bam_alignment[i]] #,seq1[bam_alignment[i][0]],seq2[bam_alignment[i][1]])]
            last_ind = i   
 
    if len(match) == k:
        return match, last_ind
    else:
        return None, None

def extract_signal_fragment(anchor1,anchor2,hdf5,ind_mult,read_len):
    """returns signal between two perfect kmer anchor points"""
    ind = min(0,ind_mult*read_len)
    all_events = hdf5['Analyses/Basecall_1D_000/BaseCalled_template/Events']
    read_start = hdf5['Analyses/Basecall_1D_000/BaseCalled_template/Events'].attrs['start_time']
    
    location = 'Analyses/Basecall_1D_000/BaseCalled_template/Model'
    shift = hdf5[location].attrs['shift']
    drift = hdf5[location].attrs['drift']
    scale = hdf5[location].attrs['scale']
    var = hdf5[location].attrs['var']    
    event_frag = [] 
    
    #print ind, anchor1, anchor2
    for event in all_events:
        ind = ind + event[6]
        if ind_mult*ind >= anchor1 and ind_mult*ind < anchor2+1:
            #print ind_mult*ind, event
            event_frag.append(event)
    #print event_frag
    sig_frag = Signal(event_frag,shift,drift,scale,var,read_start)
    return sig_frag

def extract_transitions(hdf5,loc='Analyses/Basecall_1D_000/BaseCalled_template'):
    """returns array with transition probabilities for stays, steps, skips, and double skips"""
    location = loc+"/Model"
    stay_rate = hdf5[location].attrs['stay_prob']
    step_rate = hdf5[location].attrs['step_prob']
    skip_rate = hdf5[location].attrs['skip_prob']
    dskip_rate = np.power(hdf5[location].attrs['skip_prob'],2)
    return (stay_rate,step_rate,skip_rate,dskip_rate)

def emission_prob(kmer_model,obs_mean,sig):
    """returns probability density for an event given a Gaussian model for a kmer"""
    P = scipy.stats.norm(kmer_model[0],kmer_model[1]*sig.var).pdf(obs_mean)
    return P

def viterbiish(events,sequence,transition_probs,kmer_mod,k=6,m_sequence=None):
    """implements a variation on the Viterbi algorithm where the first and last states are known"""
    nStates = len(sequence)-k+1
    nSamples = len(events.means)
    states = [sequence[i:i+k] for i in range(nStates)] 
    
    if m_sequence:
        m_states = [m_sequence[i:i+k] for i in range(nStates)] 
        final_states = m_states
    else:
        final_states = states
    
    vit = np.zeros((nSamples+1,nStates)) # initialize viterbi table
    bpt = np.zeros((nSamples+1,nStates)) # initialize the best path table
    best_path = np.zeros(nSamples+1); # initialize output vector

    #initialization 
    vit[0,0] = 1 
    
    #first phase
    for j in xrange(1,nSamples):
        for i in xrange(0,nStates): 
            results_across_i2 = [0]*max(0,i-3) + [vit[j-1,itState2]*transition_probs[i-itState2]*emission_prob(kmer_mod.model[states[i]],events.means[j],events) for itState2 in range(max(0,i-3),min(nStates,i+1))] + [0]*(nStates-i)
            bpt[j,i] = np.argmax(results_across_i2)
            vit[j,i] = results_across_i2[int(bpt[j,i])]
    
    #termination step
    #vit[-1,-1] = 1 
    bpt[-1,-1] = nStates-1 
    
    best_path[-1] = bpt[-1,-1]
    for j in range(nSamples-1,-1,-1):
        best_path[j] = bpt[j,int(best_path[j+1])]
    
    return [(x,final_states[int(x)],events.means[i]) for i,x in enumerate(best_path[1:])]

def update_signal_matrix(realigned_sig,model):
    """adds differences between measured and expected currents for kmers surrounding a realigned base to the signal matrix"""
    A_signal = [x for x in realigned_sig if len(x[1].split('M')) == 2]
    
    kmer_currents = {i:[] for i in range(1,7)}
    
    for kmer in A_signal:
        base_ind = len(kmer[1].split('M')[0])+1
        orig_kmer = 'A'.join(kmer[1].split('M'))
        kmer_currents[base_ind].append(kmer[2]-model[orig_kmer][0])
            
    kmer_current_array = [0]*6
    for base_ind in kmer_currents:
        if len(kmer_currents[base_ind]) > 0:
            kmer_current_array[base_ind-1] = np.mean(kmer_currents[base_ind])
    
    return kmer_current_array
