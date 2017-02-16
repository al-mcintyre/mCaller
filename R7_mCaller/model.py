class Model:
    """ from nanopolish:       // as per ONT documents
        scaled_states[i].level_mean = states[i].level_mean * scale + shift;
        scaled_states[i].level_stdv = states[i].level_stdv * var;
        scaled_states[i].sd_mean = states[i].sd_mean * scale_sd;
        scaled_states[i].sd_lambda = states[i].sd_lambda * var_sd;
        scaled_states[i].update_sd_stdv();"""
    def __init__(self,model_list):
        self.model = {}
        self.kmers = [k[0] for k in model_list]
        self.means = [k[1] for k in model_list]
        self.stds = [k[2] for k in model_list]
        self.model = {kmer:(mean,std) for kmer,mean,std in zip(self.kmers,self.means,self.stds)}

def extract_model(hdf5,loc='Analyses/Basecall_1D_000/BaseCalled_template'):
    """returns template strand ONT model as k-mer dict"""
    location = loc+"/Model"
    try:
        kmer_list = [(a[0],a[1],a[2]) for a in hdf5[location].value]
        mod = Model(kmer_list)
        return mod
    except KeyError:
        return None


