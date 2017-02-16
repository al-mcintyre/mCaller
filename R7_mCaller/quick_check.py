from classifier import *
import cPickle
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from collections import defaultdict
from sklearn.model_selection import train_test_split

t_signal_by_kmer = defaultdict(list)
egfi = open('ecoli_20160907_m6A.pkl','rb')
t_signals,t_labels,t_11mers = zip(*cPickle.load(egfi))
t_labels = ['m6A']*len(t_signals)
#print t_signals[:10],t_labels[:10],t_11mers[:10]
#t_prob_scores = random_forest(t_signals,[],False) 
egfi.close()

c_signal_by_kmer = defaultdict(list)
egfi = open('ecoli_20160907_A.pkl','rb')
c_signals,c_labels,c_11mers = zip(*cPickle.load(egfi))
c_labels = ['A']*len(c_signals)
#c_prob_scores = random_forest(c_signals,[],False)
egfi.close()

egfi = open('lambda_20160907_A.pkl','rb')
cl_signals,cl_labels,cl_11mers = zip(*cPickle.load(egfi))
cl_labels = ['A']*len(cl_signals)
#c_prob_scores = random_forest(c_signals,[],False)
egfi.close()


random_forest(t_signals+c_signals+cl_signals,t_labels+c_labels+cl_labels,True,list(t_11mers)+list(c_11mers)+list(cl_11mers))

"""
for signal,kmer in zip(t_signals,t_11mers):
    t_signal_by_kmer[kmer].append(signal)

for signal,kmer in zip(c_signals,c_11mers):
    c_signal_by_kmer[kmer].append(signal)

print len(set(t_11mers)),'methylated contexts and',len(set(c_11mers)),'unmethylated contexts are represented in this dataset'
kmer_intersection = set(t_11mers).intersection(set(c_11mers))

sns.set_style('white')
print len(kmer_intersection)
for kmer in list(kmer_intersection)[:20]:
    fig = plt.figure(figsize=(6,6))
    for mat,col,lab in zip([t_signal_by_kmer[kmer],c_signal_by_kmer[kmer]],['blue','red'],['m6A','A']):
            for i, entry in enumerate(mat):
                if i == 0:
                    plt.plot(entry,color=col,alpha=0.2,label=lab)
                else:
                    plt.plot(entry,color=col,alpha=0.2)
    plt.ylabel('Measured - Expected Current')
    plt.xlabel('Adenine Position in K-mer')
    plt.legend()
    #plt.show()
    plt.savefig('20160907_ecoli_measured_vs_expected_across_positions_'+kmer+'.png',dpi=500)    

c_avg_A = np.mean([x[0] for x in c_prob_scores])
c_avg_m6A = np.mean([x[1] for x in c_prob_scores])
c_m6A_probs = [x[1] for x in c_prob_scores]
t_avg_A = np.mean([x[0] for x in t_prob_scores])
t_avg_m6A = np.mean([x[1] for x in t_prob_scores])
t_m6A_probs = [x[1] for x in t_prob_scores]
print c_avg_A, c_avg_m6A
print t_avg_A, t_avg_m6A

print('unmethylated accuracy =', len([x for x in c_prob_scores if x[0] > 0.5])*100./len(c_prob_scores), len(c_prob_scores))
print('methylated accuracy =', len([x for x in t_prob_scores if x[0] < 0.5])*100./len(t_prob_scores), len(t_prob_scores))

prob_dict = {'probability':t_m6A_probs+c_m6A_probs,'base':['m6A']*len(t_m6A_probs)+['A']*len(c_m6A_probs)}
prob_db = pd.DataFrame(prob_dict)
sns.boxplot(x="base", y="probability", data=prob_db, palette="coolwarm")
sns.despine()
plt.savefig('20160907_ecoli_m6A_call_probability_rf.png',dpi=500)

fig = plt.figure(figsize=(6,6))
sns.set_style('white')
for mat,col,lab in zip([t_signals,c_signals],['blue','red'],['m6A','A']):
    for i,entry in enumerate(mat):
        if i == 0:
            plt.plot(entry,color=col,alpha=0.05,label=lab)
        else:
            plt.plot(entry,color=col,alpha=0.05)
plt.ylabel('Measured - Expected Current')
plt.xlabel('Adenine Position in K-mer')
plt.legend()
plt.savefig('20160907_ecoli_measured_vs_expected_across_positions.png',dpi=500)
"""
