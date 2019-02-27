#!/usr/bin/env python
import sys
import os
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
from sklearn.cluster import KMeans
from scipy.stats import pearsonr,percentileofscore
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score
import numpy as np
import warnings
import pandas as pd

base = 'A'
modbase = 'm6A'
base_colours = {base:'#55B196', modbase:'#B4656F'}

def plot_w_labels(klabels,labels,currents,strategy,kmer,pos,outdir,base_colours,train=False,alpha=1):
    warnings.filterwarnings("ignore", module="matplotlib")
    bin_labels = [1 if x == 'A' else 0 for x in labels]
    lstyles = {0:'-',1:'--',-1:':',2:':'}
    sns.set_style('white')
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    if train:
        ars = adjusted_rand_score(bin_labels[:-1], klabels[:-1])

    if len(set(klabels)) < 4:
        for current,label,kl in zip(currents,labels,klabels):
            plt.plot(range(1,7),current,label='{}, {}'.format(label,kl),color=base_colours[label],linestyle=lstyles[kl],alpha=alpha)

        plt.ylabel('observed-expected current (pA)')
        plt.xlabel('position in kmer')
        handles, labels = ax.get_legend_handles_labels()
        hs, ls = [],[]
        for h,l in zip(handles, labels):
            if l not in set(ls):
                ls.append(l)
                hs.append(h)
        ax.legend(hs,ls,loc='center left', bbox_to_anchor=(1, 0.5))
        title = kmer 
        if train:
            title = title + ', clustered by '+strategy+"\nAdjusted Rand Index: "+str(np.round(ars,3))
        plt.title(title)
        plt.show()
        plt.savefig(outdir+'/signals_cluster_'+str(pos)+'.pdf',dpi=500,bbox_inches='tight',transparent=True)
    if train:
        return ars

def plot_correlation_matrix(curmat,elevenmer,labels,outdir):
    plt.figure(figsize=(7,6))
    cg = sns.clustermap(curmat,metric='euclidean',xticklabels=labels,yticklabels=labels)
    plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    #sns.heatmap(curmat,xticklabels=labels,yticklabels=labels)
    plt.title(elevenmer)
    plt.show()
    plt.savefig(outdir+'correlation_matrix_'+elevenmer+'.pdf',dpi=500,transparent=True)

def plot_change_by_pos(diffs_by_context,plottype='box'):
    fig = plt.figure(figsize=(6,4))
    changes_by_position = {'position':[],'base':[],'diff':[]}
    for lab in diffs_by_context:
        for context in diffs_by_context[lab]:
            for entry in diffs_by_context[lab][context]:
                for pos,diff in enumerate(entry[:-1]):
                    changes_by_position['position'].append(pos+1)
                    changes_by_position['base'].append(lab)
                    changes_by_position['diff'].append(diff)
    dPos = pd.DataFrame(changes_by_position)
    if plottype == 'box':
        sns.boxplot(x="position", y="diff", hue="base", data=dPos, palette=[cols[base],cols[methbase]])
    elif plottype == 'violin':
        sns.violinplot(x="position",y="diff", hue="base", data=dPos, palette=[cols[base],cols[methbase]])
    sns.despine(trim=False)
    plt.xlabel('Adenine Position in 6-mer')
    plt.ylabel('Measured - Expected Current (pA)')
    plt.ylim([-20,20])
    plt.legend(title='',loc='upper center', bbox_to_anchor=(0.5, 1.05),
          ncol=3, fancybox=True)
    plt.savefig('change_by_position_box.pdf',transparent=True,dpi=500, bbox_inches='tight')

def plot_training_probabilities(prob_scores,tb):
    #prob_scores = {'m6A':[0.9,0.4,...],'A':[0.1,0.5,0.2,...]}
    sns.set_style('darkgrid')
    sns.set_palette(['#55B196','#B4656F'])
    fig = plt.figure(figsize=(3,4))
    prob_dict = {'probability':prob_scores[base]+prob_scores[modbase],'base':[base]*len(prob_scores[base])+[modbase]*len(prob_scores[modbase])}
    prob_db = pd.DataFrame(prob_dict)
    sns.boxplot(x="base", y="probability", data=prob_db)
    sns.despine()
    plt.show()
    plt.savefig('training_probability_'+tb+'_model_boxplot.pdf',transparent=True,dpi=500,bbox_inches='tight')
