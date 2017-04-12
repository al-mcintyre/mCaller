from classifier import *
import sys
import random
import cPickle
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from decimal import Decimal
from collections import defaultdict
from sklearn.ensemble import RandomForestClassifier
from scipy.stats.stats import pearsonr # use spearman instead of pearson
from scipy.stats import spearmanr 
from sklearn.metrics import roc_auc_score

def extract_data(exp,num,i,read_qualities,qual_thresh,skip_thresh,diffs_by_context,diff_fis=None):
   signals, contexts, labels, positions = [],[],[],[]
   if not diff_fis:
      diff_fis = [exp+'_R9_fwd.gm.ecoli.eventalign.diffs.'+num+'m6A',exp+'_R9_fwd.gm.ecoli.eventalign.diffs.'+num+'A']
   for diff_fi in diff_fis:
      for line in open(diff_fi,'r'):
         readname, pos, context, diffs, strand, label = line.split('\t')
         if read_qualities[readname] >= qual_thresh:
            nskips = len([float(x) for x in diffs.split(',') if float(x) == 0])
            features = [float(x) for x in diffs.split(',')]+[read_qualities[readname]]
            if len(features) == i and nskips <= skip_thresh:
               if context not in diffs_by_context[label.strip()]:
                  diffs_by_context[label.strip()][context] = [features]
               else:
                  diffs_by_context[label.strip()][context].append(features)
               signals.append(features)
               contexts.append(context)
               labels.append(label.strip())
               positions.append(pos)
      print diff_fi, len(signals)   
   for label in diffs_by_context:
      print len(diffs_by_context[label]),'contexts for',label
   print len(set(diffs_by_context[unmeth].keys())&set(diffs_by_context[meth].keys())),'overlapping contexts'
   print len(signals), 'observations' #, eg.',training_signals[:5]
   return diffs_by_context, signals, contexts, labels, positions

def plot_diffs_by_context(diffs_by_context,base,methbase,cols):
   count = 0
   sns.set_style('white')
   for context in diffs_by_context[methbase]:
      if context in diffs_by_context[base]:
         fig = plt.figure(figsize=(6,4))
         count += 1
         for label in diffs_by_context:
            col = cols[label]
            for i,entry in enumerate(diffs_by_context[label][context]):
               if i == 0:
                  plt.plot(range(1,len(entry[:-1])+1),entry[:-1],color=col,alpha=1,label=label)
               else:
                  plt.plot(range(1,len(entry[:-1])+1),entry[:-1],color=col,alpha=1)
         plt.title(context)
         plt.ylim([-8,8])
         sns.despine(trim=False)
         plt.ylabel('Measured - Expected Current (pA)')
         plt.xlabel('Adenine Position in 6-mer')
         plt.legend()
         plt.savefig('R9_measured_vs_expected_across_positions_'+context+'.pdf',dpi=500,bbox_inches='tight',transparent=True)
      if count > 10:
         break

   fig = plt.figure(figsize=(6,4))
   for lab in [meth,unmeth]:
      col = cols[lab]
      all_contexts = diffs_by_context[lab].keys()
      contexts = random.sample(all_contexts,500)
      for i,context in enumerate(contexts):
        for entry in random.sample(diffs_by_context[lab][context],1):
          if i == 0:
             plt.plot(range(1,len(entry[:-1])+1),entry[:-1],color=col,alpha=0.2,label=lab)
          else:
             plt.plot(range(1,len(entry[:-1])+1),entry[:-1],color=col,alpha=0.2)
   plt.ylabel('Measured - Expected Current (pA)')
   plt.xlabel('Adenine Position in 6-mer')
   plt.legend()
   sns.despine(trim=False)
   plt.ylim([-15,15])
   plt.savefig('R9_measured_vs_expected_across_positions.pdf',transparent=True,dpi=500,bbox_inches='tight')

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
   sns.boxplot(x="position", y="diff", hue="base", data=dPos, palette=[cols[base],cols[methbase]])
   sns.despine(trim=False)
   plt.xlabel('Adenine Position in 6-mer')
   plt.ylabel('Measured - Expected Current (pA)')
   plt.ylim([-20,20])
   plt.legend(title='',loc='upper center', bbox_to_anchor=(0.5, 1.05),
          ncol=3, fancybox=True)
   plt.savefig('R9_change_by_position_box.pdf',transparent=True,dpi=500, bbox_inches='tight')

   sys.exit(0)


unmeth = 'A'
meth = 'm6A'
roc = {}
prob_scores = {}
nums = ['4.','6.','8.']
num_feats = [5,7,9]
classifiers = ['NN','RF','LR','NBC']
base_colours = {unmeth:'#55B196', meth:'#B4656F'}

for exp in ['chiu','mason'][1:]:

   read_qualities = {x.split('\t')[0]:float(x.split('\t')[1]) for x in open(exp+'_R9_fwd_qual_list.txt','r').read().split('\n') if len(x.split('\t'))>1}

   #fo qual, qual_threshold in zip(['hq.',''],[9,0]):
   for num,i in zip(nums,num_feats):
      if num == '6.':
         quals,qthreshs = ['_hq',''],[9,0]
         skips,sthreshs = ['_sk',''],[2,0]
      else:
         quals,qthreshs = [''],[0]
         skips,sthreshs = [''],[0]
      for qual, qual_threshold in zip(quals,qthreshs):
       for skip, skip_threshold in zip(skips,sthreshs):
         print exp, num
         diffs_by_context = {unmeth:{},meth:{}}
         diffs_by_context,signals,contexts,labels,positions = extract_data(exp,num,i,read_qualities,qual_threshold,skip_threshold,diffs_by_context)
         if exp == 'chiu' and num == '6.' and qual == '' and skip == '':
            plot_diffs_by_context(diffs_by_context,unmeth,meth,base_colours)
         if qual == '' and skip == '':
            methods = classifiers[:]
         else:
            methods = ['NN']
         for clf in [x+qual+skip for x in methods]:
         #for clf in ['SVM']: #,'SVM']: #RF','SVM','NN','LR','NBC']:
           print num, clf
           if exp == 'chiu':
              sys.stdout.flush()
              model = model_signal(signals,labels,True,contexts,modelfile='model_'+str(i)+'_'+clf+'_chiu.pkl',classifier=clf)      
              sys.stdout.flush()
           if exp == 'mason':
              if clf not in prob_scores:
                 prob_scores[clf] = {}
              prob_scores[clf][num],results = model_signal(signals,labels,False,contexts,modelfile='model_'+str(i)+'_'+clf+'_chiu.pkl')

              accuracy = []
              for base in prob_scores[clf][num]:
                 if base == 'm6A':
                    ind = 1
                 else:
                    ind = 0
                 accuracy.append( len([x for x in prob_scores[clf][num][base] if x[ind] > 0.5])*100./len(prob_scores[clf][num][base]) )
                 print base, accuracy[-1]
              print 'accuracy', np.mean(accuracy) 
             
              print 'auc', roc_auc_score([1 if lab == meth else 0 for lab in labels], [x[1] for x in results])

              sns.set_style('darkgrid')
              sns.set_palette(['#55B196','#B4656F'])
              fig = plt.figure(figsize=(3,4))
              prob_dict = {'probability':[x[1] for x in prob_scores[clf][num]['A']]+[x[1] for x in prob_scores[clf][num]['m6A']],'base':['A']*len(prob_scores[clf][num]['A'])+['m6A']*len(prob_scores[clf][num]['m6A'])}
              prob_db = pd.DataFrame(prob_dict)
              sns.boxplot(x="base", y="probability", data=prob_db)
              sns.despine()
              plt.savefig('R9_mason_trained_chiu_qual.gm.'+clf+'.'+num+'pdf',transparent=True,dpi=500,bbox_inches='tight')
              plt.close()

              if clf not in roc:
                 roc[clf] = {}
              roc[clf][num] = {'tp':[],'fp':[]}

              #print prob_scores['A'][:10]
              for thresh in np.arange(0,1.1,0.1)[::-1]:
                 roc[clf][num]['tp'].append(len([x for x in prob_scores[clf][num]['m6A'] if x[0]<=thresh])*100./len(prob_scores[clf][num]['m6A']))
                 roc[clf][num]['fp'].append(len([x for x in prob_scores[clf][num]['A'] if x[0]<=thresh])*100./len(prob_scores[clf][num]['A']))
              #print roc[clf][num]['tp'], roc[clf][num]['fp'], np.max(prob_scores['m6A']), np.max(prob_scores['A']) 
         
              pos_dict = {unmeth:{},meth:{}}
              if clf == 'NN' and num == '6.':
                 roc[clf+'_pos'] = {}
                 roc[clf+'_pos'][num] = {'tp':[],'fp':[]}
                 for pos,lab,res in zip(positions,labels,results):
                    if pos not in pos_dict[lab]:
                       pos_dict[lab][pos] = []
                    if res[0] > 0.5:
                       read_score = 1
                    else:
                       read_score = 0
                    pos_dict[lab][pos].append(read_score)
                 new_pos_dict = {meth:{},unmeth:{}}
                 for lab in pos_dict:
                    for pos in pos_dict[lab]:
                       if len(pos_dict[lab][pos]) >= 15: #set coverage threshold to 15
                          new_pos_dict[lab][pos] = pos_dict[lab][pos]
                 pos_dict = new_pos_dict
                 pos_labels = [0]*len(pos_dict[meth])+[1]*len(pos_dict[unmeth])
                 pos_probs = [np.mean(pos_dict[meth][pos]) for pos in pos_dict[meth]] + [np.mean(pos_dict[unmeth][pos]) for pos in pos_dict[unmeth]] 
                 print 'auc with 15X coverage', roc_auc_score(pos_labels, pos_probs)
                 for thresh in np.arange(0,1,0.1)[::-1]:
                    roc[clf+'_pos'][num]['tp'].append(len([pos for pos in pos_dict[meth] if np.mean(pos_dict[meth][pos])<=thresh])*100./len(pos_dict[meth]))
                    roc[clf+'_pos'][num]['fp'].append(len([pos for pos in pos_dict[unmeth] if np.mean(pos_dict[unmeth][pos])<=thresh])*100./len(pos_dict[unmeth]))
                    if thresh == 0.5:
                       print 'accuracy with 15X coverage', (roc[clf+'_pos'][num]['tp'][-1]+(100-roc[clf+'_pos'][num]['fp'][-1]))/2.
                 print roc[clf+'_pos'][num]['tp']
                 print roc[clf+'_pos'][num]['fp']

if exp == 'chiu': #or exp == 'mason':
   sys.exit(0)

sns.set_style('white')
#nums = ['6.','10.','15.']
num_lab ={'4.':'7','8.':'15','':'11','6.':'11','10.':'19','15.':'29'}
clf_colours = {'NN':['#D15C9C','#B70064','#750040'],'SVM':['#E1ADE1','#C382C3','#754E75'],'NN_hq':['#E1ADE1','#C382C3','#754E75'],'RF':['#96B596','#5B8C5A','#3A5A3A'],'NBC':['#B9B9BB','#767677','#000009'],'LR':['#87BCDE','#6F9AB6','#3E5665'],'NN_sk':['#FF5C72','#FF0022','#A30016'],'NN_pos':['#A19FB6','#6C698D','#45435A']}
for clf in [x for x in roc.keys() if len(x.split('_')) == 1]:
   fig = plt.figure(figsize=(4,4))
   for num_ind,num in enumerate(nums):
      plt.plot(roc[clf][num]['fp'],roc[clf][num]['tp'],color=clf_colours[clf][num_ind],label=num_lab[num])
   plt.legend(loc=4)
   plt.tight_layout()
   plt.ylabel('True Positives')
   plt.xlabel('False Positives')
   plt.savefig('R9_mason_trained_chiu.'+clf+'.ROC.pdf',dpi=500,transparent=True,bbox_inches='tight')
   plt.close()

for num_ind,num in enumerate(nums):
   fig = plt.figure(figsize=(4,4))
   for clf in [x for x in roc.keys() if len(x.split('_')) == 1]:
      plt.plot(roc[clf][num]['fp'],roc[clf][num]['tp'],color=clf_colours[clf][num_ind],label=clf)
   plt.legend(loc=4)
   plt.tight_layout()
   plt.ylabel('True Positives')
   plt.xlabel('False Positives')
   plt.savefig('R9_mason_trained_chiu.'+num+'ROC.pdf',dpi=500,transparent=True,bbox_inches='tight')
   plt.close()

num = '6.'
num_ind = 1
plt.figure(figsize=(4,4))
for clf in ['NN','NN_hq','NN_sk','NN_pos']:
   plt.plot(roc[clf][num]['fp'],roc[clf][num]['tp'],color=clf_colours[clf][num_ind],label=clf)
plt.legend(loc=4)
plt.tight_layout()
plt.ylabel('True Positives')
plt.xlabel('False Positives')
plt.savefig('R9_mason_trained_chiu.ROC.pdf',dpi=500,transparent=True,bbox_inches='tight')
plt.close()

plt.figure(figsize=(4,4))
subsample = {}
for base in (meth,unmeth):
   subsample[base] = set(random.sample(range(len(prob_scores['RF'][num][base])),10000))
for base,col in zip([meth,unmeth],['#B4656F','#55B196']):
   yax = [x[1] for i,x in enumerate(prob_scores['RF'][num][base]) if i in subsample[base]]
   xax = [x[1] for i,x in enumerate(prob_scores['NN'][num][base]) if i in subsample[base]]
   plt.scatter(xax,yax,label=base,color=col,alpha=0.2)  
rfvalues = [x[1] for x in prob_scores['RF'][num][base] for base in [meth,unmeth]]
nnvalues = [x[1] for x in prob_scores['NN'][num][base] for base in [meth,unmeth]]
print rfvalues[:10]
spearcorr = spearmanr(rfvalues,nnvalues)
plt.title('Spearman correlation = '+str(np.round(spearcorr[0],2))+', p = ' + str(np.round(spearcorr[1],2)))
plt.ylim([0,1])
plt.xlim([0,1])
plt.legend(loc=4)
plt.xlabel('m6A probability estimates (NN)')
plt.ylabel('m6A probability estimates (RF)')
plt.savefig('NN_vs_RF_prob_scores_'+num+'png',dpi=500,transparent=True,bbox_inches='tight')

#plot fractions of methylation called by ONT and PacBio 
clf = 'NN'
i = 7
exp = 'mason'
fraction_by_pos = {}
fraction_fi = open('ecoli_fractional_m6A_positions.txt','r').read().split('\n')
for pos in fraction_fi:
   if len(pos) > 0:
      pos_info = pos.split('\t')
      fraction_by_pos[int(pos_info[0])-1] = float(pos_info[2])

read_quals = {x.split('\t')[0]:float(x.split('\t')[1]) for x in open(exp+'_R9_fwd_qual_list.txt','r').read().split('\n') if len(x.split('\t'))>1} 
fdiffs_by_context = {'fracm6A':{},unmeth:{},meth:{}}
fdiffs_by_context, fsignals, fcontexts, flabels, fpositions = extract_data('mason','6.',7,read_quals,0,0,fdiffs_by_context,diff_fis=['mason_R9_fwd.gm.ecoli.eventalign.diffs.6.fracm6A'])
frac_predictions = model_signal(fsignals,train=False,modelfile='model_'+str(i)+'_'+clf+'_chiu.pkl')

meth_by_pos = {}
for pos,meth in zip(fpositions,frac_predictions):
   if pos not in meth_by_pos:
      meth_by_pos[pos] = {'m6A':0,'tot':0}
   if meth[1] >= 0.5:
      meth_by_pos[pos]['m6A'] += 1
   meth_by_pos[pos]['tot'] += 1

#print fraction_by_pos
ont_meth = []
pb_meth = []
for pos in meth_by_pos:
   if meth_by_pos[pos]['tot']>=15:
   #if int(pos) in fraction_by_pos:
      ont_meth.append(meth_by_pos[pos]['m6A']*1./meth_by_pos[pos]['tot'])
      pb_meth.append(fraction_by_pos[int(pos)])
   #else:
   #   print int(pos)+1 
corrcoef = pearsonr(pb_meth,ont_meth)

fig = plt.figure(figsize=(4,4))
plt.scatter(pb_meth,ont_meth,color='#B70064')
plt.ylabel('Fraction methylated (ONT)')
plt.xlabel('Fraction methylated (PacBio)')
plt.ylim([0,1])
plt.xlim([0,1])
plt.title('Pearson Correlation = ' +str(np.round(corrcoef[0],2))+', p = '+'%.2E' % Decimal(corrcoef[1])) #str(np.round(corrcoef[1],2)))
print 'Pearson Correlation = ' +str(np.round(corrcoef[0],2))+', p = ',corrcoef[1]
plt.savefig('R9_fraction_methylated.png',dpi=500,transparent=False, bbox_inches='tight')
