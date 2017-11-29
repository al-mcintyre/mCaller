#!/usr/bin/env python
import sys
import numpy as np
from Bio import SeqIO
from scipy import stats
from scipy.cluster.hierarchy import fcluster,linkage
import scipy.spatial.distance as ssd
from extract_contexts import revcomp
from unsupervised_mod_detection import plot_w_labels
import os
import pandas as pd

def make_pos_set(pos_list):
    pos_set = set()
    with open(pos_list,'r') as fi:
        for line in fi:
            if len(line) > 3:
                pos_set.add(tuple(line.strip().split('\t')[:4]))
    return pos_set

def check_thresh(locus_list,mod_thresh,depth_thresh,control):
    if len(locus_list) >= depth_thresh:
        if not control and np.mean(locus_list) >= mod_thresh:
            return True
        elif control and np.mean(locus_list) < mod_thresh:
            return True
        else:
            return False

def write_gff(outfi,line_info):
    #scf7180000000004|quiver kinModCall  m4C 6644    6644    20  +   .   coverage=16;context=CGAAATCATTTCTCAGCAGGCGGTAAAAATGGCGGTTTTCG;IPDRatio=3.13;frac=1.000;fracLow=1.000;fracUp=1.000;identificationQv=7
    csome,nextpos,strand,deets = line_info
    line = '\t'.join([csome,"kinModCall","m6A",nextpos,nextpos,'10',strand,'.',deets])+'\n' #10 for all positions, giving a "log10 score" of 1 per read for motif calling with MotifMaker 
    outfi.write(line)

def ref2context(ref,pos_dict):
    ref_dict = {}
    pos_context_dict = {}
    for contig in SeqIO.parse(ref,"fasta"):
        ref_dict[contig.id] = str(contig.seq)
    for pos in pos_dict:
        if pos[0] in ref_dict:
            cx = ref_dict[pos[0]][int(pos[1])-20:int(pos[1])+21].upper()
            if pos[4] == '-':
                cx = revcomp(cx)
            pos_context_dict[pos] = cx #make adjustable length of context?
    ref_dict = {}
    return pos_context_dict

def cluster(values,context,original_labels,pos1,plot,plotdir):
    #print context, pos1
    #print values
    pdistance = ssd.pdist(values,metric='correlation')
    #print len(pd)
    dm = ssd.squareform(pdistance)
    #print dm
    link = linkage(dm,method='complete',metric='correlation')
    klabels = fcluster(link,2,'maxclust') #1,'inconsistent') #2,'maxclust')
    print klabels
    #klabels = [1 if x == 1 else 0 for x in klabels]
    #labels = ['m6A']*len(klabels)
    currents = values
    strategy = 'correlation'
    colours = {'m6A':'#B4656F','A':'#55B196'} #TODO update for other labels
    if plot:
        plot_w_labels(klabels,original_labels,currents,strategy,context+', position '+str(pos1),plotdir,colours)
    #for cluster in clusters:

def aggregate_by_pos(meth_fi,aggfi,depth_thresh,mod_thresh,pos_list,control,verbose_results,gff,ref,plot,plotdir):
    pos_dict = {}
    if verbose_results:
        pos_dict_verbose = {}

    if pos_list:
        pos_set = make_pos_set(pos_list)

    values_dict = {}
    for line in open(meth_fi,'r'):
        #try:
            #print line
            try:
                csome,read,pos,context,values,strand,label,prob = tuple(line.split('\t'))
            except: #for backwards compatibility; does not work with verbose results
                csome,read,pos,context,values,strand,label = tuple(line.split('\t'))
            nextpos = str(int(pos)+1)
            if pos_list and (csome,pos,nextpos,strand) not in pos_set:
                continue
            if (csome,pos,nextpos,context,strand) not in pos_dict:
                pos_dict[(csome,pos,nextpos,context,strand)] = []
                values_dict[(csome,pos,nextpos,context,strand)] = []
                if verbose_results:
                    pos_dict_verbose[(csome,pos,nextpos,context,strand)] = []
            if pos_list and (csome,pos,nextpos,strand) in pos_set:
                values_dict[(csome,pos,nextpos,context,strand)].append([float(v) for v in values.split(',')][:-1])
            if label[0] == 'm':
                pos_dict[(csome,pos,nextpos,context,strand)].append(1)
            else:
                pos_dict[(csome,pos,nextpos,context,strand)].append(0)
            if verbose_results:
                pos_dict_verbose[(csome,pos,nextpos,context,strand)].append(prob.strip())
        #except:
        #    pass

    if pos_list:
        for locus in values_dict:
            if plot:
                cluster(values_dict[locus],locus[3],['m6A' if x == 1 else 'A' for x in pos_dict[locus]],locus[1],plot,plotdir)
            values_df = pd.DataFrame(values_dict[locus])
            tvals = []
            pvals = []
            for i in values_df.columns: #[:-1]:
                #for j in values_df.columns[i+1:]:
                    ttest = stats.ttest_1samp(values_df[i],0)
                    #ttest = stats.ttest_rel(values_df[i],values_df[i+1])
                    #tvals.append(ttest[0])
                    pvals.append((ttest[1],ttest[0]))
            pval = (sum([-np.log10(x[0]) for x in pvals]), max([x[1] for x in pvals])) #min(pvals)
            values_dict[locus] = [np.round(x,3) for x in [pval[1],pval[0]]]

    if ref:
        context_dict = ref2context(ref,pos_dict)

    count = 0
    outfi = open(aggfi,'w')
    for locus in pos_dict.keys():
        a = (not pos_list) and check_thresh(pos_dict[locus],mod_thresh,depth_thresh,control)
        b = pos_list and (locus[0],locus[1],locus[2],locus[4]) in pos_set #and #'A' not in set(locus[4])
        if ref:
            cx = context_dict[locus]
        else:
            cx = locus[3]
        if a or b:
            count+=1 
            frac = np.mean(pos_dict[locus])
            if gff:
                deets = 'coverage='+str(len(pos_dict[locus]))+';context='+cx+';IPDRatio=5;frac='+str(frac)
                if verbose_results:
                    probs = [float(x) for x in pos_dict_verbose[locus]]
                    se_95 = 2*stats.sem(probs)
                    deets = deets + ';fracLow='+str(frac-se_95)+';fracUp='+str(frac+se_95)+';identificationQv='+str(int(100*np.mean([float(x) for x in pos_dict_verbose[locus]])))
                gff_info = (locus[0],locus[2],locus[4],deets)
                write_gff(outfi,gff_info)
            else:
                out_line = '\t'.join(list(locus)[:-1]+[str(np.mean(pos_dict[locus]))]+[locus[-1]]+[str(len(pos_dict[locus]))]) #+[str(x) for x in values_dict[locus]])
                if pos_list:
                    out_line = out_line + '\t' + '\t'.join([str(x) for x in values_dict[locus]])
                if verbose_results:
                    out_line = out_line + '\t'+','.join(pos_dict_verbose[locus])
                outfi.write(out_line+'\n')
    if not pos_list:
        if not control:
            print count, 'methylated loci found with min depth', depth_thresh, 'reads'
        else:
            print count, 'unmethylated loci found with min depth', depth_thresh, 'reads'

def main():
    #parse command line options
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Produce bed file of methylated positions based on mCaller output')
    parser.add_argument('-d','--min_read_depth',type=int,required=False,help='minimum coverage of position to determine methylation (default = 15)',default=15)
    parser.add_argument('-t','--mod_threshold',type=float,required=False,help='minimum %% of observations at a position to include in report (default = 0.5)',default=0.5)
    parser.add_argument('-f','--mCaller_file',type=str,required=True,help='the output file from mCaller to summarize')
    parser.add_argument('-p','--positions',type=str,required=False,help='~bed file of positions for which to calculate % methylated (chromosome,start,end,strand); ignores other thresholds')
    parser.add_argument('--control',action='store_true',required=False,help='take unmethylated positions as a control for motif detection')
    parser.add_argument('--gff',action='store_true',required=False,help='output PacBio-style gff instead of bed ("identificationQv" score will be average probability of methylation)')
    parser.add_argument('--ref',type=str,required=False,help='use reference fasta to output longer contexts surrounding a base, from -20 to +20')
    parser.add_argument('--plot',action='store_true',required=False,help='plot currents deviations at the positions included (not recommended for many positions)')
    parser.add_argument('--plotdir',type=str,required=False,default='mCaller_position_plots',help='output directory for plots, default=mCaller_position_plots')
    parser.add_argument('--vo',action='store_true',required=False,help='verbose output including probabilities for each position')
    parser.add_argument('-v','--version',action='store_true',required=False,help='print version')
    args = parser.parse_args()

    if args.version:
        print 'mCallerNP 0.1'
        sys.exit(0)

    assert os.path.isfile(args.mCaller_file), 'file not found at '+args.mCaller_file
    if args.positions:
        output_file = args.mCaller_file.split('.')[0]+'.methylation.positions.summary'
    elif not args.control:
        output_file = args.mCaller_file.split('.')[0]+'.methylation.summary'
    else:
        output_file = args.mCaller_file.split('.')[0]+'.methylation.control.summary'
    if args.gff:
        output_file = output_file+'.gff'
    else:
        output_file = output_file+'.bed'

    print args.mCaller_file

    aggregate_by_pos(args.mCaller_file,output_file,args.min_read_depth,args.mod_threshold,args.positions,args.control,args.vo,args.gff,args.refi,args.plot,args.plotdir)

if __name__ == "__main__":
    main()

