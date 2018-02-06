from extract_contexts import base_models

def tsv2matrix(tsvname,base):
    base_model = base_models(base)
    signals,contexts = {bm:{} for bm in base_model.values()},{bm:{} for bm in base_model.values()}
    with open(tsvname,'r') as infi:
        for line in infi:
            #ecoli1  c183b422-5dda-4a23-b732-309e8f7f331f    1794509 ATGCGMTCCAG     1.49,1.93166666667,-0.385,5.615,5.36,-0.945,15.7357504216     -       m6A
            context,sigs,strand,label = line.split('\t')[3:7]
            label = label.strip()
            twobase_model = base_model[context[int(len(context)/2):int(len(context)/2)+2]]
            if label not in signals[twobase_model]:
                signals[twobase_model][label] = []
                contexts[twobase_model][label] = []
            signals[twobase_model][label].append([float(s) for s in sigs.split(',')])
            contexts[twobase_model][label].append(context)
    return signals, contexts
