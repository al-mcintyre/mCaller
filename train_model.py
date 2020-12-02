import pickle 
import numpy as np
#from plotlib import *
from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn import svm
from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import GaussianNB

from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GroupKFold
from sklearn.model_selection import GroupShuffleSplit
from sklearn.ensemble import GradientBoostingClassifier

#make position to label dict for training
def pos2label(positions):
   pos2label_dict = {(pos.split()[0],int(pos.split()[1]),pos.split()[2]):pos.split()[3] for pos in open(positions,'r').read().split('\n') if len(pos.split()) > 1} #no longer -1 b/c switching everything to 0-based references
   return pos2label_dict

#to show best results of parameter grid search
def report(results, n_top):
   for i in range(1, n_top + 1):
       candidates = np.flatnonzero(results['rank_test_score'] == i)
       for candidate in candidates:
           print("Model with rank: {0}".format(i))
           print("Mean validation score: {0:.3f} (std: {1:.3f})".format(
                 results['mean_test_score'][candidate],
                 results['std_test_score'][candidate]))
           print("Parameters: {0}".format(results['params'][candidate]))

def train_classifier(signals,groups,modelfile,classifier='NN',plot=False): #TODO: set order of labels
    models = {}
    for twobase_model in signals:
        #print len(signals[twobase_model]['A']), len(signals[twobase_model]['m6A']) 
        #sys.exit(0)
        #print(twobase_model)
        if classifier == 'RF':
            model = RandomForestClassifier(bootstrap=True, class_weight=None, criterion='entropy',
                max_depth=10, max_features=4, max_leaf_nodes=None,
                min_impurity_split=3.93927246304e-06, min_samples_leaf=2,
                min_samples_split=3, min_weight_fraction_leaf=0.0,
                n_estimators=50, oob_score=False, random_state=None,
                verbose=0, warm_start=False)
            #param_grid = {'bootstrap':[True,False],'n_estimators':[40,50,60,100],'max_depth':[5,8,10,12],'max_features':[1,2,3,4],'min_samples_leaf':[1,2,3,10],'n_jobs':[15]}          
        elif classifier == 'NN':
            #param_grid = {'hidden_layer_sizes':[(4,4,4,4,4),(100),(100,100,100),(100,100,100,100)],'alpha':[0.001],'learning_rate':['adaptive'],'early_stopping':[False],'activation':['tanh']}
            model = MLPClassifier(hidden_layer_sizes=(100), alpha=0.001,learning_rate='adaptive',early_stopping=False,activation='tanh') #solver='lbfgs', alpha=1e-5, hidden_layer_sizes=(4), activation='tanh', random_state=1)
    
        elif classifier == 'SVM':
            #param_grid = {'kernel':['poly','rbf','sigmoid'],'degree':[3,5],'tol':[1e-3],'shrinking':[True,False],'gamma':['auto']} #,0.1,0.01,0.001]}
            model = svm.SVC(kernel='rbf',probability=True)

        elif classifier == 'LR':
            #param_grid = {'solver':['liblinear'],'multi_class':['ovr'],'penalty':['l1']}
            model = LogisticRegression(solver='liblinear',multi_class='ovr',penalty='l1')
 
        elif classifier == 'NBC':
            model = GaussianNB()

        if groups:
            gfk = GroupKFold(n_splits=5)
        else:
            gfk = 5
        
        """
        combinations = 1
        niter = 1
        for param in param_grid:
            combinations = combinations*len(param_grid[param]) 
        random_search = RandomizedSearchCV(model, param_distributions = param_grid, n_iter=combinations)
        As = signals[twobase_model]['A']
        m6As = signals[twobase_model]['m6A']
        random_search.fit(As+m6As,['A']*len(As)+['m6A']*len(m6As))
        report(random_search.cv_results_,combinations) #min(1,combinations))
        sys.exit(0)
        model=None        
        """

        num_examples = min([len(signals[twobase_model][label]) for label in signals[twobase_model]])
        labs, sigs, grps = [], [], []
        for label in signals[twobase_model]:
            labs = labs + [label]*num_examples
            sigs = sigs + signals[twobase_model][label][:num_examples]
            grps = grps + groups[twobase_model][label][:num_examples]

        print(labs[:10])
        print(sigs[:10])
        print(grps[:10])

        scores = cross_val_score(model,sigs,labs,cv=gfk,groups=grps)
        print("%s %s model scores: %s" %(classifier,twobase_model,','.join([str(s) for s in scores])))
        print("Cross validation accuracy: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))

        #scores = cross_val_score(model,[signal[:-1] for signal in sigs],labs,cv=gfk,groups=grps)
        #print('no quality',classifier, twobase_model, 'model scores:', scores)
        #print("Cross validation accuracy: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))

        model.fit(sigs,labs)
        models[twobase_model] = model

        if plot:
            model.fit(signals[twobase_model]['m6A'][:int(num_examples/2)]+signals[twobase_model]['A'][:int(num_examples/2)],['m6A']*(int(num_examples/2))+['A']*(int(num_examples/2)))
            prob_scores = {}
            prob_scores['m6A'] = [x[0] for x in model.predict_proba(signals[twobase_model]['m6A'][int(num_examples/2):])]
            prob_scores['A'] = [x[0] for x in model.predict_proba(signals[twobase_model]['A'][int(num_examples/2):])]
            plot_training_probabilities(prob_scores,twobase_model)      

    modfi = open(modelfile,'wb')
    pickle.dump(models,modfi)
    modfi.close()
    return models

