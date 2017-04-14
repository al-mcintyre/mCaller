import cPickle
import numpy as np
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

def report(results, n_top):
    for i in range(1, n_top + 1):
        candidates = np.flatnonzero(results['rank_test_score'] == i)
        for candidate in candidates:
            print("Model with rank: {0}".format(i))
            print("Mean validation score: {0:.3f} (std: {1:.3f})".format(
                  results['mean_test_score'][candidate],
                  results['std_test_score'][candidate]))
            print("Parameters: {0}".format(results['params'][candidate]))
            print("")

def train_classifier(signals,labels,groups,modelfile,classifier='NN'): #TODO: set order of labels
   classifier = classifier.split('_')[0]
   if classifier == 'RF':
      model = RandomForestClassifier(bootstrap=True, class_weight=None, criterion='entropy',
          max_depth=10, max_features=4, max_leaf_nodes=None,
          min_impurity_split=3.93927246304e-06, min_samples_leaf=2,
          min_samples_split=3, min_weight_fraction_leaf=0.0,
          n_estimators=50, n_jobs=15, oob_score=False, random_state=None,
          verbose=0, warm_start=False)
          #param_grid = {'bootstrap':[True,False],'n_estimators':[40,50,60,100],'max_depth':[5,8,10,12],'max_features':[1,2,3,4],'min_samples_leaf':[1,2,3,10],'n_jobs':[15]}          
   elif classifier == 'NN':
      #param_grid = {'hidden_layer_sizes':[(4,4,4,4.4),(100),(100,100,100),(100,100,100,100)],'alpha':[0.001],'learning_rate':['adaptive'],'early_stopping':[False],'activation':['tanh']}
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
        
   #combinations = 1
   #niter = 1
   #for param in param_grid:
   #   combinations = combinations*len(param_grid[param]) 
   #random_search = RandomizedSearchCV(model, param_distributions = param_grid, n_iter=min(1,combinations))
   #random_search.fit(signals,labels)
   #report(random_search.cv_results_,min(1,combinations))
   #model=None        

   model.fit(signals,labels)
   scores = cross_val_score(model,signals,labels,cv=gfk,groups=groups)
      
   print scores
   print("Accuracy: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))
  
   modfi = open(modelfile,'wb')
   cPickle.dump(model,modfi)
   modfi.close()
   return model

