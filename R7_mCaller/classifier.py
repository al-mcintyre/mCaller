import cPickle
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GroupKFold
from sklearn.model_selection import GroupShuffleSplit

def random_forest(signals,labels,train,groups=None):
    modfiname = 'random_forest_model.pkl'
    if train:
        #modfi = open(modfiname,'wb')
        model = RandomForestClassifier()
        gkf = GroupKFold(n_splits=5)
        scores = cross_val_score(model,signals,labels,groups=groups,cv=gkf)
        print scores
        print("Accuracy: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))
        #cPickle.dump(model,modfi)
        #modfi.close()
        return None
    else:
        modfi = open(modfiname,'rb')
        model = cPickle.load(modfi)
        probabilities = model.predict_proba(signals)
        modfi.close()
        return probabilities
        
        
