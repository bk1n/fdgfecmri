import data as d
import pandas as pd
import numpy as np

from sklearn.model_selection import cross_validate, train_test_split
from sklearn.model_selection import RepeatedKFold
from sklearn import impute
from sklearn.linear_model import LogisticRegression
from sklearn import metrics
from sklearn import svm
from sklearn.preprocessing import MinMaxScaler
from sklearn.pipeline import Pipeline


import matplotlib.pyplot as plt
 
class Model():
    def __init__(self) -> None:
        self.data = d.Data()
        
        #metrics
        f1_score = metrics.make_scorer(metrics.f1_score)
        accuracy = metrics.make_scorer(metrics.accuracy_score)
        sensitivity = metrics.make_scorer(metrics.precision_score)
        specificity = metrics.make_scorer(metrics.recall_score)
        auc = metrics.make_scorer(metrics.roc_auc_score)
        self.scorers = {'f1': f1_score, 
                        'acc': accuracy, 
                        'sens': sensitivity,
                        'spec': specificity,
                        'roc_auc': auc}
        
    def get_X_y(self, df):
        X = df.drop(['CANC', 'HIST', 'region'], axis=1)
        y = df.loc[:,'HIST']
        return X, y
    
    def knn_impute(self, X, k):
        knn = impute.KNNImputer(n_neighbors=k)
        colnames = X.columns
        X = pd.DataFrame(knn.fit_transform(X))
        X = X.set_axis(colnames, axis=1, inplace=False)
        return X

    def cross_val_scores(self, model, X, y):
        cv = RepeatedKFold(n_splits=4, n_repeats=30, random_state=1)
        return cross_validate(model, X, y, cv=cv, scoring=self.scorers)
        
    def print_scores(self, cross_val_scores):
        for k,v in cross_val_scores.items():
            if 'test' in k:
                print(f'{k}:, {v.mean():.3f}')
        print('\n')
    
if __name__ == '__main__':
    m = Model()
    m.data.remove_null_rows_without_fdg_adc()
    m.data.remove_columns_less50_null()
    #m.data.remove_calculated_variables()
    m.data.get_endo()
    m.data.remove_null_rows()
    
    m.data.stats()
    
    X,y = m.get_X_y(m.data.df)
    #X = X[['FDG_SUV', 'ADC', 'FDG_LA', 'FDG_SA']]
    X = X[['FDG_SUV_PT', 'ADC_PT', 'FDG_SUV', 'ADC', 'FDG_LA', 'FDG_SA']]
    
    # multivariate LASSO on RAW data
    range_ = np.linspace(0.1, 15, 40)
    f1 = []
    for l in range_:
        print(l)
        pipe = Pipeline([("scale", MinMaxScaler((0,1))),
                    ("logistic", LogisticRegression(random_state=5, solver='liblinear', penalty='l1', max_iter=4000, C=l))
            ])
        scores = m.cross_val_scores(pipe, X, y)
        f1.append(scores['test_f1'])
        
    f1_ = [a.mean() for a in f1]
    plt.plot(range_, f1_)

    lambda_ = range_[np.argmax(f1_)] 
    lambda_ = 8.123076923076923
    
    scale = MinMaxScaler((0,1))
    X = pd.DataFrame(scale.fit_transform(X), columns=X.columns)
    
    lr = LogisticRegression(random_state=5, solver='liblinear', penalty='l1', max_iter=4000, C=lambda_)
    lr.fit(X, y)
    print(pd.DataFrame(lr.coef_, columns=X.columns))
    
    pipe = Pipeline([("scale", MinMaxScaler((0,1))),
            ("logistic", LogisticRegression(random_state=5, solver='liblinear', penalty='l1', max_iter=4000, C=lambda_))
    ])
    scores = m.cross_val_scores(pipe, X, y)
    m.print_scores(scores)
    
    # univariate Logistic Regression w no regularisation
    f1, accuracy, sensitivity, specificity, auc = [], [], [], [], []
    cut_range = np.linspace(0,1,101)
    for cutoff in cut_range:
        X_train, X_test, y_train, y_test = train_test_split(X,y, train_size=0.7, shuffle=True)
        lr = LogisticRegression(random_state=0, solver='saga',penalty='none', max_iter=4000)
        lr.fit(X_train, y_train)
        proba_preds = lr.predict_proba(X_test)
        preds = (proba_preds[:,1]  > cutoff).astype('int') #Pr(M)
        
        f1.append(metrics.f1_score(y_test, preds))
        accuracy.append(metrics.accuracy_score(y_test, preds))
        sensitivity.append(metrics.precision_score(y_test, preds))
        specificity.append(metrics.recall_score(y_test, preds))
        auc.append(metrics.roc_auc_score(y_test, preds))