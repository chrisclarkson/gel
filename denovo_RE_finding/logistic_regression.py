import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy.special import expit
import pandas as pd
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.metrics import confusion_matrix

data=pd.read_csv('CC_AT/chr13_102160743_102162697_AAG_FGF14_stretch40_ranked_genomes.txt',sep='\t',header=None)

X=np.array(data[[4]])
y=np.array(data[[1]])
y=(y=='case').astype(float)
clf = LogisticRegression(C=1e5)
clf.fit(X, y)

plt.figure(1, figsize=(4, 3))
plt.clf()
plt.scatter(X.ravel(), y, label="example data", color="black", zorder=20)
X_test = X

def logit_pvalue(model, x):
    """ Calculate z-scores for scikit-learn LogisticRegression.
    parameters:
        model: fitted sklearn.linear_model.LogisticRegression with intercept and large C
        x:     matrix on which the model was fit
    This function uses asymtptics for maximum likelihood estimates.
    """
    p = model.predict_proba(x)
    n = len(p)
    m = len(model.coef_[0]) + 1
    coefs = np.concatenate([model.intercept_, model.coef_[0]])
    x_full = np.matrix(np.insert(np.array(x), 0, 1, axis = 1))
    ans = np.zeros((m, m))
    for i in range(n):
        ans = ans + np.dot(np.transpose(x_full[i, :]), x_full[i, :]) * p[i,1] * p[i, 0]
    vcov = np.linalg.inv(np.matrix(ans))
    se = np.sqrt(np.diag(vcov))
    t =  coefs/se  
    p = (1 - stats.norm.cdf(abs(t))) * 2
    return p

print(logit_pvalue(clf,X_test))
print(clf.coef_)
tn, fp, fn, tp=confusion_matrix(y, clf.predict(X_test)).ravel()
specificity = tn / (tn+fp)
print(specificity*tp)
sensitivity = tp / (tp+fn)
print(sensitivity)

loss = expit(X_test * clf.coef_ + clf.intercept_)
plt.plot(X_test, loss, label="Logistic Regression Model", color="red", linewidth=3)

plt.ylabel("y")
plt.xlabel("X")
plt.yticks([0, 0.5, 1])
plt.ylim(-0.25, 1.25)
plt.tight_layout()
plt.savefig('out.png')



