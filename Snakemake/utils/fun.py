from sklearn.metrics import roc_curve
from sklearn.metrics import precision_recall_curve
import sklearn.metrics as sm
from matplotlib import pyplot

def evaluate(lab,prob):
    temp = pd.DataFrame()
    temp[1]=prob
    temp[0]=1-prob
    pred = temp.idxmax(axis = 1)
    precision, recall, thresholds = sm.precision_recall_curve(lab, prob)

    p=precision.mean()
    r=recall.mean()
    a=sm.accuracy_score(lab,pred)
    f=sm.f1_score(lab, pred, average='binary')
    au=sm.roc_auc_score(lab,prob)
    l = list(zip(['presision','recall','accuracy','f score','roc'],[p,r,a,f,au]))
    print(l)

    pyplot.plot([0, 1], [0, 1], linestyle='--', label='No Skill')
    fpr, tpr, _ = roc_curve(lab, prob)
    pyplot.plot(fpr, tpr, marker='.', label='parSMURF',color='green')
    pyplot.xlabel('False Positive Rate')
    pyplot.ylabel('True Positive Rate')
    pyplot.legend()
    pyplot.show()
    #marker='.'
    y = lab#[0]
    no_skill = len(y[y==1]) / len(y)
    pyplot.plot([0, 1], [no_skill, no_skill], linestyle='--', label='No Skill')
    precision, recall, _ = precision_recall_curve(y, prob)
    pyplot.plot(recall, precision, label='parSMURF',color='green')
    pyplot.xlabel('Recall')
    pyplot.ylabel('Precision')
    pyplot.legend()
    pyplot.show()
    


