import scipy.io
import numpy as np
import matplotlib.pyplot as plt
from sklearn import svm
import lightgbm as lgb
from sklearn.metrics import accuracy_score, classification_report
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import balanced_accuracy_score


n_ctl_train = 28
n_sui_train = 20
n_ctl_test = 7
n_sui_test = 6

y_train = np.concatenate((np.zeros(n_ctl_train), np.ones(n_sui_train)))
y_test = np.concatenate((np.zeros(n_ctl_test), np.ones(n_sui_test)))

rst_folder = '/hd2/research/EEG/code/output/'

# model_name = "lgbm"
model_name = "svm"

is_PCA = True

acc_list = []

for cur_rank in range(2, 21): 

    print("current rank:", cur_rank)
    X_trains = []
    X_tests = []

    for task in ['happy', 'neutral', 'sad', 'suicide']:
        subj_coef = scipy.io.loadmat(rst_folder + 'subj_coef_' + task + '_R' + str(cur_rank) + '_ctl28_sui20_dep0_subNone.mat')
        X_train = subj_coef['subj_coef_train']
        X_test = subj_coef['subj_coef_test']

        scaler = StandardScaler()
        X_train = scaler.fit_transform(X_train)
        X_test = scaler.transform(X_test)

        X_trains.append(X_train)
        X_tests.append(X_test)

    X_train = np.concatenate(X_trains, axis=1)
    X_test = np.concatenate(X_tests, axis=1)

    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)

    # if is_PCA:
    #     n_comp_train = math.min(X_train.shape[0], X_train.shape[1])
    #     n_comp_test = math.min(X_test.shape[0], X_test.shape[1])
    #     pca = PCA(n_components=n_comp_train, svd_solver='full')
    #     X_train = pca.fit_transform(X_train)
    #     X_test = pca.

    if model_name == "svm":
        model = svm.SVC(C=1.5, kernel='rbf', gamma='auto', class_weight='balanced')
    elif model_name == "lgbm":
        params = {
            'objective': 'binary',
            'boosting_type': 'gbdt',
            'metric': 'logloss',
            'num_leaves': 16,
            'learning_rate': 0.05,
            'feature_fraction': 0.7,
            'bagging_fraction': 0.8,
            'bagging_freq': 5,
            'verbose': -1,
            'max_depth': 10,
            'min_data_in_leaf': 20,
        }

        model = lgb.LGBMClassifier(**params)

    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)

    accuracy = balanced_accuracy_score(y_test, y_pred)
    # accuracy = accuracy_score(y_test, y_pred)
    report = classification_report(y_test, y_pred)

    print("Accuracy:", accuracy)
    # print("Classification Report:\n", report)

    acc_list.append(accuracy)


# from matplotlib.ticker import MaxNLocator
import math

ax = plt.figure().gca()
# ax.yaxis.set_major_locator(MaxNLocator(integer=True))

plt.plot(range(math.floor(2), math.ceil(21)), acc_list)
plt.xticks(range(math.floor(2), math.ceil(21)))
plt.xlabel("Rank")
plt.ylabel("Accuracy")
plt.savefig("acc_vs_rank.png")
