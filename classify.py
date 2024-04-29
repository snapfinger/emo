import scipy.io
import numpy as np
import matplotlib.pyplot as plt
from sklearn import svm
from sklearn.metrics import accuracy_score, classification_report
from sklearn.preprocessing import StandardScaler


n_ctl_train = 28
n_sui_train = 20
n_ctl_test = 7
n_sui_test = 6

y_train = np.concatenate((np.zeros(n_ctl_train), np.ones(n_sui_train)))
y_test = np.concatenate((np.zeros(n_ctl_test), np.ones(n_sui_test)))

rst_folder = '/hd2/research/EEG/code/output/'

acc_list = []

for cur_rank in range(2, 11): 

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

    # Feature scaling for better performance of SVM
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)

    model = svm.SVC(kernel='rbf', gamma='scale', class_weight='balanced')
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)

    accuracy = accuracy_score(y_test, y_pred)
    report = classification_report(y_test, y_pred)

    print("Accuracy:", accuracy)
    # print("Classification Report:\n", report)

    acc_list.append(accuracy)


plt.plot(range(2, 11), acc_list)
plt.savefig("acc_vs_rank.png")
