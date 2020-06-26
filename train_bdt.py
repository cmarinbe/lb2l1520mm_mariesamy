
import uproot as ut
import numpy as np
import xgboost

from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc

import matplotlib.pyplot as plt

def plot_output(y_train, y_test, y_prob_train, y_prob_test,
                figname='overtraining.png'):
    """Plot train and test BDT output for signal and bkg.

    :param y_train: train labels
    :type y_train: array
    :param y_test: test labels
    :type y_test: array
    :param y_prob_train: train predictions
    :type y_prob_train: array
    :param y_prob_test: test predictions
    :type y_prob_test: array
    :param figname: title of output figure
    :type figname: str

    """
    fig = plt.figure()
    nbins = 50
    range = (0,1)

    # train: filled histo
    plt.hist(y_prob_train[y_train==1][:,1], bins=nbins,
             range=range,
             label="Sig Train", normed=True, #log=True,
             bottom=0,
             color="red", histtype="stepfilled", alpha=0.5)
    plt.hist(y_prob_train[y_train==0][:,1], bins=nbins,
             range=range,
             label="Bkg Train", normed=True, #log=True,
             bottom=0,
             color="blue", histtype="stepfilled", alpha=0.5)

    # test: error bars
    test_sig = y_prob_test[y_test==1][:,1]
    hist_sig, bins_sig = np.histogram(test_sig,
                                      range=range, bins=nbins,
                                      normed=True)
    cent_sig = (bins_sig[1:] + bins_sig[:-1])/2
    plt.errorbar(cent_sig, hist_sig, fmt='o',
                 color="red", label="Sig Test")

    test_bkg = y_prob_test[y_test==0][:,1]
    hist_bkg, bins_bkg = np.histogram(test_bkg,
                                      range=range, bins=nbins,
                                      normed=True, weights=test_bkg_w)
    cent_bkg = (bins_bkg[1:] + bins_bkg[:-1])/2
    plt.errorbar(cent_bkg, hist_bkg, fmt='o',
                 color="blue", label="Bkg Test")

    # plot settings
    plt.ylim(ymin=0)
    plt.xlabel("BDT output")
    plt.ylabel("Arbitrary units")
    plt.legend(loc="upper center")
    plt.savefig(figname)
    return fig


if __name__ == '__main__':
    # read the data
    sig = ut.open('/sps/lhcb/marin/RpK-fullr1r2/2017/MC/15114000/LeptonU_MagDown.root')["tuple_mmLine;1"]['DecayTree;1']
    bkg = ut.open('/sps/lhcb/marin/RpK-fullr1r2/2017/data/LeptonU_MagDown.root')["tuple_mmLine;1"]['DecayTree;1']

    # define the variables we want to use
    vars = ["Proton_IPCHI2_OWNPV", "Kaon_PT", "Lb_FD_OWNPV", "Jpsi_PT"]
    # create a pandas data frame with these variables only
    sig_array = sig.pandas.df(vars, entrystart=0, entrystop=1000).to_numpy() # use only 1000 events for this example!
    bkg_array = bkg.pandas.df(vars, entrystart=0, entrystop=1000).to_numpy() # use only 1000 events for this example!
    print("Signal shape:", sig_array.shape)
    print("Backgr shape:", bkg_array.shape)

    # merge and define signal and background labels
    X = np.concatenate((sig_array, bkg_array))
    y = np.concatenate((np.ones((sig_array.shape[0], 1)),# 1 is signal
                        np.zeros((bkg_array.shape[0], 1)))) # 0 is background

    # split data in train and test samples
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.5)

    # define a model, example: xgboost
    params_dict = {
        "n_estimators": 200,
        "learning_rate": 0.3,
        "max_depth": 3,
        }
    model = xgboost.XGBClassifier(**params_dict)

    # train it
    model.fit(X_train, y_train)

    # plot features importance
    xgboost.plot_importance(model)
    plt.savefig("vars_ranking.png")

    # make predictions
    y_prob_train = model.predict_proba(X_train)
    y_prob_test  = model.predict_proba(X_test)

    # plot train and test output
    overtr_file = os.path.join(args.out_dir, "overtraining.png")
    fig1 = plot_output(y_train, y_test, y_prob_train,
                       y_prob_test, figname=overtr_file,
                       train_w=w_train, test_w=w_test)

    # evaluate performance with area under the roc
    fpr, tpr, thresholds = roc_curve(y_train, y_prob_train[:,1])
    fnr = 1-fpr
    roc_auc = auc(tpr, fnr)
    fpr_test, tpr_test, thresholds_test = roc_curve(y_test, y_prob_test[:,1])
    fnr_test = 1-fpr_test
    roc_auc_test = auc(tpr_test, fnr_test)
    print("## Performance: ROC area")
    print("- train: {}".format(roc_auc))
    print("- test : {}".format(roc_auc_test))


# EOF
