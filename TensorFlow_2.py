import uproot as ut
import numpy as np
import xgboost
import tensorflow as tf
from tensorflow import keras

from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc

import matplotlib.pyplot as plt

Bruit = ut.open('/users/LHCb/polcherrafael/Data/Data_Bruit.root')["t"]
Signal = ut.open('/users/LHCb/polcherrafael/MC/MC_BKGCAT10.root')["t"]
data  = ut.open('/users/LHCb/polcherrafael/Data/Select_sig.root')["t"]

# define the variables we want to use
vars = ["Lb_PT", "Lb_IPCHI2_OWNPV","Lb_FDCHI2_OWNPV", "Lb_LOKI_DTF_CHI2NDOF",
"Lb_ENDVERTEX_CHI2", "Jpsi_FDCHI2_OWNPV", "Lambdastar_PT", "Lambdastar_IPCHI2_OWNPV",
"Lambdastar_ENDVERTEX_CHI2", "Proton_P",]

tf.keras.backend.set_floatx('float64') #to set the data type

# create a pandas data frame with these variables only
sig_array = Signal.pandas.df(vars).to_numpy()
bkg_array = Bruit.pandas.df(vars).to_numpy()

print("Signal shape:", sig_array.shape)
print("Backgr shape:", bkg_array.shape)

data_array = data.pandas.df(vars).to_numpy()  #to select the signal in the data at the end

# merge and define signal and background labels
X = np.concatenate((sig_array, bkg_array))
y = np.concatenate((np.ones(sig_array.shape[0]),# 1 is signal
                    np.zeros(bkg_array.shape[0]))) # 0 is background

# split data in train and test samples
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.5)

print("Train size:", X_train.shape[0])
print("Test size: ", X_test.shape[0])


# define model
model = keras.Sequential([
    keras.layers.Dense(10, activation='relu'),   #n neurones
    keras.layers.Dense(10, activation='relu'),   #n neurones Second layer
    keras.layers.Dense(1)                        # 2 outcomes bruit/signal
])

model.compile(optimizer='adam',
          loss=tf.keras.losses.BinaryCrossentropy(from_logits=True),
          metrics=['accuracy'])


# train it
model.fit(X_train, y_train, epochs=20)
print("Model has been trained")
model.summary()

# make predictions
probability_model = tf.keras.Sequential([model, 
                                     tf.keras.layers.Softmax()])

y_predict_train = probability_model.predict(X_train)
y_predict_test = probability_model.predict(X_test)

y_prob_data = probability_model.predict(data_array)
print(y_prob_data)

