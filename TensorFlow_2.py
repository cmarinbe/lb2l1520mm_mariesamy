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
sig_df = Signal.pandas.df(vars)
bkg_df = Bruit.pandas.df(vars)

data_df = data.pandas.df(vars)  #to select the signal in the data at the end

# add a target column and merge the two df
sig_df['target'] = np.ones(sig_df.shape[0])
bkg_df['target'] = np.zeros(bkg_df.shape[0])

print(sig_df.keys())
print(np.shape(sig_df['target']),np.shape(sig_df['Lb_PT']))

df = sig_df.append(bkg_df)

# split data in train and test samples
train_df, test_df = train_test_split(df, test_size=0.5)

# create datasets to train the model
target_train = train_df.pop('target').to_numpy()
target_train = np.reshape(target_train, (-1, 1))   #to have the right label shape (everything in one column)

target_test = test_df.pop('target').to_numpy()
target_test = np.reshape(target_test, (-1, 1))


print(train_df.keys()) # target is no longer there 


dataset_train = tf.data.Dataset.from_tensor_slices((train_df.to_numpy(), target_train))
dataset_test = tf.data.Dataset.from_tensor_slices((test_df.to_numpy(), target_test))


for feat, targ in dataset_train.take(5):
    print ('Features: {}, Target: {}'.format(feat, targ))


model = keras.Sequential([
    keras.layers.Dense(10, activation='relu'),   #n neurones
    keras.layers.Dense(10, activation='relu'),   #n neurones Second layer
    keras.layers.Dense(1)                        # 2 outcomes bruit/signal
])

model.compile(optimizer='adam',
          loss=tf.keras.losses.BinaryCrossentropy(from_logits=True),
          metrics=['accuracy'])


# train it
model.fit(dataset_train, epochs=20)
print("Model has been trained")
model.summary()

# make predictions
probability_model = tf.keras.Sequential([model, 
                                     tf.keras.layers.Softmax()])

y_predict_train = probability_model.predict(dataset_train)
y_predict_test = probability_model.predict(dataset_test)

y_prob_data = probability_model.predict(data_array)
print(y_prob_data)

