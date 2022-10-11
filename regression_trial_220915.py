#!/usr/bin/env python

import numpy as np
import pandas as pd
import os
import sklearn
import xgboost as xgb
from sklearn.multioutput import MultiOutputRegressor
from sklearn.model_selection import KFold
import matplotlib.pyplot as plt
from xgboost import XGBRegressor
from sklearn.metrics import mean_squared_error as MSE


#a)import data ; np.matrix
file_dir = '/Users/jinhyesu/my_project/rawData/220511_ct/swd_result_220826/'
tmp = pd.read_excel(os.path.join(file_dir,'./Vvector_raw+predicted.xlsx'),sheet_name = 'NormalizeFeature(93)')
#tmp = tmp.sample(frac=0.8, random_state=1234)
#print(tmp)

kf = KFold(n_splits = 5, shuffle = True, random_state = 2)
result = next(kf.split(tmp), None)
print(result)





train = tmp.iloc[result[0]]
test =  tmp.iloc[result[1]]

print(train, test)

var_df = train.iloc[:, 0:30]
y_df = train.iloc[:, 30:]

var_mat = np.matrix(var_df)
y_mat = np.matrix(y_df) #np.shape(y_mat)
print(var_mat, y_mat)

#b)xgboost regression
# fitting
myxgb = xgb.XGBRegressor(objective='reg:squarederror', n_estimators=100)
multioutputregressor = MultiOutputRegressor(myxgb).fit(var_mat,y_mat)

# predicting
#print(np.shape(var_mat), np.shape(y_mat))

test_var_mat = np.matrix(test.iloc[:,0:30])
test_y_mat = np.matrix(test.iloc[:,30:])
arr= np.sqrt(np.mean(np.square(multioutputregressor.predict(test_var_mat)-test_y_mat), axis=0))
np.savetxt('predicted_220915.csv',arr,delimiter=",")

np.square(multioutputregressor.predict(test_var_mat)-test_y_mat)

#importance feature
est = multioutputregressor.estimators_[0]
feature_importances = pd.DataFrame(est.feature_importances_, columns=['importance'], index=var_df.columns).sort_values('importance')
print(feature_importances)
feature_importances.plot(kind='barh')

plt.show()

#c) single regression ; target 18y
target_col = np.array(test.columns[30:])
arr_pd = pd.DataFrame(np.transpose(arr), columns=['all_rmse'], index = target_col)
arr_pd.to_excel('/Users/jinhyesu/Desktop/rmse_all_221011.xlsx')

rmse_df = pd.DataFrame()
for target in target_col:
    print(target)
    var_df = train.iloc[:, 0:30]
    y_df = train[[target]]

    var_mat = np.matrix(var_df)
    y_mat = np.matrix(y_df) #np.shape(y_mat)
    print(var_mat, y_mat)

    #b)xgboost regression
    # fitting
    myxgb = xgb.XGBRegressor(objective='reg:squarederror', n_estimators=100)
    single_regressor = myxgb.fit(var_mat,y_mat)

    # predicting
    #print(np.shape(var_mat), np.shape(y_mat))

    test_var_mat = np.matrix(test.iloc[:,0:30])
    test_y_mat = np.matrix(test[[target]])

    pred = single_regressor.predict(test_var_mat)
    rmse = np.sqrt(MSE(test_y_mat, pred))
    #arr= np.sqrt(np.mean(np.square(single_regressor.predict(test_var_mat)-test_y_mat), axis=0))

    #importance feature
    feature_importances = pd.DataFrame(myxgb.feature_importances_, columns=['importance'], index=var_df.columns).sort_values('importance')
    print(feature_importances)
    feature_importances.plot(kind='barh')
    #plt.show()
    plt.savefig(os.path.join('/Users/jinhyesu/python_pj/ct/importance_feature_221008/', '{}.png'.format(target)), format='png', dpi=500)
    plt.close()

    #generate df
    tmp_df = feature_importances
    tmp_df['target'] = np.repeat(target, len(feature_importances), axis=0)
    tmp_df['each_rmse'] = np.repeat(rmse, len(feature_importances), axis=0)
    tmp_df['all_rmse'] = np.repeat(arr_pd.loc[target][0], len(feature_importances), axis=0)

    rmse_df = pd.concat([rmse_df, tmp_df])

rmse_df.to_excel('/Users/jinhyesu/Desktop/rmse_df_221011.xlsx')
