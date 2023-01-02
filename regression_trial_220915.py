#!/usr/bin/env python
import numpy as np
import pandas as pd
import os
from sklearn.multioutput import MultiOutputRegressor
from sklearn.model_selection import KFold,GridSearchCV
from sklearn.linear_model import Ridge, Lasso, ElasticNet
import xgboost as xgb
from xgboost import XGBRegressor
from sklearn.metrics import mean_squared_error as MSE
import matplotlib.pyplot as plt
import statistics
import statsmodels.api as sm
import itertools
import time
from contextlib import redirect_stdout
import argparse

#parser = argparse.ArgumentParser()
#parser.add_argument('--dir', type=str, default='./Phase_add7/', help='phase ann path')
#parser.add_argument('--ONLY_ARMES', type=bool, default=True)
#parser.add_argument('--neglect_junkframe', type=bool, default=True,help= 'neglect 1,2,21 phases')
#parser.add_argument('--neglect_33',type=bool, default=True,help= 'neglect 33')
#args, remaining_args = parser.parse_known_args()
#dir = args.dir

#import data ; np.matrix
file_dir = '/Users/jinhyesu/my_project/rawData/220511_ct/swd_result_220826/'
tmp = pd.read_excel(os.path.join(file_dir,'./Vvector_raw+predicted.xlsx'),sheet_name = 'NormalizeFeature(93)')
#tmp = tmp.sample(frac=0.8, random_state=1234)
print(tmp.describe())
kf = KFold(n_splits = 5, shuffle = True, random_state = 2)
#result = next(kf.split(tmp), None);print(result)

#a) all data
# does it be necessary to adjust with specific x variables?
arr_fold_pd = pd.DataFrame(); y_pred_fold_pd =pd.DataFrame()
ft_import_all = pd.DataFrame(tmp.iloc[:,0:30].columns, columns=['index'])
for train_index, test_index in kf.split(tmp):
    train, test = tmp.iloc[train_index,:], tmp.iloc[test_index, :]
    print(train, test)

    var_df, y_df = train.iloc[:, 0:30], train.iloc[:, 30:]
    var_mat, y_mat = np.matrix(var_df), np.matrix(y_df) #np.shape(y_mat)
    print(var_mat, y_mat)

    #xgboost regression
    # fitting
    myxgb = xgb.XGBRegressor(objective='reg:squarederror', n_estimators=100)
    multioutputregressor = MultiOutputRegressor(myxgb).fit(var_mat,y_mat)

    #predict ; RMSE
    test_var_mat, test_y_mat = np.matrix(test.iloc[:,0:30]), np.matrix(test.iloc[:,30:])
    arr= np.sqrt(np.mean(np.square(multioutputregressor.predict(test_var_mat)-test_y_mat), axis=0))
    #np.savetxt('predicted_220915.csv',arr,delimiter=",")
    print('arr is')
    print(arr)

    #average with arr
    arr_pd = pd.DataFrame(np.transpose(arr))  
    arr_fold_pd =  pd.concat([arr_fold_pd, arr_pd], axis=1)  

    #predict ; y
    t_var_df, t_y_df = test.iloc[:, 0:30], test.iloc[:, 30:]
    y_pred = pd.DataFrame(np.mean(multioutputregressor.predict(np.matrix(t_var_df)), axis=0)) #mean by x, y, z 75
    y_pred_fold_pd = pd.concat([y_pred_fold_pd, y_pred], axis=1)

    #importance feature
    est = multioutputregressor.estimators_[0]
    feature_tmp = pd.DataFrame(est.feature_importances_, columns=['importance'])
    feature_tmp['index'] = var_df.columns
    ft_import_all = pd.merge(ft_import_all, feature_tmp, on='index', how='outer')

#all_rmse_df
arr_fold_pd.index = np.array(test.columns[30:])
arr_fold_pd['average'] = arr_fold_pd.mean(axis=1)
arr_fold_pd.to_excel('/Users/jinhyesu/Desktop/rmse_all_221213.xlsx')

#all_difference_df ; difference / predicted y^ 
#y_pred_fold_pd.index = np.array(test.columns[30:])
#y_pred_fold_pd['average'] = y_pred_fold_pd.mean(axis=1)
#y_pred_fold_pd.to_excel('/Users/jinhyesu/Desktop/diff_221213.xlsx')

#plot ; importance feature
ft_import_all['average'] = ft_import_all.mean(axis = 1)
ft_import_all = ft_import_all.sort_values('average') 
ft_import_all = ft_import_all.rename(index=ft_import_all['index'])
ft_import_all.plot(kind='barh', y='average') 
#plt.show()
plt.savefig(os.path.join('/Users/jinhyesu/python_pj/ct/importance_feature_221008/all_feature_importances'), format='png', dpi=500)
plt.close()

#b) step-wise regression ; 
# return AIC
def processSubset(X, y, feature_set):
    if(len(feature_set) == 0):
        X = sm.add_constant(X)
        feature_set += ['const']
    model = sm.OLS(y, X[list(feature_set)])
    regr = model.fit()
    AIC = regr.aic
    return {'model':regr, 'AIC':AIC}

# step for forward selection 
def forward(X, y, predictors):
    remaining_predictors = [p for p in X.columns.difference(['const']) if p not in predictors]
    tic = time.time()
    results=[]
    for p in remaining_predictors:
        results.append(processSubset(X, y, feature_set=predictors+[p]))#+['const']
    models = pd.DataFrame(results)
    
    best_model = models.loc[models['AIC'].argmin()]
    toc = time.time()
    print("Processed ",models.shape[0], "models on", len(predictors)+1, "predictors in", (toc-tic))
    print("Selected predictors:",best_model["model"].model.exog_names,"AIC: ",best_model[0])    
    return best_model

# forward selection
def forward_model(X,y):
    Fmodels = pd.DataFrame(columns=["AIC", "model"])
    tic = time.time()
    predictors = []
    
    for i in range(1, len(X.columns.difference(['const']))+1):
        Forward_result = forward(X=X,y=y,predictors=predictors)
        if i > 1:
            if Forward_result['AIC'] > Fmodel_before: # 변수를 추가하면서 AIC가 증가하면 stop
                break
        Fmodels.loc[i] = Forward_result
        predictors = Fmodels.loc[i]["model"].model.exog_names
        Fmodel_before = Fmodels.loc[i]["AIC"]
        predictors = [ k for k in predictors if k != 'const']
        toc = time.time()
        #print("Total elapsed time : ", (toc-tic), "seconds.")
        print(predictors)
        print(Fmodels['AIC'])
    return(Fmodels['model'][len(Fmodels['model'])])
#Forward_best_model = forward_model(X=var_df, y=y_df['25z'])

def backward(X,y,predictors):
    results = []
   
    for combo in itertools.combinations(predictors, len(predictors) - 1):
        results.append(processSubset(X=X, y= y,feature_set=list(combo)))#+['const']
    models = pd.DataFrame(results)
   
    best_model = models.loc[models['AIC'].argmin()]
    print('Selected predictors:',best_model['model'].model.exog_names,' AIC:',best_model[0] )
    return best_model
#BK= backward(var_df, y_df['25z'], ['5Peri', 'Gender', '1Peri', '4AP', '5Ratio', '5AP', 'AGE', '3MA', '2MA', '2Ratio', '1MA'])

def Stepwise_model(X,y):
    Stepmodels = pd.DataFrame(columns=["AIC", "model"])
    predictors = []
    Smodel_before = processSubset(X=X, y=y, feature_set=predictors)['AIC']
    predictors = []
    for i in range(1, len(X.columns) + 1):
        print(i)
        X = sm.add_constant(X)

        Forward_result = forward(X=X, y=y, predictors=predictors) # constant added
        print('forward')
        print(predictors)
        Stepmodels.loc[i] = Forward_result
        predictors = Stepmodels.loc[i]["model"].model.exog_names
        predictors = [k for k in predictors]
        Backward_result = backward(X=X, y=y, predictors=predictors)  # Check if there is anything to remove
        if Backward_result['AIC']< Forward_result['AIC']:
            Stepmodels.loc[i] = Backward_result
            predictors = Stepmodels.loc[i]["model"].model.exog_names
            Smodel_before = Stepmodels.loc[i]["AIC"]
            predictors = [ k for k in predictors]
        if Stepmodels.loc[i]['AIC']> Smodel_before:
            break
        else:
            Smodel_before = Stepmodels.loc[i]["AIC"]
        
    print(Stepmodels)
    return (Stepmodels.loc[i]['model'])

def Regular_fuc(X,y, model_name): #L: lasso, R: ridge, E:elasticNet
    # define model
    if  model_name== 'L':
        model = Lasso()
    elif model_name == 'R':
        model = Ridge()
    elif model_name == 'E':
        model = ElasticNet()
    else:
        print('check model_name!')
    # define grid
    grid = dict()
    grid['alpha'] = 10**np.linspace(10,-2,100)*0.5
    # define search
    search = GridSearchCV(model, grid, scoring='neg_mean_absolute_error', cv=[(slice(None), slice(None))], n_jobs=-1, return_train_score=True)
    # perform the search
    #results = search.fit(X, y)
    # summarize
    #print('MAE: %.3f' % results.best_score_)
    #print('Config: %s' % results.best_params_)
    return(search)

my_predictions = {}; my_pred = None; my_actual = None; my_name = None

colors = ['r', 'c', 'm', 'y', 'k', 'khaki', 'teal', 'orchid', 'sandybrown',
          'greenyellow', 'dodgerblue', 'deepskyblue', 'rosybrown', 'firebrick',
          'deeppink', 'crimson', 'salmon', 'darkred', 'olivedrab', 'olive', 
          'forestgreen', 'royalblue', 'indigo', 'navy', 'mediumpurple', 'chocolate',
          'gold', 'darkorange', 'seagreen', 'turquoise', 'steelblue', 'slategray', 
          'peru', 'midnightblue', 'slateblue', 'dimgray', 'cadetblue', 'tomato'
         ]

def plot_predictions(_name, pred, actual):
    df = pd.DataFrame({'prediction': pred, 'actual': actual})
    df = df.sort_values(by='actual').reset_index(drop=True)

    plt.figure(figsize=(11, 8))
    plt.scatter(df.index, df['prediction'], marker='x', color='r')
    plt.scatter(df.index, df['actual'], alpha=0.7, marker='o', color='black')
    plt.title(_name, fontsize=15)
    plt.legend(['prediction', 'actual'], fontsize=12)
    plt.savefig(os.path.join('/Users/jinhyesu/python_pj/ct/regression_results/221214_test/{a}/{b}scatter_pred.png'.format(a=_name, b=y)), format='png', dpi=300)
    plt.close()

def plot_coef(_name, columns, coef):
    coef_df = pd.DataFrame(list(zip(columns, coef)))
    coef_df.columns=['feature', 'coef']
    coef_df = coef_df.sort_values('coef', ascending=False).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(9, 7))
    ax.barh(np.arange(len(coef_df)), coef_df['coef'])
    idx = np.arange(len(coef_df))
    ax.set_yticks(idx)
    ax.set_yticklabels(coef_df['feature'])
    fig.tight_layout()
    plt.savefig(os.path.join('/Users/jinhyesu/python_pj/ct/regression_results/221214_test/{a}/{b}_coef.png'.format(a = _name, b= y)), format='png', dpi=300)
    plt.close()

def vis_Regular_mod_fuc(_name, model_results, y):
    
    model_results.fit(var_df, y_df[y])
    pred_y2 = model_results.predict(t_var_df)
    plot_predictions(_name, pred_y2, test[y])
    plot_coef(_name, var_df.columns,model_results.best_estimator_.coef_)

    rmse_col.append(np.sqrt(MSE(t_y_df, pred_y2)))
        
    aic_col.append(" ") #AIC lasso_results.fit(var_df, y_df['1x']).best_params_ #{'alpha': 0.005}
    tmp_var_df = pd.DataFrame({'var_name' : t_var_df.columns, 'coef' :model_results.best_estimator_.coef_})
    tmp_var_df = tmp_var_df.loc[(abs(tmp_var_df['coef'])>0.1),:]
    selected_var.append(np.array(tmp_var_df.var_name))

def selcted_modeling(_name, model_results, y=y):
    with open(os.path.join('/Users/jinhyesu/python_pj/ct/regression_results/221214_test/feature_selection/', '{a}_{b}_summary.txt'.format(a=_name, b=y)), 'w') as f:
            with redirect_stdout(f):
                model_results.summary()
    res= model_results.resid
    fig=sm.qqplot(res, fit=True, line='45') #qqplot of residual
    plt.savefig(os.path.join('/Users/jinhyesu/python_pj/ct/regression_results/221214_test/feature_selection/', '{a}_{b}_qqplot.png'.format(a=_name, b=y)), format='png', dpi=300)
    plt.close()
    pred_y = model_results.predict(var_df)
    # residual pattern ; check homoscedasticity
    fig = plt.scatter(pred_y, res, s=4)
    plt.xlabel('Fitted values')
    plt.ylabel('Residual')
    plt.savefig(os.path.join('/Users/jinhyesu/python_pj/ct/regression_results/221214_test/', '{a}_{b}_full_residual.png'.format(a=_name, b=y)), format='png', dpi=300)
    plt.close()
    t_var_df, t_y_df = test[col_tmp], test[y]
    pred_y2 = model_results.predict(t_var_df)
    # residual plot
    plt.plot(np.array(t_y_df-pred_y2), label='pred_full_features')
    plt.legend()
    plt.savefig(os.path.join('/Users/jinhyesu/python_pj/ct/regression_results/221214_test/', '{a}_{b}_pred_full.png'.format(a=_name, b=y)), format='png', dpi=300)
    plt.close()
    #MSE
    mse = MSE(t_y_df, pred_y2)
    #RMSE
    rmse = np.sqrt(mse)
    rmse_col.append(rmse)
    print(rmse)
    #etc
    index.extend(y for _ in range(5))
    fold_num.extend(fold_cur for _ in range(5))
    mod_col.extend(['full', 'step', 'lasso', 'ridge', 'elasticNet'])
    print(fitted_full_model.aic)
    aic_col.append(fitted_full_model.aic)
    selected_var.append(np.array(fitted_full_model.params.index))


print(tmp.describe())
kf = KFold(n_splits = 5, shuffle = True, random_state = 2)
index = []; mod_col = []; fold_num = []; selected_var = []
aic_col = []; rmse_col = []; fold_cur = 1

for train_index, test_index in kf.split(tmp): #fold ; 5
    train, test = tmp.iloc[train_index,:], tmp.iloc[test_index, :]
    var_df_origin, y_df = train.iloc[:, 0:30], train.iloc[:, 30:]

    for y in y_df.columns.tolist():
        print(y)
        #select specific ct's volumetry features
        col_tmp= ['Weight', 'Height', "BMI", "AGE", "Gender"]
        if(y in {'1x', '2x', '3x', '4x', '5x', '1y', '2y', '3y', '4y', '5y', '1z', '2z', '3z', '4z', '5z'}):
            col_tmp.extend(var_df_origin.filter(regex='1').columns)
        elif(y in {'6x', '7x', '8x', '9x', '10x', '6y', '7y', '8y', '9y', '10y', '6z', '7z', '8z', '9z', '10z'}):
            col_tmp.extend(var_df_origin.filter(regex='2').columns)
        elif(y in {'11x', '12x', '13x', '14x', '15x', '11y', '12y', '13y', '14y', '15y', '11z', '12z', '13z', '14z', '15z'}):
            col_tmp.extend(var_df_origin.filter(regex='3').columns)
        elif(y in {'16x', '17x', '18x', '19x', '20x', '16y', '17y', '18y', '19y', '20y', '16z', '17z', '18z', '19z', '20z'}):
            col_tmp.extend(var_df_origin.filter(regex='4').columns)
        else:
            col_tmp.extend(var_df_origin.filter(regex='5').columns)

        var_df = var_df_origin[col_tmp]

        #full model
        full_model = sm.OLS(y_df[y], var_df)
        fitted_full_model = full_model.fit()
        with open(os.path.join('/Users/jinhyesu/python_pj/ct/regression_results/221214_test/', '{}_full_model_summary.txt'.format(y)), 'w') as f:
            with redirect_stdout(f):
                fitted_full_model.summary()
        print(fitted_full_model.summary())
        res= fitted_full_model.resid
        fig=sm.qqplot(res, fit=True, line='45') #qqplot of residual
        plt.savefig(os.path.join('/Users/jinhyesu/python_pj/ct/regression_results/221214_test/', '{}_qqplot.png'.format(y)), format='png', dpi=300)
        plt.close()
        pred_y = fitted_full_model.predict(var_df)
        # residual pattern ; check homoscedasticity
        fig = plt.scatter(pred_y, res, s=4)
        plt.xlabel('Fitted values')
        plt.ylabel('Residual')
        plt.savefig(os.path.join('/Users/jinhyesu/python_pj/ct/regression_results/221214_test/', '{}_full_residual.png'.format(y)), format='png', dpi=300)
        plt.close()
        t_var_df, t_y_df = test[col_tmp], test[y]
        pred_y2 = fitted_full_model.predict(t_var_df)
        # residual plot
        plt.plot(np.array(t_y_df-pred_y2), label='pred_full_features')
        plt.legend()
        plt.savefig(os.path.join('/Users/jinhyesu/python_pj/ct/regression_results/221214_test/', '{}_pred_full.png'.format(y)), format='png', dpi=300)
        plt.close()
        #MSE
        mse = MSE(t_y_df, pred_y2)
        #RMSE
        rmse = np.sqrt(mse)
        rmse_col.append(rmse)
        print(rmse)
        #etc
        index.extend(y for _ in range(5))
        fold_num.extend(fold_cur for _ in range(5))
        mod_col.extend(['full', 'step', 'lasso', 'ridge', 'elasticNet'])
        print(fitted_full_model.aic)
        aic_col.append(fitted_full_model.aic)
        selected_var.append(np.array(fitted_full_model.params.index))

        #stepwise_model
        Stepwise_best_model = Stepwise_model(X=var_df, y=y_df[y])
        with open(os.path.join('/Users/jinhyesu/python_pj/ct/regression_results/221214_test/', '{}_step_model_summary.txt'.format(y)), 'w') as f:
            with redirect_stdout(f):
                Stepwise_best_model.summary()
        print(Stepwise_best_model.summary())
        aic_col.append(Stepwise_best_model.aic)
        print(Stepwise_best_model.aic)
        sel_var = np.array(Stepwise_best_model.params.index)
        selected_var.append(sel_var)
        #MSE
        t_var_df = sm.add_constant(t_var_df)
        pred_y2 = Stepwise_best_model.predict(t_var_df[sel_var])
        mse = MSE(t_y_df, pred_y2)
        #RMSE
        rmse = np.sqrt(mse) #print('mse:', mse); print('rmse:', rmse)
        rmse_col.append(rmse)
        print(rmse)
        res= Stepwise_best_model.resid
        fig=sm.qqplot(res, fit=True, line='45') #qqplot of residual
        plt.savefig(os.path.join('/Users/jinhyesu/python_pj/ct/regression_results/221214_test/', '{}_step_qqplot.png'.format(y)), format='png', dpi=300)
        plt.close()
        var_df = sm.add_constant(var_df)
        pred_y = Stepwise_best_model.predict(var_df[sel_var])
        # residual pattern ; check homoscedasticity
        fig = plt.scatter(pred_y, res, s=4)
        plt.xlabel('Fitted values')
        plt.ylabel('Residual')
        plt.savefig(os.path.join('/Users/jinhyesu/python_pj/ct/regression_results/221214_test/', '{}_step_residual.png'.format(y)), format='png', dpi=300)
        plt.close()
        # residual plot
        plt.plot(np.array(t_y_df-pred_y2), label='pred_step_features')
        plt.legend()
        plt.savefig(os.path.join('/Users/jinhyesu/python_pj/ct/regression_results/221214_test/', '{}_pred_step.png'.format(y)), format='png', dpi=300)
        plt.close()


        #Regularized model; lasso, ridge, elastic-net
        t_var_df, t_y_df = test[col_tmp] #RMSE

        #lasso 
        lasso_results = Regular_fuc(var_df, y_df['2x'], 'L')
        vis_Regular_mod_fuc('lasso', lasso_results, y=y)

        #ridge
        ridge_results = Regular_fuc(var_df, y_df[y], 'R')
        vis_Regular_mod_fuc('ridge', ridge_results, y=y)

        #elasticNet
        elastic_results = Regular_fuc(var_df, y_df[y], 'E') #RMSE
        vis_Regular_mod_fuc('elasticNet', elastic_results, y=y)

    fold_cur = fold_cur +1
    print(fold_cur)

fin_df = pd.DataFrame({'index':index, 'mod_col':mod_col, 'fold_num':fold_num, 'sel_var':selected_var, 'aic_col':aic_col, 'rmse_col':rmse_col})  
fin_df.to_excel('/Users/jinhyesu/Desktop/rmse_all_221227.xlsx')

def vis_score():#visalization Score as follow alpha
    plt.figure()
    plt.title('Model')
    plt.xlabel('$\\alpha$ (alpha)')
    plt.ylabel('Score')
    # plot train scores
    alphas = 10**np.linspace(10,-2,100)*0.5
    plt.semilogx(alphas, lasso_results.cv_results_['mean_train_score'], label='Train', color='navy')
    plt.semilogx(alphas, lasso_results.cv_results_['mean_test_score'], label='Test', color='darkorange')
    plt.legend()
    plt.savefig(os.path.join('/Users/jinhyesu/python_pj/ct/regression_results/221214_test/', '{}_score_', str(lasso_results.best_estimator_)[0:5], '_.png'.format(y)), format='png', dpi=300)
    plt.close()
    plt.show()

def mse_eval(name_, pred, actual):
    global my_predictions, colors, my_pred, my_actual, my_name
    
    my_name = name_
    my_pred = pred
    my_actual = actual
    
    plot_predictions(name_, pred, actual)
    
    mse = MSE(pred, actual)
    my_predictions[name_] = mse
    
    y_value = sorted(my_predictions.items(), key=lambda x: x[1], reverse=True)
    
    df = pd.DataFrame(y_value, columns=['model', 'mse'])
    print(df)
    min_ = df['mse'].min() - 10
    max_ = df['mse'].max() + 10
    
    length = len(df) / 2
    plt.show()
