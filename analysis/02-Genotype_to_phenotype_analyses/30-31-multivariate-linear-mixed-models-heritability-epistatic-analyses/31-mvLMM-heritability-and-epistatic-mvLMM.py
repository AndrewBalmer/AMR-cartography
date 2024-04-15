

### Setup - this script requires many different functions from versions of Limix and other packages
import limix
import pandas
import pandas as pd
import numpy
from numpy.random import RandomState
import matplotlib.pyplot as plt
from limix.stats import linear_kinship
from limix.her import estimate
from numpy import ones, dot, stack, eye, exp, concatenate, zeros, sqrt
from limix.vardec import VarDec
from pandas import DataFrame
from limix.qtl import scan
import math

print(limix.__version__) # This should be 3.0.4
print(pandas.__version__) # This should be 1.0.5

# Read in datasets
dummy_gen = limix.io.csv.read('S.pneumo_map_dummy_gen_relatedness_matrix.csv', verbose=False)
map_coords = limix.io.csv.read('S.pneumo_map_mvlmm_map_coords.csv', verbose=False)
mod_matrix = limix.io.csv.read('S.pneumo_map_mvlmm_mod_matrix.csv', verbose=False)
genetic_markers = limix.io.csv.read('S.pneumo_map_dummy_gen_test_markers.csv', verbose=False)
test_markers_epi = pandas.read_csv('S.pneumo_map_test_markers_incl_2nd_order_epistatic.csv')

map_coords_obs = map_coords
random = RandomState(2)


plt.scatter(map_coords['D1'], map_coords['D2'])
plt.xlabel('x-axis label')
plt.ylabel('y-axis label')
plt.title('Scatter Plot')
plt.show()

K = dummy_gen
K = K.to_numpy()
K = linear_kinship(K)
p = 2

##### Heritability estimates for map dimensions 1 and 2

V1_heritability_normal = estimate(map_coords['D1'], "normal", K, verbose=True)
V2_heritability_normal = estimate(map_coords['D2'], "normal", K, verbose=True)
print(V1_heritability_normal)
print(V2_heritability_normal)



## Per PBP variance decomposition
PBP1A = limix.io.csv.read('S.pneumo_map_PBP1A_dummy_gen_relatedness_matrix.csv', verbose=False)
PBP2B = limix.io.csv.read('S.pneumo_map_PBP2B_dummy_gen_relatedness_matrix.csv', verbose=False)
PBP2X = limix.io.csv.read('S.pneumo_map_PBP2X_dummy_gen_relatedness_matrix.csv', verbose=False)

#PBP1A = PBP1A.to_numpy()
K0 = dot(PBP1A, PBP1A.T)

#PBP2B = PBP2B.to_numpy()
K1 = dot(PBP2B, PBP2B.T)

#PBP2X = PBP2X.to_numpy()
K2 = dot(PBP2X, PBP2X.T)

### map dimension 1
Y = map_coords['D1']
x = VarDec(Y, "normal")

x.append(K0, "PBP1A")
x.append(K1, "PBP2B")
x.append(K2, "PBP2X")
x.append_iid("noise")
x.fit(verbose=True)
print(x)
x
x.plot()

variances = [c.scale for c in x._covariance]
variances = [(v / sum(variances)) * 100 for v in variances]
sum(variances)
variances # output the variances for 1a, 2b 2x and noise

### map dimension 2
Y = map_coords['D2']

map_axis_2 = VarDec(Y, "normal")

map_axis_2.append(K0, "PBP1A_Y")
map_axis_2.append(K1, "PBP2B_Y")
map_axis_2.append(K2, "PBP2X_Y")

map_axis_2.append_iid("noise")
map_axis_2.fit(verbose=True)
print(map_axis_2)
print(x)
map_axis_2
map_axis_2.plot()
variances_V2 = [c.scale for c in map_axis_2._covariance]
variances_V2 = [(v / sum(variances_V2)) * 100 for v in variances_V2]
sum(variances_V2)
variances_V2
variances # output the variances for 1a, 2b 2x and noise









### mvLMM for S pneumo
A0 = ones((p, 1))
A1 = eye(p)

print(A0)
print(A1)

mod_matrix = mod_matrix.to_numpy()

## run mvLMM without per test-variance decomposition
test = scan(G=genetic_markers, Y=map_coords, K=K, A=map_coords.cov(), verbose=False)

### run with per test variance decomposition
#a_list = list(range(0, genetic_markers.shape[1]))
a_list = list(range(0, 5))

pv_mv = ["lml0", "lml2","dof20","scale2","pv20"]
pv_mv = pandas.DataFrame(columns = pv_mv)

eff_mv_h2 = ["test", "trait", "effect_type","effect_name","env","effsize","effsize_se"]
eff_mv_h2 = pandas.DataFrame(columns = eff_mv_h2)

map_coords.cov()

for i in a_list:
    test_marker = genetic_markers.iloc[:, i].to_frame()
    K = dummy_gen
    K = K.drop(test_marker.columns[[0]], axis=1)
    K = K.to_numpy()
    K = linear_kinship(K)
    result = scan(G=test_marker, Y=map_coords, K=K, A = map_coords.cov(), verbose = False)
    pv_mv = pv_mv.append(result.stats)
    #eff_mv_h1 = eff_mv_h1.append(result.effsizes['h1'])
    eff_mv_h2 = eff_mv_h2.append(result.effsizes['h2'])

pv_mv.to_csv('mvLMM_p_values_normal_pneumo_low_freq_vars.csv')
eff_mv_h2.to_csv('mvLMM_effect_sizes_normal_pneumo_low_freq_vars.csv')






#### Random permutation GWAS
a_list = list(range(0, genetic_markers.shape[1])) # use this line to test all markers (will take a long time for the permutations)
#a_list = list(range(0, 5)) # - just test the first five markers

map_coords_full = ["D1", "D2"]
map_coords_full = pandas.DataFrame(columns = map_coords_full)

pv_mv_full = ["lml0", "lml2","dof20","scale2","pv20","repeat_p_index"]
pv_mv_full = pandas.DataFrame(columns = pv_mv_full)

eff_mv_full = ["test", "trait", "effect_type","effect_name","env","effsize","effsize_se","repeat_eff_index"]
eff_mv_full = pandas.DataFrame(columns = eff_mv_full)

# number of repeats
b_list = list(range(0, 100)) # do the full set of 100 permutations
#b_list = list(range(0, 2)) # here I just use two permutations as an example

for f in b_list:
    pv_mv = ["lml0", "lml2", "dof20", "scale2", "pv20"]
    pv_mv = pandas.DataFrame(columns=pv_mv)

    eff_mv_h2 = ["test", "trait", "effect_type", "effect_name", "env", "effsize", "effsize_se"]
    eff_mv_h2 = pandas.DataFrame(columns=eff_mv_h2)

    map_coords = map_coords_obs
    map_coords = map_coords.sample(frac=1).reset_index(drop=True)

    for i in a_list:
        test_marker = genetic_markers.iloc[:, i].to_frame()
        K = genetic_markers
        K = K.drop(test_marker.columns[[0]], axis=1)
        K = K.to_numpy()
        K = linear_kinship(K)
        result = scan(G=test_marker, Y=map_coords, K=K, A=map_coords.cov(), verbose = False)
        pv_mv = pv_mv.append(result.stats)
        # eff_mv_h1 = eff_mv_h1.append(result.effsizes['h1'])
        eff_mv_h2 = eff_mv_h2.append(result.effsizes['h2'])

    map_coords['repeat_p_index'] = pandas.DataFrame(numpy.full((map_coords.shape[0]), f + 1))
    pv_mv['repeat_p_index'] = pandas.DataFrame(numpy.full((pv_mv.shape[0]), f + 1))
    eff_mv_h2['repeat_eff_index'] = pandas.DataFrame(numpy.full((eff_mv_h2.shape[0]), f + 1))

    map_coords_full = map_coords_full.append(map_coords)
    pv_mv_full = pv_mv_full.append(pv_mv)
    eff_mv_full = eff_mv_full.append(eff_mv_h2)

minValuesObj = pv_mv.pv20.min()
pv_mv.pv20.min()

map_coords_full.to_csv('mvLMM_map_coords_normal_pneumo_random_phenotype_FWAS.csv')
pv_mv_full.to_csv('mvLMM_p_values_normal_pneumo_random_phenotype_FWAS.csv')
eff_mv_full.to_csv('mvLMM_effect_sizes_normal_pneumo_random_phenotype_FWAS.csv')





### testing for epistatic interactions

## for loop for running all epistatic interaction tests
#a_list = list(range(0, test_markers_epi.shape[1])) # test all epistatic markers
a_list = list(range(0, 3)) # just the first 3 to test


pv_mv_epi = ["lml0", "lml2","dof20","scale2","pv20"]
pv_mv_epi = pandas.DataFrame(columns = pv_mv_epi)

eff_mv_h2_epi = ["test", "trait", "effect_type","effect_name","env","effsize","effsize_se"]
eff_mv_h2_epi = pandas.DataFrame(columns = eff_mv_h2_epi)

for i in a_list:
    test_epistatic_marker = test_markers_epi.iloc[:, i].to_frame()
    colnames = test_epistatic_marker.columns[0]
    colnames = colnames.split(':')
    M = dummy_gen[colnames]
    dropped = dummy_gen.copy()
    dropped.drop(colnames, axis=1, inplace=True)
    dropped = dropped.to_numpy()
    K = linear_kinship(dropped)
    r = scan(G = test_epistatic_marker, Y = map_coords, K=K, M=M, A=mod_matrix, verbose=False)
    pv_mv_epi = pv_mv_epi.append(r.stats)
    eff_mv_h2_epi = eff_mv_h2_epi.append(r.effsizes['h2'])

pv_mv_epi.to_csv('mvLMM_p_values_epi_pneumo.csv')
eff_mv_h2_epi.to_csv('mvLMM_effect_sizes_h2_epi_pneumo.csv')




#### Random permutation GWAS - epistatic

## for loop for running all epistatic interaction tests
no_columns_to_test = 2
a_list = list(range(0, no_columns_to_test))

pv_mv_epi_full = ["lml0", "lml2","dof20","scale2","pv20"]
pv_mv_epi_full = pandas.DataFrame(columns = pv_mv_epi_full)

eff_mv_h2_epi_full = ["test", "trait", "effect_type","effect_name","env","effsize","effsize_se"]
eff_mv_h2_epi_full = pandas.DataFrame(columns = eff_mv_h2_epi_full)

p = 2
A0 = ones((p, 1))
A1 = eye(p)

b_list = list(range(0, 2))

for f in b_list:
    test_markers_sample = test_markers_epi.sample(n=no_columns_to_test, axis=1)
    pv_mv_epi = ["lml0", "lml2", "dof20", "scale2", "pv20"]
    pv_mv_epi = pandas.DataFrame(columns=pv_mv_epi)
    eff_mv_h2_epi = ["test", "trait", "effect_type", "effect_name", "env", "effsize", "effsize_se"]
    eff_mv_h2_epi = pandas.DataFrame(columns=eff_mv_h2_epi)

    map_coords = map_coords_obs
    map_coords = map_coords.sample(frac=1).reset_index(drop=True)


    for i in a_list:
        test_epistatic_marker = test_markers_sample.iloc[:, i].to_frame()
        colnames = test_epistatic_marker.columns[0]
        colnames = colnames.split(':')
        M = dummy_gen[colnames]
        # M = M.to_numpy()
        dropped = dummy_gen.copy()
        dropped.drop(colnames, axis=1, inplace=True)
        dropped = dropped.to_numpy()
        K = linear_kinship(dropped)
        result = scan(G=test_epistatic_marker, Y=map_coords, K=K, M=M,  A=map_coords.cov(), verbose = False)
        pv_mv_epi = pv_mv_epi.append(result.stats)
        eff_mv_h2_epi = eff_mv_h2_epi.append(result.effsizes['h2'])

    map_coords['repeat_p_index'] = pandas.DataFrame(numpy.full((map_coords.shape[0]), f + 1))
    pv_mv_epi['repeat_p_index'] = pandas.DataFrame(numpy.full((pv_mv_epi.shape[0]), f + 1))
    eff_mv_h2_epi['repeat_eff_index'] = pandas.DataFrame(numpy.full((eff_mv_h2_epi.shape[0]), f + 1))

    map_coords_full = map_coords_full.append(map_coords)
    pv_mv_epi_full = pv_mv_epi_full.append(pv_mv_epi)
    eff_mv_h2_epi_full = eff_mv_h2_epi_full.append(eff_mv_h2_epi)

pv_mv_epi.pv20.min()

map_coords_full.to_csv('mvLMM_map_coords_normal_pneumo_random_phenotype_EPI_FWAS.csv')
pv_mv_epi_full.to_csv('mvLMM_p_values_normal_pneumo_random_phenotype_EPI_FWAS.csv')
eff_mv_h2_epi_full.to_csv('mvLMM_effect_sizes_normal_pneumo_random_phenotype_EPI_FWAS.csv')







########### univariate LMM with per test variation decomposition
### Run LMM for each drug on the logtransformed MIC values


MIC_values = limix.io.csv.read('S.pneumo_map_mvlmm_MIC_values.csv', verbose=False)
#a_list = list(range(0, genetic_markers.shape[1]))
a_list = list(range(0, 1))

pv_uv_d1 = ["lml0", "lml2","dof20","scale2","pv20"]
pv_uv_d1 = pandas.DataFrame(columns = pv_uv_d1)

eff_uv_d1 = ["trait", "effect_type","effect_name","effsize","effsize_se"]
eff_uv_d1 = pandas.DataFrame(columns = eff_uv_d1)

b_list = list(range(0, MIC_values.shape[1]))

for f in b_list:
    drug_being_tested = MIC_values.iloc[:, f].to_frame()

    for i in a_list:
        test_marker = genetic_markers.iloc[:, i].to_frame()
        K = dummy_gen
        K = K.drop(test_marker.columns[[0]], axis=1)
        K = K.to_numpy()
        K = linear_kinship(K)
        result = scan(test_marker, drug_being_tested, "normal", K=K, M=None, verbose = False)

        pv_uv_d1 = pv_uv_d1.append(result.stats)
        eff_uv_d1 = eff_uv_d1.append(result.effsizes['h2'])


pv_uv_d1.to_csv('uniLMM_p_val_normal_MIC_pneumo.csv')
eff_uv_d1.to_csv('uniLMM_effect_normal_MIC_pneumo.csv')






### Additional inputs for genotype to phenotype prediction model
import scipy as sp
import numpy as np
import limix_legacy.modules.varianceDecomposition as var


import pylab as pl
from limix.stats import multivariate_normal as mvn
import sys
import scipy.linalg
import scipy.stats
import limix_legacy.deprecated as dlimix_legacy
import limix_legacy.deprecated.utils.preprocess as preprocess
import pdb
import time
import copy
import warnings

import pylab as pl
import scipy.stats as st

import limix_legacy.modules.qtl as qtl
import limix_legacy.io.data as data
import limix_legacy.io.genotype_reader as gr
import limix_legacy.io.phenotype_reader as phr
import limix_legacy.io.data_util as data_util
import limix_legacy.utils.preprocess as preprocess
import limix_legacy.modules.BLUP as blups_test
from y10k_prediction.train_and_test_sets import get_Itrain_Itest
from y10k_prediction.train_and_test_sets import get_CV_ind
from y10k_prediction.train_and_test_sets import select_subset
from sklearn import linear_model
import scipy.linalg as LA
import limix_legacy.modules.lmm_fast as lmm_fast
import limix_legacy.modules.BLUP as blup
from limix_legacy.utils.plot import *
from limix_legacy.stats.geno_summary import *
import limix_legacy.deprecated.io.data_util as du


## BLUP - best linear unbiased prediction of phenotypes from PBP type

def get_Itrain_Itest(N, proportion = 0.8, seed=0):
    sp.random.seed(seed)
    n1 = int(np.floor(N*proportion))
    n2 = int(N - n1)
    selected_indexes = np.concatenate((np.ones(n1, dtype=bool), np.zeros(n2, dtype=bool)))
    sp.random.shuffle(selected_indexes)
    Itrain = selected_indexes
    Itest = ~selected_indexes
    return Itrain, Itest

def get_CV_ind(Nobs, n_folds=5):
    r = sp.random.permutation(Nobs)
    Icv = np.floor(((np.ones((Nobs))*n_folds)*r)/Nobs)
    return Icv, n_folds

Icv, n_folds = get_CV_ind(map_coords.shape[0], n_folds =4)
Itrain, Itest = get_Itrain_Itest(map_coords.shape[0], proportion = 0.8, seed = 32)



def get_BLUPs(Y, K, Itrain=None, Itest=None):
    if Itrain is None:
        Itrain = np.ones(Y.shape[0], dtype=bool)
        Itest = np.ones(Y.shape[0], dtype=bool)
    m = var.VarianceDecomposition(Y[Itrain])
    m.setTestSampleSize(Itest.sum())
    m.addFixedEffect()
    m.addRandomEffect(K=K[Itrain, :][:, Itrain], Kcross=K[Itrain, :][:, Itest])
    m.addRandomEffect(is_noise=True)
    m.optimize()
    blups = m.predictPhenos()
    return blups

def get_BLUPs_with_confidence(Y, K, Itrain, Itest):
    m = var.VarianceDecomposition(Y[Itrain])
    m.setTestSampleSize(Itest.sum())
    m.addFixedEffect()
    m.addRandomEffect(K=K[Itrain, :][:, Itrain], Kcross=K[Itrain, :][:, Itest])
    m.addRandomEffect(is_noise=True)
    m.optimize()
    blups = m.predictPhenos()
    varcomps = m.getVarianceComps().ravel()
    Sigma = varcomps[0]*K + varcomps[1]*np.eye(Y.shape[0])
    Sigma_train_inv = np.linalg.inv(Sigma[Itrain, :][:, Itrain])
    var_predictive = Sigma[Itest, :][:, Itest] - np.dot(np.dot(Sigma[Itest, :][:, Itrain], Sigma_train_inv), Sigma[Itrain, :][:, Itest])
    return {"pred": blups, "predictive_sd": np.sqrt(np.diag(var_predictive))}


LABID = limix.io.csv.read('S.pneumo_map_mvlmm_LABID.csv', verbose=False)


## for loop for running all epistatic interaction tests
Icv, n_folds = get_CV_ind(map_coords.shape[0], n_folds =4)
a_list = list(range(1, n_folds))

Itest = (Icv==0)
Itrain = (Icv!=0)

blups_output_with_conf = get_BLUPs_with_confidence(map_coords['D1'], K, Itrain=Itrain, Itest= Itest)
blups_output_with_conf_V2 = get_BLUPs_with_confidence(map_coords['D2'], K, Itrain=Itrain, Itest= Itest)

LABID_test = pandas.DataFrame(np.array(LABID['LABID'][Itest]))
fold_1_test_V1 = pandas.DataFrame(np.array(map_coords['D1'][Itest]))
fold_1_test_V2 = pandas.DataFrame(np.array(map_coords['D2'][Itest]))
predictions_V1 = pandas.DataFrame(blups_output_with_conf['pred'])
SE_V1 = pandas.DataFrame(blups_output_with_conf['predictive_sd'])
predictions_V2 = pandas.DataFrame(blups_output_with_conf_V2['pred'])
SE_V2 = pandas.DataFrame(blups_output_with_conf_V2['predictive_sd'])
fold_index = pandas.DataFrame(np.full((fold_1_test_V1.shape[0]), 1))
fold_blups = pd.concat([LABID_test, fold_1_test_V1, fold_1_test_V2, predictions_V1, SE_V1,predictions_V2, SE_V2, fold_index], axis=1)


for i in a_list:
    Itest = (Icv==i)
    Itrain = (Icv!=i)

    blups_output_with_conf = get_BLUPs_with_confidence(map_coords['D1'], K, Itrain=Itrain, Itest= Itest)
    blups_output_with_conf_V2 = get_BLUPs_with_confidence(map_coords['D2'], K, Itrain=Itrain, Itest= Itest)

    LABID_test = pandas.DataFrame(np.array(LABID['LABID'][Itest]))
    fold_1_test_V1 = pandas.DataFrame(np.array(map_coords['D1'][Itest]))
    fold_1_test_V2 = pandas.DataFrame(np.array(map_coords['D2'][Itest]))
    predictions_V1 = pandas.DataFrame(blups_output_with_conf['pred'])
    SE_V1 = pandas.DataFrame(blups_output_with_conf['predictive_sd'])
    predictions_V2 = pandas.DataFrame(blups_output_with_conf_V2['pred'])
    SE_V2 = pandas.DataFrame(blups_output_with_conf_V2['predictive_sd'])
    fold_index = pandas.DataFrame(np.full((fold_1_test_V1.shape[0]), i+1))
    current_fold_blups = pd.concat([LABID_test,fold_1_test_V1, fold_1_test_V2, predictions_V1, SE_V1,predictions_V2, SE_V2, fold_index], axis=1)
    fold_blups = fold_blups.append(current_fold_blups)


fold_blups = fold_blups.set_axis(['LABID','D1_true_value', 'D2_true_value', 'D1_prediction', 'D1_prediction_SE', 'D2_prediction', 'D2_prediction_SE', 'Fold'], axis=1, inplace=False)


fold_blups.to_csv('blups_predictions_4fold_S.pneumo.csv')
Icv = pandas.DataFrame(Icv)
Icv.to_csv('Fold_index_for_S.pneumo.csv')
