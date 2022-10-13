#! /usr/bin/env python
#Writed by LiHongxiang on 26/4/2022
#The input duplicate_file don't need to be logged, the script can transfrom
#输入的testing duplicate_file不再需要自己进行log处理
#Add the function which can correspond the dup ancestor gene id with the results: classificaions and predictions
#The dupgene_list format is:Dup1\tDup2\tAncestor, header is needed
import numpy as np
import pandas as pd
import random
import math
import tensorflow
from keras.models import Sequential
from keras.layers import Dense
from keras import regularizers
from keras.models import load_model

def GenerateTrainingData(m, Nk, singlecopy_filename, duplicate_filename, training_prefix):
    #生成随机数据根据的参数
    #sub类型会用到，取0.5表示parent和child对合并表达量贡献相同
    propPvc = 0.5
    #logThetaMin表达下限
    logThetaMin = -4
    #logThetaMax表达上限
    logThetaMax = 4
    logAlphaMin = 0
    logAlphaMax = 3
    logSigmaSqMin = -2
    logSigmaSqMax = 3

    scenarios = ['cons','neoparent','neochild','sub','spec']
    tissues = [i for i in range(1,m+1)]
    training_data = np.zeros((5*Nk,2+8*m))
    #得到duplicate file中的分化时间比：tPC
    set_of_tPC = []
    with open(duplicate_filename) as f:
        for lines in f.readlines()[1:]:
            line = lines.strip().split()
            tPC = float(line[0])/float(line[1])
            set_of_tPC.append(tPC)
    #为每种类型的功能化生成符合其性质的随机数据，每一行都生成一组数据，每一类Nk行，共5类,每次运行生成的数据都不一样
    for sindex in range(len(scenarios)):
        scenario = scenarios[sindex]
        for i in range(Nk):
            tPCA = 1
            tPC = float(random.sample(set_of_tPC,1)[0])
            alpha = 10 ** np.random.uniform(logAlphaMin,logAlphaMax,(m))
            sigmaSq = 10 ** np.random.uniform(logSigmaSqMin,logSigmaSqMax,(m))

            thetaP = []
            thetaC = []
            thetaA = []
    #生成随机的参数
            #cons的P、C、A三者要求一致
            if scenario == 'cons':
                thetaA = np.random.uniform(logThetaMin,logThetaMax,(m))
                thetaP = thetaA
                thetaC = thetaA
            #neoparent要求P、A有差别
            elif scenario == 'neoparent':
                thetaA = np.random.uniform(logThetaMin,logThetaMax,(m))
                thetaP = np.random.uniform(logThetaMin,logThetaMax,(m))
                thetaC = thetaA
            #neochild要求P、C有差别
            elif scenario == 'neochild':
                thetaA = np.random.uniform(logThetaMin,logThetaMax,(m))
                thetaC = np.random.uniform(logThetaMin,logThetaMax,(m))
                thetaP = thetaA
            #sub要求P、C二者加权相加后和A一致
            elif scenario == 'sub':
                thetaP = np.random.uniform(logThetaMin,logThetaMax,(m))
                thetaC = np.random.uniform(logThetaMin,logThetaMax,(m))
                thetaA = propPvc * thetaP + (1-propPvc) * thetaC
            #spec要求P、C、A都不一致
            elif scenario == 'spec':
                thetaA = np.random.uniform(logThetaMin, logThetaMax, (m))
                thetaP = np.random.uniform(logThetaMin, logThetaMax, (m))
                thetaC = np.random.uniform(logThetaMin, logThetaMax, (m))
            #生成随机的表达数据
            expression_vec = []
            mu = [0,0,0]
            CovMat = np.zeros((3,3))
            #根据alpha，tPC来模拟进化过程中分化时间不同，而造成的表达量差异大小不同，主要想模拟的是tPC分化时间造成的影响，alpha就是一个约束参数，通过对theta的加权，模拟到真实的进化场景
            #tPC取值0到1，如果tPC=0.9，这就说明parent和child很早就分开了，他们之间的差别就会大些，因为tPC大，所以‘1 - math.exp(-alpha[j] * tPC)’的值就接近1，thetaP、C起主导作用，‘math.exp(-alpha[j] * tPC)’值接近0，thetaA参与的影响很小
            #如果tPC很小，比如0.000005，说明parent和child刚分开，那么他们即使有新功能产生，他们的表达量差别也不会特别大，所以‘math.exp(-alpha[j] * tPC)’的值会接近1，thetaA加权之后在表达量中占主要部分，这样P和C的表达量的差异就会缩小，符合进化中的实况
            #math.exp(-alpha[j] * tPC) * thetaA[j]相当于在表达量中加入一部分稳定的值，置于这个稳定的值占多少比例，由tPC决定
            for j in range(m):
                mu[0] = (1 - math.exp(-alpha[j] * tPC)) * thetaP[j] + math.exp(-alpha[j] * tPC) * thetaA[j]
                mu[1] = (1 - math.exp(-alpha[j] * tPC)) * thetaC[j] + math.exp(-alpha[j] * tPC) * thetaA[j]
                mu[2] = thetaA[j]
                CovMat[0, 0] = sigmaSq[j] / (2 * alpha[j])
                CovMat[1, 1] = CovMat[0, 0]
                CovMat[2, 2] = CovMat[0, 0]
                CovMat[0, 1] = math.exp(-2 * alpha[j] * tPC) * sigmaSq[j] / (2 * alpha[j])
                CovMat[1, 0] = CovMat[0, 1]
                CovMat[0, 2] = math.exp(-2 * alpha[j] * tPCA) * sigmaSq[j] / (2 * alpha[j])
                CovMat[2, 0] = CovMat[0, 2]
                CovMat[1, 2] = CovMat[0, 2]
                CovMat[2, 1] = CovMat[0, 2]
                for exp in np.random.multivariate_normal(mean=mu, cov=CovMat, size= 1).reshape(3):
                    expression_vec.append(exp)
                #np.random.multivariate_normal以mean（多维分布的均值）和cov（协方差矩阵），size指定矩阵的shape,生成一个随机正态矩阵
            #用数据填充数据框
            rowIndex = sindex*Nk + i
            training_data[rowIndex, 0] = sindex+1
            training_data[rowIndex, 1] = tPC
            for j in range(3*m):
                training_data[rowIndex, 2 + j] = expression_vec[j]
            for j in range(m):
                training_data[rowIndex, 2 + 3*m + j] = thetaP[j]
            for j in range(m):
                training_data[rowIndex, 2 + 4*m + j] = thetaC[j]
            for j in range(m):
                training_data[rowIndex, 2 + 5*m + j] = thetaA[j]
            for j in range(m):
                training_data[rowIndex, 2 + 6*m + j] = math.log10(alpha[j])
            for j in range(m):
                training_data[rowIndex, 2 + 7*m + j] = math.log10(sigmaSq[j])
    #制作每列的标签
    column_labels = ["Class","tPC"]
    for j in range(1,m+1):
        column_labels.append("eP" + str(j))
        column_labels.append("eC" + str(j))
        column_labels.append("eA" + str(j))
    for j in range(1,m+1):
        column_labels.append("ThetaP" + str(j))
    for j in range(1,m+1):
        column_labels.append("ThetaC" + str(j))
    for j in range(1,m+1):
        column_labels.append("ThetaA" + str(j))
    for j in range(1,m+1):
        column_labels.append("LogAlpha" + str(j))
    for j in range(1,m+1):
        column_labels.append("LogSigmaSq" + str(j))
    #输出到文件
    training_data_out = open(training_prefix+".data",'w')
    print('\t'.join(column_labels),file=training_data_out)
    for j in range(training_data.shape[0]):
        for k in range(training_data.shape[1]):
            if k == 0:
                print(str(int(training_data[j, k])), end='\t', file=training_data_out)
            else:
                print(str(training_data[j,k]),end='\t',file=training_data_out)
        print('\n',end='',file=training_data_out)
    training_data_out.close()
    #调用其他函数，生成其他文件
    GenerateFeatures(m, singlecopy_filename, training_prefix + ".data", training_prefix + ".features")
    GenerateClassifierResponse(training_prefix + ".data", training_prefix + ".classes")
    GeneratePredictorResponse(training_prefix + ".data", training_prefix + ".responses")

def GenerateFeatures(m, singlecopy_filename, input_filename, feature_filename):
    minexp = 10 ** -4
    errorexp = 10 ** -5
    # 读取singledata
    single = pd.read_csv(singlecopy_filename, sep='\t', index_col=False)
    # 将pandas数据框转为np数组，便于运算
    single = np.array(single)
    # Transform data to log scale, accounting for expression of 0
    for i in range(single.shape[0]):
        for j in range(single.shape[1]):
            single[i, j] = math.log10(single[i, j] + minexp + np.random.normal(0, errorexp, 1)[0])

    # Euclidean distance between single-copy expression profiles
    eS1S2dist = np.sqrt(np.sum((single[:, :m] - single[:, m:2 * m]) ** 2, axis=1))
    maxS1S2dist = max(abs(eS1S2dist))
    # 计算S1和S2的相关性
    eS1S2cor = []
    for i in range(single.shape[0]):
        eS1S2cor.append(pd.DataFrame({'S1': single[i, :m], 'S2': single[i, m:2 * m]}).corr().iloc[0, 1])
    eS1S2cor = np.array(eS1S2cor)

    features = pd.read_csv(input_filename, sep='\t', index_col=False)
    eP = ["eP" + str(i) for i in range(1, m + 1)]
    eC = ["eC" + str(i) for i in range(1, m + 1)]
    eA = ["eA" + str(i) for i in range(1, m + 1)]
    tmplist = ["tPC"]
    tmplist.extend(eP)
    tmplist.extend(eC)
    tmplist.extend(eA)
    features = features.reindex(columns=tmplist)

    eSumPC = features.copy().loc[:, eP].values + features.copy().loc[:, eC].values
    eSumPC_labels = ["eSumPC" + str(i) for i in range(1, m + 1)]
    eSumPC = pd.DataFrame(eSumPC, columns=eSumPC_labels)
    features = features.join(eSumPC)
    # Euclidean distance between P and C
    eDistPvC = np.sqrt(np.sum((features.loc[:, eP].values - features.loc[:, eC].values) ** 2, axis=1))
    # Euclidean distance between P and A
    eDistPvA = np.sqrt(np.sum((features.loc[:, eP].values - features.loc[:, eA].values) ** 2, axis=1))
    # Euclidean distance between C and A
    eDistCvA = np.sqrt(np.sum((features.loc[:, eC].values - features.loc[:, eA].values) ** 2, axis=1))
    # Euclidean distance between PC and A
    eDistPCvA = np.sqrt(np.sum((features.loc[:, eSumPC_labels].values - features.loc[:, eA].values) ** 2, axis=1))
    features["DistPvC"] = eDistPvC
    features["DistPvA"] = eDistPvA
    features["DistCvA"] = eDistCvA
    features["DistPCvA"] = eDistPCvA
    features["PBS_P"] = (features["DistPvC"] + features["DistPvA"] - features["DistCvA"]) / 2
    features["PBS_C"] = (features["DistPvC"] + features["DistCvA"] - features["DistPvA"]) / 2
    features["PBS_At"] = (features["DistPvA"] + features["DistCvA"] - features["DistPvC"]) / 2
    #对数据进行各种转换，相互运算，生成可以充分展示数据的各种属性的各种计算后的数据，这些数据有助于机器学习更精准的识别
    RankDistPvC = []
    for i in features['DistPvC']:
        if eS1S2dist[eS1S2dist < i].size > 0:
            RankDistPvC.append(np.mean(eS1S2dist < i))
        else:
            RankDistPvC.append(i)
    features["RankDistPvC"] = RankDistPvC
    RankDistPvA = []
    for i in features['DistPvA']:
        if eS1S2dist[eS1S2dist < i].size > 0:
            RankDistPvA.append(np.mean(eS1S2dist < i))
        else:
            RankDistPvA.append(i)
    features["RankDistPvA"] = RankDistPvA
    RankDistCvA = []
    for i in features['DistCvA']:
        if eS1S2dist[eS1S2dist < i].size > 0:
            RankDistCvA.append(np.mean(eS1S2dist < i))
        else:
            RankDistCvA.append(i)
    features["RankDistCvA"] = RankDistCvA
    RankDistPCvA = []
    for i in features['DistPCvA']:
        if eS1S2dist[eS1S2dist < i].size > 0:
            RankDistPCvA.append(np.mean(eS1S2dist < i))
        else:
            RankDistPCvA.append(i)
    features["RankDistPCvA"] = RankDistPCvA
    # When center and absolute are both FALSE, the moment is simply sum(x ^ order) / length(x).
    # m1
    DistPvC_m1 = []
    for i in features['DistPvC']:
        DistPvC_m1.append(np.mean((eS1S2dist - i) / maxS1S2dist))
    features['DistPvC_m1'] = DistPvC_m1
    DistPvA_m1 = []
    for i in features['DistPvA']:
        DistPvA_m1.append(np.mean((eS1S2dist - i) / maxS1S2dist))
    features['DistPvA_m1'] = DistPvA_m1
    DistCvA_m1 = []
    for i in features['DistCvA']:
        DistCvA_m1.append(np.mean((eS1S2dist - i) / maxS1S2dist))
    features['DistCvA_m1'] = DistCvA_m1
    DistPCvA_m1 = []
    for i in features['DistPCvA']:
        DistPCvA_m1.append(np.mean((eS1S2dist - i) / maxS1S2dist))
    features['DistPCvA_m1'] = DistPCvA_m1
    # m2
    DistPvC_m2 = []
    for i in features['DistPvC']:
        DistPvC_m2.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 2))
    features['DistPvC_m2'] = DistPvC_m2
    DistPvA_m2 = []
    for i in features['DistPvA']:
        DistPvA_m2.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 2))
    features['DistPvA_m2'] = DistPvA_m2
    DistCvA_m2 = []
    for i in features['DistCvA']:
        DistCvA_m2.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 2))
    features['DistCvA_m2'] = DistCvA_m2
    DistPCvA_m2 = []
    for i in features['DistPCvA']:
        DistPCvA_m2.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 2))
    features['DistPCvA_m2'] = DistPCvA_m2
    # m3
    DistPvC_m3 = []
    for i in features['DistPvC']:
        DistPvC_m3.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 3))
    features['DistPvC_m3'] = DistPvC_m3
    DistPvA_m3 = []
    for i in features['DistPvA']:
        DistPvA_m3.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 3))
    features['DistPvA_m3'] = DistPvA_m3
    DistCvA_m3 = []
    for i in features['DistCvA']:
        DistCvA_m3.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 3))
    features['DistCvA_m3'] = DistCvA_m3
    DistPCvA_m3 = []
    for i in features['DistPCvA']:
        DistPCvA_m3.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 3))
    features['DistPCvA_m3'] = DistPCvA_m3
    # m4
    DistPvC_m4 = []
    for i in features['DistPvC']:
        DistPvC_m4.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 4))
    features['DistPvC_m4'] = DistPvC_m4
    DistPvA_m4 = []
    for i in features['DistPvA']:
        DistPvA_m4.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 4))
    features['DistPvA_m4'] = DistPvA_m4
    DistCvA_m4 = []
    for i in features['DistCvA']:
        DistCvA_m4.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 4))
    features['DistCvA_m4'] = DistCvA_m4
    DistPCvA_m4 = []
    for i in features['DistPCvA']:
        DistPCvA_m4.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 4))
    features['DistPCvA_m4'] = DistPCvA_m4
    # m5
    DistPvC_m5 = []
    for i in features['DistPvC']:
        DistPvC_m5.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 5))
    features['DistPvC_m5'] = DistPvC_m5
    DistPvA_m5 = []
    for i in features['DistPvA']:
        DistPvA_m5.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 5))
    features['DistPvA_m5'] = DistPvA_m5
    DistCvA_m5 = []
    for i in features['DistCvA']:
        DistCvA_m5.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 5))
    features['DistCvA_m5'] = DistCvA_m5
    DistPCvA_m5 = []
    for i in features['DistPCvA']:
        DistPCvA_m5.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 5))
    features['DistPCvA_m5'] = DistPCvA_m5
    # m6
    DistPvC_m6 = []
    for i in features['DistPvC']:
        DistPvC_m6.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 6))
    features['DistPvC_m6'] = DistPvC_m6
    DistPvA_m6 = []
    for i in features['DistPvA']:
        DistPvA_m6.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 6))
    features['DistPvA_m6'] = DistPvA_m6
    DistCvA_m6 = []
    for i in features['DistCvA']:
        DistCvA_m6.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 6))
    features['DistCvA_m6'] = DistCvA_m6
    DistPCvA_m6 = []
    for i in features['DistPCvA']:
        DistPCvA_m6.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 6))
    features['DistPCvA_m6'] = DistPCvA_m6
    # m7
    DistPvC_m7 = []
    for i in features['DistPvC']:
        DistPvC_m7.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 7))
    features['DistPvC_m7'] = DistPvC_m7
    DistPvA_m7 = []
    for i in features['DistPvA']:
        DistPvA_m7.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 7))
    features['DistPvA_m7'] = DistPvA_m7
    DistCvA_m7 = []
    for i in features['DistCvA']:
        DistCvA_m7.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 7))
    features['DistCvA_m7'] = DistCvA_m7
    DistPCvA_m7 = []
    for i in features['DistPCvA']:
        DistPCvA_m7.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 7))
    features['DistPCvA_m7'] = DistPCvA_m7
    # m8
    DistPvC_m8 = []
    for i in features['DistPvC']:
        DistPvC_m8.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 8))
    features['DistPvC_m8'] = DistPvC_m8
    DistPvA_m8 = []
    for i in features['DistPvA']:
        DistPvA_m8.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 8))
    features['DistPvA_m8'] = DistPvA_m8
    DistCvA_m8 = []
    for i in features['DistCvA']:
        DistCvA_m8.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 8))
    features['DistCvA_m8'] = DistCvA_m8
    DistPCvA_m8 = []
    for i in features['DistPCvA']:
        DistPCvA_m8.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 8))
    features['DistPCvA_m8'] = DistPCvA_m8

    eCorPvC = []
    eCorPvA = []
    eCorCvA = []
    eCorPCvA = []

    # Pearson correlations
    for i in range(features.shape[0]):
        eCorPvC.append(pd.DataFrame({'eP': np.array(features.iloc[i, features.columns.str.startswith('eP')]),'eC': np.array(features.iloc[i, features.columns.str.startswith('eC')])}).corr().iloc[0, 1])
    features['CorPvC'] = eCorPvC
    for i in range(features.shape[0]):
        eCorPvA.append(pd.DataFrame({'eP': np.array(features.iloc[i, features.columns.str.startswith('eP')]),'eA': np.array(features.iloc[i, features.columns.str.startswith('eA')])}).corr().iloc[0, 1])
    features['CorPvA'] = eCorPvA
    for i in range(features.shape[0]):
        eCorCvA.append(pd.DataFrame({'eC': np.array(features.iloc[i, features.columns.str.startswith('eC')]),'eA': np.array(features.iloc[i, features.columns.str.startswith('eA')])}).corr().iloc[0, 1])
    features['CorCvA'] = eCorCvA
    for i in range(features.shape[0]):
        eCorPCvA.append(pd.DataFrame({'ePC': np.array(features.iloc[i, features.columns.str.startswith('eSumPC')]),'eA': np.array(features.iloc[i, features.columns.str.startswith('eA')])}).corr().iloc[0, 1])
    features['CorPCvA'] = eCorPCvA

    RankCorPvC = []
    for i in features['CorPvC']:
        if eS1S2cor[eS1S2cor < i].size > 0:
            RankCorPvC.append(np.mean(eS1S2cor < i))
        else:
            RankCorPvC.append(i)
    features["RankCorPvC"] = RankCorPvC

    RankCorPvA = []
    for i in features['CorPvA']:
        if eS1S2cor[eS1S2cor < i].size > 0:
            RankCorPvA.append(np.mean(eS1S2cor < i))
        else:
            RankCorPvA.append(i)
    features["RankCorPvA"] = RankCorPvA

    RankCorCvA = []
    for i in features['CorCvA']:
        if eS1S2cor[eS1S2cor < i].size > 0:
            RankCorCvA.append(np.mean(eS1S2cor < i))
        else:
            RankCorCvA.append(i)
    features["RankCorCvA"] = RankCorCvA

    RankCorPCvA = []
    for i in features['CorPCvA']:
        if eS1S2cor[eS1S2cor < i].size > 0:
            RankCorPCvA.append(np.mean(eS1S2cor < i))
        else:
            RankCorPCvA.append(i)
    features["RankCorPCvA"] = RankCorPCvA

    # m1
    CorPvC_m1 = []
    for i in features['CorPvC']:
        CorPvC_m1.append(np.mean(eS1S2cor - i))
    features['CorPvC_m1'] = CorPvC_m1
    CorPvA_m1 = []
    for i in features['CorPvA']:
        CorPvA_m1.append(np.mean(eS1S2cor - i))
    features['CorPvA_m1'] = CorPvA_m1
    CorCvA_m1 = []
    for i in features['CorCvA']:
        CorCvA_m1.append(np.mean(eS1S2cor - i))
    features['CorCvA_m1'] = CorCvA_m1
    CorPCvA_m1 = []
    for i in features['CorPCvA']:
        CorPCvA_m1.append(np.mean(eS1S2cor - i))
    features['CorPCvA_m1'] = CorPCvA_m1
    # m2
    CorPvC_m2 = []
    for i in features['CorPvC']:
        CorPvC_m2.append(np.mean((eS1S2cor - i) ** 2))
    features['CorPvC_m2'] = CorPvC_m2
    CorPvA_m2 = []
    for i in features['CorPvA']:
        CorPvA_m2.append(np.mean((eS1S2cor - i) ** 2))
    features['CorPvA_m2'] = CorPvA_m2
    CorCvA_m2 = []
    for i in features['CorCvA']:
        CorCvA_m2.append(np.mean((eS1S2cor - i) ** 2))
    features['CorCvA_m2'] = CorCvA_m2
    CorPCvA_m2 = []
    for i in features['CorPCvA']:
        CorPCvA_m2.append(np.mean((eS1S2cor - i) ** 2))
    features['CorPCvA_m2'] = CorPCvA_m2
    # m3
    CorPvC_m3 = []
    for i in features['CorPvC']:
        CorPvC_m3.append(np.mean((eS1S2cor - i) ** 3))
    features['CorPvC_m3'] = CorPvC_m3
    CorPvA_m3 = []
    for i in features['CorPvA']:
        CorPvA_m3.append(np.mean((eS1S2cor - i) ** 3))
    features['CorPvA_m3'] = CorPvA_m3
    CorCvA_m3 = []
    for i in features['CorCvA']:
        CorCvA_m3.append(np.mean((eS1S2cor - i) ** 3))
    features['CorCvA_m3'] = CorCvA_m3
    CorPCvA_m3 = []
    for i in features['CorPCvA']:
        CorPCvA_m3.append(np.mean((eS1S2cor - i) ** 3))
    features['CorPCvA_m3'] = CorPCvA_m3
    # m4
    CorPvC_m4 = []
    for i in features['CorPvC']:
        CorPvC_m4.append(np.mean((eS1S2cor - i) ** 4))
    features['CorPvC_m4'] = CorPvC_m4
    CorPvA_m4 = []
    for i in features['CorPvA']:
        CorPvA_m4.append(np.mean((eS1S2cor - i) ** 4))
    features['CorPvA_m4'] = CorPvA_m4
    CorCvA_m4 = []
    for i in features['CorCvA']:
        CorCvA_m4.append(np.mean((eS1S2cor - i) ** 4))
    features['CorCvA_m4'] = CorCvA_m4
    CorPCvA_m4 = []
    for i in features['CorPCvA']:
        CorPCvA_m4.append(np.mean((eS1S2cor - i) ** 4))
    features['CorPCvA_m4'] = CorPCvA_m4
    # m5
    CorPvC_m5 = []
    for i in features['CorPvC']:
        CorPvC_m5.append(np.mean((eS1S2cor - i) ** 5))
    features['CorPvC_m5'] = CorPvC_m5
    CorPvA_m5 = []
    for i in features['CorPvA']:
        CorPvA_m5.append(np.mean((eS1S2cor - i) ** 5))
    features['CorPvA_m5'] = CorPvA_m5
    CorCvA_m5 = []
    for i in features['CorCvA']:
        CorCvA_m5.append(np.mean((eS1S2cor - i) ** 5))
    features['CorCvA_m5'] = CorCvA_m5
    CorPCvA_m5 = []
    for i in features['CorPCvA']:
        CorPCvA_m5.append(np.mean((eS1S2cor - i) ** 5))
    features['CorPCvA_m5'] = CorPCvA_m5
    # m6
    CorPvC_m6 = []
    for i in features['CorPvC']:
        CorPvC_m6.append(np.mean((eS1S2cor - i) ** 6))
    features['CorPvC_m6'] = CorPvC_m6
    CorPvA_m6 = []
    for i in features['CorPvA']:
        CorPvA_m6.append(np.mean((eS1S2cor - i) ** 6))
    features['CorPvA_m6'] = CorPvA_m6
    CorCvA_m6 = []
    for i in features['CorCvA']:
        CorCvA_m6.append(np.mean((eS1S2cor - i) ** 6))
    features['CorCvA_m6'] = CorCvA_m6
    CorPCvA_m6 = []
    for i in features['CorPCvA']:
        CorPCvA_m6.append(np.mean((eS1S2cor - i) ** 6))
    features['CorPCvA_m6'] = CorPCvA_m6
    # m7
    CorPvC_m7 = []
    for i in features['CorPvC']:
        CorPvC_m7.append(np.mean((eS1S2cor - i) ** 7))
    features['CorPvC_m7'] = CorPvC_m7
    CorPvA_m7 = []
    for i in features['CorPvA']:
        CorPvA_m7.append(np.mean((eS1S2cor - i) ** 7))
    features['CorPvA_m7'] = CorPvA_m7
    CorCvA_m7 = []
    for i in features['CorCvA']:
        CorCvA_m7.append(np.mean((eS1S2cor - i) ** 7))
    features['CorCvA_m7'] = CorCvA_m7
    CorPCvA_m7 = []
    for i in features['CorPCvA']:
        CorPCvA_m7.append(np.mean((eS1S2cor - i) ** 7))
    features['CorPCvA_m7'] = CorPCvA_m7
    # m8
    CorPvC_m8 = []
    for i in features['CorPvC']:
        CorPvC_m8.append(np.mean((eS1S2cor - i) ** 8))
    features['CorPvC_m8'] = CorPvC_m8
    CorPvA_m8 = []
    for i in features['CorPvA']:
        CorPvA_m8.append(np.mean((eS1S2cor - i) ** 8))
    features['CorPvA_m8'] = CorPvA_m8
    CorCvA_m8 = []
    for i in features['CorCvA']:
        CorCvA_m8.append(np.mean((eS1S2cor - i) ** 8))
    features['CorCvA_m8'] = CorCvA_m8
    CorPCvA_m8 = []
    for i in features['CorPCvA']:
        CorPCvA_m8.append(np.mean((eS1S2cor - i) ** 8))
    features['CorPCvA_m8'] = CorPCvA_m8

    # print features
    features_out = open(feature_filename, 'w')
    print('\t'.join(features.columns), file=features_out)
    for j in range(features.shape[0]):
        for k in range(features.shape[1]):
            print(str(features.iloc[j, k]), end='\t', file=features_out)
        print('\n', end='', file=features_out)
    features_out.close()

def GenerateTestingFeatures(m, singlecopy_filename, duplicate_filename, testing_prefix):
    minexp = 10**-4
    errorexp = 10**-5
    #读取singledata
    single = pd.read_csv(singlecopy_filename, sep='\t', index_col=False)
    # 将pandas数据框转为np数组，便于运算
    single = np.array(single)
    #Transform data to log scale, accounting for expression of 0
    for i in range(single.shape[0]):
        for j in range(single.shape[1]):
            single[i,j] = math.log10(single[i,j] + minexp + np.random.normal(0, errorexp, 1)[0])

    # Euclidean distance between single-copy expression profiles
    eS1S2dist = np.sqrt(np.sum((single[:,:m] - single[:,m:2*m]) ** 2, axis = 1))
    maxS1S2dist = max(abs(eS1S2dist))
    #计算S1和S2的相关性
    eS1S2cor = []
    for i in range(single.shape[0]):
        eS1S2cor.append(pd.DataFrame({'S1':single[i,:m],'S2':single[i,m:2*m]}).corr().iloc[0,1])
    eS1S2cor = np.array(eS1S2cor)

    features = pd.read_csv(duplicate_filename, sep='\t', index_col=False)
    eP = ["eP" + str(i) for i in range(1, m + 1)]
    eC = ["eC" + str(i) for i in range(1, m + 1)]
    eA = ["eA" + str(i) for i in range(1, m + 1)]
    tmplist = ["tPC"]
    tmplist.extend(eP)
    tmplist.extend(eC)
    tmplist.extend(eA)
    if 'tPC' in features.columns:
        features = features.reindex(columns = tmplist)
        # Transform data to log scale, accounting for expression of 0
        for i in range(features.shape[0]):
            for j in range(1,features.shape[1]):
                features.iloc[i, j] = math.log10(features.iloc[i, j] + minexp + np.random.normal(0, errorexp, 1)[0])
    else:
        features['tPC'] = features['TPC'] / features['TPCA']
        features = features.reindex(columns = tmplist)
        # Transform data to log scale, accounting for expression of 0
        for i in range(features.shape[0]):
            for j in range(1,features.shape[1]):
                features.iloc[i, j] = math.log10(features.iloc[i, j] + minexp + np.random.normal(0, errorexp, 1)[0])

    eSumPC = features.copy().loc[:,eP].values + features.copy().loc[:,eC].values
    eSumPC_labels = ["eSumPC"+str(i) for i in range(1,m+1)]
    eSumPC = pd.DataFrame(eSumPC,columns=eSumPC_labels)
    features = features.join(eSumPC)
    # Euclidean distance between P and C
    eDistPvC = np.sqrt(np.sum((features.loc[:,eP].values - features.loc[:,eC].values) ** 2, axis = 1))
    # Euclidean distance between P and A
    eDistPvA = np.sqrt(np.sum((features.loc[:,eP].values - features.loc[:,eA].values) ** 2, axis = 1))
    # Euclidean distance between C and A
    eDistCvA = np.sqrt(np.sum((features.loc[:,eC].values - features.loc[:,eA].values) ** 2, axis = 1))
    # Euclidean distance between PC and A
    eDistPCvA = np.sqrt(np.sum((features.loc[:,eSumPC_labels].values - features.loc[:,eA].values) ** 2, axis = 1))
    features["DistPvC"] = eDistPvC
    features["DistPvA"] = eDistPvA
    features["DistCvA"] = eDistCvA
    features["DistPCvA"] = eDistPCvA
    features["PBS_P"] = (features["DistPvC"] + features["DistPvA"] - features["DistCvA"]) / 2
    features["PBS_C"] = (features["DistPvC"] + features["DistCvA"] - features["DistPvA"]) / 2
    features["PBS_At"] = (features["DistPvA"] + features["DistCvA"] - features["DistPvC"]) / 2
    RankDistPvC = []
    for i in features['DistPvC']:
        if eS1S2dist[eS1S2dist < i].size > 0:
            RankDistPvC.append(np.mean(eS1S2dist < i))
        else:
            RankDistPvC.append(i)
    features["RankDistPvC"] = RankDistPvC
    RankDistPvA = []
    for i in features['DistPvA']:
        if eS1S2dist[eS1S2dist < i].size > 0:
            RankDistPvA.append(np.mean(eS1S2dist < i))
        else:
            RankDistPvA.append(i)
    features["RankDistPvA"] = RankDistPvA
    RankDistCvA = []
    for i in features['DistCvA']:
        if eS1S2dist[eS1S2dist < i].size > 0:
            RankDistCvA.append(np.mean(eS1S2dist < i))
        else:
            RankDistCvA.append(i)
    features["RankDistCvA"] = RankDistCvA
    RankDistPCvA = []
    for i in features['DistPCvA']:
        if eS1S2dist[eS1S2dist < i].size > 0:
            RankDistPCvA.append(np.mean(eS1S2dist < i))
        else:
            RankDistPCvA.append(i)
    features["RankDistPCvA"] = RankDistPCvA
#When center and absolute are both FALSE, the moment is simply sum(x ^ order) / length(x).
    # m1
    DistPvC_m1 = []
    for i in features['DistPvC']:
        DistPvC_m1.append(np.mean((eS1S2dist - i) / maxS1S2dist))
    features['DistPvC_m1'] = DistPvC_m1
    DistPvA_m1 = []
    for i in features['DistPvA']:
        DistPvA_m1.append(np.mean((eS1S2dist - i) / maxS1S2dist))
    features['DistPvA_m1'] = DistPvA_m1
    DistCvA_m1 = []
    for i in features['DistCvA']:
        DistCvA_m1.append(np.mean((eS1S2dist - i) / maxS1S2dist))
    features['DistCvA_m1'] = DistCvA_m1
    DistPCvA_m1 = []
    for i in features['DistPCvA']:
        DistPCvA_m1.append(np.mean((eS1S2dist - i) / maxS1S2dist))
    features['DistPCvA_m1'] = DistPCvA_m1
    # m2
    DistPvC_m2 = []
    for i in features['DistPvC']:
        DistPvC_m2.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 2))
    features['DistPvC_m2'] = DistPvC_m2
    DistPvA_m2 = []
    for i in features['DistPvA']:
        DistPvA_m2.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 2))
    features['DistPvA_m2'] = DistPvA_m2
    DistCvA_m2 = []
    for i in features['DistCvA']:
        DistCvA_m2.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 2))
    features['DistCvA_m2'] = DistCvA_m2
    DistPCvA_m2 = []
    for i in features['DistPCvA']:
        DistPCvA_m2.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 2))
    features['DistPCvA_m2'] = DistPCvA_m2
    # m3
    DistPvC_m3 = []
    for i in features['DistPvC']:
        DistPvC_m3.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 3))
    features['DistPvC_m3'] = DistPvC_m3
    DistPvA_m3 = []
    for i in features['DistPvA']:
        DistPvA_m3.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 3))
    features['DistPvA_m3'] = DistPvA_m3
    DistCvA_m3 = []
    for i in features['DistCvA']:
        DistCvA_m3.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 3))
    features['DistCvA_m3'] = DistCvA_m3
    DistPCvA_m3 = []
    for i in features['DistPCvA']:
        DistPCvA_m3.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 3))
    features['DistPCvA_m3'] = DistPCvA_m3
    # m4
    DistPvC_m4 = []
    for i in features['DistPvC']:
        DistPvC_m4.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 4))
    features['DistPvC_m4'] = DistPvC_m4
    DistPvA_m4 = []
    for i in features['DistPvA']:
        DistPvA_m4.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 4))
    features['DistPvA_m4'] = DistPvA_m4
    DistCvA_m4 = []
    for i in features['DistCvA']:
        DistCvA_m4.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 4))
    features['DistCvA_m4'] = DistCvA_m4
    DistPCvA_m4 = []
    for i in features['DistPCvA']:
        DistPCvA_m4.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 4))
    features['DistPCvA_m4'] = DistPCvA_m4
    # m5
    DistPvC_m5 = []
    for i in features['DistPvC']:
        DistPvC_m5.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 5))
    features['DistPvC_m5'] = DistPvC_m5
    DistPvA_m5 = []
    for i in features['DistPvA']:
        DistPvA_m5.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 5))
    features['DistPvA_m5'] = DistPvA_m5
    DistCvA_m5 = []
    for i in features['DistCvA']:
        DistCvA_m5.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 5))
    features['DistCvA_m5'] = DistCvA_m5
    DistPCvA_m5 = []
    for i in features['DistPCvA']:
        DistPCvA_m5.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 5))
    features['DistPCvA_m5'] = DistPCvA_m5
    # m6
    DistPvC_m6 = []
    for i in features['DistPvC']:
        DistPvC_m6.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 6))
    features['DistPvC_m6'] = DistPvC_m6
    DistPvA_m6 = []
    for i in features['DistPvA']:
        DistPvA_m6.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 6))
    features['DistPvA_m6'] = DistPvA_m6
    DistCvA_m6 = []
    for i in features['DistCvA']:
        DistCvA_m6.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 6))
    features['DistCvA_m6'] = DistCvA_m6
    DistPCvA_m6 = []
    for i in features['DistPCvA']:
        DistPCvA_m6.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 6))
    features['DistPCvA_m6'] = DistPCvA_m6
    # m7
    DistPvC_m7 = []
    for i in features['DistPvC']:
        DistPvC_m7.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 7))
    features['DistPvC_m7'] = DistPvC_m7
    DistPvA_m7 = []
    for i in features['DistPvA']:
        DistPvA_m7.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 7))
    features['DistPvA_m7'] = DistPvA_m7
    DistCvA_m7 = []
    for i in features['DistCvA']:
        DistCvA_m7.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 7))
    features['DistCvA_m7'] = DistCvA_m7
    DistPCvA_m7 = []
    for i in features['DistPCvA']:
        DistPCvA_m7.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 7))
    features['DistPCvA_m7'] = DistPCvA_m7
    # m8
    DistPvC_m8 = []
    for i in features['DistPvC']:
        DistPvC_m8.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 8))
    features['DistPvC_m8'] = DistPvC_m8
    DistPvA_m8 = []
    for i in features['DistPvA']:
        DistPvA_m8.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 8))
    features['DistPvA_m8'] = DistPvA_m8
    DistCvA_m8 = []
    for i in features['DistCvA']:
        DistCvA_m8.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 8))
    features['DistCvA_m8'] = DistCvA_m8
    DistPCvA_m8 = []
    for i in features['DistPCvA']:
        DistPCvA_m8.append(np.mean(((eS1S2dist - i) / maxS1S2dist) ** 8))
    features['DistPCvA_m8'] = DistPCvA_m8

    eCorPvC = []
    eCorPvA = []
    eCorCvA = []
    eCorPCvA = []

#Pearson correlations
    for i in range(features.shape[0]):
            eCorPvC.append(pd.DataFrame({'eP': np.array(features.iloc[i, features.columns.str.startswith('eP')]), 'eC': np.array(features.iloc[i, features.columns.str.startswith('eC')])}).corr().iloc[0, 1])
    features['CorPvC'] = eCorPvC
    for i in range(features.shape[0]):
            eCorPvA.append(pd.DataFrame({'eP': np.array(features.iloc[i, features.columns.str.startswith('eP')]), 'eA': np.array(features.iloc[i, features.columns.str.startswith('eA')])}).corr().iloc[0, 1])
    features['CorPvA'] = eCorPvA
    for i in range(features.shape[0]):
            eCorCvA.append(pd.DataFrame({'eC': np.array(features.iloc[i, features.columns.str.startswith('eC')]), 'eA': np.array(features.iloc[i, features.columns.str.startswith('eA')])}).corr().iloc[0, 1])
    features['CorCvA'] = eCorCvA
    for i in range(features.shape[0]):
            eCorPCvA.append(pd.DataFrame({'ePC': np.array(features.iloc[i, features.columns.str.startswith('eSumPC')]), 'eA': np.array(features.iloc[i, features.columns.str.startswith('eA')])}).corr().iloc[0, 1])
    features['CorPCvA'] = eCorPCvA

    RankCorPvC = []
    for i in features['CorPvC']:
        if eS1S2cor[eS1S2cor < i].size > 0:
            RankCorPvC.append(np.mean(eS1S2cor < i))
        else:
            RankCorPvC.append(i)
    features["RankCorPvC"] = RankCorPvC

    RankCorPvA = []
    for i in features['CorPvA']:
        if eS1S2cor[eS1S2cor < i].size > 0:
            RankCorPvA.append(np.mean(eS1S2cor < i))
        else:
            RankCorPvA.append(i)
    features["RankCorPvA"] = RankCorPvA

    RankCorCvA = []
    for i in features['CorCvA']:
        if eS1S2cor[eS1S2cor < i].size > 0:
            RankCorCvA.append(np.mean(eS1S2cor < i))
        else:
            RankCorCvA.append(i)
    features["RankCorCvA"] = RankCorCvA

    RankCorPCvA = []
    for i in features['CorPCvA']:
        if eS1S2cor[eS1S2cor < i].size > 0:
            RankCorPCvA.append(np.mean(eS1S2cor < i))
        else:
            RankCorPCvA.append(i)
    features["RankCorPCvA"] = RankCorPCvA

    # m1
    CorPvC_m1 = []
    for i in features['CorPvC']:
        CorPvC_m1.append(np.mean(eS1S2cor - i))
    features['CorPvC_m1'] = CorPvC_m1
    CorPvA_m1 = []
    for i in features['CorPvA']:
        CorPvA_m1.append(np.mean(eS1S2cor - i))
    features['CorPvA_m1'] = CorPvA_m1
    CorCvA_m1 = []
    for i in features['CorCvA']:
        CorCvA_m1.append(np.mean(eS1S2cor - i))
    features['CorCvA_m1'] = CorCvA_m1
    CorPCvA_m1 = []
    for i in features['CorPCvA']:
        CorPCvA_m1.append(np.mean(eS1S2cor - i))
    features['CorPCvA_m1'] = CorPCvA_m1
    # m2
    CorPvC_m2 = []
    for i in features['CorPvC']:
        CorPvC_m2.append(np.mean((eS1S2cor - i) ** 2))
    features['CorPvC_m2'] = CorPvC_m2
    CorPvA_m2 = []
    for i in features['CorPvA']:
        CorPvA_m2.append(np.mean((eS1S2cor - i) ** 2))
    features['CorPvA_m2'] = CorPvA_m2
    CorCvA_m2 = []
    for i in features['CorCvA']:
        CorCvA_m2.append(np.mean((eS1S2cor - i) ** 2))
    features['CorCvA_m2'] = CorCvA_m2
    CorPCvA_m2 = []
    for i in features['CorPCvA']:
        CorPCvA_m2.append(np.mean((eS1S2cor - i) ** 2))
    features['CorPCvA_m2'] = CorPCvA_m2
    # m3
    CorPvC_m3 = []
    for i in features['CorPvC']:
        CorPvC_m3.append(np.mean((eS1S2cor - i) ** 3))
    features['CorPvC_m3'] = CorPvC_m3
    CorPvA_m3 = []
    for i in features['CorPvA']:
        CorPvA_m3.append(np.mean((eS1S2cor - i) ** 3))
    features['CorPvA_m3'] = CorPvA_m3
    CorCvA_m3 = []
    for i in features['CorCvA']:
        CorCvA_m3.append(np.mean((eS1S2cor - i) ** 3))
    features['CorCvA_m3'] = CorCvA_m3
    CorPCvA_m3 = []
    for i in features['CorPCvA']:
        CorPCvA_m3.append(np.mean((eS1S2cor - i) ** 3))
    features['CorPCvA_m3'] = CorPCvA_m3
    # m4
    CorPvC_m4 = []
    for i in features['CorPvC']:
        CorPvC_m4.append(np.mean((eS1S2cor - i) ** 4))
    features['CorPvC_m4'] = CorPvC_m4
    CorPvA_m4 = []
    for i in features['CorPvA']:
        CorPvA_m4.append(np.mean((eS1S2cor - i) ** 4))
    features['CorPvA_m4'] = CorPvA_m4
    CorCvA_m4 = []
    for i in features['CorCvA']:
        CorCvA_m4.append(np.mean((eS1S2cor - i) ** 4))
    features['CorCvA_m4'] = CorCvA_m4
    CorPCvA_m4 = []
    for i in features['CorPCvA']:
        CorPCvA_m4.append(np.mean((eS1S2cor - i) ** 4))
    features['CorPCvA_m4'] = CorPCvA_m4
    # m5
    CorPvC_m5 = []
    for i in features['CorPvC']:
        CorPvC_m5.append(np.mean((eS1S2cor - i) ** 5))
    features['CorPvC_m5'] = CorPvC_m5
    CorPvA_m5 = []
    for i in features['CorPvA']:
        CorPvA_m5.append(np.mean((eS1S2cor - i) ** 5))
    features['CorPvA_m5'] = CorPvA_m5
    CorCvA_m5 = []
    for i in features['CorCvA']:
        CorCvA_m5.append(np.mean((eS1S2cor - i) ** 5))
    features['CorCvA_m5'] = CorCvA_m5
    CorPCvA_m5 = []
    for i in features['CorPCvA']:
        CorPCvA_m5.append(np.mean((eS1S2cor - i) ** 5))
    features['CorPCvA_m5'] = CorPCvA_m5
    # m6
    CorPvC_m6 = []
    for i in features['CorPvC']:
        CorPvC_m6.append(np.mean((eS1S2cor - i) ** 6))
    features['CorPvC_m6'] = CorPvC_m6
    CorPvA_m6 = []
    for i in features['CorPvA']:
        CorPvA_m6.append(np.mean((eS1S2cor - i) ** 6))
    features['CorPvA_m6'] = CorPvA_m6
    CorCvA_m6 = []
    for i in features['CorCvA']:
        CorCvA_m6.append(np.mean((eS1S2cor - i) ** 6))
    features['CorCvA_m6'] = CorCvA_m6
    CorPCvA_m6 = []
    for i in features['CorPCvA']:
        CorPCvA_m6.append(np.mean((eS1S2cor - i) ** 6))
    features['CorPCvA_m6'] = CorPCvA_m6
    # m7
    CorPvC_m7 = []
    for i in features['CorPvC']:
        CorPvC_m7.append(np.mean((eS1S2cor - i) ** 7))
    features['CorPvC_m7'] = CorPvC_m7
    CorPvA_m7 = []
    for i in features['CorPvA']:
        CorPvA_m7.append(np.mean((eS1S2cor - i) ** 7))
    features['CorPvA_m7'] = CorPvA_m7
    CorCvA_m7 = []
    for i in features['CorCvA']:
        CorCvA_m7.append(np.mean((eS1S2cor - i) ** 7))
    features['CorCvA_m7'] = CorCvA_m7
    CorPCvA_m7 = []
    for i in features['CorPCvA']:
        CorPCvA_m7.append(np.mean((eS1S2cor - i) ** 7))
    features['CorPCvA_m7'] = CorPCvA_m7
    # m8
    CorPvC_m8 = []
    for i in features['CorPvC']:
        CorPvC_m8.append(np.mean((eS1S2cor - i) ** 8))
    features['CorPvC_m8'] = CorPvC_m8
    CorPvA_m8 = []
    for i in features['CorPvA']:
        CorPvA_m8.append(np.mean((eS1S2cor - i) ** 8))
    features['CorPvA_m8'] = CorPvA_m8
    CorCvA_m8 = []
    for i in features['CorCvA']:
        CorCvA_m8.append(np.mean((eS1S2cor - i) ** 8))
    features['CorCvA_m8'] = CorCvA_m8
    CorPCvA_m8 = []
    for i in features['CorPCvA']:
        CorPCvA_m8.append(np.mean((eS1S2cor - i) ** 8))
    features['CorPCvA_m8'] = CorPCvA_m8

    #print features
    features_out = open(testing_prefix+'.features','w')
    print('\t'.join(features.columns),file=features_out)
    for j in range(features.shape[0]):
        for k in range(features.shape[1]):
            print(str(features.iloc[j, k]), end='\t', file=features_out)
        print('\n', end='', file=features_out)
    features_out.close()

def GenerateClassifierResponse(input_filename, response_filename):
    training_data = pd.read_csv(input_filename, sep='\t', index_col=False)
    response_out = open(response_filename,'w')
    print('cons', 'neoparent', 'neochild', 'sub', 'spec', sep='\t', file=response_out)
    response = np.zeros([training_data.iloc[:,0].size,5], dtype='int')
    for i in range(response.shape[0]):
        response[i, int(training_data.iloc[i,0])-1] = 1
    for j in range(response.shape[0]):
        for k in range(response.shape[1]):
            print(response[j,k],end='\t',file=response_out)
        print('\n',end='',file=response_out)
    response_out.close()

def GeneratePredictorResponse(input_filename, response_filename):
    training_data = pd.read_csv(input_filename, sep='\t', index_col=False)
    response_out = open(response_filename,'w')
    theta = training_data.loc[:, training_data.columns.str.startswith('Theta')]
    logalpha = training_data.loc[:, training_data.columns.str.startswith('LogAlpha')]
    logsigma = training_data.loc[:, training_data.columns.str.startswith('LogSigma')]
    response = theta.join(logalpha.join(logsigma))
    print('\t'.join(response.columns), file=response_out)
    for j in range(response.shape[0]):
        for k in range(response.shape[1]):
            print(response.iloc[j,k], end = '\t', file=response_out)
        print('\n', end='', file=response_out)
    response_out.close()

def ClassifierCV(m, batchsize, num_epochs, log_lambda_min, log_lambda_max, num_lambda, gamma_min, gamma_max, num_gamma, training_prefix):
    CV = 5
    lambdas = 10 ** np.linspace(log_lambda_min, log_lambda_max, num_lambda)
    gammas = np.linspace(gamma_min, gamma_max, num_gamma)

    X = pd.read_csv(training_prefix+'.features', sep = '\t', index_col=False)
    Y = pd.read_csv(training_prefix+'.classes', sep = '\t', index_col=False)

    # standardize the input for training（标准化）
    X_means = X.mean(axis=0)
    #计算每一列的平均值
    X_sds = X.std(axis=0)
    #计算矩阵每一列的标准差
    for j in range(X.shape[1]):
        X.iloc[:,j] = (X.iloc[:,j] - X_means[j]) / X_sds[j]

    # CV-fold cross validation
    val_loss = np.zeros([num_gamma, num_lambda])

    # Randomly choose balanced training/validation sets per fold
    #80%训练，20%预测
    #分成5类取随机是为了保证每一类中都均匀的从1-5随机，如果合在一起取，就不能保证每一类都有固定的1-5的值
    set_of_value = [1,2,3,4,5]

    foldid_cons = []
    while len(foldid_cons) < X.shape[0] / 5:
        foldid_cons = foldid_cons + random.sample(set_of_value,5)
    foldid_cons = foldid_cons[:int(X.shape[0] / 5)]

    foldid_neoparent = []
    while len(foldid_neoparent) < X.shape[0] / 5:
        foldid_neoparent = foldid_neoparent + random.sample(set_of_value,5)
    foldid_neoparent = foldid_neoparent[:int(X.shape[0] / 5)]

    foldid_neochild = []
    while len(foldid_neochild) < X.shape[0] / 5:
        foldid_neochild = foldid_neochild + random.sample(set_of_value,5)
    foldid_neochild = foldid_neochild[:int(X.shape[0] / 5)]

    foldid_sub = []
    while len(foldid_sub) < X.shape[0] / 5:
        foldid_sub = foldid_sub + random.sample(set_of_value,5)
    foldid_sub = foldid_sub[:int(X.shape[0] / 5)]

    foldid_spec = []
    while len(foldid_spec) < X.shape[0] / 5:
        foldid_spec = foldid_spec + random.sample(set_of_value, 5)
    foldid_spec = foldid_spec[:int(X.shape[0] / 5)]

    foldid = np.array(foldid_cons + foldid_neoparent +foldid_neochild + foldid_sub + foldid_spec)

    # Perform K-fold CV, where K = CV
    for curr_fold in range(1,CV+1):
        Xval = X.copy().iloc[foldid == curr_fold, :]
        Xtrain = X.copy().iloc[foldid != curr_fold, :]
        Yval = Y.copy().iloc[foldid == curr_fold, :]
        Ytrain = Y.copy().iloc[foldid != curr_fold, :]

        # standardize the input for train and val based on train
        temp_means = Xtrain.mean(axis=0)
        temp_sds = Xtrain.std(axis=0)
        for j in range(Xtrain.shape[1]):
            Xtrain.iloc[:,j] = (Xtrain.iloc[:,j] - temp_means[j]) / temp_sds[j]
            Xval.iloc[:,j] = (Xval.iloc[:,j] - temp_means[j]) / temp_sds[j]
        #参数两两组合，每对参数都构建模型进行训练
        for i in range(num_gamma):
            for j in range(num_lambda):
                model = Sequential()
                model.add(Dense(units = 256,
                                activation = 'relu',
                                kernel_regularizer = regularizers.l1_l2(l1 = gammas[i] * lambdas[j], l2 = (1 - gammas[i]) * lambdas[j]),
                                input_shape = (Xtrain.shape[1],)))
                model.add(Dense(units = 128,
                                activation = 'relu',
                                kernel_regularizer = regularizers.l1_l2(l1 = gammas[i] * lambdas[j], l2 = (1 - gammas[i]) * lambdas[j])))
                model.add(Dense(units = 5,
                                activation = 'softmax',
                                kernel_regularizer = regularizers.l1_l2(l1 = gammas[i] * lambdas[j], l2 = (1 - gammas[i]) * lambdas[j])))

                model.compile(loss = 'categorical_crossentropy',
                              optimizer = 'Adam',
                              metrics = ['categorical_crossentropy'])

                history = model.fit(Xtrain, Ytrain, epochs = num_epochs, batch_size = batchsize, validation_data = (Xval, Yval), verbose = 0)

                val_loss[i,j] = val_loss[i,j] + pd.DataFrame(history.history).loc[:,'val_categorical_crossentropy'].min()

    val_loss = val_loss / CV
    gamma_opt = gammas[np.where(val_loss == np.min(val_loss))[0][0]]
    lambda_opt = lambdas[np.where(val_loss == np.min(val_loss))[1][0]]
    cv_results = []
    cv_results.append(str(np.min(val_loss)))
    cv_results.append(str(gamma_opt))
    cv_results.append(str(lambda_opt))

    classifier_cv_out = open(training_prefix+'.classifier_cv', 'w')
    print('Loss','Gamma','Lambda',sep='\t',file=classifier_cv_out)
    print('\t'.join(cv_results),file=classifier_cv_out)
    classifier_cv_out.close()

    ClassifierFit(m,batchsize,num_epochs,training_prefix)

def PredictorCV(m, batchsize, num_epochs, log_lambda_min, log_lambda_max, num_lambda, gamma_min, gamma_max, num_gamma, training_prefix):
    CV = 5
    lambdas = 10 ** np.linspace(log_lambda_min, log_lambda_max, num_lambda)
    gammas = np.linspace(gamma_min, gamma_max, num_gamma)

    X = pd.read_csv(training_prefix+'.features', sep = '\t', index_col=False)
    Y = pd.read_csv(training_prefix+'.responses', sep = '\t', index_col=False)

    # standardize the input for training（标准化）
    X_means = X.mean(axis=0)
    X_sds = X.std(axis=0)
    Y_means = Y.mean(axis=0)
    Y_sds = Y.std(axis=0)

    for j in range(X.shape[1]):
        X.iloc[:,j] = (X.iloc[:,j] - X_means[j]) / X_sds[j]
    for j in range(Y.shape[1]):
        Y.iloc[:,j] = (Y.iloc[:,j] - Y_means[j]) / Y_sds[j]

    # CV-fold cross validation
    val_loss = np.zeros([num_gamma, num_lambda])

    # Randomly choose balanced training/validation sets per fold
    set_of_value = [1,2,3,4,5]

    foldid_cons = []
    while len(foldid_cons) < X.shape[0] / 5:
        foldid_cons = foldid_cons + random.sample(set_of_value,5)
    foldid_cons = foldid_cons[:int(X.shape[0] / 5)]

    foldid_neoparent = []
    while len(foldid_neoparent) < X.shape[0] / 5:
        foldid_neoparent = foldid_neoparent + random.sample(set_of_value,5)
    foldid_neoparent = foldid_neoparent[:int(X.shape[0] / 5)]

    foldid_neochild = []
    while len(foldid_neochild) < X.shape[0] / 5:
        foldid_neochild = foldid_neochild + random.sample(set_of_value,5)
    foldid_neochild = foldid_neochild[:int(X.shape[0] / 5)]

    foldid_sub = []
    while len(foldid_sub) < X.shape[0] / 5:
        foldid_sub = foldid_sub + random.sample(set_of_value,5)
    foldid_sub = foldid_sub[:int(X.shape[0] / 5)]

    foldid_spec = []
    while len(foldid_spec) < X.shape[0] / 5:
        foldid_spec = foldid_spec + random.sample(set_of_value, 5)
    foldid_spec = foldid_spec[:int(X.shape[0] / 5)]

    foldid = np.array(foldid_cons + foldid_neoparent +foldid_neochild + foldid_sub + foldid_spec)

    # Perform K-fold CV, where K = CV
    for curr_fold in range(1,CV+1):
        Xval = X.copy().iloc[foldid == curr_fold, :]
        Xtrain = X.copy().iloc[foldid != curr_fold, :]
        Yval = Y.copy().iloc[foldid == curr_fold, :]
        Ytrain = Y.copy().iloc[foldid != curr_fold, :]

        # standardize the input for train and val based on train
        temp_Xmeans = Xtrain.mean(axis=0)
        temp_Xsds = Xtrain.std(axis=0)
        temp_Ymeans = Ytrain.mean(axis=0)
        temp_Ysds = Ytrain.std(axis=0)

        for j in range(Xtrain.shape[1]):
            Xtrain.iloc[:,j] = (Xtrain.iloc[:,j] - temp_Xmeans[j]) / temp_Xsds[j]
            Xval.iloc[:,j] = (Xval.iloc[:,j] - temp_Xmeans[j]) / temp_Xsds[j]
        for j in range(Ytrain.shape[1]):
            Ytrain.iloc[:,j] = (Ytrain.iloc[:,j] - temp_Ymeans[j]) / temp_Ysds[j]
            Yval.iloc[:,j] = (Yval.iloc[:,j] - temp_Ymeans[j]) / temp_Ysds[j]

        for i in range(num_gamma):
            for j in range(num_lambda):
                model = Sequential()
                model.add(Dense(units = 256,
                                activation = 'relu',
                                kernel_regularizer = regularizers.l1_l2(l1 = gammas[i] * lambdas[j], l2 = (1 - gammas[i]) * lambdas[j]),
                                input_shape = (Xtrain.shape[1],)))
                model.add(Dense(units = 128,
                                activation = 'relu',
                                kernel_regularizer = regularizers.l1_l2(l1 = gammas[i] * lambdas[j], l2 = (1 - gammas[i]) * lambdas[j])))
                model.add(Dense(units = Ytrain.shape[1],
                                activation = 'linear',
                                kernel_regularizer = regularizers.l1_l2(l1 = gammas[i] * lambdas[j], l2 = (1 - gammas[i]) * lambdas[j])))

                model.compile(loss = 'mean_squared_error',
                              optimizer = 'Adam',
                              metrics = ['mean_squared_error'])

                history = model.fit(Xtrain, Ytrain, epochs = num_epochs, batch_size = batchsize, validation_data = (Xval, Yval), verbose = 0)

                val_loss[i, j] = val_loss[i, j] + pd.DataFrame(history.history).loc[:,'val_mean_squared_error'].min()

    val_loss = val_loss / CV
    #where 用来取得val_loss中的最小值的索引，gamma_opt,取最小值行索引的gamma值，lambda_opt,取最小值的列索引的lambda值
    gamma_opt = gammas[np.where(val_loss == np.min(val_loss))[0][0]]
    lambda_opt = lambdas[np.where(val_loss == np.min(val_loss))[1][0]]
    cv_results = []
    cv_results.append(str(np.min(val_loss)))
    cv_results.append(str(gamma_opt))
    cv_results.append(str(lambda_opt))

    predictor_cv_out = open(training_prefix+'.predictor_cv', 'w')
    print('Loss','Gamma','Lambda', sep='\t', file=predictor_cv_out)
    print('\t'.join(cv_results), file=predictor_cv_out)
    predictor_cv_out.close()

    PredictorFit(m,batchsize,num_epochs,training_prefix)

def ClassifierFit(m,batchsize,num_epochs,training_prefix):
    cv_results = pd.read_csv(training_prefix+'.classifier_cv',sep='\t', index_col=False)
    gamma_opt = cv_results.loc[0,'Gamma']
    lambda_opt = cv_results.loc[0,'Lambda']

    X = pd.read_csv(training_prefix+'.features',sep='\t',index_col=False)
    Y = pd.read_csv(training_prefix+'.classes',sep='\t',index_col=False)

    # standardize the input for training
    X_means = X.mean(axis=0)
    X_sds = X.std(axis=0)
    for j in range(X.shape[1]):
        X.iloc[:, j] = (X.iloc[:, j] - X_means[j]) / X_sds[j]
    #make DataFrame
    std_params = pd.DataFrame({"Xmeans":X_means,"Xsds":X_sds})
    #save to file
    std_params.to_csv(training_prefix+'.X_stdparams',sep='\t',index=False)

    model = Sequential()
    model.add(Dense(units=256,
                    activation='relu',
                    kernel_regularizer=regularizers.l1_l2(l1=gamma_opt * lambda_opt, l2=(1 - gamma_opt) * lambda_opt),
                    input_shape=(X.shape[1],)))
    model.add(Dense(units=128,
                    activation='relu',
                    kernel_regularizer=regularizers.l1_l2(l1=gamma_opt * lambda_opt, l2=(1 - gamma_opt) * lambda_opt)))
    model.add(Dense(units=5,
                    activation='softmax',
                    kernel_regularizer=regularizers.l1_l2(l1=gamma_opt * lambda_opt, l2=(1 - gamma_opt) * lambda_opt)))

    model.compile(loss='categorical_crossentropy',
                  optimizer='Adam',
                  metrics=['categorical_crossentropy','accuracy'])

    history = model.fit(X, Y, epochs=num_epochs, batch_size=batchsize,verbose=0)

    model.save(training_prefix+'.classifier.hdf5')

def PredictorFit(m,batchsize,num_epochs,training_prefix):
    cv_results = pd.read_csv(training_prefix + '.predictor_cv', sep='\t', index_col=False)
    gamma_opt = cv_results.loc[0, 'Gamma']
    lambda_opt = cv_results.loc[0, 'Lambda']

    X = pd.read_csv(training_prefix+'.features',sep='\t',index_col=False)
    Y = pd.read_csv(training_prefix+'.responses',sep='\t',index_col=False)

    # standardize the input for training
    X_means = X.mean(axis=0)
    X_sds = X.std(axis=0)
    Y_means = Y.mean(axis=0)
    Y_sds = Y.std(axis=0)
    for j in range(X.shape[1]):
        X.iloc[:, j] = (X.iloc[:, j] - X_means[j]) / X_sds[j]
    for j in range(Y.shape[1]):
        Y.iloc[:, j] = (Y.iloc[:, j] - Y_means[j]) / Y_sds[j]

    std_params = pd.DataFrame({"Xmeans": X_means, "Xsds": X_sds})
    std_params.to_csv(training_prefix + '.X_stdparams', sep='\t', index=False)
    std_params = pd.DataFrame({"Ymeans": Y_means, "Ysds": Y_sds})
    std_params.to_csv(training_prefix + '.Y_stdparams', sep='\t', index=False)

    model = Sequential()
    model.add(Dense(units=256,
                    activation='relu',
                    kernel_regularizer=regularizers.l1_l2(l1=gamma_opt * lambda_opt, l2=(1 - gamma_opt) * lambda_opt),
                    input_shape=(X.shape[1],)))
    model.add(Dense(units=128,
                    activation='relu',
                    kernel_regularizer=regularizers.l1_l2(l1=gamma_opt * lambda_opt, l2=(1 - gamma_opt) * lambda_opt)))
    model.add(Dense(units=Y.shape[1],
                    activation='linear',
                    kernel_regularizer=regularizers.l1_l2(l1=gamma_opt * lambda_opt, l2=(1 - gamma_opt) * lambda_opt)))

    model.compile(loss='mean_squared_error',
                  optimizer='Adam',
                  metrics=['mean_squared_error'])

    history = model.fit(X, Y, epochs=num_epochs, batch_size=batchsize, verbose=0)

    model.save(training_prefix + '.predictor.hdf5')

def CLOUDClassify(training_prefix, testing_prefix):
    X = pd.read_csv(testing_prefix + '.features', sep='\t', index_col=False)
    # standardize the input for training
    std_params = pd.read_csv(training_prefix+'.X_stdparams', sep='\t', index_col=False)
    X_means = std_params['Xmeans']
    X_sds = std_params['Xsds']
    for j in range(X.shape[1]):
        X.iloc[:, j] = (X.iloc[:, j] - X_means[j]) / X_sds[j]
    #loading model
    model = load_model(training_prefix+'.classifier.hdf5')
    #starting predicting
    #predict_classes was deprecated in tensorflow 2.6, so 'predict_classes' was replaced 'predict() and np.argmax()'
    #Yest_num = model.predict_classes(X)
    Yest_proba = model.predict(X)
    Yest_num = np.argmax(Yest_proba, axis=1)
    #writing outfile:classifications
    classifications_out = open(testing_prefix+'.classifications','w')
    print('Classifications',file=classifications_out)
    for i in Yest_num:
        if i == 0:
            print('Conservation',file=classifications_out)
        elif i == 1:
            print('Neofunctionalization(Parent)',file=classifications_out)
        elif i == 2:
            print('Neofunctionalization(Child)',file=classifications_out)
        elif i == 3:
            print('Subfunctionalization',file=classifications_out)
        else:
            print('Specialization',file=classifications_out)
    classifications_out.close()

def CLOUDPredict(training_prefix, testing_prefix):
    X = pd.read_csv(testing_prefix + '.features', sep='\t', index_col=False)
    # standardize the input for training
    std_params = pd.read_csv(training_prefix + '.X_stdparams', sep='\t', index_col=False)
    X_means = std_params['Xmeans']
    X_sds = std_params['Xsds']
    for j in range(X.shape[1]):
        X.iloc[:, j] = (X.iloc[:, j] - X_means[j]) / X_sds[j]

    std_params = pd.read_csv(training_prefix + '.Y_stdparams', sep='\t', index_col=False)
    Y_means = std_params['Ymeans']
    Y_sds = std_params['Ysds']

    # loading model
    model = load_model(training_prefix + '.predictor.hdf5')

    #starting predicting
    Yest_std = model.predict(X)
    Yest = Yest_std.copy()

    for j in range(Yest.shape[1]):
        Yest[:, j] = Yest_std[:, j] * Y_sds[j] + Y_means[j]

    Y = pd.read_csv(training_prefix+'.responses', sep = '\t', index_col=False)
    Yest = pd.DataFrame(Yest, columns = Y.columns)
    Yest.to_csv(testing_prefix+'.predictions', sep='\t', index=False)

def Format_results(testing_prefix, dupgene_id_list):
    classifications = pd.read_csv(testing_prefix+'.classifications',index_col=False)
    predictions = pd.read_csv(testing_prefix+'.predictions', sep = '\t', index_col=False)
    dup_list = pd.read_csv(dupgene_id_list, sep = '\t', index_col=False)
    classifications_formatted = dup_list.join(classifications)
    predictions_formatted = dup_list.join(predictions)
    classifications_formatted.to_csv(testing_prefix+'.final_classifications',sep='\t', index=False)
    predictions_formatted.to_csv(testing_prefix+'.final_predictions',sep='\t', index=False)
    d = {}
    for i in classifications.values:
        if i[0] in d:
            d[i[0]] += 1
        else:
            d[i[0]] = 1
    d = pd.Series(d)
    d.to_csv(testing_prefix+'.classifications_count',sep='\t',header=False)
    









