#! /usr/bin/env python
#Writed by LiHongxiang on 4/26/2022
#tensorflow = 2.7.0
#python = 3.10.4
#keras = 2.7.0
#h5py = 3.6.0
#numpy = 1.22.3
#pandas = 1.4.2
#conda install tensorflow=2.7.0 -c conda-forge
#conda install pandas -c conda-forge
import my_CLOUD as c
import argparse
parser = argparse.ArgumentParser(description="Run CLOUD !")
parser.add_argument("-S", help="SingleData")
parser.add_argument("-D", help="Duplicate Data")

parser.add_argument("-m", type=int, help="Condition numbers of expression data ")
parser.add_argument("-Nk", type=int, help="Number of training data for each class")

parser.add_argument("-T","--training_prefix",help="Training file prefix")
parser.add_argument("-E","--testing_prefix",help="Testing file prefix")

parser.add_argument("-B","--batchsize", type=int, help="Batch size")
parser.add_argument("-N","--num_epochs", type=int, help="Number of epochs")

parser.add_argument("-LMin","--log_lambda_min", default=-5, help="Log lambda min")
parser.add_argument("-LMax","--log_lambda_max", default=-1, help="Log lambda max")
parser.add_argument("-NL","--num_lambda", type=int, default=5, help="Number of lambda")

parser.add_argument("-GMin","--gamma_min", default=0, help="Gamma min")
parser.add_argument("-GMax","--gamma_max", default=1,help="Gamma max")
parser.add_argument("-NG","--num_gamma", default=3, type=int, help="Number of gamma")

parser.add_argument("-DL","--dup_list",help="Dupgene id list,format:Dup1\tDup2\tAncestor, header is needed")

parser.add_argument("-GenerateTrainingData", type=bool, default=0, help="If the parameter is set to 1, it will only run GenerateTrainingData! -m, -Nk, -S, -D, -T are required!")
parser.add_argument("-ClassifierCV", type=bool, default=0, help="If the parameter is set to 1, it will only run ClassifierCV! -m, -B, -N, -T are required! If -LMin, -LMax, -NL, -GMin, -GMax, -NG not setted, they were setted as default values:-5, -1, 5, 0, 1,3.")
parser.add_argument("-PredictorCV", type=bool, default=0, help="If the parameter is set to 1, it will only run PredictorCV! -m, -B, -N, -T are required! If -LMin, -LMax, -NL, -GMin, -GMax, -NG not setted, they were setted as default values:-5, -1, 5, 0, 1,3.")
parser.add_argument("-GenerateTestingFeatures", default=0, type=bool, help="If the parameter is set to 1, it will only run GenerateTestingFeatures! -m, -S, -D, -E are required!")
parser.add_argument("-CLOUDClassify", type=bool, default=0, help="If the parameter is set to 1, it will only run CLOUDClassify! -T, -E are required!")
parser.add_argument("-CLOUDPredict", type=bool, default=0, help="If the parameter is set to 1, it will only run CLOUDPredict! -T, -E are required!")
parser.add_argument("-Format_results", type=bool, default=0, help="If the parameter is set to 1, it will only run Format_results! -E, -DL are required!")
args = parser.parse_args()

m = args.m
Nk = args.Nk
batch_size = args.batchsize
num_epochs = args.num_epochs

log_lambda_max = float(args.log_lambda_max)
log_lambda_min = float(args.log_lambda_min)
num_lambda = args.num_lambda

gamma_min = float(args.gamma_min)
gamma_max = float(args.gamma_max)
num_gamma = args.num_gamma

if args.GenerateTrainingData:
    print('\033[36mStarting GenerateTrainingData......\n\033[0m')
    c.GenerateTrainingData(m, Nk, args.S, args.D, args.training_prefix)
    print('\033[34mGenerateTrainingData Done!\n\033[0m')
elif args.ClassifierCV:
    print('\033[36mStarting ClassifierCV......\n\033[0m')
    c.ClassifierCV(m, batch_size, num_epochs, log_lambda_min, log_lambda_max, num_lambda, gamma_min, gamma_max, num_gamma, args.training_prefix)
    print('\033[34mClassifierCV Done!\n\033[0m')
elif args.PredictorCV:
    print('\033[36mStarting PredictorCV......\n\033[0m')
    c.PredictorCV(m, batch_size, num_epochs, log_lambda_min, log_lambda_max, num_lambda, gamma_min, gamma_max, num_gamma, args.training_prefix)
    print('\033[34mPredictorCV Done!\n\033[0m')
elif args.GenerateTestingFeatures:
    print('\033[36mStarting GenerateTestingFeatures......\n\033[0m')
    c.GenerateTestingFeatures(m, args.S, args.D, args.testing_prefix)
    print('\033[34mGenerateTestingFeatures Done!\n\033[0m')
elif args.CLOUDClassify:
    print('\033[36mStarting CLOUDClassify......\n\033[0m')
    c.CLOUDClassify(args.training_prefix, args.testing_prefix)
    print('\033[34mCLOUDClassify Done!\n\033[0m')
elif args.CLOUDPredict:
    print('\033[36mStarting CLOUDPredict......\n\033[0m')
    c.CLOUDPredict(args.training_prefix, args.testing_prefix)
    print('\033[34mCLOUDPredict Done!\n\033[0m')
elif args.Format_results:
    print('\033[36mStarting Format_results......\n\033[0m')
    c.Format_results(args.testing_prefix, args.dup_list)
    print('\033[34mFormat_results Done!\n\033[0m')
else:
    print('\033[36mStarting GenerateTrainingData......\n\033[0m')
    c.GenerateTrainingData(m, Nk, args.S, args.D, args.training_prefix)
    print('\033[34mGenerateTrainingData Done!\n\033[0m')

    print('\033[36mStarting ClassifierCV......\n\033[0m')
    c.ClassifierCV(m, batch_size, num_epochs, log_lambda_min, log_lambda_max, num_lambda, gamma_min, gamma_max, num_gamma, args.training_prefix)
    print('\033[34mClassifierCV Done!\n\033[0m')

    print('\033[36mStarting PredictorCV......\n\033[0m')
    c.PredictorCV(m, batch_size, num_epochs, log_lambda_min, log_lambda_max, num_lambda, gamma_min, gamma_max, num_gamma, args.training_prefix)
    print('\033[34mPredictorCV Done!\n\033[0m')

    print('\033[36mStarting GenerateTestingFeatures......\n\033[0m')
    c.GenerateTestingFeatures(m, args.S, args.D, args.testing_prefix)
    print('\033[34mGenerateTestingFeatures Done!\n\033[0m')

    print('\033[36mStarting CLOUDClassify......\n\033[0m')
    c.CLOUDClassify(args.training_prefix, args.testing_prefix)
    print('\033[34mCLOUDClassify Done!\n\033[0m')

    print('\033[36mStarting CLOUDPredict......\n\033[0m')
    c.CLOUDPredict(args.training_prefix, args.testing_prefix)
    print('\033[34mCLOUDPredict Done!\n\033[0m')

    print('\033[36mStarting Format_results......\n\033[0m')
    c.Format_results(args.testing_prefix, args.dup_list)
    print('\033[34mFormat_results Done!\n\033[0m')
