# python_CLOUD
python_CLOUD is a python version adapted from the R implementation of DeGiorgio and Assis (2020) multilayer neural networks for classifying repetitive gene retention mechanisms and predicting their evolutionary parameters from gene expression data in two species.

	
---------------	
Getting started
---------------
Before you are able to use the python_CLOUD , you will need to have the numpy, pandas, random, math， tensorflow and keras packages installed in your python environment. These packages can be installed with the following commands in python:

  pip install numpy
  pip install pandas
  pip install random
  pip install math
  pip install tensorflow
  pip install keras

Our python version is 3.10.4, numpy version is 1.22.3, pandas version is 1.4.2, TensorFlow version is 2.7.0, Keras version is 2.7.0.
The python_CLOUD package comes with the script my_CLOUD.py and run_CLOUD.py. 
An ExampleFiles directory containing example files to help you get started.
  
The python_CLOUD package can be used by:

python run_CLOUD.py -h

Note that this command assumes that you are in the directory of run_cloud.py script, and that this directory contains the my_CLOUD script.

To run python_CLOUD, you will need to provide two input files containing expression data for single-copy and duplicate genes. The formats of these files are described in the next two sections, and example files are provided in the ExampleFiles directory.

----------------------------------
Single-copy gene input file format
----------------------------------
A tab-delimited file with N+1 rows and 2m columns, where N is the number of genes and m is the number of conditions (e.g., tissues, ages, etc.) for expression data. The first row is a header.

Each row (rows 2 to N+1) is a single-copy gene, and each column is an absolute expression value at a particular condition in a particular gene. The first m columns are the gene expression values for the m conditions in the first species, and the next m columns are the gene expression values for the m conditions in the second species. Specifically, columns j and m+j represent the gene expression values at condition j (j = 1, 2, ..., m) in species 1 and species 2, respectively.

The header takes the form:

  eS1t1 eS1t2 ... eS1tm eS2t1 eS2t2 ... eS2tm,

where S1 and S2 denote species 1 and 2, respectively, and t1, t2, ..., tm denote conditions 1, 2, ..., m, respectively.
  
The ExampleFiles directory contains a file called SingleData, which illustates this format for N=5133 single-copy genes at m=6 tissues in two species.

--------------------------------
Duplicate gene input file format
--------------------------------
A tab-delimited file with N+1 rows and 3m+2 columns, where N is the number of genes and m is the number of conditions (e.g., tissues, ages, etc.) for expression data. The first row is a header.

Each row (rows 2 to N+1) is a duplicate gene, and each column is an absolute expression value at a particular condition in a particular gene. The first column is the divergence time between the parent and child genes on each gene tree, the second column is the divergence time between the parent and child gene with the ancestral gene on each gene tree, the next m columns (columns 3 to m+2) are the gene expression values for the m conditions in the parent gene, the next m columns (columns m+3 to 2m+2) are the gene expression values for the m conditions in the child gene, and the last m columns (columns 2m+3 to 3m+2) are the gene expression values for the m conditions in the ancestral gene. Specifically, columns j+2, m+j+2, and 2m+j+2 represent the gene expression values at condition j (j = 1, 2, ..., m) in the parent, child, and ancestral genes, respectively.

The header takes the form:

TPC TPCA  eP1 eP2 ... ePm eC1 eC2 ... eCm eA1 eA2 ... eAm

where TPC is the divergence time between the parent and child genes on the gene tree, TPCA is the divergence time between the parent and child genes with the ancestral gene on the gene tree, and where eP1 to ePm denote parent gene expression values for conditions 1 to m, eC1 to eCm denote child gene expression values for conditions 1 to m, and eA1 to eAm denote ancestral gene expression values for conditions 1 to m.

The ExampleFiles directory contains a file called DuplicateData, which illustates this format for N=100 duplicate genes at m=6 tissues.


-------------------------------------------
Generating training data from an OU process
-------------------------------------------
Given input single-copy and duplicate gene datasets as described in "Single-copy gene input file format" and "Duplicate gene input file format", respectively, training data can be generated with the command:

  python run_CLOUD.py -GenerateTrainingData 1 -m 6 -Nk 100 -S singlecopy_filename -D duplicate_filename -T training_prefix
  
  eg. python run_CLOUD.py -GenerateTrainingData 1 -m 6 -Nk 100 -S SingleData -D DuplicateData -T training
  
where m is the number of conditions, Nk is the number of training observations for each of the five classes, singlecopy_filename is the dataset of single-copy genes described in "Single-copy gene input file format", duplicate_filename is the dataset of duplicate genes described in "Duplicate gene input file format", and training_prefix is the prefix to all files output by this function.

-----------------------------
Training the python_CLOUD classifier
-----------------------------
Given training data outputted by the 'GenerateTrainingData' function as described in "Generating training dta from an OU processes", we can perform five-fold cross-validation to identify appropriate hyper parameters and fit the python_CLOUD classifier with two hidden layers (see DeGiorgio and Assis 2020) under optimal hyper parameters lambda and gamma with the command:
  
  python run_CLOUD.py -ClassifierCV 1 -m 6 -B batchsize -N num_epochs -T training_prefix -LMin log_lambda_min -LMax log_lambda_max -NL num_lambda -GMin gamma_min -GMax gamma_max -NG num_gamma 
  If -LMin, -LMax, -NL, -GMin, -GMax, -NG not setted, they were setted as default values:-5, -1, 5, 0, 1,3.
  
  eg. python run_CLOUD.py -ClassifierCV 1 -m 6 -B 50 -N 50 -T training
  
where m is the number of conditions, batchsize is the number of training observations used in each epoch, num_epochs is the number of training epochs, hyper parameter lambda is drawn from log10(lambda) in interval [log_lambda_min, log_lambda_max] for num_lambda evenly spaced points, hyper parameter gamma is drawn from interval [gamma_min, gamma_max] for num_gamma evenly spaced points, and training_prefix is the prefix to all files outputted by this function.

The ClassifierCV function outputs the set of optimal hyper parameters chosen through cross-validation to the file training_prefix.classifier_cv containing, the means and standard deviations of the input features used for standardizing the input during training to the file training_prefix.X_stdprams, and the fitted model in TensorFlow format to the file training_prefix.classifier.hdf5.

----------------------------
Training the CLOUD predictor
----------------------------
Given training data outputted by the 'GenerateTrainingData' function as described in "Generating training dta from an OU processes", we can perform five-fold cross-validation to identify appropriate hyper parameters and fit the python_CLOUD predictor with two hidden layers (see DeGiorgio and Assis 2020) under optimal hyper parameters lambda and gamma with the command:

  python run_CLOUD.py -PredictorCV 1 -m 6 -B batchsize -N num_epochs -T training_prefix -LMin log_lambda_min -LMax log_lambda_max -NL num_lambda -GMin gamma_min -GMax gamma_max -NG num_gamma 
  If -LMin, -LMax, -NL, -GMin, -GMax, -NG not setted, they were setted as default values:-5, -1, 5, 0, 1,3.
  
  eg. python run_CLOUD.py -PredictorCV 1 -m 6 -B 50 -N 50 -T training

where m is the number of conditions, batchsize is the number of training observations used in each epoch, num_epochs is the number of training epochs, hyper parameter lambda is drawn from log10(lambda) in interval [log_lambda_min, log_lambda_max] for num_lambda evenly spaced points, hyper parameter gamma is drawn from interval [gamma_min, gamma_max] for num_gamma evenly spaced points, and training_prefix is the prefix to all files outputted by this function.

The PredictorCV function outputs the set of optimal hyper parameters chosen through cross-validation to the file training_prefix.predictor_cv, the means and standard deviations of the input features used for standardizing the input during training to the file training_prefix.X_stdprams, the means and standard deviations of the responses used for standardizing the responses during training to the file training_prefix.Y_stdprams, and the fitted model in TensorFlow format to the file training_prefix.predictor.hdf5.

-------------------------------------------
Generating testing data
-------------------------------------------
Given input single-copy and duplicate gene datasets as described in "Single-copy gene input file format" and "Duplicate gene input file format", respectively, testing data can be generated with the command:

  python run_CLOUD.py -GenerateTestingData 1 -m 6 -S singlecopy_filename -D duplicate_filename -E testing_prefix
  
  eg. python run_CLOUD.py -GenerateTestingData 1 -m 6 -S SingleData -D DuplicateData -E testing
  
where m is the number of conditions, singlecopy_filename is the dataset of single-copy genes described in "Single-copy gene input file format", duplicate_filename is the dataset of duplicate genes described in "Duplicate gene input file format", and testing_prefix is the prefix to all files output by this function.

-------------------------------
Performing test classifications
-------------------------------
Given a trained classifier from the ClassifierCV function as described in "Training the CLOUD classifier" and an input test dataset with features outputted by the GenerateTestingData function described in "Generating testing data", we can perform classification of the five classes of duplicate gene retention mechanisms on the test dataset with the command:

  python run_CLOUD.py -CLOUDClassify 1 -T training_prefix -E testing_prefix
  
  eg. python run_CLOUD.py -CLOUDClassify 1 -T training -E testing
  
where training_prefix is the prefix to training files outputted by ClasifierCV, and testing_prefix is the prefix to the testing files ouputted by GenerateTestingData applied to test data.

**NOTE: The Expression Values for the Test DataSet will be automatically normalized, just enter the original value.

The CLOUDCLassify function outputs predicted classes on each row for each of the duplicate genes in the test dataset to the file test_prefix.classifications.

---------------------------
Performing test predictions
---------------------------
Given a trained predictor from the PredictorCV function as described in "Training the CLOUD predictor" and an input test dataset with features outputted by the GenerateTestingData function described in "Generating testing data", we can perform prediction of the 5m evolutionary model parameters (thetaP, thetaC, thetaA, alpha, and sigmaSq for each of the m conditions) on the test dataset with the command:

  python run_CLOUD.py -CLOUDPredict 1 -T training_prefix -E testing_prefix
  
  eg. python run_CLOUD.py -CLOUDPredict 1 -T training -E testing

where training_prefix is the prefix to training files outputted by PredictorCV, and testing_prefix is the prefix to the testing files outputted by GenerateTestingData applied to test data.

**NOTE: The Expression Values for the Test DataSet will be automatically normalized, just enter the original value.

The CLOUDPredictor function outputs the 5m predicted evolutionary parameters on each row for each of the duplicate genes in the test dataset to the file test_prefix.predictions. 


----------------------------
Example application of CLOUD
----------------------------
Within the R environment, set the working directory to the directory containing both the CLOUD.r script and the subdirectory ExampleFiles containing the example files. 

Load the functions of the CLOUD package by typing the command:
  
  source("CLOUD.r")

Next, generate a training dataset of 100 observations for each of the five classes at six tissues based on the single-copy and duplicate gene data in the ExampleFiles directory, and store the training data with prefix Training, by typing the command:

  GenerateTrainingData(6, 100, "ExampleFiles/SingleData", "ExampleFiles/DuplicateData", "Training")

The above operation will output the files Training.data, Training.features, Training.classes, and Training.responses.
  
To train a CLOUD classifier on this training datset using five-fold cross-validation, with hyper parameters log(lambda) in {-5, -4, -3, -2, -1} and gamma in {0, 0.5, 1}, assuming a batch size of 50 observations per epoch and trained for 50 epochs, type the command:

  ClassifierCV(6, 50, 50, -5, -1, 5, 0, 1, 3, "Training")
 
The above operation will output the files Training.classifier_cv, Training.X_stdparams, and Training.classifier.hdf5.

To train a CLOUD predictor on this training datset using five-fold cross-validation, with hyper parameters log(lambda) in {-5, -4, -3, -2, -1} and gamma in {0, 0.5, 1}, assuming a batch size of 50 observations per epoch and trained for 50 epochs, type the command:

  PredictorCV(6, 50, 50, -5, -1, 5, 0, 1, 3, "Training")
 
The above operation will ouput the files Training.predictor_cv, Training.X_stdparams, Training.Y_stdparams, and Training.predictor.hdf5.
 
Then, generate a test dataset of 100 observations for each of the five classes at six tissues based on the single-copy and duplicate gene data located in the ExampleFiles directory, and store the testing data with prefix Testing, by typing the command:

  GenerateTrainingData(6, 100, "ExampleFiles/SingleData", "ExampleFiles/DuplicateData", "Testing")

The above operation will output the files Testing.data, Testing.features, Testing.classes, and Testing.responses. 

Finally, to perform classification and predictions on these test data, type the commands:
  
  CLOUDClassify("Training", "Testing")
  CLOUDPredict("Training", "Testing")
  
The two above operations will output the files Testing.classifications and Testing.predictions.
