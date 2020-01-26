
1. Input data
The demo input data are provided in the folder 'data'. The subfolder 'goterm_cv' includes data for evulating the performance of IsoResolve with cross validation (cv). The subfolder 'goterm_traintest' includes training data for building models and test data for evaluating the performance of IsoResolve.


2. How to run IsoResolve
We provide two scripts to show how to run IsoResolve
2.1 To test IsoResolve by cross validation, run the following command from command line:
python run_isoresolve_cv.py
This command will excecute CV on the provided data.


2.2 To test IsoResolve by on an independent test dataset, run the following command from command line:
python run_isoresolve_traintest.py
This command will first build a model on the training data and then make predictions on the test data.




