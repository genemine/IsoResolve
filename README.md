# IsoResolve


<img width="450" height="150" src="https://github.com/genemine/iread/blob/master/images/ir1.png"/>


## 1 Description
IsoResolve is a computational approach for isoform function prediction by leveraging the information of gene function prediction models with domain adaptation (DA). IsoResolve treats gene- level and isoform-level features as source and target domain, respectively. It employs DA to project the two domains to a latent variable (LV) space in such a way that the LVs projected from the gene and isoform domain features are of approximately the same distribution, enabling that the gene domain information can be leveraged for predicting isoform functions.


## 2 Input data
The demo input data are provided in the folder 'data'. The subfolder 'goterm_cv' includes data for evulating the performance of IsoResolve with cross validation (cv). The subfolder 'goterm_traintest' includes training data for building models and test data for evaluating the performance of IsoResolve.

## 3 Download

* **software**

The source codes are tested and work on both MacOS and Linux operating system. They are freely available for non-commercial use.<br>

| **Version** | **Changes** |
| - | - |
| [iREAD_0.5.0.zip](https://github.com/genemine/iread/raw/master/history_version/iREAD_0.5.0.zip) |  |

## 4. Usage
We provide two scripts to show how to run IsoResolve
### 4.1 To test IsoResolve by cross validation, run the following command from command line:
python run_isoresolve_cv.py
This command will excecute CV on the provided data.


### 4.2 To test IsoResolve by on an independent test dataset, run the following command from command line:
python run_isoresolve_traintest.py
This command will first build a model on the training data and then make predictions on the test data.

# 5. Contact
If any questions, please do not hesitate to contact me at:
<br>
Hongdong Li `hongdong@csu.edu.cn`
<br>
Jianxin Wang `jxwang@csu.edu.cn`


