import os,sys,glob
import pandas  as pd
import numpy as np
from numpy import array
import scipy

#from imblearn.over_sampling import SMOTE
from collections import Counter
from sklearn.model_selection import StratifiedKFold
from sklearn import metrics
from sklearn.metrics import precision_recall_curve


#script,inputdir,fold,A,lmd,outdir=sys.argv

#if not os.path.exists(outdir):
#	os.makedirs(outdir)

###############
# sub-routines

def dipls(X,Y,Xs,Xt,A,lmd):
	# X, Y, Xs, Xt are numpy matrix
	#
	varX=np.square(X).sum()
	varY=np.square(Y).sum()
	n=X.shape[0]
	pvar=X.shape[1]
	#
	ns=Xs.shape[0]
	nt=Xt.shape[0]
	#
	T =np.matrix(np.zeros((n,A)))
	Ts=np.matrix(np.zeros((ns,A)))
	Tt=np.matrix(np.zeros((nt,A)))
	P =np.matrix(np.zeros((pvar,A)))
	Ps=np.matrix(np.zeros((pvar,A)))
	Pt=np.matrix(np.zeros((pvar,A)))
	W =np.matrix(np.zeros((pvar,A)))
	q=np.matrix(np.zeros((1,A)))
	B=np.matrix(np.zeros((pvar,A)))

	#
	for i in range(A):
		part1=Xs.T*Xs/(ns-1) 
		part2=Xt.T*Xt/(nt-1)

		y2=Y.T*Y
		y2=y2[0,0]
		coef=lmd/(2*y2)
		Q=np.matrix(np.eye(pvar))+coef*(part1-part2)

		#
		wt= Y.T*X*(Q.I)/y2
		w=wt.T
		w=w/np.linalg.norm(w)
		#
		#
		t=X*w
		ts=Xs*w
		tt=Xt*w
		#
		pt=(t.T*t).I*t.T*X
		p=pt.T
		#
		pst=(ts.T*ts).I*ts.T*Xs
		ps=pst.T
		#
		ptt=(tt.T*tt).I*tt.T*Xt
		pt=ptt.T
		#
		qa=(t.T*t).I*Y.T*t
		qa=qa[0,0]
		X=X-t*p.T
		Xs=Xs-ts*ps.T
		Xt=Xt-tt*pt.T
		Y=Y-qa*t
		# store
		T[:,i] =t	
		Ts[:,i]=ts
		Tt[:,i]=tt
		P[:,i] =p
		Ps[:,i]=ps
		Pt[:,i]=pt
		W[:,i]=w
		q[:,i]=qa

	for k in range(1,A+1):
		Wstar=W[:,0:k]*((P[:,0:k].T*W[:,0:k]).I)
		bk=Wstar*q[:,0:k].T
		B[:,k-1]=bk
	#
	return(B)

def dipls_pred(X,b):
	ypred=X*b
	
	return(ypred)


def find_sig_mig(genelist):
	counter=Counter(genelist)
	allgenes=np.array(list(counter.keys()))
	sizes=np.array(list(counter.values()))
	sig=list(allgenes[sizes==1])
	mig=list(allgenes[sizes>1])
	ret={}
	ret['sig']=sig
	ret['mig']=mig
	return(ret)

def cal_auc_auprc(label,scores):
	fpr, tpr, thresholds = metrics.roc_curve(label, scores)
	precision, recall, thresholds = precision_recall_curve(label,scores)
	auc   = metrics.auc(fpr,tpr)
	auprc = metrics.auc(recall,precision)

	ret={}
	ret['auc']=auc
	ret['auprc']=auprc

	return(ret)


def model_eval(scores,iso2gene,posigene):
	niso=iso2gene.shape[0]
	genes=list(set(iso2gene[:,1]))
	ngenes=len(genes)
	
	# AUC for sig
	tmp=find_sig_mig(list(iso2gene[:,1]))
	sig=tmp['sig']
	mig=tmp['mig']
	
	#scores
	genescores={}
	for i in range(niso):
		iiso=iso2gene[i,0]
		igene=iso2gene[i,1]
		iscore=scores[i]
		if igene in genescores:
			if iscore > genescores[igene]:
				genescores[igene]=iscore
		else:
			genescores[igene]=iscore

	# prepare label and scores for AUC
	label=np.zeros(ngenes,dtype='int')
	scores=np.zeros(ngenes)
	geneclass=np.ones(ngenes)
	combined=np.zeros(shape=(ngenes,3))
	i=0

	rownames=[]
	for igene in genescores.keys():
		scores[i]=genescores[igene]
		if igene in posigene:
			label[i]=1
		if igene in mig:
			geneclass[i]=2

		rownames.append(igene)
		i=i+1
	# combine results
	combined[:,0]=scores
	combined[:,1]=label
	combined[:,2]=geneclass
	ksig=(geneclass==1)
	kmig=(geneclass==2)

	scores_sig=combined[ksig,0]
	label_sig=combined[ksig,1]


	scores_mig=combined[kmig,0]
	label_mig=combined[kmig,1]

	combined=pd.DataFrame(combined,index=rownames)


	# AUC and AUPRC
	res0=cal_auc_auprc(label, scores)
	res1=cal_auc_auprc(label_sig, scores_sig)
	res2=cal_auc_auprc(label_mig, scores_mig)


	ret={}
	ret['auc']      =round(res0['auc'],4)
	ret['auprc']    =round(res0['auprc'],4)
	ret['auc_sig']  =round(res1['auc'],4)
	ret['auprc_sig']=round(res1['auprc'],4)
	ret['auc_mig']  =round(res2['auc'],4)
	ret['auprc_mig']=round(res2['auprc'],4)
	ret['gene_score']   =combined
	return(ret)





#def dipls_paraopt(X,y,Xs,Xt,Xnames,Xsnames,Xtnames,fold=3,Xt_iso2gene,Xt_posigene):
##
#	niso=X.shape[0]
#	isogroup = np.array([s % fold for s in range(niso)])
#	# para ranges
#	AS=[5,8]
#	LMD=[0.001,0.01,0.1,1,10,100]
#	Abest=0
#	lmdbest=0
#	aucbest=0

	

#	for A in AS:
#		for lmd in LMD:
#			print(A,lmd)
#			# pred scores
#			pred=np.zeros((niso,2))
#			for i in range(fold):
#				train_index  = np.where(isogroup!=i)[0]
#				test_index   = np.where(isogroup==i)[0]
#
#				Xtrain=X[train_index,:]
#				ytrain=y[train_index,:]
		
#				Xtest=X[test_index,:]
#				ytest=y[test_index,:]

				# find corresponding Xs and Xt
#				trainnames=Xnames[train_index]
#				ksels=~pd.Series(Xsnames).isin(trainnames)
#				kselt=~pd.Series(Xtnames).isin(trainnames)
#				Xstrain=Xs[ksels,:]
#				Xttrain=Xt[kselt,:]

#				# build models
#				beta=dipls(Xtrain,ytrain,Xstrain,Xttrain,A,lmd)
				# make predictions
#				fold_scores=dipls_pred(Xtest,beta)
				# store
#				pred[test_index,:]=np.column_stack((ytest,fold_scores))

			# AUC
            # result=model_eval(predscore,iso2gene,posigene)
            
            
            
#			tmp=cal_auc_auprc(pred[sig_index:,0],pred[sig_index:,1])
#			thisauc=tmp['auc']
			
#			if thisauc > aucbest:
#				Abest=A
#				lmdbest=lmd
#				aucbest=thisauc


#	return(Abest,lmdbest,aucbest)
    
    
    


def partitionisoform4cv(iso2gene,fold=3):
    niso=iso2gene.shape[0]
    ugenes=list(set(iso2gene[:,1]))
    ngenes=len(ugenes)
    ggroup = np.array([s % fold for s in range(ngenes)])
    gene2group=np.column_stack((ugenes,ggroup))
    
    
    iso2group=np.zeros((niso,1))
    for i in range(niso):
        igene=iso2gene[i,1]
        ki=np.where(gene2group[:,0]==igene)[0][0]
        igroup=gene2group[ki,1]
        iso2group[i,0]=float(igroup)
        
    iso2group=iso2group[:,0]
    
    return(iso2group)
    


def dipls_paraoptimize(X,y,Xs,Xt,Xnames,Xsnames,Xtnames,iso2gene,posigene,fold=3):
	niso=Xt.shape[0]
	#isogroup = np.array([s % fold for s in range(niso)])
	isogroup=partitionisoform4cv(iso2gene,fold=fold)
	# para ranges
	AS=[1,3,5,7,10]
	LMD=[0.001,0.01,0.1,1,10,100]
	Abest=0
	lmdbest=0
	aucbest=0

	for A in AS:
		for lmd in LMD:
			#print(A,lmd)
			# pred scores
			pred=np.zeros((niso,1))
			for i in range(fold):
				train_index  = np.where(isogroup!=i)[0]
				test_index   = np.where(isogroup==i)[0]
				iso2gene_train=iso2gene[train_index,:]
				iso2gene_test=iso2gene[test_index,:]
				
				traingenes=list(set(iso2gene_train[:,1]))
				testgenes =list(set(iso2gene_test[:,1]))
				# target domain
				Xt_train_i=Xt[train_index,:]
				Xt_test_i =Xt[test_index,:]
                
                		# source domain
				ksels=pd.Series(Xsnames).isin(traingenes)
				Xsi=Xs[ksels,:]

                		# combined domain
				ksel_Xy=~pd.Series(Xnames).isin(testgenes)
				Xi=X[ksel_Xy,:]
				yi=y[ksel_Xy,:]
				
				# build models
				beta=dipls(Xi,yi,Xsi,Xt_train_i,A,lmd)
				# make predictions
				fold_scores=dipls_pred(Xt_test_i,beta)
				# store
				pred[test_index,:]=fold_scores

			# AUC
			tmp_result=model_eval(pred,iso2gene,posigene)
			#tmp=cal_auc_auprc(pred[sig_index:,0],pred[sig_index:,1])
			thisauc=tmp_result['auc']

			if thisauc > aucbest:
				Abest=A
				lmdbest=lmd
				aucbest=thisauc


	return(Abest,lmdbest,aucbest)





def isoresolve_cv(inputdir,fold,A=15,lmd=5,optimize=False):
	posigene=[]
	predscore=np.zeros((0,A))
	iso2gene=np.zeros((0,2),dtype=object)
	for i in range(fold):
		trainisofile  =inputdir+'/train_iso_'+ str(i)+'.tsv'
		traingenefile =inputdir+'/train_gene_'+str(i)+'.tsv'
		trainlabelfile=inputdir+'/train_label_'+str(i)+'.label'

		testisofile  =inputdir+'/test_iso_'+ str(i)+'.tsv'
		testgenefile =inputdir+'/test_gene_'+str(i)+'.tsv'
		testlabelfile=inputdir+'/test_label_'+str(i)+'.label'

		train=read_data(trainisofile,traingenefile,trainlabelfile)
		test =read_data(testisofile,testgenefile,testlabelfile)

		# build a model
		Xs=train['Xs']         # source samples
		ys=train['ys']
		Xt=train['Xt']         # target samples
		Xtl=train['Xtl']       # labelled target samples
		ytl=train['ytl']
		Xsnames=train['Xsnames']
		Xtnames=train['Xtnames']
		Xtlnames=train['Xtlnames']
		train_iso2gene=train['iso2gene']
		train_posigene=train['posigene']
        
        
		X=np.vstack((Xs,Xtl))  # supervised manner
		y=np.vstack((ys,ytl))  # supervised manner
		Xnames=np.concatenate((Xsnames,Xtlnames))

		#if optimize:
		#	A,lmd,aucbest=dipls_paraoptimize(X,y,Xs,Xt,Xnames,Xsnames,Xtnames,train_iso2gene,train_posigene,fold)
		#	print('Fold ',i,':',A,lmd,aucbest)

		beta=dipls(X,y,Xs,Xt,A,lmd)

		# make predictions
		Xtest=test['Xt']
		thisiso2gene=test['iso2gene']
		thisposigene=test['posigene']

		# make predictions
		ypred=dipls_pred(Xtest,beta)
		
		# store
		predscore=np.vstack((predscore,ypred))
	
	
		iso2gene =np.vstack((iso2gene,thisiso2gene))
		posigene=list(set(posigene+thisposigene))
		

	# AUC
	result={}
	result['auc']=0
	auc_nlv=[]
	for i in range(A):
		result_i=model_eval(predscore[:,i],iso2gene,posigene)
		pscore_i=np.column_stack((iso2gene,predscore[:,i]))
		result_i['isoscore']=pscore_i
		auc_nlv.append(result_i['auc'])
		if result_i['auc'] > result['auc']:
			result=result_i
			result['Aopt']=i+1		
		
	result['auc_nlv']=auc_nlv

	return(result)



def isoresolve_train_test(inputdir,fold,A=15,lmd=5):
    posigene=[]
    predscore=np.zeros((0,A))
    iso2gene=np.zeros((0,2),dtype=object)
    #for i in range(fold):
    trainisofile  =inputdir+'/train_iso.tsv'
    traingenefile =inputdir+'/train_gene.tsv'
    trainlabelfile=inputdir+'/train_label.label'

    testisofile  =inputdir+'/test_iso.tsv'
    testgenefile =inputdir+'/test_gene.tsv'
    testlabelfile=inputdir+'/test_label.label'

    train=read_data(trainisofile,traingenefile,trainlabelfile)
    test =read_data(testisofile,testgenefile,testlabelfile)

    # build a model
    Xs=train['Xs']         # source samples
    ys=train['ys']
    Xt=train['Xt']         # target samples
    Xtl=train['Xtl']       # labelled target samples
    ytl=train['ytl']
    Xsnames=train['Xsnames']
    Xtnames=train['Xtnames']
    Xtlnames=train['Xtlnames']
    train_iso2gene=train['iso2gene']
    train_posigene=train['posigene']
        
        
    X=np.vstack((Xs,Xtl))  # supervised manner
    y=np.vstack((ys,ytl))  # supervised manner
    Xnames=np.concatenate((Xsnames,Xtlnames))

    beta=dipls(X,y,Xs,Xt,A,lmd)

    # make predictions
    Xtest=test['Xt']
    thisiso2gene=test['iso2gene']
    thisposigene=test['posigene']

    # make predictions
    ypred=dipls_pred(Xtest,beta)
		
    # store
    predscore=np.vstack((predscore,ypred))
	
	
    iso2gene =np.vstack((iso2gene,thisiso2gene))
    posigene=list(set(posigene+thisposigene))
		

    # AUC
    result={}
    result['auc']=0
    auc_nlv=[]
    for i in range(A):
        result_i=model_eval(predscore[:,i],iso2gene,posigene)
        pscore_i=np.column_stack((iso2gene,predscore[:,i]))
        result_i['isoscore']=pscore_i
        auc_nlv.append(result_i['auc'])
        if result_i['auc'] > result['auc']:
            result=result_i
            result['Aopt']=i+1
    result['auc_nlv']=auc_nlv
    return(result)

def read_data(isofile,genefile,labelfile):
	data=pd.read_csv(isofile,header=None,sep='\t')
	Xt=np.matrix(data.iloc[:,2:])


def read_data(isofile,genefile,labelfile):
	data=pd.read_csv(isofile,header=None,sep='\t')
	Xt=np.matrix(data.iloc[:,2:])
	Xtnames=np.array(data.values[:,0])
	iso2gene=data.values[:,0:2]
	tgenes=data.values[:,1]
	nt=Xt.shape[0]



	# read source domain data
	data=pd.read_csv(genefile,header=None,sep='\t')
	Xs=np.matrix(data.iloc[:,1:])
	sgenes=data.values[:,0]
	Xsnames=sgenes
	ns=Xs.shape[0]
	nsample=Xs.shape[1]

	# get source response
	gl=pd.read_csv(labelfile,header=None,sep='\t').values
	ys=np.matrix(np.zeros((ns,1)))
	i=0
	for igene in sgenes:
		k=np.where(gl[:,0]==igene)
		ys[i,0]=gl[k,2]
		i=i+1
	
	# find posi genes
	kposi=np.where(gl[:,2]==1)
	posigene=list(np.unique(gl[kposi][:,0]))
	# get sig label: labled target samples
	ret=find_sig_mig(iso2gene[:,1])
	sig=ret['sig']
	nsig=len(sig)
	Xtl=np.zeros((nsig,nsample))
	Xtlnames=np.array(sig)
	ytl=np.zeros((nsig,1))
	ksig=[]
	j=0
	for isig in sig:
		k=np.where(iso2gene[:,1]==isig)
		k=k[0]
		k=k[0]
		Xtl[j,:]=Xt[k,:]

		k1=np.where(gl[:,0]==isig)
		k1=k1[0]
		k1=k1[0]
		ytl[j,0]=gl[k1,2]
		j=j+1


	ret={}

	Xs_expand=np.column_stack((Xs,np.ones((Xs.shape[0],1))))
	Xt_expand=np.column_stack((Xt,np.ones((Xt.shape[0],1))))
	Xtl_expand=np.column_stack((Xtl,np.ones((Xtl.shape[0],1))))
	
	ret['Xs']=Xs_expand
	ret['ys']=ys
	ret['Xt']=Xt_expand
	ret['Xtl']=Xtl_expand
	ret['ytl']=ytl
	ret['Xsnames']=Xsnames
	ret['Xtnames']=Xtnames
	ret['Xtlnames']=Xtlnames
	ret['posigene']=posigene
	ret['iso2gene']=iso2gene

	return(ret)
# end of subroutines

# main
# basic settings: users can modify paras here

##############################
###### NO changes below ######
##############################
