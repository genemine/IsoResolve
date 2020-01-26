import glob,sys,os
from isoresolve_lib import *



# input data
inputdir='data/goterm_cv'
outdir='output_cv/'
if not os.path.exists(outdir):
    os.makedirs(outdir)
fold=5
A=20
lmd=0.01

# MAIN
term=inputdir.split('/')[-1]
ret=isoresolve_cv(inputdir,fold,A=A,lmd=lmd)

outfile=outdir+'/'+term+'.isoresolve.auc.auprc'
iso_scorefile=outdir+'/'+term+'.isoresolve.iso.score'
gene_scorefile=outdir+'/'+term+'.isoresolve.gene.score'
gene_scorefile=outdir+'/'+term+'.isoresolve.gene.score'

f=open(outfile,'w')

f.write('AUC\t'+str(ret['auc'])+'\n')
#f.write('AUC_sig\t'+str(ret['auc_sig'])+'\n')
#f.write('AUC_mig\t'+str(ret['auc_mig'])+'\n')
f.write('AUPRC\t'+str(ret['auprc'])+'\n')
#f.write('AUPRC_sig\t'+str(ret['auprc_sig'])+'\n')
#f.write('AUPRC_mig\t'+str(ret['auprc_mig'])+'\n')
f.write('Aopt\t'+str(ret['Aopt'])+'\n')
f.write('AUC_nlv\t'+str(ret['auc_nlv'])+'\n')
f.close()

print('AUC:'+str(ret['auc']))

isoscore=pd.DataFrame(ret['isoscore'])
isoscore.to_csv(iso_scorefile,index=None,header=None,sep='\t')


genescore=pd.DataFrame(ret['gene_score'])
genescore.to_csv(gene_scorefile,header=None,sep='\t')

