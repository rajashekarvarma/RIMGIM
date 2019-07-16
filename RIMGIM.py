# from rpy2.robjects.packages import importr

# # do the following _only the first time_, to install the package seqLogo
# base = importr('base')
# base.source("http://www.bioconductor.org/biocLite.R")
# biocinstaller = importr("BiocInstaller")
# biocinstaller.biocLite("ALL")


# from rpy2.robjects import r
# r_src = """
# source("http://www.bioconductor.org/biocLite.R")
# biocLite("limma")
# """
# r(r_src)

import GEOparse # Python package to upload a geo data
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats import multitest
import seaborn as sns
import matplotlib.pyplot as plt


############### Fetch Agilent data ###############
print('\n\n',"******...Hi Welcome to RIPGEO...******",'\n\n')
GSE_ID = input('Please enter your GSE ID (ex:GSE62893): ')

print('\n',"Provided GSE ID: ",GSE_ID)

print('\n',"Intitating data extraction...",'\n\n')

gse = GEOparse.get_GEO(geo=GSE_ID, destdir="./")
plt_name=[]
# print(gse.gpls)
for pl_name, pl in gse.gpls.items():
    plt_name.append(pl_name)
plt_name=''.join(plt_name)

print("Platform Name:", plt_name)

pivoted_control_samples = gse.pivot_samples('VALUE')
# print(pivoted_control_samples.head())

######## Filter probes that are not expressed worst 25% genes are filtered out

pivoted_control_samples_average = pivoted_control_samples.median(axis=1)
# print("Number of probes before filtering: ", len(pivoted_control_samples_average))
expression_threshold = pivoted_control_samples_average.quantile(0.25)
expressed_probes = pivoted_control_samples_average[pivoted_control_samples_average >= expression_threshold].index.tolist()
# print("Number of probes above threshold: ", len(expressed_probes))

samples = gse.pivot_samples("VALUE").ix[expressed_probes]
# print(samples.head())
print(gse.gpls[plt_name].table.head())

######## Annotate matrix table

samples_annotated = samples.reset_index().merge(gse.gpls[plt_name].table, left_on='ID_REF', right_on="ID").set_index('ID_REF')

# print(samples_annotated.head())
del samples_annotated["ID"]
# print(samples_annotated.head())
samples_annotated = samples_annotated.dropna(subset=["GB_ACC"])
samples_annotated = samples_annotated[~samples_annotated.GB_ACC.str.contains("///")]
samples_annotated = samples_annotated.groupby("GB_ACC").median()
# print(samples_annotated.index)

print('\n','Column names from the matrix: ',samples_annotated.columns)

######## Extract matrix data to a csv file
exprs = []
gsmNames = []
metadata = {}

for gsm_name, gsm in gse.gsms.items():
    # print(gsm.metadata['type'][0])
    if gsm.metadata['type'][0]=='RNA':
        # Expression data
        if len(gsm.table)>0:
            tmp = gsm.table['VALUE']
            # print(tmp)
            tmp.index = gsm.table['ID_REF']
            gsmNames.append(gsm_name)
            if len(exprs)==0:
                exprs = tmp.to_frame()
            else:
                exprs = pd.concat([exprs,tmp.to_frame()],axis=1)

print('\n','Extracting metadata...','\n')

######## extract metadata to csv file

for gsm_name, gsm in gse.gsms.items():
    if gsm.metadata['type'][0]=='RNA':
                for key,value in gsm.metadata.items():
                # print(key)
                # print(value)
                    if (key=='characteristics_ch1' or key=='characteristics_ch2') and (len([i for i in value if i!=''])>1 or value[0].find(': ')!=-1):
                        # print(value)
                        tmpVal = 0
                        for tmp in value:
                            splitUp = [i.strip() for i in tmp.split(':')]
                            # print(splitUp)
                            if len(splitUp)==2:
                                if not splitUp[0] in metadata:
                                    metadata[splitUp[0]] = {}
                                metadata[splitUp[0]][gsm_name] = splitUp[1]
                            else:
                                if not key in metadata:
                                    metadata[key] = {}
                                metadata[key][gsm_name] = splitUp[0]
                    else:
                        if not key in metadata:
                            metadata[key] = {}
                        if len(value)==1:
                            metadata[key][gsm_name] = ' '.join([j.replace(',',' ') for j in value])

# Write expression data matrix to file
exprs.columns = gsmNames
with open(GSE_ID+'exprs.csv','w') as outFile:
    exprs.to_csv(outFile)

# Write metadata matrix to file
with open(GSE_ID+'metadata.csv','w') as outFile:
    outFile.write('Metadata,'+','.join(gsmNames))
    for key in metadata:
        tmp = [key]
        for gsm_name in gsmNames:
            if gsm_name in metadata[key]:
                tmp.append(metadata[key][gsm_name])
            else:
                tmp.append('NA')
        outFile.write('\n'+','.join(tmp))
print('\n','Data matrix and metadata for',GSE_ID,'have been written to',GSE_ID+'exprs.csv',GSE_ID+'metadata.csv @ cwd','\n')

######## select control and test sample columns

samples = samples_annotated.astype(float)

######################## use this for illumina

# import pandas as pd
# import glob
# import numpy as np
# from scipy import stats
# from statsmodels.stats import multitest
# import seaborn as sns
# import matplotlib.pyplot as plt

# c_fo = input('Please input control sample files folder name(ex:control):')
# t_fo = input('Please input test sample files folder name(ex:test):')

# cpath ="*.gz"
# # tpath =t_fo+"/*.gz"
# samples = pd.DataFrame()
# # samples_annotated = pd.DataFrame()
# for i, filename in enumerate(glob.glob(cpath)):
#     # print(filename)
#     colV = filename.strip('.txt.gz')
#     # colV=colV.replace('ill2\\','')
#     samples2 = pd.DataFrame()

#     df = pd.read_csv(filename, compression='gzip',sep='\t')
#     if 'GeneID' not in samples.columns:
#         samples = df[['GeneID','RPKM']]
#         samples.columns=['GeneID',colV]
#         # print(samples.head())
#     else:
#         samples2 = df[['GeneID','RPKM']]
#         samples2.columns=['GeneID',colV]
#         # print(samples2.head())

#         samples = pd.merge(samples,samples2, on='GeneID', )

# print("\n",'List of file names read to dataframe:','\n')
# samples.set_index("GeneID", inplace = True) 
# print(samples.head())
########################################################################
samples.to_csv('matrix.csv',)
    # print(colV)
samples_annotated = samples
import rpy2.robjects as robjects 
from rpy2.robjects.packages import importr 
ALL = importr('ALL') 
limma = importr('limma')
exprs = robjects.r['exprs']
summary = robjects.r['summary']
matrix = robjects.r['as.matrix']
new = robjects.r['new']
robjects.r('data("ALL")')
data = robjects.globalenv['ALL'] 
featureNames = robjects.r['featureNames'] 
ExpressionSet = robjects.r['ExpressionSet']
character = robjects.r['as.character']
pas = robjects.r['paste']
fac = robjects.r['as.factor']
mmax = robjects.r['model.matrix']





from rpy2.robjects import pandas2ri
pandas2ri.activate()

#print(samples.head())
samples_log2=samples.apply(np.log2)
# print(samples_log2.head())

import rpy2.robjects as ro

# print(asda)
R = ro.r

control_sample = input('Please enter column numbers of control samples (ex:0,2,4): ')
control_samples = control_sample.split(',')
control_samples = [int(i) for i in control_samples]
Test_sample =  input('Please enter column numbers of test samples (ex:3,5,7): ')
Test_samples = Test_sample.split(',')
Test_samples = [int(i) for i in Test_samples]

f_samples = pd.DataFrame() 
f_samples2 = samples_log2.iloc[:,Test_samples] 
f_samples1 = samples_log2.iloc[:,control_samples]
f_samples = f_samples1.join(f_samples2)


r_samples = pandas2ri.py2ri(f_samples)

a_labels=control_samples + Test_samples
# print(a_labels)
a_lab=list()
for i, t in enumerate(a_labels):
    if i < len(control_samples):
        a_lab.append("1")
    else:
        a_lab.append("0")
sml = pas("G",a_lab,sep="")
fl = fac(sml)
# print(fl)
R.assign('fl',fl)
R.assign('r_samples',r_samples)
R('design<-model.matrix(~ fl + 0,r_samples)')
R('print(design)')
R('colnames(design)<-levels(fl)')
R('print(design)')
R('fit<-lmFit(r_samples,design)')
R('cont.matrix <- makeContrasts(G1-G0, levels=design)')
R('fit2 <- contrasts.fit(fit, cont.matrix)')
R('fit2 <- eBayes(fit2, 0.01)')
R('tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)')
R('tT <- subset(tT, select=c("adj.P.Val","P.Value","t","B","logFC"))')
# R('print(tT)')

R('write.csv(tT, file="F_matrix.csv", sep=",")')
