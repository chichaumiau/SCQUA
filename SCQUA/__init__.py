#  spikein_utils.py
#  Single Cell Sequencing Quality Assessment: scqua
#  
#  Copyright 2018 Chichau Miau <zmiao@ebi.ac.uk>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

import pandas as pd
from glob import iglob
import click
from sklearn.linear_model import LogisticRegression
import numpy as np
import scipy
import seaborn as sns
import statsmodels.api as sm
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt
import matplotlib

def get_ERCC():
	ercc = pd.read_table('https://raw.githubusercontent.com/Teichlab/readquant/master/readquant/ERCC.tsv', index_col=1)
	ercc = np.log(ercc['concentration in Mix 1 (attomoles/ul)'])
	return(ercc)

def get_SIRV():
	sirv = pd.read_csv('https://raw.githubusercontent.com/chichaumiau/scRNA_metadata/master/SIRV_E2/SIRV_concentration.csv', index_col=1)
	sirv = np.log(sirv['E2 molarity [fmoles/µl]']*1000)
	return(sirv)

def get_detection_limit(spike, quant, det_threshold=0.1):
	X = spike[:, None]
	y = quant[spike.index] >= det_threshold
	if y.sum() < 8:
		return np.inf
	lr = LogisticRegression(solver='liblinear', fit_intercept=True)
	lr.fit(X, y)
	midpoint = -lr.intercept_ / lr.coef_[0]
	return np.exp(midpoint[0])

def get_accuracy(ercc, quant, det_threshold=0.1):
	y = np.log(quant[ercc.index]) \
			.replace([np.inf, -np.inf], np.nan) \
			.dropna()
	if (y >= np.log(det_threshold)).sum() < 8:
		return -np.inf
	correlation = y.corr(ercc, method='pearson')
	return correlation

def get_phn(cts_file,tpm_file,phn_file, ercc, sirv, spike):
	cts = pd.read_csv(cts_file, index_col=0)
	tpm = pd.read_csv(tpm_file, index_col=0)
	phn = pd.read_csv(phn_file, index_col=0)
	df = get_result(tpm, ercc, sirv, spike)
	phn = pd.concat([phn,df,cts.loc[cts.index.str.startswith("ENS")].T], axis=1)
	phn["Total_counts"] = cts.loc[cts.index.str.startswith("ENS")].sum()
	return(phn)

def get_result(tpm, ercc=None, sirv=None, spike=None):
	df = pd.DataFrame()
	for col in tpm.columns:
		quant = tpm[col]
	
		qc_data = pd.Series()
		if not spike is None:
			qc_data['detection_limit'] = get_detection_limit(spike, quant)
			qc_data['accuracy'] = get_accuracy(spike, quant)
		
		if not ercc is None:
			qc_data['detection_limit_ERCC'] = get_detection_limit(ercc, quant)
			qc_data['accuracy_ERCC'] = get_accuracy(ercc, quant)
		
		if not sirv is None:
			qc_data['detection_limit_SIRV'] = get_detection_limit(sirv, quant)
			qc_data['accuracy_SIRV'] = get_accuracy(sirv, quant)
	
		try:
			qc_data['ERCC_content'] = quant[ercc.index].sum()
			quant = quant.drop(ercc.index)
			quant = quant / quant.sum() * 1e6
		except ValueError:
			# ERCCs not present
			pass
	
		try:
			qc_data['SIRV_content'] = quant[sirv.index].sum()
			quant = quant.drop(sirv.index)
			quant = quant / quant.sum() * 1e6
		except ValueError:
			# ERCCs not present
			pass
	
		qc_data['n_genes'] = (quant.loc[quant.index.str.startswith("ENS")] > 1.).sum()
		df[col] = qc_data
	
	return(df.T)

ercc = get_ERCC()
sirv = get_SIRV()
spike = pd.concat([ercc,sirv])

def get_result_ad(ad, ercc=None, sirv=None, spike=None):
	df = pd.DataFrame()
	for id in ad.obs_names:
		quant = pd.Series(ad[id,:].X, index=ad.var_names)
		qc_data = pd.Series()
		if not spike is None:
			qc_data['detection_limit'] = get_detection_limit(spike, quant)
			qc_data['accuracy'] = get_accuracy(spike, quant)
		
		if not ercc is None:
			qc_data['detection_limit_ERCC'] = get_detection_limit(ercc, quant)
			qc_data['accuracy_ERCC'] = get_accuracy(ercc, quant)
	
		try:
			qc_data['ERCC_content'] = quant[ercc.index].sum()
			quant = quant.drop(ercc[ercc.index.isin(quant.index)].index)
			quant = quant / quant.sum() * 1e6
		except ValueError:
			# ERCCs not present
			pass
		df[id] = qc_data
	return(df.T)

def plot_Fig1D(df1, df2, \
			   protocol = "SMARTer", \
				key2 = 'Protocol2', \
				key3 = 'BGISEQ-500', \
				key4 = 'Fullname',\
				xlabelsize = None, ylabelsize = None, titlesize = None):
	from sklearn.decomposition import PCA
	from sklearn.preprocessing import scale	
	
	df = pd.concat([df1,df2], axis = 0)
	exprs = np.log1p(df.T.loc[df.columns.str.startswith("ENS")].astype(np.float32))
	pca = PCA(n_components=50)
	pca_res = pca.fit_transform(scale(exprs.T, 1))

	df["PC1"] = pca_res[:,0]
	df["PC2"] = pca_res[:,1]

	g = df.groupby('Protocol2')
	colors = ["#00cc99","#e0115f"]
	for i, group in enumerate(g.indices):
		tmp = g.indices[group]
		plt.scatter(df.PC1.iloc[tmp], df.PC2.iloc[tmp], label=group.replace(protocol,""), lw = 0, c=colors[i])
	for i in df[df[key2].str.startswith(key3)][key4]:
		j = df[(df[key2].str.startswith(key3)) & (df[key4] == i)][['PC1','PC2']]
		k = df[(~df[key2].str.startswith(key3)) & (df[key4] == i)][['PC1','PC2']]
		plt.plot([j["PC1"],k["PC1"]],[j["PC2"],k["PC2"]],'k-.')
	
	plt.xlabel("PC1 (%.3f)"%pca.explained_variance_ratio_[0], fontsize = xlabelsize)
	plt.ylabel("PC2 (%.3f)"%pca.explained_variance_ratio_[1], fontsize = ylabelsize)
	if not titlesize is None:
		plt.title(protocol, fontsize = titlesize)
	else:
		plt.title(protocol)
	plt.legend(scatterpoints=3, bbox_to_anchor=(1.03, 1), borderaxespad=0.);

def plot_Fig1E(df, \
	protocol = "SMARTer",\
	key1 = 'Batch',\
	key2 = 'Protocol2',\
	key3 = 'Cell',\
	key4 = 'Fullname',\
	lim = 210,\
	cutoff = 1000, \
	xlabel = None, ylabel=None,\
	xlabelsize = None, ylabelsize=None, titlesize = None, legendsize = None):

	dfx = df[df[key1] == protocol]
	dfx = dfx[dfx.index.str.find('e6')>0]
	df1 = dfx[dfx[key2].str.startswith('BGISEQ-500')].sort_values(by=key3)
	df2 = dfx[dfx[key2].str.startswith('HiSeq-2000')].sort_values(by=key3)
	df1 = df1[df1[key4].isin(df2[key4])]
	df2 = df2[df2[key4].isin(df1[key4])]

	xx1 = df1.detection_limit[(df2.detection_limit <cutoff).tolist()]
	xx2 = df2.detection_limit[(df2.detection_limit <cutoff).tolist()]

	slope, intercept, r_value, p_value, std_err = \
		scipy.stats.linregress(xx1, xx2)
	plt.plot(xx1, xx2, 'k.', label =  "R=%.2f"%r_value)

	plt.xlim(0,lim)
	plt.ylim(0,lim)
	
	if xlabel is None: xlabel = key2
	if ylabel is None: ylabel = key1
	if xlabelsize is None:
		plt.xlabel(xlabel)
	else:
		plt.xlabel(xlabel, fontsize=xlabelsize)
	if ylabelsize is None:
		plt.ylabel(ylabel)
	else:
		plt.ylabel(ylabel, fontsize=ylabelsize)
	if not titlesize is None:
		plt.title(protocol, fontsize = titlesize)
	else:
		plt.title(protocol)
	if legendsize is None:
		plt.legend(loc='upper right')
	else:
		plt.legend(loc='upper right', fontsize=legendsize)
	return(df1, df2)
	
def fit_sensitivity(df, fun = 'np.log10(detection_limit)',key1 = "detection_limit", key2 = "n_counts", key3 = "protocol", \
					names = None, xlabel = None, ylabel = None,\
					xscale = 'log', yscale=None, \
					xlim = None, ylim = None, \
					xlabelsize = None, ylabelsize = None, \
					title = None, titlesize = None, \
					colors = None, colordots = False, \
					save = None):
	formula = '%s ~ np.power(np.log10(%s), 2) + np.log10(%s) + C(%s) + 1'%(fun,key2,key2,key3)
	print(formula)
	mod = smf.ols(formula=formula, data=df)
	res = mod.fit()
	print(res.summary())
	print(res.params.sort_values())
	read_satureation = 10 ** (res.params['np.log10(%s)'%key2] / (-2 * res.params['np.power(np.log10(%s), 2)'%key2]))
	print("read_satureation: %d"%read_satureation)
	
	normalised_sensitivity = pd.Series()
	for protocol in df[key3].unique():
		normalised_sensitivity[protocol] = 10 ** res.predict(pd.DataFrame({key2: [1e6], key3: [protocol]}))

	if names is None:
		names = normalised_sensitivity.index
	if colors is None:
		colors = [matplotlib.colors.rgb2hex(plt.get_cmap("Paired")(i)) for i in range(len(names))]
	
	xx=np.logspace(3, np.log10(read_satureation), 50)
	pdf = pd.DataFrame({key2: xx})
	
	fig, ax = plt.subplots()
	if not xscale is None:
		plt.xscale(xscale)
	if not yscale is None:
		plt.yscale(yscale)
	
	plt.scatter(df[key2], df[key1], c='#BBBBBB', edgecolor='none', s=50, label=None, rasterized=True);

	for i, protocol in enumerate(names):
		pdf[key3] = protocol
		if fun.startswith("np.log"):
			yy = 10 ** res.predict(pdf)
		else:
			yy = res.predict(pdf)
		xx2 = np.concatenate((xx, np.array([1e8])))
		yy2 = np.append(yy,yy.iloc[-1])
		plt.plot(xx2, yy2, label=protocol, color = colors[i])
		if colordots:
			plt.scatter(df.loc[df[key3] == protocol][key2], df.loc[df[key3] == protocol][key1], c=colors[i], edgecolor='none', s=50, label=None, rasterized=True);
		
	plt.axvline(xx[-1], linestyle='--', c='r', lw=1)
	for x in [1e4, 1e5, 1e6]:
		plt.axvline(x, linestyle='--', c='grey', lw=1)

	lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5));

	sns.despine()
	if xlabel is None: xlabel = key2
	if ylabel is None: ylabel = key1
	if xlabelsize is None:
		plt.xlabel(xlabel)
	else:
		plt.xlabel(xlabel, fontsize=xlabelsize)
	if ylabelsize is None:
		plt.ylabel(ylabel)
	else:
		plt.ylabel(ylabel, fontsize=ylabelsize)
	
	if not xlim is None:
		plt.xlim(xlim)
	if not ylim is None:
		plt.ylim(ylim)
	if not title is None:
		if titlesize is None:
			plt.title(title)
		else:
			plt.title(title, fontsize = titlesize)
	if not save is None:
		plt.savefig(save)
	
def plot_jitter(df, key1 = "detection_limit", key2 = "lane", ylog = False, xlabel = None, ylabel = None):
	ax = sns.violinplot(x=key2,y=key1,data=df,inner=None, color=".9", cut=0)
	if ylog: ax.set_yscale("log", nonposy='clip')
	ax = sns.stripplot(x=key2,y=key1,data=df, jitter=0.3)
	ax.set_xticklabels(ax.get_xticklabels(),rotation=90)
	if xlabel is None: xlabel = key2
	if ylabel is None: ylabel = key1
	plt.xlabel(xlabel)
	plt.ylabel(ylabel);
	return(ax)

	
import GPy
import sys
import struct
import gzip

import sys,os
from glob import iglob
import numpy as np
import pandas as pd
from tqdm import tqdm


def read_kallisto(sample_path):
    ''' Function for reading a Kallisto quantification result.
    Returns
    -------
    A pandas.Series with the expression values in the sample.
    '''
    quant_file = sample_path + '/abundance.tsv'
    df = pd.read_table(quant_file, engine='c',
                                   usecols=['target_id', 'tpm','est_counts'],
                                   index_col=0,
                                   dtype={'target_id': np.str, 'tpm': np.float64, 'est_counts':np.float64 })

    df = df.rename(columns={'tpm': 'TPM', 'est_counts':'counts','target_id':'Name'})
    return df



def read_salmon(sample_path, tp ='gene'):
    ''' Function for reading a Salmon quantification result.

    Parameters
    ----------
    isoforms : bool, default False
        Whether to parse isoform level expression or gene level expression.

    Returns
    -------
    A pandas.Series with the expression values in the sample.
    '''
    read_kwargs = {
        'engine': 'c',
        'usecols': ['Name', 'TPM','NumReads'],
        'index_col': 0,
        'dtype': {'Name': np.str, 'TPM': np.float64, 'NumReads': np.float64}
    }
    if tp == 'gene':
        df = pd.read_table('%s/quant.genes.sf'%sample_path, **read_kwargs)
    else:
        df = pd.read_table('%s/quant.sf'%sample_path, **read_kwargs)
    df = df.rename(columns={'NumReads': 'counts'})
    return df

def read_salmon_qc(sample_path, flen_lim=(100, 100)):
    ''' Parse technical quality control data from a Salmon quantification
    result.

    Parameters
    ----------
    flen_lim, tuple (int start, int end), default (100, 100)
        How many bases to remove from start and end of fragment length
        distribution when calculating the robust mode. This is too see if
        things roughly worked out even if the max FLD Salmon parameter was
        set too small.

    Returns
    -------
    A pandas.Series with technical information from the Salmon results for
    the sample.
    '''
    flen_dist = np.fromfile('%s/libParams/flenDist.txt'%sample_path, sep='\t')
    global_fl_mode = flen_dist.argmax()
    robust_fl_mode = flen_dist[flen_lim[0]:-flen_lim[1]].argmax() + flen_lim[0]

    qc_data = pd.read_json('%s/aux_info/meta_info.json'%sample_path, typ='series')  
    qc_data = qc_data[['num_processed', 'num_mapped', 'percent_mapped']]
    qc_data['global_fl_mode'] = global_fl_mode
    qc_data['robust_fl_mode'] = robust_fl_mode

    return qc_data

def read_quants(pattern='salmon/*_salmon_out', tool='salmon', **kwargs):
    ''' Read quantification results from every directory matching the glob
    in pattern.
    Parameters
    ----------
    tool, str, default 'salmon'
        The quantification tool used to generate the results. Currently
        supports 'salmon', 'sailfish', 'kallisto', and 'cufflinks'.
    **kwargs,
        kwargs are passed on to the tool specific sample parser. See documentation
        for individual parsers for details.
    Returns
    -------
    A pandas.DataFrame where columns are samples, rows are genes, and cells
    contain the expression value.
    '''
    sample_readers = {
        'salmon': read_salmon,
	'kallisto':read_kallisto
    }

    quant_reader = sample_readers[tool]

    tpm_quants = pd.DataFrame()
    cts_quants = pd.DataFrame()
    for sample_path in tqdm(iglob(pattern)):
        if os.path.isdir(sample_path):
            #print sample_path
            if not os.path.isfile('%s/quant.genes.sf'%sample_path):
                print("Not found %s"%sample_path)
                #os.system('rm %s -rf'%sample_path)
                continue
            sample_quant = quant_reader(sample_path, **kwargs)
            
            id=sample_path.split('/')[-1]
            tpm_quants[id] = sample_quant['TPM']
            cts_quants[id] = sample_quant['counts']

    return tpm_quants,cts_quants

def read_qcs(pattern='salmon/*_salmon_out', tool='salmon', **kwargs):
    ''' Read technical quality control data results from every directory
    matching the glob in pattern.
    Parameters
    ----------
    tool, str, default 'salmon'
        The quantification tool used to generate the results. Currently
        supports 'salmon' and 'sailfish'.
    **kwargs,
        kwargs are passed on to the tool specific sample parser. See documentation
        for individual parsers for details.
    Returns
    -------
    A pandas.DataFrame where rows are samples, and columns are technical
    features extrated from the tool results.
    '''
    sample_readers = {
        'salmon': read_salmon_qc,
    }

    qc_reader = sample_readers[tool]

    QCs = pd.DataFrame()
    for sample_path in tqdm(iglob(pattern)):
        if os.path.isdir(sample_path):
            if not os.path.isfile('%s/quant.genes.sf'%sample_path):
                #os.system('rm %s -rf'%sample_path)
                print('Not found: %s'%sample_path)
                continue
            sample_qc = qc_reader(sample_path, **kwargs)
            id=sample_path.split('/')[-1]
            QCs[id] = sample_qc

    return QCs.T

class PosModel:
	def __init__(self):
		self.models = {}

	def from_file(self, fn):
		f = gzip.open(fn)
		b = f.read()

		offset = 0
		uint_size = struct.calcsize('I')
		double_size = struct.calcsize('d')

		num_models = struct.unpack_from('I', b[offset:])[0]
		offset += uint_size

		length_bins = []
		for i in range(0, num_models):
			length_bins.append(struct.unpack_from('I', b[offset:])[0])
			offset += uint_size

		models = []
		for i in range(0, num_models):
			model_bins = struct.unpack_from('I', b[offset:])[0]
			offset += uint_size
			model = list(struct.unpack_from('d'*model_bins, b[offset:]))
			offset += model_bins * double_size
			models.append(model)
		f.close()
		self.models = dict(zip(length_bins, models))
	
	
	
	
	
	
	
