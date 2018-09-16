import os, sys, gzip, glob, string
import numpy as np
import pandas as pd
from shutil import copyfile
from subprocess import call
from collections import Counter
from wordcloud import WordCloud
from geopy.geocoders import ArcGIS
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
import matplotlib
matplotlib.use('Agg') # Stop matplotlib from opening a display window when generating images
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.pyplot import cm
from matplotlib import rcParams
from adjustText import adjust_text
rcParams.update({'figure.autolayout': True}) # Auto-resize plot

''' Main function for importing the input files, calling the ancestry c++ code, preprocessing the data, and calling the functions that create our images '''
def prepare_data():
	files = list(['/23andme.txt','/refpanel_public_ref.gz','/gwas_ref.tsv','/23andme.txt.gz'])
	[copyfile(sys.argv[index+1], files[index]) for index in range(len(sys.argv)-1)]
	with open(files[0], 'rb') as f_in, gzip.GzipFile(files[3], 'wb') as gzip_file:
		gzip_file.write(f_in.read())
	call(["/home/ancestry/src/ancestry", "-i", files[1], "-b", "2", "-ranN", "100000", "-g23", files[3], "-o", "/home/ancestry/test"])
	ref = pd.read_csv(files[2], sep='\t', index_col = [0]).rename(index=str, columns={'MAPPED_TRAIT':'trait'})
	ref = ref[pd.to_numeric(ref['P-VALUE'], errors='coerce').notnull()]
	ref = ref[pd.to_numeric(ref['RISK ALLELE FREQUENCY'], errors='coerce').notnull()]
	ref = ref.loc[(ref['P-VALUE'].astype(float) < 0.1) & (ref['RISK ALLELE FREQUENCY'].astype(float) < 0.25)][list(ref)]
	data23 = pd.read_csv(files[0], sep='\t', names = ['SNPS','chrom','pos','geno'],  skiprows = 20)
	ref['trait'] = ref['trait'].map(lambda trait: ','.join(str(trait).split(', ')))
	ref = ref.assign(rallele = [variant.split('-')[1] if '-' in variant else '?' for variant in ref['STRONGEST SNP-RISK ALLELE'].values])
	ref = ref.assign(variant = [str(row['SNPS']) + ':' +str(row['rallele']) for index, row in ref.iterrows()])
	data23 = data23.assign(variant1 = [str(row['SNPS']) + ':' + str(row['geno'])[0] for _, row in data23.iterrows()])
	data23 = data23.assign(variant2 = [str(row['SNPS']) + ':' + (str(row['geno'])[1] if len(str(row['geno'])) > 1 else '') for _, row in data23.iterrows()])
	same1 = ref.loc[ref['variant'].isin(data23['variant1']) | ref['variant'].isin(data23['variant2'])][list(ref)]
	same2 = data23.loc[data23['variant1'].isin(ref['variant']) | data23['variant2'].isin(ref['variant'])][list(data23)]
	c_dict = {'X': 23, 'Y' : 24, 'MT':25}
	filter_chrom = lambda c: float(c) if (type(c) != str or (type(c) == str and unicode(c, 'utf-8').isnumeric())) else c_dict[c]
	same2['chrom'] = same2['chrom'].map(filter_chrom)
	chrom_and_pos = np.asarray([same2.loc[same2['SNPS'] == snp][['pos','chrom']].values.reshape((2,)) for snp in same1['SNPS']])
	clusters = KMeans(n_clusters=150).fit_predict(chrom_and_pos)
	same1 = same1.assign(clusters = clusters)
	traits_per_snp = same1.drop_duplicates(["clusters", "trait"]).groupby('clusters')['trait'].apply(list)
	create_ancestry_image() # create the ancestry map
	create_wordcloud_image(same1) # create the word-cloud
	create_cluster_image(dict(Counter(same1['trait'].values)),traits_per_snp) # create the t-sne trait-association picture

''' Function for creating your ancestry map.'''
def create_ancestry_image(): 
	anc_loc = {'ASHKENAZI': 'Poland', 'WEURASIA': 'Europe', 'BALOCHI-MAKRANI-BRAHUI': 'Pakistan', 'INDPAK': 'India', 'BANTUKENYA': 'Kenya', 'EAFRICA': 'Somalia', 'AFRICA': 'Africa', 'BANTUNIGERIA': 'Nigeria', 'WAFRICA': 'Ghana', 'BIAKA': 'Congo', 'CAFRICA': 'Chad', 'CAMBODIA-THAI': 'Cambodia', 'SEASIA': 'Phillipines', 'EASIA': 'Indonesia', 'CSAMERICA': 'Panama', 'AMERICAS': 'Nicaragua', 'CYPRUS-MALTA-SICILY': 'Sicily', 'EMED': 'Lebanon', 'SBALKANS': 'Serbia', 'ITALY': 'Itay', 'SWEUROPE': 'Sweden', 'EASTSIBERIA': 'Siberia', 'NEASIA': 'South Korea', 'FINNISH': 'Finland', 'NEEUROPE': 'England', 'GAMBIA': 'Gambia', 'GUJARAT': 'Gujarat', 'GUJARAT_PATEL': 'Gujarat', 'HADZA': 'Tanzania', 'HAZARA-UYGUR-UZBEK': 'Mongolia', 'CASIA': 'Kazakhstan', 'JAPAN-KOREA': 'Japan', 'KALASH': 'Pakistan', 'MENDE': 'Sierra Leone', 'NAFRICA': 'Tunisia', 'NCASIA': 'Kazakhstan', 'NEAREAST': 'Turkey', 'NEUROPE': 'Denmark', 'NGANASAN': 'Siberia', 'OCEANIA': 'Australia', 'PATHAN-SINDHI-BURUSHO': 'Pakistan', 'SAFRICA': 'South Africa', 'SAMERICA': 'Brazil', 'SARDINIA': 'Sardinia', 'SSASIA': 'Phillipines', 'BENGALI': 'India', 'TAIWAN': 'Taiwan', 'TUBALAR': 'Russia', 'TURK-IRAN-CAUCASUS': 'Iran'}
	anc_prop,anc_coord,anc_count = {}, {}, {}
	for anc_prop_fn in glob.glob('/home/ancestry/test*.Q'):
		with open(anc_prop_fn,'r') as f_in:
			for entry in f_in.readlines():
				anc_prop[entry.split(' ')[0]]  = anc_prop[entry.split(' ')[0]] + float(entry.split(' ')[1]) if entry.split(' ')[0] in anc_prop else float(entry.split(' ')[1])
				anc_count[entry.split(' ')[0]] = anc_count[entry.split(' ')[0]] + 1 if entry.split(' ')[0] in anc_count else 1
	arc_gis = ArcGIS(user_agent="ancestry_visualiser",timeout=3)
	for anc_group in anc_prop.keys():
		anc_prop[anc_group] /= anc_count[anc_group]
		anc_coord[anc_group] = arc_gis.geocode(anc_loc[anc_group]) if anc_prop[anc_group] > 0.01 else None
	map = Basemap()
	map.drawcoastlines()
	map.drawmapboundary(fill_color='aqua')
	map.fillcontinents(color='coral',lake_color='aqua')
	longi, lati = map([loc.longitude if loc != None else -200 for loc in anc_coord.values()], [loc.latitude if loc != None else -200 for loc in anc_coord.values()])
	s = np.asarray([500*anc_prop[key] if anc_prop[key] != None else 0 for key in anc_coord.keys()])
	map.scatter(longi, lati, s = s, marker='o',color='r')
	for index,anc_group in enumerate(anc_coord.keys()):
		rank = np.argsort(np.asarray(list(anc_prop.values())))[index]
		plt.text(longi[index],lati[index],str(rank+1),fontsize=12,fontweight='bold')
		plt.text(5, 200 + 15*rank, str(rank+1) + ') ' + anc_group + ' ' + str(anc_prop[anc_group]), fontsize=12,fontweight='bold')
	plt.savefig('./ancestry_map.png',bbox_inches="tight", dpi=400)

''' Function for creating a word cloud of your unique SNP traits '''
def create_wordcloud_image(same1): 
	plt.figure()
	frequencies = Counter([item for sublist in same1[['SNPS','trait']].groupby('SNPS')['trait'].apply(list).values for item in sublist])
	plt.imshow(WordCloud(max_font_size=30).generate_from_frequencies(frequencies), interpolation="bilinear")
	plt.axis("off")
	plt.savefig('./wordcloud.png', dpi=400)

''' Function for creating a cluster image of your traits '''
def create_cluster_image(trait_dict,snp_traits):
	top = [tup[0] for tup in Counter('#'.join(['#'.join(t) for t in snp_traits.values]).split('#')).most_common(20)]
	ohe_traits = np.asarray([[1 if t in ts else 0 for t in trait_dict.keys()] for ts in snp_traits])
	occurrences = [np.nonzero(ohe_traits[:,list(trait_dict.keys()).index(trait)]) for trait in top]
	tsne_values = TSNE(n_components=2).fit_transform(ohe_traits)*5
	plt.figure()
	txts = []
	c = [color for color in iter(cm.rainbow(np.linspace(0,1,20)))]
	for indx in range(20):
		xy = (tsne_values[occurrences[indx],0],tsne_values[occurrences[indx],1])
		txts.append(plt.text(np.mean(xy[0]),np.mean(xy[1]), '\n'.join(top[indx].split(' ')), fontsize=2, alpha=0.8, color=c[indx]))
		plt.scatter(xy[0], xy[1],marker=(2, 2, indx*360/(20*2)), s = 40, alpha=0.5, c=c[indx])
	adjust_text(txts)
	plt.axis("off")
	plt.savefig('./clusters.png',dpi=400)

if __name__ == "__main__":
	prepare_data()