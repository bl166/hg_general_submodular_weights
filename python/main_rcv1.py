import numpy as np
from scipy.io import savemat
import os
import data

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

#%%
class_0, class_1 = ['C22'], ['C23']
categories = sorted(class_0 + class_1)

sparsity_b1, sparsity_b2 = 0.1, 0.002
sparsity_b1, sparsity_b2 = 0.01, 0.001
num_words = 30

save_to = f'../matlab/data/rcv1/data.mat'
print(save_to)

#%%
dataset = data.TextRCV1(data_dir='./data/RCV1', subset='all', categories=categories)
dataset.clean_text(num='substitute')
dataset.vectorize(stop_words='english')
dataset.remove_short_documents(nwords=50, vocab='full')

dataset.remove_encoded_images()
dataset.remove_frequent_words(sparsity_b1=sparsity_b1, sparsity_b2=sparsity_b2)
dataset.keep_top_words(num_words, 0)
dataset.remove_short_documents(nwords=5, vocab='selected')

dataset.compute_tfidf()
dataset.data_info(show_classes=True)

tfidf = dataset.tfidf.astype(float).T.toarray()  # size: (num of words) x (num of documents)
print('max/min tfidf', np.max(tfidf), np.min(tfidf[tfidf > 0]))

card = np.sum(tfidf > 0, 1)
print('max/min edge cardinality', max(card), min(card))

#%%
index2class = {i: dataset.class_names[i] for i in range(len(dataset.class_names))}
true_classes = np.array([int(index2class[i] in class_1) for i in dataset.labels])

#%%
data_mat = {'R': tfidf, 'y': true_classes}
savemat(save_to, data_mat)
