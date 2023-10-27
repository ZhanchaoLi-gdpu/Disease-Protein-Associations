# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 11:31:30 2022

@author: Administrator
"""

import networkx as nx
from karateclub import DeepWalk    
import scipy.io as scio
import numpy as np

G=nx.read_weighted_edgelist('HumanPDNetwork.txt')
#G=nx.read_edgelist('HumanPPINetwork3.txt')
G = nx.convert_node_labels_to_integers(G, first_label=0, ordering='default')
#len(G.edges)
#len(G.nodes)

###############################################################################
#DeepWalk
modelDeepWalk = DeepWalk(walk_number=500, walk_length=80, dimensions=128, workers=20, window_size=5, epochs=1, learning_rate=0.05, min_count=1, seed=42)
#DeepWalk(walk_number: int = 10, walk_length: int = 80, dimensions: int = 128, workers: int = 4, window_size: int = 5, epochs: int = 1, learning_rate: float = 0.05, min_count: int = 1, seed: int = 42)
modelDeepWalk.fit(G)    
DeepWalk_embedding = modelDeepWalk.get_embedding()
scio.savemat('DeepWalk_embedding128.mat',{'DeepWalk_embedding': DeepWalk_embedding})
#DeepWalk_embedding = scio.loadmat('DeepWalk_embedding.mat')
np.save("DeepWalk_embedding128.npy",DeepWalk_embedding)
#DeepWalk_embedding = np.load("DeepWalk_embedding.npy")

