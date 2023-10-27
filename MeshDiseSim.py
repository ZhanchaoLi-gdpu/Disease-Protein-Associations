# -*- coding: utf-8 -*-
"""
Created on Sun May 22 21:30:00 2022

@author: ly
"""

import ssmpy
import numpy as np
import scipy.io as scio

#ssmpy.create_semantic_base("mesh.nt", "mesh.db", "http://id.nlm.nih.gov/mesh/", "http://id.nlm.nih.gov/mesh/vocab#broaderDescriptor", '')
ssmpy.semantic_base("mesh.db")

MeSHName = []
for i in range (1,100000):
    temp = ssmpy.get_name(i)
    MeSHName.append(temp)
 
MeSHName = np.array(MeSHName) 
MeSHName = np.unique(MeSHName)
MeSHName = np.delete(MeSHName,0)
del i
del temp

DiseSimeResnik = np.zeros((len(MeSHName),len(MeSHName)))
DiseSimeLin = np.zeros((len(MeSHName),len(MeSHName)))
DiseSimeJiang = np.zeros((len(MeSHName),len(MeSHName)))

for i in range (0,30263):
    for j in range (i+1,30264):
        e1 = ssmpy.get_id(MeSHName[i,])
        e2 = ssmpy.get_id(MeSHName[j,])
        
        temp1 = ssmpy.ssm_resnik(e1,e2)
        DiseSimeResnik[i,j] = temp1
        temp2 = ssmpy.ssm_lin(e1,e2)
        DiseSimeLin[i,j] = temp2
        temp3 = ssmpy.ssm_jiang_conrath(e1,e2)
        DiseSimeJiang[i,j] = temp3
        
scio.savemat('DiseSimeResnik.mat', {'DiseSimeResnik': DiseSimeResnik}) 
scio.savemat('DiseSimeLin.mat', {'DiseSimeLin': DiseSimeLin}) 
scio.savemat('DiseSimeJiang.mat', {'DiseSimeJiang': DiseSimeJiang})
scio.savemat('MeSHName.mat', {'MeSHName': MeSHName})
  
np.save("DiseSimeResnik.npy",DiseSimeResnik)        
np.save("DiseSimeLin.npy",DiseSimeLin) 
np.save("DiseSimeJiang.npy",DiseSimeJiang) 
np.save("MeSHName.npy",MeSHName) 
















