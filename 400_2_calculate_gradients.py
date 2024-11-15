import warnings
warnings.simplefilter('ignore')
import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
import os 
import subprocess 



sub=['sub-010005'	, 'sub-010009'	, 'sub-010011'	, 'sub-010014'	, 'sub-010015'	, 'sub-010016'	, 'sub-010017'	, 'sub-010018'	, 'sub-010020'	, 'sub-010021'	, 'sub-010022'	, 'sub-010023'	, 'sub-010024'	, 'sub-010027'	, 'sub-010029'	, 'sub-010030'	, 'sub-010031'	, 'sub-010033'	, 'sub-010034'	, 'sub-010035'	, 'sub-010038'	, 'sub-010040'	, 'sub-010042'	, 'sub-010052'	, 'sub-010056'	, 'sub-010058'	, 'sub-010060'	, 'sub-010061'	, 'sub-010062'	, 'sub-010064'	, 'sub-010065'	, 'sub-010067'	, 'sub-010068'	, 'sub-010069'	, 'sub-010071'	, 'sub-010072'	, 'sub-010073'	, 'sub-010074'	, 'sub-010075'	, 'sub-010076'	, 'sub-010077'	, 'sub-010078'	, 'sub-010079'	, 'sub-010080'	, 'sub-010082'	, 'sub-010084'	, 'sub-010094'	, 'sub-010096'	, 'sub-010097'	, 'sub-010098'	, 'sub-010099'	, 'sub-010101'	, 'sub-010102'	, 'sub-010103'	, 'sub-010105'	, 'sub-010106'	, 'sub-010107'	, 'sub-010108'	, 'sub-010111'	, 'sub-010112'	, 'sub-010113'	, 'sub-010114'	, 'sub-010115'	, 'sub-010116'	, 'sub-010117'	, 'sub-010118'	, 'sub-010120'	, 'sub-010121'	, 'sub-010122'	, 'sub-010123'	, 'sub-010124'	, 'sub-010125'	, 'sub-010126'	, 'sub-010127'	, 'sub-010128'	, 'sub-010129'	, 'sub-010130'	, 'sub-010131'	, 'sub-010132'	, 'sub-010135'	, 'sub-010136'	, 'sub-010137'	, 'sub-010138'	, 'sub-010141'	, 'sub-010142'	, 'sub-010143'	, 'sub-010144'	, 'sub-010145'	, 'sub-010147'	, 'sub-010151'	, 'sub-010154'	, 'sub-010155'	, 'sub-010157'	, 'sub-010158'	, 'sub-010159'	, 'sub-010161'	, 'sub-010162'	, 'sub-010163'	, 'sub-010164'	, 'sub-010165'	, 'sub-010167'	, 'sub-010168'	, 'sub-010171'	, 'sub-010173'	, 'sub-010174'	, 'sub-010176'	, 'sub-010185'	, 'sub-010186'	, 'sub-010187'	, 'sub-010188'	, 'sub-010189'	, 'sub-010190'	, 'sub-010191'	, 'sub-010192'	, 'sub-010194'	, 'sub-010195'	, 'sub-010196'	, 'sub-010200'	, 'sub-010201'	, 'sub-010203'	, 'sub-010204'	, 'sub-010205'	, 'sub-010206'	, 'sub-010208'	, 'sub-010209'	, 'sub-010210'	, 'sub-010211'	, 'sub-010212'	, 'sub-010213'	, 'sub-010216'	, 'sub-010217'	, 'sub-010218'	, 'sub-010220'	, 'sub-010221'	, 'sub-010224'	, 'sub-010225'	, 'sub-010226'	, 'sub-010227'	, 'sub-010228'	, 'sub-010229'	, 'sub-010230'	, 'sub-010231'	, 'sub-010232'	, 'sub-010233']

#load individual connectivity matrices
conn_matrix_ind_fisher= np.load('all_lsd/gradients_400_parcels/conn_matrix_ind_fisher.npy')
#load group connectivity matrix
conn_matrix_group_fisher= np.load('all_lsd/gradients_400_parcels/conn_matrix_group_fisher.npy')


from brainspace.gradient import GradientMaps
from brainspace.datasets import load_group_fc, load_parcellation, load_conte69
from brainspace.plotting import plot_hemispheres
from brainspace.utils.parcellation import map_to_labels


#set number of components
NOcomp = 10

#set kernel
kernel = 'normalized_angle'

#set approach
approach = 'pca'

#set interactive plotting
plt.ion()

#Load 400 Schaefer parcellation
labeling = load_parcellation('schaefer', scale=400, join=True)
#and load the conte69 hemisphere surfaces
surf_lh, surf_rh = load_conte69()
mask = labeling != 0

#load HCP data from BrainSpace and calculate group level gradients
conn_matrix = load_group_fc('schaefer', scale=400)
ref_gradient = GradientMaps(n_components=NOcomp, kernel= kernel, approach = approach)
ref_gradient.fit(conn_matrix)

#np.savetxt('all_lsd/gradients_400_parcels/HCP_pca.csv', ref_gradient.gradients_[0], delimiter =",")

#calculate my group level gradients
group_matrix = conn_matrix_group_fisher
group_gradient = GradientMaps(n_components=NOcomp, kernel= kernel, approach = approach, alignment='procrustes')
#align my group level gradients to the HCP gradients
group_gradient.fit([group_matrix], reference = ref_gradient.gradients_)
np.savetxt('all_lsd/gradients_400_parcels/group_grad_pca_HCP_aligned.csv', group_gradient.aligned_[0], delimiter =",")

#calculate individual gradients aligned to the group level gradient gradients
  
#For aligning all subjects together
newarr = []
for i in range(len(conn_matrix_ind_fisher)):
    newarr.append(conn_matrix_ind_fisher[i])

ind_gradient = GradientMaps(n_components=NOcomp, kernel= kernel, approach = approach, alignment='procrustes')
ind_gradient.fit(newarr[0:144],reference = group_gradient.aligned_[0])

#save newly aligned ind gradients as csv
ind_grad=[]
for s in range(len(sub)) :
    ind_grad = ind_gradient.aligned_[s]
    np.savetxt('all_lsd/gradients_400_parcels/grad_pca_group_aligned_' +str(sub[s])+'.csv', ind_grad, delimiter =",")
   




