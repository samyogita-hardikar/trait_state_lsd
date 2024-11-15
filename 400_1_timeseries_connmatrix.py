#python 

import numpy as np

from nilearn.datasets import fetch_atlas_schaefer_2018
dataset = fetch_atlas_schaefer_2018(n_rois=400)
atlas_filename = dataset.maps
labels = dataset.labels

from nilearn.input_data import NiftiLabelsMasker
masker = NiftiLabelsMasker(labels_img=atlas_filename, standardize=True)


sub=['sub-010005'	, 'sub-010009'	, 'sub-010011'	, 'sub-010014'	, 'sub-010015'	, 'sub-010016'	, 'sub-010017'	, 'sub-010018'	, 'sub-010020'	, 'sub-010021'	, 'sub-010022'	, 'sub-010023'	, 'sub-010024'	, 'sub-010027'	, 'sub-010029'	, 'sub-010030'	, 'sub-010031'	, 'sub-010033'	, 'sub-010034'	, 'sub-010035'	, 'sub-010038'	, 'sub-010040'	, 'sub-010042'	, 'sub-010052'	, 'sub-010056'	, 'sub-010058'	, 'sub-010060'	, 'sub-010061'	, 'sub-010062'	, 'sub-010064'	, 'sub-010065'	, 'sub-010067'	, 'sub-010068'	, 'sub-010069'	, 'sub-010071'	, 'sub-010072'	, 'sub-010073'	, 'sub-010074'	, 'sub-010075'	, 'sub-010076'	, 'sub-010077'	, 'sub-010078'	, 'sub-010079'	, 'sub-010080'	, 'sub-010082'	, 'sub-010084'	, 'sub-010094'	, 'sub-010096'	, 'sub-010097'	, 'sub-010098'	, 'sub-010099'	, 'sub-010101'	, 'sub-010102'	, 'sub-010103'	, 'sub-010105'	, 'sub-010106'	, 'sub-010107'	, 'sub-010108'	, 'sub-010111'	, 'sub-010112'	, 'sub-010113'	, 'sub-010114'	, 'sub-010115'	, 'sub-010116'	, 'sub-010117'	, 'sub-010118'	, 'sub-010120'	, 'sub-010121'	, 'sub-010122'	, 'sub-010123'	, 'sub-010124'	, 'sub-010125'	, 'sub-010126'	, 'sub-010127'	, 'sub-010128'	, 'sub-010129'	, 'sub-010130'	, 'sub-010131'	, 'sub-010132'	, 'sub-010135'	, 'sub-010136'	, 'sub-010137'	, 'sub-010138'	, 'sub-010141'	, 'sub-010142'	, 'sub-010143'	, 'sub-010144'	, 'sub-010145'	, 'sub-010147'	, 'sub-010151'	, 'sub-010154'	, 'sub-010155'	, 'sub-010157'	, 'sub-010158'	, 'sub-010159'	, 'sub-010161'	, 'sub-010162'	, 'sub-010163'	, 'sub-010164'	, 'sub-010165'	, 'sub-010167'	, 'sub-010168'	, 'sub-010171'	, 'sub-010173'	, 'sub-010174'	, 'sub-010176'	, 'sub-010185'	, 'sub-010186'	, 'sub-010187'	, 'sub-010188'	, 'sub-010189'	, 'sub-010190'	, 'sub-010191'	, 'sub-010192'	, 'sub-010194'	, 'sub-010195'	, 'sub-010196'	, 'sub-010200'	, 'sub-010201'	, 'sub-010203'	, 'sub-010204'	, 'sub-010205'	, 'sub-010206'	, 'sub-010208'	, 'sub-010209'	, 'sub-010210'	, 'sub-010211'	, 'sub-010212'	, 'sub-010213'	, 'sub-010216'	, 'sub-010217'	, 'sub-010218'	, 'sub-010220'	, 'sub-010221'	, 'sub-010224'	, 'sub-010225'	, 'sub-010226'	, 'sub-010227'	, 'sub-010228'	, 'sub-010229'	, 'sub-010230'	, 'sub-010231'	, 'sub-010232'	, 'sub-010233']

filename= 'all_lsd/derivatives/%s/func/%s_ses-02_task-rest_acq-%s_MNI2mm.nii.gz'
outfile = 'all_lsd/gradients_400_parcels/timeseries/timeseries_%s.npy'
group = []
time_series = []
for s in sub:
    try:
        ts_all  =  []
        fmri_filenames = [filename%(s, s, 'AP_run-01'), filename%(s, s, 'AP_run-02'), filename%(s, s, 'PA_run-01'), filename%(s, s, 'PA_run-02')]
        for f in fmri_filenames:
            try:
                ts = masker.fit_transform(f)
                ts_all.append(ts)
            except:
                pass
        time_series = np.concatenate(ts_all)
        np.save(outfile%(s), time_series)
        group.append(time_series)
        np.save('all_lsd/gradients_400_parcels/timeseries/timeseries_group.npy', group) 
    except:
        pass




from nilearn.connectome import ConnectivityMeasure
correlation_measure = ConnectivityMeasure(kind='correlation')
conn_matrix_ind = correlation_measure.fit_transform([group][0])
conn_matrix_ind_fisher = np.arctanh([conn_matrix_ind][0])
conn_matrix_ind_fisher[np.isposinf([conn_matrix_ind_fisher][0])]=1


np.save('all_lsd/gradients_400_parcels/conn_matrix_ind.npy', conn_matrix_ind)
np.save('all_lsd/gradients_400_parcels/conn_matrix_ind_fisher.npy', conn_matrix_ind_fisher)

#correlation matrix across all subjects

conn_matrix_group = correlation_measure.mean_
np.save('all_lsd/gradients_400_parcels/conn_matrix_group', conn_matrix_group)

conn_matrix_group_fisher = np.mean(conn_matrix_ind_fisher, axis=0)
conn_matrix_group_fisher[np.isposinf(conn_matrix_group_fisher)]=1
np.save('all_lsd/gradients_400_parcels/conn_matrix_group_fisher.npy', conn_matrix_group_fisher)





