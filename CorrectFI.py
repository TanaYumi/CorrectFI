#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 20:26:22 2018

@author: tanakayumiko


This script corrects the Fluorescent Intensity (FI) of sequencial time-lapse images

imput  : FI_Data.csv, M_Roi_of_analyzewell.csv, M_Roi_of_emptywell.csv
output : corrected_FI_Data.csv

methods
----------
* read_file
* calc_cr_coef
* calc_cr_FI

Please excecute these methods sequentialy .

Example
----------
import CorrectFI as cf

wd = '~/data/'
FI_Data, empwell_Data = cf.read_file(wd, 'ROI_Data.csv', 'PositiveWellIndexes_empty.csv')
cr_coef_Data = cf.calc_cr_coef(FI_Data,empwell_Data)
cr_FI_Data = cf.calc_cr_FI(cr_coef_Data,FI_Data,medfilt=3)
cf.plot(cr_FI_Data)

"""    

import pandas as pd
import scipy.signal as ss
import matplotlib.pyplot as plt


def read_file(wd, FI_Data_file, empwell_file):
    
    """
    Read the FI_file and empty_well_file
    
    parameters
    ----------
    wd           : string
                   Path of working_directory
    
    FI_Data_file : string
                   filename of FI_file
    
    empwell_file : string
                   Filename of empty_well_file
              
    Returns
    ----------
    pd.DataFrame(FI_Data), pd.DataFrame(empwell_Data)
    
    """
        
    FI_Data = pd.read_csv(wd + FI_Data_file, index_col = 0)
    empwell_Data = pd.read_csv(wd + empwell_file, index_col=[0,1])
    
    return FI_Data, empwell_Data


    
def calc_cr_coef(FI_Data, empwell_Data, chm_end=False, roi_list=[1,2,3,4], channel='MeanIL13'):
    
    """
    Calculate the correction coefficient using the FI of empty wells
    
    
    parameters
    ----------
    FI_Data : pd.DataFrame
    
    empwell_Data : pd.DataFrame
    
    chm_end : False or nested list, default False
              If you calculate the back noise of each chamber separately,
              imput first and last Multipoint in each chamber, like [[s1,l1],[s2,l2]]
    
    roi_list: list
              Roi numbers
              
    channel : string, default 'MeanIL13'
              Channel name
              
    Returns
    ----------
    pd.DataFrame (correction_coefficients of each chamber and each roi)
    
    """
    
    #chm_end = [[1,80],[81,138]]
    
    #pivot table of FI_Data
    piv_FI_Data = pd.pivot_table(FI_Data,index=['ND.M','RoiID'],columns='ND.T',values=channel)
    
    #crop FI_Data of empty wells
    empty_FI_Data = piv_FI_Data.loc[empwell_Data.index]
    
    cr_coef_Data = pd.DataFrame()
    
    if chm_end ==False:
        chm_end = [[piv_FI_Data.index.min()[0], piv_FI_Data.index.max()[0]]]
    else:
        pass
    
    #chamber forloop
    d = {}
    for chm_n in range(len(chm_end)):
        
        #crop FI_Data of selected chamber
        crop_empty_FI_Data = empty_FI_Data.loc[chm_end[chm_n][0]:chm_end[chm_n][1],:]
        
        #calculate trim_mean of each time and each roi(ave_r_Data)
        ave_r_Data = pd.DataFrame(columns=crop_empty_FI_Data.columns)
        for roi in roi_list:
            crop_Data = crop_empty_FI_Data.xs(roi,level='RoiID')
            q_Data = crop_Data.quantile(q=[0.1, 0.9])
            ave_r_Data.loc[roi] = crop_Data[crop_Data>q_Data.loc[0.1]][crop_Data<q_Data.loc[0.9]].mean()
            
        #calculate trim_mean of each roi(ave_r_t_Data)
        t_q_Data = ave_r_Data.quantile(q=[0.1,0.9], axis=1)
        ave_r_t_Se = pd.Series()
        for roi in roi_list:
            crop_ave_Data = ave_r_Data.loc[roi]
            ave_r_t_Se.loc[roi] = \
            crop_ave_Data[crop_ave_Data>t_q_Data.loc[0.1][roi]][crop_ave_Data<t_q_Data.loc[0.9][roi]].mean()
        
        #calculate correction coefficient
        temp_cr_coef_Data = ave_r_Data.apply(lambda x: ave_r_t_Se/x)
        temp_cr_coef_Data.index.name = 'RoiID'
        
        #to concat dataframes of different chamber
        d[chm_n] = temp_cr_coef_Data
        
    cr_coef_Data = pd.concat(d.values(), keys=d.keys())
    
    return cr_coef_Data



def calc_cr_FI(cr_coef_Data, FI_Data, medfilt=False, subtmean=True, chm_end=False, roi_list=[1,2,3,4], channel='MeanIL13'):
    
    """
    calculate the corrected FI using the correction_coefficient
    
    Parameter
    ----------
    
    cr_coef_Data: pd.DataFrame
                   correction coefficient
                   
    FI_Data     : pd.DataFrame
                  FI_Data to correct(mandatory: ND.M, RoiID, ND.T, FIvalue)
               
    medfilt     : default is False, 
                  If you apply, imput kernel_size of median filter
                 
    subtmean    : True or False, default True
                  If True, corrected FI was subtracted by trimmean of each multipoint
                     
    chm_end     : False or nested list, default False
                  If you calculate the back noise of each chamber separately,
                  imput first and last Multipoint in each chamber, like [[s1,l1],[s2,l2]]
    
    roi_list    : list, default [1,2,3,4]
                  Roi numbers
    
    channel     : string, default 'MeanIL13'
                  Channel name
              
    Returns
    ----------
    pd.DataFrame (corrected FI_Data)
    
    """
    
    #make pivot table of FI_Data
    if len(FI_Data.index.names) == 1:
        crop_FI_Data = FI_Data[['ND.M','RoiID','MeanIL13','ND.T']]
        piv_FI_Data = pd.pivot_table(crop_FI_Data,index=['ND.M','RoiID'],columns='ND.T', values=channel)
        
    elif len(FI_Data.index.names) == 2:
        piv_FI_Data = FI_Data
    
    else:
        print('FI_Data must be pd.DataFrame containing [ND.M, RoiID, ND.T, FIvalue]')
        raise TypeError
    
    #chamber
    if chm_end ==False:
        chm_end = [[piv_FI_Data.index.min()[0], piv_FI_Data.index.max()[0]]]
    else:
        pass
    
    #chamber forloop
    cr_coef_FI_Data = pd.DataFrame()
    for chm_n in range(len(chm_end)):
        crop_cr_coef_Data = cr_coef_Data.loc[chm_n,:]
        
        for roi in roi_list:
            crop_piv_FI_Data = piv_FI_Data.loc[chm_end[chm_n][0]:chm_end[chm_n][1],:]
            crop_crFI_Data = crop_piv_FI_Data.xs(roi,level=("RoiID")) * crop_cr_coef_Data.loc[roi]
            crop_crFI_Data['RoiID'] = [roi]*len(crop_crFI_Data)
            cr_coef_FI_Data = pd.concat([cr_coef_FI_Data, crop_crFI_Data])
    
    reset_crFI_Data = cr_coef_FI_Data.reset_index()
    piv_crFI_Data = pd.pivot_table(reset_crFI_Data,index=['ND.M','RoiID'])
    
    #subtract FI of first time point
    crFI_Data = piv_crFI_Data.apply(lambda x: x-piv_crFI_Data[1])
    
    #apply subtmean
    if subtmean == False:
        pass
    else:
        sum_crFI_Data = crFI_Data.groupby('ND.M').sum()-crFI_Data.groupby('ND.M').max()
        crFI_Data = crFI_Data.groupby('ND.M').apply(lambda x:x-sum_crFI_Data.div(len(roi_list)-1))

    #apply median filter
    if medfilt == False:
        pass
    else:
        crFI_Data = crFI_Data.apply(lambda x: ss.medfilt(x,kernel_size=medfilt), axis=1)
        
    return crFI_Data


    
def plot(FI_Data, select=False, channel='MeanIL13'):
    
    """
    plot FI_Data
    
    parameter
    ----------
    FI_Data : pd.DataFrame
              FI_Data to plot(mandatory: ND.M, RoiID, ND.T, FIvalue)
    
    select  : list or pd.DataFrame.index, default False
              If you want to plot FI of selected wells, imput mupliindex of the wells
              otherwise, imput nested list, like [[m1,r1],[m2,r2],...]
              
    channel : string, default 'MeanIL13'
              Channel name
    
    Returns
    ----------
    """
    
    if len(FI_Data.index.names) == 1:
        crop_FI_Data = FI_Data[['ND.M','RoiID','MeanIL13','ND.T']]
        piv_FI_Data = pd.pivot_table(crop_FI_Data,index=['ND.M','RoiID'],columns='ND.T', values=channel)
        
    elif len(FI_Data.index.names) == 2:
        piv_FI_Data = FI_Data
    
    else:
        print('FI_Data must be pd.DataFrame containing [ND.M, RoiID, ND.T, FIvalue]')
        raise TypeError
    
    plt.figure(figsize=(6,6))
    
    if select == False:
        for i in piv_FI_Data.index:
            plt.plot(piv_FI_Data.loc[i])
            
    elif type(select) == list:
        for m,r in select:
            plt.plot(piv_FI_Data.loc[m,r])
    
    else:
        for i in select:
            plt.plot(piv_FI_Data.loc[i])
            
    plt.show()
    
    
 
