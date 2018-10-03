# FI-corrector-for-timelapse-images
This script corrects the Fluorescent Intensity (FI) of sequential time-lapse images.

imput  : FI_Data.csv, M_Roi_of_analyzewell.csv, M_Roi_of_emptywell.csv
output : corrected_FI_Data.csv

Methods
----------
* read_file
* calc_cr_coef
* calc_cr_FI

Please excecute these methods sequentialy .

*read_file | Read the FI_file and empty_well_file
    
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
    
    
    
*calc_cr_coef | Calculate the correction coefficient using the FI of empty wells
    
    
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
    
    

*calc_cr_FI | calculate the corrected FI using the correction_coefficient
    
    parameter
    ----------
    
    cr_coef_Data: pd.DataFrame
                   correction coefficient
                   
    FI_Data     : pd.DataFrame
                  FI_Data to correct(mandatory: ND.M, RoiID, ND.T, FIvalue)
               
    medfilt     : default is False, 
                  If you apply, imput kernel_size of median filter
                 
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
    

    
*plot | plot FI_Data
    
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
    ----


Example
----------
import CorrectFI as cf

wd = '~/data/'
FI_Data, empwell_Data = cf.read_file(wd, 'ROI_Data.csv', 'PositiveWellIndexes_empty.csv')
cr_coef_Data = cf.calc_cr_coef(FI_Data,empwell_Data)
cr_FI_Data = cf.calc_cr_FI(cr_coef_Data,FI_Data,medfilt=3)
cf.plot(cr_FI_Data)
