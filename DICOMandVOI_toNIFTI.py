# DICOMandVOI_toNIFTI.py
#   Created 8/4/22 by Katie Merriman
#   Searches through designated folder for DICOM files of all patients on csv list
#   Converts DICOMs to Nifti, determines DICOM type, and saves to designated folder
#   Resamples ADC, B0, and DCE based on T2, if ADC and B0 exist
#   If ADC NIfTI with same resampling exists, high B can be calculated with calculate_highB.py after running this code

# Requirements
#   pydicom
#   SimpleITK
#   Pandas

# Inputs
#   1: csv file with column of MRNs to convert labeled as "MRN"
#           (MRNs can be JUST the MRN or in MRN_DateOfMRI format)
#   2: path to source folder with DICOMS
#   3: path to folder where converted files and csv reports should be saved
#   example usage:
#   python highB_resample.py "/path/to/cvsFileName.csv" "/path/to/source/folder" "/path/to/save/folder"

# Outputs
#   Nifti files saved as T2, ADC, DCE, or highB in folders created for each patient within designated save folder


import os
#import sys
import pydicom
from pydicom import dcmread
from pydicom.errors import InvalidDicomError
import SimpleITK as sitk
import csv
#import glob
import pandas as pd
import numpy as np
import glob
from skimage import draw


class DICOMtoNIFTI():
    def __init__(self):
        pathway = 0

        if pathway:
            self.csv_file = r'T:\MIP\Katie_Merriman\ZoeKatieCollab\Dicoms2Convert.csv'
            self.sourceFolder = r'T:\MRIClinical\surgery_cases'
            self.saveFolder = r'T:\MIP\Katie_Merriman\ZoeKatieCollab\convertedDICOM2'
        else:
            self.csv_file = 'Mdrive_mount/MIP/Katie_Merriman/ZoeKatieCollab/Dicoms2Convert.csv'
            self.sourceFolder = 'Mdrive_mount/MRIClinical/surgery_cases'
            self.saveFolder = 'Mdrive_mount/MIP/Katie_Merriman/ZoeKatieCollab/convertedDICOM2'

        self.DICOMconversions = []
        self.VOIconversions = []
        self.error_data = []

    def startConversion(self):
        patient = []

        df_csv = pd.read_csv(self.csv_file, sep=',', header=0)
        for rows, file_i in df_csv.iterrows():
            p = (str(file_i['MRN']))
            if len(p)<7:
                p ='0'+ p
            p = p + '*'
            patientPath = glob.glob(os.path.join(self.sourceFolder, p))
            if patientPath:
                patientPath = patientPath[0]
                patientID = os.path.basename(patientPath)
                patient.append([patientID, patientPath, 0])

        for i in range(0, len(patient)):
            patient[i][2] = self.sortDICOMS(patient[i])
            #patient[i][2] = convertedFiles[0]

    def sortDICOMS(self, p):
        success = 0
        DICOMsuccess = 0
        VOIsuccess = 0
        dicomFoldersList = []
        dicomSeriesList = []
        dicomProtocolList = []
        voi_list = []
        suffix_list = ['.stl', '.xml', '.csv', '.xlsx', '.doc', '.txt', '.jpg', '.png']
        ignore_list = ['DICOMDIR', 'LOCKFILE', 'VERSION']

        print('searching ',p[0],' for DICOMs and VOIs')
        for root, dirs, files in os.walk(p[1]):
            for name in files:

                ## check for and separate VOIs
                if name.endswith('.voi'):
                    if name.endswith('bt.voi'):
                        if 'gg' in name:
                            continue
                        if 'urethra' in name:
                            continue
                        else:
                            filePath = os.path.join(root, name)
                            print('voi found: ', filePath)
                            voi_list.append(filePath)

                ## ignore other common non-DICOM formats
                elif name.endswith(tuple(suffix_list)):
                    continue

                else:
                    if name not in ignore_list:
                        filePath = os.path.join(root, name)
                        try:
                            ds = dcmread(filePath)
                        except IOError:
                            # print(f'No such file')
                            continue
                        except InvalidDicomError:
                            # print(f'Invalid Dicom file')
                            continue
                        dicomString = filePath[:-(len(name) + 1)]
                        if dicomString not in dicomFoldersList:
                            if 'delete' in dicomString:
                                 break
                            dicomFoldersList.append(dicomString)
                            dicomSeries = ds.ProtocolName.replace('/', '-')
                            dicomSeries = dicomSeries.replace(" ", "_")
                            # print(f'DICOM found...')
                            ADClist = ['Apparent Diffusion Coefficient', 'adc', 'ADC', 'dWIP', 'dSSh', 'dReg']
                            if ('T2' in dicomSeries or 't2' in dicomSeries):
                                dicomSeriesType = 'T2'
                                print('T2 found')
                            elif (any([substring in dicomString for substring in ADClist])) or (
                            any([substring in dicomSeries for substring in ADClist])):
                                dicomSeriesType = 'ADC'
                                print('ADC found')
                            elif ('extracted' in filePath):
                                dicomSeriesType = 'B0_extracted'
                                print('B0 found')
                            else:
                                series_descript = ds.SeriesDescription
                                if any([substring in series_descript for substring in ADClist]) or (
                                        dicomString.endswith('ADC') or dicomString.endswith('adc')):
                                    dicomSeriesType = 'ADC'
                                    print('ADC found')
                                elif ('T2' in series_descript or 't2' in series_descript) or (
                                        dicomString.endswith('T2') or dicomString.endswith('t2')):
                                    dicomSeriesType = 'T2'
                                    print('T2 found')
                                elif (dicomString.endswith('DCE') or dicomString.endswith('dce')):
                                    dicomSeriesType = 'DCE'
                                    print('DCE found')
                                else:
                                    dicomSeriesType = 'highB'
                                    print('highB found')
                            dicomSeriesList.append(dicomSeries)
                            dicomProtocolList.append(dicomSeriesType)



        if dicomFoldersList:
            print('converting files')
            success = self.DICOMconvert(dicomFoldersList, dicomSeriesList, dicomProtocolList, voi_list, p)

        return success

    def DICOMconvert(self, dicomFoldersList, dicomSeriesList, dicomProtocolList, voi_list, p):
        success = 0
        MRN = p[0]
        adc_img = []
        highB_img = []
        B0_img = []
        t2_img = []
        niftiInfo = []
        convertedDICOMs = []
        convertedVOIfiles = []
        convertedVOInames = []

        # create path and file name for NIfTI
        patientPath = os.path.join(self.saveFolder, MRN)
        try:
            os.mkdir(patientPath)
        except FileNotFoundError:
            print(f'Error! Invalid path!')
        except FileExistsError:
            print(f'Patient folder already created')
        suffix = '.nii.gz'

        print(f'Converting DICOM to NIfTI...')  # let user know code is running
        # read, convert, and save NIfTIs
        reader = sitk.ImageSeriesReader()
        for refPath, DCMseries, dicomProtocol in zip(dicomFoldersList, dicomSeriesList, dicomProtocolList):

            # read original dicom
            #print(f'Reading images...')
            dicom_names = reader.GetGDCMSeriesFileNames(refPath)
            niftiInfo.append([MRN, DCMseries, dicomProtocol, patientPath])

            # convert to NIfTI
            reader.SetFileNames(dicom_names)
            #print(f'Converting files...')
            imageFile = os.path.join(patientPath, dicomProtocol + suffix)
            if dicomProtocol == 'ADC':
                adc_img = reader.Execute()
                print('     ADC ready to resample')
                niiADC = imageFile
            elif dicomProtocol == 'highB':
                highB_img = reader.Execute()
                print('     highB ready to resample')
                niihighB = imageFile
            elif dicomProtocol == 'B0_extracted':
                B0_img = reader.Execute()
                print('      B0 ready to resample')
                niiB0 = imageFile
            elif dicomProtocol == 'T2':
                t2_img = reader.Execute()
                print(f'     Writing T2 images...')
                sitk.WriteImage(t2_img, imageFile)
                convertedDICOMs.append('T2')
                success = 1
            else:
                image = reader.Execute()
                print(f'Saving ',dicomProtocol,' image')
                sitk.WriteImage(image, imageFile)
                convertedDICOMs.append('DCE')
                success = 1

        # Resample
        if t2_img:
            # set the resample filter based on T2 header info
            Filter = sitk.ResampleImageFilter()
            Filter.SetReferenceImage(t2_img)
            Filter.SetOutputDirection(t2_img.GetDirection())
            Filter.SetOutputOrigin(t2_img.GetOrigin())
            Filter.SetOutputSpacing(t2_img.GetSpacing())
            if adc_img:
                ## Resample ADC
                adc_resamp = Filter.Execute(adc_img)
                adc_resamp.CopyInformation(t2_img)

                # make sure all header info matches
                for meta_elem in t2_img.GetMetaDataKeys():
                    adc_resamp.SetMetaData(meta_elem, t2_img.GetMetaData(meta_elem))

                sitk.WriteImage(adc_resamp, niiADC)
                convertedDICOMs.append('ADC')
                print('saving resampled ADC')
            else:
                print('No ADC')
            if highB_img:
                ## Resample high B
                highB_resamp = Filter.Execute(highB_img)
                highB_resamp.CopyInformation(t2_img)

                # make sure all header info matches
                for meta_elem in t2_img.GetMetaDataKeys():
                    highB_resamp.SetMetaData(meta_elem, t2_img.GetMetaData(meta_elem))

                sitk.WriteImage(highB_resamp, niihighB)
                convertedDICOMs.append('highB')
                print('saving resampled highB')
            else:
                print('No highB')
            if B0_img:
                ## Resample high B
                B0_resamp = Filter.Execute(B0_img)
                B0_resamp.CopyInformation(t2_img)

                # make sure all header info matches
                for meta_elem in t2_img.GetMetaDataKeys():
                    B0_resamp.SetMetaData(meta_elem, t2_img.GetMetaData(meta_elem))

                sitk.WriteImage(highB_resamp, niiB0)
                convertedDICOMs.append('B0')
                print('saving resampled highB')
            else:
                print('No B0')
            if voi_list:
                for file in voi_list:
                    img_array = sitk.GetArrayFromImage(t2_img)
                    img_array = np.swapaxes(img_array, 2, 0)

                    # iterate over mask and update empty array with mask
                    numpy_mask = np.empty(img_array.shape)
                    mask_dict = self.mask_coord_dict(patient_id=MRN, file=file,
                                                     img_shape=(img_array.shape[0], img_array.shape[1]))
                    if mask_dict:
                        for key in mask_dict.keys():
                            numpy_mask[:, :, int(key)] = mask_dict[key]

                        numpy_mask = np.swapaxes(numpy_mask, 2, 0)
                        img_out = sitk.GetImageFromArray(numpy_mask)

                        # need to save as nifti
                        img_out.CopyInformation(t2_img)
                        for meta_elem in t2_img.GetMetaDataKeys():
                            img_out.SetMetaData(meta_elem, t2_img.GetMetaData(meta_elem))
                        niftiPath = os.path.basename(file)
                        niftiPath = str.replace(niftiPath, '.voi', '.nii.gz')
                        niftiPath = os.path.join(self.saveFolder, MRN, niftiPath)
                        sitk.WriteImage(img_out, niftiPath)
                        convertedVOIfiles.append(file)
                        VOIname = os.path.basename(file)
                        convertedVOInames.append(VOIname)

                        success = 2
                    else:
                        self.error_data.append([MRN, 'No keys in mask_dict', file])
        else:
            print('No T2')
        self.DICOMconversions.append([MRN,convertedDICOMs, convertedVOIfiles, convertedVOInames])



        return success

    def mask_coord_dict(self,patient_id='',file='',img_shape=()):
        '''
        creates a dictionary where keys are slice number and values are a mask (value 1) for area
        contained within .voi polygon segmentation
        :param patient_dir: root for directory to each patient
        :param type: types of file (wp,tz,urethra,PIRADS)
        :return: dictionary where keys are slice number, values are mask
        '''

        # define path to voi file
        #voi_path=os.path.join(self.dicom_folder, patient_id, file)
        voi_path = file

        #read in .voi file as pandas df
        pd_df = pd.read_fwf(voi_path)

        # use get_ROI_slice_loc to find location of each segment
        dict=self.get_ROI_slice_loc(voi_path)

        output_dict={}
        if dict:
            for slice in dict.keys():
                values=dict[slice]
                select_val=list(range(values[1],values[2]))
                specific_part=pd_df.iloc[select_val,:]
                split_df = specific_part.join(specific_part['MIPAV VOI FILE'].str.split(' ', 1, expand=True).rename(columns={0: "X", 1: "Y"})).drop(['MIPAV VOI FILE'], axis=1)
                X_coord=np.array(split_df['X'].tolist(),dtype=float).astype(int)
                Y_coord=np.array(split_df['Y'].tolist(),dtype=float).astype(int)
                mask=self.poly2mask(vertex_row_coords=X_coord, vertex_col_coords=Y_coord, shape=img_shape)
                output_dict[slice]=mask

        return(output_dict)

    def get_ROI_slice_loc(self,path):
        '''
        selects each slice number and the location of starting coord and end coord
        :return: dict of {slice number:(tuple of start location, end location)}

        '''

        pd_df=pd.read_fwf(path)

        #get the name of the file
        filename=path.split(os.sep)[-1].split('.')[0]

        #initialize empty list and empty dictionary
        slice_num_list=[]
        last_line=[]
        loc_dict={}

        #find the location of the last line -->
        for line in range(len(pd_df)):
            line_specific=pd_df.iloc[line,:]
            as_list=line_specific.str.split(r"\t")[0]
            if "# slice number" in as_list: #find location of all #slice numbers
                slice_num_list.append(line)
            if '# unique ID of the VOI' in as_list:
                last_line.append(line)

        if len(slice_num_list) < 1:
            return None
        else:
            for i in range(len(slice_num_list)):
                # for all values except the last value
                if i<(len(slice_num_list)-1):
                    loc=slice_num_list[i]
                    line_specific=pd_df.iloc[loc,:]
                    slice_num=line_specific.str.split(r"\t")[0][0]
                    start=slice_num_list[i]+3
                    end=slice_num_list[i+1]-1
                    loc_dict.update({slice_num:(filename,start,end)})

                #for the last value
                if i == (len(slice_num_list) - 1):
                    loc = slice_num_list[i]
                    line_specific=pd_df.iloc[loc,:]
                    slice_num=line_specific.str.split(r"\t")[0][0]
                    start=slice_num_list[i]+3
                    end=(last_line[0]-1)
                loc_dict.update({slice_num: (filename, start, end)})

        return (loc_dict)

    def poly2mask(self,vertex_row_coords, vertex_col_coords, shape):
        ''''''
        fill_row_coords, fill_col_coords = draw.polygon(vertex_row_coords, vertex_col_coords, shape)
        mask = np.zeros(shape, dtype=int)
        mask[fill_row_coords, fill_col_coords] = 1
        return mask

    def create_csv_files(self):
        nifti_cvsFileName = os.path.join(self.saveFolder, 'ConvertedDICOMs.csv')
        niftiHeader = ['MRN', 'Converted DICOMs', 'Converted VOIs']
        with open(nifti_cvsFileName, 'w', newline="") as file2write:
            csvwriter = csv.writer(file2write)
            csvwriter.writerow(niftiHeader)
            csvwriter.writerows(self.DICOMconversions)


if __name__ == '__main__':
    c = DICOMtoNIFTI()
    c.startConversion()
    c.create_csv_files()


