# HemorrhageDataPreprocessing.py
#   Created 8/8/22 by Katie Merriman
#   Searches through designated folder for DICOM files of all patients on csv list
#   Determines DICOM orientation, sorts DICOMS, converts to NIfTI, and saves to designated folder

# Requirements
#   pydicom
#   SimpleITK
#   csv
#   pandas
#   shutil

# Inputs (currently in class definition)
#   1: csv file with column of MRNs to convert labeled as "MRN"
#   2: path to source directory with all patients' DICOMS
#           - converted files will be saved in each patient folder and csv report will be saved in top level folder

# Outputs
#   1: Sorted DICOMs saved in folders within each DICOM series folder
#   2: Converted NIfTIs saved in folders in each patient folder
#   3: Csv report saved in top level source folder

import os
from pydicom import dcmread
from pydicom.errors import InvalidDicomError
import SimpleITK as sitk
import csv
import pandas as pd
import shutil



class DICOMconverter():
    def __init__(self):
        pathways = 1
            # I just do this so I can easily switch between the versions of the filepaths needed to debug on my computer vs run on Lambda.
            # if on lambda, I just change "pathways" to 0 instead of commenting/uncommenting multiple lines

        if pathways:
            self.csv_file = r'T:\MIP\Robert Huang\2022_06_01\other\hemorrhage\patients2convert.csv'
            self.patientFolder = r'T:\MIP\Robert Huang\2022_06_01\other\hemorrhage'

        else:
            self.csv_file = 'Mdrive_mount/MIP/Robert Huang/2022_06_01/other/hemorrhage/patients2convert.csv'
            self.patientFolder = 'Mdrive_mount/MIP/Robert Huang/2022_06_01/other/hemorrhage'


        self.DICOMconversions = []

    def convertPatients(self):
        patient = []

        patientFolders = os.listdir(self.patientFolder)
        df_csv = pd.read_csv(self.csv_file, sep=',', header=0)
        for rows, file_i in df_csv.iterrows():
            folderList = []
            p = (str(file_i['MRN']))  # get patient MRNs
            if len(p)<7:              # adds leading 0 if needed
                p ='0'+ p

            #find and store path to patient folder with MRN
            # p = p + '*'               # allows search regardless of what follows MRN in folder names
            folderList = [x for x in patientFolders if p in x]
            for patientID in folderList:
                patientPath = os.path.join(self.patientFolder, patientID)
                patient.append([patientID, patientPath, 0])  # 0 will be updated with number of successful conversions below


        for i in range(0, len(patient)):
            patient[i][2] = self.sortDICOMS(patient[i]) # stores number of successful conversions

        self.create_csv_files(patient)
        print('conversion complete')

    def sortDICOMS(self, p):
        success = 0
        dicomFoldersList = []
        suffix_list = ['.voi', '.stl', '.xml', '.csv', '.xlsx', '.doc', '.txt', '.jpg', '.png']
        ignore_list = ['DICOMDIR', 'LOCKFILE', 'VERSION']
        print('searching ', p[0], ' for DICOMs')

        for root, dirs, files in os.walk(p[1]):
            for name in files:

                ## ignore common non-DICOM files
                if name.endswith(tuple(suffix_list)):
                    continue

                ## look for DICOMS and sort by orientation
                else:
                    if name not in ignore_list:
                        filePath = os.path.join(root, name)

                        # check if file is DICOM, get orientation
                        try:
                            ds = dcmread(filePath)
                            orthog = ds.ImageOrientationPatient
                            if abs(orthog[1]) > .5:
                                orientation = 'sagittal'
                            elif abs(orthog[4]) > .5:
                                orientation = 'axial'
                            else:
                                orientation = 'coronal'
                        except IOError:
                            # print(f'No such file')
                            break
                        except InvalidDicomError:
                            # print(f'Invalid Dicom file')
                            break

                        # create and save a new orientation-based folder within DICOM folder
                        dicomFolder = os.path.dirname(filePath)
                        orthogPath = os.path.join(dicomFolder, orientation + '_DICOM')
                        if not os.path.exists(orthogPath):
                            try:
                                os.mkdir(orthogPath)
                                dicomFoldersList.append([orthogPath, orientation])  # add to list of DICOM folders for later convertion to NIfTI
                            except FileNotFoundError:
                                print(f'Error! Invalid path!')


                        # copy DICOM file to orientation-based folder
                        shutil.copy(filePath, orthogPath)

        if dicomFoldersList:
            print('converting files')
            success = self.DICOMconvert(dicomFoldersList, p)

        return success

    def DICOMconvert(self, dicomFoldersList, p):
        success = 0
        c_copy = 1
        s_copy = 1
        a_copy = 1

        reader = sitk.ImageSeriesReader()
        for folder in dicomFoldersList:
            # create path and file name for NIfTI
            savePath = os.path.join(p[1], folder[1] + '_NIfTI')
            try:
                os.mkdir(savePath)
                incr = 0  # file is first of that orientation being saved
            except FileNotFoundError:
                print(f'Error! Invalid path!')
            except FileExistsError:
                incr = 1  # another file of same orientation exists - need to increment filename
            if folder[1] == 'sagittal':
                s_copy = s_copy + incr
                copy = s_copy
            elif folder[1] == 'coronal':
                c_copy = c_copy + incr
                copy = c_copy
            else:
                a_copy = a_copy + incr
                copy = a_copy

            print(f'Converting', folder[1], copy, ' DICOM to NIfTI...')  # let user know code is running
            # read, convert, and save NIfTIs
            f = folder[0]
            dicom_names = reader.GetGDCMSeriesFileNames(f)
            reader.SetFileNames(dicom_names)
            imageFile = os.path.join(savePath, folder[1] + str(copy) + '.nii.gz')
            image = reader.Execute()
            sitk.WriteImage(image, imageFile)
            success = success + 1
            self.DICOMconversions.append(imageFile)

        return success


    def create_csv_files(self, list):
        nifti_cvsFileName = os.path.join(self.patientFolder, 'DICOMconversions.csv')
        niftiHeader = ['MRN', 'Filepath', '#DICOMS converted']
        with open(nifti_cvsFileName, 'w', newline="") as file2write:
            csvwriter = csv.writer(file2write)
            csvwriter.writerow(niftiHeader)
            csvwriter.writerows(list)


if __name__ == '__main__':
    c = DICOMconverter()
    c.convertPatients()


