import pandas as pd
import numpy as np
import h5py
import datetime
import math
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from scipy import stats


def spec_fit(x, c0, c1, c2, c3, w1, w2, f1):
    '''
    A function for fitting peaks to spectra.
    This fit function was copy-pasted in as text from David's fit function on 5/11/21.
    There are two water peaks (free water and bound water) modeled as gaussian peaks w1 and w2.
    There is the f1 fatty acid peak.
    There is a three-term polynomial function.

    :param x: input data
    :param c0: polynomial term 1
    :param c1: polynomial term 2
    :param c2: polynomial term 3
    :param c3: polynomial term 4
    :param w1: water peak term 1
    :param w2: water peak term 2
    :param f1: fatty acid peak
    :return:
    '''

    f_of_x = c0+c1*(x-1700)+c2*(x-1700)**2+c3*(x-1700)**3 + \
             w1*(0.0751747*np.exp(-(x-1903.82)**2/26.4725**2)+0.225213*np.exp(-(x-1942.21)**2/48.8781**2) +
                 0.005*np.exp(-(x-1779.71)**2/32.1869**2))/7.715 + \
             w2*(0.0280945*np.exp(-(x-1913.6)**2/25.0449**2)+0.103527*np.exp(-(x-1949.5)**2/52.2024**2))/3.07 + \
             f1*(0.31*np.exp(-(x-(1730-24))**2/17**2)+np.exp(-(x-1730)**2/17**2)+0.39 *
                 np.exp(-(x-(1730+31))**2/17**2))/25.484
    return f_of_x

# This version had a four-term polynomial, which David didn't want to use for now.
# def blah_spec_fit(x, c0, c1, c2, c3, c4, w1, w2, f1):
#     f_of_x = c0+c1*(x-1700)+c2*(x-1700)**2+c3*(x-1700)**3+c4*(x-1700)**4 + \
#              w1*(0.0751747*np.exp(-(x-1903.82)**2/26.4725**2)+0.225213*np.exp(-(x-1942.21)**2/48.8781**2) +
#                  0.005*np.exp(-(x-1779.71)**2/32.1869**2))/7.715 + \
#              w2*(0.0280945*np.exp(-(x-1913.6)**2/25.0449**2)+0.103527*np.exp(-(x-1949.5)**2/52.2024**2))/3.07 + \
#              f1*(0.31*np.exp(-(x-(1730-24))**2/17**2)+np.exp(-(x-1730)**2/17**2)+0.39 *
#                  np.exp(-(x-(1730+31))**2/17**2))/25.484
#     return f_of_x


def model_fit(df, wave_array):
    '''
    Performs spec_fit function on the input spectra.

    Input parameters:
    df is dataframe containing only the spectra, without extraneous columns
    wave_array is np.array containing the wavelengths of spectra

    Return values:
    modelparams is np.array of parameters output from curve_fit
    modelcovar is np.array of covariance array output from curve_fit
    modeled_spectra is np.array of the modeled spectra using the curve_fit parameters
    residual_spectra is np.array of the original spectrum - the modeled spectrum

    '''
    number = df.shape[0]
    modelparams = np.empty((number, 7), dtype=float)
    modelcovar = np.empty((number, 7, 7), dtype=float)
    modeled_spectra = np.empty((number, wave_array.shape[0]), dtype=float)
    residual_spectra = np.empty((number, wave_array.shape[0]), dtype=float)
    for i in range(number):
        row = df.iloc[i, :]
        modeled, pcov = curve_fit(spec_fit, wave_array, row)
        modelparams[i, :] = modeled
        modelcovar[i, ::] = pcov
        modeled_spectra[i, :] = spec_fit(wave_array, *modeled)
        residual_spectra[i, :] = row - modeled_spectra[i, :]
    return modelparams, modelcovar, modeled_spectra, residual_spectra


def extract_farray(modelparams):
    '''
    Convenience function to create an np.array containing only the f1 values

    Input:  the modelparams np.array created by model_fit().
    '''
    return modelparams[:, 6]


def create_modelparams_df(modelparams):
    columns = ['c0', 'c1', 'c2', 'c3', 'w1', 'w2', 'f1']
    modelparams_df = pd.DataFrame(data=modelparams, columns=columns)
    return modelparams_df


def extract_f_sd_err(modelcovar):
    '''
    Convenience function to create an np.array containing the standard deviation
       of the covariance matrix for the f1 term.

    Input:  the modelcovar np.array created by model_fit().

    '''
    f_sd_err = np.empty(modelcovar.shape[0])
    for i in range(modelcovar.shape[0]):
        f_sd_err[i] = np.diag(modelcovar[i, :, :])[6]
    return f_sd_err


def high_f1_numbers(farray, threshold):
    '''
    Returns np.array of the index positions of spectra with f1 values greater than the
       threshold.  Use this with iloc to find the appropriate rows in the dataframe.

    Input:  np.array output from extract_farray() and desired threshold value.
    '''
    return np.where(farray > threshold)[0]


def get_visible_wavelength_vector(file_path, calibration_path_str):
    with h5py.File(file_path, 'r') as h5_file:
        wavelength_vector = h5_file[calibration_path].attrs['spec1_wavelengths_vector'][:]
        return wavelength_vector


def get_visible_insertion_absorbance_depth_df(file_path, insertion_path_str):
    with h5py.File(file_path, 'r') as h5_file:
        absorbance_depth = h5_file[f"{insertion_path_str}/derived/absorbance_depth"][:]
        calibration_path_str = insertion_path_str[:17]
        wavelength_vector = h5_file[calibration_path_str].attrs['spec1_wavelengths_vector'][:]
        length = wavelength_vector.shape[0]
        columns = list(np.arange(0, length, 1))
        columns.append('force')
        columns.append('depth')
        dataframe = pd.DataFrame(data=absorbance_depth, columns=columns)
        return dataframe


def create_insertions_date_dataframe(spreadsheet, date_selection):
    df = pd.read_csv(spreadsheet)
    selected_df = df.loc[df['date'] == date_selection]
    return selected_df









if __name__ == "__main__":
    path_name = "/Users/linda/OneDrive/Documents/S4_mine_p/Projects/Data_collected/"
    spreadsheet = 'data/nirone_misc_insertions.csv'
    date_selection = '4/15/21'
    f1_threshold = 0.06

    selected_date_df = create_insertions_date_dataframe(spreadsheet, date_selection)

    high_f1_df = pd.DataFrame()
    all_values_df = pd.DataFrame()

    for index in selected_date_df.index.unique():

        # Determine variables for getting the appropriate data
        file_name = selected_date_df['file_name'][index]
        file = path_name + file_name
        calibration_path = selected_date_df['session'][index] + '/' + selected_date_df['calibration'][index]
        insertion_path = calibration_path + '/' + selected_date_df['insertion'][index]

        # Get the wavelength vector
        waves = get_visible_wavelength_vector(file, calibration_path)

        # Get the absorbance-force-depth dataset and put in DataFrame format
        absorbance_depth_df = get_visible_insertion_absorbance_depth_df(file, insertion_path)
        absorbance_depth_df['file_name'] = file_name
        absorbance_depth_df['insertion'] = insertion_path

        # The model input is only the absorbances part of the dataset
        absorbances = absorbance_depth_df.iloc[:, :-4]

        # Run the model
        modelparams, modelcovar, modeled_spectra, residual_spectra = model_fit(absorbances, waves)

        # Calculate the standard deviation of the covariance matrix for the f1 term
        sd_err = extract_f_sd_err(modelcovar)
        absorbance_depth_df['sd_err'] = sd_err

        # Get the f1 array from the model parameters
        farray = extract_farray(modelparams)
        # absorbance_depth_df['f1_value'] = farray
        absorbance_depth_df[['c0', 'c1', 'c2', 'c3', 'w1', 'w2', 'f1']] = modelparams

        # Determine which of the f1 numbers meet the desired threshold
        ins_high_f1s = high_f1_numbers(farray, f1_threshold)

        # Find the spectra that meet the high f1 threshold
        high_spectra = absorbance_depth_df.iloc[ins_high_f1s, :].copy()

        # Append the spectra to the dataframes
        high_f1_df = pd.concat([high_f1_df, high_spectra], ignore_index=True, sort=False)
        all_values_df = pd.concat([all_values_df, absorbance_depth_df], ignore_index=True, sort=False)

    print(absorbance_depth_df.shape)
    print(absorbance_depth_df.iloc[0, :])
    print(high_f1_df[['insertion', 'depth', 'f1']])
    print()
    print(all_values_df[['insertion', 'depth', 'sd_err']])






















