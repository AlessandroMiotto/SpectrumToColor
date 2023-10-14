import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps


def main():
    illuminant_file = 'd65.csv'
    spectra_name = 'chlorophyll_b.txt'

    illuminant = np.loadtxt('illuminant/'+illuminant_file, delimiter=',')
    spectra = np.loadtxt('spectra/'+spectra_name, delimiter='\t')
    spectra_normalization(spectra)

    #print('Illuminant: ', illuminant_file)
    #print('Molecule: ', spectra_name)

    # computing max density
    density = 1.0 * gamut_edge(illuminant, spectra)
    #print('Fictious density: ', density)

    #XYZ_val = XYZ(illuminant, spectra, density)
    #RGB_val = XYZ_to_sRGB(XYZ_val)

    #print("\nXYZ: ", XYZ_val)
    #print("RGB: ", RGB_val*255)

    
# ---------------------------------------------------------------------------------

# transform XYZ to sRGB color space with D65 illuminant
# http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
def XYZ_to_sRGB(XYZ):
    matrix = np.array([[3.2404542, -1.5371385, -0.4985314],
                       [-0.9692660,  1.8760108,  0.0415560],
                       [0.0556434, -0.2040259, 1.0572252]])

    RGB = np.matmul(matrix, XYZ)

    for i in range(3):
        if RGB[i] > 0.0031308:
            RGB[i] = 1.055 * RGB[i]**(1/2.4) - 0.055
        else:
            RGB[i] = 12.92 * RGB[i]

        if RGB[i] < 0:
            RGB[i] = 0
        if RGB[i] > 1:
            RGB[i] = 1

    return RGB

# calculate XYZ color space from a absorption spectra
def XYZ(illuminant, spectra, density):
    # ILLUMINANT DATA
    l_data = illuminant[:, 0]    # wavelenghts of illuminant data
    d65_data = illuminant[:, 1]  # relative spectral power distribution

    # SPECTRAL DATA
    l_spectra_data = spectra[:, 0]           # wavelenghts of spectra data
    absorption_spectra_data = spectra[:, 1]  # absorption coefficients

    # Slice array to have same size
    # Find wavelenghts domain
    l_max = min(l_data.max(), l_spectra_data.max())
    l_min = max(l_data.min(), l_spectra_data.min())
    # Create a boolean mask to slice arrays
    mask_illuminant = (l_data >= l_min) & (l_data <= l_max)
    mask_spectra = (l_spectra_data >= l_min) & (l_spectra_data <= l_max)
    # Slice arrays of illuminant and spectra
    d65 = d65_data[mask_illuminant]
    l = l_data[mask_illuminant]
    absorption_spectra = absorption_spectra_data[mask_spectra]
    l_spectra = l_spectra_data[mask_spectra]

    # interpolate array to match values
    absorption_spectra = np.interp(l,l_spectra,absorption_spectra)

    # calculation of color matching functions
    xyz = xyz_functions(l)
    x = xyz[0]
    y = xyz[1]
    z = xyz[2]

    N = simps(d65*y) # normalization constant

    # XYZ integrals (transmissive case)
    X = np.trapz(np.exp(-absorption_spectra*density)*d65*x)/N
    Y = np.trapz(np.exp(-absorption_spectra*density)*d65*y)/N
    Z = np.trapz(np.exp(-absorption_spectra*density)*d65*z)/N

    XYZ = np.array([X,Y,Z])

    return XYZ


# CREATION OF xyz COLOR MATCHING FUNCTIONS
def xyz_functions(l):
    x = np.zeros(len(l))
    y = np.zeros(len(l))
    z = np.zeros(len(l))

    # gaussian function
    def g(x, mu, sigma_1, sigma_2):
        if x < mu:
            return np.exp(-0.5*(x-mu)**2/sigma_1**2)
        else:
            return np.exp(-0.5*(x-mu)**2/sigma_2**2)

    # calculation of xyz using experimental data about the observer chromatic response
    for i in range(len(l)):
        x[i] = 1.056*g(l[i], 599.8, 37.9, 31.0) + 0.362*g(l[i], 442.0, 16.0, 26.7) - 0.065*g(l[i], 501.1, 20.4, 26.2)
        y[i] = 0.821*g(l[i], 568.8, 46.9, 40.5) + 0.286*g(l[i], 530.9, 16.3, 31.1)
        z[i] = 1.217*g(l[i], 437.0, 11.8, 36.0) + 0.681*g(l[i], 459.0, 26.0, 13.8)

    xyz = np.array([x, y, z])

    return xyz

# Rough normalization of spectral data
def spectra_normalization(spectra):
    mask = (spectra[:,0] >= 300) & (spectra[:,0] <= 830)
    norm = simps(spectra[:,1][mask])
    spectra[:,1] = spectra[:,1] / norm


def gamut_edge(illuminant, spectra):
    density = 0.0
    while True:
        density += 10.0
        XYZ_val = XYZ(illuminant, spectra, density)
        RGB_val = XYZ_to_sRGB(XYZ_val)
        print(str(RGB_val[0]*255)+'\t'+str(RGB_val[1]*255)+'\t'+str(RGB_val[2]*255)+'\n')
        if RGB_val[0] == 0 or RGB_val[1] == 0 or RGB_val[2] == 0:
            break

    return density


if __name__ == "__main__":
    main()