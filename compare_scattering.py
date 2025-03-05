#! python3.11
"""
Comparing the scattering of two calculations
and outputs the percentage difference.

Patrick Wang
patrick.wang@chem.ox.ac.uk

05/03/2025
"""

import numpy as np
import matplotlib.pyplot as plt
import colorcet as cc
import seaborn as sns
import sys
import scipy.integrate as spi

sns.set_palette(cc.glasbey_bw)

def integrate(y, dx):
    #  return np.sum(y) * dx
    return spi.trapezoid(y, dx=dx)

def calc_pdf(data):
    """
    Calculate the PDF of the data using the sine transform.

    Parameters
    ----------
    data : np.ndarray
        The data to calculate the PDF of. Can be directly read from a file using np.loadtxt.
        Atomic units are assumed.

    Returns
    -------
    r : np.ndarray
        The r values of the PDF.
    pdf : np.ndarray
        The PDF values.
    V_ee_int : np.ndarray
        The integrand of the electron-electron interaction energy.
    
    """
    q = data[:, 0]
    dq = data[1, 0] - data[0, 0]
    no_points = len(q)
    r = np.linspace(0, 1 / dq, no_points)
    pdf = np.zeros_like(r)
    for i, rr in enumerate(r):
        integrand = data[:, 1] * np.sin(rr * q) * q
        pdf[i] = 1 / np.pi * rr * integrate(integrand, dq)
    integral = integrate(pdf, r[1] - r[0])
    print(f"Integral = {integral:10.8f}")

    V_ee_int = np.zeros_like(pdf)
    V_ee_int[1:] = pdf[1:] / r[1:]
    V_ee = integrate(V_ee_int, r[1] - r[0])
    print(f"V_ee     = {V_ee:10.8f}")
    return r, pdf, V_ee_int


def read_file(file):
    return np.loadtxt(file)

def sine_transform(data):
    return np.fft.fft(data)

def main():
    '''
    Reads two scattering files and plots the percentage difference between the two.
    Also plots the PDFs and percentage difference of the PDFs.
    '''
    
    args = sys.argv[1:]
    file1 = args[0]
    file2 = args[1]
    file1_label = args[2]
    file2_label = args[3]

    data1 = read_file(file1)
    data2 = read_file(file2)
    pdf1 = calc_pdf(data1)
    pdf2 = calc_pdf(data2)
    r1, pdf1, V_ee_int1 = pdf1
    r2, pdf2, V_ee_int2 = pdf2

    fig, ax = plt.subplots(2, 2, figsize=(12, 12))
    ax[0,1].plot(r1, pdf1, label=file1_label)
    ax[0,1].plot(r2, pdf2, label=file2_label)
    ax[0,1].set_xlabel("r")
    ax[0,1].set_ylabel("PDF")
    ax[0,1].legend()
    ax[0,1].set_title("PDF")
    ax[0,0].plot(data1[:, 0], data1[:, 1], label=file1_label)
    ax[0,0].plot(data2[:, 0], data2[:, 1], label=file2_label)
    ax[0,0].set_xlabel("Scattering")
    ax[0,0].set_ylabel("Intensity")
    ax[0,0].legend()
    ax[0,0].set_xlabel('q / a.u.')

    # Percentage difference plots
    diff = np.abs(pdf1 - pdf2) / pdf1 * 100
    ax[1,1].plot(r1, diff)
    ax[1,1].set_xlabel("r")
    ax[1,1].set_ylabel("Percentage difference")
    ax[1,1].set_title("Percentage difference in PDF")
    
    diff_scat = np.abs(data1[:, 1] - data2[:, 1]) / data1[:, 1] * 100
    ax[1,0].plot(data1[:, 0], diff_scat)
    ax[1,0].set_xlabel("Scattering")
    ax[1,0].set_ylabel("Percentage difference")
    ax[1,0].set_title("Percentage difference in scattering")
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()


    




    
