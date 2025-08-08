import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftshift
from scipy.signal import find_peaks
#from radar_system import Radar, Target, Signal, Filter, RadarPlotter, CFAR
from radar_system import *

B = float(input('Unesite bandwidth ( GHz): ')) * 1e9  # Bandwidth u Hz
N = int(input('Unesite broj chirp signala: '))  # Broj chirp signala
Tc = float(input('Unesite period chirp signala (us): ')) *1e-6
number_of_targets = int(input('Unesite broj meta: '))

P_t = float(input('Unesite snagu transmitovanog signala (dBm): '))
G = float(input('Unesite Gain antene (dBi): '))
rcs = float(input('Unesite Radar Cross Section (m^2): '))

SNR_min = 20
objects = []
signals = []
reflected = []
radar = Radar(B, N, Tc, P_t, G, rcs)
radar.print_radar_limitations()
# Unos razdaljine i brzine za svaku metu

for i in range(number_of_targets):
    R = float(input(f'Unesi razdaljinu za objekat {i+1} (u metrima): '))
    V = float(input(f'Unesi brzinu za objekat {i+1} (u m/s): '))
    objects.append(Target(R, V))
    signals.append(Signal(radar,objects[i]))
    reflected.append(signals[i].generate_signal())


x_combined_matrix = np.zeros((len(signals[0].t), radar.N))
x_final_matrix = np.zeros((signals[0].new_length(), N), dtype=complex)
X_f_complex = np.zeros((signals[0].new_length(), N), dtype=complex)
filtering  = Filter(radar)

for i  in range(N):
    for j in range(number_of_targets):
        x_combined_matrix[:,i] = x_combined_matrix[:,i] + reflected[j][:,i]

for i in range(N):
    x_final_matrix[:,i] =Radar.Radar_VGA_gain* (filtering.filter_signal(signals[0].x_cos* x_combined_matrix[:,i] + 1j*signals[0].x_sin * x_combined_matrix[:,i]))[::radar.fs_factor()]
    x_final_matrix[:, i] = Signal.generate_quantized_signal(x_final_matrix[:,i])
    X_f_complex[:,i] = fft(x_final_matrix[:, i])

f = np.linspace(0, radar.Fs,signals[0].new_length())
L = (radar.Tc * Radar.c * f/2) / (2 * radar.B)


plotting =RadarPlotter()
plotting.plot_grafics(L,(np.abs(X_f_complex[:, N//2])), 'Razdaljina između objekta i radara','Daljina (m)','Amplituda[dB]')

cfar = CFAR(number_of_ref_cells = 12, number_of_gap_cells = 8, bias=12)
threshold_values, targets_only ,peaks = cfar.calculate_cfar_threshold((np.abs(X_f_complex[:, N//2])))
CFAR.plot_threshold(L,np.abs( X_f_complex[:, N//2]), threshold_values)

#DOPPLER FFT

range_doppler_map = RangeDoppler(radar)
range_doppler_map.plot_range_doppler(L, X_f_complex)
print("hello")
print("mama ti se kupa gola ")
print("i tebi tata")


