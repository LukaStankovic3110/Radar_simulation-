import numpy as np
from scipy.signal import butter, filtfilt, find_peaks
from scipy.fft import fft, fftshift
import matplotlib.pyplot as plt
class Target:
    c = 3e8
    def __init__(self, R, V):
        self.R = R
        self.V = V
    def update_distance(self, i, Tc):
        return 2 * (self.R - i * self.V * Tc) / Target.c



class Radar:
    c = 3e8
    lambda_ = c / 59.25e9  # Talasna dužina u m (npr., 10 GHz signal)
    T = 290  # Temperatura
    k = 1.38e-23
    NF = 10 ** (12 / 10)
    SNR_min = 20
    max_value_adc = 0.9 # Maksimalni napon (V)
    min_value_adc = -0.9 # Minimalni napon (V)
    bit_depth = 12
    Radar_VGA_gain = 1e3
    Fs_adc_max = 20e6




    def __init__(self, B, N, Tc ,P_t, G, rcs):
        self.P_t =  10**(P_t / 10) * 10**-3
        self.G = 10**(G / 10)
        self.rcs = rcs
        self.B = B  # Bandwidth
        self.f0 = (60e9 - self.B / 2)
        self.N = N  # Broj chirp signala
        self.Tc = Tc  # Period chirp signala
        self.Lmax_pocetno  = ((self.rcs * self.P_t * G**2 * Radar.lambda_**2 * N * Tc) / (10**(Radar.SNR_min / 10) * (4 * np.pi)**3 * Radar.k * Radar.T * Radar.NF))**(1/4)
        self.Fs_pocetno = 2 * 2 * self.f0 + self.B  # Frekvencija uzorkovanja
        self.fbmax_pocetno = (self.B * 2 * self.Lmax_pocetno) / (self.Tc * self.c)  # Maksimalna frekvencija
        self.Fs = (2.5 * self.fbmax_pocetno) if (2.5 * self.fbmax_pocetno)<Radar.Fs_adc_max else Radar.Fs_adc_max
        self.fbmax = self.Fs/2.5
        self.Lmax = (self.fbmax * self.Tc * self.c) / (2 * self.B)

    def print_radar_limitations(self):
        print(f'Rezolucija distance je: {Radar.c / (2 * self.B)} m')
        print(f'Maksimalna brzina je: {Radar.c / (4 * (60e9 - self.B / 2) * self.Tc)} m/s')
        print(f'Rezolucija brzine je: {Radar.c / (2 * (60e9 - self.B / 2) * self.N * self.Tc)} m/s')
        print(f'Maksimalna distanca koja ce se prikazivati je: {self.Lmax} m')

    def received_signal_amplitude(self, R):
        return  np.sqrt(2*(self.P_t * self.G ** 2 * Radar.lambda_ ** 2 * self.rcs) / ((4 * np.pi) ** 2 * R ** 4))

    def transmited_signal_amplitude(self):
        return  np.sqrt(2* self.P_t)

    def calculate_snr(self, R):
        snr_db = 10 * np.log10((self.rcs * self.P_t * self.G ** 2 * Radar.lambda_ ** 2 * self.N * self.Tc) / ((4 * np.pi) ** 3 * R ** 4 * Radar.k * Radar.T * Radar.NF))
        print(f"Izracunat SNR za distancu {R} je : {snr_db} ")

    def noise_amplitude(self):
        return Radar.k * Radar.T * (self.Fs/2)

    def fs_factor(self):
        return round(self.Fs_pocetno/self.Fs)

    def calculate_threshold(self, objects):
        # Početna vrednost praga
        prag = float('inf')

        # Iteriramo kroz sve objekte i računamo prag za svakog
        for obj in objects:
            # Proračun praga za trenutni objekat
            prag_current = 0.7 * ((1 / 2) * self.Fs_pocetno * self.received_signal_amplitude(obj.R) * self.transmited_signal_amplitude() * self.Tc / self.fs_factor())

            # Ažuriramo minimum praga
            prag =Radar.Radar_VGA_gain* min(prag, prag_current)

        # Ispisujemo rezultat
        print(f'Prag razlikovanja signala i suma je: {20 * np.log10(prag) }')

        return prag





class Signal:
    def __init__(self, radar, target):
        self.radar = radar
        self.target = target  # Lista objekata
        self.t = np.arange(0, self.radar.Tc, 1 / self.radar.Fs_pocetno)
        self.noise = np.sqrt(self.radar.noise_amplitude()) * np.random.randn(len(self.t))
        self.x_cos = self.radar.transmited_signal_amplitude() * np.cos(2 * np.pi * (self.radar.f0 + self.radar.B * self.t / self.radar.Tc) * self.t)
        self.x_sin = self.radar.transmited_signal_amplitude() * np.sin(2 * np.pi * (self.radar.f0 + self.radar.B * self.t / self.radar.Tc) * self.t)
        self.window = np.hanning(len(self.x_cos))


    def generate_signal(self):

        x_matrix = np.zeros((len(self.x_cos), self.radar.N))

        for i in range(self.radar.N):
            x_matrix[:, i] = self.radar.received_signal_amplitude(self.target.R) * np.cos(2 * np.pi * (self.radar.f0 + self.radar.B * (self.t - self.target.update_distance( i, self.radar.Tc)) / self.radar.Tc) * (self.t - self.target.update_distance( i, self.radar.Tc)))
            x_matrix[:, i] +=  self.noise
            x_matrix[:, i] *= self.window


        return  x_matrix

    def new_length(self):
        new_length = round(len(self.x_cos) / self.radar.fs_factor())
        if new_length < len(self.x_cos[::self.radar.fs_factor()]):
            new_length += 1
        return new_length

    @staticmethod
    def generate_quantized_signal( signal):
        levels = 2 ** Radar.bit_depth  # Broj nivoa u ADC-u
        quantization_step = (Radar.max_value_adc - Radar.min_value_adc) / (levels - 1)  # Korak kvantizacije
        quantized_signal = np.round((signal - Radar.min_value_adc) / quantization_step) * quantization_step + Radar.min_value_adc
        return  quantized_signal

class Filter:
    def __init__(self, radar):
        self.lowpass_cutoff = radar.f0*0.5
        self.nyquist = 0.5 * radar.Fs_pocetno
        self.normal_cutoff = self.lowpass_cutoff / self.nyquist

    def filter_signal(self, signal):
        [b,a]= butter(2, self.normal_cutoff, btype='low')
        return filtfilt(b,a,signal)




class RadarPlotter:
    def __init__(self):
        pass

    def plot_grafics(self, L, x, title, xlabel, ylabel):
        plt.figure()
        plt.plot(L, 20 * np.log10(x))
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.grid(True)  # uključivanje interaktivnog režima
        plt.show()


class CFAR:
    def __init__(self,number_of_ref_cells,number_of_gap_cells,bias) :
        self.number_of_ref_cells = number_of_ref_cells
        self.number_of_gap_cells = number_of_gap_cells
        self.bias = bias

    def calculate_cfar_threshold(self, signal_fft):
        N = signal_fft.size
        cfar_values = np.zeros(signal_fft.shape)
        for center_index in range(
                self.number_of_gap_cells +  self.number_of_ref_cells, N - (self.number_of_gap_cells +  self.number_of_ref_cells)
        ):
            min_index = center_index - (self.number_of_gap_cells +  self.number_of_ref_cells)
            min_guard = center_index - self.number_of_gap_cells
            max_index = center_index + (self.number_of_gap_cells +  self.number_of_ref_cells) + 1
            max_guard = center_index + self.number_of_gap_cells + 1

            lower_nearby = signal_fft[min_index:min_guard]
            upper_nearby = signal_fft[max_guard:max_index]

            lower_mean = np.mean(lower_nearby)
            upper_mean = np.mean(upper_nearby)

            mean = np.mean(np.concatenate((lower_nearby, upper_nearby)))

            output = mean * self.bias
            cfar_values[center_index] = output

        first_non_zero_from_start_index = np.argmax(cfar_values != 0)
        cfar_values[0:first_non_zero_from_start_index] = cfar_values[first_non_zero_from_start_index]
        first_non_zero_from_end_index = np.argmax(cfar_values[::-1] != 0)
        first_non_zero_from_end_index = len(cfar_values) - first_non_zero_from_end_index - 1
        cfar_values[first_non_zero_from_end_index+1:] = cfar_values[first_non_zero_from_start_index]



        targets_only = np.copy(signal_fft)

        targets_only[np.where(signal_fft < cfar_values)] = np.ma.masked
        # Koristi find_peaks za detekciju lokalnih maksimuma
        peaks = find_peaks(targets_only,distance=5)

        return cfar_values, targets_only ,peaks

    @staticmethod
    def plot_threshold(L, signal, threshold_values):
        plt.figure()
        plt.plot(L, 20 * np.log10(signal))
        plt.plot(L, 20 * np.log10(threshold_values))
        plt.grid()
        plt.show()

class RangeDoppler:
    def __init__(self, radar):
        self.radar = radar
        self.freq_axis = np.linspace(-np.pi, np.pi, radar.N)
        self.x_axis = -(Radar.c * self.freq_axis) / (4 * np.pi * radar.f0 * radar.Tc)

    def generate_doppler_fft(self ,X_f_complex):
        """
        Generisanje Doppler FFT i priprema podataka za Range-Doppler spektar.
        """
        # Kreiranje niza sa kompleksnim brojevima
        doppler_fft = np.zeros((self.radar.N, len(X_f_complex[:, 0])), dtype=complex)

        # Računanje FFT za svaki red matrice X_f_complex
        for i in range(len(X_f_complex[:, 0])):
            doppler_fft[:, i] = fftshift(fft(X_f_complex[i, :]))

        return doppler_fft

    def plot_range_doppler(self, L, X_f_complex):

        doppler_fft= self.generate_doppler_fft(X_f_complex)

        # Postavljanje granica za X i Y osovine
        extent = [-self.x_axis.min(), -self.x_axis.max(), L.min(), L.max()]

        # Kreiranje figure i prikazivanje Range-Doppler spektra
        fig, ax = plt.subplots(figsize=(8, 8))
        range_doppler_plot = ax.imshow(
            20 * np.log10(np.abs(doppler_fft.T)),  # Prikaz u dB skali
            aspect='auto',  # Omogućava automatsko prilagođavanje proporcija slike
            extent=extent,  # Postavljanje opsega osovina (X: brzina, Y: range)
            origin='lower',  # Početak koordinatnog sistema u donjem levom kutu
            cmap='viridis',  # Colormap (moguće je koristiti i 'jet', 'gray', ...)
            vmax=25,  # Maksimalna vrednost za colormap
            vmin=-25  # Minimalna vrednost za colormap
        )

        # Dodavanje boje (colorbar)
        fig.colorbar(range_doppler_plot)

        # Postavljanje osovina i naslova
        ax.set_title("Range Doppler Spectrum", fontsize=24)
        ax.set_xlabel("Velocity (m/s)", fontsize=22)
        ax.set_ylabel("Range (m)", fontsize=22)

        # Prikazivanje grafika
        plt.show()