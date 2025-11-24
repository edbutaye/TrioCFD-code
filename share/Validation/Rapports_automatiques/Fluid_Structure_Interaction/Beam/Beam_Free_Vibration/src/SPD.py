import numpy as np


# -----------------------------
# Paramètres fichier
# -----------------------------
file1 = "Poutre_Displacement1D_plane_0.out"

# -----------------------------
# Lecture signal
# -----------------------------
t = np.loadtxt(file1, usecols=[0])
x = np.loadtxt(file1, usecols=[1])

# -----------------------------
# Vérifier si le signal est uniformément échantillonné
dt = np.diff(t)
if np.max(dt) - np.min(dt) > 1e-12:
    print("Attention : le signal n'est pas parfaitement uniforme, interpolation recommandée")
    # Ici on interpole le signal pour uniformiser
    t_uniform = np.linspace(t[0], t[-1], len(t))
    x_uniform = np.interp(t_uniform, t, x)
else:
    t_uniform = t
    x_uniform = x

# -----------------------------
# FFT
# -----------------------------
N = len(x_uniform)
fs = 1.0 / np.median(np.diff(t_uniform))  # fréquence d'échantillonnage
Xf = np.fft.fft(x_uniform)
freqs = np.fft.fftfreq(N, d=1/fs)

# On prend uniquement les fréquences positives
pos_mask = freqs > 0
freqs = freqs[pos_mask]
power = np.abs(Xf[pos_mask])**2

# Fréquence dominante
idx_max = np.argmax(power)
freq_dom = freqs[idx_max]

freq= round(freq_dom,2)
fichier = open("Freq.txt", "w")       
fichier.write(str(freq))        
fichier.close()


