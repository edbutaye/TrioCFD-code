import numpy as np
import matplotlib.pyplot as plt

# ============================
# Geometrical and material data
# ============================
cylinder_radius = 5e-3
cylinder_length = 0.250
solid_density = 7850
young_modulus = 2e11
fluid_density = 1000
kinematic_viscosity = 0.884e-6
initial_const_A = 0. 
initial_const_B = 1E-5


# ============================================================
# THEORETICAL REFERENCE
# ------------------------------------------------------------
# "Fluid-Induced Vibration Frequency and Damping of a Coaxial
#  Cylinder in a Quiescent Viscous Medium: Theoretical and
#  Numerical Predictions", JAM 2025, Lagrange & Puscas
#
# Governing displacement form for a clamped-free cylinder:
#
# Tip displacement (at x = L):
#    u_tip(t) = exp(-σ*t) * [ A*sin(ω*t) + B*cos(ω*t) ]
#
# Displacement at any x ∈ [0, L]:
#    u(x,t) = exp(-σ*t) * [ A*sin(ω*t) + B*cos(ω*t) ] * W(x)
#
# where W(x) is the first mode shape (clamped-free beam):
#    W(x) = (1/N) * [ sinh(β*x/L) - sin(β*x/L)
#             - h * (cosh(β*x/L) - cos(β*x/L)) ]
# ============================================================

# ============================================================
# NUMERICAL REFERENCES
# ------------------------------------------------------------
# 1) "Optimizing Coupled Fluid-Structure Simulations for
#     Nuclear-Relevant Geometries"
#    2024, Vivaldi & Ricciardi,
#    Journal of Pressure Vessel Technology
#
# 2) "Assessment of an Euler-Bernoulli Beam Model Coupled to CFD
#     in order to Perform Fluid-Structure Simulations"
#    2022, Vivaldi & Ricciardi, PVP2022-81755
# ============================================================




# ============================
# Mode shape function
# ============================
def mode_shape(x, L):
    N = -2.7
    beta = 1.875
    h = 1.3622
    return (1 / N) * (
        np.sinh(beta * x / L)
        - np.sin(beta * x / L)
        - h * (np.cosh(beta * x / L) - np.cos(beta * x / L))
    )

# ============================
# Displacement computation
# ============================
def compute_displacement(freq_hz, log_dec_percent, time, A, B, x, L):
    omega = freq_hz * 2 * np.pi
    sigma = omega * log_dec_percent / 100
    W = mode_shape(x, L)
    return np.exp(-sigma * time) * (A * np.sin(omega * time) + B * np.cos(omega * time)) * W
    
    
# ============================
# Velocity computation
# ============================
def compute_velocity(freq_hz, log_dec_percent, time, A, B, x, L):
    omega = freq_hz * 2 * np.pi
    sigma = omega * log_dec_percent / 100
    W = mode_shape(x, L)
    
    exp = np.exp(-sigma * time)
    sin = np.sin(omega * time)
    cos = np.cos(omega * time)
    coeff_sin = -(B * omega + sigma * A)
    coeff_cos =  (A * omega - sigma * B)
    
    return exp * (coeff_sin * sin + coeff_cos * cos) * W
    
    # ============================
# Velocity acceleration
# ============================
def compute_acceleration(freq_hz, log_dec_percent, time, A, B, x, L):
    omega = freq_hz * 2 * np.pi
    sigma = omega * log_dec_percent / 100
    W = mode_shape(x, L)
    
    exp = np.exp(-sigma * time)
    sin = np.sin(omega * time)
    cos = np.cos(omega * time)

    C_sin = A*(sigma**2 - omega**2) + 2*sigma*B*omega
    C_cos = B*(sigma**2 - omega**2) - 2*sigma*A*omega

    return exp * (C_sin*sin + C_cos*cos) * W

# ============================
# Plotting functions
# ============================
def plot_modal_deformation(beam_length, filename):
    x = np.linspace(0, beam_length, num=1000)
    w = mode_shape(x, beam_length)

    fig, ax = plt.subplots()
    ax.plot(x, w, linewidth=2, color='black', label="Mode shape")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("Normalized displacement")
    ax.set_title("First deformation mode of a clamped-free beam")
    ax.legend()
    fig.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.show()


def plot_displacements(medium_name, freq_theoretical, logdec_theoretical,     
                       A, B, x, L, filename, out):
    """Generate and save displacement comparison plot at a given x."""
    time = np.loadtxt(filename, usecols=[0])
    trio = np.loadtxt(filename, usecols=[1])

    disp_theory = compute_displacement(freq_theoretical, logdec_theoretical, time, A, B, x, L)



    fig, ax = plt.subplots()
    ax.plot(time, disp_theory, linewidth=2, color='black', label="Theory")
    ax.plot(time, trio, linewidth=2, linestyle='dashed', color='red', label="Trio")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Displacement [m]")
    ax.set_title(f"Displacement ")
    ax.legend()
    ax.set_ylim(-1.2 * np.max(np.abs(disp_theory)), 1.2 * np.max(np.abs(disp_theory)))
    fig.tight_layout()
    plt.savefig(out, dpi=300)


# =============================
def plot_velocity(medium_name, freq_theoretical, logdec_theoretical,     
                       A, B, x, L, filename, out):
    """Generate and save velocity comparison plot at a given x."""
    time = np.loadtxt(filename, usecols=[0])
    trio = np.loadtxt(filename, usecols=[1])

    disp_theory = compute_velocity(freq_theoretical, logdec_theoretical, time, A, B, x, L)



    fig, ax = plt.subplots()
    ax.plot(time, disp_theory, linewidth=2, color='black', label="Theory")
    ax.plot(time, trio, linewidth=2, linestyle='dashed', color='red', label="Trio")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Velocity [m/s]")
    ax.set_title(f"Velocity ")
    ax.legend()
    ax.set_ylim(-1.2 * np.max(np.abs(disp_theory)), 1.2 * np.max(np.abs(disp_theory)))
    fig.tight_layout()
    plt.savefig(out, dpi=300)


# =============================
def plot_acceleration(medium_name, freq_theoretical, logdec_theoretical,     
                       A, B, x, L, filename, out):
    """Generate and save acceleration comparison plot at a given x."""
    time = np.loadtxt(filename, usecols=[0])
    trio = np.loadtxt(filename, usecols=[1])

    disp_theory = compute_acceleration(freq_theoretical, logdec_theoretical, time, A, B, x, L)



    fig, ax = plt.subplots()
    ax.plot(time, disp_theory, linewidth=2, color='black', label="Theory")
    ax.plot(time, trio, linewidth=2, linestyle='dashed', color='red', label="Trio")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Acceleration [m/s2]")
    ax.set_title(f"Acceleration ")
    ax.legend()
    ax.set_ylim(-1.2 * np.max(np.abs(disp_theory)), 1.2 * np.max(np.abs(disp_theory)))
    fig.tight_layout()
    plt.savefig(out, dpi=300)


# =============================
# ----------- AIR -------------
# =============================
plot_displacements(
    medium_name="Air",
    freq_theoretical=113.0,
    logdec_theoretical=0.0,
    A=initial_const_A,
    B=initial_const_B,
    x=cylinder_length,         # tip
    L=cylinder_length,
    filename = "Air/Poutre_Displacement1D_plane_0.out",
    out='Displacement_air.png'
)

plot_velocity(
    medium_name="Air",
    freq_theoretical=113.0,
    logdec_theoretical=0.0,
    A=initial_const_A,
    B=initial_const_B,
    x=cylinder_length,         # tip
    L=cylinder_length,
    filename = "Air/Poutre_Velocity1D_plane_0.out",
    out='Velocity_air.png'
)

plot_acceleration(
    medium_name="Air",
    freq_theoretical=113.0,
    logdec_theoretical=0.0,
    A=initial_const_A,
    B=initial_const_B,
    x=cylinder_length,         # tip
    L=cylinder_length,
    filename = "Air/Poutre_Acceleration1D_plane_0.out",
    out='Acceleration_air.png'
)

# =============================
# ----------- WATER -----------
# =============================
plot_displacements(
    medium_name="Water",
    freq_theoretical=107.24,
    logdec_theoretical=0.1,
    A=initial_const_A,
    B=initial_const_B,
    x=cylinder_length,         # tip
    L=cylinder_length,
    filename = "Water/Poutre_Displacement1D_plane_0.out",
    out='Displacement_water.png'
)

plot_velocity(
    medium_name="Water",
    freq_theoretical=107.24,
    logdec_theoretical=0.1,
    A=initial_const_A,
    B=initial_const_B,
    x=cylinder_length,         # tip
    L=cylinder_length,
    filename = "Water/Poutre_Velocity1D_plane_0.out",
    out='Velocity_water.png'
)

plot_acceleration(
    medium_name="Water",
    freq_theoretical=107.24,
    logdec_theoretical=0.1,
    A=initial_const_A,
    B=initial_const_B,
    x=cylinder_length,         # tip
    L=cylinder_length,
    filename = "Water/Poutre_Acceleration1D_plane_0.out",
    out='Acceleration_water.png'
)


