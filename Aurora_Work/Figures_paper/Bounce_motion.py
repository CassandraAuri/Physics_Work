import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Constants
q_e = -1.6e-19  # Electron charge (C)
q_p = 1.6e-19   # Proton charge (C)
m_e = 9.11e-31  # Electron mass (kg)
m_p = 1.67e-27  # Proton mass (kg)
B0 = 1e-4       # Magnetic field strength at equator (T)
R_E = 6.371e6   # Earth radius (m)
v_perp = 1e6    # Perpendicular velocity (m/s)
v_parallel = 1e6 # Parallel velocity (m/s)

# Time settings
time_steps = 5000
dt = 1e-1

def dipole_field(x, z):
    r = np.sqrt(x**2 + z**2)
    Bx = -3 * x * z / r**5
    Bz = (2 * z**2 - x**2) / r**5
    return Bx, Bz

def update_position(x, z, vx, vz, q, m):
    Bx, Bz = dipole_field(x, z)
    B_mag = np.sqrt(Bx**2 + Bz**2)
    
    # Lorentz force
    ax = q * (vz * Bz - vx * B_mag) / m
    az = q * (-vx * Bx) / m
    
    # Update velocity
    vx += ax * dt
    vz += az * dt
    
    # Update position
    x += vx * dt
    z += vz * dt
    
    return x, z, vx, vz

def simulate_particle(q, m):
    x, z = 5 * R_E, 5 * R_E  # Initial position (start near equator)
    vx, vz = 0, v_parallel
    
    # Initial perpendicular velocity
    Bx, Bz = dipole_field(x, z)
    B_mag = np.sqrt(Bx**2 + Bz**2)
    vx = v_perp * Bz / B_mag
    vz = v_perp * Bx / B_mag
    
    positions = []
    for _ in range(time_steps):
        x, z, vx, vz = update_position(x, z, vx, vz, q, m)
        positions.append((x / R_E, z / R_E))  # Normalize to Earth radii
    return np.array(positions)

# Simulate electron and proton
electron_positions = simulate_particle(q_e, m_e)
proton_positions = simulate_particle(q_p, m_p)
print(proton_positions)

# Create the animation
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_xlim(-10, 10)  # Normalize to Earth radii
ax.set_ylim(-10, 10)  # Normalize to Earth radii
ax.set_aspect('equal')
ax.set_xlabel("x (Earth radii)")
ax.set_ylabel("z (Earth radii)")

# Draw Earth
earth = plt.Circle((0, 0), 1, color='blue', fill=True, alpha=0.3)  # Radius = 1 Earth radius
ax.add_artist(earth)

# Initialize plot elements
electron_line, = ax.plot([], [], lw=2, label="Electron", color='red')
proton_line, = ax.plot([], [], lw=2, label="Proton", color='green')

# Update function for animation
def update(frame):
    electron_line.set_data(electron_positions[:frame, 0], electron_positions[:frame, 1])
    proton_line.set_data(proton_positions[:frame, 0], proton_positions[:frame, 1])
    return electron_line, proton_line

# Create the animation
ani = FuncAnimation(fig, update, frames=time_steps, interval=20, blit=True)
ax.legend()
plt.show()
