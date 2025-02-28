import numpy as np
import matplotlib.pyplot as plt

# FROM CHAT GPT

# Define range of lambda*h values (real and imaginary parts)
real_values = np.linspace(-5, 2, 400)  # Real axis
imag_values = np.linspace(-3, 3, 400)  # Imaginary axis

# Create a meshgrid for complex numbers
Re, Im = np.meshgrid(real_values, imag_values)
Lambda_h = Re + 1j * Im  # Complex lambda*h

# Compute stability function Q(lambda*h)
Q = 1 / (1 - Lambda_h)

# Compute magnitude |Q|
Q_magnitude = np.abs(Q)

# Plot contour levels where |Q| = 1 (stability boundary)
plt.figure(figsize=(8, 6))
contour = plt.contour(Re, Im, Q_magnitude, levels=[1], colors='red', linewidths=2)
plt.clabel(contour, fmt="|Q|=1", colors='red')

# Additional contours to visualize |Q|
plt.contourf(Re, Im, Q_magnitude, levels=np.linspace(0, 2, 30), cmap='coolwarm')

# Labels and title
plt.xlabel(r'Re$(\lambda h)$')
plt.ylabel(r'Im$(\lambda h)$')
plt.title(r'Stability Region of Implicit Euler: $Q(\lambda h) = \frac{1}{1 - \lambda h}$')
plt.colorbar(label=r'$|Q(\lambda h)|$')
plt.axhline(0, color='black', linewidth=0.5)
plt.axvline(0, color='black', linewidth=0.5)

# Show the plot
plt.show()