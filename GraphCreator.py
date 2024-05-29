import matplotlib.pyplot as plt
import numpy as np

def f(x):
    a = np.where(x < 45, 0.01, 0.003)
    return a * (x - 45)**2 + 36

x_values = np.linspace(0, 80) 

y_values = f(x_values)

plt.plot(x_values, y_values)
plt.title('Graph representing the relationship between speed and fuel efficiency')
plt.xlabel('Speed (Mph)')
plt.ylabel('Fuel Efficiency (Mpg)')
plt.grid(True)
plt.show()