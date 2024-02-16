
from toolbox import secante, bissection
import matplotlib.pyplot as plt

def f(x):
    return x**2 - 2*x - 3

# Plotting Graph
x = [i-5 for i in range(11)]
y = [f(x-5) for x in range(11)]

plt.plot(x, y)
plt.grid(True)
plt.show()

# Calculating roots of f
x, statut = bissection(f, 4, 1, 0.01)
print(x, statut)
