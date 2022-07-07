import numpy as np
import matplotlib.pyplot as plt
plt.style.use("ggplot")

top_data_f = open("top_data.txt", "r")
xs = np.array(list(map(float, top_data_f.readline().split())))
Z = np.array(list(map(float, top_data_f.readline().split())))

plt.plot(xs, Z)
plt.show()