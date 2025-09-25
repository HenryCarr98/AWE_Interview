import sys
import matplotlib.pyplot as plt

xs = []
ys = []
with open(sys.argv[1]) as f:
    lines = f.readlines()
    for line in lines:
        x, y = [float(v) for v in line.split(",")]
        xs.append(x)
        ys.append(y)

plt.plot(xs, ys)
plt.show()
