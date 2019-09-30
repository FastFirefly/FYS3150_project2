from  matplotlib import pyplot as plt
import numpy as np
from scipy.stats import norm

fout = "output_d_R"
figfile = fout+".jpeg"
data = np.loadtxt(fout)
data2 = np.loadtxt("output_d_A")
N = sum(1 for line in open(fout))
min_i = 100
h = 10/(N+1);
hh = h*h;
rho = np.zeros(200)
print(data2[3][3])
for i in range(0, 200):
	rho[i] = (i+1)*h;
	if data2[i][i] < min_i:
		min_i = data2[i][i]
		k = i

vec = data[:, k]
mag =  np.dot(vec,vec)
vec = vec/mag

plt.plot(rho, (abs(vec))**2, label=r"$n=%i$"%k)

fout = "output_b_R"
figfile = fout+".jpeg"
data = np.loadtxt(fout)
data2 = np.loadtxt("output_b_A")
min_i = 100
h = 10/(N+1);
hh = h*h;
rho = np.zeros(200)
print(data2[3][3])
for i in range(0, 200):
	rho[i] = (i+1)*h;
	if data2[i][i] < min_i:
		min_i = data2[i][i]
		k = i

vec = data[:, k]
mag =  np.dot(vec,vec)
vec = vec/mag

plt.plot(rho, (abs(vec))**2, label=r"$n=%i$"%k)

plt.xlabel(r'$\rho$')
plt.ylabel(r'$|\psi(\rho)|^2$')
plt.legend([r"With $\rho^2$", "Without"])
plt.title(r"Quantum dots 3D, one electron")
plt.savefig("2b_and_d")