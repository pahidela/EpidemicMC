import matplotlib.pyplot as plt

filename = "SF_N500_g27_mu9.txt"
with open(filename, "r") as f:
    lines = f.readlines()
    beta = [float(x.replace("\n","").split(" ")[0]) for x in lines]
    rho = [float(x.replace("\n","").split(" ")[1]) for x in lines]
    mmca = [float(x.replace("\n","").split(" ")[3]) for x in lines]
plt.plot(beta, mmca, '-g', label = r'MMCA')
plt.plot(beta, rho, '.g', label = r'MC')
plt.xlabel(r"$\beta$")
plt.ylabel(r"$\rho$")
plt.title(r'$\mu$ = 0.9')
plt.legend()
plt.show()
