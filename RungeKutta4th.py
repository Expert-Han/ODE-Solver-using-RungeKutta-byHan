import numpy as np
from matplotlib import pyplot as plt

# The system dynamics is given as that of the pdf file.
def plant(x, u):
    dx[0] = x[1]
    dx[1] = 980 * x[0] - 2.8 * x[2]
    dx[2] = -100 * x[2] + 100 * u
    return dx

def rk5(x, u, T):
    k1=plant(x,u)*T
    k2=plant(x+k1*0.5,u)*T
    k3=plant(x+k2*0.5,u)*T
    k4=plant(x+k3,u)*T
    dx = x + ((k1+k4)/6+(k2+k3)/3)
    return dx

n = 3
K = [[365.46613439, 21.77699344, 0.49931265]]
x0 = np.array([0.1, -0.1, -0.1], dtype=np.float64)
u0 = np.matmul(K,x0.T)
dx = np.array([0,0,0], dtype=np.float64)

T = 0.001
tf = 1
sam = int(tf / T)
tspan = np.linspace(0, tf, sam + 1)

xs = len(tspan)
x = np.array([x0], dtype=np.float64)
u_sig = np.array([u0], dtype=np.float64)

for i in range(0, xs - 1):
    u = np.matmul(K, x[i].T)
    x_next = rk5(x[i], u, T)
    x = np.vstack((x,x_next))
    if i >= 1:
        u_sig = np.vstack((u_sig,u))

u = np.matmul(K, x[i+1].T)
u_sig = np.vstack((u_sig,u))

plt.figure()
plt.plot(tspan, x[:, 0], label="x1")
plt.plot(tspan, x[:, 1], label="x2")
plt.plot(tspan, x[:, 2], label="x3")
plt.grid()
plt.xlabel("Time")
plt.ylabel("State response")
plt.legend()
plt.show()

plt.figure()
plt.plot(tspan, u_sig[:, 0], label="u")
plt.grid()
plt.xlabel("Time")
plt.ylabel("Control signal")
plt.legend()
plt.show()
