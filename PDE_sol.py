import numpy as np
import matplotlib.pyplot as plt
import time

# П06-1. Решить уравнение теплопроводности (u_{th}(x,t)- точное решение)

def  he2_step(x, t, u, h, dt):
    nn = len(u)
    q_int = lambda x, t:  (x + 2*t**2*(np.tanh(x*t)))/(np.cosh(x*t)**2)
    x1 = 0.0;  x2 =  1.0; 
    hs = h*h  
    aL, bL, fL = 1.0, 1.0, 1 + t
    aR, bR, fR = 1.0, 0.0, 1 + t / np.cosh(t)**2
    U1= u[1];    UN2= u[nn-2]
     # Реализация явной схемы для внутр.точек
    unew = np.zeros(nn)    
    # for j in range(1,nn-1):
    #     unew[j] = u[j] + dt*( (u[j+1] - 2.0*u[j] + u[j-1])/hs + q_int(x[j],t))

    unew[1:nn-1] = u[1:nn-1] + dt*((u[2:nn] - 2.0*u[1:nn-1] +  u[0:nn-2])/hs  + q_int(x[1:nn-1], t))

    # Граничные условия 1-го рода:
    #unew[0] = fL
    #unew[nn-1]= fR

    # Граничные условия 2-го или 3-го родов:
    uL= U1 - 2.0*h*(fL - bL*u[0]) / aL
    # unew[0] = u[0] + dt*( (U1 - 2.0*u[0] + uL)/hs + q_int(x1,t))
    uR= UN2 + 2.0*h*(fR - bR*u[nn-1]) / aR
    unew[nn-1] = u[nn-1] + dt*( (uR - 2.0*u[nn-1] + UN2)/hs + q_int(x2,t))

    # ***************
    return unew

# Точное решение
u_th = lambda x, t: x + np.tanh(x*t)

# Входные параметры
h = 0.01
x = np.arange(0, 1, h)
t0 = 0.0  # начальное время
tend = 5 # конечное время
dt = 5.0e-5
nout = 1000   # nout*h - шаг вывода результатов
nstep = int(np.round(tend / dt))  # полное количество шагов

u0 =  x   # начальное условие

tm0 = time.time()
t = t0;  u = np.array(u0, copy = True)
# Инициализация данных для вывода
tdata = t;   udata = np.array(u, copy = True)
# Интегрирование НУТ по t = t0:tend
for ks in range(nstep):
# ***** Введите ваш код: *****
    u = he2_step(x, t, u, h, dt)
    t = t + dt

# ***************
    # Если iter кратно nout, то сохранять результаты
    if (ks +1) % nout == 0:
        tdata = np.append(tdata, t)
        udata = np.vstack((udata, u))

tm1 = time.time()

# Построить графики функций координаты и скорости
# ***** Введите ваш код: *****
# nt, nx = udata.shape
plt.plot(x, udata[0, :], 'o', x, u_th(x, 0))
# plt.plot(x, udata[-1, :], 'x', x, u_th(x, tend))
# ***************
print(tm1-tm0)
plt.show()
