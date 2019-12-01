#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np


# In[2]:


E1 = 50
E2 = 25
E3 = 0

R1 = 10
R2 = 3
R3 = 8
R4 = 6
R5 = 7
R6 = 0

L1 = 0.15
L2 = 0.03
L3 = 0
L4 = 0.8
L5 = 0.09
L6 = 0.1

# * 10^-6
C1 = 36
C2 = 20
C3 = 15
C4 = 0
C5 = 10
C6 = 9

f = 50
w = 2 * 3.14 * f


# In[3]:


Xl_mtrx = np.zeros(6)
Xc_mtrx = np.zeros(6)
X_mtrx = np.zeros(6, dtype=complex)
Z_mtrx = np.zeros(6, dtype=complex)
kirchhoff_mtrx = np.zeros((6, 6), dtype=complex)


# In[4]:


C_mtrx = [C1, C2, C3, C4, C5, C6]
L_mtrx = [L1, L2, L3, L4, L5, L6]
R_mtrx = [R1, R2, R3, R4, R5, R6]
E_mtrx = [0, 0, 0, E1, -E2, (E1 + E2)]


# In[5]:


for x in range(0, 6):
    Xl_mtrx[x] = w * L_mtrx[x]


# In[6]:


for x in range(0, 6):
    if C_mtrx[x] != 0:
        Xc_mtrx[x] = ((1 / (w * C_mtrx[x])) * 1000000)
    else:
            Xc_mtrx[x] = 0


# In[7]:


for x in range(0, 6):
    X_mtrx[x] = (Xl_mtrx[x] - Xc_mtrx[x]) * 1j


# In[8]:


for x in range(0, 6):
    Z_mtrx[x] = R_mtrx[x] + X_mtrx[x]


# In[9]:


kirchhoff_mtrx = [[1, -1, -1, 0, 0, 0],
                  [0, 0, 1, 1, 0, -1],
                  [0, 1, 0, -1, -1, 0],
                  [Z_mtrx[0], 0, Z_mtrx[2], -Z_mtrx[3], Z_mtrx[4], 0],
                  [0, 0, 0, -Z_mtrx[3], Z_mtrx[4], -Z_mtrx[5]],
                  [Z_mtrx[0], 0, Z_mtrx[2], 0, 0, Z_mtrx[5]]]


# In[10]:


I_mtrx = np.linalg.solve(kirchhoff_mtrx, E_mtrx)


# In[14]:


# Не понятно... Тоже самое, но цифры другие (https://media1.tenor.com/images/1f5acd2e50e75ad23bc9fd96aad1db73/tenor.gif?itemid=15569009)
I_mtrx = np.matmul(np.linalg.pinv(kirchhoff_mtrx), E_mtrx)


# In[11]:


P_mtrx = R_mtrx * ((np.abs(I_mtrx) ** 2) / np.sqrt(2))


# In[12]:


Ql_mtrx = Xl_mtrx * ((np.abs(I_mtrx) ** 2) / np.sqrt(2))


# In[13]:


Qc_mtrx = Xc_mtrx * ((np.abs(I_mtrx) ** 2) / np.sqrt(2))


# In[14]:


Q_mtrx = Ql_mtrx - Qc_mtrx


# In[15]:


S_mtrx = np.sqrt((Q_mtrx ** 2) + (P_mtrx ** 2))


# In[16]:


Cos_mtrx = P_mtrx / S_mtrx


# In[25]:


print(' ///ИСХОДНЫЕ ДАННЫЕ\\\ ')
print(" ", " Частота сети:", f, sep='\n')
print(" ", " Активное сопротивление:", R_mtrx, sep='\n')
print(" ", " Индуктивность:", L_mtrx, sep='\n')
print(" ", " Емкость:", C_mtrx, ' ', sep='\n')
print(' \\\ИСХОДНЫЕ ДАННЫЕ/// \n\n\n')

print(" ", "Результаты: ")
print(" ", " Индуктивное сопротивление:", Xl_mtrx, sep='\n')
print(" ", " Eмкостное сопротивление:", Xc_mtrx, sep='\n')
print(" ", " Полное сопротивление:", Z_mtrx, sep='\n')
print(" ", " ЭДС:", E_mtrx, sep='\n')
print(" ", " Коэффициенты Киркхоффа:", kirchhoff_mtrx, sep='\n')
print(" ", " Токи:", I_mtrx, sep='\n')
print(" ", " Реактивная мощность на индуктивности:", Ql_mtrx, sep='\n')
print(" ", " Реактивная мощность на емкости:", Qc_mtrx, sep='\n')
print(" ", " Активная мощность:", P_mtrx, sep='\n')
print(" ", " Реактивная мощность:", Q_mtrx, sep='\n')
print(" ", " Полная мощность:", S_mtrx, sep='\n')
print(" ", " Коэффициенты мощности:", Cos_mtrx, sep='\n')


# In[ ]:





# In[ ]:




