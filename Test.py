#!/usr/bin/env python
# coding: utf-8

# In[1]:


#import data_analysis.ipynb 
import import_ipynb
from Data_Analysis import *


# In[2]:


x = [1, 2 , 2.5, 1.8, 1.2, 3.6]
errore_x = [0.5, 1, 0.1, 0.5, 0.4, 0.3]
a = compatibilità_list(x,errore_x)
print(a)


# In[3]:


velocità = [0.12501857418816506, 0.3377796936499026 , 0.9786705738924276 , 5.173305742369384, 5.261613816678994 , 0.1405670474694918, 0.6069139638764814, 10.511562718990888 , 0.22527052576531698,  1.5775275935868251, 2.671602770833734]
#errore_velocità = [0.00006409509095846773, 0.00023353098399002224, 0.00044779853853757436 , 0.020596272008780083, 0.026006174551478824 , 0.00007262138712827439, 0.0005581732544674304 , 0.034082021804363004 , 0.00029703316783571374, 0.0018420509717760717 , 0.015661958100380263]
errore_velocità = [0.000036670628445924206, 0.00010164401284265719, 0.0003669667339277374, 0.02766855314564073,  0.026006174551478824, 0.000014444489758834329, 0.0003406037386506163, 0.034082021804363004, 0.00009246674574684704, 0.0008918906094775717, 0.003744536276137677]
v_inf = [0.13047902831277786, 0.3612326547392903, 1.0924069791328825 , 6.139864188406879 , 6.919876735275608, 0.14707625288413445 , 0.663315298393586, 12.723557449194109 , 0.23840358287407704, 1.7977149912447228, 3.1079169706569005 ]
velocità2 = [0.12501857418816506, 0.3377796936499026 , 0.9786705738924276 , 5.173305742369384, 0.1405670474694918, 0.6069139638764814, 0.22527052576531698,  1.5775275935868251, 2.671602770833734]
#errore_velocità2 = [0.00006409509095846773, 0.00023353098399002224, 0.00044779853853757436 , 0.020596272008780083, 0.00007262138712827439, 0.0005581732544674304 , 0.00029703316783571374, 0.0018420509717760717 , 0.015661958100380263]
errore_velocità2 = [0.000036670628445924206, 0.00010164401284265719, 0.0003669667339277374, 0.02766855314564073, 0.000014444489758834329, 0.0003406037386506163, 0.00009246674574684704, 0.0008918906094775717, 0.003744536276137677]
v_inf2 = [0.13047902831277786, 0.3612326547392903, 1.0924069791328825 , 6.139864188406879 , 0.14707625288413445 , 0.663315298393586, 0.23840358287407704, 1.7977149912447228, 3.1079169706569005 ]
         
diametri = [0.15, 0.2381, 0.3969, 0.635, 0.635, 0.1588, 0.3175, 0.7144, 0.20, 0.4763, 0.5556] #cm
diametri2 = [0.15, 0.2381, 0.3969, 0.635, 0.1588, 0.3175, 0.20, 0.4763, 0.5556] #cm
raggi_quadro = []
raggi_quadro2 = []
for i in range(len(diametri)):
    raggi_quadro.append(pow((diametri[i]/2),2))
for i in range(len(diametri2)):
    raggi_quadro2.append(pow((diametri2[i]/2),2))
print(len(velocità), len(errore_velocità), len(v_inf), len(velocità2), len(errore_velocità2), len(v_inf2))
print(raggi_quadro)


# In[4]:


h = [10, 15, 20, 25, 30, 35, 40, 45] #cm
data = [29.404999999999998, 44.175, 59.092, 73.72800000000001, 88.546, 103.334, 118.21800000000002, 133.07]
errore_data = 0.1*m.sqrt(2)/m.sqrt(24)
g = interpolazione_no_errore(h,data,errore_data)
print(g)


# In[5]:


h = [10, 15, 20, 25, 30, 35, 40, 45] #cm
data = [29.404999999999998, 44.175, 59.092, 73.72800000000001, 88.546, 103.334, 118.21800000000002, 133.07]
errore_data = 0.1*m.sqrt(2)/m.sqrt(24)
andamento_temporale(2.96, -0.22, data, errore_data, h)


# In[6]:


z = interpolazione(velocità2, errore_velocità2, raggi_quadro2, None)
layout = Layout()
traccia1 = traccia(velocità2, errore_velocità2, raggi_quadro2)
retta1 = retta(z[0], z[2], velocità2, errore_velocità2, raggi_quadro2) 
x = [traccia1, retta1]


# In[7]:


g = Grafico(x,layout)


# # Prova per Relazione

# In[8]:


x = [1, 2 , 2.5, 1.8, 1.2, 3.6]
errore_x = [0.5, 1, 0.1, 0.5, 0.4, 0.3]
a = compatibilità_list(x,errore_x)
Prova = "Prova.txt"
scrittura(Prova,a)


# In[9]:


h = [10, 15, 20, 25, 30, 35, 40, 45] #cm
data = [29.404999999999998, 44.175, 59.092, 73.72800000000001, 88.546, 103.334, 118.21800000000002, 133.07]
errore_data = 0.1*m.sqrt(2)/m.sqrt(24)
andamento_temporale_relazione(2.96, -0.22, data, errore_data, h)


# In[10]:


z = interpolazione(velocità2, errore_velocità2, raggi_quadro2, None, False)
layout = Layout_relazione()
traccia1 = traccia_relazione(velocità2, errore_velocità2, raggi_quadro2)
retta1 = retta_relazione(z[0], z[2], velocità2, errore_velocità2, raggi_quadro2) 
p = [traccia1, retta1]
G = Grafico(p,layout)


# In[11]:


#Write on -- #True at the end
z = interpolazione(velocità2, errore_velocità2, raggi_quadro2, None, True)
scrittura_aggiunta(Prova, z)


# In[12]:


r = coefficiente_Pearson(raggi_quadro2, velocità2)
print(r)
R = coefficiente_Pearson(raggi_quadro2, velocità2, True)
scrittura_aggiunta(Prova, R)


# In[21]:


print("ci"\r"bello")

