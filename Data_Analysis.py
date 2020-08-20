#!/usr/bin/env python
# coding: utf-8

# In[1]:


import plotly as py
import numpy as np
import math as m
py.offline.init_notebook_mode(connected=True)
import plotly.graph_objs as go


# # Compatibilità

# In[2]:


def compatibilità(x, errore_x, y, errore_y):
    num = abs(x-y)
    den = m.sqrt(errore_x**2 + errore_y**2)
    r = num/den
    if r >= 0 and r < 1:
        print ("Buona, r = ", round(r,3))
    elif r > 1 and r < 2:
        print("Sufficiente r = ",  round(r,3))
    elif r > 2 and r < 3:
        print ("Scarsa r = ",  round(r,3))
    else:
        print ("Incompatibilità r = ",  round(r,3))
    return r


# In[3]:


def compatibilità_list(x:list, errore_x:list):
    i, j = 0, 0
    while i < len(x)-1:
        print("La compatibilità tra l'elemento ", i, " e ", i+1, " è: ")
        compatibilità(x[i],errore_x[i], x[i+1],errore_x[i+1])
        i += 1
    while j < len(x)-2:
        print("La compatibilità tra l'elemento ", j, " e ", j+1, " è: ")
        compatibilità(x[j], errore_x[j], x[j+2], errore_x[j+2])
        j += 1


# 
# # Media Ponderata

# In[4]:


def media_ponderata(x:list,errore_x:list):
    num = 0
    den = 0
    for i in range(0,len(x)):
        num += x[i]/pow(errore_x[i],2)
        den += 1/pow(errore_x[i],2)
    media = num/den
    errore_media = 1/m.sqrt(den)
    if len(x) == 2:
        r = compatibiltà(x[0],errore_x[0], x[1], errore_x[1])
        if r <3 :
            print("Poiché c'è compatibilità ne facciamo la media ponderata:")
    print("La media ponderata vale: ", media)
    print("La sue incertezza è:", errore_media)
    return media, errore_media


# # Coefficiente di Pearson

# In[5]:


def coefficiente_Pearson(x:list, y:list):
    sum1, sum2, sum_quadro_1, sum_quadro_2 = 0, 0, 0, 0
    for i in range(len(x)):
        sum1 += x[i]
    media_x = sum1/len(x)
    for j in range(len(y)):
        sum2 += y[j]
    media_y = sum2/len(y)
    for k in range(len(x)):
        sum_quadro_1 += pow(x[k]-media_x,2) 
        sum_quadro_2 += pow(y[k]-media_y,2)
    numeratore = 0
    denominatore = m.sqrt(sum_quadro_1)*m.sqrt(sum_quadro_2)
    for l in range(len(x)):
        numeratore += (x[l]-media_x)*(y[l]-media_y)
    rho = numeratore/denominatore
    print("L'indice di correlazione di Pearson vale: $\rho = $", rho)
    return rho


# # INTERPOLAZIONE LINEARE  $y = a + bx $
# 

# In[6]:


def interpolazione(y:list, errore_y:list, x:list, errore_x = None):
    uno_sigma_quadro1, x_su_sigma_quadro1, x_per_y_su_sigma_quadro1, y_su_sigma_quadro1, x_quadro_su_sigma_quadro1 = 0, 0, 0, 0, 0

    for i in range (0, len(x)):
        uno_sigma_quadro1 += 1/pow(errore_y[i],2)
        x_su_sigma_quadro1 += x[i]/pow(errore_y[i],2)    
        x_quadro_su_sigma_quadro1  += pow(x[i],2)/pow(errore_y[i],2) 
        x_per_y_su_sigma_quadro1 += x[i]*y[i]/pow(errore_y[i],2)
        y_su_sigma_quadro1 += y[i]/pow(errore_y[i],2)

    delta = uno_sigma_quadro1*x_quadro_su_sigma_quadro1-pow(x_su_sigma_quadro1,2)
    b = 1/delta * (uno_sigma_quadro1*x_per_y_su_sigma_quadro1 - x_su_sigma_quadro1*y_su_sigma_quadro1)
    a = 1/delta * (x_quadro_su_sigma_quadro1*y_su_sigma_quadro1 - x_su_sigma_quadro1*x_per_y_su_sigma_quadro1)

    print("Il coefficiente angolare 'b' vale: ", b)
    errore_b = m.sqrt(1/delta * uno_sigma_quadro1)
    print("L'incertezza sul coefficiente angolare vale: ", errore_b)
    errore_a = m.sqrt(1/delta * x_quadro_su_sigma_quadro1)
    print("L'intercetta 'a' della retta vale: ", a)
    print("L'incertezza sull'intercetta vale: ", errore_a)
    return b, errore_b, a, errore_a


# In[7]:


# velocità = [0.12501857418816506, 0.3377796936499026 , 0.9786705738924276 , 5.173305742369384, 5.261613816678994 , 0.1405670474694918, 0.6069139638764814, 10.511562718990888 , 0.22527052576531698,  1.5775275935868251, 2.671602770833734]
# #errore_velocità = [0.00006409509095846773, 0.00023353098399002224, 0.00044779853853757436 , 0.020596272008780083, 0.026006174551478824 , 0.00007262138712827439, 0.0005581732544674304 , 0.034082021804363004 , 0.00029703316783571374, 0.0018420509717760717 , 0.015661958100380263]
# errore_velocità = [0.000036670628445924206, 0.00010164401284265719, 0.0003669667339277374, 0.02766855314564073,  0.026006174551478824, 0.000014444489758834329, 0.0003406037386506163, 0.034082021804363004, 0.00009246674574684704, 0.0008918906094775717, 0.003744536276137677]
# v_inf = [0.13047902831277786, 0.3612326547392903, 1.0924069791328825 , 6.139864188406879 , 6.919876735275608, 0.14707625288413445 , 0.663315298393586, 12.723557449194109 , 0.23840358287407704, 1.7977149912447228, 3.1079169706569005 ]
# velocità2 = [0.12501857418816506, 0.3377796936499026 , 0.9786705738924276 , 5.173305742369384, 0.1405670474694918, 0.6069139638764814, 0.22527052576531698,  1.5775275935868251, 2.671602770833734]
# #errore_velocità2 = [0.00006409509095846773, 0.00023353098399002224, 0.00044779853853757436 , 0.020596272008780083, 0.00007262138712827439, 0.0005581732544674304 , 0.00029703316783571374, 0.0018420509717760717 , 0.015661958100380263]
# errore_velocità2 = [0.000036670628445924206, 0.00010164401284265719, 0.0003669667339277374, 0.02766855314564073, 0.000014444489758834329, 0.0003406037386506163, 0.00009246674574684704, 0.0008918906094775717, 0.003744536276137677]
# v_inf2 = [0.13047902831277786, 0.3612326547392903, 1.0924069791328825 , 6.139864188406879 , 0.14707625288413445 , 0.663315298393586, 0.23840358287407704, 1.7977149912447228, 3.1079169706569005 ]
         
# diametri = [0.15, 0.2381, 0.3969, 0.635, 0.635, 0.1588, 0.3175, 0.7144, 0.20, 0.4763, 0.5556] #cm
# diametri2 = [0.15, 0.2381, 0.3969, 0.635, 0.1588, 0.3175, 0.20, 0.4763, 0.5556] #cm
# raggi_quadro = []
# raggi_quadro2 = []
# for i in range(len(diametri)):
#     raggi_quadro.append(pow((diametri[i]/2),2))
# for i in range(len(diametri2)):
#     raggi_quadro2.append(pow((diametri2[i]/2),2))
# print(len(velocità), len(errore_velocità), len(v_inf), len(velocità2), len(errore_velocità2), len(v_inf2))
# print(raggi_quadro)


# In[8]:


#interpolazione(velocità2, errore_velocità2, raggi_quadro2)


# In[9]:


def interpolazione_no_errore(x:list, y:list, errore_y):
    N = len(y)
    x_quadri, x_singoli, y_singoli, x_per_y = 0, 0, 0, 0
    for i in range (0,len(y)):
        x_quadri += pow(x[i],2)
        x_singoli += x[i]
        y_singoli += y[i]
        x_per_y += y[i]*x[i]

    delta = N*x_quadri -pow(x_singoli,2)
    a = (x_quadri*y_singoli-x_singoli*x_per_y)/delta
    b = (N*x_per_y-x_singoli*y_singoli)/delta
    errore_a = errore_y*m.sqrt(x_quadri/delta)
    errore_b = errore_y*m.sqrt(N/delta)
    print("Questo è a, l'intercetta della retta: ", a,"\nQuesto il suo errore:", errore_a, "\nQuesto è b, il coefficiente angolare, cioè la velocità 1/v: ", b, "\nQuesto il suo errore:", errore_b)
    return (b, errore_b, a, errore_a)


# In[10]:


# h = [10, 15, 20, 25, 30, 35, 40, 45] #cm
# data = [29.404999999999998, 44.175, 59.092, 73.72800000000001, 88.546, 103.334, 118.21800000000002, 133.07]
# errore_data = 0.1*m.sqrt(2)/m.sqrt(24)
# interpolazione_no_errore(h,data,errore_data)


# 
# # Grafico

# In[11]:


def grafico(b, a, y, errore_y, x, errore_x = None):
    
    layout = go.Layout(
        title= str(input("Inserisci il nome del grafico ") ),
        yaxis=dict(
                title = str("$"+str(input("Inserisci il titolo dell'asse y ") ) + " [" +str(input("Inserisci l'unità di misura ")+ "]" + "$"))
        ),
        xaxis=dict(
            title= str("$" + str(input("Inserisci il titolo dell'asse x ") ) + " [" +str(input("Inserisci l'unità di misura "))+ "]" +"$"),
        )
    )
    i = int(input("Scegli di quanto vuoi arrotondare b: ") )
    j = int(input("Scegli di quanto vuoi arrotondare a: ") )
    traccia = go.Scatter(
        x = x,
        y= y,
        mode='markers',
        name= "$" + str(input("Inserisci la legenda ") ) +"$",
        showlegend=True,
       # line = dict(
        #    shape='spline'
       # ),
        error_y=dict(
                type='data',
                array = errore_y,
        ),
        error_x=dict(
                type='data',
                array = errore_x,
        ),
    )
#     traccia2 = go.Scatter(
#         x = raggi_quadro2,
#         y= velocità2,
#         mode='markers',
#         name='$v_i, i \in [1,9]$',
#         showlegend=True,
#        # line = dict(
#         #    shape='spline'
#        # )

#             error_y=dict(
#                 type='data',
#                 array = errore_velocità2,
#             )
#     )
    t = np.linspace(0., max(x), 1000)
    lalla = a + b*t
    if a >= 0:
        retta = go.Scatter(
            x = t,
            y = lalla,
            mode = 'lines',
            name = str("y = " + str(round(b,i))+ "x +" + str(round(a,j)) + str( input("Inserisci l'unità di misura: ")))
            #name = '$y = (25.594x - 0.02081) cm/s $',
    #             line = dict(
    #             shape='spline',
    #             color = 'orange'
    #        )  
        )
    else:
        retta = go.Scatter(
            x = t,
            y = lalla,
            mode = 'lines',
            name = str("$y = " + "("+str(round(b,i))+ "x" + str(round(a,j)) +")" + str( input("Inserisci l'unità di misura: "))+"$")
            #name = '$y = (25.594x - 0.02081) cm/s $',
    #             line = dict(
    #             shape='spline',
    #             color = 'orange'
    #        )  
        )
    fig = go.Figure(data=[traccia, retta], layout = layout)
    py.offline.iplot(fig)


# In[ ]:





# In[12]:


#z = interpolazione(velocità2, errore_velocità2, raggi_quadro2, None)


# In[13]:


#grafico(z[0], z[2], velocità2, errore_velocità2, raggi_quadro2)


# In[14]:


def chi_square(b, a, y:list, errore_y:list):
    #vincoli = N-2 -- caso interpolazione
    chi_quadro = 0
    for i in range(0,len(x)):
        chi_quadro += pow((y[i]-(b*x[i]+a)),2)/pow(errore_y[i],2)
    print("Il chi quadro vale:", chi_quadro)
    print("Il numero di DEF nel caso di 2 vincoli è: ", len(y)-2)
    return chi_quadro


# In[15]:


def errore_posteriori(b, a, y:list, x:list):
    somma = 0
    for i in range(0, len(y)):
        somma += pow(a+b*x[i]-y[i],2)
    errore_a_posteriori = somma/(len(y)-2)
    print("L'errore a posteriori dell'interpolazione vale: ", errore_a_posteriori)
    return errore_a_posteriori
    
    


# # Andamento temporale

# In[16]:


def andamento_temporale(b, a, y:list, errore_y, x:list):
    scarti = []
    #vincoli = N-2 == 9
    for i in range(0,len(y)):
        scarti.append((y[i]-(b*x[i]+a))/errore_y)
    layout = go.Layout(
        title= "",
        yaxis=dict(
                title = 'Scarti'
        ),
        xaxis=dict(
            title= str("$" + str(input("Inserisci il titolo dell'asse x: ")) + "$")
        )
    )
    l = np.linspace(0,99,100)
    traccia = go.Scatter(
        x = l,
        y = scarti,
        mode = 'lines',
        name = 'Scarti',
        #y = 7.98x -1.1
        line = dict(
        #shape='spline',
        color = str(input("Inserisci un colore (default #EF553B): "))
                   ),
        )
#     print(traccia)   
#     print(type(traccia))
    l1 = np.linspace(0., len(y), 1000)
    retta = go.Scatter(
        x = l1,
        y = [0]*len(l1),
        mode = 'lines',
        name = "$y = 0$",
        line = dict(
            shape='spline',
            color = 'firebrick',
            width=5
            )
        )
   
    #fig = go.Figure(data=[traccia1, traccia2, traccia3, traccia4, traccia5, traccia6, traccia7, traccia8, traccia9, retta], layout = layout)
    fig = go.Figure(data=[traccia, retta], layout=layout)
    fig.show()


# In[17]:


# h = [10, 15, 20, 25, 30, 35, 40, 45] #cm
# data = [29.404999999999998, 44.175, 59.092, 73.72800000000001, 88.546, 103.334, 118.21800000000002, 133.07]
# errore_data = 0.1*m.sqrt(2)/m.sqrt(24)
# andamento_temporale(2.96, -0.22, data, errore_data, h)


# In[23]:


def traccia(y, errore_y, x, errore_x = None):
    #l = np.linspace(0,99,100)
    traccia1 = go.Scatter(
        x = x,
        y = y,
        mode = str(input("Inserisci il mode per la traccia: lines, markers, markers+lines ")),
        name= "$" + str(input("Inserisci la legenda della traccia ") ) +"$",
        showlegend=True,
        error_y=dict(
                type='data',
                array = errore_y,
        ),
        error_x=dict(
                type='data',
                array = errore_x,
        ),
    )
    return traccia1
def retta(b, a, y, errore_y, x, errore_x = None):
    i = int(input("Scegli di quanto vuoi arrotondare b: ") )
    j = int(input("Scegli di quanto vuoi arrotondare a: ") )
    t = np.linspace(0., max(x), 1000)
    lalla = a + b*t
    if a >= 0:
        retta1 = go.Scatter(
            x = t,
            y = lalla,
            mode = 'lines',
            name = str("y = " + str(round(b,i))+ "x +" + str(round(a,j)) + str( input("Inserisci l'unità di misura: ")))
            #name = '$y = (25.594x - 0.02081) cm/s $',
    #             line = dict(
    #             shape='spline',
    #             color = 'orange'
    #        )  
        )
    else:
        retta1 = go.Scatter(
            x = t,
            y = lalla,
            mode = 'lines',
            name = str("$y = " + "("+str(round(b,i))+ "x" + str(round(a,j)) +")" + str( input("Inserisci l'unità di misura della legenda: "))+"$")
            #name = '$y = (25.594x - 0.02081) cm/s $',
    #             line = dict(
    #             shape='spline',
    #             color = 'orange'
    #        )  
        )
        return retta1
def Layout():
    layout = go.Layout(
        title= str(input("Inserisci il nome del grafico ") ),
        yaxis=dict(
                title = str("$"+str(input("Inserisci il titolo dell'asse y ") ) + " [" +str(input("Inserisci l'unità di misura ")+ "]" + "$"))
        ),
        xaxis=dict(
            title= str("$" + str(input("Inserisci il titolo dell'asse x ") ) + " [" +str(input("Inserisci l'unità di misura "))+ "]" +"$"),
        )
    )
    return layout
    


# In[25]:


def Grafico(x:list, layout):
    fig = go.Figure(data=x, layout=layout)
    fig.show()


# In[20]:


# layout = Layout()
# traccia1 = traccia(velocità2, errore_velocità2, raggi_quadro2)
# retta1 = retta(z[0], z[2], velocità2, errore_velocità2, raggi_quadro2) 
# x = [traccia1, retta1]
# Grafico(x)

