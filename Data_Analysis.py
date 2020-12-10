#!/usr/bin/env python
# coding: utf-8

# In[1]:


import plotly as py
import numpy as np
import math as m
from scipy.odr import *
py.offline.init_notebook_mode(connected=True)
import plotly.graph_objs as go


# # Compatibilità

# In[2]:


def compatibilità(x, errore_x, y, errore_y, docs:bool = False):
    """Ritorna il valore compatibilità o una stringa a seconda che docs sia False o True; Compatibilità di default
    """
    num = abs(x-y)
    den = m.sqrt(errore_x**2 + errore_y**2)
    r = num/den
    output = ""
    if r >= 0 and r < 1:
        output += f"Buona, r = {round(r,3)}\n"
    elif r > 1 and r < 2:
        output += f"Sufficiente r = {round(r,3)}\n"
    elif r > 2 and r < 3:
        output += f"Scarsa r = {round(r,3)}\n"
    else:
        output += f"Incompatibilità r = {round(r,3)}\n"
    if docs:
        return output
    return r


# In[3]:


def compatibilità_list(x:list, errore_x:list, docs:bool = False):
    """Ritorna una stringa
    """
    output = ""
    for i, val1 in enumerate(x):
        for j, val2 in enumerate(x):
            if i == j or i-j > 0:
                pass
            else:
                output += f"La compatibilità tra l'elemento {i+1} e {j+1} è:"
                comp = compatibilità(val1,errore_x[i], val2,errore_x[j], True)
                output += comp
                output += " "
                print("La compatibilità tra l'elemento", {i+1}, "e ", {j+1}, " è:", comp)
    # while i < len(x)-1:
    #     output += f"La compatibilità tra l'elemento {i+1} e {i+2} è: \n"
    #     output += compatibilità(x[i],errore_x[i], x[i+1],errore_x[i+1])
    #     i += 1
    # while j < len(x)-2:
    #     output += f"La compatibilità tra l'elemento {j+1} e {j+2} è: \n"
    #     output += compatibilità(x[j], errore_x[j], x[j+2], errore_x[j+2])
    #     j += 1
    if docs: 
        return output


# 
# # Media Ponderata

# In[4]:


def media_ponderata(x:list, errore_x:list, docs:bool = False):
    """Ritorna un tuple: media, errore_media, stringa per il documento
    """
    num = 0
    den = 0
    output = ""
    for i in range(0,len(x)):
        num += x[i]/pow(errore_x[i],2)
        den += 1/pow(errore_x[i],2)
    media = num/den
    errore_media = 1/m.sqrt(den)
    if len(x) == 2:
        r = compatibilità(x[0],errore_x[0], x[1], errore_x[1])
        if r <3 :
            output += f"Poiché c'è compatibilità ne facciamo la media ponderata:\n"
    output += f"La media ponderata vale: {media}\n"
    output += f"La sue incertezza è: {errore_media}\n"
    print("La media ponderata vale: ", media)
    print("La sue incertezza è: ", errore_media)
    if docs:
        return output
    return media, errore_media


# # Coefficiente di Pearson

# In[5]:


def coefficiente_Pearson(x:list, y:list, docs:bool = False):
    sum1, sum2, sum_quadro_1, sum_quadro_2 = 0, 0, 0, 0
    output = ""
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
    output += f"L'indice di correlazione di Pearson vale: $\rho = $ {rho}\n"
    print("L'indice di correlazione di Pearson vale: rho = ", rho)
    if docs:
        return output
    return rho


# # Covarianza

# In[6]:


def covarianza(x:list, y:list):
    N = len(x)
    N_covariance = 0
    x_mean = sum(x)/N
    y_mean = sum(y)/N
    for i in range(N):
        N_covariance += (x[i] - x_mean)*(y[i] - y_mean)
    return N_covariance/N


# # INTERPOLAZIONE LINEARE  $y = a + bx $
# 

# In[7]:


def interpolazione(y:list, errore_y:list, x:list, errore_x = None, docs:bool = False):
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
    output = ""
    print("Il coefficiente angolare 'b' vale: ", b)
    errore_b = m.sqrt(1/delta * uno_sigma_quadro1)
    print("L'incertezza sul coefficiente angolare vale: ", errore_b)
    errore_a = m.sqrt(1/delta * x_quadro_su_sigma_quadro1)
    print("L'intercetta 'a' della retta vale: ", a)
    print("L'incertezza sull'intercetta vale: ", errore_a)
    output += f"Il coefficiente angolare 'b' vale: {b}\n"
    output += f"L'incertezza sul coefficiente angolare vale:: {errore_b}\n"
    output += f"IL'intercetta 'a' della retta vale: {a}\n"
    output += f"L'incertezza sull'intercetta vale: {errore_a}\n"
    if docs:
        return output
    return b, errore_b, a, errore_a


# In[8]:


def interpolazione_no_errore(x:list, y:list, errore_y, docs:bool = False):
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
    print("Questo è a, l'intercetta della retta: ", a,"\nQuesto il suo errore:", errore_a, "\nQuesto è b, il coefficiente angolare", b, "\nQuesto il suo errore:", errore_b)
    output = ""
    output += f"Questo è a, l'intercetta della retta: {a}\nQuesto il suo errore: {errore_a} \nQuesto è b, il coefficiente angolare {b} \nQuesto il suo errore: {errore_b}"
    if docs:
        return output
    return (b, errore_b, a, errore_a)


# 
# # Grafico

# In[9]:


def grafico(b, a, y, errore_y, x, errore_x = None):
    
    layout = go.Layout(
        title= "",
        yaxis=dict(
                title = "y"
        ),
        xaxis=dict(
            title= "x",
        )
    )
    traccia = go.Scatter(
        x = x,
        y= y,
        mode='markers',
        name= "Traccia",
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
            name = "Retta"
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
            name = "Retta"
            #name = '$y = (25.594x - 0.02081) cm/s $',
    #             line = dict(
    #             shape='spline',
    #             color = 'orange'
    #        )  
        )
    fig = go.Figure(data=[traccia, retta], layout = layout)
    py.offline.iplot(fig)


# In[10]:


def grafico_relazione(b, a, y, errore_y, x, errore_x = None):
    
    layout = go.Layout(
        title= str(input("Inserisci il nome del grafico ") ),
        yaxis=dict(
                title = str("$"+str(input("Inserisci il titolo dell'asse y ") ) + " [" + str(input("Inserisci l'unità di misura ")) + "]" + "$")
        ),
        xaxis=dict(
            title= str("$" + str(input("Inserisci il titolo dell'asse x ") ) + " [" +str(input("Inserisci l'unità di misura ")) + "]" +"$"),
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
            name = str("$y = " + "(" + str(round(b,i))+ "x +" + str(round(a,j))+ ")" + str(input("Inserisci l'unità di misura ")) +"$" )
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
            name = str("$y = " + "("+str(round(b,i))+ "x" + str(round(a,j)) +")" + str(input("Inserisci l'unità di misura ")) +"$")
            #name = '$y = (25.594x - 0.02081) cm/s $',
    #             line = dict(
    #             shape='spline',
    #             color = 'orange'
    #        )  
        )
    fig = go.Figure(data=[traccia, retta], layout = layout)
    py.offline.iplot(fig)


# In[ ]:





# In[11]:


#z = interpolazione(velocità2, errore_velocità2, raggi_quadro2, None)


# In[12]:


#grafico(z[0], z[2], velocità2, errore_velocità2, raggi_quadro2)


# In[13]:


def chi_square(b, a, y:list, errore_y:list, x:list, docs:bool = False):
    #vincoli = N-2 -- caso interpolazione
    chi_quadro = 0
    output = ""
    for i in range(0,len(x)):
        chi_quadro += pow((y[i]-(b*x[i]+a)),2)/pow(errore_y[i],2)
    print("Il chi quadro vale:", chi_quadro)
    print("Il numero di DOF nel caso di 2 vincoli è: ", len(y)-2)
    output += f"Il chi quadro vale: {chi_quadro}"
    output += f"Il numero di DOF nel caso di 2 vincoli è: {len(y)-2}"
    if docs:
        return output
    return chi_quadro


# In[14]:


def errore_posteriori(b, a, y:list, x:list, docs:bool = False):
    somma = 0
    output = ""
    for i in range(0, len(y)):
        somma += pow(a+b*x[i]-y[i],2)
    errore_a_posteriori = somma/(len(y)-2)
    print("L'errore a posteriori dell'interpolazione vale: ", errore_a_posteriori)
    output += f"L'errore a posteriori dell'interpolazione vale: {errore_a_posteriori}"
    if docs:
        return output
    return errore_a_posteriori


# In[15]:


def t_Student(rho, N):
    errore_rho = m.sqrt((1-rho**2)/(N-2))
    t = rho/errore_rho
    print("La variabile di Student vale t =", t, "con ", N-2, "gradi di libertà")


# # Andamento temporale

# In[16]:


def andamento_temporale_relazione(b, a, y:list, errore_y, x:list, x_axis_name = "Occorrenze" ,colore = "#EF553B"):
    scarti = []
    #vincoli = N-2 == 9
    for i in range(0,len(y)):
        scarti.append((y[i]-(b*x[i]+a))/errore_y[i])
    layout = go.Layout(
        title= "",
        yaxis=dict(
                title = 'Scarti'
        ),
        xaxis=dict(
            title= str("$" + x_axis_name + "$")
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
#         color = str(input("Inserisci un colore (default #EF553B): "))
        color = str(colore)
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


def andamento_temporale(b, a, y:list, errore_y:list, x:list):
    scarti = []
    #vincoli = N-2 == 9
    for i in range(0,len(y)):
        scarti.append((y[i]-(b*x[i]+a))/errore_y[i])
    layout = go.Layout(
        title= "",
        yaxis=dict(
                title = 'Scarti'
        ),
        xaxis=dict(
            title= "x"
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
        color = '#EF553B'
                   ),
        )
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
    fig = go.Figure(data=[traccia, retta], layout=layout)
    fig.show()


# In[18]:


def traccia(y, errore_y, x, errore_x = None, modo = 'markers'):
    #l = np.linspace(0,99,100)
    traccia1 = go.Scatter(
        x = x,
        y = y,
        mode = str(modo),
        name= "$x_i$",
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
    max_x, min_x = max(x), min(x)
    t_max = max_x if max_x >= abs(min_x) else min_x
    t = np.linspace(-t_max, t_max, 1000) 
    lalla = a + b*t
    retta1 = go.Scatter(
        x = t,
        y = lalla,
        mode = 'lines',
        name = "Retta"
        #name = '$y = (25.594x - 0.02081) cm/s $',
#             line = dict(
#             shape='spline',
#             color = 'orange'
#        )  
    )

    return retta1

def Layout():
    layout = go.Layout(
        title= "",
        yaxis=dict(
                title = ""
        ),
        xaxis=dict(
            title= "",
        )
    )
    return layout
    


# In[ ]:





# In[19]:


def traccia_relazione(y, errore_y, x, errore_x = None):
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

def retta_relazione(b:float, a:float, y, errore_y, x, errore_x = None):
    i = int(input("Scegli di quanto vuoi arrotondare b: ") )
    j = int(input("Scegli di quanto vuoi arrotondare a: ") )
    max_x, min_x = max(x), min(x)
    t_max = max_x if max_x >= abs(min_x) else min_x
    t = np.linspace(-t_max, t_max, 1000)
    lalla = float(round(a,j)) + float( round(b,i))*t
    if a >= 0:
        retta1 = go.Scatter(
            x = t,
            y = lalla,
            mode = 'lines',
            name = str("$y = " + "("+str(float(round(b,i)))+ "x +" + str(float(round(a,j))) +")" + str( input("Inserisci l'unità di misura della legenda: "))+"$")
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
            name = str("$y = " + "("+str(float(round(b,i)))+ "x" + str(float(round(a,j))) +")" + str( input("Inserisci l'unità di misura della legenda: "))+"$")
            #name = '$y = (25.594x - 0.02081) cm/s $',
    #             line = dict(
    #             shape='spline',
    #             color = 'orange'
    #        )  
        )
    return retta1

def Layout_relazione():
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
    


# In[20]:


def Grafico(x:list, layout):
    fig = go.Figure(data=x, layout=layout)
    fig.show()


# # Fit Lineare

# # $y = bx + a$

# In[21]:


def f(B, x):
    '''Linear function y = m*x + b'''
    # B is a vector of the parameters.
    # x is an array of the current x values.
    # x is in the same format as the x passed to Data or RealData.
    #
    # Return an array in the same format as y passed to Data or RealData.
    return B[0]*x + B[1]
def linear_fit(x, y, sigma_x, sigma_y, b, a):
    mydata = RealData(x, y, sx=sigma_x, sy=sigma_y)
    linear = Model(f)
    myodr = ODR(mydata, linear, beta0=[b, a])
    myoutput = myodr.run()
    myoutput.pprint()
    return myoutput.beta[0], myoutput.sd_beta[0], myoutput.beta[1], myoutput.sd_beta[1]


# In[ ]:





# # Fit Parabolico

# # $y = ax^2 + bx + c$

# In[22]:


def g(B, x):
    '''quadratic function y = a*x**2 + b*x + c'''
    # B is a vector of the parameters.
    # x is an array of the current x values.
    # x is in the same format as the x passed to Data or RealData.
    #
    # Return an array in the same format as y passed to Data or RealData.
    return B[0]*x**2 + B[1]*x + B[2]
def parabolic_fit(x, y, sigma_x, sigma_y, a, b, c):
    mydata = RealData(x, y, sx=sigma_x, sy=sigma_y)
    quadratic = Model(g)
    myodr = ODR(mydata, quadratic, beta0=[a, b, c])
    myoutput = myodr.run()
    myoutput.pprint()
    return myoutput.beta[0], myoutput.sd_beta[0], myoutput.beta[1], myoutput.sd_beta[1], myoutput.beta[2], myoutput.sd_beta[2]


# In[23]:


def max_quadratic(a, b, c):
    x_max = -b/(2*a)
    y_max = a*x_max**2 + b*x_max +c
    return (y_max, x_max)


# In[ ]:





# # Scrittura su file

# In[24]:


def scrittura(nome_file, contenuto):
    ciao = nome_file
    f = open(ciao, "w")
    f.write(contenuto)
    #f.read()
    #f.readlines()
    f.close()
    f = open(nome_file, "a")
    a_capo = "\n"
    f.write(a_capo)
    f.close()
def scrittura_aggiunta(nome_file, contenuto):
    f = open(nome_file, "a")
    f.write(contenuto)
    a_capo = "\n"
    f.write(a_capo)
    f.close()


# # Importare da Excel

# In[25]:


import openpyxl
def excel_import(nome_file:str, start, stop, col):
    #excel_document = openpyxl.load_workbook('Guidovia.xlsx')
    #sheet = excel_document.get_sheet_by_name('Foglio1')
    nome_file += '.xlsx'
    wb = openpyxl.load_workbook(nome_file)
    sheets = wb.sheetnames
    ws = wb[sheets[0]]
    vettore = [] #array contenete le medie delle misure effettuate della variabile diretta.
    for i in range(start, stop + 1):
        casella = ws.cell(row=i, column=col)
        vettore.append(casella.value) #aggiunge elementi al vettore (array)
    return vettore


# # Leggi file

# In[26]:


def text(nome_file:str):
    h = open(nome_file, 'r')
    # Reading from the file
    content = h.readlines()
    # print(content)
    h.close()
    n = 0
    x = []
    y = []
    for line in content:
        if n==0:
            N = line[0]
        if n>0:
            words = line.split()
            if len(words):
                x.append(words[0])
                y.append(words[1])
        n+=1
    return(x,y)

