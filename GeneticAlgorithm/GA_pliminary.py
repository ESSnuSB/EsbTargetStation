import numpy as np
1 of 266

￼

￼

Print all

In new window

GA_my_Version

￼

vidhya hariharan <vidhya.ph.96@gmail.com>

Nov 16, 2023, 3:04 PM (4 days ago)

￼

￼

to vidhya.thara.hariharan

￼

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 09:52:56 2023

@author: vidhya
"""

import seaborn as sns
from scipy.optimize import curve_fit
import math
import numpy as np
import pandas as pd
import os
import random as rn
import subprocess as sp
import matplotlib.pyplot as plt
from scipy.stats import moyal
from deap import base
from deap import creator
from deap import tools

from numpy import linspace


from matplotlib.ticker import NullFormatter, MaxNLocator

import tensorflow as tf
import tensorflow_probability as tfp

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.optimize import curve_fit

from matplotlib import pylab
from pylab import *

tfd = tfp.distributions

if os.path.exists('./ga.txt'):
  os.remove("ga.txt")
  print("ga.txt file removed")
else:
  fn=open("ga.txt","w")
  fn.close()

# The provided Fitness class is an abstract class that needs a weights attribute in order to be functional.
# A maximizing fitness has positive weights. The following line creates, in the creator,
# a ready to use single objective maximizing fitness named FitnessMax.
# The create() function takes at least two arguments, a name for the newly created class and a base class.
# Any subsequent argument becomes an attribute of the class.
# weights attribute must be a tuple so that multi-objective and single objective fitnesses can be treated the same way.
# The first individual created will be a simple list containing floats.
# In order to produce this kind of individual, we need to create an Individual class,
# using the creator, that will inherit from the standard list type and have a fitness attribute.

creator.create("FitnessMax", base.Fitness, weights=(1.0,))
# creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
creator.create("Individual", list, fitness=creator.FitnessMax)

toolbox = base.Toolbox()


# Attribute generator
#  define 'attr_bool' to be an attribute ('gene')
#  which corresponds to integers sampled uniformly
#  from the range [0,1] (i.e. 0 or 1 with equal
#  robability)
#### register(alias, method[, argument[, ...]])
# random.randint(start, stop)

toolbox.register("attr_bool", rn.randint, 0, 100)

# Structure initializers
#     define 'individual' to be an individual
#     consisting of 100 'attr_bool' elements ('genes')
#### deap.tools.initRepeat(container, func, n)
#        container – The type to put in the data from func.
#        func – The function that will be called n times to fill the container.
#        n – The number of times to repeat func.
#        Returns:    An instance of the container filled with data from func.

toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_bool, 9)

# define the population to be a list of individuals
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

# The first registered function attr_bool redirects to the random.randint() function
# with its arguments fixed to sample numbers from the given range 0-100.
# The second registered function individual is a shortcut to the initRepeat() function,
# with its container argument set to the creator.Individual class and its generator argument
# to the toolbox.atts_bool() alias.
#[container – The type to put in the data from func.# generator – A function returning an iterable (list, tuple, …), the content of this iterable will fill the container.]
# Calling toolbox.individual() will call initRepeat() with the fixed arguments and return a complete creator.
# Individual composed of a permutation with a minimizing single objective fitness attribute.

count=0
g=0
strvar="test_{}_{}"

def fxy(x):

  L01=63.98
  L02=46.8
  L03=60.3
  L04=47.5
  L05=1.08
  R00=27.2
  h00=3.3
  h01=22.4
  zz0=6.8
  Ii0=350.0
 
 
  # L01=0
  # L02=0
  # L03=0
  # L04=0
  # L05=0
  # R00=0
  # h00=0
  # h01=0
  # zz0=0
  # Ii0=0  

  alpha1=0.5+2.0*x[0]/100.
  alpha2=0.5+2.0*x[1]/100.
  alpha3=0.5+2.0*x[2]/100.
  alpha4=0.5+2.0*x[3]/100.
  alpha5=0.5+2.0*x[4]/100.
  alpha6=0.5+2.0*x[5]/100.
  alpha7=0.5+2.0*x[6]/100.
  alpha8=0.5+2.0*x[7]/100.  # 1.0
  alpha9=0.5+2.0*x[8]/100.
 
  # alpha1=1
  # alpha2=1
  # alpha3=1
  # alpha4=1
  # alpha5=1
  # alpha6=1
  # alpha7=1
  # alpha8=1  # 1.0
  # alpha9=1

  L1=alpha1*L01
  L2=alpha2*L02
  L3=alpha3*L03
  L4=alpha4*L04
  L5=L05
  R0=alpha5*R00
  h0=h00
  h1=alpha6*h01
  zoff = alpha7*zz0
  I0 = alpha8*Ii0
 
  param =[L1, L2, L3, L4, L5, R0, h0, h1, zoff, I0]
  aux=[alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8, alpha9]

  # print(L01, L02, L03, L04, L05, R00, h00, h01, zz0, Ii0)
  # print(x)
  print("Horn paramrter List:", L1, L2, L3, L4, L5, R0, h0, h1, zoff, I0)
  # print("Horn paramrter List:", L01, L02, L03, L04, L05, R00, h00, h01, zz0, Ii0)

  return param

# the goal ('fitness') function to be maximized
def evalOneMax(individual):
  global count
  global g
  count=count+1  
  print("File:",count)
  print("Individual",individual)
  auxvar=fxy(individual)
  # auxvar1=np.zeros(len(auxvar))
  # auxvar1[0]=auxvar[0]
 
  # mydata = pd.read_csv("./test_{}_{}.txt".format(g,count), sep="  ", header=None,engine='python')
  # print(mydata)
 
 
  # with open('/home/vidhya/Documents/ds_und_code/20230406_EsbGA_fl/EsbGA_fl/test_{}_{}.txt', 'r', encoding='utf-8') as f:
  #     lines = f.readlines()

    # print(lines)

  mydata = pd.read_csv("./test _{}_{}.txt".format(g,count), delim_whitespace=True,header=None,engine='python')
  mydata.columns=["idp0","idpp0","xx0","yy0","zz0","px0","py0","pz0","ww0","ppar","apar"]
 
 

 
  npi=np.size(mydata.idp0.values)
  idp=mydata.idp0.values
  idpp=mydata.idpp0.values
  xx=mydata.xx0.values
  yy=mydata.yy0.values
  zz=mydata.zz0.values
  px=mydata.px0.values
  py=mydata.py0.values
  pz=mydata.pz0.values
  ww=mydata.ww0.values

  weight0=1.e2/(4.e4)*2.16e23;
  Ldet=1.e7;

  th0=0.05;
  dr=Ldet*np.sin(th0);
  area=2*np.arccos(-1)*Ldet**2*(1-np.cos(th0))*1.e-4;

  weight=weight0/area;
  wnu0=ww*weight

  ausx=xx+px/pz*(Ldet-zz);
  ausy=yy+py/pz*(Ldet-zz);
  # print("ausx",ausx)
  # print("ausy",ausy)

  ausr=np.sqrt(ausx**2+ausy**2);

  kk=np.sqrt(px**2+py**2+pz**2);

  kpp=[]
  kpp1=[]
  kpp2=[]
  kpp3=[]
  kpp4=[]
  npp1=0
  npp2=0
  npp3=0
  npp4=0
  knp=[]

 
  npp=0

  xp=[]
  xn=[]

  wpp=[]
  wnp=[]
  nnp=0
  xp=[]
  xo=[]
  yo=[]
  zo=[]
  x1=[]
  y1=[]
  x2=[]
  y2=[]
  x3=[]
  y3=[]
  x4=[]
  y4=[]
 
  x11=[]
  y11=[]
  x22=[]
  y22=[]
  x33=[]
  y33=[]
  x44=[]
  y44=[]  
 
  z1=[]
  z2=[]
  z3=[]
  z4=[]
  n1=0
  n11=0
  x2=[]
  y2=[]
  z2=[]
  # print("IdP",idpp)
  print(npi)
  # print("MaxZ:",max(zz),"MinZ:",min(zz),"LenZ",len(zz))

 #relevant for neutrino studies only
  kppc1=[]
  nppc1=0
  kppc=[]
  nppc=0



########################################## NOTE: This study(neutrino) uses daughter particles hence the first column(idp0)
##########################################       However the for pion flux we use second column(idpp0) since we are not interested in the daughter particles
  for i in range(0,npi):
      if(idpp[i]==13):              #pi+
          kpp.append(kk[i])
          npp+=1
          xo.append(xx[i])
          yo.append(yy[i])
          zo.append(zz[i])
         
     
         
          if(zz[i]>=255 and zz[i]<505):# & zz[i]<355):
          # if(zz[i]==1255):    
              kpp1.append(kk[i])
              npp1+=1              
              if((xx[i]>0 and xx[i]<=200)and(yy[i]>0 and yy[i]<=200)):
                  x1.append(xx[i])
                  y1.append(yy[i])
                  z1.append(zz[i])
                  if((xx[i]>75 and xx[i]<=125)and(yy[i]>75 and yy[i]<=125)):
                      x11.append(xx[i])
                      y11.append(yy[i])

          if(zz[i]>=505 and zz[i]<755):# & zz[i]<355):
          # if(zz[i]==1255):    
              kpp2.append(kk[i])
              npp2+=1              
              if((xx[i]>0 and xx[i]<=200)and(yy[i]>0 and yy[i]<=200)):
                  x2.append(xx[i])
                  y2.append(yy[i])
                  z2.append(zz[i])
                  if((xx[i]>75 and xx[i]<=125)and(yy[i]>75 and yy[i]<=125)):
                      x22.append(xx[i])
                      y22.append(yy[i])
                 
          if(zz[i]>=755 and zz[i]<1005):# & zz[i]<355):
          # if(zz[i]==1255):    
              kpp3.append(kk[i])
              npp3+=1              
              if((xx[i]>0 and xx[i]<=200)and(yy[i]>0 and yy[i]<=200)):
                  x3.append(xx[i])
                  y3.append(yy[i])
                  z3.append(zz[i])
                  if((xx[i]>75 and xx[i]<=125)and(yy[i]>75 and yy[i]<=125)):
                      x33.append(xx[i])
                      y33.append(yy[i])
               
           # if((xx[i]>0 and xx[i]<=200)):
               # xp.append(xx[i])
          if(zz[i]>=1005 and zz[i]<1255):# & zz[i]<355):
              kpp4.append(kk[i])
              npp4+=1
               # if((xx[i]>0 and xx[i]<=200)and(yy[i]>0 and yy[i]<=200)):
              x4.append(xx[i])
              y4.append(yy[i])
              z4.append(zz[i])
              if((xx[i]>75 and xx[i]<=125)and(yy[i]>75 and yy[i]<=125)):
                  x44.append(xx[i])
                  y44.append(yy[i])

             
  print("**************************************")
  print("Total:",npp)
  print("**************************************")
  print("Number of pi+ 1255>z>255:",npp1)
  print("RATIO:",npp1/npp)
  print("**************************************")
  print("Number of pi+ 2255>z>1255:",npp2)
  print("RATIO",npp2/npp)
 
  fig = plt.figure(figsize = (12,9))
  plt.hist(kpp,label="Total $\pi^+$,  $\Sigma \pi^+$:%d" %(npp),bins=100,color='navy')
  plt.hist(kpp1,label="focused $\pi^+$ 2.5m , $\Sigma \pi^+$:%d" %(npp1),bins=100,alpha=0.75,color="palevioletred")
  plt.hist(kpp2,label="focused $\pi^+$ 2.5-5.0m,$\Sigma \pi^+$:%d" %(npp2),bins=100,alpha=0.75,color="yellow")
  plt.hist(kpp3,label="focused $\pi^+$ 5.0-7.5m, $\Sigma \pi^+$:%d" %(npp3),bins=100,alpha=0.75,color="violet")
  plt.hist(kpp4,label="focused $\pi^+$ 7.5-10.0m, $\Sigma \pi^+$:%d" %(npp4),bins=100,alpha=0.75,color="green")

  #  # plt.hist(kppc,label="",bins=100)
  #  # plt.hist(kppc1,label="",bins=100,alpha=0.7)
  plt.legend(loc="upper right", fontsize=15)
  # plt.text(0.95, 10000, auxvar, fontsize = 5)
  # plt.text(1.5, 7500,npp1, fontsize = 11)
  # plt.text(0.75, 7500, auxvar[:5], fontsize = 12)
  plt.xlabel('Energy [GeV]',fontsize = 18)
  plt.ylabel('$\pi^{+}$flux',fontsize = 18)
  plt.xticks(fontsize=18)
  plt.yticks(fontsize=18)
  # plt.title(npp1)
  plt.show()
 
 
  fig = plt.figure(figsize = (12,9))
  def gaussian(x, mean, amplitude, standard_deviation):
      return amplitude * np.exp( - (x - mean)**2 / (2*standard_deviation ** 2))
  x = x11
 
  # plt.hist(x1,)
  bin_heights, bin_borders, _ = plt.hist(x1,histtype='bar',color='red' ,bins=200,label="$\mu$:%f, $\sigma$:%f " % (np.mean(x1), np.std(x1)),alpha=0.4)
  bin_heights, bin_borders, _ = plt.hist(x, bins=30,alpha=0.01)
  bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2
        # popt, _ = curve_fit(gaussian, bin_centers, bin_heights, p0=[1., 0., 1.])
  popt, _ = curve_fit(gaussian, bin_centers, bin_heights, p0=[75., 100., 125.])
       #popt returns mean, amplitude and standard deviation                                                                                                                                                        
       #_ returns the corresponding errors
  x_interval_for_fit = np.linspace(bin_borders[0], bin_borders[-1], 10000)
  # x_interval_for_fit1 = np.linspace(bin_borders1[0], bin_borders1[-1], 10000)
  plt.plot(x_interval_for_fit, gaussian(x_interval_for_fit, *popt),label="$\mu$:%f , $\sigma$:%f " % (popt[0], popt[2]),linewidth=3.5)
  plt.legend(loc='upper right',fontsize = 15)
  plt.xlabel('x [cm]',fontsize = 18)
  plt.ylabel('entries',fontsize = 18)
  plt.xticks(fontsize=18)
  plt.yticks(fontsize=18)
       # plt.text(0.65, 2.9, r'$\Sigma \pi^{+}:$',fontsize = 12)
       # plt.text(0.69, 2.9, npp,fontsize = 12 )
       # plt.title("$\Sigma \pi^{+}$:",npp)
       # plt.xlim([-2,8])
     # plt.title(kpp)  
  plt.show()
 
 
  fig = plt.figure(figsize = (12,9))
  def gaussian(x, mean, amplitude, standard_deviation):
      return amplitude * np.exp( - (x - mean)**2 / (2*standard_deviation ** 2))
  x = y11
      # plt.hist(x1)
   # bin_heights, bin_borders, _ = plt.hist(x, bins=40,density='True',label="$\mu$:%f, $\sigma$:%f " % (np.mean(y11), np.std(y11)),alpha=0.7)
  bin_heights, bin_borders, _ = plt.hist(y1, bins=200,histtype='bar',color='red',label="$\mu$:%f, $\sigma$:%f " % (np.mean(y1), np.std(y1)),alpha=0.4)
  bin_heights, bin_borders, _ = plt.hist(x, bins=30,alpha=0.01)
  bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2
  popt, _ = curve_fit(gaussian, bin_centers, bin_heights, p0=[75., 100., 125.])
  x_interval_for_fit = np.linspace(bin_borders[0], bin_borders[-1], 10000)
  plt.plot(x_interval_for_fit, gaussian(x_interval_for_fit, *popt),label="$\mu$:%f , $\sigma$:%f " % (popt[0], popt[2]),linewidth=3.5)
  plt.legend(loc='upper right',fontsize = 15)
  plt.xlabel('y [cm]',fontsize = 18)
  plt.ylabel('entries',fontsize = 18)
  plt.xticks(fontsize=18)
  plt.yticks(fontsize=18)
       # plt.text(0.65, 2.9, r'$\Sigma \pi^{+}:$',fontsize = 12)
       # plt.text(0.69, 2.9, npp,fontsize = 12 )
       # plt.title("$\Sigma \pi^{+}$:",npp)
       # plt.xlim([-2,8])
     # plt.title(kpp)  
  plt.show()
 
 
  # plt.figure(figsize=(18, 9))
  # hist3 = plt.hexbin(x1, y1, gridsize=30, cmap='Greens', mincnt=10,alpha=0.02)
  # hist2 = plt.hexbin(x2, y2, gridsize=30, cmap='Oranges', mincnt=10,alpha=0.02)
  # hist1 = plt.hexbin(x1, y1, gridsize=30, cmap='Blues', mincnt=10,alpha=0.02)
 
  # plt.colorbar(hist3, orientation='vertical')
  # plt.colorbar(hist2, orientation='vertical')
  # plt.colorbar(hist1, orientation='vertical')
  #  # plt.xlabel("x [cm]",fontsize = 12)
  #  # plt.ylabel("y [cm]",fontsize = 12)
  #  # plt.xticks(fontsize=12)
  #  # plt.yticks(fontsize=12)
  # plt.show() 
 
  # fig = plt.subplots(figsize =(12, 9))
  # plt.hist2d(z,x,bins=(180,180),cmap=('plasma'))
  # plt.colorbar()
  # plt.xlabel("z[cm]",fontsize = 12)
  # plt.ylabel("x [cm]",fontsize = 12)
  # plt.xticks(fontsize=12)
  # plt.yticks(fontsize=12)
  # plt.show()
  # plt.close()
 
  #spherical plots 
  # def ellipse(ra,rb,ang,x0,y0,Nb=100):
  #     xpos,ypos=x0,y0
  #     radm,radn=ra,rb
  #     an=ang
  #     co,si=np.cos(an),np.sin(an)
  #     the=linspace(0,2*np.pi,Nb)
 
  #     X=radm*np.cos(the)*co-si*radn*np.sin(the)+xpos
  #     Y=radm*np.cos(the)*si+co*radn*np.sin(the)+ypos
  #     return X,Y
   
  #    # Define the x and y data
  #    # For example just using random numbers
  # x = x1
  # y = y1
   
  #    # Set up default x and y limits
  # xlims = [min(x),max(x)]  
  # ylims = [min(y),max(y)]
   
  #    # Set up your x and y labels
  # xlabel = '$\mathrm{x [cm]}$'
  # ylabel = '$\mathrm{y [cm]}$'
   
  #    # Define the locations for the axes
  # left, width = 0.12, 0.55
  # bottom, height = 0.12, 0.55
  # bottom_h = left_h = left+width+0.02
   
  #    # Set up the geometry of the three plots
  # rect_temperature = [left, bottom, width, height] # dimensions of temp plot
  # rect_histx = [left, bottom_h, width, 0.25] # dimensions of x-histogram
  # rect_histy = [left_h, bottom, 0.25, height] # dimensions of y-histogram
   
  #    # Set up the size of the figure
  # fig = plt.figure(1, figsize=(9.5,9))
   
  #    # Make the three plots
  # axTemperature = plt.axes(rect_temperature) # temperature plot
  # axHistx = plt.axes(rect_histx) # x histogram
  # axHisty = plt.axes(rect_histy) # y histogram
   
  #    # Remove the inner axes numbers of the histograms
  # nullfmt = NullFormatter()
  # axHistx.xaxis.set_major_formatter(nullfmt)
  # axHisty.yaxis.set_major_formatter(nullfmt)
   
  #    # Find the min/max of the data
  # xmin = min(xlims)
  # xmax = max(xlims)
  # ymin = min(ylims)
  # ymax = max(y)
   
  #    # Make the 'main' temperature plot
  #    # Define the number of bins
  # nxbins = 50
  # nybins = 50
  # nbins = 100
   
  # xbins = linspace(start = xmin, stop = xmax, num = nxbins)
  # ybins = linspace(start = ymin, stop = ymax, num = nybins)
  # xcenter = (xbins[0:-1]+xbins[1:])/2.0
  # ycenter = (ybins[0:-1]+ybins[1:])/2.0
  # aspectratio = 1.0*(xmax - 0)/(1.0*ymax - 0)
   
  # H, xedges,yedges = np.histogram2d(y,x,bins=(ybins,xbins))
  # X = xcenter
  # Y = ycenter
  # Z = H
   
  #    # Plot the temperature data
  # cax = (axTemperature.imshow(H, extent=[xmin,xmax,ymin,ymax],
  #        interpolation='nearest', origin='lower',aspect=aspectratio))
   
  #    # Plot the temperature plot contours
  # contourcolor = 'white'
  # xcenter = np.mean(x)
  # ycenter = np.mean(y)
  # ra = np.std(x)
  # rb = np.std(y)
  # ang = 0
   

  #    #Plot the axes labels
  # axTemperature.set_xlabel(xlabel,fontsize=25)
  # axTemperature.set_ylabel(ylabel,fontsize=25)
   
  #    #Make the tickmarks pretty
  # ticklabels = axTemperature.get_xticklabels()
  # for label in ticklabels:
  #     label.set_fontsize(18)
  #     label.set_family('serif')
   
  # ticklabels = axTemperature.get_yticklabels()
  # for label in ticklabels:
  #     label.set_fontsize(18)
  #     label.set_family('serif')
   
  #    #Set up the plot limits
  # axTemperature.set_xlim(xlims)
  # axTemperature.set_ylim(ylims)
   
  #    #Set up the histogram bins
  # xbins = np.arange(xmin, xmax, (xmax-xmin)/nbins)
  # ybins = np.arange(ymin, ymax, (ymax-ymin)/nbins)
   
  #     #Plot the histograms
  # axHistx.hist(x, bins=xbins, color = 'blue',histtype='bar',ec="red")
  # axHisty.hist(y, bins=ybins, orientation='horizontal', color = 'red',histtype="bar",ec='blue')
 
 
  #    #Make the tickmarks pretty
  # ticklabels = axHistx.get_yticklabels()
  # for label in ticklabels:
  #     label.set_fontsize(12)
  #     label.set_family('serif')
   
  #    #Make the tickmarks pretty
  # ticklabels = axHisty.get_xticklabels()
  # for label in ticklabels:
  #     label.set_fontsize(12)
  #     label.set_family('serif')
   
  #    #Cool trick that changes the number of tickmarks for the histogram axes
  # axHisty.xaxis.set_major_locator(MaxNLocator(4))
  # axHistx.yaxis.set_major_locator(MaxNLocator(4))

  #    #Show the plot
  # plt.draw()
   

