# -*- coding: utf-8 -*-
"""
Created on Tue Apr 15 12:57:17 2025

@author: herve
"""

import numpy as np
from tqdm import tqdm
import pandas as pd
from scipy.stats.stats import pearsonr
import matplotlib.pyplot as plt

plt.close("all")

def pearsonr_2D(y, x):
    upper = np.sum((x - np.mean(x)) * (y - np.mean(y, axis=1)[:,None]), axis=1)
    lower = np.sqrt(np.sum(np.power(x - np.mean(x), 2)) * np.sum(np.power(y - np.mean(y, axis=1)[:,None], 2), axis=1))
    rho = upper / lower
    return rho
class Multivariant_Analysis:
    
    def __init__(self,path,ustar,part,time_ini,time_fin,wind_dir,field,coords=["X (m)", "Y (m)"],color_map ="seismic"):
        
        self.path=path
        self.ustar=ustar
        self.part=part
        self.time_ini=time_ini
        self.time_fin=time_fin
        self.wind_dir=wind_dir
        self.field=field
        
        self.color_map = color_map
        
        hours=np.arange(self.time_ini, self.time_fin + 3600,3600)
        self.hours=hours
        
        
        return
    def load_single_hour(self,hour,variable,norm=1):
        
        df=pd.read_csv(self.path+self.ustar+self.wind_dir+str(hour)+"_"+self.field+self.part+".csv")
        
        df.sort_values(by="X (m)",inplace=True,ignore_index=True)
        
        data=np.array(df[variable])
        
        self.data=data.T*norm
       
        X = np.array(df["X (m)"])
        self.X =X
        
        Y = np.array(df["Y (m)"])
        self.Y = Y
        
        area=np.array(df["Area: Magnitude (m^2)"])
        self.area=area
        
        return
    
    
    def load_data(self,variable,norm=1):
        

        matrix=[]
        
        for hour in tqdm(range(len(self.hours))):
            
            try:
               
                self.load_single_hour(self.hours[hour], variable,norm=norm)
                
                matrix+=[self.data]
                
                # if hour==0:
                    
                #     X = np.array(df["X (m)"])
                #     self.X =X
                    
                #     Y = np.array(df["Y (m)"])
                #     self.Y = Y
                    
                #     area=np.array(df["Area: Magnitude (m^2)"])
                #     self.area=area
                
            except FileNotFoundError:
                print("File Not Found for %i seconds" %self.hours[hour])
                
        self.matrix=np.array(matrix).T
        
        return 
    
    def regression_map(self,Index,variable,standarized=False,norm=1,alpha=0.1,sig="test-t",pp=100):
        
        #datos del campo (contaminante, velocidad, temperature, tke)
        self.load_data(variable,norm=norm)
        Data=self.matrix #dimensiones espacio x tiempo
        
        if standarized == True:
            Data = ( Data - Data.mean(axis=1,keepdims=True) ) / ( Data.std(axis=1, keepdims=True) )
        else:
            Data = Data - Data.mean(axis=1,keepdims=True)
        print(Data.shape) 
        
        try:
            [ns,nt]=Data.shape # ns=espacio, nt=tiempo #en este caso data es un campo
        except ValueError:
            # si Data es un índice
            ns=1
            nt=len(Data)
            Data = np.array([Data])
            
        cor=np.ma.empty([ns,]) #aquí voy a ir guardando los coeficientes de correlación para cada punto del espacio entre el campo y el índice
        Pvalue=np.ma.empty([ns,]) #lo mismo de arriba pero para los p-values
        reg=np.dot(Data,Index)/(nt-1) #matriz de covarianzas C=Y*X^T
        
        for nn in range(ns): 
            bb=pearsonr(Data[nn,:],Index) #bb es el coeficiente de correlacion de pearson entre el indice y el valor del campo en cada punto del espacio
            cor[nn]=bb[0] #guardo bb en un vector
            Pvalue[nn]=bb[1] #guardo el p-valor que me dice si bb es significativo en otro vector
        
        if sig == 'test-t':
            cor_sig=np.ma.masked_where(Pvalue>alpha,cor)
            reg_sig=np.ma.masked_where(Pvalue>alpha,reg)
            
        if sig == 'MonteCarlo':
            corp = np.ma.empty([ns,pp])
            for p in range(pp):
                corp[:,p] = pearsonr_2D(Data,np.random.permutation(Index))
                # aquí uso la función pearsonr_2D y me ahorro un bucle en ns
            
            for nn in range(ns): 
                hcor = np.count_nonzero((cor[nn]>0)&(corp[nn,:]<cor[nn])|(cor[nn]<0)&(corp[nn,:]>cor[nn]))
                # nivel de confianza
                Pvalue[nn] = hcor/pp
                
            cor_sig = np.ma.masked_where(Pvalue<(1-alpha),cor)
            reg_sig = np.ma.masked_where(Pvalue<(1-alpha),reg)
        
        
        self.cor_sig=cor_sig
        self.reg_sig=reg_sig
        self.reg=reg
        self.cor=cor
        self.Pvalue=Pvalue
            
            
        return
    
    def plot_regression_map(self,Index,variable,norm=1,standarized=False,alpha=0.1,sig="test-t",pp=100,area_factor=1.5):
        
        self.regression_map(Index, variable,norm=norm,standarized=standarized,alpha=alpha,sig="test-t",pp=100)
        
        self.plot(self.reg_sig, str(self.ustar)+str(self.wind_dir) +" Significative Regression Map for \n" + variable)
        self.plot(self.cor_sig, str(self.ustar)+str(self.wind_dir) +" Significative Correlation Map for \n" + variable)
        self.plot(self.reg, str(self.ustar)+str(self.wind_dir) +" Regression Map for \n" + variable)
        self.plot(self.cor, str(self.ustar)+str(self.wind_dir) +" Correlation Map for \n" + variable)

             
        
        
        return
    
    def principal_components(self,Index,variable,norm=1):
        
        #datos del campo (contaminante, velocidad, temperature, tke)
        self.load_data(variable,norm=norm)
        Data=self.matrix #dimensiones espacio x tiempo
        Data = Data - Data.mean(axis=1,keepdims=True) #anomalo
        
        #calculo la matriz de covarianza del campo anómalo
        C=np.dot(Data,Data.T)
        eof,d,eof2=np.linalg.svd(C) #eof son los autovalores, d son los autovectores
        
        lambdas_normalizados=(d)/sum(d)
        
        #calculo las componentes principales
        PC=np.dot(Data.T,eof[:,:3]) # el campo es de dimensiones nsxnt y eof es de dimensiones nsx3y 
        PCs=(PC-np.mean(PC))/np.std(PC) #estandarizo las componentes principales
        
        self.lambdas=lambdas_normalizados
        self.eof=eof
        self.eof2=eof
        self.PC = PCs
        
        fig=plt.figure()
        ax1=fig.add_subplot(121)
        ax1.set_ylabel("Explained Variance")
        ax1.set_xlabel("Eigenvalue")
        ax2=fig.add_subplot(122)
        ax2.set_ylabel("Principal Component")
        ax2.set_xlabel("Time")
       
        ax1.plot(lambdas_normalizados[:5],"o",color='Blue')
        ax2.plot(PCs[:,0],"-o",label= "First PC")
        ax2.plot(PCs[:,1],"-o",label= "Second PC")
        ax2.plot(PCs[:,2],"-o",label= "Third PC")
        ax2twin=ax2.twinx()
        ax2twin.plot(Index,label="Index")
        
        fig.legend()
        
        return
    
    def plot_principal_component(self,Index,variable,number,norm=1,standarized=False,alpha=0.1,sig="test-t",pp=100):
        
        self.principal_components(Index, variable,norm=norm)
        PCs=self.PC[:,number]
        self.regression_map(np.transpose(PCs), variable,standarized=standarized,norm=norm,
                            alpha=alpha,sig=sig,pp=pp)
        
        
        self.plot(self.reg_sig, str(self.ustar)+str(self.wind_dir) +"PC number" +str(int(number+1)) +" Significative Regression Map for \n" + variable)
        self.plot(self.cor_sig, str(self.ustar)+str(self.wind_dir) +"PC number" +str(int(number +1))+" Significative Correlation Map for \n" + variable)
        self.plot(self.reg, str(self.ustar)+str(self.wind_dir) +"PC number" +str(int(number)+1)+" Regression Map for \n" + variable)
        self.plot(self.cor, str(self.ustar)+str(self.wind_dir) +"PC number" +str(int(number)+1)+" Correlation Map for \n" + variable)
        
        return
        
        
        
        
        
        
        
        
    
    def plot(self,data,title,area_factor=1.5,bounds=None):
        
        if bounds != None:
            vmin,vmax=bounds[0],bounds[1]
        
        else:
            vmin,vmax=data.min(),data.max()
        
        fig=plt.figure(figsize=(area_factor*6.4,area_factor*4.8))
        fig.suptitle(title,fontsize=10*area_factor)
        ax=fig.add_subplot(111)
        ax.set_aspect('equal', adjustable='box')
        im=ax.scatter(self.X,self.Y,s=area_factor*self.area**2,c=data,marker="h",vmin=vmin,vmax=vmax,cmap=self.color_map,edgecolor=None)
        ax.set_xlabel("X (m) ",fontsize=area_factor*10)
        ax.set_ylabel("Y (m) ",fontsize=area_factor*10)
        ax.set_facecolor('gray')
        
        cbar=fig.colorbar(im,location="right")
        ax.tick_params(axis='both', which='both', labelsize=10*area_factor)
        cbar.ax.tick_params(labelsize=10*area_factor) 
        fig.text(0.4,0.2,"South",transform=ax.transAxes,color="white",fontsize=10*area_factor)
        fig.text(0.4,0.75,"North",transform=ax.transAxes,color="white",fontsize=10*area_factor)
        fig.text(0.65,0.5,"East",transform=ax.transAxes,color="white",fontsize=10*area_factor)
        fig.text(0.2,0.5,"West",transform=ax.transAxes,color="white",fontsize=10*area_factor)
        
        #plt.show()
        
        return
        
        
        

        


if __name__ == "__main__":
    
    #path_to_folder=r"C:/TFM_DAVIDH_2025/Tablas/" #PC CIEMAT
    
    path_to_folder=r"C:/Users/herve/Documents/TFM_DAVIDH_2025/Tablas/" #Asus
   
    campo="PS_T_V_" #passive scalar
    part="threshold_plane_z=3m" 
    con=Multivariant_Analysis(path=path_to_folder, ustar="ustar0.33_", part=part,
                         time_ini=3600, time_fin= 176400, wind_dir="NO_",
                         field=campo,color_map="seismic")
    
  
    
    #Indice
    bc=pd.read_csv(path_to_folder+"condiciones_contorno.csv",sep=";")
    indice = 'R.Dir (Wh/m2)'
    indice_anomalo = bc[indice] - bc[indice].mean()
    indice_standarizado = indice_anomalo / bc[indice].std()  
    
    #Regression and correlations maps
    con.plot_regression_map(indice_anomalo, "PS2", standarized=False, alpha=0.1)  
    
    #Principal Component Analysis
    
    #con.principal_components(indice_anomalo,"PS2")
   # con.plot_principal_component(indice_anomalo, "Vorticity: Magnitude (/s)", 0)
    