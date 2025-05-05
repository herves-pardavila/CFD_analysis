# -*- coding: utf-8 -*-
"""
Created on Tue Mar 25 16:13:06 2025

@author: e6040
"""

import pandas as pd
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

plt.close("all")
class Temporal_Analysis:
    
    def __init__(self,path,ustar,wind_dir,field,part,coords):
        
        self.path=path
        self.ustar=ustar
        self.wind_dir=wind_dir
        self.field=field
        self.part=part
        self.coords=coords
        self.time=np.arange(3600, 3600*50,3600)
        self.hours=["7","8","9","10","11","12","13","14","15","16","17",
                    "18","19","20","21","22","23","24","1","2","3","4","5",
                    "6b","7b","8b","9b","10b","11b","12b","13b","14b","15b",
                    "16b","17b","18b","19b","20b","21b","22b","23b","24b","1b",
                    "2b","3b","4b","5b","6c"]
        
        return
    
    
    def surface_average(self,dataset, field,printing=False):
        """
        Calcula la media de un campo escalar en un plano, ponderando con el área
        de cada celda

        Parameters
        ----------
        dataset : xarray.Dataset
            Los datos de una tabla exportada de Star-CCM pasados a xarray dataset
        field : string
            Nombre del campo escalar

        Returns
        -------
        surf_avg : xarray.DataArray
            Media superficial del campo escalar

        """
        
        total_area = dataset["Area: Magnitude (m^2)"].sum()
        
        surface_sum = ( dataset[field] * dataset["Area: Magnitude (m^2)"] ).sum()
        
        surf_avg = surface_sum / total_area
        
        self.average=surf_avg
        
        
        return
    
    def surface_std1(self,dataset, field,printing=False):
        """
        Calcula la desviación estándar de un campo escalar en un plano, ponderando 
        con el área de cada celda

        Parameters
        ----------
        dataset : xarray.Dataset
            Los datos de una tabla exportada de Star-CCM pasados a xarray dataset
        field : string
            Nombre del campo escalar

        Returns
        -------
        surface_std : xarray.DataArray
           Desviación estándar superficial del campo escalar

        """
        
        total_area = dataset["Area: Magnitude (m^2)"].sum()
        
        surf_avg=self.average
       
        surface_sum = ( dataset["Area: Magnitude (m^2)"] * ( dataset[field] - surf_avg )**2 ).sum()
        
        surf_std = np.sqrt( surface_sum / total_area )
        
        self.standard_deviation=surf_std
        
        return
    
    def load_boundary_conditions(self,table_name="condiciones_contorno.csv",variable="T_AIRE"):
        
        
        df=pd.read_csv(self.path+"condiciones_contorno.csv",sep=";")
        
        self.boundary_condition=np.array(df[variable])
        return

    def load_data(self,variable,norm=1):
        
        means=np.zeros(len(self.time))
        means[:]=np.nan
        stds=np.copy(means)
        variables = [variable] +  ["Area: Magnitude (m^2)"]
        for hour in range(len(self.time)):
            
            #print(self.time[hour])
            
            try:
                df=pd.read_csv(self.path+self.ustar+self.wind_dir+str(self.time[hour]) +"_"+self.field + self.part + ".csv" )
                
                ds = xr.Dataset(
                    {var: ("punto", df[var].values) for var in variables},  # Variables escalares
                    coords={coord: ("punto", df[coord].values) for coord in self.coords},)
                
                
                self.surface_average(ds,variable)
                means[hour] = self.average*norm
             
                self.surface_std1(ds,variable)
                stds[hour]= self.standard_deviation*norm
            except FileNotFoundError:
                print("File Not Found for %i seconds" %self.time[hour])
                
                
        self.means=means
        self.stds=stds
        
    
    
    def plot_twin_axes_means(self,variables_left,variables_right,left_axis_name,right_axis_name,x_axis_name):
        
        fig=plt.figure()
        ax=fig.add_subplot(111)
        ax.set_ylabel(left_axis_name)
        ax.set_xlabel(x_axis_name)
        
        for var in variables_left:
            self.load_data(var)
            ax.plot(self.time,self.means,"-o",label=var,lw=3)
        
        if variables_right != None:
            
            axtwin=ax.twinx()
            axtwin.set_ylabel(right_axis_name)
            
            for var in variables_right:
                self.load_data(var)
                axtwin.plot(self.time,self.means,"-o",label=var,lw=3,ls="dashed")
       
        fig.legend()
                
            
            

        return   
        
                
        

if __name__== "__main__":
    
    #path_to_folder=r"C:/TFM_DAVIDH_2025/Tablas/"
    path_to_folder="../Tablas/"
    
    rho_air=1.18415 #kg/m3
    Qv=9.041e3 #kg/s
    A= 9.041e3 #m2
    
    eje1="Temperature (K)"
    eje2="Turbulent Kinetic Energy (J/kg)"
    
    
    print("u*=0.33 y NO")
    #simulacion NO u*= 0.33
    clase=Temporal_Analysis(path=path_to_folder, ustar="ustar0.33_", 
                            wind_dir="NO_", field="PS_T_V_", part="threshold_plane_z=3m",
                            coords=["X (m)", "Y (m)"])
    
    clase.load_data(eje1,norm=1)
    ps_NO_033=clase.means
    ps_NO_033_sd=clase.stds
    clase.load_data(eje2)
    tke_NO_033=clase.means
    
    
    
    print("u*=0.44 y NO")
    #simulacion NO u*= 0.44
    clase=Temporal_Analysis(path=path_to_folder, ustar="ustar0.44_", 
                            wind_dir="NO_", field="PS_T_V_", part="threshold_plane_z=3m",
                            coords=["X (m)", "Y (m)"])
    
    
    
    
    
    
    clase.load_data(eje1,norm=1)
    ps_NO_044=clase.means
    ps_NO_044_sd=clase.stds
    clase.load_data(eje2)
    tke_NO_044=clase.means
    
    
    print("u*=0.33 y NE")
    #simulacion NE u*= 0.33
    clase=Temporal_Analysis(path=path_to_folder, ustar="ustar0.33_", 
                            wind_dir="NE_", field="PS_T_V_", part="threshold_plane_z=3m",
                            coords=["X (m)", "Y (m)"])
    
    clase.load_data(eje1,norm=1)
    ps_NE_033=clase.means
    ps_NE_033_sd=clase.stds
    clase.load_data(eje2)
    tke_NE_033=clase.means
    
    
    
    print("u*=0.44 y NE")
    #simulacion NE u*= 0.44
    clase=Temporal_Analysis(path=path_to_folder, ustar="ustar0.44_", 
                            wind_dir="NE_", field="PS_T_V_", part="threshold_plane_z=3m",
                            coords=["X (m)", "Y (m)"])
    
    clase.load_data(eje1,norm=1)
    ps_NE_044=clase.means
    ps_NE_044_sd=clase.stds
    clase.load_data(eje2)
    tke_NE_044=clase.means
    
    
    clase.load_boundary_conditions(variable=["T_AIRE", "R.Dir (Wh/m2)"])
    tair=clase.boundary_condition[:,0]
    rdir=clase.boundary_condition[:,1]
    
    
    fig=plt.figure(figsize=(10,5))
    
    ax=fig.add_subplot(111)
    ax.set_ylabel(eje1)
    # ax.errorbar(clase.time,ps_NO_033,ps_NO_033_sd,fmt="-o",label="u*=0.33 NO",lw=3,capsize=10,capthick=2)
    # ax.errorbar(clase.time,ps_NO_044,ps_NO_044_sd,fmt="-o",label="u*=0.44 NO",lw=3,capsize=10,capthick=2)
    # ax.errorbar(clase.time,ps_NE_033,ps_NE_033_sd,fmt="-o",label="u*=0.33 NE",lw=3,capsize=10,capthick=2)
    # ax.errorbar(clase.time,ps_NE_044,ps_NE_044_sd,fmt="-o",label="u*=0.33 NE",lw=3,capsize=10,capthick=2)
    ax.plot(clase.time,ps_NO_033,"-o",label="u*=0.33 NO",lw=3)
    ax.plot(clase.time,ps_NO_044,"-o",label="u*=0.44 NO",lw=3)
    ax.plot(clase.time,ps_NE_033,"-o",label="u*=0.33 NE",lw=3)
    ax.plot(clase.time,ps_NE_044,"-o",label="u*=0.44 NE",lw=3)

    
    axtwin=ax.twinx()
    
    # axtwin.plot(clase.time,tke_NO_033,"-o", lw=2, ls="dashed")
    # axtwin.plot(clase.time,tke_NO_044,"-o", lw=2, ls="dashed")
    # axtwin.plot(clase.time,tke_NE_033,"-o", lw=2, ls="dashed")
    # axtwin.plot(clase.time,tke_NE_044,"-o", lw=2, ls="dashed")
    ax.plot(clase.time,tair + 273.15,"-o", lw=2, ls="dashed",color="black",label="Taire")
    axtwin.set_ylabel("T aire")
    fig.legend()
    ax.set_xlabel("Tiempo (s)")
    #ax.set_ylim(280,310)

    
    # fig2=plt.figure()
    # ax2=fig2.add_subplot(111)
    # ax2.set_xlabel(eje2)
    # ax2.set_ylabel(eje1)
    
    
    # ax2.plot(tke_NO_033,ps_NO_033,"-o",label="u*=0.33 NO")
    # ax2.plot(tke_NO_044,ps_NO_044,"-o",label="u*=0.44 NO")
    # ax2.plot(tke_NE_033,ps_NE_033,"-o",label="u*=0.33 NE")
    # ax2.plot(tke_NE_044,ps_NE_044,"-o",label="u*=0.44 NE")
    
    # fig2.legend()
    
    
    
    

    
    

    
    

    
    
    
    
    
    
    
    
    
    
    
    