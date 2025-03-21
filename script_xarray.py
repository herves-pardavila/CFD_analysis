# -*- coding: utf-8 -*-
"""
Editor de Spyder

Script para representar el ciclo diurno de algunas variables escalares.

Autor David Hervés Pardavila 6 de Marzo del 2025
"""


import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt




def surface_average(dataset, field,printing=False):
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
    
    if printing == True:
        print("Surface average para %s = " %field , surf_avg)
    
    return surf_avg

def surface_std1(dataset, field,printing=False):
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
    
    surf_avg=surface_average(dataset, field)
    
    surface_sum = ( dataset["Area: Magnitude (m^2)"] * ( dataset[field] - surf_avg )**2 ).sum()
    
    surf_std = np.sqrt( surface_sum / total_area )
    
    if printing == True:
        print("Surface standard devitation para %s = " %field, surf_std)
    
    return surf_std

def surface_std2(dataset, field, printing=False):
    """
    Calcula la desviación estándar de un campo escalar en un plano, ponderando 
    con el área de cada celda. Esta función lo calcula teniendo en cuenta la 
    siguiente expresión para la varianza
    
    Var[x] = E[x^2] - E[x]^2
        

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
    #primero calculo E[x^2]
    
    total_area = dataset["Area: Magnitude (m^2)"].sum()
    
    surface_sum2 = ( ( dataset[field] ** 2 ) * dataset["Area: Magnitude (m^2)"] ).sum()
    
    surf_avg2 = surface_sum2 / total_area
    #print("djkfajkdas", surf_avg2)
    
    #segundo calculo la media E[x]
    
    surf_avg=surface_average(dataset, field)
    
    
    #ahora calculo s como E[x^2] - E[x]^2
    
    surf_variance = surf_avg2 - surf_avg**2
    
    surf_std = np.sqrt(surf_variance)
   
    if printing == True: 
       print("Surface standard devitation para %s = " %field, surf_std)
    
    return surf_std


if __name__ == "__main__":
    
    path_to_file=r"C:\TFM_DAVIDH_2025/Tablas/"

    ustar="ustar0.33_"
    viento="NO_"
    tiempo="7200_" # 3600 son las 7 am
    campo="PS_T_V_" #passive scalar
    part="threshold_plane_z=3m"
    
    # Cargar CSV en un DataFrame
    df = pd.read_csv(path_to_file+ustar+viento+tiempo+campo+part+".csv")
    df=df.drop(columns="Z (m)")
    # Suponiendo que las columnas están en orden: temperatura, presión, humedad, X, Y, Z
    variables = ["PS1", "PS2", "PS3", "Area: Magnitude (m^2)"]
    coords = ["X (m)", "Y (m)"]

    # Crear Dataset de xarray
    ds = xr.Dataset(
        {var: ("punto", df[var].values) for var in variables},  # Variables escalares
        coords={coord: ("punto", df[coord].values) for coord in coords},  # Coordenadas
        )

    print(ds)
    
    X = ds["X (m)"].values
    Y = ds["Y (m)"].values
    area=ds["Area: Magnitude (m^2)"]
    ps1 = ds["PS1"].values  # Puedes cambiar por "presion" o "humedad"
    
    
    #ds_grid = ds.set_index(punto=["X (m)", "Y (m)"]).unstack()
    #ds_grid.PS1.plot(cmap="coolwarm")
    #plt.title("Mapa de Temperatura")
    #plt.show()
        

    # Crear el gráfico de dispersión
    fig=plt.figure()
    ax=fig.add_subplot(111)
    im=ax.scatter(X,Y,s=area,c=ps1,marker="h",cmap="coolwarm",edgecolor=None)
    ax.set_aspect('equal', adjustable='box')
    fig.colorbar(im,location="top")
    fig.savefig("imagen_scatter_xarray.png")
    plt.show()

    
    #calculo el promedio ponderando con el area de cada celda
    s_avg = surface_average(ds, "PS1")
    s_std = surface_std1(ds, "PS1")
    s_std = surface_std2(ds, "PS1")



    
   