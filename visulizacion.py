import os
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
from ipywidgets import widgets, Output
from IPython.display import display as dp
import matplotlib.patches as patches
import pygmt
import xarray as xr
from pygmt import datasets
from pyproj import CRS, Transformer
import imageio.v2 as imageio
import math

output = Output() # usados para los Loading...
output0 = Output()  #usado para las representaciones en PyGMT
output1 = Output()  #usado para seleccionar opciones en las carpetas
output2 = Output()  #usado para los plots y algunos widgets
output3 = Output()  #usado para los sliders de las variables temporales
output4 = Output() #usado para botones de variables temporales

def display(path="./simulaciones"):
    """
    Escanea el directorio dado y genera botones para los archivos .nc
    y las carpetas, asociándolos a las funciones elegirdeuna y elegirarch.
    """
    print("Please, select a simulation output:")
    
    botones = []
    
    for nombre in sorted(os.listdir(path)):
        ruta_completa = os.path.join(path, nombre)

        if os.path.isfile(ruta_completa) and nombre.endswith(".nc"):
            boton = widgets.Button(description=f"{nombre[:-3]} Output", layout=widgets.Layout(width='300px', height='30px'))
            boton.on_click(lambda b, ruta=nombre: elegirdeuna(ruta))
            botones.append(boton)

        elif os.path.isdir(ruta_completa):
            boton = widgets.Button(description=f"[F] {nombre} Output", layout=widgets.Layout(width='300px', height='30px'))
            boton.on_click(lambda b, ruta=nombre: elegirarch(ruta))
            botones.append(boton)

    dp(widgets.HBox(botones))

    dp(output)
    dp(output0)
    dp(output1)
    dp(output2)
    dp(output3)
    dp(output4)
    
'''
def display():
    """
    Crea botones para seleccionar una carpeta y muestra el resultado.
    """
    print("Please, select a simulation output:")

    boton1 = widgets.Button(description="1 - Japón Output", layout=widgets.Layout(width='300px', height='30px'))
    boton2 = widgets.Button(description="2 - Cádiz Output", layout=widgets.Layout(width='300px', height='30px'))
    boton3 = widgets.Button(description="3 - Rules Output", layout=widgets.Layout(width='300px', height='30px'))
    boton4 = widgets.Button(description="4 - Trial Output", layout=widgets.Layout(width='300px', height='30px'))

    boton1.on_click(lambda b: elegirdeuna("japon7.nc"))  
    boton2.on_click(lambda b: elegirarch("Cadiz"))
    boton3.on_click(lambda b: elegirdeuna("Rules/rules.nc"))
    boton4.on_click(lambda b: elegirarch("prueba2")) 

    dp(widgets.HBox([boton1, boton2, boton3, boton4]))

    dp(output)
    dp(output0)
    dp(output1)
    dp(output2)
    dp(output3)
    dp(output4)
    
'''


def elegirdeuna(nombre):
    output0.clear_output()
    output1.clear_output()
    output2.clear_output() 
    output3.clear_output() 
    output4.clear_output()
    
    archivo = "./simulaciones/"+nombre 
    op = ""
    with output:
        print("Loading...")
    with output0:
            
        try:
            dataset = nc.Dataset(archivo, mode="r")
        except FileNotFoundError:
            print("Error: El archivo no existe.")
            return

        if 'lat' in dataset.variables and 'lon' in dataset.variables or 'y' in dataset.variables and 'x' in dataset.variables:
            keyy = 'lat' if 'lat' in dataset.variables else 'y'
            keyx = 'lon' if 'lon' in dataset.variables else 'x'
            x = dataset.variables[keyx][:]
            y = dataset.variables[keyy][:]

            fig = pygmt.Figure()

            if x.max() > 400 or y.max() > 400:
                fig.basemap(
                    region=[x.min(), x.max(), y.min(), y.max()], 
                    projection="x1:100000",
                    frame=True
                )
            else:
                fig.basemap(
                    region=[x.min(), x.max(), y.min(), y.max()],
                    projection="M6i",
                    frame=True
                )
    
            
            if 'original_bathy' in dataset.variables or 'bathymetry' in dataset.variables:
                key = 'original_bathy' if 'original_bathy' in dataset.variables else 'bathymetry'
                z = -dataset.variables[key][:]

            
                y
                grid = xr.DataArray(
                    z,
                    dims=["lat", "lon"],
                    coords={"lat": y, "lon": x}
                )
        
                fig.grdimage(grid=grid, shading="+a270+nt1", cmap="geo")
                fig.colorbar(frame=["a", "x+l'Elevación (m)'"])
            else : 
                fig.coast(
                    land="gray",
                    water="skyblue",
                    shorelines=True
                )
        else : 
            print("No hay variables lat y lon")
            
        fig.show()
        output.clear_output()

    elegirop(nombre)
        

def elegirarch(nombre):
    with output:
        print("Loading...")
    with output2:
        output2.clear_output()  
        carpeta = os.path.join("./simulaciones", nombre)

        archivos = os.listdir(carpeta)
        dataset_info = []
        
        for archivo in archivos:
            ruta_archivo = os.path.join(carpeta, archivo)
            nomvar = os.path.join(nombre,archivo)

            try:
                dataset = nc.Dataset(ruta_archivo, mode="r")
                #print(len(dataset.variables['lon'][:])*len(dataset.variables['lat'][:])*len(dataset.variables['time'][:]))
                if 'lat' in dataset.variables and 'lon' in dataset.variables or 'y' in dataset.variables and 'x' in dataset.variables:
                    keyy = 'lat' if 'lat' in dataset.variables else 'y'
                    keyx = 'lon' if 'lon' in dataset.variables else 'x'
                    x = dataset.variables[keyx][:]
                    y = dataset.variables[keyy][:]
        
                    if x.max()>400 or y.max() > 400:
                        coord_type = "UTM"
                        print("Coordenadas detectadas en UTM. Convirtiendo a lat/lon...")
            
                        # Ajustar zona UTM y hemisferio según corresponda
                        utm_zone = 30  # <-- ajusta esto según tu zona
                        hemisferio_sur = False  # cambia a True si estás en el hemisferio sur
            
                        crs_utm = CRS.from_proj4(f"+proj=utm +zone={utm_zone} +datum=WGS84 +units=m +{'south' if hemisferio_sur else 'north'}")
                        crs_wgs84 = CRS.from_epsg(4326)
                        transformer = Transformer.from_crs(crs_utm, crs_wgs84)
            
                        x,y = transformer.transform(x, y)
                        print("Conversión completa.")  
                    
                    if len(x) == 1 or len(y) == 1:
                        x_min, x_max = np.nan, x
                        y_min, y_max = np.nan, y
                        area = 1  #poner por debajo de 1/1000000 si quieres que en la numeracion aparezcan antes los recuadros
                    else:         
                        x_min, x_max = x.min(), x.max()
                        y_min, y_max = y.min(), y.max()
                        area = (x_max - x_min) * (y_max - y_min)  
                    
                    dataset_info.append((nomvar, x_min, x_max, y_min, y_max, area)) 
                """  
                quitar comillas para ver cuantos datos tiene cada archivo  
                
                if "lon" in dataset.variables and "lat" in dataset.variables and "time" in dataset.variables:    
                    print(len(dataset.variables['lon'][:])*len(dataset.variables['lat'][:])*len(dataset.variables['time'][:]))
                """   
                dataset.close()  
            except Exception:
                pass
                
        
        
        dataset_info.sort(key=lambda item: item[5], reverse=True)
        multiarch(dataset_info)
        
# esta funcion debe cambiar por la de gmt 
def multiarch(dataset_info):
    output0.clear_output() 
    output1.clear_output() 
    output2.clear_output() 
    output3.clear_output() 
    output4.clear_output()
    
    if dataset_info:
        ruta_archivo_top = "./simulaciones/"+dataset_info[0][0]  
        dataset = nc.Dataset(ruta_archivo_top, mode="r")  
        
        #x = dataset.variables['lon'][:]
        #y = dataset.variables['lat'][:]
        
        if 'lat' in dataset.variables and 'lon' in dataset.variables or 'y' in dataset.variables and 'x' in dataset.variables:
            keyy = 'lat' if 'lat' in dataset.variables else 'y'
            keyx = 'lon' if 'lon' in dataset.variables else 'x'
            x = dataset.variables[keyx][:]
            y = dataset.variables[keyy][:]

        if x.max() > 400 or y.max() > 400:
            # Si son arrays 1D, los convertimos a malla si es necesario
            if x.ndim == 1 and y.ndim == 1:
                x, y = np.meshgrid(x, y)
    
            # Flatten para tener vectores del mismo tamaño
            x_flat = x.flatten()
            y_flat = y.flatten()
    
            
            coord_type = "UTM"
            print("Coordenadas detectadas en UTM. Convirtiendo a lat/lon...")

            utm_zone = 30  # ⚠️ Ajusta según tu zona
            hemisferio_sur = False

            crs_utm = CRS.from_proj4(
                f"+proj=utm +zone={utm_zone} +datum=WGS84 +units=m +{'south' if hemisferio_sur else 'north'}"
            )
            crs_wgs84 = CRS.from_epsg(4326)
            transformer = Transformer.from_crs(crs_utm, crs_wgs84)

            x, y = transformer.transform(x_flat, y_flat)
            print("Conversión completa.") 
        
        a = 10
        l = int(a*len(x)/len(y))
        if len(x) != 1 and len(y) != 1:
            with output0:
                output0.clear_output()
                
                fig = pygmt.Figure()
                
                fig.basemap(
                    region=[x.min(), x.max(), y.min(), y.max()],
                    #projection=f"M{str(a)}c/{str(l)}c",
                    projection="M6i",
                    frame=True
                )
                if 'original_bathy' in dataset.variables:
                    z = -dataset.variables['original_bathy'][:]
            
                    # Crear grid de tipo xarray
                    grid = xr.DataArray(
                        z,
                        dims=["lat", "lon"],
                        coords={"lat": y, "lon": x}
                    )
            
                    fig.grdimage(grid=grid, shading="+a270+nt1", cmap="geo")
                    fig.colorbar(frame=["a", "x+l'Elevación (m)'"])
                else : 
                    fig.coast(
                        land="gray",
                        water="skyblue",
                        shorelines=True
                    )
                # Dibujar los rectángulos
                i = 0
                for datos in dataset_info[1:]:
                    x_min2, x_max2, y_min2, y_max2 = datos[1:5]
                    
                    # Dibujar el contorno del rectángulo
                    fig.plot(
                        x=[x_min2, x_max2, x_max2, x_min2, x_min2],
                        y=[y_min2, y_min2, y_max2, y_max2, y_min2],
                        pen="1p,black,--"
                    )
                    
                    # Añadir número encima del rectángulo
                    fig.text(
                        text=str(i+1),
                        x=x_min2,
                        y=y_max2,
                        font="10p,Helvetica-Bold,white",
                        justify="LT",
                        fill="black"
                    )
            
                    i += 1
                    
                fig.show()
                output.clear_output()
                
            with output1:
                output1.clear_output()  
                
                botones = []
                for i in range(0,len(dataset_info)):
                    boton = widgets.Button(description=str(i))
                    boton.on_click(lambda b,op = i: elegirdemult(dataset_info[op][0]))  
                    botones.append(boton)
                dp(widgets.HBox(botones))
        else:
            sectemporal_boya(dataset.variables['lo que sea'][:], dataset.variables['time'][:], 'boya')
       
def elegirdemult(nombre):
    output0.clear_output()
    output2.clear_output() 
    output3.clear_output() 
    output4.clear_output()
    # en esta funcion debe ir lo de gmt
    with output:
        print("Loading...")
    archivo = "/home/vedia/hysea-project/venv/simulaciones/"+nombre 
    op = ""
    #output2.clear_output()
    
    with output0:        
        
        try:
            dataset = nc.Dataset(archivo, mode="r")
        except FileNotFoundError:
            print("Error: El archivo no existe.")
            return

        if 'lat' in dataset.variables and 'lon' in dataset.variables or 'y' in dataset.variables and 'x' in dataset.variables:
            keyy = 'lat' if 'lat' in dataset.variables else 'y'
            keyx = 'lon' if 'lon' in dataset.variables else 'x'
            x = dataset.variables[keyx][:]
            y = dataset.variables[keyy][:]

            if x.max()>400 or y.max() > 400:
                coord_type = "UTM"
                print("Coordenadas detectadas en UTM. Convirtiendo a lat/lon...")
    
                # Ajustar zona UTM y hemisferio según corresponda
                utm_zone = 30  # <-- ajusta esto según tu zona
                hemisferio_sur = False  # cambia a True si estás en el hemisferio sur
    
                crs_utm = CRS.from_proj4(f"+proj=utm +zone={utm_zone} +datum=WGS84 +units=m +{'south' if hemisferio_sur else 'north'}")
                crs_wgs84 = CRS.from_epsg(4326)
                transformer = Transformer.from_crs(crs_utm, crs_wgs84)
    
                x,y = transformer.transform(x, y)
                print("Conversión completa.")  
        
        fig = pygmt.Figure()
        fig.basemap(
            region=[x.min(), x.max(), y.min(), y.max()],
            projection="M6i",
            frame=True
        )
        if 'original_bathy' in dataset.variables:
            z = -dataset.variables['original_bathy'][:]
    
            # Crear grid de tipo xarray
            grid = xr.DataArray(
                z,
                dims=["lat", "lon"],
                coords={"lat": y, "lon": x}
            )
    
            fig.grdimage(grid=grid, shading="+a270+nt1", cmap="geo")
            fig.colorbar(frame=["a", "x+l'Elevación (m)'"])
        else : 
            fig.coast(
                land="gray",
                water="skyblue",
                shorelines=True
            )
        fig.show()
        output.clear_output()

    elegirop(nombre)

    
def elegirop(nombre):
    # En esta función elegimos una de las variables del archivo NetCDF
    with output:
        print("Loading...")

    archivo = "/home/vedia/hysea-project/venv/simulaciones/" + nombre

    with output2:
        output2.clear_output()
        output3.clear_output()
        output4.clear_output()

        try:
            dataset = nc.Dataset(archivo, mode="r")
        except FileNotFoundError:
            print("Error: El archivo no existe.")
            return

        print("Select output variable:")

        botones = []

        for key in dataset.variables:
            if len(dataset.variables[key].shape) > 1:
                # Caso especial: arrival_times
                if key == "arrival_times":
                    boton = widgets.Button(description="arrival_times")
                    boton.on_click(lambda b, ds=dataset: procesar_arrival_times(ds))
                else:
                    boton = widgets.Button(description=key)
                    boton.on_click(lambda b, key=key, ds=dataset: elegirfun(key, ds))
                botones.append(boton)

        # Si están eta y original_bathy, añade botón flood
        if 'eta' in dataset.variables and 'original_bathy' in dataset.variables:
            boton = widgets.Button(description='flood')
            boton.on_click(lambda b, op1='eta', op2='original_bathy': inundacion(dataset, op1, op2))
            botones.append(boton)

        # Botón para operar variables
        boton = widgets.Button(description="Operate variables")
        boton.on_click(lambda b: operarvar(nombre))
        botones.append(boton)

        output.clear_output()
        dp(widgets.VBox(botones))



def elegirfun(op,dataset):
    """Función para manejar la selección de variables dentro del NetCDF."""
    #en funcion del tamaño de la variable elegida 3d o 2d ejecutaremos un camino u otro
    with output:
        print("Loading...")
    output0.clear_output()
    if len(dataset.variables[op].shape) == 3:
        variable = dataset.variables[op]
        
        if 'lat' in dataset.variables and 'lon' in dataset.variables and "time" in dataset.variables or 'y' in dataset.variables and 'x' in dataset.variables and "time" in dataset.variables:
            keyy = 'lat' if 'lat' in dataset.variables else 'y'
            keyx = 'lon' if 'lon' in dataset.variables else 'x'
            x = dataset.variables[keyx][:]
            y = dataset.variables[keyy][:]
            t = dataset.variables["time"][:]
            coord = True
        else:
            x = np.arange(variable.shape[2])  
            y = np.arange(variable.shape[1])
            t = np.arange(variable.shape[0])
            coord = False
        menu3d(variable ,op ,x ,y ,t ,coord)
        
        
    elif len(dataset.variables[op].shape) == 2:
        variable = dataset.variables[op][:]
        if 'lat' in dataset.variables and 'lon' in dataset.variables or 'y' in dataset.variables and 'x' in dataset.variables:
            keyy = 'lat' if 'lat' in dataset.variables else 'y'
            keyx = 'lon' if 'lon' in dataset.variables else 'x'
            x = dataset.variables[keyx][:]
            y = dataset.variables[keyy][:]
            coord = True
        else:
            x = np.arange(variable.shape[1])  
            y = np.arange(variable.shape[0])
            coord = False
        menu2d (variable, op, x, y,coord)


def menu2d (variable, op, x, y,coord):
    output0.clear_output()
    output2.clear_output()
    output3.clear_output()
    output4.clear_output()
    
    mostmap2d(variable, op, x, y,coord)
    with output4:
        
        
        
        #boton1 = widgets.Button(description= "Make a zoom",layout=widgets.Layout(width='300px', height='30px'))
        boton2 = widgets.Button(description= "Make transect (recomendable hacer un zoom primero)",
                                layout=widgets.Layout(width='300px', height='30px'))       

        #boton1.on_click(lambda b: zoom2d(variable, x, y, op, coord))
        boton2.on_click(lambda b: seccion2d (variable, x, y, op,coord))
            
        dp(widgets.VBox([boton2]))


def menu3d(variable, op, x, y,t ,coord):
    output0.clear_output()
    output2.clear_output()
    output3.clear_output()
    output4.clear_output()
    mostmap3d(variable, op, x, y, t, coord)
    with output4:
        
        
        
        #boton1 = widgets.Button(description= "Make a zoom",layout=widgets.Layout(width='300px', height='30px'))
        boton2 = widgets.Button(description= "Extract time series",layout=widgets.Layout(width='300px', height='30px'))  
        boton3 = widgets.Button(description= "Make a transect ",layout=widgets.Layout(width='300px', height='30px'))

        #boton1.on_click(lambda b: zoom3d(variable, op, x, y, t, coord))
        boton2.on_click(lambda b: sectemporal(variable, op, x, y, t, coord))
        boton3.on_click(lambda b: seccion3d(variable, op, x, y, t, coord))
            
        dp(widgets.VBox([boton2,boton3]))
 
def menu_inun(eta_ajust ,bathy ,op ,x ,y ,t ,coord):
    output0.clear_output()
    output2.clear_output()
    output3.clear_output()
    output4.clear_output()
    mostrar_inun(eta_ajust ,bathy ,op ,x ,y ,t ,coord)
    '''
    with output4:
        
        boton1 = widgets.Button(description= "Make a zoom",layout=widgets.Layout(width='300px', height='30px'))
        
        boton1.on_click(lambda b: zoom_inun(eta_ajust ,bathy ,op ,x ,y ,t ,coord))
            
        dp(widgets.VBox([boton1]))
     '''


def mostmap2d(variable, op, x, y,coord):
    #en este funcion mostramos las variables 2d 
        # Graficar la variable
    with output2:
        
        if (len(x)*len(y)>4000000):
            #print(len(x),len(y))
            x ,y , variable = comprimir(x, y, variable, 4)
            #print(len(x),len(y))
        
        a = 10
        l = int(a*len(x)/len(y))
        #print( x.min(), x.max(), y.min(), y.max())  
        if op == 'original_bathy' or op == 'deformed_bathy' or op == 'bathymetry':
            vmin = -np.nanmin(variable)
            vmax = -np.nanmax(variable)
            sig=-1
        else:
            vmin = np.nanmin(variable)
            vmax = np.nanmax(variable)
            sig=1

        escala = np.linspace(int(vmin), int(vmax),50)
        escala_ent = [round(num) for num in escala]
        
        slider1 = widgets.SelectionRangeSlider(
            options=escala_ent,  
            index=(0, len(escala) - 1),  
            description='Scale:',
            disabled=False,
            layout=widgets.Layout(width="700px")
        )
        boton = widgets.Button(description="Save figure")
        fig_container = {"fig": None}

        
        def actualizar(val):
            with output3:
                output3.clear_output()
                valmin,valmax = val
                fig, ax = plt.subplots(figsize=(l, a))
                cax = ax.imshow(sig * variable, cmap="jet", origin="lower", extent=[x.min(), x.max(), y.min(), y.max()], vmin=valmin, vmax=valmax)
                fig.colorbar(cax, label=f"Valores de {op}")
                ax.contour(x, y, variable, levels=[0], colors='black', linewidths=1)
                ax.set_title(f"Variable '{op}'")
                ax.set_xlabel("Longitud" if coord else "X")
                ax.set_ylabel("Latitud" if coord else "Y")
                plt.show()
                fig_container["fig"] = fig

        def guardar_plot(b):
            fig = fig_container["fig"]
            if fig:
                fig.savefig("plots/mapa_variable_2d.png")
                with output3:
                    print("Imagen guardada como 'mapa_variable_2d.png' en la carpeta plots.")
            else:
                with output3:
                    print("Primero genera el gráfico antes de guardar.")
                    
        boton.on_click(guardar_plot)
        widgets.interactive(actualizar,val=slider1)
                
        output.clear_output()
        dp(widgets.VBox([slider1,boton])) 


def mostmap3d(variable, op, x, y, t, coord):
    #en esta funcion mostramos las variables 3d junto con un slider para poder ver la evolucion temporal
    
    with output2:
        output2.clear_output()
        
        if (len(x)*len(y)>4000000):
            #print(len(x),len(y))
            x ,y , variable = comprimir3d(x, y, variable, 4)
            #print(len(x),len(y))
        
        num_capas = len(t)
        if variable.shape[0]>20:
            vmin = np.nanmin(variable[0:20])
            vmax = np.nanmax(variable[0:20])
        else:
            vmin = np.nanmin(variable)
            vmax = np.nanmax(variable)
        slider = widgets.IntSlider(min=0, max=num_capas - 1, step=1, value=0, description="Time Step",
            layout=widgets.Layout(width="700px"))
        
        escala = np.linspace(int(vmin), int(vmax),50)
        escala_ent = [round(num) for num in escala]
        
        slider1 = widgets.SelectionRangeSlider(
            options=escala_ent,  # Convertir el array en una lista para que sea compatible con options
            index=(0, len(escala) - 1),  # Establecer el valor inicial del rango (por ejemplo, de vmin a vmax)
            description='Scale:',
            disabled=False,
            layout=widgets.Layout(width="700px")
        )
        boton1 = widgets.Button(description="Save Figure")
        boton2 = widgets.Button(description="Save animated GIF")
        fig_container = {"fig": None}
        
        def actualizar(capa,val):
            with output3:
                output3.clear_output()  
                a = 8
                l = int(a*len(x)/len(y))
                
                valmin,valmax = val

                fig, ax = plt.subplots(figsize=(l, a))
                cax = ax.imshow(variable[capa], cmap="jet", origin="lower", extent=[x.min(), x.max(), y.min(), y.max()], vmin=valmin, vmax=valmax)
                fig.colorbar(cax, label=f"Valores de {op}")
                ax.set_title(f"Variable '{op}' - t = {int(t[capa])} s")
                ax.set_xlabel("Longitud" if coord else "X")
                ax.set_ylabel("Latitud" if coord else "Y")
                plt.show()
                fig_container["fig"] = fig

        def guardar_plot(b):
            fig = fig_container["fig"]
            if fig:
                fig.savefig("plots/mapa_variable_3d.png")
                with output3:
                    print("Imagen guardada como 'mapa_variable_3d.png' en la carpeta plots.")

        def guardar_gif(b):
            with output3:
                output3.clear_output()
                print("Generando imágenes para el GIF...")
        
                ruta = "frames_gif"
                os.makedirs(ruta, exist_ok=True)
                archivos = []
        
                valmin, valmax = slider1.value
        
                a = 8
                l = int(a*len(x)/len(y))
        
                for i in range(len(t)):
                    fig, ax = plt.subplots(figsize=(l, a))
                    cax = ax.imshow(variable[i], cmap="jet", origin="lower",
                                    extent=[x.min(), x.max(), y.min(), y.max()],
                                    vmin=valmin, vmax=valmax)
                    fig.colorbar(cax, label=f"Valores de {op}")
                    ax.set_title(f"Variable '{op}' - t = {int(t[i])} s")
                    ax.set_xlabel("Longitud" if coord else "X")
                    ax.set_ylabel("Latitud" if coord else "Y")
                    
                    nombre_archivo = os.path.join(ruta, f"frame_{i:03d}.png")
                    fig.savefig(nombre_archivo)
                    archivos.append(nombre_archivo)
                    plt.close(fig)
        
                print("Creando GIF...")
                imagenes = [imageio.imread(archivo) for archivo in archivos]
                imageio.mimsave("plots/animacion_variable.gif", imagenes, duration=0.01)
        
                print("GIF guardado en la carpeta 'plots' como 'animacion_variable.gif'")
        
                
                for archivo in os.listdir(ruta):
                    archivo_path = os.path.join(ruta, archivo)
                    if os.path.isfile(archivo_path):
                        os.remove(archivo_path)
                
                os.rmdir(ruta)
                    
        boton1.on_click(guardar_plot)
        boton2.on_click(guardar_gif)
        widgets.interactive(actualizar, capa=slider, val=slider1)
                
        dp(widgets.VBox([slider,slider1,boton1,boton2])) 
        output.clear_output()
        

def mostrar_inun(eta_ajust ,bathy ,op ,x ,y ,t ,coord):
     with output2:
        
        num_capas = len(t)

        vmin = -1
        vmax = 2
        slider = widgets.IntSlider(min=0, max=num_capas - 1, step=1, value=0, description="Time Step",
                                    layout=widgets.Layout(width="700px"))
        
        
        boton1 = widgets.Button(description="Guardar imagen")
        boton2 = widgets.Button(description="Guardar GIF")
        fig_container = {"fig": None}
        
        def actualizar(capa):
            with output3:
                output3.clear_output()  
                a = 8
                l = int(a*len(x)/len(y))
                
                fig, ax = plt.subplots(figsize=(l, a))
                cax = ax.imshow(eta_ajust[capa], cmap="jet", origin="lower", extent=[x.min(), x.max(), y.min(), y.max()],vmin=vmin, vmax=vmax)
                fig.colorbar(cax, label=f"Valores de {op}")
                ax.contour(x, y, bathy, levels=[0], colors='black', linewidths=1)
                ax.set_title(f"Variable '{op}' - t = {int(t[capa])} s")
                ax.set_xlabel("Longitud" if coord else "X")
                ax.set_ylabel("Latitud" if coord else "Y")

                plt.show()

                fig_container["fig"] = fig

        def guardar_plot(b):
            fig = fig_container["fig"]
            if fig:
                fig.savefig("plots/mapa_variable_3d.png")
                with output3:
                    print("Imagen guardada como 'mapa_variable_3d.png' en la carpeta plots.")

        def guardar_gif(b):
            with output3:
                output3.clear_output()
                print("Generando imágenes para el GIF...")
        
                ruta = "frames_gif"
                os.makedirs(ruta, exist_ok=True)
                archivos = []
        
                a = 8
                l = int(a*len(x)/len(y))
        
                for i in range(len(t)):
                    fig, ax = plt.subplots(figsize=(l, a))
                    cax = ax.imshow(eta_ajust[i], cmap="jet", origin="lower", extent=[x.min(), x.max(), y.min(), y.max()],vmin=vmin, vmax=vmax)
                    fig.colorbar(cax, label=f"Valores de {op}")
                    ax.contour(x, y, bathy, levels=[0], colors='black', linewidths=1)
                    ax.set_title(f"Variable '{op}' - t = {int(t[i])} s")
                    ax.set_xlabel("Longitud" if coord else "X")
                    ax.set_ylabel("Latitud" if coord else "Y")
                    
                    nombre_archivo = os.path.join(ruta, f"frame_{i:03d}.png")
                    fig.savefig(nombre_archivo)
                    archivos.append(nombre_archivo)
                    plt.close(fig)
        
                print("Creando GIF...")
                imagenes = [imageio.imread(archivo) for archivo in archivos]
                imageio.mimsave("plots/animacion_variable.gif", imagenes, duration=6)
        
                print("GIF guardado en la carpeta 'plots' como 'animacion_variable.gif'")
        
                # Limpiar imágenes temporales
                for archivo in archivos:
                    os.remove(archivo)
                os.rmdir(ruta)
                    
        boton1.on_click(guardar_plot)
        boton2.on_click(guardar_gif)
        
        widgets.interactive(actualizar, capa=slider)
        
        dp(widgets.VBox([slider,boton1,boton2]))  
        output.clear_output()

'''
def zoom2d(variable, x, y, op,coord):
    #en esta funcion para hacer zoom, enseñamos la cuadricula de posibilidades 
    #y mostramos los botones con los que se ejecuta la funcion que hace el zoom
    
    with output2:
        output0.clear_output()
        output2.clear_output()
        output3.clear_output()
        output4.clear_output()
        
        n = 2  # Dividir la imagen en n x n secciones

        a = 10
        l = int(a*len(x)/len(y))
        
        fig, ax = plt.subplots(figsize=(l, a))
        img = ax.imshow(-variable, cmap="jet", origin="lower", 
                        extent=[x.min(), x.max(), y.min(), y.max()], aspect="auto")
        plt.colorbar(img, ax=ax, label=f"Valores de {op}")
        
        xticks = np.linspace(x.min(), x.max(), n+1)
        yticks = np.linspace(y.min(), y.max(), n+1)
        plt.contour(x, y, variable, levels=[0], colors='black', linewidths=1)
        ax.set_xticks(xticks)
        ax.set_yticks(yticks)
        ax.grid(True, linestyle="--", linewidth=0.7, color="white")
        
        # Numerar las celdas
        for i in range(n):
            for j in range(n):
                cx = (xticks[j] + xticks[j+1]) / 2
                cy = (yticks[i] + yticks[i+1]) / 2
                ax.text(cx, cy, f"{i*n + j + 1}", fontsize=14, ha="center", va="center", color="red", weight="bold")
        
        plt.show()
        botones = []
        for i in range(0,n*n):
            boton = widgets.Button(description=str(i+1))
            boton.on_click(lambda b, num=i+1: manejar_recuadro_zoom2d(num, variable, x, y, op, n))  
            botones.append(boton)
        dp(widgets.HBox(botones)) 
              

def manejar_recuadro_zoom2d(recuadro, variable, x, y, op, n):
    #en esta funcion calculamos que partes del array se deben de mostrar en funcion de lo elegido en los botones
    # y llamamos a menu2d para mostrarlo
    
    with output2:
        output2.clear_output()

        coordy = (recuadro - 1) // n  
        coordx = (recuadro - 1) % n 

        yini = (variable.shape[0] // n) * coordy
        yfin = (variable.shape[0] // n) * (coordy + 1) - 1
        xini = (variable.shape[1] // n) * coordx
        xfin = (variable.shape[1] // n) * (coordx + 1) - 1

        menu2d (variable[yini:yfin, xini:xfin], op, x[xini:xfin], y[yini:yfin],True)


def zoom_inun(eta_ajust ,bathy ,op ,x ,y ,t ,coord):
    #en esta funcion para hacer zoom, enseñamos la cuadricula de posibilidades 
    #podemos usar un slider para ver la evolucion temporal
    #y mostramos los botones con los que se ejecuta la funcion que hace el zoom
    
    
    with output2:
        output0.clear_output()
        output2.clear_output()
        output3.clear_output()
        output4.clear_output()
        
        num_capas = len(t)
        n = 2 
        vmin = -1
        vmax = 2
        
        slider = widgets.IntSlider(min=0, max=num_capas - 1, step=1, value=0, description="Time step",
    layout=widgets.Layout(width="700px"))
        
        
        def actualizar(capa):
            with output3:
                output3.clear_output()  
                  
                a = 10
                l = int(a*len(x)/len(y))
            
                fig, ax = plt.subplots(figsize=(l, a))
                img = ax.imshow(eta_ajust[capa], cmap="jet", origin="lower", 
                                extent=[x.min(), x.max(), y.min(), y.max()],vmin=vmin, vmax=vmax, aspect="auto")
                plt.colorbar(img, ax=ax, label=f"Valores de {op}")
                plt.contour(x, y, bathy, levels=[0], colors='black', linewidths=1)
                xticks = np.linspace(x.min(), x.max(), n+1)
                yticks = np.linspace(y.min(), y.max(), n+1)
                ax.set_xticks(xticks)
                ax.set_yticks(yticks)
                ax.grid(True, linestyle="--", linewidth=0.7, color="black")
                
               
                for i in range(n):
                    for j in range(n):
                        cx = (xticks[j] + xticks[j+1]) / 2
                        cy = (yticks[i] + yticks[i+1]) / 2
                        ax.text(cx, cy, f"{i*n + j + 1}", fontsize=14, ha="center", va="center", color="red", weight="bold")
                
                plt.show()

        widgets.interactive(actualizar, capa=slider)
                
        dp(widgets.VBox([slider]))  
        
        botones = []
        for i in range(0,n*n):
            boton = widgets.Button(description=str(i+1))
            boton.on_click(lambda b, num=i+1: manejar_recuadro_zoom_inun(num, eta_ajust, bathy, x, y, t, op, n))  
            botones.append(boton)
        dp(widgets.HBox(botones)) 
        

def manejar_recuadro_zoom_inun(recuadro, variable, bathy, x, y, t, op, n):
    #en esta funcion hacemos los calculos para obtener las posiciones a mostrar
    
    with output2:
        output2.clear_output()
        print("se esta ejecutando la funcion ")
        coordy = (recuadro - 1) // n  
        coordx = (recuadro - 1) % n 
        
        yini = (variable.shape[1] // n) * coordy
        yfin = (variable.shape[1] // n) * (coordy + 1) - 1
        xini = (variable.shape[2] // n) * coordx
        xfin = (variable.shape[2] // n) * (coordx + 1) - 1
        
      
        menu_inun (variable[:,yini:yfin, xini:xfin],bathy[yini:yfin, xini:xfin], op, x[xini:xfin], y[yini:yfin], t,True)

def zoom3d(variable, op, x, y, t, coord):
    #en esta funcion para hacer zoom, enseñamos la cuadricula de posibilidades 
    #podemos usar un slider para ver la evolucion temporal
    #y mostramos los botones con los que se ejecuta la funcion que hace el zoom
    
    
    with output2:
        output0.clear_output()
        output2.clear_output()
        output3.clear_output()
        output4.clear_output()
        
        num_capas = len(t)
        n = 2 
        vmin = -1
        vmax = 2
        
        slider = widgets.IntSlider(min=0, max=num_capas - 1, step=1, value=0, description="Time Step",
    layout=widgets.Layout(width="700px"))
        
        
        def actualizar(capa):
            with output3:
                output3.clear_output()  
                  
                a = 10
                l = int(a*len(x)/len(y))
            
                fig, ax = plt.subplots(figsize=(l, a))
                img = ax.imshow(variable[capa], cmap="jet", origin="lower", 
                                extent=[x.min(), x.max(), y.min(), y.max()],vmin=vmin, vmax=vmax, aspect="auto")
                plt.colorbar(img, ax=ax, label=f"Valores de {op}")
                
                xticks = np.linspace(x.min(), x.max(), n+1)
                yticks = np.linspace(y.min(), y.max(), n+1)
                ax.set_xticks(xticks)
                ax.set_yticks(yticks)
                ax.grid(True, linestyle="--", linewidth=0.7, color="white")
                
               
                for i in range(n):
                    for j in range(n):
                        cx = (xticks[j] + xticks[j+1]) / 2
                        cy = (yticks[i] + yticks[i+1]) / 2
                        ax.text(cx, cy, f"{i*n + j + 1}", fontsize=14, ha="center", va="center", color="red", weight="bold")
                
                plt.show()

        widgets.interactive(actualizar, capa=slider)
                
        dp(widgets.VBox([slider]))  
        
        botones = []
        for i in range(0,n*n):
            boton = widgets.Button(description=str(i+1))
            boton.on_click(lambda b, num=i+1: manejar_recuadro_zoom3d(num, variable, x, y, t, op, n))  
            botones.append(boton)
        dp(widgets.HBox(botones)) 
        

def manejar_recuadro_zoom3d(recuadro, variable, x, y, t, op, n):
    #en esta funcion hacemos los calculos para obtener las posiciones a mostrar
    
    with output2:
        output2.clear_output()
        
        coordy = (recuadro - 1) // n  
        coordx = (recuadro - 1) % n 
        
        yini = (variable.shape[1] // n) * coordy
        yfin = (variable.shape[1] // n) * (coordy + 1) - 1
        xini = (variable.shape[2] // n) * coordx
        xfin = (variable.shape[2] // n) * (coordx + 1) - 1

        menu3d (variable[:,yini:yfin, xini:xfin], op, x[xini:xfin], y[yini:yfin], t,True)
'''
    
def seccion2d(variable, x, y, op, coord):
    # Esta función muestra una cuadricula con botones para elegir una sección,
    # y añade campos de texto para introducir coordenadas manualmente

    global val1
    global val2
    val1 = None
    val2 = None

    with output2:
        output0.clear_output()
        output2.clear_output()
        output3.clear_output()
        output4.clear_output()
        
        n = 5  # Dividir la imagen en n x n secciones

        a = 10
        l = int(a * len(x) / len(y))

        fig, ax = plt.subplots(figsize=(l, a))
        img = ax.imshow(-variable, cmap="jet", origin="lower", 
                        extent=[x.min(), x.max(), y.min(), y.max()], aspect="auto")
        plt.colorbar(img, ax=ax, label=f"Valores de {op}")

        xticks = np.linspace(x.min(), x.max(), n + 1)
        yticks = np.linspace(y.min(), y.max(), n + 1)
        ax.set_xticks(xticks)
        ax.set_yticks(yticks)
        ax.grid(True, linestyle="--", linewidth=0.7, color="white")
        plt.contour(x, y, variable, levels=[0], colors='black', linewidths=1)

        # Numerar las celdas
        for i in range(n):
            for j in range(n):
                cx = (xticks[j] + xticks[j + 1]) / 2
                cy = (yticks[i] + yticks[i + 1]) / 2
                ax.text(cx, cy, f"{i * n + j + 1}", fontsize=14, ha="center", va="center", color="red", weight="bold")

        plt.show()

        # Crear botones por cuadrícula
        botones = []
        for i in range(0, n * n):
            boton = widgets.Button(description=str(i + 1), layout=widgets.Layout(width='40px', height='30px'))
            boton.on_click(lambda b, num=i + 1: act_valores_seccion2d(num, variable, x, y, op, n))
            botones.append(boton)

        tabla = widgets.GridBox(
            botones,
            layout=widgets.Layout(
                grid_template_columns="repeat(5, 60px)",
                grid_gap="0px",
                width="1200px",
            )
        )

        # Campos de texto para coordenadas manuales
        x1 = widgets.FloatText(description="X1:", layout=widgets.Layout(width="250px"))
        y1 = widgets.FloatText(description="Y1:", layout=widgets.Layout(width="250px"))
        x2 = widgets.FloatText(description="X2:", layout=widgets.Layout(width="250px"))
        y2 = widgets.FloatText(description="Y2:", layout=widgets.Layout(width="250px"))
        boton_coord = widgets.Button(description="Continuar", button_style='success')

        mensaje_salida = widgets.Output()

        def verificar_coordenadas(b):
            mensaje_salida.clear_output()
            x_val1 = x1.value
            y_val1 = y1.value
            x_val2 = x2.value
            y_val2 = y2.value

            dentro_x1 = x.min() <= x_val1 <= x.max()
            dentro_y1 = y.min() <= y_val1 <= y.max()
            dentro_x2 = x.min() <= x_val2 <= x.max()
            dentro_y2 = y.min() <= y_val2 <= y.max()

            with mensaje_salida:
                if dentro_x1 and dentro_y1 and dentro_x2 and dentro_y2:
                    print(f"Coordenadas válidas: X1 = {x_val1}, Y1 = {y_val1} y X2 = {x_val2}, Y1 = {y_val2}")
                    manejar_coordenadas_seccion2d(x_val1, y_val1, x_val2, y_val2, variable, x, y, op)
                    # Aquí podrías continuar con el procesamiento, como definir val1 o val2
                else:
                    print("❌ Coordenadas fuera del dominio.")
                    if not dentro_x1:
                        print(f" - X1={x_val1} está fuera del rango ({x.min():.2f}, {x.max():.2f})")
                    if not dentro_y1:
                        print(f" - Y1={y_val1} está fuera del rango ({y.min():.2f}, {y.max():.2f})")
                    if not dentro_x2:
                        print(f" - X2={x_val2} está fuera del rango ({x.min():.2f}, {x.max():.2f})")
                    if not dentro_y2:
                        print(f" - Y2={y_val2} está fuera del rango ({y.min():.2f}, {y.max():.2f})")

        boton_coord.on_click(verificar_coordenadas)

        entrada_manual = widgets.HBox([x1,y1,x2,y2, boton_coord])
        dp(widgets.VBox([tabla, entrada_manual, mensaje_salida]))

        
def act_valores_seccion2d(num, variable, x, y, op, n):
    #en esta funcion recibimos los dos valores seleccionados por los botones
    #de la funcion anterior y llamamos a la funcion para representarlos
    
    global val1
    global val2
    with output2:
        if val1 == None and val2 == None:
            val1 = num
            
        elif val2 == None:
            val2 = num 
            
            manejar_recuadro_seccion2d(val1, val2, variable, x, y, op, n)

        
def manejar_recuadro_seccion2d(recuadro1, recuadro2, variable, x, y, op, n):
    # Calculamos la sección a mostrar y añadimos una línea negra en el corte 
    
    with output2:
        output2.clear_output()
        
        
        coordy1 = (recuadro1 - 1) // n  
        coordx1 = (recuadro1 - 1) % n 
        coy1 = (variable.shape[0]//n)*(2*coordy1+1)//2
        cox1 = (variable.shape[1]//n)*(2*coordx1+1)//2
    
        coordy2 = (recuadro2 - 1) // n  
        coordx2 = (recuadro2 - 1) % n 
        coy2 = (variable.shape[0]//n)*(2*coordy2+1)//2
        cox2 = (variable.shape[1]//n)*(2*coordx2+1)//2


        a = 10
        l = int(a*len(x)/len(y))
        
        fig, ax = plt.subplots(figsize=(l, a))
        img = ax.imshow(variable, cmap="jet", origin="lower",
                        extent=[x.min(), x.max(), y.min(), y.max()], aspect="auto")
        plt.colorbar(img, ax=ax, label=f"Valores de {op}")
        plt.contour(x, y, variable, levels=[0], colors='black', linewidths=1)
       
        ax.plot([x[cox1], x[cox2]], [y[coy1], y[coy2]], 'k-', linewidth=2, label="Corte seleccionado")

        ax.set_xlabel("Eje X")
        ax.set_ylabel("Eje Y")
        ax.legend()
        plt.title("Mapa con la sección seleccionada")
        plt.show()
        
        
        fx = []
        if np.abs(cox2 - cox1) > np.abs(coy2 - coy1):  
            x_range = range(min(cox1, cox2), max(cox1, cox2))  
            for i in x_range:
                y_interp = int((i - cox1) * (coy2 - coy1) / (cox2 - cox1) + coy1)
                fx.append(variable[y_interp, i])
            a = np.array(list(x_range))
        else: 
            y_range = range(min(coy1, coy2), max(coy1, coy2))  
            for i in y_range:
                x_interp = int((i - coy1) * (cox2 - cox1) / (coy2 - coy1) + cox1)
                fx.append(variable[i, x_interp])
            a = np.array(list(y_range))
        
        
        plt.figure(figsize=(7, 5))
        plt.plot(a, fx, linestyle='-', color='blue')
        plt.xlabel("Coordenada")
        plt.ylabel(f"Valores de {op}")
        plt.title("Perfil seleccionado")
        plt.grid()
        plt.show()

       
        boton1 = widgets.Button(description="Continuar")
        boton1.on_click(lambda b: menu2d(variable, op, x, y, True))
        dp(widgets.HBox([boton1]))

def manejar_coordenadas_seccion2d(x1, y1, x2, y2, variable, x, y, op):
    with output2:
        output2.clear_output()
        output3.clear_output()
        output4.clear_output()
        

        # Verificar que las coordenadas estén dentro del dominio
        if not (x.min() <= x1 <= x.max() and x.min() <= x2 <= x.max() and 
                y.min() <= y1 <= y.max() and y.min() <= y2 <= y.max()):
            print("❌ Las coordenadas están fuera del dominio.")
            return

        # Crear grilla regular a partir de x, y
        xx, yy = np.meshgrid(x, y)

        a = 10
        l = int(a * len(x) / len(y))

        # Mostrar mapa con línea negra del corte
        fig, ax = plt.subplots(figsize=(l, a))
        img = ax.imshow(-variable, cmap="jet", origin="lower",
                        extent=[x.min(), x.max(), y.min(), y.max()], aspect="auto")
        plt.colorbar(img, ax=ax, label=f"Valores de {op}")
        plt.contour(x, y, variable, levels=[0], colors='black', linewidths=1)

        ax.plot([x1, x2], [y1, y2], 'k-', linewidth=2, label="Corte seleccionado")
        ax.set_xlabel("Eje X")
        ax.set_ylabel("Eje Y")
        ax.legend()
        plt.title("Mapa con la sección seleccionada")
        plt.show()

        n_points = 200
        xi = np.linspace(x1, x2, n_points)
        yi = np.linspace(y1, y2, n_points)

        # Convertir coordenadas a índices más cercanos en x, y
        def coord_a_indice(valor, vector):
            return np.argmin(np.abs(vector - valor))

        valores = []
        distancias = []

        for i in range(n_points):
            ix = coord_a_indice(xi[i], x)
            iy = coord_a_indice(yi[i], y)
            val = -variable[iy, ix]  # cuidado con el orden: [y, x]
            valores.append(val)
            d = np.sqrt((xi[i] - x1)**2 + (yi[i] - y1)**2)
            distancias.append(d)

        plt.figure(figsize=(7, 5))
        plt.plot(distancias, valores, linestyle='-', color='blue')
        plt.xlabel("Distancia a lo largo del corte")
        plt.ylabel(f"Valores de {op}")
        plt.title("Perfil seleccionado")
        plt.grid()
        plt.show()

        # Botón para continuar
        boton1 = widgets.Button(description="Continuar")
        boton1.on_click(lambda b: menu2d(variable, op, x, y, True))
        dp(widgets.HBox([boton1]))



def sectemporal(variable, op, x, y, t, coord):
    output0.clear_output()
    output2.clear_output()
    output3.clear_output()
    output4.clear_output()

    with output2:
        num_capas = len(t)
        n = 5

        slider = widgets.IntSlider(min=0, max=num_capas - 1, step=1, value=0, description="Time Step")

        def actualizar(capa):
            with output3:
                output3.clear_output()
                a = 10
                l = int(a * len(x) / len(y))

                fig, ax = plt.subplots(figsize=(l, a))
                img = ax.imshow(variable[capa], cmap="jet", origin="lower",
                                extent=[x.min(), x.max(), y.min(), y.max()], aspect="auto")
                plt.colorbar(img, ax=ax, label=f"Valores de {op}")

                xticks = np.linspace(x.min(), x.max(), n + 1)
                yticks = np.linspace(y.min(), y.max(), n + 1)
                ax.set_xticks(xticks)
                ax.set_yticks(yticks)
                ax.grid(True, linestyle="--", linewidth=0.7, color="white")

                for i in range(n):
                    for j in range(n):
                        cx = (xticks[j] + xticks[j + 1]) / 2
                        cy = (yticks[i] + yticks[i + 1]) / 2
                        ax.text(cx, cy, f"{i * n + j + 1}", fontsize=14, ha="center", va="center", color="red", weight="bold")

                plt.show()

        interactivo = widgets.interactive_output(actualizar, {"capa": slider})

        # ------------------ Sección de botones ------------------
        botones = []
        for i in range(n * n):
            boton = widgets.Button(description=str(i + 1))
            boton.on_click(lambda b, num=i + 1: manejar_recuadro_sectemporal(num, variable, x, y, t, op, n, coord))
            botones.append(boton)

        print("Elige una celda para obtener una secuencia temporal")

        tabla = widgets.GridBox(
            botones,
            layout=widgets.Layout(
                grid_template_columns="repeat(5, 150px)",
                grid_gap="0px",
                width="3000px",
            )
        )

        # ------------------ Entrada manual ------------------
        x1 = widgets.FloatText(description="X:", layout=widgets.Layout(width="250px"))
        y1 = widgets.FloatText(description="Y:", layout=widgets.Layout(width="250px"))
        boton_coord = widgets.Button(description="Continuar", button_style='success')
        mensaje_salida = widgets.Output()

        def verificar_coordenadas(b):
            mensaje_salida.clear_output()
            x_val1 = x1.value
            y_val1 = y1.value

            dentro_x1 = x.min() <= x_val1 <= x.max()
            dentro_y1 = y.min() <= y_val1 <= y.max()

            with mensaje_salida:
                if dentro_x1 and dentro_y1:
                    print(f"Coordenadas válidas: X = {x_val1}, Y = {y_val1}")
                    manejar_coordenada_sectemporal(x_val1, y_val1, variable, x, y, t, op, coord)
                else:
                    print("❌ Coordenadas fuera del dominio.")
                    if not dentro_x1:
                        print(f" - X1={x_val1} está fuera del rango ({x.min():.2f}, {x.max():.2f})")
                    if not dentro_y1:
                        print(f" - Y1={y_val1} está fuera del rango ({y.min():.2f}, {y.max():.2f})")

        boton_coord.on_click(verificar_coordenadas)
        entrada_manual = widgets.HBox([x1, y1, boton_coord])

        # Mostrar todo junto
        dp(widgets.VBox([slider, interactivo, tabla, entrada_manual, mensaje_salida]))




def manejar_recuadro_sectemporal(recuadro, variable, x, y, t, op, n, coord):
   with output3:
        output3.clear_output()
        
        coordy = (recuadro - 1) // n  
        coordx = (recuadro - 1) % n 
        coy = (variable.shape[1]//n)*(2*coordy+1)//2
        cox = (variable.shape[2]//n)*(2*coordx+1)//2

        plt.figure(figsize=(8, 4))
        plt.plot(t, variable[:, coy, cox], marker="o", linestyle="-", color="b")
        plt.xlabel("Tiempo")
        plt.ylabel(f"Valores de {op}")
        plt.title(f"Evolución temporal en la celda {recuadro}")
        plt.grid()
        plt.show()

        boton1 = widgets.Button(description= "continuar")
        boton1.on_click(lambda b: menu3d(variable, op, x, y, t, coord))
        dp(widgets.HBox([boton1]))

def manejar_coordenada_sectemporal(xval, yval, variable, x, y, t, op, coord):
    with output3:
        output3.clear_output()

        # Convertir coordenadas a índices más cercanos en la malla
        def coord_a_indice(valor, vector):
            return np.argmin(np.abs(vector - valor))

        cox = coord_a_indice(xval, x)
        coy = coord_a_indice(yval, y)

        # Graficar evolución temporal en esa coordenada
        plt.figure(figsize=(8, 4))
        plt.plot(t, variable[:, coy, cox], marker="o", linestyle="-", color="b")
        plt.xlabel("Tiempo")
        plt.ylabel(f"Valores de {op}")
        plt.title(f"Evolución temporal en ({xval:.3f}, {yval:.3f})")
        plt.grid()
        plt.show()

        # Botón para volver
        boton1 = widgets.Button(description="Continuar")
        boton1.on_click(lambda b: menu3d(variable, op, x, y, t, coord))
        dp(widgets.HBox([boton1]))


def sectemporal_boya(variable, t, op):
   with output3:
        output3.clear_output()
        
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(t, variable, marker="o", linestyle="-", color="b")
        ax.set_xlabel("Tiempo")
        ax.set_ylabel(f"Valores de la {op} [m]")
        ax.set_title("Evolución temporal en la boya seleccionada")
        ax.grid()

        plt.show()
       
        boton = widgets.Button(description="Guardar imagen")
        boton_guardar.on_click(guardar_plot)
        dp(boton)
        def guardar_plot(b):
            fig.savefig("grafico_boya.png")
            print("Imagen guardada como 'grafico_boya.png'")

        
        

        

def seccion3d(variable, op, x, y, t, coord):
    
    output0.clear_output()
    output2.clear_output()
    output3.clear_output()
    output4.clear_output()
    with output2:
        
        num_capas = len(t)
        n = 5 
        
        slider = widgets.IntSlider(min=0, max=num_capas - 1, step=1, value=0, description="Time Step",
    layout=widgets.Layout(width="700px"))
        
        
        def actualizar(capa):
            with output3:
                output3.clear_output()   
                a = 10
                l = int(a*len(x)/len(y))
                
                fig, ax = plt.subplots(figsize=(l, a))
                img = ax.imshow(variable[capa], cmap="jet", origin="lower", 
                                extent=[x.min(), x.max(), y.min(), y.max()], aspect="auto")
                plt.colorbar(img, ax=ax, label=f"Valores de {op}")
                
                xticks = np.linspace(x.min(), x.max(), n+1)
                yticks = np.linspace(y.min(), y.max(), n+1)
                ax.set_xticks(xticks)
                ax.set_yticks(yticks)
                ax.grid(True, linestyle="--", linewidth=0.7, color="white")
                
               
                for i in range(n):
                    for j in range(n):
                        cx = (xticks[j] + xticks[j+1]) / 2
                        cy = (yticks[i] + yticks[i+1]) / 2
                        ax.text(cx, cy, f"{i*n + j + 1}", fontsize=14, ha="center", va="center", color="red", weight="bold")
                
                plt.show()

        widgets.interactive(actualizar, capa=slider)        
                
        dp(widgets.VBox([slider]))  

        boton1 = widgets.Button(description= "grabar tiempo")
        boton1.on_click(lambda b: selec_seccion3d(variable[slider.value], op, x, y, t, coord))
        dp(widgets.HBox([boton1]))

        
def selec_seccion3d(variable, op, x, y, t, coord):
    #en esta funcion mostramos una cuadricula con las posibilidades para elegir 
    #punto inicial y punto final en nuestra seccion y mostramos una cuadricula de 
    #botones para seleccionar estos valores 
    
    global val1
    global val2
    val1 = None
    val2 = None
    
    with output2:
        output2.clear_output()
        
        n = 5  # Dividir la imagen en n x n secciones
        a = 10
        l = int(a*len(x)/len(y))
            
        fig, ax = plt.subplots(figsize=(l, a))
        img = ax.imshow(variable, cmap="jet", origin="lower", 
                        extent=[x.min(), x.max(), y.min(), y.max()], aspect="auto")
        plt.colorbar(img, ax=ax, label=f"Valores de {op}")
        
        xticks = np.linspace(x.min(), x.max(), n+1)
        yticks = np.linspace(y.min(), y.max(), n+1)
        ax.set_xticks(xticks)
        ax.set_yticks(yticks)
        ax.grid(True, linestyle="--", linewidth=0.7, color="white")
        
        # Numerar las celdas
        for i in range(n):
            for j in range(n):
                cx = (xticks[j] + xticks[j+1]) / 2
                cy = (yticks[i] + yticks[i+1]) / 2
                ax.text(cx, cy, f"{i*n + j + 1}", fontsize=14, ha="center", va="center", color="red", weight="bold")
        
        plt.show()
        botones = []
        for i in range(0,n*n):
            boton = widgets.Button(description=str(i+1),layout=widgets.Layout(width='40px', height='30px'))
            boton.on_click(lambda b, num=i+1: act_valores_seccion2d(num, variable, x, y, op, n))  
            botones.append(boton)
        tabla = widgets.GridBox(
            botones,
            layout=widgets.Layout(
                grid_template_columns="repeat(5, 50px)",  # 3 columnas de 200px
                grid_gap="0px",  # Espacio entre botones
                width="1200px",  # Ancho total del contenedor
            )
        )
# Campos de texto para coordenadas manuales
        x1 = widgets.FloatText(description="X1:", layout=widgets.Layout(width="250px"))
        y1 = widgets.FloatText(description="Y1:", layout=widgets.Layout(width="250px"))
        x2 = widgets.FloatText(description="X2:", layout=widgets.Layout(width="250px"))
        y2 = widgets.FloatText(description="Y2:", layout=widgets.Layout(width="250px"))
        boton_coord = widgets.Button(description="Continuar", button_style='success')

        mensaje_salida = widgets.Output()

        def verificar_coordenadas(b):
            mensaje_salida.clear_output()
            x_val1 = x1.value
            y_val1 = y1.value
            x_val2 = x2.value
            y_val2 = y2.value

            dentro_x1 = x.min() <= x_val1 <= x.max()
            dentro_y1 = y.min() <= y_val1 <= y.max()
            dentro_x2 = x.min() <= x_val2 <= x.max()
            dentro_y2 = y.min() <= y_val2 <= y.max()

            with mensaje_salida:
                if dentro_x1 and dentro_y1 and dentro_x2 and dentro_y2:
                    print(f"Coordenadas válidas: X1 = {x_val1}, Y1 = {y_val1} y X2 = {x_val2}, Y1 = {y_val2}")
                    manejar_coordenadas_seccion2d(x_val1, y_val1, x_val2, y_val2, variable, x, y, op)
                    # Aquí podrías continuar con el procesamiento, como definir val1 o val2
                else:
                    print("❌ Coordenadas fuera del dominio.")
                    if not dentro_x1:
                        print(f" - X1={x_val1} está fuera del rango ({x.min():.2f}, {x.max():.2f})")
                    if not dentro_y1:
                        print(f" - Y1={y_val1} está fuera del rango ({y.min():.2f}, {y.max():.2f})")
                    if not dentro_x2:
                        print(f" - X2={x_val2} está fuera del rango ({x.min():.2f}, {x.max():.2f})")
                    if not dentro_y2:
                        print(f" - Y2={y_val2} está fuera del rango ({y.min():.2f}, {y.max():.2f})")

        boton_coord.on_click(verificar_coordenadas)

        entrada_manual = widgets.HBox([x1,y1,x2,y2, boton_coord])
        dp(widgets.VBox([tabla, entrada_manual, mensaje_salida]))

def procesar_arrival_times(dataset):
    with output:
        print("Loading...")
    output0.clear_output()

    op = 'arrival_times'
    variable = dataset.variables[op][:]
    if 'lat' in dataset.variables and 'lon' in dataset.variables or 'y' in dataset.variables and 'x' in dataset.variables:
        keyy = 'lat' if 'lat' in dataset.variables else 'y'
        keyx = 'lon' if 'lon' in dataset.variables else 'x'
        x = dataset.variables[keyx][:]
        y = dataset.variables[keyy][:]
        coord = True
    else:
        x = np.arange(variable.shape[1])  
        y = np.arange(variable.shape[0])
        coord = False
    menu2d_arrival_times (variable, op, x, y,coord)

def menu2d_arrival_times(variable, op, x, y, coord):
    output0.clear_output()
    output2.clear_output()
    output3.clear_output()
    output4.clear_output()

    with output2:
        arrival_time = np.where(variable == -1, np.nan, variable)
        arrival_time = arrival_time / 60  # convertir a minutos

        if (len(x)*len(y) > 4000000):
            x, y, arrival_time = comprimir(x, y, arrival_time, 4)

        vmax = np.nanmax(arrival_time)
        vmin = np.nanmin(arrival_time[arrival_time > 0])  # Evitar valores cero o negativos
        vmin = 0 if np.isnan(vmin) else vmin

        escala = np.linspace(vmin, vmax, 50)
        escala_ent = [round(num) for num in escala]

        slider1 = widgets.SelectionRangeSlider(
            options=escala_ent,
            index=(0, len(escala_ent) - 1),
            description='Scale:',
            disabled=False,
            layout=widgets.Layout(width="700px")
        )

        boton = widgets.Button(description="Save figure")
        fig_container = {"fig": None}

        def actualizar(val):
            with output3:
                output3.clear_output()
                valmin, valmax = val
                interval = 5
                max_level = math.ceil(valmax / interval) * interval
                levels = np.arange(0, max_level + interval, interval)

                fig, ax = plt.subplots(figsize=(12, 8))
                extent = [x.min(), x.max(), y.min(), y.max()]
                cmap = plt.cm.get_cmap('coolwarm')
                im = ax.imshow(-arrival_time, origin='lower', extent=extent, cmap=cmap,
                               vmin=-valmax, vmax=-valmin, aspect='auto')

                cs = ax.contour(arrival_time, levels=levels, colors='k', linewidths=0.5,
                                extent=extent, origin='lower')
                ax.clabel(cs, inline=True, fontsize=8, fmt='%d')

                ax.set_xlabel('Longitude' if coord else 'X')
                ax.set_ylabel('Latitude' if coord else 'Y')
                ax.set_title(f"Tsunami Arrival Time ({op})", fontsize=14)
                ax.grid(True, linestyle='--', alpha=0.5)

                cbar = fig.colorbar(im, ax=ax, orientation='horizontal', pad=0.05)
                cbar.set_ticks([-valmax, -((valmax + valmin) / 2), -valmin])
                cbar.set_ticklabels([f"{valmax:.0f}", f"{(valmax + valmin)/2:.0f}", f"{valmin:.0f}"])
                cbar.set_label('Tsunami Arrival Time (min)', fontsize=12)

                plt.show()
                fig_container["fig"] = fig

        def guardar_plot(b):
            fig = fig_container["fig"]
            if fig:
                fig.savefig("plots/mapa_variable_2d.png")
                with output3:
                    print("Imagen guardada como 'mapa_variable_2d.png' en la carpeta plots.")
            else:
                with output3:
                    print("Primero genera el gráfico antes de guardar.")

        boton.on_click(guardar_plot)
        widgets.interactive(actualizar, val=slider1)

        output.clear_output()
        dp(widgets.VBox([slider1, boton]))




def inundacion (dataset,op1,op2):
    output0.clear_output()
    output2.clear_output()
    output3.clear_output()
    output4.clear_output()
    with  output: 
        print("Loading...")
    with output2:
        xmax = dataset.variables[op1].shape[2]
        xmin = 0
        ymax = dataset.variables[op1].shape[1]
        ymin = 0
        tmax = dataset.variables[op1].shape[0]
        tmin = 0
        print(xmax*ymax*tmax)
        if xmax*ymax*tmax > 200000000:
            tmax = 150000000//(xmax*ymax)
        eta = dataset.variables[op1][tmin:tmax,ymin:ymax,xmin:xmax]
        bathy = dataset.variables[op2][ymin:ymax,xmin:xmax]
        
        if eta.ndim != 3 or bathy.ndim != 2:
            print("Error: eta debe ser 3D (time, lat, lon) y bathy 2D (lat, lon)")
            return
        
        
        mask = bathy < 0
        eta_ajust = np.where(mask, eta - bathy, np.nan)
        eta_ajust[eta_ajust < 0] = np.nan
        
        
        if "lon" in dataset.variables and "lat" in dataset.variables and "time" in dataset.variables:
            x = dataset.variables["lon"][xmin:xmax]
            y = dataset.variables["lat"][ymin:ymax]
            t = dataset.variables["time"][tmin:tmax]
            coord = True
        else:
            x = np.arange(eta_ajust.shape[2])
            y = np.arange(eta_ajust.shape[1])
            t = np.arange(eta_ajust.shape[0])
            coord = False
        menu_inun(eta_ajust ,bathy ,"inundacion" ,x ,y ,t ,coord)

selvars = []

def operarvar(nombre):
    selvars.clear()
    with output:
        output.clear_output()
        print("Cargando...")

    archivo = "./simulaciones/" + nombre 

    with output2:
        output0.clear_output()
        output2.clear_output()
        output3.clear_output()
        output4.clear_output()

        try:
            dataset = nc.Dataset(archivo, mode="r")
        except FileNotFoundError:
            print("Error: El archivo no existe.")
            return

        print("Select variables to operate:")    

        botones = []
        for key in dataset.variables:
            if len(dataset.variables[key].shape) > 1:
                toggle = widgets.ToggleButton(description=key, value=False)
                toggle.observe(lambda change, key=key, dataset=dataset, nombre=nombre: cambtog(change, key, dataset, nombre), names='value')
                botones.append(toggle)

        output.clear_output()    
        dp(widgets.VBox(botones)) 

# Controlador de selección de variables (toggle)
def cambtog(change, key, dataset, nombre):
    if change['name'] == 'value':
        if change['new']:
            if len(selvars) < 2:
                selvars.append((key, dataset))
            else:
                # No más de 2 variables
                change['owner'].value = False
                with output4:
                    output4.clear_output()
                    output3.clear_output()
                    print("Solo se pueden seleccionar 2 variables.")
        else:
            # Quitar si se desactiva
            selvars[:] = [item for item in selvars if item[0] != key]

        # Mostrar variables seleccionadas
        with output3:
            output3.clear_output()
            if selvars:
                print("Seleccionadas:")
                for var, _ in selvars:
                    print(f" - {var}")
            if len(selvars) == 2:
                print("Listo para operar con estas dos variables.")
                mostrar_operadores(nombre)

# Muestra los botones de operación + - *
def mostrar_operadores(nombre):
    with output4:
        output4.clear_output()
        operador_suma = widgets.Button(description='+', layout=widgets.Layout(width='50px'))
        operador_resta = widgets.Button(description='-', layout=widgets.Layout(width='50px'))
        operador_multi = widgets.Button(description='*', layout=widgets.Layout(width='50px'))

        operador_suma.on_click(lambda b: operar_con_dos_variables(selvars[0], selvars[1], '+', nombre))
        operador_resta.on_click(lambda b: operar_con_dos_variables(selvars[0], selvars[1], '-', nombre))
        operador_multi.on_click(lambda b: operar_con_dos_variables(selvars[0], selvars[1], '*', nombre))

        dp(widgets.HBox([operador_suma, operador_resta, operador_multi]))

# Realiza la operación y visualiza resultado
def operar_con_dos_variables(var1, var2, operador, nombre):
    with output:
        output.clear_output()
        print("Loading...")

    with output4:
        output2.clear_output()
        output3.clear_output()
        output4.clear_output()
        print(f"Operando: {var1[0]} {operador} {var2[0]}")
        if 'lat' in var1[1].variables and 'lon' in var1[1].variables or 'y' in var1[1].variables and 'x' in var1[1].variables:
            keyy = 'lat' if 'lat' in var1[1].variables else 'y'
            keyx = 'lon' if 'lon' in var1[1].variables else 'x'
            x = var1[1].variables[keyx][:]
            y = var1[1].variables[keyy][:]
        try:
            data1 = var1[1].variables[var1[0]][:]
            data2 = var2[1].variables[var2[0]][:]

            
            
            if data1.ndim == 2 and data2.ndim == 2:
                d3 = False
            else:
                d3 = True
                t = var2[1].variables['time']
                
            # Si una es 2D y otra 3D, convertir la 2D en 3D
            if data1.ndim == 2 and data2.ndim == 3:
                d3 = True
                xmax = var2[1].variables[var2[0]].shape[2]
                xmin = 0
                ymax = var2[1].variables[var2[0]].shape[1]
                ymin = 0
                tmax = var2[1].variables[var2[0]].shape[0]
                tmin = 0
                print(xmax*ymax*tmax)
                if xmax*ymax*tmax > 200000000:
                    tmax = 150000000//(xmax*ymax)
                data2 = var2[1].variables[var2[0]][tmin:tmax,ymin:ymax,xmin:xmax]
                t = var2[1].variables['time']
                
            elif data1.ndim == 3 and data2.ndim == 2:
                d3 = True
                xmax = var1[1].variables[var1[0]].shape[2]
                xmin = 0
                ymax = var1[1].variables[var1[0]].shape[1]
                ymin = 0
                tmax = var1[1].variables[var1[0]].shape[0]
                tmin = 0
                print(xmax*ymax*tmax)
                if xmax*ymax*tmax > 200000000:
                    tmax = 150000000//(xmax*ymax)
                data2 = var1[1].variables[var1[0]][tmin:tmax,ymin:ymax,xmin:xmax]
                t = var1[1].variables['time']
                
            elif data1.ndim != data2.ndim:
                print("Error: Las dimensiones no son compatibles.")
                return

            # Operación
            if operador == '+':
                result = data1 + data2
            elif operador == '-':
                result = data1 - data2
            elif operador == '*':
                result = data1 * data2
            else:
                print("Operador no válido.")
                return

            print("Resultado calculado. Mostrando visualización...")

            
            if d3:
                mostmap3d(result, f"Operación de: {var1[0]} {operador} {var2[0]}", x, y,t, True)
            else:
                mostmap2d(result, f"Operación de: {var1[0]} {operador} {var2[0]}", x, y, True)

        except Exception as e:
            print(f"Error en la operación: {e}")
        
def comprimir(x, y, u, bloque=4):
    with output:
        print("Comprimiendo array...")
        ny, nx = u.shape
        nx_c = nx - (nx % bloque)
        ny_c = ny - (ny % bloque)
    
        u_c = u[:ny_c, :nx_c]
        x_c = x[:nx_c]
        y_c = y[:ny_c]
    
        u_red = u_c.reshape(ny_c // bloque, bloque, nx_c // bloque, bloque)
        u_red = u_red.mean(axis=(1, 3)) 
        x_red = x_c.reshape(-1, bloque).mean(axis=1)
        y_red = y_c.reshape(-1, bloque).mean(axis=1)
    
        return x_red, y_red, u_red

def comprimir3d(x, y, u, bloque=4):
    with output:
        print("Comprimiendo array...")
        nt, ny, nx = u.shape
        nx_c = nx - (nx % bloque)
        ny_c = ny - (ny % bloque)
    
        u_c = u[:,:ny_c, :nx_c]
        x_c = x[:nx_c]
        y_c = y[:ny_c]
        
        u_red = u_c.reshape(nt, ny_c // bloque, bloque, nx_c // bloque, bloque)
        u_red = u_red.mean(axis=(2, 4)) 
        
        x_red = x_c.reshape(-1, bloque).mean(axis=1)
        y_red = y_c.reshape(-1, bloque).mean(axis=1)
    
        return x_red, y_red, u_red
