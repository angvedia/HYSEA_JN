import os
import ipywidgets as widgets
from IPython.display import display, clear_output
from datetime import datetime

def agregar_ayuda_contextual(widget, texto_ayuda):
    ayuda_widget = widgets.HTML(value='')

    def mostrar_ayuda(change):
        if change['new']:  # Solo muestra si hay algo escrito
            ayuda_widget.value = f'<span style="color: gray;">{texto_ayuda}</span>'
        else:
            ayuda_widget.value = ''

    widget.observe(mostrar_ayuda, names='value')

    return widgets.VBox([widget, ayuda_widget])



def funcion():
   
    name = widgets.Text(description='Name:', placeholder='Simulation name')
    name_con_ayuda = agregar_ayuda_contextual(name, 'Introduce un nombre para la simulación.')
    
    nameout = widgets.Text(
        description='Name NETCDF file:',
        placeholder='Name for Output',
        tooltip='Nombre del archivo netCDF que se generará'
    )
    carpeta_malla = "topobatimetrias"
    archivos_malla = sorted([f for f in os.listdir(carpeta_malla) if os.path.isfile(os.path.join(carpeta_malla, f))])
    opciones_bath = archivos_malla if archivos_malla else ['(sin archivos)']
    
    Bath_file = widgets.Dropdown(description='Bathymetry file:', options=opciones_bath)

    initstate = widgets.Dropdown(
        options=[
            ('0: Sea surface displacement from file', 0),
            ('1: Standard Okada', 1),
            ('2: Standard Okada from file', 2),
            ('3: Triangular Okada', 3),
            ('4: Triangular Okada from file', 4),
            ('5: Sea floor deformation from file', 5
            ('6: Sea floor dynamic deformation from file', 6),
            ('7: Gaussian', 7)
        ],
        description='Initstate:'
    )
    
    
    initstate_fields_vbox = widgets.VBox()

    carpeta_mallaesf = "topobathymetry initialization"
    archivos_mallaesf = sorted([f for f in os.listdir(carpeta_mallaesf) if os.path.isfile(os.path.join(carpeta_mallaesf, f))])
    opciones_bathesf = archivos_mallaesf if archivos_mallaesf else ['pepe']
    init_file = widgets.Dropdown(description='Init file:',options=opciones_bathesf)
    
    kajiura = widgets.Checkbox(description='Kajiura:',value=False)
    ref_depth = widgets.FloatText(description='Ref. depth (m):')
    num_faults = widgets.IntText(description='Num. faults:')
    faults_info = widgets.Textarea(description = 'Data for every fault')
    faults_info_con_ayuda = agregar_ayuda_contextual(faults_info, 'Use one line per fault, each line must contain the following values separated by an space: Time(sec), Lon-epicenter, Lat-epicenter, Depth-hypocenter, Fault-length(km), Fault-width(km), Strike, Dip, Rake, Slip(m)')
    faults_info_con_ayuda2 = agregar_ayuda_contextual(faults_info, 'For every fault, a line containing:Time(sec), File with the accumulated deformation')
    faults_info_con_ayuda3 = agregar_ayuda_contextual(faults_info, 'Use one line per fault, each line must contain the following values separated by an space: Time(sec), Lon_v1, Lat_v1, Depth_v1(km), Lon_v2, Lat_v2, Depth_v2(km), Lon_v3, Lat_v3, Depth_v3(km), Rake, Slip(m)')
    faults_info_con_ayuda4 = agregar_ayuda_contextual(faults_info, 'A line containing: File with the dynamic deformation')
    
    carpeta_mallafwf = "files with faults"
    archivos_mallafwf = sorted([f for f in os.listdir(carpeta_mallafwf) if os.path.isfile(os.path.join(carpeta_mallafwf, f))])
    opciones_bathfwf = archivos_mallafwf if archivos_mallafwf else ['hola']
    okada_file = widgets.Dropdown(description='Fault file:',options=opciones_bathfwf)
    
    gaussian_params = widgets.Text(description='Gaussian data:')
    gaussian_params_con_ayuda = agregar_ayuda_contextual(gaussian_params, 'Use one line per fault, each line must contain the following values separated by an space: Time(sec), Lon_center, Lat_center, Height(m), Sigma(km)')
    crop_type = widgets.Dropdown(
        options=[('No', 0), ('Relative', 1), ('Absolute', 2), ('Window', 3)],
        description='Cropping:'
    )
    crop_value = widgets.Text(description='Crop param:')
    ayuda_crop = widgets.HTML(value='')  # ayuda vacía al principio
    crop_value_con_ayuda = widgets.VBox([crop_value, ayuda_crop])  # combinado
    
    # Función que actualiza el texto de ayuda según la opción elegida
    def actualizar_ayuda(change):
        value = change["new"]
        if value == 0:
            text = ""
        elif value == 1:
            text = "Percentage [0,1]"
        elif value == 2:
            text = "Threshold value (m)"
        elif value == 3:
            text = "Lon_center Lat_center Radius(km)"
        else:
            text = ""
        ayuda_crop.value = f'<span style="color: gray;">{text}</span>' if text else ''
    
    # Inicializar ayuda
    actualizar_ayuda({"new": crop_type.value})
    
    # Vincular el observador al dropdown
    crop_type.observe(actualizar_ayuda, names="value")
    
    # Campos para cada initstate (ejemplo para 0, 1 y 7)
    def campos_para_estado(valor):
        if valor == 0:
            return [init_file]
        elif valor == 1:
            return [kajiura, ref_depth, num_faults,faults_info_con_ayuda,okada_file, crop_type, crop_value_con_ayuda]
        elif valor == 2:
            return [kajiura, ref_depth,okada_file, crop_type, crop_value_con_ayuda]
        elif valor == 3:
            return [kajiura, ref_depth,num_faults, faults_info_con_ayuda3, crop_type, crop_value_con_ayuda]
        elif valor == 4:
            return [kajiura, ref_depth,okada_file, crop_type, crop_value_con_ayuda]
        elif valor == 5:
            return [kajiura, ref_depth,num_faults,faults_info_con_ayuda2]
        elif valor == 6:
            return [kajiura, ref_depth, faults_info_con_ayuda4]
        elif valor == 7:
            return [num_faults, gaussian_params_con_ayuda, crop_type, crop_value_con_ayuda]
        else:
            return []

    def actualizar_campos(change):
        nuevo_valor = change['new']
        nuevos_campos = campos_para_estado(nuevo_valor)
        initstate_fields_vbox.children = nuevos_campos
    
    initstate.observe(actualizar_campos, names='value')
    
    actualizar_campos({'new': initstate.value})

    
    opciones_netcdf = [
        "eta","maximum eta","velocities","maximum velocities",
        "modulus of velocity","maximum modulus of velocity","maximum modulus of mass flow",
        "momentum_flux","maximum momentum flux","arrival times"
    ]
    checkboxes = {opcion: widgets.Checkbox(value=False, description=opcion) for opcion in opciones_netcdf}
    checkboxes_vbox = widgets.VBox(list(checkboxes.values()))

    num_levels = widgets.BoundedIntText(
        value=1, min=1, max=10, step=1,
        description='Number of levels:',
        style={'description_width': 'initial'}
    )

    # Contenedor para los campos dinámicos
    campos_condicionales_vbox = widgets.VBox([])

    contenedores_submallas = []  

    def actualizar_campos_condicionales(change):
        n = change["new"]
        campos = []
        contenedores_submallas.clear()
    
        if n <= 1:
            campos_condicionales_vbox.children = []
            return
    
        # Opciones de ejemplo para desplegables y checkboxes
        
        carpeta_mallabf = "topobatimetrias"
        archivos_mallabf = sorted([f for f in os.listdir(carpeta_mallabf) if os.path.isfile(os.path.join(carpeta_mallabf, f))])
        opciones_bathbf = archivos_mallabf if archivos_mallabf else ['']
        
        carpeta_mallaesf = "topobathymetry initialization"
        archivos_mallaesf = sorted([f for f in os.listdir(carpeta_mallaesf) if os.path.isfile(os.path.join(carpeta_mallaesf, f))])
        opciones_bathesf = archivos_mallaesf if archivos_mallaesf else ['']
        
        opciones_check = ["eta","maximum eta","velocities","maximum velocities",
        "modulus of velocity","maximum modulus of velocity","maximum modulus of mass flow",
        "momentum_flux","maximum momentum flux","arrival times"]
    
        for i in range(n - 1):
            refinement = widgets.FloatText(
                description=f'Refinement ratio {i+1}:',
                placeholder=f'Refinement ratio for level {i+1}'
            )
            submallas = widgets.IntText(
                description=f'Submallas {i+1}:',
                placeholder=f'Número de submallas nivel {i+1}',
                value=1  # Valor inicial
            )
    
            # Contenedor dinámico donde se insertarán las opciones por cada submalla
            contenedor_submallas_dinamico = widgets.VBox([])
    
            def generar_campos_submalla(change, contenedor=contenedor_submallas_dinamico):
                total = change['new']
                widgets_por_submalla = []
    
                for j in range(total):
                    
                    dropdown1 = widgets.Dropdown(
                        description=f'bathymetry file [{i+1}][{j+1}]:',
                        options=opciones_bathbf
                    )
                    dropdown2 = widgets.Dropdown(
                        description=f'Initial state [{i+1}][{j+1}]:]:',
                        options=opciones_bathesf
                    )
                    salida_nombre = widgets.Text(
                        description=f'Output name [{j+1}]:',
                        placeholder='Nombre de salida'
                    )
                    checkboxes = [widgets.Checkbox(description=label) for label in opciones_check]
                    checks_vbox = widgets.VBox(checkboxes)
    
                    bloque_submalla = widgets.VBox([
                        widgets.HTML(f"<h4>Submalla {j+1} (nivel {i+1}):</h4>"),
                        dropdown1,
                        dropdown2,
                        salida_nombre,
                        checks_vbox
                    ])
                    widgets_por_submalla.append(bloque_submalla)
    
                contenedor.children = widgets_por_submalla
    
            submallas.observe(generar_campos_submalla, names='value')
            generar_campos_submalla({'new': submallas.value})
    
            grupo_nivel = widgets.VBox([
                refinement,
                submallas,
                contenedor_submallas_dinamico
            ])
    
            contenedores_submallas.append(grupo_nivel)
            campos.append(grupo_nivel)

        campos_condicionales_vbox.children = campos


    num_levels.observe(actualizar_campos_condicionales, names='value')
    actualizar_campos_condicionales({'new': num_levels.value})


    Upper_border_condition = widgets.IntText(
        description='Upper border condition:',
        value=0
    )
    Upper_border_condition_con_ayuda = agregar_ayuda_contextual(Upper_border_condition, '(1: open, -1: wall)')
    
    Lower_border_condition = widgets.IntText(
        description='Lower border condition:',
        value=0
    )
    Lower_border_condition_con_ayuda = agregar_ayuda_contextual(Lower_border_condition, '(1: open, -1: wall)')
    
    Left_border_condition = widgets.IntText(
        description='Left border condition:',
        value=0
    )
    Left_border_condition_con_ayuda = agregar_ayuda_contextual(Left_border_condition, '(1: open, -1: wall)')
    
    Right_border_condition = widgets.IntText(
        description='Right border condition:',
        value=0
    )
    Right_border_condition_con_ayuda = agregar_ayuda_contextual(Right_border_condition, '(1: open, -1: wall)')
    
    Simulation_time = widgets.IntText(
        description='Simulation time:',
        value=1
    )
    
    Saving_time_netCDFfiles = widgets.IntText(
        description='Saving time netCDFfiles:',
        value=1
    )
    Saving_time_netCDFfiles_con_ayuda = agregar_ayuda_contextual(Saving_time_netCDFfiles, '(sec) (-1: do not save)')
    
    Read_points = widgets.Checkbox(
        description='Read points:',
        value=False
    )
    Show_advanced_options = widgets.Checkbox(
        description='Show advanced options :',
        value=False
    )

    points_box = widgets.VBox([])

    File_points = widgets.Text(
        description='File points:',
        placeholder='Text files specifying the points for series outputs',
        layout=widgets.Layout(width='50%')
    )
    Saving_time_time_series = widgets.FloatText(description='Saving time time series:', value=20.0)

    def toggle_read_points(change):
        if change['new']:
            points_box.children = [
                widgets.HTML("<h4>Read points:</h4>"),
                File_points,
                Saving_time_time_series
            ]
        else:
            points_box.children = []
    
    Read_points.observe(toggle_read_points, names='value')
    toggle_read_points({'new': Read_points.value})

    
    advanced_numerics_box = widgets.VBox([])
    
    cfl = widgets.FloatText(description='CFL:', value=0.5)
    epsilon_h = widgets.FloatText(description='Epsilon h:', value=0.01)
    threshold_2swaf = widgets.FloatText(description='Threshold 2sWAF:', value=20.0)
    stability_coeffs = widgets.Text(
        description='Stability coeffs:',
        placeholder='e.g. 1.0 0.5 0.25',
        layout=widgets.Layout(width='50%')
    )
    
    def toggle_advanced_options(change):
        if change['new']:
            advanced_numerics_box.children = [
                widgets.HTML("<h4>Advanced Numerical Options:</h4>"),
                cfl,
                epsilon_h,
                threshold_2swaf,
                stability_coeffs
            ]
        else:
            advanced_numerics_box.children = []
    
    Show_advanced_options.observe(toggle_advanced_options, names='value')
    toggle_advanced_options({'new': Show_advanced_options.value})
    
    extra_fields_vbox3 = widgets.VBox([
        widgets.HTML("<h3>Simulation and Saving Parameters:</h3>"),
        Upper_border_condition_con_ayuda,
        Lower_border_condition_con_ayuda,
        Left_border_condition_con_ayuda,
        Right_border_condition_con_ayuda,
        Simulation_time,
        Saving_time_netCDFfiles_con_ayuda,
        Read_points,
        Show_advanced_options 
    ])
    
    friction_type = widgets.BoundedIntText(value=0, min=0, max=2, description='Friction type:')
    friction_type_con_ayuda = agregar_ayuda_contextual(friction_type, '(0: Fixed \n1: Variable friction specifying the frictions of the level 0 mesh \n2: Variable friction specifying the frictions of all the submeshes)')

    bottom_friction_coeff = widgets.FloatText(
        description='Bottom friction coeff:',
        value=0.03
    )

    bottom_friction_file = widgets.Text(
        description='Friction file:',
        placeholder='Enter friction file name'
    )

    max_velocity = widgets.FloatText(
        description='Max water velocity:',
        value=100.0
    )

    typical_length = widgets.FloatText(
        description='Typical length (m):',
        value=100000.0
    )

    typical_depth = widgets.FloatText(
        description='Typical depth (m):',
        value=1000.0
    )

    threshold_arrival = widgets.FloatText(
        description='Threshold arrival:',
        value=0.0
    )

    extra_fields_vbox4 = widgets.VBox([
        widgets.HTML("<h3>Friction and Physical Parameters:</h3>"),
        friction_type_con_ayuda,
        bottom_friction_coeff,
        bottom_friction_file,
        max_velocity,
        typical_length,
        typical_depth,
        threshold_arrival
    ])

    # Botón y salida
    guardar_btn = widgets.Button(description='Guardar', button_style='success', icon='save')
    salida = widgets.Output()

    init_map = {
        'Sea surface displacement from file': 0,
        'Standard Okada': 1,
        'Standard Okada from file': 2,
        'Triangular Okada': 3,
        'Triangular Okada from file': 4,
        'Sea floor deformation from file': 5,
        'Sea floor dynamic deformation from file': 6,
        'Gaussian': 7
    }

    
    def guardar_datos(b):
        with salida:
            clear_output()
            seleccionados = [op for op, cb in checkboxes.items() if cb.value]  # checkboxes GLOBAL
    
            niveles_info = []
    
            for nivel in contenedores_submallas:
                children = nivel.children
                refinement_ratio = children[0].value  # FloatText
                num_submallas = children[1].value     # IntText
                submallas_contenedor = children[2]    # VBox con submallas
            
                submallas_info = []
            
                for submalla in submallas_contenedor.children:
                    # Cada submalla es un VBox con 5 hijos
                    dropdown1 = submalla.children[1].value
                    dropdown2 = submalla.children[2].value
                    output_name = submalla.children[3].value
            
                    checkboxes_submalla = submalla.children[4].children  # Nuevo nombre
                    checks_seleccionados = [cb.description for cb in checkboxes_submalla if cb.value]
            
                    submallas_info.append({
                        "dropdown1": dropdown1,
                        "dropdown2": dropdown2,
                        "output_name": output_name,
                        "checks": checks_seleccionados
                    })
            
                niveles_info.append({
                    "refinement_ratio": refinement_ratio,
                    "num_submallas": num_submallas,
                    "submallas": submallas_info
                })
    
            datos = {
                "name": name.value,
                "bath_file": Bath_file.value,
                "initstate": initstate.value,
                "init_file": init_file.value,
                "okada_file": okada_file.value,
                "gaussian_params": gaussian_params.value,
                "kajiura": kajiura.value,
                "ref_depth": ref_depth.value,
                "num_faults": num_faults.value,
                "faults_info": faults_info.value,
                "crop_type": crop_type.value,
                "crop_value": crop_value.value,
                "nameout": nameout.value,
                "netcdf_outputs": seleccionados,
                "num_levels": num_levels.value,
                "Upper_border_condition": Upper_border_condition.value,
                "Lower_border_condition": Lower_border_condition.value,
                "Left_border_condition": Left_border_condition.value,
                "Right_border_condition": Right_border_condition.value,
                "Simulation_time": Simulation_time.value,
                "Saving_time_netCDFfiles": Saving_time_netCDFfiles.value,
                "File_points": File_points.value,
                "Saving_time_time_series": Saving_time_time_series.value,
                "cfl": cfl.value,
                "epsilon_h": epsilon_h.value,
                "threshold_2swaf": threshold_2swaf.value,
                "stability_coeffs": stability_coeffs.value,
                "friction_type": friction_type.value,
                "bottom_friction_coeff": bottom_friction_coeff.value,
                "bottom_friction_file": bottom_friction_file.value,
                "max_velocity": max_velocity.value,
                "typical_length": typical_length.value,
                "typical_depth": typical_depth.value,
                "threshold_arrival": threshold_arrival.value,
                "fecha": datetime.now().isoformat(),
                "niveles": niveles_info
            }
            
            ruta = "fichero entrada hysea/datos_formateados.txt"
            os.makedirs(os.path.dirname(ruta), exist_ok=True)
            
            with open(ruta, "w", encoding="utf-8") as f:
                # Bathymetry block
                f.write(f"{datos['name']}\n")  # Bathymetry name
                f.write(f"{datos['bath_file']}\n")  # Bathymetry file
            
                init_code = datos["initstate"]
                f.write(f"{init_code}\n")

                if init_code in [0, 5, 6] and not datos["init_file"]:
                    print("⚠️ Falta especificar 'init_file'.")
                    return
                
                if init_code == 0:
                    f.write(f"{datos['init_file']}\n")
                elif init_code == 1:
                    f.write(f"{int(datos['kajiura'])}\n")
                    if int(datos['kajiura']) == 1:
                        f.write(f"{datos['ref_depth']}\n")
                    f.write(f"{datos['num_faults']}\n")
                    f.write(f"{datos['faults_info']}\n")
                    f.write(f"{datos['okada_file']}\n")
                    f.write(f"{datos['crop_type']}\n")
                    if (datos['crop_type'] != 0):
                        f.write(f"{datos['crop_value']}\n")
                elif init_code == 2:
                    f.write(f"{int(datos['kajiura'])}\n")
                    if int(datos['kajiura']) == 1:
                        f.write(f"{datos['ref_depth']}\n")
                    f.write(f"{datos['okada_file']}\n")
                    f.write(f"{datos['crop_type']}\n")
                    if (datos['crop_type'] != 0):
                        f.write(f"{datos['crop_value']}\n")
                elif init_code == 3:
                    f.write(f"{int(datos['kajiura'])}\n")
                    if int(datos['kajiura']) == 1:
                        f.write(f"{datos['ref_depth']}\n")
                    f.write(f"{datos['num_faults']}\n")
                    f.write(f"{datos['faults_info']}\n")
                    f.write(f"{datos['crop_type']}\n")
                    if (datos['crop_type'] != 0):
                        f.write(f"{datos['crop_value']}\n")
                elif init_code == 4:
                    f.write(f"{int(datos['kajiura'])}\n")
                    if int(datos['kajiura']) == 1:
                        f.write(f"{datos['ref_depth']}\n")
                    f.write(f"{datos['okada_file']}\n")
                    f.write(f"{datos['crop_type']}\n")
                    if (datos['crop_type'] != 0):
                        f.write(f"{datos['crop_value']}\n")
                elif init_code == 5:
                    f.write(f"{int(datos['kajiura'])}\n")
                    if int(datos['kajiura']) == 1:
                        f.write(f"{datos['ref_depth']}\n")
                    f.write(f"{datos['num_faults']}\n")
                    f.write(f"{datos['faults_info']}\n")
                elif init_code == 6:
                    f.write(f"{int(datos['kajiura'])}\n")
                    if int(datos['kajiura']) == 1:
                        f.write(f"{datos['ref_depth']}\n")
                    f.write("1\n")
                    f.write(f"{datos['faults_info']}\n")
                elif init_code == 7:
                    f.write(f"{datos['num_faults']}\n")
                    f.write(f"{datos['gaussian_params']}\n")
                    f.write(f"{datos['crop_type']}\n")
                    if (datos['crop_type'] != 0):
                        f.write(f"{datos['crop_value']}\n")
                
                # NetCDF outputs for level 0
                f.write(f"{datos['nameout']}\n")
                
                orden_var = [
                    "eta", "maximum eta", "velocities", "maximum velocities",
                    "modulus of velocity", "maximum modulus of velocity",
                    "maximum modulus of mass flow", "momentum_flux",
                    "maximum momentum flux", "arrival times"
                ]
                linea_netcdf = " ".join(["1" if var in datos["netcdf_outputs"] else "0" for var in orden_var])
                f.write(f"{linea_netcdf}\n")
            
                # Levels
                f.write(f"{datos['num_levels']}\n")
            
                # Submesh information if levels > 1
                for nivel in datos["niveles"]:
                    f.write(f"{nivel['refinement_ratio']}\n")
                    f.write(f"{nivel['num_submallas']}\n")
                    for sub in nivel["submallas"]:
                        f.write(f"{sub['dropdown1']}\n")  # Bathymetry file
                        f.write(f"{sub['dropdown2']}\n")  # Initial state file
                        f.write(f"{sub['output_name']}\n")
                        linea_submesh = " ".join(["1" if v in sub["checks"] else "0" for v in orden_var])
                        f.write(f"{linea_submesh}\n")
            
                # Border & simulation params
                f.write(f"{datos['Upper_border_condition']}\n")
                f.write(f"{datos['Lower_border_condition']}\n")
                f.write(f"{datos['Left_border_condition']}\n")
                f.write(f"{datos['Right_border_condition']}\n")
                f.write(f"{datos['Simulation_time']}\n")
                f.write(f"{datos['Saving_time_netCDFfiles']}\n")
                
                # Read points
                if datos["File_points"]:
                    f.write("1\n")
                    f.write(f"{datos['File_points']}\n")
                    f.write(f"{datos['Saving_time_time_series']}\n")
                else:
                    f.write("0\n")
            
                # Advanced numerics
                f.write(f"{datos['cfl']}\n")
                f.write(f"{datos['epsilon_h']}\n")
                f.write(f"{datos['threshold_2swaf']}\n")
                f.write(f"{datos['stability_coeffs'] or '0.2'}\n")
            
                # Friction
                f.write(f"{datos['friction_type']}\n")
                if datos['friction_type'] == 0:
                    f.write(f"{datos['bottom_friction_coeff']}\n")
                elif datos['friction_type'] == 1:
                    f.write(f"{datos['bottom_friction_file']}\n")
                else:
                    f.write(f"{datos['bottom_friction_file']}\n")  # Placeholder
            
                # Física
                f.write(f"{datos['max_velocity']}\n")
                f.write(f"{datos['typical_length']}\n")
                f.write(f"{datos['typical_depth']}\n")
            
                # Arrival threshold (only if activated)
                if "arrival times" in datos["netcdf_outputs"]:
                    f.write(f"{datos['threshold_arrival']}\n")
            
                
            print("✅ Datos guardados correctamente en:", ruta)


    guardar_btn.on_click(guardar_datos)

    formulario = widgets.VBox([
        widgets.HTML("<h1>Formulario interactivo</h1>"),
        name_con_ayuda,
        Bath_file,
        widgets.HTML("<h3>Configuración de estado inicial</h3>"),
        initstate,
        initstate_fields_vbox,
        widgets.HTML("<h3>NETCDF Outputs:</h3>"),
        nameout,
        checkboxes_vbox,
        num_levels,
        campos_condicionales_vbox,
        extra_fields_vbox3,
        points_box,
        advanced_numerics_box,
        extra_fields_vbox4,  
        guardar_btn,
        salida
    ])
    display(formulario)

    
    
