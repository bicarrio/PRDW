# -*- coding: cp1252 -*-
##ultima version Julio 2015.
##  redeficion de funciones e interaccion con archivo 'project'
##  limpiza general del archivo.
##  no es necesario definir los triangulos vecinos.
##  no es necesario definir los segmentos de los hoyos. ANUGA lo hace solo.
##  incorpora definicion de subregion, para usar diferente fricción.
##Benjamín Carrion PRDW

# import project
import numpy as np
# from Scientific.IO.NetCDF import NetCDFFile
import netCDF4 as netcdf

from anuga.file.sww import Write_sww
from anuga.config import netcdf_mode_w, netcdf_mode_r, netcdf_float
from anuga.coordinate_transforms.geo_reference import Geo_reference

##DEFINICION DE FUNCIONES ------------------------------------------------------
def lee_vertices_y_triangulos(malla_mike, bathy_anuga = 'bathy.csv'):
    """Lee los vertices y triangulos desde la malla de MIKE.
       Retorna:
               vertices     [['x0', 'y0'], ['x1', 'x2'], ...]
               tags         ['tag0', 'tag1', ...]
               tri_v0       ['v00', 'v01', ...]
               tri_v1       ['v10', 'v11', ...]
               tri_v2       ['v20', 'v21', ...]
    """
    
    fid = open(malla_mike)
    print('Leyendo vertices y triangulos de '+ malla_mike +' ...')
    
    #1. lee los vertices -------------------------------------------------------
    vertices = []
    tags = []
    
    line = fid.readline()
    fields = line.split()
    num_vertices = int(fields[2])
    print('    Numero de vertices: ' + str(num_vertices) + '.')    
    
    fout = open(bathy_anuga, 'w')
    fout.write('easting,northing,elevation\n')

    for i in range(num_vertices):
        line = fid.readline()
        fields = line.split()        
        vertice = [fields[1], fields[2]]
        vertices.append(vertice)
        tag = equivalencia_tag(int(fields[4]))
        tags.append(tag)
        
        #imprime bathy en formato XYZ
        fout.write(fields[1] + ',' + fields[2] + ',' + fields[3] + '\n')
    fout.close()    #cierra archivo XYZ

    #2. lee los triangulos -----------------------------------------------------
    tri_v0 = []
    tri_v1 = []
    tri_v2 = []
        
    line = fid.readline()
    fields = line.split()
    num_triangulos = int(fields[0])
    print('    Numero de triangulos: ' + str(num_triangulos) + '.')
          
    for i in range(num_triangulos):
        line = fid.readline()
        fields = line.split()

        vertice0 = int(fields[1])-1
        vertice1 = int(fields[2])-1
        vertice2 = int(fields[3])-1
        
        tri_v0.append(repr(vertice0))
        tri_v1.append(repr(vertice1))
        tri_v2.append(repr(vertice2))
    fid.close()     #cierra malla mike

    return [vertices, tags, tri_v0, tri_v1, tri_v2]


def lee_poligonos_mike(bounding_polygon, *inner_polygons):
    """
    Lee los polígonos de contorno e interiores (si existen)
    y los trans forma a segmentos de la malla de anuga.         
                                                                
    Los polígonos deben estar en formato MIKE.                  
    El polígono exterior puede tener varios segmentos           
    La conversión de tags para los segmentos es la siguiente:   
        1: interior                                                 
        2: water            
        3: land
        4: subregion
    """

    vertices_malla = []
    segmentos_malla = []
    
    fid = open(bounding_polygon,'r')
    print('    Poligono actual: ' + bounding_polygon + '.')
    line = fid.readline()
    fields = line.split()

    conector1 = int(fields[2])
    vertices_malla.append([fields[0], fields[1]])
    
    tag1 = int(fields[4])
    tag_ini = tag1

    vertice = 0
    while 1:
        line = fid.readline()
        if not line:
            break
       
        fields = line.split()
        conector2 = int(fields[2])
        tag2 = int(fields[4])

        if conector2 == 1:
            vertice = vertice + 1
            vertices_malla.append([fields[0], fields[1]])
            segmentos_malla.append([repr(vertice-1), repr(vertice), equivalencia_tag(min(tag1,tag2))])
            tag1 = tag2
            
        else:
            if vertices_malla[0] == [fields[0], fields[1]]:
                segmentos_malla.append([repr(vertice), repr(0), equivalencia_tag(min(tag_ini,tag2))])
            
    fid.close()

    #esquisita recursión para los otros polígonos
    for poligono in inner_polygons:
        [vertices_dummy, segmentos_dummy] = lee_poligonos_mike(poligono)
        for i in range(len(vertices_dummy)):
            vertices_malla.append(vertices_dummy[i])
            segmentos_malla.append([repr(int(segmentos_dummy[i][0])+vertice+1),
                                    repr(int(segmentos_dummy[i][1])+vertice+1),
                                    segmentos_dummy[i][2]])        
    
    return [vertices_malla, segmentos_malla]


def lee_segmentos_y_tags(tags, tri_v0, tri_v1, tri_v2, rango = ['water', 'land']):
    """segmentos =[#vertex1, #vertex2, tag] """

    print('Leyendo segmentos interiores:...')
    segmentos = []
    num_vertices = len(tags)
    vertices_ok = range(num_vertices)

    #iterar solo para "water" y "land"
    for tag in rango:
        #encuentra el primero y lo marca como usado
        index = findall(tags, tag)     #todos los nodos de ese tag  
        extremos = sub_busca_primero(tag, tags, tri_v0, tri_v1, tri_v2)
        if tag == 'water':
            extremos_water = extremos
        elif tag == 'land':
            extremos_land = extremos
        
        if extremos:
            primero = extremos[0]
        else:
            primero = index[0]

        vertices_ok[primero] = -1
        #subrutina recursiva
        [segmentos, vertices_ok] = sub_segmentos_y_tags(tag, primero, vertices_ok,
                                                        segmentos, tags,
                                                        tri_v0, tri_v1, tri_v2)

    #cierra los segmentos entre 'water' y 'land'
##    for i in range(len(tri_v0)):
##        if 
##
##    for i in range(len(extremos_land)):
##        land = extremos_land[i]
##        tri_l = argwhere(triangulos[:, 0:3] == land)[:,0]
##        for j in range(len(extremos_water)):
##            water = int(extremos_water[j])
##            tri_w = argwhere(triangulos[:, 0:3] == water)[:,0]
##
##            for ii in range(len(tri_l)):
##                for jj in range(len(tri_w)):
##                    if tri_l[ii] == tri_w[jj]:
##                        segmentos.append([repr(land), repr(water), 'land'])
##                        print('agrega segmentos entre tierra y mar')        
                                        
    return segmentos


def prepare_timeboundary(filename):
    """Convert benchmark 2 time series to NetCDF tms file.
    This is a 'throw-away' code taylor made for files like
    'Benchmark_2_input.txt' from the LWRU2 benchmark
    """
    
    import numpy as num
    from Scientific.IO.NetCDF import NetCDFFile
    from anuga.config import netcdf_float

    print 'Creating', filename

    # Read the ascii (.txt) version of this file
    fid = open(filename[:-4] + '.txt')

    # Skip first line
    line = fid.readline()

    # Read remaining lines
    lines = fid.readlines()
    fid.close()


    N = len(lines)
    T = num.zeros(N, num.float)  #Time
    Q = num.zeros(N, num.float)  #Values

    for i, line in enumerate(lines):
        fields = line.split()

        T[i] = float(fields[0])
        Q[i] = float(fields[1])


    # Create tms NetCDF file

    fid = netcdf.Dataset(filename, 'w')
    fid.institution = 'PRDW'
    fid.description = 'Input wave'
    fid.starttime = 0.0
    fid.createDimension('number_of_timesteps', len(T))
    fid.createVariable('time', netcdf_float, ('number_of_timesteps',))
    fid.variables['time'][:] = T

    fid.createVariable('stage', netcdf_float, ('number_of_timesteps',))
    fid.variables['stage'][:] = Q[:]

    fid.createVariable('xmomentum', netcdf_float, ('number_of_timesteps',))
    fid.variables['xmomentum'][:] = 0.0

    fid.createVariable('ymomentum', netcdf_float, ('number_of_timesteps',))
    fid.variables['ymomentum'][:] = 0.0

    fid.close()

#Fuciones auxiliares ----------------------------------------------------------
def sub_segmentos_y_tags(tag, primero, vertices_ok, segmentos, tags, tri_v0, tri_v1, tri_v2):

    print('    Empezando bucle. Tag: ' + tag)
    print('        primero: ' + repr(primero))
    vertex1 = primero
    index = findall(tags, tag)      #todos los nodos de ese tag    

    seguir = True
    while seguir:
        #elige el sigueinte vecino logico.
        vecinos = sub_busca_vecinos(vertex1, tags, tri_v0, tri_v1, tri_v2)
        vecinos_ok = []
        for i in vecinos:
            if vertices_ok[i] != -1:
                vecinos_ok.append(i)
                  
        if len(vecinos_ok) == 1:
            vertex2 = vecinos_ok[0]
            
        elif len(vecinos_ok) == 2:
            for i in range(2):
                if len(sub_busca_vecinos(vecinos_ok[i], tags, tri_v0, tri_v1, tri_v2))== 2:
                    vertex2 = vecinos_ok[i]

        elif len(vecinos_ok) == 0:            
            #se usaron todos. debiera estar "primero" aca.
            if vecinos.count(primero):
                print('        cerrando bucle')
                vertex2 = primero
            else:
                nuevo_primero = []                
                for i in index:
                    if vertices_ok[i] != -1:
                        nuevo_primero = i
                        break
                if nuevo_primero:
                    print('        nuevo bucle')
                    [segmentos, vertices_ok] = sub_segmentos_y_tags(tag, nuevo_primero,
                                                                    vertices_ok, segmentos,
                                                                    tags, tri_v0, tri_v1, tri_v2)
                else:
                    print('        todos los puntos usados')
                    seguir = False                    
        else:
            print('error: demasiados vecinos')
            print('vertice: '+ repr(vertex1))
            break

        if vertex1 != vertex2:
            segmentos.append([repr(vertex1), repr(vertex2), tag])
            vertices_ok[vertex2] = -1
            vertex1 = vertex2

    return [segmentos, vertices_ok]
              

def sub_busca_primero(tag, tags, tri_v0, tri_v1, tri_v2):
    """busca los extremos de un segmento.
       de no existir, entrega []"""

    primero = []
    index = findall(tags, tag)

    for i in index:
        vecinos = sub_busca_vecinos(i, tags, tri_v0, tri_v1, tri_v2)
        if len(vecinos) == 1:
            primero.append(int(i))

    return primero


def sub_busca_vecinos(punto, tags, tri_v0, tri_v1, tri_v2):
    """indica los vertices vecinos a "punto" que tengan el mismo tag
       vecinos = [vecino1, vecino2, ...]
       """
    
    vecinos = []
    tag = tags[punto]
    
    ind0 = findall(tri_v0, repr(punto))
    for i in ind0:
        punto2 = int(tri_v1[i])
        tag2 = tags[punto2]
        if tag2 == tag and not vecinos.count(punto2):
            vecinos.append(punto2)

        punto2 = int(tri_v2[i])
        tag2 = tags[punto2]
        if tag2 == tag and not vecinos.count(punto2):
            vecinos.append(punto2)
           
    ind1 = findall(tri_v1, repr(punto))
    for i in ind1:
        punto2 = int(tri_v0[i])
        tag2 = tags[punto2]
        if tag2 == tag and not vecinos.count(punto2):
            vecinos.append(punto2)

        punto2 = int(tri_v2[i]) 
        tag2 = tags[punto2]
        if tag2 == tag and not vecinos.count(punto2):
            vecinos.append(punto2)
            
    ind2 = findall(tri_v2, repr(punto))
    for i in ind2:
        punto2 = int(tri_v0[i])
        tag2 = tags[punto2]
        if tag2 == tag and not vecinos.count(punto2):
            vecinos.append(punto2)

        punto2 = int(tri_v1[i])
        tag2 = tags[punto2]
        if tag2 == tag and not vecinos.count(punto2):
            vecinos.append(punto2)
    
    return vecinos


def findall(lista, valor):
    """entrega los indices donde "lista[indices] == valor"""
    i = -1
    indeces = []
    try:
        while 1:
            i = lista.index(valor, i+1)
            indeces.append(i)            
    except ValueError:
        pass

    return indeces

    
def equivalencia_tag(numero_mike):
    """cambia de numero a texto"""
    if numero_mike == 3:
        tag_anuga = 'land'
    elif numero_mike == 2:
        tag_anuga = 'water'
    elif numero_mike == 1:
        tag_anuga = 'interior'
    elif numero_mike == 4:
        tag_anuga = 'subregion'
    else:
        tag_anuga = ''

    return tag_anuga


#Fucion que escribe la malla -----------------------------------------------------------
def escribe_malla_anuga(malla_anuga, vertices, tri_v0, tri_v1, tri_v2,
                        segmentos, vertices_malla, segmentos_malla, regiones):
    """
    Escribe archivo de salida.
    Malla en formato anuga .tsh
    """
    
    print('Escribiendo archivos para ANUGA:...')

    num_vertices = len(vertices)
    num_triangulos = len(tri_v0)
    num_segmentos = len(segmentos)
    num_vert_malla = len(vertices_malla)
    num_seg_malla = len(segmentos_malla)
    num_regiones = len(regiones)
    
    fout = open(malla_anuga, 'w')
    fout.write(repr(num_vertices) + ' 0 # <# of verts> <# of vert attributes>, next lines <vertex #> <x> <y> [attributes] ...Triangulation Vertices...' + '\n')

    for vertice in range(num_vertices):
        fout.write(repr(vertice)+ ' ')
        for i in range(2):
            fout.write(vertices[vertice][i] + ' ')
        fout.write('\n')

    fout.write('# attribute column titles ...Triangulation Vertex Titles... \n')
    fout.write(repr(num_triangulos) + " # <# of triangles>, next lines <triangle #> [<vertex #>] [<neigbouring triangle #>] [attribute of region] ...Triangulation Triangles..." + '\n')
    for triangulo in range(num_triangulos):
        fout.write(repr(triangulo) + ' ')
        fout.write(tri_v0[triangulo] + ' ')
        fout.write(tri_v1[triangulo] + ' ')
        fout.write(tri_v2[triangulo] + ' -1 -1 -1\n')

    fout.write(repr(num_segmentos) + ' # <# of segments>, next lines <segment #> <vertex #>  <vertex #> [boundary tag] ...Triangulation Segments...' + '\n')
    for seg in range(num_segmentos):
        fout.write(repr(seg) + ' ')
        for i in range(2):
            fout.write(segmentos[seg][i] + ' ')
        fout.write(segmentos[seg][2] + '\n')

    fout.write(repr(num_vert_malla) + ' 0 # <# of verts> <# of vert attributes>, next lines <vertex #> <x> <y> [attributes] ...Mesh Vertices...' + '\n')
    for vert in range(num_vert_malla):
        fout.write(repr(vert) + ' ')
        for i in range(2):
            fout.write(vertices_malla[vert][i] + ' ')
        fout.write('\n')

    fout.write(repr(num_seg_malla) + ' # <# of segments>, next lines <segment #> <vertex #>  <vertex #> [boundary tag] ...Mesh Segments...' + '\n')
    for seg in range(num_seg_malla):
        fout.write(repr(seg) + ' ')
        for i in range(2):
            fout.write(segmentos_malla[seg][i] + ' ')
        fout.write(segmentos_malla[seg][2] + '\n')

    fout.write('0 # <# of holes>, next lines <Hole #> <x> <y> ...Mesh Holes...\n') 
    fout.write(repr(num_regiones) + ' # <# of regions>, next lines <Region #> <x> <y> <tag>...Mesh Regions...\n')
    for reg in range(num_regiones):
        fout.write(repr(reg) + ' ')
        for i in range(3):
            fout.write(regiones[reg][i] + ' ')
        fout.write('\n')
##    fout.write('0 0.00001 0.00001\n')
    fout.write(repr(num_regiones) + ' # <# of regions>, next lines <Region #> [Max Area] ...Mesh Regions...\n')
    for reg in range(num_regiones):
        fout.write(repr(reg) + ' ')
        fout.write(regiones[reg][3] + '\n')
##    fout.write('0 0.1\n')
    fout.write('#geo reference\n')
    fout.write('-1\n')
    fout.write('0.0\n')
    fout.write('0.0\n')        
        
    fout.close()
    print('Terminado.')
    #-------------------------------------------------------------

def create_SWW_input(lista_puntos = 'lista_puntos.txt', time_series = 'time_series.txt',
                     sww_file_out = 'boundary_fixed.sww', num_params = 3):
                     
    tipo_datos = np.float32()
    #-------------------------------------------------
    #Lee archivos con posiciones y series de tiempo
    header = 3
    fid = open(lista_puntos,'r')
    line = fid.readline()
    fields = line.split()
    nx = int(fields[1])
    line = fid.readline()
    fields = line.split()
    ny = int(fields[1])
    line = fid.readline()
    fields = line.split()
    dt = float(fields[1])
    fid.close()

    points_utm = np.loadtxt(lista_puntos, dtype = tipo_datos, skiprows = header)[:,0:2]
    number_of_points = len(points_utm)
    elevation = np.loadtxt(lista_puntos, dtype = tipo_datos, skiprows = header)[:,2]

    #lee el tamaño del tiempo, y lo hace relativo a 0
    dummy = np.loadtxt(time_series, dtype = tipo_datos, skiprows = header, delimiter = '\t',
                       usecols = (1,2))
    number_of_times = len(dummy)
    times = np.arange(0,number_of_times*dt,dt)

    #determina número de datos válidos
    fid = open(time_series,'r')
    line = fid.readline()
    line = fid.readline()
    fields = line.split()
    num_datos = (len(fields)-1)/(3*num_params)
    fid.close()

    index = []
    for i in range(1, num_datos*3, 3):
        index.append(int(fields[i][:-1]))

    index_tupla = tuple(np.array(index)-1)

    stage_aux = np.loadtxt(time_series, dtype = tipo_datos, skiprows = header, delimiter = '\t',
                       usecols = tuple(range(1,num_datos+1)))

    xmom_aux = np.loadtxt(time_series, dtype = tipo_datos, skiprows = header, delimiter = '\t',
                      usecols = tuple(range(num_datos+1,2*num_datos+1)))

    ymom_aux = np.loadtxt(time_series, dtype = tipo_datos, skiprows = header, delimiter = '\t',
                      usecols = tuple(range(2*num_datos+1,3*num_datos+1)))

    #crea
    stage = np.zeros((number_of_times,number_of_points),dtype = tipo_datos)
    xmom = np.zeros((number_of_times,number_of_points),dtype = tipo_datos)
    ymom = np.zeros((number_of_times,number_of_points),dtype = tipo_datos)

    stage[:,index_tupla] = stage_aux
    xmom[:,index_tupla] = xmom_aux
    ymom[:,index_tupla] = ymom_aux


    #triangulacion a mano
    vertices = np.zeros((nx, ny))
    i = 0
    for k in range(ny):
        for l in range(nx):
            vertices[l, k] = i
            i += 1

    volumes = []
    for l in range(nx-1):
        for k in range(ny-1):
            v1 = vertices[l, k+1]
            v2 = vertices[l, k]
            v3 = vertices[l+1, k+1]
            v4 = vertices[l+1, k]

            volumes.append([v1, v2, v3])
            volumes.append([v4, v3, v2])

    volumes = np.array(volumes, np.int)

    ##Escribe archivo SWW
    # outfile = NetCDFFile(sww_file_out, netcdf_mode_w)
    outfile = netcdf.Dataset(sww_file_out, netcdf_mode_w)
    new_origin = Geo_reference(-1, 0, 0)
    number_of_volumes = len(volumes)

    sww = Write_sww(['elevation'], ['stage', 'xmomentum', 'ymomentum'])

    sww.store_header(outfile, times, number_of_volumes, number_of_points,
                     description='SWW creado en PRDW',
                     sww_precision = netcdf_float,
                     verbose=True)

    sww.store_triangulation(outfile, points_utm, volumes,
                            new_origin = new_origin,
                            verbose = True)       

    sww.store_static_quantities(outfile, elevation = elevation, verbose = True)

    for i in range(len(times)):
        sww.store_quantities(outfile, sww_precision = netcdf_float, slice_index = i,
                             verbose = True,
                             stage = stage[i,:],
                             xmomentum = xmom[i,:],
                             ymomentum = ymom[i,:])
        
    outfile.close()

        


                     
                     
                     
                     