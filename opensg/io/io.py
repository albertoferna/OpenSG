"""Input/output utilities for OpenSG.

This module provides functions for reading and writing mesh data in various formats,
including YAML and GMSH formats.
"""

import yaml
import numpy as np
import pyvista as pv


def load_yaml(yaml_file):
    """Load mesh data from a YAML file.

    Parameters
    ----------
    yaml_file : str
        Path to the YAML file containing mesh data

    Returns
    -------
    dict
        Dictionary containing the mesh data with keys:
        - nodes: list of node coordinates
        - elements: list of element definitions
        - sets: dictionary of element and node sets
        - materials: dictionary of material definitions
        - sections: dictionary of section definitions
        - elementOrientations: list of element orientation matrices
    """
    with open(yaml_file, "r") as file:
        mesh_data = yaml.load(file, Loader=yaml.CLoader)
    return mesh_data


def write_yaml(data, yaml_file):
    """Write data to a YAML file.

    Parameters
    ----------
    data : dict
        Data to write to the file
    yaml_file : str
        Path to the output YAML file

    Returns
    -------
    None
    """
    with open(yaml_file, "w") as file:
        yaml.dump(data, file, Loader=yaml.CLoader)
    return


# TODO write a function to validate mesh data schema
def validate_mesh_data(mesh_data):
    """Validate the structure and content of mesh data.

    This function is not yet implemented. It will check that the mesh data
    contains all required fields and that they have the correct structure.

    Parameters
    ----------
    mesh_data : dict
        Dictionary containing mesh data to validate

    Returns
    -------
    bool
        True if the data is valid, False otherwise

    Raises
    ------
    NotImplementedError
        This function is not yet implemented
    """
    raise NotImplementedError("Mesh data validation is not yet implemented")

def write_mesh(filename, mesh_data, format='VTK'):
    """
    Writes mesh data to a file in the specified format.
    Parameters
    ----------
    mesh_data : dict
        Data to write to the file
    format : str
       Format of the output file (default: VTK)
    Returns
    -------
    None
    """
    # check that nodes is a list of float triplets and not a string
    old_format = False
    if isinstance(mesh_data['nodes'][0][0], str):
        # this needs preprocessing. Assuming old data format
        points = []
        old_format = True
        for node in mesh_data['nodes']:
            points.append([float(_) for _ in node[0].split()])
    else:
        points = np.array(mesh_data['nodes'])
    # we need to know the type of element. 4 noded can be quad or tetra.
    # This is slow and it could be vectorize, but let's check all elements:
    cell_types = []
    cells = []
    for element in mesh_data['elements']:
        # check that is a list of integers and not a string
        if old_format:
            # assume it is in the old format and it starts at 1
            elem = [(int(_) - 1) for _ in element[0].split()]
        else:
            elem = element
        if len(elem) == 4:
            # check normals
            p1, p2, p3, p4 = points[elem]
            v1 = p2 - p1
            v2 = p4 - p1
            v3 = p4 - p3
            v4 = p2 - p3
            n1 = np.cross(v1, v2)
            n2 = np.cross(v3, v4)
            check = np.dot(n1/np.linalg.norm(n1), n2/np.linalg.norm(n2))
            if abs(1-check) < 1e-3:
                cell_types.append(pv.CellType.QUAD)
                cells.append([4, *elem])
            else:
                cell_types.append(pv.CellType.TETRA)
                cells.append([4, *elem])
        elif len(elem) == 3:
            cell_types.append(pv.CellType.TRIANGLE)
            cells.append([3, *elem])
        elif len(elem) == 8:
            cell_types.append(pv.CellType.HEXAHEDRON)
            cells.append([8, *elem])
        elif len(elem) == 6:
            cell_types.append(pv.CellType.WEDGE)
            cells.append([6, *elem])
        elif len(elem) == 5:
            cell_types.append(pv.CellType.PYRAMID)
            cells.append([5, *elem])
        else:
            raise ValueError(f"Invalid number of vertices in element {elem}")
    # check if old format was improperly used
    if np.min(cells) < 0:
        raise ValueError('Old format assumed but new format mesh (Node IDs start at 0)')
    grid = pv.UnstructuredGrid(cells, cell_types, points)
    grid['elementOrientations'] = mesh_data['elementOrientations']
    grid['elem_set'] = np.zeros(grid.n_cells, dtype=int)
    grid.cell_data['elem_set'] = 99999
    for i, element_set in enumerate(mesh_data['sets']['element']):
        if old_format:
            labels = np.array(element_set['labels']) - 1
        else:
            labels = np.array(element_set['labels'])
        grid['elem_set'][labels] = i
    grid.save(f'{filename}.vtu')
    return grid


# def write_mesh(filename, blade_mesh):
#     """Write mesh data to a GMSH format file.

#     Parameters
#     ----------
#     filename : str
#         Path to the output mesh file
#     blade_mesh : BladeMesh
#         BladeMesh object containing the mesh data to write

#     Returns
#     -------
#     None

#     Notes
#     -----
#     The mesh is written in GMSH 2.2 format with the following sections:
#     - MeshFormat: version and type information
#     - Nodes: node coordinates
#     - Elements: element definitions with tags
#     """
#     mesh_file = open(filename, 'w')

#     mesh_file.write('$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n')
#     newNumNds = np.max(ndNewLabs)
#     mesh_file.write(str(newNumNds) + '\n')

#     for i, nd in enumerate(nodes):
#         lab = ndNewLabs[i]
#         if(lab > -1):
#             ln = [str(lab),str(nd[2]),str(nd[0]),str(nd[1])]
#         #  ln = [str(lab),str(nd[0]),str(nd[1]),str(nd[2])]
#             mesh_file.write(' '.join(ln) + '\n')

#     mesh_file.write('$EndNodes\n$Elements\n')

#     newNumEls = np.max(elNewLabs)
#     mesh_file.write(str(newNumEls) + '\n')

#     for i, el in enumerate(elements):
#         lab = elNewLabs[i]
#         if(lab > -1):
#             ln = [str(lab)]
#             if(el[3] == -1):
#                 ln.append('2')
#             else:
#                 ln.append('3')
#             ln.append('2')
#             ln.append(str(elLayID[i]+1))
#             ln.append(str(elLayID[i]+1))
#             for nd in el:
#                 if(nd > -1):
#                     ln.append(str(ndNewLabs[nd]))
#             mesh_file.write(' '.join(ln) + '\n')
#     mesh_file.write('$EndElements\n')

#     mesh_file.close()
