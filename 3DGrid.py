import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np

def read_nodes(filename):
    with open(filename) as f:
        nodes = [tuple(map(float, line.strip().split())) for line in f]
    return nodes

def read_elements(filename, zero_indexed=True):
    with open(filename) as f:
        elements = [list(map(int, line.strip().split())) for line in f]
        if not zero_indexed:
            elements = [[i - 1 for i in elem] for elem in elements]
    return elements

def plot_3d_mesh(nodes, elements):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Определение граней гексаэдра
    faces = [
        [0, 1, 3, 2], [4, 5, 7, 6],
        [0, 1, 5, 4], [2, 3, 7, 6],
        [0, 2, 6, 4], [1, 3, 7, 5]
    ]

    for elem in elements:
        verts = [nodes[i] for i in elem]
        for face in faces:
            square = [verts[idx] for idx in face]
            poly = Poly3DCollection(
                [square], 
                alpha=0.3, 
                edgecolor='k',
                facecolor='pink'  
            )
            ax.add_collection3d(poly)

    # Отрисовка узлов
    nodes_np = np.array(nodes)
    ax.scatter(nodes_np[:, 0], nodes_np[:, 1], nodes_np[:, 2], color='r', s=5)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title("3D Сетка конечных элементов")
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    nodes = read_nodes("processFiles/nodes_out.txt")
    elements = read_elements("processFiles/elements_out.txt", zero_indexed=True)
    plot_3d_mesh(nodes, elements)