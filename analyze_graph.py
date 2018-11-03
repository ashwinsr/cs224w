import numpy as np
import matplotlib.pyplot as plt

GRAPH_WEIGHT_MATRIX_FILE = 'weights.csv'
GRAPH_NODE_LABELS_FILE = 'labels.csv'

# The weight matrix specifies the edge weight in a graph
# with M_{i,j} representing the weight of the edge from
# node i to node j
def load_weight_matrix():
    return np.loadtxt(GRAPH_WEIGHT_MATRIX_FILE, delimiter=',')

def plot_hist(values, title, xlabel, output_filename):
    plt.clf()
    plt.hist(values, bins=20)

    plt.xlabel(xlabel)
    plt.ylabel('Count')
    plt.title(title)

    plt.savefig('figures/' + output_filename, format='png')

def heatmap_weights(M):
    plt.clf()
    plt.imshow(M, cmap='hot')

    plt.xlabel('To Node')
    plt.ylabel('From Node')
    plt.title('Heatmap of edge weights')

    plt.savefig('figures/heatmap_weights.png', format='png')

def compute_degree_distributions(M):
    # Compute in degree distribution
    in_degrees = np.sum(M, axis=0)
    plot_hist(in_degrees, 'In Degree Distribution', 'In Degree', 'in_degrees.png')

    # Print out degree distribution
    out_degrees = np.sum(M, axis=1).T
    plot_hist(out_degrees, 'Out Degree Distribution', 'Out Degree', 'out_degrees.png')

M = load_weight_matrix()
heatmap_weights(M)
compute_degree_distributions(M)
