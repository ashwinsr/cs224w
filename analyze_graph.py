import numpy as np
import matplotlib.pyplot as plt

GRAPH_WEIGHT_MATRIX_FILE = 'weights.csv'
GRAPH_NODE_LABELS_FILE = 'labels.csv'

EPSILON = 1e-8

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

# Compute the pagerank vector given a matrix M
# Weights in M are intended to represent transition
# probabilities such that M_{i, j} indicates weight
# from node i to node j (opposite of typical page rank)
#
# This implies that rows sum to 1.
def validate_row_sum(M):
    row_sum_deviation = np.sum(M, axis=1) - 1

    if np.max(row_sum_deviation) >= EPSILON:
        raise Exception("Row sum of transition matrix is not 1")

def compute_pagerank_eigenvector(M):
    validate_row_sum(M)

    w, v = np.linalg.eig(M.T)
    return v[:,np.argmax(w)].T

def compute_pagerank_power(M):
    validate_row_sum(M)

    n = M.shape[0]
    ranks = np.matrix([1.0/n] * n).T
    
    # Power iterate
    while True:
        ranks_new = M.T.dot(ranks)

        delta = np.max(np.absolute(ranks_new - ranks))
        if delta < EPSILON:
            break

        ranks = ranks_new

    # Normalize and return
    ranks = ranks / np.linalg.norm(ranks)
    return ranks.T

M = load_weight_matrix()
heatmap_weights(M)
compute_degree_distributions(M)
print compute_pagerank_eigenvector(M)
print compute_pagerank_power(M)
