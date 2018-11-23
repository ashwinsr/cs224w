import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

GRAPH_WEIGHT_MATRIX_FILE = 'cds/WeeklySTD_CDS_Network_GVD_p8_keep2500.csv'

EPSILON = 1e-8
EDGE_WEIGHT_THRESHOLD = 0.06

def load_country_list():
    with open(GRAPH_WEIGHT_MATRIX_FILE) as f:
        countries = f.readline().strip().replace('"','').replace('.', ' ').split(',')
    return countries

# The weight matrix specifies the edge weight in a graph
# with M_{i,j} representing the weight of the edge to
# node i from node j.
# Rows sum to 100 in original file and 1 in returned matrix
def load_weight_matrix():
    graph = np.loadtxt(GRAPH_WEIGHT_MATRIX_FILE, delimiter=',', skiprows=1)
    return graph/100

# Plot the high weight edges only
def plot_network(M):
    M = M * (M > EDGE_WEIGHT_THRESHOLD)

    G = nx.Graph()
    for i in range(len(M)):
        for j in range(len(M)):
            if M[i][j] > 0:
                G.add_edge(countries[i], countries[j], weight=M[i][j])

    pos = nx.spring_layout(G)
    nx.draw_networkx_nodes(G, pos=pos)
    nx.draw_networkx_edges(G, pos=pos)
    nx.draw_networkx_labels(G, pos=pos, font_size=9, font_color='b')

    #plt.show()

def plot_hist(values, title, xlabel, output_filename):
    plt.clf()
    plt.hist(values, bins=20)
    plt.gcf().subplots_adjust(bottom=0.15)

    plt.xlabel(xlabel)
    plt.ylabel('Count')
    plt.title(title)

    plt.savefig('figures/' + output_filename, format='png')

def plot_bar(x, y, xlabel, ylabel, title, output_filename):
    # Sort ascending
    y, x = zip(*sorted(zip(y,x)))

    plt.clf()
    
    n = range(len(x))
    plt.bar(n, y, align='center')
    plt.xticks(n, x, rotation='vertical')
    plt.gcf().subplots_adjust(bottom=0.25)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)

    plt.savefig('figures/' + output_filename, format='png')

def heatmap_weights(M):
    plt.clf()
    plt.imshow(M, cmap='hot')

    plt.xlabel('From Node')
    plt.ylabel('To Node')
    plt.title('Heatmap of edge weights')

    plt.gca().set_xticks(range(0,30,1))
    plt.gca().set_yticks(range(0,30,1))

    plt.gca().set_xticklabels(countries, rotation='vertical')
    plt.gca().set_yticklabels(countries)
    plt.gcf().subplots_adjust(bottom=0.25)

    #plt.show()

def compute_degree_distributions(M):
    # Compute in degree distribution
    in_degrees = np.sum(M, axis=1).T
    plot_hist(in_degrees, 'In Degree Distribution', 'In Degree', 'in_degrees_dist.png')
    plot_bar(countries, in_degrees, 'Country', 'In Degree', 'In Degree', 'in_degrees.png')

    # Compute out degree distribution
    out_degrees = np.sum(M, axis=0).T
    plot_hist(out_degrees, 'Out Degree Distribution', 'Out Degree', 'out_degrees_dist.png')
    plot_bar(countries, out_degrees, 'Country', 'Out Degree', 'Out Degree', 'out_degrees.png')

def validate_row_sum(M):
    row_sum_deviation = np.sum(M, axis=1) - 1

    if np.max(np.absolute(row_sum_deviation)) >= EPSILON:
        raise Exception("Row sum of transition matrix is not 1")

def compute_pagerank_eigenvector(M):
    w, v = np.linalg.eig(M)
    pagerank = v[:,np.argmax(w)].T

    plot_bar(countries, pagerank, 'Country', 'Pagerank Score', 'Pagerank Scores', 'pagerank.png')

def compute_pagerank_power(M):
    n = M.shape[0]
    ranks = np.matrix([1.0/n] * n).T
    
    # Power iterate
    while True:
        ranks_new = M.dot(ranks)

        delta = np.max(np.absolute(ranks_new - ranks))
        if delta < EPSILON:
            break
        #print delta

        ranks = ranks_new

    # Normalize and return
    ranks = ranks / np.linalg.norm(ranks)

    plot_bar(countries, ranks, 'Country', 'Pagerank Score', 'Pagerank Scores', 'pagerank.png')
    return ranks

countries = load_country_list()
M = load_weight_matrix()
np.fill_diagonal(M, 0.0)
plot_network(M)

heatmap_weights(M)
compute_degree_distributions(M)
compute_pagerank_power(M) #is this correct? Is directionality of edges correct?
#compute_pagerank_eigenvector(M)

# Use network deconvolution to get rid of second and third order relations. Look at lecture 4
# Compute motif intensity and coherence
# Node centrality
# Clustering coefficient
# The stochastic block model thing
# Role detection with recursive features (lecture 5)
# Community detection with Louvain algorithm (lecture 6)
