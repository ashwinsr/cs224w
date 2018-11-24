import numpy as np
import networkx as nx
from networkx.algorithms.community.kernighan_lin import kernighan_lin_bisection
import matplotlib.pyplot as plt

GRAPH_WEIGHT_MATRIX_FILE = 'cds/WeeklySTD_CDS_Network_GVD_p8_keep2500.csv'
EPSILON = 1e-8

# Edge weight threshold only used for
# plotting of high weight edges
EDGE_WEIGHT_THRESHOLD = 0.06

# Mapping from countries to continents
CONTINENT_MAP = {'Brazil': 'South America', 'Chile': 'South America', 'China': 'Asia', 'Colombia': 'South America', 'Indonesia': 'Asia', 'Malaysia': 'Asia', 'Mexico': 'North America', 'Panama': 'North America', 'Peru': 'South America', 'Philippines': 'Asia', 'Russia': 'Asia', 'South Africa': 'Africa', 'South Korea': 'Asia', 'Thailand': 'Asia', 'Turkey': 'Asia', 'Argentina': 'South America', 'Austria': 'Europe', 'Belgium': 'Europe', 'Bulgaria': 'Europe', 'Croatia': 'Europe', 'France': 'Europe', 'Germany': 'Europe', 'Hungary': 'Europe', 'Italy': 'Europe', 'Poland': 'Europe', 'Portugal': 'Europe', 'Romania': 'Europe', 'Slovakia': 'Europe', 'Spain': 'Europe', 'Venezuela': 'South America'}

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

def get_nx_network(M):
    G = nx.DiGraph()
    for i in range(len(M)):
        for j in range(len(M)):
            if M[i][j] > 0:
                G.add_edge(countries[j], countries[i], weight=M[i][j])
    return G

def get_undirected_nx_network(M):
    G = nx.Graph()
    for i in range(len(M)):
        for j in range(i, len(M)):
            if M[i][j] > 0:
                G.add_edge(countries[j], countries[i], weight=(M[i][j]+M[j][i])/2.0)
    return G

# Plot the high weight edges only
def plot_network(M):
    M = M * (M > EDGE_WEIGHT_THRESHOLD)
    G = get_nx_network(M)

    pos = nx.spring_layout(G)
    nx.draw_networkx_nodes(G, pos=pos)
    nx.draw_networkx_edges(G, pos=pos, arrows=True)
    nx.draw_networkx_labels(G, pos=pos, font_size=9, font_color='b')

    #plt.show()

    return G

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

        ranks = ranks_new

    # Normalize and return
    ranks = ranks / np.linalg.norm(ranks)

    plot_bar(countries, ranks, 'Country', 'Pagerank Score', 'Pagerank Scores', 'pagerank.png')
    return ranks

def output_weighted_edge_graph(nx_graph):
    for edge in nx_graph.edges(data=True):
        print edge[0].replace(' ','-'), edge[1].replace(' ','-'), edge[2]['weight']

def convert_community_to_continent(community):
    return sorted([CONTINENT_MAP[c] for c in community])

def get_initial_partition(size):
    # Split the countries array into two partitions,
    # one of which is of size 'size'
    return (countries[:size], countries[size:])

def detect_communities(M):
    # Uses the Kernighan Lin bisection for
    # community detection
    nx_graph = get_undirected_nx_network(M)

    # Try out partitions of all sizes
    for size in range(1, len(M)):
        partition = get_initial_partition(size)
        partition = kernighan_lin_bisection(nx_graph, partition=partition, max_iter=100, weight='weight')
        print "------"
        print "Size:", size
        print "Size a:", len(partition[0])
        print "Size b:", len(partition[1])
        print sorted(partition[0])
        print sorted(partition[1])
        print convert_community_to_continent(partition[0])
        print convert_community_to_continent(partition[1])

def spectral_detect_communities(M):
    # Use spectral clustering to detect communities
    w, v = np.linalg.eig(M)

    # Use the second largest eigenvector to detect communities
    community_split = v[:,1]
    assert all(np.isreal(community_split))
    community_split = community_split.real

    community_a = []
    community_b = []
    for i in range(len(community_split)):
        if community_split[i] > 0:
            community_a.append(countries[i])
        else:
            community_b.append(countries[i])

    print "Size a:", len(community_a)
    print "Size b:", len(community_b)
    print sorted(community_a)
    print sorted(community_b)
    print convert_community_to_continent(community_a)
    print convert_community_to_continent(community_b)

########## BEGIN ANALYSIS ##########

countries = load_country_list()
M = load_weight_matrix()

# We don't care about self influence
np.fill_diagonal(M, 0.0)

# Visualize basic network properties
plot_network(M)
heatmap_weights(M)
compute_degree_distributions(M)

# The highest page rank scores will be nodes that
# are INFLUENCED a lot
compute_pagerank_power(M)

# Find some communities
#detect_communities(M)
spectral_detect_communities(M)

# Make a random graph null model
# Compute motif intensity and coherence

# Role extraction

########## END ANALYSIS ##########
