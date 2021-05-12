import networkx as nx

def printNetwork(G, filename):
    m = nx.linalg.graphmatrix.adjacency_matrix(G)
    with open(filename + ".txt", "w") as f:
        f.write("{}\n".format(G.number_of_nodes()))
        f.write("{}\n".format(G.number_of_edges()))
        for edge in G.edges:
            f.write("{} {}\n".format(edge[0], edge[1]))

    with open(filename + "AM.txt", "w") as f:
        string = ""
        for i in range(m.shape[0]):
            for j in range(m.shape[1]):
                string += str(int(m[i,j]))
            string += "\n"
        f.write(string)

#Barabasi Albert networks
G = nx.generators.random_graphs.barabasi_albert_graph(500,3)
G = nx.Graph(G)
printNetwork(G, "BA_N500_m3")
nx.write_pajek(G, "BA_N500_m3.net")

G = nx.generators.random_graphs.barabasi_albert_graph(2000,2)
G = nx.Graph(G)
printNetwork(G, "BA_N2000_m2")
nx.write_pajek(G, "BA_N2000_m2.net")

#Erdos Renyi networks
G = nx.generators.random_graphs.erdos_renyi_graph(500, 0.012)
G = nx.Graph(G)
printNetwork(G, "ER_N500_p012")
nx.write_pajek(G, "ER_N500_p012.net")

G = nx.generators.random_graphs.erdos_renyi_graph(1000, .002)
G = nx.Graph(G)
printNetwork(G, "ER_N1000_p002")
nx.write_pajek(G, "ER_N1000_p002.net")

#Scale free networks
s = nx.utils.powerlaw_sequence(500, 2.7)
G = nx.expected_degree_graph(s, selfloops=False)
G = nx.Graph(G)
printNetwork(G, "SF_N500_g27")
nx.write_pajek(G, "SF_N500_g27.net")

s = nx.utils.powerlaw_sequence(2000, 2.5)
G = nx.expected_degree_graph(s, selfloops=False)
G = nx.Graph(G)
printNetwork(G, "SF_N2000_g25")
nx.write_pajek(G, "SF_N2000_g25.net")
