import rdkit
import networkx as nx

def graph_from_mol(mol):
    """
    Create a networkx.graph based on a molecule structure.

    :param rdkit.mol mol: Structure

    :rtype: networkx.graph
    :return: A graph based on the structure, where atoms are graph nodes,
        and bonds are added as edges between nodes.
    """

    graph = nx.Graph()
    for atom in mol.GetAtoms():
        graph.add_node(atom.GetIdx(), atomic_num=atom.GetAtomicNum())

    for bond in mol.GetBonds():
        graph.add_edge(bond.GetBeginAtomIdx(),
                       bond.GetEndAtomIdx(),
                       weight=bond.GetBondTypeAsDouble())

    return graph
