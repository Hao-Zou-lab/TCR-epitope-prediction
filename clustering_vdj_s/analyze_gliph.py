import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from collections import Counter
from Levenshtein import distance

def build_gliph_graph(df_input, netw, label_known='known', label_unknown='unknown', return_full=False):
    df_input['nodes'] = df_input['CDR3b'] + '-' + df_input['source']
    cdr_known = df_input[df_input['source'] == label_known]['CDR3b'].values.tolist()
    cdr_unknown = df_input[df_input['source'] == label_unknown]['CDR3b'].values.tolist()
    cdr_overlap = list(set(cdr_known) & set(cdr_unknown))
    #read network (cdr connections) data
    node1 = [l.split()[0] for l in netw]
    node2 = [l.split()[1] for l in netw]
    link_type = [l.split()[2] for l in netw]
    #rename node names, add connection type
    netw_new = []
    for ndx, node in enumerate(node1):
        if node in cdr_overlap and node2[ndx] not in cdr_overlap:
            if node2[ndx] in cdr_known:
                netw_new.append([node + '-' + label_unknown, node2[ndx] + '-' + label_known, link_type[ndx]])
                netw_new.append([node + '-' + label_known, node2[ndx] + '-' + label_known, link_type[ndx]])
            else:
                netw_new.append([node + '-' + label_known, node2[ndx] + '-' + label_unknown, link_type[ndx]])
                netw_new.append([node + '-' + label_unknown, node2[ndx] + '-' + label_unknown, link_type[ndx]])
        elif node2[ndx] in cdr_overlap and node not in cdr_overlap:
            if node in cdr_known:
                netw_new.append([node + '-' + label_known, node2[ndx] + '-' + label_unknown, link_type[ndx]])
                netw_new.append([node + '-' + label_known, node2[ndx] + '-' + label_known, link_type[ndx]])
            else:
                netw_new.append([node + '-' + label_unknown, node2[ndx] + '-' + label_known, link_type[ndx]])
                netw_new.append([node + '-' + label_unknown, node2[ndx] + '-' + label_unknown, link_type[ndx]])
        elif node in cdr_overlap and node2[ndx] in cdr_overlap:
            if node == node2[ndx]:
                netw_new.append([node + '-' + label_known, node2[ndx] + '-' + label_unknown, link_type[ndx]])
            else:
                netw_new.append([node + '-' + label_known, node2[ndx] + '-' + label_unknown, link_type[ndx]])
                netw_new.append([node + '-' + label_unknown, node2[ndx] + '-' + label_known, link_type[ndx]])
        else:
            if node in cdr_known:
                if node2[ndx] in cdr_known:
                    netw_new.append([node + '-' + label_known, node2[ndx] + '-' + label_known, link_type[ndx]])
                else:
                    netw_new.append([node + '-' + label_known, node2[ndx] + '-' + label_unknown, link_type[ndx]])
            else:
                if node2[ndx] in cdr_known:
                    netw_new.append([node + '-' + label_unknown, node2[ndx] + '-' + label_known, link_type[ndx]])
                else:
                    netw_new.append([node + '-' + label_unknown, node2[ndx] + '-' + label_unknown, link_type[ndx]])
    for cdr in cdr_overlap:
        netw_new.append([cdr + '-' + label_known, cdr + '-' + label_unknown, 'global'])
    #generate networkx graph from edges
    df_netw = pd.DataFrame(netw_new, columns=['node1', 'node2', 'type'])
    df_netw['label'] = ['hetero' if df_netw['node1'][i].split('-')[1] != df_netw['node2'][i].split('-')[1] else 'homo' for i in df_netw.index]
    cdr_graph = nx.from_pandas_edgelist(df_netw, 'node1', 'node2', edge_attr=['type', 'label'])
    ##assign epitope information for epitope-known cdr nodes
    #for n in cdr_graph.nodes:
    #    if 'unknown' not in n:
    #        cdr = n.split('-')[0]
    #        epitope = df_input[(df_input['CDR3b'] == cdr) & (df_input['source'] == label_known)]['epitope'].values[0]
    #        cdr_graph.add_node(n, epitope=epitope)
    df_netw_hetero = df_netw[df_netw['label'] == 'hetero']
    cdr_graph_trimed = nx.from_pandas_edgelist(df_netw_hetero, 'node1', 'node2', edge_attr=['type', 'label'])
    #assign epitope information for epitope-known cdr nodes
    for n in cdr_graph_trimed.nodes:
        if f'-{label_known}' in n:
            cdr = n.split('-')[0]
            epitope = df_input[(df_input['CDR3b'] == cdr) & (df_input['source'] == label_known)]['epitope'].values[0]
            cdr_graph_trimed.add_node(n, epitope=epitope)
    #predict epitope for epitope-unknown cdrs based on epi-known cdrs in the same cluster
    for n in cdr_graph_trimed.nodes:
        candidates = []
        if f'-{label_unknown}' not in n:
            continue
        for nbr, attr in cdr_graph_trimed.adj[n].items():
            if f'-{label_known}' in nbr:
                epitope = cdr_graph_trimed.nodes[nbr]['epitope']
                cdr_dist = distance(n.split('-')[0], nbr.split('-')[0])
                link_type = 0 if attr['type'] == 'local' else 1 if attr['type'] == 'global' else 2
                candidates.append((epitope, cdr_dist, link_type))
        candidates = sorted(candidates, key=lambda x: (x[2], x[1]))
        cdr_graph_trimed.add_node(n, epitope=candidates[0][0])
    if return_full:
        return cdr_graph, cdr_graph_trimed
    else:
        return cdr_graph_trimed


if __name__ == "__main__":
    gliph_input = 'tcr_gliph_input.tsv'
    network = gliph_input + '-clone-network.txt'
    label_unknown = 'Cell'
    label_known = 'vdj'
    df_input = pd.read_csv(gliph_input, sep='\t')
    with open(network, 'r') as f:
        netw = [l.strip() for l in f.readlines() if len(l.strip()) > 0]
    cdr_graph, cdr_graph_trimed = build_gliph_graph(df_input, netw, label_known, label_unknown, return_full=True)
    predictions = [(n.split('-')[0], cdr_graph_trimed.nodes[n]['epitope']) for n in cdr_graph_trimed.nodes if f'-{label_unknown}' in n]
    df_predictions = pd.DataFrame(predictions, columns=['CDR3b', 'epitope'])
    df_highly_expressed = df_input[(df_input['source'] == label_unknown) & (df_input['TCR_clone.aafreq'] > 100)]
    df_highly_expressed = df_highly_expressed[df_highly_expressed['CDR3b'].isin(df_predictions['CDR3b'].values)].reset_index(drop=True)
    df_highly_expressed.sort_values(by='TCR_clone.aafreq', inplace=True, ascending=False)
    df_clusters = pd.DataFrame()
    for ndx, cdr in enumerate(df_highly_expressed['CDR3b'].values):
        cdrs_linked = [n for n, attr in cdr_graph.adj[cdr + f'-{label_unknown}'].items()]
        for n in cdrs_linked:
            if f'-{label_known}' in n:
                for n2 in cdr_graph.adj[n]:
                    if f'-{label_unknown}' in n2:
                        cdrs_linked.append(n2)
        cdrs_linked = list(set([n.split('-')[0] for n in cdrs_linked])) + [cdr]
        df_tmp = df_input[df_input['CDR3b'].isin(cdrs_linked)].sort_values(by='TCR_clone.aafreq', ascending=False)
        df_tmp['ClusterId'] = ndx + 1
        df_clusters = pd.concat([df_clusters, df_tmp])
    df_clusters.to_csv('clusters_highly_expressed.tsv', sep='\t')
