import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from collections import Counter
from Levenshtein import distance

def analyze_gliph(df_input, netw, label_known='known', label_unknown='unknown', return_full=False):
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
    gliph_input = 'input_gliph.tsv'
    network = gliph_input + '-clone-network.txt'
    label_unknown = 'unknown'
    label_known = 'known'
    df_input = pd.read_csv(gliph_input, sep='\t')
    with open(network, 'r') as f:
        netw = [l.strip() for l in f.readlines() if len(l.strip()) > 0]
    cdr_graph_trimed = analyze_gliph(df_input, netw, label_known, label_unknown)
    predictions = [(n.split('-')[0], cdr_graph_trimed.nodes[n]['epitope']) for n in cdr_graph_trimed.nodes if f'-{label_unknown}' in n]
    df_input = pd.read_csv(gliph_input, sep='\t')
    df_results = pd.DataFrame(predictions, columns=['CDR3b', 'epitope-pred'])
    true_values = [df_input[df_input['CDR3b'] == n[0]]['epitope'].values[0] for n in predictions]
    df_results['epitope-true'] = true_values
    df_results.to_csv('output_gliph_pred.tsv', sep='\t', index=False)
    #calculate true positive rate and false positive rate
    tcr_pos = pd.read_csv('../../tcr_pmhc_pos.csv')
    tcr_neg = pd.read_csv('../../tcr_pmhc_neg.csv')
    ture_pos = df_results[(df_results['CDR3b'].isin(tcr_pos['CDR3B'])) & (df_results['epitope-pred'] == df_results['epitope-true'])].shape[0]
    false_pos = df_results[(df_results['CDR3b'].isin(tcr_neg['CDR3B']))].shape[0] + \
                df_results[(df_results['CDR3b'].isin(tcr_pos['CDR3B'])) & (df_results['epitope-pred'] != df_results['epitope-true'])].shape[0]
    true_neg = df_input[df_input['CDR3b'].isin(tcr_neg['CDR3B'])].shape[0] - df_results[(df_results['CDR3b'].isin(tcr_neg['CDR3B']))].shape[0]
    false_neg = df_input[df_input['CDR3b'].isin(tcr_pos['CDR3B'])].shape[0] - ture_pos
    tpr = ture_pos / (ture_pos + false_neg)
    fpr = false_pos / (false_pos + true_neg)
    precision = ture_pos / (ture_pos + false_pos)
    recall = ture_pos / (ture_pos + false_neg)
    df_rate = pd.DataFrame([[tpr, fpr, precision, recall]], columns=['tpr', 'fpr', 'precision', 'recall'])
    df_rate.to_csv('tpr_fpr_gliph.tsv', sep='\t', index=False)
    print(f'true positive rate: {tpr}\nfalse positive rate: {fpr}')
