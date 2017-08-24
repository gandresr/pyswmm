from pyswmm import Simulation, Nodes, Links
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from numpy import linalg as LA
from collections import defaultdict

# This script can be useful to estimate 
# the parameters of the virtual tanks model 
# for a urban drainage system modelled in 
# SWMM using pyswmm

def save_data_for_estimation(filename, inputs, groups={}):
    '''Saves the data required for the parameter estimation
    algorithm in Matlab (i.e., the ssest function)

    The SWMM model is approximated to a linear state-space 
    representation, ruled by the VT model dynamics:
    
    \dot{v} = q_in - q_out
    q_out = k*v

    in which runoff is the input of the system. 
    Hence, it is necessary to determine the discharge coefficients
    k for each virtual tank in the system, i.e., the coefficient
    associated to the total outflow and the volume of each tank.
    Also, the volumes (states) and the runoff signals (inputs) from
    SWMM.

    The ssest Matlab function will be used to estimate the
    unknown entries of A wich are related to discharge coefficients
    of flows between tanks.

    The function saves the required information in .csv files.

    Args:
        filename: path to the .inp SWMM file
        inputs: list of nodes with lateral inflows (i.e., the nodes
            where the input signals get into the system)
        groups: dict of groups of retention structures
            groups = {'group_id' : 
                            {'tanks': [tank1, tank2, ...],
                             'outflows': [tank1, tank2, ...]}, ...}
    Returns:
        G: graph representing interactions among states (i.e., virtual tanks
            volumes)
    '''

    G = graph(filename, groups)
    
    for node in G.nodes():
        G.node[node]['volume'] = []
        G.node[node]['outflow'] = []
    
    # Collect simulation results    
    volumes = defaultdict(list)
    inflows = defaultdict(list) # Only for nodes
    outflows = defaultdict(list) # Only for links
    runoffs = defaultdict(list) # Only for input nodes
    link_ids = []
    node_ids = []
    # Extract link and node IDs
    with Simulation(filename) as sim:
        link_ids = [l.linkid for l in Links(sim)]
        node_ids = [n.nodeid for n in Nodes(sim)]
        print(sim.flow_units)
    # Extract results
    with Simulation(filename) as sim:
        links = Links(sim)
        nodes = Nodes(sim)
        for step in sim:
            for link in link_ids:
                l = links[link]
                volumes[link].append(l.volume)
                outflows[link].append(l.flow)
            for node in node_ids:
                n = nodes[node]
                volumes[node].append(n.volume)
                inflows[node].append(n.total_inflow)
                if node in inputs:    
                    runoffs[node].append(n.lateral_inflow)
               
    # Updates volumes and outflows in G to estimate K
    for node in G.nodes():
        volume = []
        outflow = []
        # IMPORTANT
        # When the graph is created it is guaranteed
        #   that only storage nodes and conduit links
        #   are part of the states interaction graph.
        
        # Handle groups
        tlist = G.node[node]['tanks']
        olist = G.node[node]['outflows']
        if len(tlist) > 0:
            for t in tlist:
                # Update volume
                if len(volume) > 0:
                    volume += np.array(volumes[t])
                else:
                    volume = np.array(volumes[t])
            
            # Update outflow 
            for ot in olist:
                if ot in inflows.keys(): # it is a node   
                    dv = np.diff(volumes[ot])
                    dv = np.append(dv, dv[-1]) # to operate with proper dimensions
                    if len(outflow) > 0:
                        outflow += inflows[ot] - dv
                    else:
                        outflow = inflows[ot] - dv
                if ot in outflows.keys(): # it is a link
                    if len(outflow) > 0:
                        outflow += np.array(outflows[ot])
                    else:
                        outflow = np.array(outflows[ot])
        else:
            if node in inflows.keys(): # it is a node
                dv = np.diff(volumes[node])
                dv = np.append(dv, dv[-1]) # to operate with proper dimensions
                outflow = inflows[node] - dv
            if node in outflows.keys(): # it is a link
                outflow = np.array(outflows[node])
            volume = np.array(volumes[node])
        G.node[node]['volume'] = volume
        G.node[node]['outflow'] = outflow
    
    # Estimation of K
    K = {}
    for node in G.nodes():
        K[node] = estimate_K(G.node[node]['volume'], G.node[node]['outflow'])
        print(node)
    
    # TODO confirm that input values are the proper ones
    plt.plot(runoffs['2'])
    plt.show()
    return K

def estimate_K(V, Q, plot=True):
    '''Returns the value of the constant k for the virtual
    tanks (VT) model, i.e.,

    q_out = k * V

    Args:
        Q: outflow vector
    
        V: volume vector
    Returns:
        k: constant of the VT model
        
    
    '''
    Q = np.array(Q)
    V = np.array(V)
    
    def fcost(k):
        return LA.norm(Q - k*V, 2)        
    
    cons = [{'type': 'ineq', 'fun': lambda x:  1-x},
            {'type': 'ineq', 'fun': lambda x:  x}]

    k0 = 0
    k = minimize(fcost, k0, constraints=cons)
    if plot:
        plt.plot(k.x*V)
        plt.plot(Q)
        plt.show()
    return k.x

def graph(filename, groups={}):
    '''Returns the graph representation of the states interaction
    in the drainage network
    
    Args:
        filename: path to the .inp SWMM file
        groups: dict of groups of retention structures
            groups = {'group_id' : 
                            {'tanks': [tank1, tank2, ...],
                             'outflows': [tank1, tank2, ...]}, ...}
    Returns:
        G: digraph representing interactions between states in the
            drainage network
    '''
    G = nx.DiGraph()
    
    def del_connect(n):
        '''Deletes a node and connects its upstream nodes
        with its downstream nodes
        '''
        Gr = G.reverse()
        upstream = Gr.neighbors(n)
        downstream = G.neighbors(n)
        G.remove_node(n)
        for l1 in upstream:
            for l2 in downstream:
                G.add_edge(l1, l2)
    
    def replace(n, new):
        Gr = G.reverse()
        upstream = Gr.neighbors(n)
        downstream = G.neighbors(n)
        for u in upstream:
            G.add_edge(u, new)
        for d in downstream:
            G.add_edge(new, d)
        G.remove_node(n)
    
    
    upstream = {}; downstream = {}
    with Simulation(filename) as sim:
        nodes = Nodes(sim)
        links = Links(sim)
       
        for link in links:
            n = link.connections # (upstream, downstream)
            if nodes[n[0]].is_storage():
                try:
                    downstream[n[0]+'A'].append(link.linkid)
                except:
                    downstream[n[0]+'A'] = [link.linkid]
                try:
                    upstream[n[0]+'A'].append(n[0])
                except:
                    upstream[n[0]+'A'] = [n[0]]    
            else:
                try:
                    downstream[n[0]].append(link.linkid)
                except:
                    downstream[n[0]] = [link.linkid]
            if nodes[n[1]].is_storage():
                try:
                    downstream[n[1]+'B'].append(n[1])
                except:
                    downstream[n[1]+'B'] = [n[1]]    
                try:
                    upstream[n[1]+'B'].append(link.linkid)
                except:
                    upstream[n[1]+'B'] = [link.linkid]
            else:
                try:
                    upstream[n[1]].append(link.linkid)
                except:
                    upstream[n[1]] = [link.linkid]
        
        for n in upstream.keys():
            try:
                for l1 in upstream[n]:
                    for l2 in downstream[n]:
                        G.add_edge(l1, l2)
            except:
                continue

        for link in G.nodes():
            try:
                if not links[link].is_conduit():
                    del_connect(link)
            except:
                continue
        
        for link in G.nodes():
            G.node[link]['tanks'] = []
            G.node[link]['outflows'] = []
        if groups:
            for group in groups.keys():
                tlist = groups[group]['tanks'] # Tanks list
                olist =  groups[group]['outflows'] # Tanks for outflows list
                for i in range(len(tlist)-1):
                    del_connect(tlist[i])
                replace(tlist[-1], group)
                G.node[group]['tanks'] = tlist
                G.node[group]['outflows'] = olist
    return G


# TODO there is a problem with 4states_simple when there is a storage tank
fname = 'C:/Users/ga.riano949/Documents/GitHub/mpc/Paper ACC 2018/models/4 states/4states_simple.inp'
g = {'T2' : 
        {
            'tanks':
                ['T2A', 'T2B'], 
            'outflows': 
                ['T2B']
        }
    }
# G = graph(fname, g)
# K = save_data_for_estimation(fname, ('2','3',), groups=g)
K = save_data_for_estimation(fname, ('2','3',))
print(K)