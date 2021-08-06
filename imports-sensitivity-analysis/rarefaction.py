"""
Perform rarefaction sensitivity analysis on African CoV data
Ancestral state reconstructions are performed using max parsimony

@author: david rasmussen
"""

from pathlib import Path
import numpy as np
import pandas as pd
import dendropy
from ete3 import Tree
from DateTimeUtils import date2FloatYear
import seaborn as sns
from matplotlib import pyplot as plt


"""
    Reconstruct ancestral states using Sankoff's max parsimony algorithm (see Felsenstein p15-18)
    Parsimony scores are computed assuming all transitions have a cost of one.
    Tip/ancestral features are input and output as dictionaries
    Internal nodes are labeled as n<X> where X is an int determined by the position of the node in a pre-order traversal
"""
def reconstruct_MP(tree,feature_dic):
    
    tree, tree_times = add_tree_times(tree)
    
    state_set = set([feature_dic[node.name] for node in tree.traverse() if node.is_leaf()])
    states = len(state_set)
    
    "Add a state2int dict map so this works with non-integer data types"
    state2int = {state:index for index, state in enumerate(state_set)}
    int2state = {index:state for index, state in enumerate(state_set)}
    
    "Post-order traversal to compute costs in terms of required state transitions"
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            costs = [np.inf]*states
            tip_state = state2int[feature_dic[node.name]] # feature_dic[node.name]
            costs[tip_state] = 0
            node.add_features(costs = costs)
        else:
            costs = [0]*states
            for i in range(states):
                child_costs = []
                for child in node.children:
                    temp_costs = [0]*states
                    for j in range(states):
                        temp_costs[j] = child.costs[j]
                        if i != j:
                            temp_costs[j] = temp_costs[j] + 1 # add cost for transitioning between i and j
                    child_costs.append(temp_costs)
                costs[i] = sum([min(c) for c in child_costs])
            node.add_features(costs = costs)
    
    "Pre-order traversal to select anc states based on least cost parsimony score"
    anc_dic = {}
    internal_cntr = 0 # counter for internal nodes encountered
    
    changes = 0
    times = []
    origins = []
    destinations = []
    
    for node in tree.traverse("preorder"):
        costs = node.costs
        if node.is_root():
            root_state = costs.index(min(costs)) # or np.argmin(node.costs)
            node.add_features(state = root_state)
            anc_dic['root'] = root_state
        else:
            parent_state = node.up.state
            least_cost = min(costs)
            least_cost_state = costs.index(least_cost)
            if parent_state == least_cost_state:
                anc_state = parent_state
            else:
                parent_st_cost = costs[parent_state]
                if parent_st_cost < (least_cost+1):
                    anc_state = parent_state # if parent state costs less than transitioning to least cost state
                else:
                    anc_state = least_cost_state
            node.add_features(state = anc_state)
            if node.is_leaf():
                anc_dic[node.name] = anc_state
            else:
                name = 'n' + str(internal_cntr)
                anc_dic[name] = anc_state
                internal_cntr += 1
            
            "Record direction of movements"
            if anc_state != parent_state:
                changes += 1
                times.append(node.time)
                origins.append(int2state[parent_state])
                destinations.append(int2state[anc_state])
    
    "Covert from integers back to original states"
    for k,v in anc_dic.items():
        anc_dic[k] = int2state[v]
        
    print("Total number of state changes: " + str(changes))
    df = pd.DataFrame({'EventTime' : times})
    df['Origin'] = origins
    df['Destination'] = destinations

    "Should return tree and anc_dic"
    return tree, df

"""
    Remove quotes surrounding taxa names in newick tree file
"""
def scrub_quotes(tree_file):
    
    tree_file = 'subsampled_temp.tre'
    tree_file_clean = tree_file
    f = open(tree_file)
    line = f.readline()
    out = ''
    while line:
        out += line.replace('\'','') # replace quotes in taxa names
        line = f.readline()
    f.close()
    nwk=open(tree_file_clean,"w")
    nwk.write(out + '\n')
    nwk.close()

"""
    Add node times/heights to ETE tree
"""    
def add_tree_times(tree):
        
    tree_times = []
    for i in tree.traverse():
        if i.is_root():
            i.add_features(time=0)
            tree_times.append(i.time)
        else:
            i.add_features(time=i.up.time + i.dist)
            tree_times.append(i.time)
    return tree, tree_times
    
      
"Set directories"
base_dir = Path(__file__).parent
#align_dir = Path.home() / 'Desktop' / 'msa_0314' # if different from base
tree_dir = base_dir / '2021-05-03_treetime'

"Set input files"
tree_file = str(tree_dir / "timetree.tre")
meta_file = str(base_dir / "metadata.tsv")

"Get metadata for global samples"
meta_df = pd.read_csv(meta_file,sep="\t",index_col='strain')

"Remove samples with unknown sampling dates"
clean_df = meta_df[meta_df['dates'].notna()]
meta_df_len = len(meta_df.index)
clean_df_len = len(clean_df.index)
if clean_df_len < meta_df_len:
    na_count = meta_df_len - clean_df_len
    meta_df = clean_df
    print('Removed ' + str(na_count) + ' samples with unknown (NA) dates') 
      
"Get sample dates"
meta_df['sample_datetime'] = pd.to_datetime(meta_df['dates'].astype('str')) # convert to pandas datetime
bins_dt = pd.date_range(start='2020-01-01',end='2021-06-01', freq='1M')
bin_labels = bins_dt.strftime('%b-%Y').tolist()[:-1]
meta_df['sample_month'] = pd.cut(meta_df['sample_datetime'], bins_dt, labels=bin_labels)
sample_month_counts = meta_df['sample_month'].value_counts()

"Load in global tree"
taxa = dendropy.TaxonNamespace()
global_tree = dendropy.Tree.get(file=open(tree_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa) 
if not global_tree.is_rooted:
    print('WARNING: Global tree is not rooted') # Should be rooted, this is just to check

"Get initial sample counts by region"
region_counts = meta_df['region'].value_counts()
print("Initial sample counts by region")
print(region_counts)
print()

"""
    Add external/internal classifier variable
    Unknown (NA's) will be treated as external
"""
externalMap = {'Asia':'External',
               'Europe':'External',
               'South America':'External',
               'Oceania':'External',
               'North America':'External',
               'Africa':'Internal'}
meta_df['region_type'] = meta_df['region'].apply(lambda x: externalMap[x] if isinstance(x, str) else 'External')
region_type_counts = meta_df['region_type'].value_counts()
print("External/internal counts")
print(region_type_counts)
print()

"Split meta_df by internal/external samples"
internal_df = meta_df[meta_df['region_type'] == 'Internal']
internal_taxa = internal_df.index.tolist() 
external_df = meta_df[meta_df['region_type'] == 'External']

"Main rarefaction params"
subsampling_loc = 'internal' # internal or external
sampling_fractions = np.arange(0.1,1.01,0.1) # sampling fracs to use in rarefaction analysis
boot_reps = 10 # number of replicates to run at each sampling frac

"Lists for storing results"
s_fractions = []
bootstrap_rep = []
imports = []
exports = []

"Run rarefaction"
for s in sampling_fractions:
    
    print("Starting analysis with sampling fraction = " + f'{s:.2f}')
    print()
    
    for r in range(boot_reps):
        
        print("Bootstrap rep = " + str(r))
        print()

        "Subsample a given fraction of internal or external samples"
        if subsampling_loc == 'internal':
            sub_df = internal_df.sample(frac=s, axis=0)
            merged_df = pd.concat([external_df, sub_df], axis=0) # Merge external df with subsampled external df 
        else: #external
            sub_df = external_df.sample(frac=s, axis=0)
            merged_df = pd.concat([internal_df, sub_df], axis=0) # Merge internal df with subsampled external df 
        
        subsampled_taxa = merged_df.index.tolist() 
        
        "Extract subtree with non-sampled taxa pruned"
        subsampled_tree_file = 'subsampled_temp.tre'
        taxa_to_retain = set([taxon for taxon in global_tree.taxon_namespace if taxon.label in subsampled_taxa])
        print("Pruning subsampled tree")
        subsampled_tree = global_tree.extract_tree_with_taxa(taxa=taxa_to_retain)
        subsampled_tree.write(path=subsampled_tree_file,schema='newick',suppress_annotations=True,suppress_rooting=True)
        
        "Get subsampled tree"
        scrub_quotes(subsampled_tree_file) # remove extra quotes dendropy adds to taxa names
        tree = Tree(subsampled_tree_file, format=1)
        
        "Run MP ancestral state reconstruction for all features"
        ancestral_tree_file = 'subsampled_temp_ancestral.tre'
        feature_dic = merged_df['region_type'].to_dict()
        print("Running ancestral state reconstuction")
        tree, changes_df = reconstruct_MP(tree,feature_dic)
        #tree.write(format=1, outfile=ancestral_tree_file)
        
        "Find imports"
        imports_df = changes_df[(changes_df['Origin'] == 'External') & (changes_df['Destination'] == 'Internal')]
        import_count = len(imports_df.index)
        print("Found " + str(import_count) + " imports")
        
        "Find exports"
        exports_df = changes_df[(changes_df['Origin'] == 'Internal') & (changes_df['Destination'] == 'External')]
        export_count = len(exports_df.index)
        print("Found " + str(export_count) + " exports")
        print()
        
        "Add results to output lists"
        s_fractions.append(s)
        bootstrap_rep.append(r)
        imports.append(import_count)
        exports.append(export_count)

df = pd.DataFrame({'SamplingFraction':s_fractions,'Imports':imports,'Exports':exports})
df.to_csv('rarefaction-results-internal.csv')

"""
    Plot introductions by sampling fraction
"""
sns.set(style="darkgrid")
fig, axs = plt.subplots(figsize=(6,4))
sns.lineplot(x='SamplingFraction',y='Imports',data=df, linewidth=2, dashes=False)

axs.set_xlabel('External sampling fraction', fontsize=14)
axs.set_ylabel('Imports into Africa', fontsize=14)
plt.tight_layout
plt.savefig('cov-imports-bySampleFrac-internal.png', dpi=300, bbox_inches="tight")
