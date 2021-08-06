"""
Perform rarefaction sensitivity analysis on African CoV data
Ancestral state reconstructions are performed using max parsimony

To do:
        -Try using time tree
        -Record times of change events
        -Plot import/export times
        -Subsample by time, capping at maximum 2020

@authors: david rasmussen
"""

from pathlib import Path
import numpy as np
import pandas as pd
import dendropy
from ete3 import Tree
from DateTimeUtils import date2FloatYear
from DateTimeUtils import floatYear2Date

import seaborn as sns
from matplotlib import pyplot as plt


"""
    Reconstruct ancestral states using Sankoff's max parsimony algorithm (see Felsenstein p15-18)
    Parsimony scores are computed assuming all transitions have a cost of one.
    Tip/ancestral features are input and output as dictionaries
    Internal nodes are labeled as n<X> where X is an int determined by the position of the node in a pre-order traversal
"""
def reconstruct_MP(tree,feature_dic):
    
    #tree, tree_times = add_tree_times(tree)
    
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
    Remove quotes surrounding taxa names
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

def make_map(meta_df):
    
    "Get country to region map"
    country2RegionMap = {}
    for sample, row in meta_df.iterrows():
        country = row['country']
        if country not in country2RegionMap:
            country2RegionMap[country] = row['region']
    pd.DataFrame.from_dict(data=country2RegionMap, orient='index').to_csv('country2RegionMap.csv', header=False)
      
"Set directories"
base_dir = Path(__file__).parent
#tree_dir = base_dir / '2021-05-03_treetime' # old time tree
tree_dir = base_dir / 'newtimetree' # Eduan's newer tree

"Set input files"
#tree_file = str(base_dir / "Africa_COVID19_2021.04.28.mafft.aln.fasta.treefile")
tree_file = str(tree_dir / "timetree.nwk")
meta_file = str(tree_dir / "metadata.tsv")

"Get metadata for global samples"
meta_df = pd.read_csv(meta_file,sep="\t",index_col='strain')

"Get country to region map"
#make_map(meta_df) # make country-to-region map
country_df = pd.read_csv('country2RegionMap.csv', index_col=0)
country2Region = {index:row['Region'] for index, row in country_df.iterrows()}

"Remove samples with unknown sampling dates"
clean_df = meta_df[meta_df['dates'].notna()]
meta_df_len = len(meta_df.index)
clean_df_len = len(clean_df.index)
if clean_df_len < meta_df_len:
    na_count = meta_df_len - clean_df_len
    meta_df = clean_df
    print('Removed ' + str(na_count) + ' samples with unknown (NA) dates') 
      
"Get sample dates"
#meta_df['sample_datetime'] = pd.to_datetime(meta_df['dates'].astype('str')) # convert to pandas datetime
meta_df['sample_datetime'] = pd.to_datetime(floatYear2Date(meta_df['dates'].tolist())) # convert to pandas datetime
final_sample_date = meta_df['sample_datetime'].max()
bins_dt = pd.date_range(start='2019-11-01',end='2021-02-01', freq='1M')
months = bins_dt.strftime('%b-%Y').tolist()[1:]
meta_df['sample_month'] = pd.cut(meta_df['sample_datetime'], bins_dt, labels=months)
sample_month_counts = meta_df['sample_month'].value_counts()

"Load in global tree"
taxa = dendropy.TaxonNamespace()
global_tree = dendropy.Tree.get(file=open(tree_file, 'r'), schema="newick", rooting="default-rooted", preserve_underscores=True, taxon_namespace=taxa) 
if not global_tree.is_rooted:
    print('WARNING: Global tree is not rooted') # Should be rooted, this is just to check

"Run checks making sure all samples in meta_df are in tree and all samples in tree are in meta_df"
meta_set = set(meta_df.index.tolist())
tree_set = set([taxon.label for taxon in global_tree.taxon_namespace])
missing_in_tree = meta_set.difference(tree_set) # samples present in meta_set but not in tree_set
if missing_in_tree:
    print('WARNING: There are ' + str(len(missing_in_tree)) + ' present in meta data but not in the tree taxa')
missing_in_meta = tree_set.difference(meta_set) # samples present in tree_set but not in meta_set
if missing_in_meta:
    print('WARNING: There are ' + str(len(missing_in_meta)) + ' present in tree taxa but not in the meta data')

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
internal_df = meta_df[meta_df['region_type'] == 'Internal'] # african samples
external_df = meta_df[meta_df['region_type'] == 'External'] # non-african samples

internal_sample_month_counts = internal_df['sample_month'].value_counts()

"Main sensitivity params"
boot_reps = 1
max_sample_count = 400 # max threshold sample count for each month -- set to np.inf for non-uniform sampling
rec_state = 'country' # state to reconstruct (region, region_type, country)

"Lists for storing results"
boot_dfs = []

print("Starting sensitivity analysis")
print()

for r in range(boot_reps):
    
    print("Bootstrap rep = " + str(r))
    print()

    "Subsample a given number of internal (African) samples each month"
    sub_internal_df = pd.DataFrame(columns=meta_df.columns)
    for m in months:
        month_df = internal_df[internal_df['sample_month'] == m]
        month_count = len(month_df.index)
        sample_count = min(month_count,max_sample_count)
        sub_month_df = month_df.sample(n=sample_count, axis=0)
        sub_internal_df = pd.concat([sub_internal_df, sub_month_df], axis=0)
    
    #sample_month_counts = sub_internal_df['sample_month'].value_counts()
    
    "Merge external df with subsampled external df"
    merged_df = pd.concat([external_df, sub_internal_df], axis=0)
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
    tree, tree_times = add_tree_times(tree) # tree times are heights/distances from root
    
    "Run MP ancestral state reconstruction for all features"
    ancestral_tree_file = 'subsampled_temp_ancestral.tre'
    feature_dic = merged_df[rec_state].to_dict()
    print("Running ancestral state reconstuction")
    tree, changes_df = reconstruct_MP(tree,feature_dic)
    #tree.write(format=1, outfile=ancestral_tree_file)
    
    boot_dfs.append(changes_df)

"""
    Plot introduction times by month
"""

df = boot_dfs[0]
df['Origin_Region'] = df['Origin'].apply(lambda x: country2Region[x])
df['Destination_Region'] = df['Destination'].apply(lambda x: country2Region[x])

"Convert float times back to date times"
final_sample_float = date2FloatYear(final_sample_date)
abs_root_time = final_sample_float - max(tree_times)
abs_times = list(abs_root_time + df['EventTime'].to_numpy())
#abs_times = list(final_sample_float - df['EventTime'].to_numpy())
df['event_datetime'] = pd.to_datetime(floatYear2Date(abs_times))
df['event_month'] = pd.cut(df['event_datetime'], bins_dt, labels=months)

internal_event_df = df[(df['Origin_Region'] == 'Africa') & (df['Destination_Region'] == 'Africa')]
internal_event_df['IntroType'] = 'Internal'
external_event_df = df[(df['Origin_Region'] != 'Africa') & (df['Destination_Region'] == 'Africa')]
external_event_df['IntroType'] = 'External'

internal_month_counts = internal_event_df['event_month'].value_counts()
external_month_counts = external_event_df['event_month'].value_counts()
proportion_external = external_month_counts / (internal_month_counts + external_month_counts)

"Plot temporal distribution of sample times by month"
sns.set(style="darkgrid")
fig, axs = plt.subplots(figsize=(6,4))
sns.barplot(x=proportion_external.index, y=proportion_external.values, label="Total", color="b")
axs.set_ylim([0,1.0])
axs.set_xlabel('Time', fontsize=12)
axs.set_ylabel('Proportion external', fontsize=12)
axs.set_title('Uniform sampling')
plt.xticks(rotation = 90)
plt.tight_layout
plt.savefig('cov-introductions-proportionExternal-byTime-uniform-newTimeTree.png', dpi=300, bbox_inches="tight")

