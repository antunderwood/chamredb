# Superseded by the `analysis` package (python -m analysis.run), which reproduces this script's behaviour (see analysis/tests/README.md) with the
# decisions of ADR 0006. Kept as the source of the characterisation test. It needs networkx below 3.6 and reads the bundled graph.
import os
import json
import re
from scipy import stats

import networkx as nx
from networkx.readwrite import json_graph

from strsimpy import NormalizedLevenshtein

import numpy as np

import seaborn as sns
from matplotlib import pyplot as plt

import pandas as pd

from chamredb.functions.graph_visualisation_functions import __best_target_node_match, __node_target_info
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None) 
pd.set_option('display.max_colwidth', None)

json_path=os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..", "chamredb", "data", "graph.json"
        )
with open(json_path) as json_file:
    G=json_graph.node_link_graph(json.load(json_file))


def clean_name(name,database,remove_alleles=False):
    # remove species names from card entries
    # card pattern with species
    species_pattern_one=r"^[A-Z][a-z]+ [a-z]+ (.+) conferring resistance to .+$"
    species_pattern_two=r"^(.+) from [A-Z][a-z]+ [a-z]$"
    species_pattern_three=r"^[A-Z][a-z]+ [a-z]+ (.+)$"


    if database == 'card' and re.search(species_pattern_one, name):
        cleaned_name=re.search(species_pattern_one, name).group(1)
    elif database == 'card' and re.search(species_pattern_two, name):
        cleaned_name=re.search(species_pattern_one, name).group(1)
    elif database == 'card' and re.search(species_pattern_three, name):
        cleaned_name=re.search(species_pattern_three, name).group(1)
    else:
        cleaned_name=name

    # lower case
    cleaned_name=cleaned_name.lower()
    # remove bla prefix
    cleaned_name=cleaned_name.replace('bla', '')
    # remove parentheses
    cleaned_name=cleaned_name.replace('(', '')
    cleaned_name=cleaned_name.replace(')', '')
    # remove hyphens
    cleaned_name=cleaned_name.replace('-', '')
    
    if remove_alleles:
        allele_pattern=r"^([A-Za-z]+)[0-9]+[A-Za-z0-9.]*$"
        if re.search(allele_pattern, cleaned_name):
            cleaned_name=re.search(allele_pattern, cleaned_name).group(1)

    return(cleaned_name)


# calculate normalized levenshtein similarities and make dataframe
normalized_levenshtein=NormalizedLevenshtein()
distance_data=[]
blaPAO_count=0
for source_node in G.nodes():
    source_node_targets=__node_target_info(G, source_node)
    for target_database in source_node_targets:
        target_node_info=__best_target_node_match(source_node_targets[target_database], 0, 0)
        source_name=G.nodes[source_node]['name']
        source_database=G.nodes[source_node]['database']
        target_name=target_node_info['node']['name']
        target_database=target_node_info['node']['database']
        identity=target_node_info['edge_data']['identity']

        # blaPAO-N and blaPDC matches are the source of many low name similarities. Therefore skip
        if source_name == 'blaPAO' or target_name == 'blaPAO':
            blaPAO_count += 1
            continue
        cleaned_source_name=clean_name(source_name, source_database)
        cleaned_target_name=clean_name(target_name, target_database)
        normalized_similarity=normalized_levenshtein.similarity(cleaned_source_name, cleaned_target_name)

        distance_data.append(
            {
                'source_name': source_name,
                'cleaned_source_name': cleaned_source_name,
                'source_database': source_database,
                'target_name': target_name,
                'cleaned_target_name': cleaned_target_name,
                'target_database': target_database,
                'sequence_identity': identity,
                'name_similarity': normalized_similarity,
            }
        )
distance_dataframe=pd.DataFrame(distance_data)
print(f"Number of blaPAO matches={blaPAO_count}")


linear_fit=np.polyfit(
    distance_dataframe['sequence_identity'],
    distance_dataframe['name_similarity'],
    1
)

fitted_data=(
    distance_dataframe
    .assign(
        predicted_similarity=lambda x: np.polyval(linear_fit, x['sequence_identity']),
        similarity_difference=lambda x: abs(x['name_similarity'] - x['predicted_similarity'])
    )
)
# # get quantiles
# lower_quantile=fitted_data['similarity_difference'].quantile(0.025)
# upper_quantile=fitted_data['similarity_difference'].quantile(0.975)

# outlier_data=fitted_data.query(
#     'similarity_difference < @lower_quantile or similarity_difference > @upper_quantile'
# )
# The predicted name similarity for a sequence identity of 0.95 is 0.69
sequence_identity_threshold=0.95
predicted_name_similarity=np.polyval(linear_fit, sequence_identity_threshold)



high_seq_id_similarity=fitted_data.query('sequence_identity >= @sequence_identity_threshold')
low_name_similarity=fitted_data.query('sequence_identity >= @sequence_identity_threshold and name_similarity <= @predicted_name_similarity')


low_name_similarity_path=os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "low_name_similarity.csv"
)
low_name_similarity.to_csv(low_name_similarity_path, index=False)


def clean_name_in_dataframe(df, node_type):
    """
    clean_name_in_dataframe clean database column

    Args:
        df (pandas dataframe): input dataframe
        node_type (string): source or target

    Returns:
        pandas dataframe: dataframe with name cleaned
    """
    name_column=f"{node_type}_name"
    database_column=f"{node_type}_database"
    cleaned_name_column=f"cleaner_{node_type}_name"
    df=df.copy()
    df[cleaned_name_column]=df.apply(lambda x: clean_name(x[name_column], x[database_column], remove_alleles=True), axis=1)
    return df

def calculate_name_similarity(df):
    df=df.copy()
    df['name_similarity_without_alleles']=df.apply(
        lambda x: normalized_levenshtein.similarity(
            x['cleaner_source_name'],
            x['cleaner_target_name']
        ),
        axis=1
    )
    return df

# ignore alleles
distance_df_ignoring_alleles=(
    distance_dataframe
    .pipe(clean_name_in_dataframe, 'source')
    .pipe(clean_name_in_dataframe, 'target')
    .pipe(calculate_name_similarity)
)

linear_fit_ignoring_alleles=np.polyfit(
    distance_df_ignoring_alleles['sequence_identity'],
    distance_df_ignoring_alleles['name_similarity_without_alleles'],
    1
)

distance_df_ignoring_alleles=(
    distance_df_ignoring_alleles
    .assign(
        predicted_similarity=lambda x: np.polyval(linear_fit_ignoring_alleles, x['sequence_identity']),
        similarity_difference=lambda x: abs(x['name_similarity_without_alleles'] - x['predicted_similarity'])
    )
)

# The predicted name similarity for a sequence identity of 0.95 is 0.69
sequence_identity_threshold=0.95
predicted_name_similarity_without_alleles=np.polyval(linear_fit, sequence_identity_threshold)

high_seq_id_similarity_without_alleles=distance_df_ignoring_alleles.query('sequence_identity >= @sequence_identity_threshold')
low_name_similarity_without_alleles=distance_df_ignoring_alleles.query('sequence_identity >= @sequence_identity_threshold and name_similarity_without_alleles <= @predicted_name_similarity_without_alleles')

low_name_similarity_without_alleles_path=os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "low_name_similarity_without_alleles.csv"
)
low_name_similarity_without_alleles.to_csv(low_name_similarity_without_alleles_path , index=False)

####################### Plotting ########################
fig=plt.figure(figsize=(15,17))
fig.subplots_adjust(top=0.9)
fig.suptitle(
    'Comparison of sequence identity with name smilarity of matches within the graph',
    fontsize='xx-large'
)
distance_distribution_plot_ax=plt.subplot2grid((2, 2), (0, 0), rowspan=1, colspan=1)
distance_distribution_plot_ax.set_title('A: Distribution of normalised Levenshtein\ndistances between match names')

distance_distribution_plot=sns.histplot(
    ax=distance_distribution_plot_ax,
    x="name_similarity",
    data=distance_dataframe,
    binwidth=0.05
)

identity_distribution_plot_ax=plt.subplot2grid((2, 2), (0, 1), rowspan=1, colspan=1)
identity_distribution_plot_ax.set_title('B: Distribution of sequence\nidentities between match')

identity_distribution_plot=sns.histplot(
    ax=identity_distribution_plot_ax,
    x="sequence_identity",
    data=distance_dataframe,
    binwidth=0.05
)

# get coeffs of linear fit
slope, intercept, r_value, p_value, std_err = stats.linregress(distance_dataframe['sequence_identity'],distance_dataframe['name_similarity'])

identity_vs_distance_plot_ax=plt.subplot2grid((2, 2), (1, 0), rowspan=1, colspan=2)
identity_vs_distance_plot_ax.set_title('C: Normalised Levenshtein distances\nvs sequence identities')

identity_vs_distance_plot=sns.regplot(
    ax=identity_vs_distance_plot_ax,
    x="sequence_identity",
    y="name_similarity",
    scatter_kws={"s": 10, 'alpha': 0.1},
    line_kws={"color": "red"},
    label="y={0:.2f}x+{1:.2f}. R^2={2:.2f}".format(slope,intercept,r_value),
    data=distance_dataframe,
    marker="+").legend(loc="best", fontsize=15)

identity_vs_distance_plot_ax.set_xlim(0,1)
identity_vs_distance_plot_ax.set_ylim(0,1)

figure_path=os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "seq_id_and_name_sim_analysis.png"
)
fig.savefig(figure_path)

# distribution of name smilarities
fig=plt.figure(figsize=(15,17))

fig.suptitle(
    "Analysis of distributions of name similarities and differences\nbetween observed and predicted name similarities",
    fontsize='xx-large'
)
########## violin plot 1 ###############
name_similarity_plot_ax=plt.subplot2grid((2, 4), (0, 0), rowspan=1, colspan=2)
name_similarity_plot_ax.set_title(f'A: Dist of name similarities including alleles\nfrom matches with >= 0.95 id')

sns.boxplot(
    ax=name_similarity_plot_ax,
    y="name_similarity",
    data=fitted_data.query('sequence_identity >= 0.95'),
    showfliers=False,
    showbox=False,
    whis=[2.5,97.5]
)
sns.violinplot(
    ax=name_similarity_plot_ax,
    y="name_similarity",
    cut=0.0,
    data=fitted_data.query('sequence_identity >= 0.95')
)

########## violin plot 2 ###############
name_similarity_wo_alleles_plot_ax=plt.subplot2grid((2, 4), (0, 2), rowspan=1, colspan=2)
name_similarity_wo_alleles_plot_ax.set_title(f'B: Dist of name similarities without alleles\nfrom matches with >= 0.95 id')

sns.boxplot(
    ax=name_similarity_wo_alleles_plot_ax,
    y="name_similarity_without_alleles",
    data=distance_df_ignoring_alleles.query('sequence_identity >= 0.95'),
    showfliers=False,
    showbox=False,
    whis=[2.5,97.5]
)
sns.violinplot(
    ax=name_similarity_wo_alleles_plot_ax,
    y="name_similarity_without_alleles",
    cut=0.0,
    data=distance_df_ignoring_alleles.query('sequence_identity >= 0.95')
)
########## violin plot 3 ###############
similarity_difference_plot_ax=plt.subplot2grid((2, 4), (1, 0), rowspan=1, colspan=1)
similarity_difference_plot_ax.set_title(
f"""C: Dist of differences from 
linear regression predictions incl 
alleles from matches with >= 0.95 id""")

sns.boxplot(
    ax=similarity_difference_plot_ax,
    y="similarity_difference",
    data=high_seq_id_similarity,
    showfliers=False,
    showbox=False,
    whis=[2.5,97.5]
)
sns.violinplot(
    ax=similarity_difference_plot_ax,
    y="similarity_difference",
    cut=0.0,
    data=high_seq_id_similarity
)

########## violin plot 4 ###############
low_similarity_difference_plot_ax=plt.subplot2grid((2, 4), (1, 1), rowspan=1, colspan=1)
low_similarity_difference_plot_ax.set_title(
f"""D: Dist of differences from linear 
regr. predictions incl. alleles from matches 
with >= 0.95 id + name similarity <= {round(predicted_name_similarity,2)}""")

sns.boxplot(
    ax=low_similarity_difference_plot_ax,
    y="similarity_difference",
    data=low_name_similarity,
    showfliers=False,
    showbox=False,
    whis=[2.5,97.5]
)
sns.violinplot(
    ax=low_similarity_difference_plot_ax,
    y="similarity_difference",
    cut=0.0,
    data=low_name_similarity
)

########## violin plot 5 ###############
similarity_difference_wo_alleles_plot_ax=plt.subplot2grid((2, 4), (1, 2), rowspan=1, colspan=1)
similarity_difference_wo_alleles_plot_ax.set_title(
f"""E: Dist of differences from linear 
regr. predictions without 
alleles from matches with >= 0.95 id"""
)

sns.boxplot(
    ax=similarity_difference_wo_alleles_plot_ax,
    y="similarity_difference",
    data=high_seq_id_similarity_without_alleles,
    showfliers=False,
    showbox=False,
    whis=[2.5,97.5]
)
sns.violinplot(
    ax=similarity_difference_wo_alleles_plot_ax,
    y="similarity_difference",
    cut=0.0,
    data=high_seq_id_similarity_without_alleles
)

########## violin plot 6 ###############
low_similarity_difference_wo_alleles_plot_ax=plt.subplot2grid((2, 4), (1, 3), rowspan=1, colspan=1)
low_similarity_difference_wo_alleles_plot_ax.set_title(
f"""F: Dist of differences from linear regr.
predictions without alleles from matches
with >= 0.95 id + name similarity <= {round(predicted_name_similarity_without_alleles,2)}"""
)

sns.boxplot(
    ax=low_similarity_difference_wo_alleles_plot_ax,
    y="similarity_difference",
    data=low_name_similarity_without_alleles,
    showfliers=False,
    showbox=False,
    whis=[2.5,97.5]
)
sns.violinplot(
    ax=low_similarity_difference_wo_alleles_plot_ax,
    y="similarity_difference",
    cut=0.0,
    data=low_name_similarity_without_alleles
)

plt.tight_layout()
fig.subplots_adjust(top=0.9, right=0.95)
figure_path=os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "name_similarity_distribution_analysis.png"
)
fig.savefig(figure_path)
