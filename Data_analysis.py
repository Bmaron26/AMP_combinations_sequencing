

# # Analysis of Mutations and Profile Similarity in Evolved Strains
#
# This script analyzes mutation data from strains evolved under different
# antimicrobial peptide (AMP) treatments (single and combinations).
# It assesses:
# 1. Differences in the number of mutations acquired.
# 2. Similarity (Dice index) of mutation profiles within and between treatment groups.
# Statistical significance is often assessed using non-parametric tests
# (Mann-Whitney U, Kruskal-Wallis) and permutation tests.

# ## 1. Setup: Imports and Configuration

import pandas as pd
import numpy as np
import scipy.stats as stats
from scipy.stats import mannwhitneyu, kruskal
from scikit_posthocs import posthoc_dunn
from sklearn.metrics import pairwise_distances
from itertools import combinations
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
from matplotlib.colors import LinearSegmentedColormap

# Configure plots
plt.style.use('seaborn-v0_8-ticks') # Use seaborn style
sns.set_style('ticks')
sns.set_context('paper', font_scale=1.5) # Adjust font scale as needed

# Ignore common warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning) # Can be more specific if needed

# --- Configuration ---
# Input files
MUTATION_COUNTS_FILE = 'mutations_per_strain_data_new_paul.xlsx'
MUTATION_PROFILES_FILE = 'mutations_and_resistance_correlation_data_new_paul.xlsx'
PROFILES_SHEET_NAME = 'resistance'

# Output files / settings
OUTPUT_DIR = '.' # Save outputs in the current directory
NUM_PERMUTATIONS = 2000 # Number of permutations for significance testing
FDR_ALPHA = 0.05 # Significance level for FDR corrected p-values (e.g., Dunn's test)

# Define treatments and color palette consistently
TREATMENTS = ['Temp', 'Mel', 'Pex', 'Temp_Mel', 'Pex_Temp', 'Pex_Mel']
SINGLE_AMPS = ['Pex', 'Temp', 'Mel']
COMB_AMPS = ['Temp_Mel', 'Pex_Temp', 'Pex_Mel']

# Consistent color palette
# Start with a base palette
base_palette = sns.color_palette('deep', n_colors=len(TREATMENTS))
COLOR_PALETTE = {treatment: color for treatment, color in zip(TREATMENTS, base_palette)}
# Customize specific combinations if desired (as in original script)
COLOR_PALETTE['Pex_Mel'] = 'crimson'
COLOR_PALETTE['Temp_Mel'] = 'gold'
# Add gray for general boxplot outlines etc. if needed
COLOR_PALETTE['outline'] = 'gray'
COLOR_PALETTE['boxplot_single'] = 'skyblue'
COLOR_PALETTE['boxplot_comb'] = 'lightyellow'
COLOR_PALETTE['boxplot_within'] = 'white'


# Set pandas display options
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)

print("Setup Complete.")

# ## 2. Load and Prepare Data

print(f"Loading mutation counts from: {MUTATION_COUNTS_FILE}")
df_mut_counts = pd.read_excel(MUTATION_COUNTS_FILE)
# Basic validation
required_cols_counts = ['treatment', 'single_comb', 'mutations']
if not all(col in df_mut_counts.columns for col in required_cols_counts):
    raise ValueError(f"Mutation counts file must contain columns: {required_cols_counts}")
print(f"Loaded {len(df_mut_counts)} records for mutation counts.")

print(f"\nLoading mutation profiles from: {MUTATION_PROFILES_FILE}, sheet: {PROFILES_SHEET_NAME}")
# Load, transpose, and select relevant columns for similarity analysis
# Assuming index_col=0 is the strain identifier before transposing
df_profiles_raw = pd.read_excel(MUTATION_PROFILES_FILE, sheet_name=PROFILES_SHEET_NAME, index_col=0).T
# Assuming the first 3 columns are metadata not needed for similarity calc (as per original iloc[:,3:])
# Adjust this slice if your data structure is different
df_dor_similarity = df_profiles_raw.iloc[:, 3:].copy()
# Ensure data is numeric (binary 0/1 for Dice)
df_dor_similarity = df_dor_similarity.astype(int)

# Add treatment column based on index (strain name assumed to be 'TreatmentTypeX')
df_dor_similarity['treatment'] = [idx[:-1] for idx in df_dor_similarity.index]
print(f"Loaded mutation profiles for {len(df_dor_similarity)} strains.")
print(f"Unique treatments found in profiles: {df_dor_similarity['treatment'].unique()}")

# Save the prepared similarity data matrix if needed (optional)
# df_dor_similarity.to_csv(f"{OUTPUT_DIR}/prepared_mutation_profiles_matrix.csv")


# ## 3. Analysis: Mutation Counts per Strain

print("\n--- Starting Mutation Count Analysis ---")

# ### 3.1. Single vs. Combination Treatment Mutation Counts

single_group_mut = df_mut_counts[df_mut_counts['single_comb'] == 's']['mutations']
combination_group_mut = df_mut_counts[df_mut_counts['single_comb'] == 'c']['mutations']

# Perform Mann-Whitney U test
u_statistic, p_value_mw = mannwhitneyu(single_group_mut, combination_group_mut, alternative='two-sided') # Explicitly two-sided

print("\nMann-Whitney U Test (Single vs. Combination Mutations):")
print(f'  U statistic: {u_statistic:.2f}')
print(f'  P-value: {p_value_mw:.4f}')

# Create Boxplot with Swarmplot (Single vs. Combination)
plt.figure(figsize=(8, 6))
sns.boxplot(x='single_comb', y='mutations', data=df_mut_counts,
            width=0.4, fill=False, color=COLOR_PALETTE['outline'], showfliers=False)
sns.swarmplot(data=df_mut_counts, x='single_comb', y='mutations', hue='treatment',
              palette=COLOR_PALETTE, edgecolor='black', linewidth=0.5, legend=True)

plt.xlabel('Treatment Type', fontsize=18)
plt.ylabel('Mutations per Population', fontsize=18)
plt.xticks([0, 1], ['Single', 'Combination'], fontsize=16)
plt.yticks(fontsize=14)
plt.legend(title='Treatment', bbox_to_anchor=(1.05, 1), loc='upper left', edgecolor='white')
plt.tight_layout(rect=[0, 0, 0.85, 1]) # Adjust layout for legend

plot_filename_s_vs_c = f"{OUTPUT_DIR}/mutations_single_vs_combination_boxplot_swarm.png"
plt.savefig(plot_filename_s_vs_c, dpi=300, bbox_inches='tight')
print(f"Saved plot: {plot_filename_s_vs_c}")
plt.close() # Close figure to free memory

# ### 3.2. Mutation Counts Across All Treatments

print("\nKruskal-Wallis H-Test (Mutations Across All Treatments):")
# Prepare data for Kruskal-Wallis
samples_mut = [group['mutations'].values for name, group in df_mut_counts.groupby('treatment')]

# Perform Kruskal-Wallis H-test
h_statistic, p_value_kw = stats.kruskal(*samples_mut)
print(f'  H-statistic: {h_statistic:.4f}')
print(f'  P-value: {p_value_kw:.4f}')

# Perform Dunn's post-hoc test if Kruskal-Wallis is significant
if p_value_kw < FDR_ALPHA:
    print("\nPerforming Dunn's post-hoc test (p-value adjustment: fdr_bh):")
    dunn_results = posthoc_dunn(df_mut_counts, val_col='mutations', group_col='treatment', p_adjust='fdr_bh')
    # Save Dunn's test results
    dunn_filename = f"{OUTPUT_DIR}/dunn_test_results_mutations.csv"
    dunn_results.to_csv(dunn_filename)
    print(f"  Saved Dunn's test results to: {dunn_filename}")
    print("  Dunn's Test Results (P-values):")
    print(dunn_results)
else:
    print("\nKruskal-Wallis test not significant; Dunn's test not performed.")

# Create Boxplot with Stripplot (All Treatments)
plt.figure(figsize=(10, 7)) # Adjusted size for rotation
sns.boxplot(x='treatment', y='mutations', data=df_mut_counts,
            width=0.6, palette=COLOR_PALETTE, showfliers=False, order=TREATMENTS) # Use consistent order
sns.stripplot(x='treatment', y='mutations', data=df_mut_counts,
              color='black', size=5, alpha=0.6, jitter=True, order=TREATMENTS) # Use consistent order

plt.xlabel('Treatment', fontsize=18)
plt.ylabel('Mutations per Population', fontsize=18)
plt.xticks(rotation=45, ha='right', fontsize=14) # Rotate for readability
plt.yticks(fontsize=14)
plt.tight_layout()

plot_filename_all_treat = f"{OUTPUT_DIR}/mutations_across_treatments_boxplot_strip.png"
plt.savefig(plot_filename_all_treat, dpi=300)
print(f"\nSaved plot: {plot_filename_all_treat}")
plt.close()

print("\n--- Mutation Count Analysis Complete ---")


# ## 4. Analysis: Mutation Profile Similarity

print("\n--- Starting Mutation Profile Similarity Analysis ---")

# ### 4.1. Helper Function: Calculate Pairwise Similarity

def get_pair_similarity(df_profiles):
    """
    Calculates pairwise Dice similarity between strains based on mutation profiles.

    Args:
        df_profiles (pd.DataFrame): DataFrame where rows are strains, columns are
                                   mutation features (e.g., 0/1). Must contain a
                                   'treatment' column derived from the index.

    Returns:
        pd.DataFrame: DataFrame containing pairwise similarity information with columns:
                      ['strain1', 'strain2', 'similarity', 'treatment1', 'treatment2',
                       'within_between', 'n1', 'n2', 'AMP11', 'AMP12', 'AMP21', 'AMP22']
                      where n1/n2 are number of AMPs in treatment1/2, and AMPxx are
                      the component AMPs (assuming max 2).
    """
    if 'treatment' not in df_profiles.columns:
        raise ValueError("Input DataFrame must have a 'treatment' column.")

    # Prepare data matrix (numeric profiles only)
    df_matrix = df_profiles.drop('treatment', axis=1)

    # Calculate Dice similarity (1 - Dice distance)
    # Ensure data is boolean or binary for dice metric
    dist_matrix = pairwise_distances(df_matrix.values.astype(bool), metric='dice')
    similarity_matrix = pd.DataFrame(1 - dist_matrix, columns=df_matrix.index, index=df_matrix.index)

    data = []
    strains = df_matrix.index
    treatments_dict = df_profiles['treatment'].to_dict()

    for i, (s1, s2) in enumerate(combinations(strains, 2)):
        sim = similarity_matrix.loc[s1, s2]
        t1 = treatments_dict[s1]
        t2 = treatments_dict[s2]
        within_between = 'w' if t1 == t2 else 'b'

        amps1 = t1.split('_')
        n1 = len(amps1)
        amp11 = amps1[0]
        amp12 = amps1[0] if n1 == 1 else amps1[1] # Handles single AMP case

        amps2 = t2.split('_')
        n2 = len(amps2)
        amp21 = amps2[0]
        amp22 = amps2[0] if n2 == 1 else amps2[1] # Handles single AMP case

        # Original column names were slightly confusing (AMP11 vs AMP21) - let's clarify
        # AMP11, AMP12 = components of treatment 1
        # AMP21, AMP22 = components of treatment 2
        row = [s1, s2, sim, t1, t2, within_between, n1, n2, amp11, amp12, amp21, amp22]
        data.append(row)

    # Adjusting column names for clarity based on definition
    columns = ['strain1', 'strain2', 'similarity', 'treatment1', 'treatment2',
               'within_between', 'n1', 'n2', 'AMP1_1', 'AMP1_2', 'AMP2_1', 'AMP2_2']
    return pd.DataFrame(data, columns=columns)

# Calculate all pairwise similarities
print("\nCalculating all pairwise similarities...")
# Use the data prepared in step 2
data_similarity_all_pairs = get_pair_similarity(df_dor_similarity)
similarity_filename = f"{OUTPUT_DIR}/pairwise_similarity_all_pairs.csv"
data_similarity_all_pairs.to_csv(similarity_filename, index=False)
print(f"Calculated {len(data_similarity_all_pairs)} pairwise similarities.")
print(f"Saved all pairwise similarities to: {similarity_filename}")

# Filter for within-treatment pairs for specific analyses
df_within = data_similarity_all_pairs.query('within_between=="w"').copy()
# The 'treatment' column is implicitly present via treatment1/treatment2 which are equal here.
# We can add it explicitly if preferred for plotting etc.
df_within['treatment'] = df_within['treatment1']

within_similarity_filename = f"{OUTPUT_DIR}/pairwise_similarity_within_treatment.csv"
df_within.to_csv(within_similarity_filename, index=False)
print(f"Extracted {len(df_within)} within-treatment pairs.")
print(f"Saved within-treatment similarities to: {within_similarity_filename}")


# ### 4.2. Helper Functions: Statistics for Permutation Tests

# Statistic: Difference in mean within-treatment similarity (Combination vs Single)
def calculate_stat_single_vs_comb_within_sim(df_profiles):
    data_sim = get_pair_similarity(df_profiles)
    df_w = data_sim.query('within_between=="w"')
    if df_w.empty: return 0 # Handle empty case
    # Group by n1 (number of components in treatment1) and calculate mean similarity
    sim_means = df_w.groupby('n1')['similarity'].mean()
    mean_comb = sim_means.get(2, 0) # Default to 0 if no combination treatments
    mean_single = sim_means.get(1, 0) # Default to 0 if no single treatments
    return mean_comb - mean_single

# Statistic: Kruskal-Wallis H-statistic for within-treatment similarities across groups
def calculate_stat_kruskal_within_sim(df_profiles):
    data_sim = get_pair_similarity(df_profiles)
    df_w = data_sim.query('within_between=="w"')
    if df_w.empty or df_w['treatment1'].nunique() < 2: return 0 # Need at least 2 groups
    # Prepare list of similarity values for each treatment group
    sims_by_treatment = [group['similarity'].values
                         for name, group in df_w.groupby('treatment1')]
    if len(sims_by_treatment) < 2: return 0 # Need at least 2 groups for KW test
    # Check if all groups have data
    sims_by_treatment = [s for s in sims_by_treatment if len(s) > 0]
    if len(sims_by_treatment) < 2: return 0

    h_stat, _ = kruskal(*sims_by_treatment)
    return h_stat if not np.isnan(h_stat) else 0 # Handle potential NaN result

# Statistic: Difference in mean similarity (Within vs Between treatments)
def calculate_stat_within_vs_between_sim(df_profiles):
    data_sim = get_pair_similarity(df_profiles)
    if data_sim.empty: return 0
    sim_means = data_sim.groupby('within_between')['similarity'].mean()
    mean_within = sim_means.get('w', 0)
    mean_between = sim_means.get('b', 0)
    return mean_within - mean_between

# Statistic: Difference in mean within-treatment similarity (Specific Single AMP vs Others)
def calculate_stat_specific_amp_vs_others_within_sim(df_profiles, target_amp):
    data_sim = get_pair_similarity(df_profiles)
    df_w = data_sim.query('within_between=="w"')
    if df_w.empty: return 0

    amp_mean = df_w.query('treatment1 == @target_amp')['similarity'].mean()
    # Ensure not_amp_mean calculation avoids division by zero if only target_amp exists
    other_sims = df_w.query('treatment1 != @target_amp')['similarity']
    not_amp_mean = other_sims.mean() if not other_sims.empty else 0

    # Handle cases where one group might be missing -> return 0 difference
    if pd.isna(amp_mean) or pd.isna(not_amp_mean):
         return 0
    return amp_mean - not_amp_mean

# Statistic: Diff between mean sim WITHIN a comb-treatment vs. mean sim BETWEEN that comb and a specific single-AMP treatment
def calculate_stat_comb_within_vs_comb_between_single(df_profiles, amp_comb, amp_single):
    data_sim = get_pair_similarity(df_profiles)
    if data_sim.empty: return 0

    # Similarity WITHIN the combination treatment
    within_comb_sims = data_sim.query('within_between == "w" and treatment1 == @amp_comb')['similarity']
    mean_within_comb = within_comb_sims.mean() if not within_comb_sims.empty else 0

    # Similarity BETWEEN the combination and the specific single AMP treatment
    between_comb_single_sims = data_sim.query(
        'within_between == "b" and '+
        f'((treatment1 == "{amp_comb}" and treatment2 == "{amp_single}") or ' +
        f'(treatment1 == "{amp_single}" and treatment2 == "{amp_comb}"))'
        )['similarity']
    mean_between_comb_single = between_comb_single_sims.mean() if not between_comb_single_sims.empty else 0

    if pd.isna(mean_within_comb) or pd.isna(mean_between_comb_single):
        return 0
    return mean_within_comb - mean_between_comb_single


# ### 4.3. Helper Function: Run Permutation Test

def run_permutation_test(df_profiles, statistic_func, num_permutations=NUM_PERMUTATIONS, **kwargs):
    """
    Runs a permutation test by shuffling the index (strain labels).

    Args:
        df_profiles (pd.DataFrame): The mutation profile data (rows=strains).
        statistic_func (callable): Function to calculate the test statistic.
                                  It should accept df_profiles and potentially **kwargs.
        num_permutations (int): Number of permutations to run.
        **kwargs: Additional keyword arguments passed to statistic_func.

    Returns:
        tuple: (observed_statistic, permuted_statistics, p_value)
    """
    print(f"  Running permutation test ({num_permutations} permutations)...")
    observed_statistic = statistic_func(df_profiles, **kwargs)

    permuted_statistics = np.zeros(num_permutations)
    original_index = df_profiles.index.copy()

    for i in range(num_permutations):
        if (i+1) % (num_permutations // 10) == 0: # Print progress
             print(f"    Permutation {i+1}/{num_permutations}")
        shuffled_df = df_profiles.copy()
        # Shuffle index (breaks link between profile and original treatment label)
        permuted_index = np.random.permutation(original_index)
        shuffled_df.index = permuted_index
        # The 'treatment' column derived from the *original* index remains,
        # but it's now associated with a shuffled *profile*. We need to
        # re-derive the treatment column based on the *new* (shuffled) index
        # for the statistic function to work correctly.
        # However, the get_pair_similarity function internally derives treatments
        # from the index it receives. So, just shuffling the index of the input
        # dataframe IS the correct way to shuffle the profile-treatment link.
        permuted_statistics[i] = statistic_func(shuffled_df, **kwargs)

    # Calculate two-tailed p-value
    p_value = np.mean(np.abs(permuted_statistics) >= np.abs(observed_statistic))

    print(f"  Observed statistic: {observed_statistic:.4f}")
    print(f"  Permutation P-value: {p_value:.4f}")
    return observed_statistic, permuted_statistics, p_value

def plot_permutation_histogram(observed_statistic, permuted_statistics, title, filename):
    """Plots histogram of permuted statistics with observed value lines."""
    plt.figure(figsize=(6, 4))
    plt.hist(permuted_statistics, bins=50, alpha=0.7, label='Permuted Statistics')
    ylim = plt.ylim()
    plt.plot([observed_statistic, observed_statistic], ylim, 'r-', lw=2, label='Observed Statistic')
    plt.plot([-observed_statistic, -observed_statistic], ylim, 'r--', lw=2, label='-Observed Statistic')
    plt.ylim(ylim)
    plt.xlabel("Statistic Value")
    plt.ylabel("Frequency")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    print(f"Saved permutation histogram: {filename}")
    plt.close()


# ### 4.4. Permutation Test: Within-Treatment Similarity (Single vs. Combination)

print("\nPermutation Test: Within-Treatment Similarity (Combination vs Single)")
obs_s_vs_c, perm_s_vs_c, p_s_vs_c = run_permutation_test(
    df_dor_similarity,
    calculate_stat_single_vs_comb_within_sim
)
hist_filename = f"{OUTPUT_DIR}/perm_hist_sim_single_vs_comb.png"
plot_permutation_histogram(obs_s_vs_c, perm_s_vs_c, "Permutation Dist: Sim (Comb - Single)", hist_filename)


# ### 4.5. Visualize Within-Treatment Similarity (Single vs. Combination)

# Using the pre-calculated `df_within` DataFrame
plt.figure(figsize=(10, 6))
# Boxplot base
sns.boxplot(data=df_within, y='similarity', x='n1', color=COLOR_PALETTE['boxplot_within'],
            width=0.7, showfliers=False)
# Swarmplot overlay
sns.swarmplot(data=df_within, y='similarity', x='n1', hue='treatment',
              edgecolor='black', linewidth=0.5, palette=COLOR_PALETTE, legend=True)

plt.xlabel('Treatment Type', fontsize=16)
plt.ylabel('Within-Treatment Similarity (Dice)', fontsize=16)
plt.xticks([0, 1], ['Single', 'Combination'], fontsize=14)
plt.yticks(fontsize=14)
plt.legend(title='Treatment', loc='center left', bbox_to_anchor=(1.05, 0.5), edgecolor='white')
plt.tight_layout(rect=[0, 0, 0.85, 1]) # Adjust layout for legend

plot_filename_sim_s_vs_c = f"{OUTPUT_DIR}/similarity_within_single_vs_comb_boxplot_swarm.png"
plt.savefig(plot_filename_sim_s_vs_c, dpi=300)
print(f"\nSaved plot: {plot_filename_sim_s_vs_c}")
plt.close()


# ### 4.6. Permutation Test: Equality of Within-Treatment Similarities Across Treatments

print("\nPermutation Test: Equality of Within-Treatment Similarities (Kruskal-Wallis H)")
obs_kw, perm_kw, p_kw = run_permutation_test(
    df_dor_similarity,
    calculate_stat_kruskal_within_sim
)
hist_filename_kw = f"{OUTPUT_DIR}/perm_hist_sim_kruskal_within.png"
plot_permutation_histogram(obs_kw, perm_kw, "Permutation Dist: Sim (Kruskal-Wallis H)", hist_filename_kw)


# ### 4.7. Visualize Within-Treatment Similarity (Per Treatment Group)

# Boxplot with Swarmplot
plt.figure(figsize=(12, 7)) # Adjusted size
# Boxplot base
sns.boxplot(data=df_within, y='similarity', x='treatment', color=COLOR_PALETTE['boxplot_within'],
            order=TREATMENTS, showfliers=False) # Use consistent order
# Swarmplot overlay
sns.swarmplot(data=df_within, y='similarity', x='treatment', hue='treatment',
              edgecolor='black', linewidth=0.5, palette=COLOR_PALETTE, legend=False, # Legend redundant here
              order=TREATMENTS) # Use consistent order

plt.xlabel('Treatment', fontsize=18)
plt.ylabel('Within-Treatment Similarity (Dice)', fontsize=18)
plt.xticks(rotation=45, ha='right', fontsize=14)
plt.yticks(fontsize=14)
plt.tight_layout()

plot_filename_sim_per_treat = f"{OUTPUT_DIR}/similarity_within_per_treatment_boxplot_swarm.png"
plt.savefig(plot_filename_sim_per_treat, dpi=300)
print(f"\nSaved plot: {plot_filename_sim_per_treat}")
plt.close()

# Heatmaps per treatment group
print("\nGenerating within-treatment similarity heatmaps...")
fig, axes = plt.subplots(2, 3, figsize=(24, 16), sharey=True, sharex=True) # Shared axes might be tricky if strain names differ
axes = axes.flatten()

# Custom colormap from original script
colors_heatmap = ['#CDF599','#007098'] # Light green to dark blue
cmap_name = 'custom_div'
custom_cmap = LinearSegmentedColormap.from_list(cmap_name, colors_heatmap, N=256)

# Determine global vmin/vmax across all within-treatment similarities
sim_min = df_within['similarity'].min()
sim_max = df_within['similarity'].max()

cbar_ax = fig.add_axes([.91, .3, .03, .4]) # Position for color bar

for i, treatment in enumerate(TREATMENTS):
    ax = axes[i]
    # Filter data for the current treatment
    data_treat = df_within[df_within['treatment'] == treatment]

    if not data_treat.empty:
        # Pivot requires unique index/column pairs. Strain pairs are unique in df_within.
        # We need to create a symmetric matrix.
        strains_in_treat = pd.unique(data_treat[['strain1', 'strain2']].values.ravel('K'))
        matrix = pd.DataFrame(np.nan, index=strains_in_treat, columns=strains_in_treat)
        # Fill from pairwise data
        for _, row in data_treat.iterrows():
            matrix.loc[row['strain1'], row['strain2']] = row['similarity']
            matrix.loc[row['strain2'], row['strain1']] = row['similarity']
        # Fill diagonal with 1 (similarity of a strain with itself)
        np.fill_diagonal(matrix.values, 1.0)

        # Create heatmap
        sns.heatmap(matrix, ax=ax, cbar=i == 0, cbar_ax=None if i else cbar_ax,
                    vmin=0, vmax=1, # Fixed scale 0-1 for Dice similarity
                    yticklabels=True, xticklabels=True,
                    cmap=custom_cmap, square=True) # Use the custom colormap

        ax.set_title(treatment, fontsize = 24, pad=10)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha='right')
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0, ha='right')
        ax.set_xlabel('')
        ax.set_ylabel('')
    else:
        ax.set_title(f"{treatment}\n(No within-pairs)", fontsize=20)
        ax.axis('off') # Hide empty axes

# Configure color bar
cbar_ax.set_ylabel('Similarity (Dice)', rotation=-90, va="bottom", fontsize=16)
cbar_ax.tick_params(labelsize=14)

plt.tight_layout(h_pad= 2, rect=[0, 0, .9, 1]) # Adjust layout for colorbar
heatmap_filename = f"{OUTPUT_DIR}/similarities_within_treatments_heatmaps.png"
plt.savefig(heatmap_filename, dpi=300)
print(f"Saved heatmap figure: {heatmap_filename}")
plt.close()

# Calculate and save average/median within-treatment similarity stats
print("\nCalculating average/median within-treatment similarities:")
similarity_stats_within = df_within.groupby('treatment')['similarity'].agg(['mean', 'median'])
similarity_stats_within.columns = ['Average Similarity', 'Median Similarity']
print(similarity_stats_within)
stats_filename = f"{OUTPUT_DIR}/similarity_stats_within_treatments.csv"
similarity_stats_within.to_csv(stats_filename)
print(f"Saved within-treatment similarity stats to: {stats_filename}")


# ### 4.8. Permutation Test: Within vs. Between Treatment Similarity (Overall)

print("\nPermutation Test: Overall Within vs. Between Treatment Similarity")
obs_w_vs_b, perm_w_vs_b, p_w_vs_b = run_permutation_test(
    df_dor_similarity,
    calculate_stat_within_vs_between_sim
)
hist_filename_w_vs_b = f"{OUTPUT_DIR}/perm_hist_sim_within_vs_between.png"
plot_permutation_histogram(obs_w_vs_b, perm_w_vs_b, "Permutation Dist: Sim (Within - Between)", hist_filename_w_vs_b)


# ### 4.9. Permutation Tests: Specific Treatment Comparisons

# Test 1: Within-similarity of specific single AMP vs Others
print("\nPermutation Tests: Within-Similarity (Single AMP vs Others)")
for amp in SINGLE_AMPS:
    print(f"\n-- Testing AMP: {amp} --")
    obs_amp, perm_amp, p_amp = run_permutation_test(
        df_dor_similarity,
        calculate_stat_specific_amp_vs_others_within_sim,
        target_amp=amp
    )
    hist_filename_amp = f"{OUTPUT_DIR}/perm_hist_sim_within_{amp}_vs_others.png"
    plot_permutation_histogram(obs_amp, perm_amp, f"Permutation Dist: Sim ({amp} vs Others)", hist_filename_amp)


# Test 2: Within-similarity of Combination vs. Between-similarity (Comb vs Single Component)
print("\nPermutation Tests: Within-Comb vs. Between-(Comb-Single)")
for amp_comb in COMB_AMPS:
    # Get components of the combination
    components = amp_comb.split('_')
    for amp_single in components: # Test against each component
        if amp_single in SINGLE_AMPS: # Ensure component is a valid single AMP treatment
             print(f"\n-- Testing: {amp_comb} (Within) vs. {amp_comb}-{amp_single} (Between) --")
             obs_comb, perm_comb, p_comb = run_permutation_test(
                 df_dor_similarity,
                 calculate_stat_comb_within_vs_comb_between_single,
                 amp_comb=amp_comb,
                 amp_single=amp_single
             )
             hist_filename_comb = f"{OUTPUT_DIR}/perm_hist_sim_{amp_comb}_vs_{amp_single}.png"
             plot_permutation_histogram(obs_comb, perm_comb,
                                        f"Permutation Dist:\nSim (Within {amp_comb}) - (Between {amp_comb}-{amp_single})",
                                        hist_filename_comb)


# ### 4.10. Between-Treatment Similarity Matrix and Clustermap

print("\nGenerating between-treatment similarity clustermap...")

# Calculate average similarity between each pair of treatments
# Use the full similarity data (within and between pairs)
avg_similarity_matrix = data_similarity_all_pairs.groupby(['treatment1', 'treatment2'])['similarity'].mean().unstack()

# Fill NaN and ensure symmetry (Avg(A,B) should ideally be same as Avg(B,A))
# This fills based on the assumption that if (A,B) exists, (B,A) should be the same.
# If only one exists, it uses that value. If neither exists, remains NaN (though unlikely with combinations).
avg_similarity_matrix = avg_similarity_matrix.combine_first(avg_similarity_matrix.T)

# Fill diagonal with within-treatment averages (calculated earlier)
within_avg_sim = similarity_stats_within['Average Similarity']
for treat in avg_similarity_matrix.index:
    if treat in within_avg_sim.index:
        avg_similarity_matrix.loc[treat, treat] = within_avg_sim[treat]
    else:
        avg_similarity_matrix.loc[treat, treat] = np.nan # Or 1.0 if preferred diagonal

# Reorder matrix columns/rows consistently based on TREATMENTS list for interpretability
ordered_treatments = [t for t in TREATMENTS if t in avg_similarity_matrix.index] # Keep only treatments present
avg_similarity_matrix = avg_similarity_matrix.loc[ordered_treatments, ordered_treatments]

print("\nAverage Similarity Matrix (Between/Within Treatments):")
print(avg_similarity_matrix)
matrix_filename = f"{OUTPUT_DIR}/average_similarity_matrix_treatments.csv"
avg_similarity_matrix.to_csv(matrix_filename)
print(f"Saved average similarity matrix to: {matrix_filename}")

# Create the clustermap
# Define custom colormap from original script
colors_clustermap = ['#D0F5ED', '#3B62B5'] # Light teal to dark blue
cmap_name_cluster = 'custom_cluster'
custom_cmap_cluster = LinearSegmentedColormap.from_list(cmap_name_cluster, colors_clustermap, N=256)

try:
    g = sns.clustermap(avg_similarity_matrix.fillna(0), # Fill NaNs for clustering, e.g., with 0
                       cmap=custom_cmap_cluster,
                       annot=True, fmt=".3f",
                       linewidths=.5,
                       col_cluster=True, row_cluster=True, # Cluster both rows and columns
                       cbar_pos=None, # Remove default color bar if adding manually or not needed
                       figsize=(12, 12))

    # Adjust plot aesthetics
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right')
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.ax_heatmap.set_xlabel("Treatment", fontsize=14)
    g.ax_heatmap.set_ylabel("Treatment", fontsize=14)
    g.fig.suptitle('Clustermap of Average Treatment Similarities (Dice)', y=1.02, fontsize=18)


    clustermap_filename = f"{OUTPUT_DIR}/clustermap_average_treatment_similarities.png"
    plt.savefig(clustermap_filename, dpi=300, bbox_inches='tight')
    print(f"\nSaved clustermap: {clustermap_filename}")
    plt.close(g.fig) # Close the figure associated with the clustermap

except Exception as e:
    print(f"\nCould not generate clustermap. Error: {e}")
    print("This might happen if the similarity matrix has too many NaNs or lacks variance.")


print("\n--- Mutation Profile Similarity Analysis Complete ---")
print("\n--- Script Finished ---")