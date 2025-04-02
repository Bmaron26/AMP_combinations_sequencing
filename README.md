# Analysis Code for AMP Treatment Effects on Bacterial Mutations and Profile Similarity

## Description

This repository contains the Python code and associated data files used for the analysis presented in the manuscript:

**Uncovering the Genetic Basis of Staphylococcus aureus Resistance to Single Antimicrobial Peptides and Their Combinations**

(Accepted to iScience)

The script `Data_analysis.py` performs statistical analyses and generates plots to investigate the effects of different antimicrobial peptide (AMP) treatments (single vs. combination) on the number of mutations and the similarity of mutation profiles in evolved bacterial populations.

## Citation

**To cite this code/repository:**

Please cite the specific version archived on Zenodo:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15125182.svg)](https://doi.org/10.5281/zenodo.15125182)


**To cite the associated manuscript:**
Maron, B., et al. 2025. iScience (final update will be added after publishing)


## Contents

* `Data_analysis.py`: The main Python script performing all analyses and generating outputs.
* `requirements.txt`: A list of required Python packages and their versions needed to run the script.
* `mutations_per_strain_data_new_paul.xlsx`: Input data file containing mutation counts per strain, treatment type, and single/combination classification.
* `mutations_and_resistance_correlation_data_new_paul.xlsx`: Input data file containing the detailed mutation profiles (presence/absence) for each strain, used for similarity calculations.
* `README.md`: This explanatory file.
* `LICENSE`: Contains the software license information.

## System Requirements

* **Operating System:** Tested on Windows 10/11. Should be compatible with macOS and Linux.
* **Python:** Python 3.9 or higher recommended (tested with Python 3.11).

## Installation

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/Bmaron26/AMP_combinations_sequencing
    cd AMP_combinations_sequencing
    ```

2.  **Create a Virtual Environment (Recommended):**
    It's good practice to use a virtual environment to manage dependencies for a specific project.
    ```bash
    python -m venv venv
    ```
    Activate the environment:
    * On Windows:
        ```bash
        .\venv\Scripts\activate
        ```
    * On macOS/Linux:
        ```bash
        source venv/bin/activate
        ```

3.  **Install Dependencies:**
    Install all required packages using pip and the `requirements.txt` file:
    ```bash
    pip install -r requirements.txt
    ```

## Usage

1.  Ensure you are in the root directory of the cloned repository (e.g., `AMP_combinations_sequencing`).
2.  If you created a virtual environment, make sure it is activated.
3.  Run the main analysis script from your terminal:
    ```bash
    python Data_analysis.py
    ```
4.  The script will print status messages to the console as it loads data, performs analyses (Mann-Whitney U, Kruskal-Wallis, Dunn's test, similarity calculations, permutation tests), and saves output files.

## Output

The script will generate several output files in the same directory where it is run (the repository root). These include:

* **Plots (.png format):** Visualizations of mutation counts and similarities. Key plots include:
    * `mutations_single_vs_combination_boxplot_swarm.png`
    * `mutations_across_treatments_boxplot_strip.png`
    * `similarity_within_single_vs_comb_boxplot_swarm.png`
    * `similarity_within_per_treatment_boxplot_swarm.png`
    * `similarities_within_treatments_heatmaps.png`
    * `clustermap_average_treatment_similarities.png`
    * Permutation test histograms (e.g., `perm_hist_sim_*.png`)

* **Data / Results (.csv format):** Processed data and statistical test results. Key files include:
    * `dunn_test_results_mutations.csv`
    * `pairwise_similarity_all_pairs.csv`
    * `pairwise_similarity_within_treatment.csv`
    * `similarity_stats_within_treatments.csv`
    * `average_similarity_matrix_treatments.csv`


## License

This project is licensed under the 'MIT license'. See the `LICENSE` file for details.

## Contact

For questions about the code or analysis, please contact Bar Maron at bar.maron@mail.huji.ac.il or open an issue in this GitHub repository.
