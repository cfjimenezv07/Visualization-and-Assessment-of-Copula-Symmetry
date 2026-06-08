<br />
<div align="center">
  <a href="https://github.com/cfjimenezv07/Visualization-and-Assessment-of-Copula-Symmetry">
    <img src="UoY.png" alt="York Logo" height="150">
    <img src="Kaustlogo.png" alt="KAUST Logo" height="150">
  </a>

<h3 align="center">Visualization and Assessment of Copula Symmetry</h3>
</div>

## Abstract
<p align="justify">
Visualization and assessment of copula structures are crucial for accurately understanding and modeling the dependencies in multivariate data analysis. In this article, we introduce an innovative method that employs functional boxplots and rank-based testing procedures to evaluate copula symmetries. This approach is specifically designed to assess key characteristics such as reflection symmetry, radial symmetry, and joint symmetry. We first construct test functions for each specific property and then investigate the asymptotic properties of their empirical estimators. We demonstrate that the functional boxplot of these sample test functions serves as an informative visualization tool of a given copula structure, effectively measuring the departure from zero of the test function. Furthermore, we introduce a nonparametric testing procedure to assess the significance of deviations from symmetry, ensuring the accuracy and reliability of our visualization method. Through extensive simulation studies involving various copula models, we demonstrate the effectiveness of our testing approach. Finally, we apply our visualization and testing techniques to three real-world datasets: a nutritional habits survey with five variables, stock price data for the five top companies in the NASDAQ-100 stock index, and two major stock indices, the US S&P500 and German DAX. Supplementary materials for this article are available online.
</p>

### Main Results
The R script files in the `R Code` folder should be used in the following order:

#### 1. Setup, Data Construction, and Auxiliary Utilities
* **Aux_sym_test.R**: Contains auxiliary subroutines to simulate symmetric copulas, construct functional data matrices for functional boxplots, and calculate modified band depth test statistics ($W$) for reflection, radial, and joint symmetry.
* **get_data.R**: Formulates the continuous difference functions needed to construct empirical curves for reflection, radial, and joint copula symmetry testing.
* **nutrients.txt**: Raw text dataset containing nutritional intake measurements utilized in the empirical application.
* **Wind_Dataset1.rds**: A serialized R data file containing the wind data used for the directional empirical case study.

#### 2. Main Simulation & Hypothesis Testing Core
* **copula_simulation.R**: Runs Monte Carlo simulation loops to generate synthetic copula distributions and analyze their baseline properties.
* **HT.R**: Implements the main Hypothesis Testing core functions used to formally calculate the copula symmetry test statistics.

#### 3. Empirical Case Studies & Multi-population Tests
* **DA_nutrition.R**: Conducts the empirical data analysis on the nutritional intake dataset to check for multivariate copula symmetry.
* **DA_wind.R**: Executes the empirical data analysis on the wind speed/direction data to evaluate directional dependence structures.
* **Table3.R**: Master script dedicated to computing and replicating the simulation or empirical results shown in Table 3 of the paper.
* **Table4_joint.R**: Computes the joint symmetry results and statistics required to replicate the first portion of Table 4.
* **Table4_radial.R**: Computes the radial symmetry test results and statistics required to replicate the second portion of Table 4.

#### 4. Diagnostic Visualization, Functional Plots, & Replication
* **copula_structure_visualization.R**: Implements the 2D and 3D empirical copula surface plotting routines to visually inspect symmetry deviations.
* **simulate_plots.R**: Simulates different types of copulas (symmetric and asymmetric) and creates diagnostic visualization plots using `ggplot2`.
* **fb_visualization.R**: Generates functional boxplots (fb) using color-coded p-values to visualize and analyze the distributions of the copula surfaces.
* **tfplot.R**: Visualizes the obtained functional test bands, color-coding the central region in red when rejected ($p\text{-value} < \alpha$) or green when not rejected.
* **Figure1.r**: Contains the specific plotting code to reproduce and generate Figure 1 in the published paper.

## How to Cite
If you use this code or methodology in your research, please cite our paper:

> Jiménez-Varón, C. F., Lee, H., Genton, M. G., & Sun, Y. (2025). Visualization and Assessment of Copula Symmetry. *Journal of Computational and Graphical Statistics*, 34(3), 1140–1152. https://doi.org/10.1080/10618600.2024.2432978

**BibTeX:**
```bibtex
@article{jimenezvaron2025visualization,
  author    = {Jim{\'e}nez-Var{\'o}n, Cristian F. and Lee, Hyoungsub and Genton, Marc G. and Sun, Ying},
  title     = {Visualization and Assessment of Copula Symmetry},
  journal   = {Journal of Computational and Graphical Statistics},
  volume    = {34},
  number    = {3},
  pages     = {1140--1152},
  year      = {2025},
  doi       = {10.1080/10618600.2024.2432978},
  url       = {[https://doi.org/10.1080/10618600.2024.2432978](https://doi.org/10.1080/10618600.2024.2432978)}
}
