# coxkl: Integrating External Risk Models with Time-to-Event Data


**Key Features:**

*   **Integrates External Risk Scores:** Seamlessly incorporate predictions from *any* external risk model.
*   **Handles Population Heterogeneity:** Accounts for differences between the population used to train the external model and your study population.
*   **Efficient Computation:** Designed to be computationally efficient, even with high-dimensional data.
*   **Improved Accuracy:**  Demonstrated to improve estimation efficiency and prediction accuracy in survival models.
*   **Cross-Validation:** Includes built-in cross-validation for optimal tuning of the data integration parameter.


# Installation

**Note:** *This package is still in its early stages of development, so please don't hesitate to report any problems you may experience.* 

The package only works for R 4.1.0+.

You can install 'coxkl' via CRAN or github:

    install.packages("coxkl")

    #or
    require("devtools")
    require("remotes")
    remotes::install_github("UM-KevinHe/coxkl")

We recommand to start with <a href="https://um-kevinhe.github.io/coxkl/articles/coxkl.html#quick-start" target="_blank">tutorial</a>, as it provides an overview of the package's usage, including preprocessing, model training, selection of penalization parameters, and post-estimation procedures.


## Detailed tutorial

For detailed tutorial and model paramter explaination, please go to <a href="https://um-kevinhe.github.io/coxkl/index.html" target="_blank">here</a>.

## Getting Help:

If you encounter any problems or bugs, please contact us at: [lfluo\@umich.edu](mailto:lfluo@umich.edu){.email},
[xhliuu\@umich.edu](mailto:xhliuu@umich.edu){.email},
[kevinhe\@umich.edu](mailto:kevinhe@umich.edu){.email}.

## Contributing

We welcome contributions to the **coxkl** package. Please see our [CONTRIBUTING.md](https://github.com/UM-KevinHe/coxkl/blob/main/.github/CONTRIBUTING.md) file for detailed guidelines of how to contribute.

# References

  \[1\] Wang, D., Ye, W., Zhu, J., Xu, G., Tang, W., Zawistowski, M., Fritsche, L. G., & He, K. (2023). Incorporating external risk information with the Cox model under population heterogeneity: Applications to trans-ancestry polygenic hazard scores. *arXiv preprint arXiv:2302.11123*.
  
  \[2\] Luo, L., Taylor, J. M. G., Wang, D., & He, K. (2024). Flexible Deep Learning Techniques for Cox Models with Data Integration
  
  
