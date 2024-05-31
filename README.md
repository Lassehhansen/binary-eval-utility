# From Development to Deployment: Addressing The Downstream Utility Effects of Machine Learning Healthcare Performance

This repository contains the code and analysis for the thesis addressing the gap in ML model evaluations by providing a framework that guides models from development through deployment, emphasizing how evaluation metrics impact downstream utility for various subpopulations under subpopulation shifts.

![infographic_thesis](https://github.com/Lassehhansen/binary-eval-utility/assets/54820693/ee0332c5-58f2-4302-9f56-d51e80fd9597)

The thesis explores the perceived superiority of the discrimination metric AUPRC under class imbalance by reviewing a large volume of machine learning papers. A large-scale simulation is conducted to show that AUPRC is not universally preferred in such settings and to provide a framework to understand the downstream utility effects of model evaluation metric choices. This involves generating over 800,000 model scores and examining the sensitivity of AUROC and AUPRC to subpopulation shifts. Subsequently, the impact of calibration on downstream utility under subpopulation shifts is analyzed through a case study on severity of illness prediction using the GOSSIS model, trained on more than 300,000 individuals from the US, Australia, and New Zealand. It demonstrates how miscalibration in subgroups leads to lower utility in different deployment settings.

## Repository Structure

- **decision-curve-analysis**: Contains code to reproduce the analysis for downstream utility effects of choosing AUPRC versus AUROC in deployment settings.
- **expected-cost-analysis**: Contains code for the expected cost analysis, which is not included in the thesis but follows the same process as the decision-curve analysis.
- **figures**: Contains figures from the analysis.
- **gossis-analysis**: Contains analysis of the downstream utility effects of miscalibration on utility downstream. The GOSSIS dataset can be found [here](https://physionet.org/content/gossis-1-eicu/1.0.0/).
- **measure-sensitivity**: Contains the analysis and plots provided for the sensitivity of AUPRC and AUROC to subpopulation shifts.
- **simulation**: Contains code to reproduce the simulations.

## Simulation Details

To simulate risk scores between 0 and 1, the beta distribution is used due to its flexibility and ability to generate model probabilities. The simulation involves generating 847,800 distinct sets of positive and negative disease distributions along with their corresponding labels for both subpopulations A and B. The parameters for each distribution are determined based on expected values and confidence levels. The simulated probabilities and labels from both subpopulations are combined to form a unified dataset, representing the total population. Each combination of parameters represents a distinct model, resulting in a comprehensive exploration of the parameter space.

## Analysis and Findings

- **Decision Curve Analysis (DCA)**: Demonstrates how popular evaluation frameworks like DCA fail to provide an appropriate deployment threshold under subgroup miscalibration.
- **Expected Cost Minimization (ECM)**: Introduces an ECM framework to set optimal thresholds for different subgroups, addressing the shortcomings of DCA and ensuring a more nuanced handling of false positive and false negative costs.

Collectively, these findings demonstrate that the evaluation frameworks set forth by the TRIPOD + AI statement and DCA are insufficient for evaluating individuals at the margins, under imperfect calibration, and in cases of subpopulation shifts. This emphasizes the importance of understanding subgroup structures and the costs and benefits algorithms will have on these structures in decision policies to ensure that algorithms do not unnecessarily harm protected groups.
