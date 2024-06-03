# From Development to Deployment: Addressing The Downstream Utility Effects of Machine Learning Healthcare Performance

This repository contains the code and analysis for the thesis addressing the gap in ML model evaluations by providing a framework that guides models from development through deployment, emphasizing how evaluation metrics impact downstream utility for various subpopulations under subpopulation shifts.

![infographic_thesis](https://github.com/Lassehhansen/binary-eval-utility/assets/54820693/ee0332c5-58f2-4302-9f56-d51e80fd9597)

**Figure**: This infographic illustrates the ML development and deployment evaluation cycle as described in the thesis. The process begins with problem definition and data choice, followed by training and tuning models on split data. Models are evaluated based on discriminative accuracy (AUROC/AUPRC) and calibration. The thesis highlights the gap between model development and deployment, emphasizing downstream utility effects of evaluation choices on clinical utility. In deployment, models chosen based on AUPRC can lower utility for minority subgroups if thresholds are not appropriately set, while calibration issues can lead to lowered utility.

**Thesis Objective**

The thesis that this code is produced for, explores the perceived superiority of the discrimination metric AUPRC under class imbalance by reviewing a large volume of machine learning papers. A large-scale simulation is conducted to show that AUPRC is not universally preferred in such settings and to provide a framework to understand the downstream utility effects of model evaluation metric choices. This involves generating over 800,000 model scores and examining the sensitivity of AUROC and AUPRC to subpopulation shifts. 

Subsequently, the impact of calibration on downstream utility under subpopulation shifts is analyzed through a case study on severity of illness prediction using the GOSSIS model, trained on more than 300,000 individuals from the US, Australia, and New Zealand. It demonstrates how miscalibration in subgroups leads to lower utility in different deployment settings.

## Simulation Details

A large-scale simulation of more than 800,000 model scores for a binary classification task was conducted. The simulations aimed to determine whether two popular model evaluation metrics—Area Under the Precision-Recall Curve (AUPRC) and Area Under the Receiver Operating Characteristic Curve (AUROC)—display sensitivities to subpopulation shifts and whether the choice of metric results in differing clinical utility for two subpopulations under model deployment.

To simulate risk scores between 0 and 1, the beta distribution was used due to its flexibility and ability to generate model probabilities. The simulation involved generating 847,800 distinct sets of positive and negative disease distributions along with their corresponding labels for both subpopulations A and B. The parameters for each distribution were determined based on expected values and confidence levels. The simulated probabilities and labels from both subpopulations were combined to form a unified dataset, representing the total population. Each combination of parameters represents a distinct model, resulting in a comprehensive exploration of the parameter space.

- **Parameter Selection**: The expected value for positive cases (EV Positive) ranged between 0.51 and 0.80, and negative cases (EV Negative) between 0.20 and 0.49.
- **Prevalence Levels - Minority/Majority Group (A/B)**: Prevalence levels for subgroups A/B were simulated at 0.05, 0.20, 0.35, and 0.50.
- **Attribute Ratio**: The attribute ratio, reflecting the relative size of subgroup A within the dataset, was varied at 5\%, 15\%, and 25\%.
- **Confidence Parameters**: The beta parameter was simulated at three levels of confidence: 0.10, 2.05, and 4.00.

### Repository Structure for Simulation

- **simulation**: Contains code to reproduce the simulations.
- **measure-sensitivity**: Contains the analysis and plots provided for the sensitivity of AUPRC and AUROC to subpopulation shifts.
- **decision-curve-analysis**: Contains code to reproduce the analysis for downstream utility effects of choosing AUPRC versus AUROC in deployment settings.
- **expected-cost-analysis**: Contains code for the expected cost analysis, which is not included in the thesis but follows the same process as the decision-curve analysis.
- **figures**: Contains figures from the analysis.

## GOSSIS Analysis

The thesis also focuses on understanding the effect of calibration on downstream utility under subpopulation shifts through a case study on severity of illness prediction using the GOSSIS model.

### GOSSIS Analysis Details

The GOSSIS model is trained on more than 300,000 individuals from the US, Australia, and New Zealand. The analysis demonstrates how miscalibration in subgroups leads to lower utility in different deployment settings. The dataset can be downloaded from [PhysioNet](https://physionet.org/content/gossis-1-eicu/1.0.0/).

### Repository Structure for GOSSIS Analysis

- **gossis-analysis**: Contains the analysis of the downstream utility effects of miscalibration on utility downstream.

## Analysis and Findings

The simulations showed that AUPRC displayed sensitivity to subpopulation shifts and that prevalence differences between subgroups undermined performance in deployment settings for minority subpopulations. At the same time, AUROC remained stable across varying subpopulation shifts. Nonetheless, 424 papers claimed AUPRC’s superiority over AUROC without substantial evidence, often citing misinterpreted sources. Furthermore, the deployment analysis revealed the implications of not considering people at the margins in model evaluations from development through deployment. Popular evaluation frameworks like decision curve analysis (DCA), which assumes perfect calibration, were ineffective for setting deployment thresholds when miscalibration from subpopulation shifts occurred, as this undermined fairness. Instead, an expected cost minimization (ECM) framework was introduced to set optimal thresholds for different subgroups to address these shortcomings.

Collectively, these findings demonstrate that the evaluation frameworks set forth by the TRIPOD + AI statement and DCA are insufficient for evaluating individuals at the margins, under imperfect calibration, and in cases of subpopulation shifts. This emphasizes the importance of understanding subgroup structures and the costs and benefits algorithms will have on these structures in decision policies to ensure that algorithms do not unnecessarily harm protected groups.
