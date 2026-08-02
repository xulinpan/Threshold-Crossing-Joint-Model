# Peer-Reviewed Journal Recommendations for the GLW CML Paper

Prepared: 2026-07-09

Basis: recommendations use the current academic evaluation of the paper plus current journal scope pages checked on 2026-07-09. The paper is strongest as an applied biostatistics and real-world CML molecular monitoring workflow, not as a purely clinical discovery paper.

## Recommended Submission Strategy

### First-choice target: BMC Medical Research Methodology

**Recommendation:** strongest fit if the manuscript emphasizes the statistical and data-structure contribution.

**Why it fits**

- The journal publishes methodological approaches to healthcare research, including data analysis, statistics, and modeling.
- The GLW paper's best contribution is an interval-aware joint modeling framework for real-world clinical monitoring data.
- The paper can be framed as a solution to a common healthcare research problem: repeated biomarkers, irregular visits, interval-observed clinical endpoints, and incomplete covariates.

**Best manuscript framing**

> A reproducible joint longitudinal-interval modeling framework for deep molecular response under irregular real-world CML monitoring.

**What to strengthen before submission**

- Run the full Stan model with posterior summaries.
- Add simulation or resampling checks showing why exact-time analysis is biased or less appropriate when DMR is interval-observed.
- Include sensitivity analyses for assay-floor handling and visit-gap terms.

**Risk**

- If the paper remains mostly descriptive EDA plus a proposed model, reviewers may ask for a fuller methodological evaluation.

## Strong Alternative Targets

### 2. BMC Bioinformatics

**Recommendation:** good fit if the manuscript is framed as computational/statistical modeling of biological data.

**Why it fits**

- The journal considers computational models, tools, statistical methods, machine learning, and AI for modeling and analysis of biological data.
- The GLW workflow includes cleaned model-ready biological/clinical molecular monitoring data, Stan model code, R benchmark scripts, and reproducible outputs.

**Best manuscript framing**

> A reproducible computational workflow and joint model for longitudinal BCR-ABL1/log-MRD monitoring in CML.

**What to strengthen before submission**

- Emphasize code availability and reproducibility.
- Present the model as a reusable workflow, not only a single-cohort analysis.
- Add a clear benchmark comparison: Kaplan-Meier/exact-time model vs interval-censored Weibull vs joint model.

**Risk**

- Reviewers may expect a more generalizable algorithm, software package, or multi-dataset validation.

### 3. BMC Medical Informatics and Decision Making

**Recommendation:** good fit if the paper is positioned around clinical monitoring, decision support, and real-world information systems.

**Why it fits**

- The journal covers health information technologies, decision-making, healthcare information systems, machine learning, and modeling.
- The GLW paper can be framed as converting routine molecular monitoring records into a decision-ready response model.

**Best manuscript framing**

> A privacy-preserving real-world monitoring workflow for estimating deep molecular response from irregular CML follow-up data.

**What to strengthen before submission**

- Connect model outputs to clinical decision-making, such as identifying patients likely to reach DMR or become candidates for TKI discontinuation discussions.
- Add an applied decision-use section with predicted interval DMR probabilities.

**Risk**

- If the manuscript stays purely statistical and does not connect to decision support, the fit is weaker than BMC Medical Research Methodology.

### 4. PLOS ONE

**Recommendation:** strong practical fallback if the paper is methodologically sound but not positioned as high novelty.

**Why it fits**

- PLOS ONE has a broad scope and evaluates manuscripts based on scientific validity, strong methodology, and ethics rather than perceived significance.
- The paper's novelty is moderate but the reproducible workflow, EDA, data cleaning, references, and model specification can form a sound applied study.

**Best manuscript framing**

> Reproducible modeling of interval-observed deep molecular response in real-world chronic myeloid leukemia monitoring.

**What to strengthen before submission**

- Make the paper complete and transparent: full code, model output, sensitivity analyses, and clear data privacy handling.
- Avoid overclaiming clinical discovery.

**Risk**

- Reviewers may still ask whether the sample size and lack of external validation support the conclusions.

### 5. Scientific Reports

**Recommendation:** possible if full model results are strong and the paper is written as a rigorous applied medical-data modeling study.

**Why it fits**

- Scientific Reports publishes original research across medicine, biomedical and clinical sciences, and engineering.
- The journal has broad reach and can accommodate interdisciplinary applied modeling papers.

**Best manuscript framing**

> Real-world longitudinal molecular monitoring reveals interval-observed DMR dynamics in CML.

**What to strengthen before submission**

- Complete the Bayesian model fitting.
- Add clinically interpretable predicted DMR probability curves.
- Add validation through bootstrap, cross-validation, or an external/temporal holdout if feasible.

**Risk**

- The current single-cohort dataset may be viewed as too narrow unless the analysis is very complete and the contribution is clearly original.

## Oncology-Specific Targets

### 6. Frontiers in Oncology, Hematologic Malignancies section

**Recommendation:** reasonable oncology-specific target if the clinical hematology interpretation is strengthened.

**Why it fits**

- Frontiers in Oncology has a Hematologic Malignancies section and welcomes clinical research across cancer domains.
- The GLW data are a real CML cohort rather than public-database-only analysis, which helps with validation expectations.

**Best manuscript framing**

> Interval-observed deep molecular response dynamics in real-world CML monitoring.

**What to strengthen before submission**

- Add a stronger clinical hematology discussion: DMR, MR4.5, CMR, TKI discontinuation eligibility, and practical monitoring implications.
- Include complete model results and sensitivity analyses.

**Risk**

- If framed mainly as computational analysis without enough clinical validation, it may be considered less suitable.

### 7. BMC Cancer

**Recommendation:** possible but not the best first choice.

**Why it fits**

- BMC Cancer covers cancer diagnosis, treatment, clinical research, biomarkers, and computational biology.
- The CML topic is clinically relevant, and the dataset is real rather than a public-database-only analysis.

**Why it is weaker**

- BMC Cancer explicitly cautions that diagnostic/prognostic marker studies should be validated and that purely computational analyses without validation may not be considered.
- The current paper is stronger as a methods/workflow paper than as a validated cancer biomarker study.

**Best use**

- Submit here only if the paper has complete model results and a clinically oriented validation/sensitivity section.

### 8. Cancers

**Recommendation:** optional oncology open-access target if the goal is a broader oncology readership and fast applied dissemination.

**Why it fits**

- The journal publishes oncology research and asks for detailed, reproducible results.
- The GLW workflow is reproducible and clinically relevant to CML monitoring.

**Risk**

- The paper should be clinically complete and transparent about the single-cohort limitation.

## Aspirational Methodology Targets

### 9. Statistical Methods in Medical Research

**Recommendation:** aspirational, not a first submission in the current state.

**Fit**

- Strong only if the manuscript develops or evaluates a genuinely general statistical method.

**Required strengthening**

- Add simulation studies.
- Compare against existing joint modeling and interval-censoring approaches.
- Demonstrate method performance under irregular monitoring, assay-floor censoring, and small samples.

### 10. Statistics in Medicine

**Recommendation:** aspirational, similar to Statistical Methods in Medical Research.

**Fit**

- Suitable only after a substantial methods contribution beyond applying an existing joint model to one cohort.

## Ranked Recommendation

| Rank | Journal | Fit | Best angle | Submit now? |
|---:|---|---|---|---|
| 1 | BMC Medical Research Methodology | Very strong | Healthcare methods, statistics, interval-aware modeling | After full model fitting and sensitivity analyses |
| 2 | BMC Bioinformatics | Strong | Computational/statistical workflow for molecular monitoring data | After emphasizing reusable code/workflow |
| 3 | BMC Medical Informatics and Decision Making | Strong | Real-world monitoring and decision-ready modeling | After adding clinical decision-use framing |
| 4 | PLOS ONE | Strong fallback | Scientifically valid reproducible applied study | Yes, if conclusions are conservative |
| 5 | Scientific Reports | Moderate to strong | Broad biomedical modeling study | After stronger validation/results |
| 6 | Frontiers in Oncology | Moderate | Hematologic malignancy clinical application | After strengthening clinical interpretation |
| 7 | BMC Cancer | Moderate | Cancer clinical/biomarker application | Only with stronger validation |
| 8 | Cancers | Moderate | Oncology applied data/modeling paper | Optional |
| 9 | Statistical Methods in Medical Research | Stretch | New statistical method | Only after simulation/method expansion |
| 10 | Statistics in Medicine | Stretch | Medical statistics method paper | Only after simulation/method expansion |

## Final Recommendation

The best first target is **BMC Medical Research Methodology** if the authors complete the Bayesian joint model fitting and add sensitivity/simulation evidence. The best practical fallback is **PLOS ONE** because the study is reproducible, scientifically valid, and moderate in novelty. If the paper is rewritten for computational biology with reusable code and model workflow emphasized, **BMC Bioinformatics** becomes a strong alternative. If the goal is a hematology/oncology audience, **Frontiers in Oncology** is preferable to BMC Cancer because the current manuscript is more applied-modeling than biomarker-validation.

## Source Links Checked

- BMC Medical Research Methodology: https://bmcmedresmethodol.biomedcentral.com/about
- BMC Bioinformatics: https://bmcbioinformatics.biomedcentral.com/about
- BMC Medical Informatics and Decision Making: https://bmcmedinformdecismak.biomedcentral.com/about
- PLOS ONE journal information: https://journals.plos.org/plosone/s/journal-information
- Scientific Reports about page: https://www.nature.com/srep/about
- Frontiers in Oncology about page: https://www.frontiersin.org/journals/oncology/about
- BMC Cancer: https://bmccancer.biomedcentral.com/about
- Cancers aims and scope: https://www.mdpi.com/journal/cancers/about

