# HE-WSI-lungTILs
Computerized tumor-infiltrating lymphocytes density score predicts survival of patients with resectable lung adenocarcinoma: a multi-institutional study


Codes:
The R code:https://github.com/panxipeng/HE-WSI-lungTILs
nuclei segmentation and classification: https://github.com/vqdang/hover_net
Tissue segmentation: https://github.com/ChuHan89/WSSS-Tissue
HistoQC: https://github.com/choosehappy/HistoQC

Datasets:
Camelyon Challenge dataset: https://camelyon17.grand-challenge.org
Lung adenocarcinoma whole slide images dataset: https://portal.gdc.cancer.gov/
Lung adenocarcinoma RNA transcripts dataset: https://portal.gdc.cancer.gov/

Abstract
Background: The presence of a high abundance of tumor-infiltrating lymphocytes (TILs) has a positive impact on the prognosis of patients with lung adenocarcinoma (LUAD). This study aimed to develop and validate an artificial intelligence-driven pathological scoring system for the assessment of TILs on hematoxylin and eosin (H&E)-stained whole-slide images (WSIs) of patients with resectable LUAD. 
Methods: Deep learning-based methods were applied to identify cancer epithelium, cancer stroma, and lymphocytes. The densities of lymphocytes in cancer epithelium (DLCE) and cancer stroma (DLCS) were calculated; a risk score (named WELL score) was built through linear WEighting of DLCE score and DLCS score. Association between WELL score and patient outcome was explored in 793 patients with stage I-III LUAD in four cohorts (D1, n=276; V1, n=139; V2, n=115; V3, n=263). A discovery cohort was used to optimize WELL scores in terms of overall survival (OS) and disease-free survival (DFS) independently. The median value of WELL scores in the discovery cohort was selected as the cut-off to stratify patients into high- and low-risk groups. Univariable and multivariable Cox regression analysis were performed to evaluate the association of WELL-score and clinical risk factors with OS and DFS. A prediction model (full model) integrating the WELL score and clinical risk factors was constructed, and the performance of which was then evaluated and compared with the clinical model. CIBERSORT analysis based on RNA transcripts was performed to explore the lymphocyte infiltration landscape distinction between high-risk and low-risk groups determined by the WELL score.
Findings: WELL score was an independent prognostic factor in terms of OS in D1 (hazard ratio [HR], 1.95; 95% confidence interval [CI], 1.12???3.41; p=0.019) and V1 (HR, 3.65; 95% CI, 1.61???8.27; p=0.002), and similar trend could also be found in V2 and V3. WELL score was also an independent prognostic factor in terms of DFS in D1 (HR, 1.68; 95% CI, 1.14???2.47; p=0.008) and V2 (HR, 3.66; 95% CI, 1.53???8.76; p=0.004), respectively. The full model demonstrated better discrimination performance than the clinical model in four cohorts (C-index, D1, 0.751 vs. 0.723; V1, 0.622 vs. 0.596; V2, 0.683 vs. 0.650; V3, 0.596 vs. 0.586). The fractions of plasma cells, CD8+ T cells, and activated CD4+ memory T cells in the low-risk group were statistically significantly higher than those in the high-risk group (p<0.05).
Interpretation: This study developed and validated an artificial intelligence-based WSIs analysis workflow and TILs scoring system that could support clinical decision-making and promote stratification for patients with resectable LUAD.
