# paper-prot-atlas-gbm
Activation level of 64 gene signatures in glioblastoma

This repository includes the XGBoost models and the R scripts to infer the activation status of 64 gene signatures in GBM, as published in this manuscript:

**Topographic mapping of the glioblastoma proteome reveals a triple axis model of intra-tumoral heterogeneity** <br/>
*Lam KHB, Djuric U, Leon AJ, Hui W, Lee SCE, Batruch I, Faust K, Koritzinsky M, Richer M, Diamandis P* **(under review)**

Execute the following command to analyze the sample data:
```
Rscript predictXGB.R --data sample/sample_prot.rds \
                     --models data/list_xgb-gbm_64_signatures-prot-v01.RData \
                     --out_csv predictXGB_sample_rna.csv
```

These models are also part of our online tool: (https://cancerhub.shinyapps.io/prot-atlas-gbm/)




