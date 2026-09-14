# fdgpet-dwmri

## Overview
Codebase for "Diagnostic performance of quantitative measures from [18F]FDG PET/CT, [18F]FEC PET/CT, and DW-MRI in the detection of lymph node metastases in endometrial and cervical cancer: data from the MAPPING study" (doi: [10.1007/s00259-025-07587-3](https://doi.org/10.1007/s00259-025-07587-3)).

## Fetching data
Data available upon reasonable request to corresponding author.
Create a `.env` file in the project root with the Windows path to the data folder, e.g:

```
DATA_PATH="C:\path\to\data"
```

Then run `./fetch.sh` to sync it into `./data/`.
