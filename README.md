
---

# 🧬 Population Sampling and Representation Pipelines for Human Pangenome Reference Consortium

*This projects aims to extend the human pangenome reference panel improving the coverage of human genetic diversity.*

---

## 📋 Table of Contents

* [Overview](#overview)
* [Requirements](#requirements)
* [Installation](#installation)
* [Usage](#usage)

---

## 🧩 Overview

To extend the HPRC referece, the Population Sampling and Representation working group has designed 2 computational strategies that aim to maximize the coverage of human genetic diversity. Each of the two strategies focus on a specific aspect of genetic diversity. As we recruit more individuals to the HPRC reference, the PSR will continue to explore principled approaches to define and maximize genetic diversity in the refence.

---

## ✨ Requirements

Principled approaches we design require the following characteristics:

* Simple and Transparent
* Reproducible
* Scalable and compatible with a wide range of genetic variant data.

---

## ⚙️ Installation

The following dependencies are required for the algorithm to work:

* Python v3.5 or newer
* Numpy v1.10 or newer
* PLINK v1.9
* Pandas
* Seaborn for the plotting notebook

<!--
> ```bash
# Clone the repository
git clone https://github.com/roohy/hprc_psr_eval.git
cd hprc_psr_eval

# (Optional) Create and activate a virtual environment
python -m venv venv
source venv/bin/activate  # on Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```
-->

---

## Usage

To run MaxVar, first, you need to have the following files ready.

```
main_directory/
├── OOS/
├── Phase1/
├── OOS_phase1_removed/
├── selection/
└── basic_res/

```

* The "OOS" directory should include per-chromosome PLINK binary files containing the sequencing data for all individuals that are not part of the reference.
* The "Phase1" directory contains per-chromosome binary files containing the sequencing data for the previous traunch of selected samples.
* The "OOS_phase1_removed" directory contains the per-chromosome vcf files containing sequencing data for OOS individuals with all variants present in the baseline reference panel removed.
* The "basic_res" directory should be empty.
* In the greed_selection.sh file, you need to set the address so that it points to a full freeze of [1000GP data](http://www.internationalgenome.org).

You can then run the greedy_selection.sh shell script:

```bash
./greedy_selection.sh 
```

Once the script is finished running, you can then run the extract_stat.py script to generate summary statistics of the MaxVar selection.

```bash
python extract_stats.py [chromosome_number] 
```

The script should be run once per chromosome. The python notebooks provide examples of how to use the resulting summary stat data. 

