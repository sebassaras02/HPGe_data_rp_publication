# HPGe detector Calibration and Data Analysis


This repository contains the GEANT4 applications built, experimental data, and jupyter notebooks employed for the data analysis.

**Table of Contents**
* **e_results folder**: has the experimental data for detector calibration.

* **s_results folder**: has the experimental data for detector calibration obtained from the GEANT4 application.

* **hpge_data_treat.ipynb** is the data treatment, data analisis, and data visualization of the results.

* **d_funciones.py** has functions and classes used to modularize the code.

### Data 

----
- The experimental data come from several experiments done in the Nuclear Science Department at Escuela Politécnica Nacional, Quito - Ecuador.
- The simulated data come from simulations run on the GEANT4 application. 
### HPGe detector
---
**1. What is a HPGe detector?**
A HPGe detector is a Hyper Pure Germanium Detector. This kind of detector is used in gamma spectroscopy to identify different isotopes in a sample. 

**2. How is calibrated an HPGe detector?**
The HPGe detector is calibrated through standars with a known activity. It is counted how many particles is the detector able to register effectively. 
The calibration is done with the Full Energy Peak Efficiency (FEPE). 

$$ FEPE = \frac{AP}{E*P} $$

Where:
- AP: is the area of a photopeak for a energy window
- E: is the activity of the radioactive source
- P: is the probability of the emission

How a HPGe detector looks like:

![P-type-400](https://user-images.githubusercontent.com/82113558/189552768-dd2538f3-13c4-468e-a219-9d8dfbe7f3d9.jpg)
                
### FlowChart of the Data Treatment
---
```mermaid 
  graph TD
    A[Collect experimental data] --> B(Collect simulated data)
    B --> C{FEPE calibration for experimental data}
    C --> D[Is the calibration succesfull]
    D -->|Yes| E[FEPE calibration for simulated data]
    D -->|No| F[Try all again]
    E --> G[Scale simulated data]
    G --> H[Compare simulated and experimental results]

```

### Contacts:
---
- [Twitter](https://twitter.com/sarasti_seb)
- [LinkedIn](https://linkedin.com/in/sebastiansarasti)
- [ResearchGate](https://www.researchgate.net/profile/Sebastian-Sarasti-2)
