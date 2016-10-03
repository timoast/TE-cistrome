## Transposable element transciption factor bindings

Aim: to identify TF binding sites within Arabidopsis TEs, in order to make testable predictions about the environmental conditions leading to TE transcritional activation.

### Methods

Certain TEs in the Arabidopsis are expected to contain peaks, but these are absent in the published peaks from O'Malley et al. In particular the *ONSEN* transposons should have a HSFA1 peak (which is now well documented), but don't. Instead these regions have a complete lack of mapped reads, which I believe is due to the discarding of multimapping reads during alignment. To address this issue I will need to download and process the data again myself, and either allow multimapping reads, or map to a single representative of each TE family.

#### Step 1: Download data

```bash
sh DataAquisition/download_from_sra.sh -f RawData/SraAccList.txt
```

#### Step 2: Map

```bash
sh ProcessingCode/map.sh
```

#### Step 3: Call peaks

```bash
sh ProcessingCode/call_peaks.sh
```

#### Step 4: Find conserved bindings in TE families


### References

O'Malley RC, Huang S-SC, Song L, Lewsey MG, Bartlett A, Nery JR, et al. Cistrome and Epicistrome Features Shape the Regulatory DNA Landscape. Cell. 2016.[doi:10.1016/j.cell.2016.04.038](http://dx.doi.org/10.1016/j.cell.2016.04.038)
