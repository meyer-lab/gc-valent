# 4/23: Plate 1

Gates tentatively determined relative to shapes and percentages outlined in Farhat et al. Figure S3. Here, magnitudes have been transformed by log10, but do not appear consistent with S3.


```python
import os
import string
import numpy as np
import pandas as pd
from FlowCytometryTools import ThresholdGate, PolyGate
from matplotlib import pyplot as plt
from matplotlib import cm
from flow import importF
```


```python
path_ = os.path.abspath('')
```


```python
sample1A, _ = importF(path_ + "/data/PBMC-rec-quant/04-23/Plate 1/Plate 1 - Panel 1 IL2R/", "A")
sample2B, _ = importF(path_ + "/data/PBMC-rec-quant/04-23/Plate 1/Plate 1 - Panel 2 IL2R/", "B")
sample3C, _ = importF(path_ + "/data/PBMC-rec-quant/04-23/Plate 1/Plate 1 - Panel 3 IL2R/", "C")
```


```python
def combineWells(samples, channels_):
    """Accepts sample array returned from importF, and array of channels, returns transformed combined well data"""
    combinedSamples = samples[0]
    for sample in samples[1:]:
        combinedSamples.data = combinedSamples.data.append(sample.data)
    t_combinedSamples = combinedSamples.transform('tlog', channels=channels_)
    return t_combinedSamples
```


```python
panel1 = combineWells(sample1A, ['VL6-H', 'VL4-H', 'BL1-H', 'VL1-H', 'BL3-H'])
panel2 = combineWells(sample2B, ['VL4-H', 'BL3-H', 'FSC-H'])
panel3 = combineWells(sample3C, ['VL6-H', 'VL4-H', 'BL3-H'])
```

## C) T helper and T regulatory cells


```python
cd3cd4_gate = PolyGate([(5.6e03, 5.5e03), (7.9e03, 6.0e03), (7.9e03, 7.2e03), (5.5e03, 7.2e03)], ('VL4-H', 'VL6-H'), region='in', name='cd3cd4')
_ = panel1.plot(['VL4-H','VL6-H'], gates=cd3cd4_gate, cmap=cm.jet) #CD3, CD4
plt.title("Singlet Lymphocytes")
plt.xlabel("CD3")
plt.ylabel("CD4")
print("Percent gated as CD3+CD4+: ", panel1.gate(cd3cd4_gate).counts/panel1.counts*100., "%")
```

The above shape and percentage are consistent with S3(C).


```python
samplecd3cd4 = panel1.gate(cd3cd4_gate)
thelp_gate = PolyGate([(3.1e03, 4.5e03), (3.5e03, 4.5e03), (6.4e03, 6.2e03), (6.4e03, 6.8e03), (3.1e03, 6.8e03)], ('VL1-H', 'BL1-H'), region='in', name='thelp')
treg_gate = PolyGate([(4.9e03, 4e03), (5.6e03, 5.6e03), (6.8e03, 6.3e03), (8e03, 3.8e03), (8e03, 2.0e3), (4.9e03, 2.0e3)], ('VL1-H', 'BL1-H'), region='in', name='treg')
_ = samplecd3cd4.plot(['VL1-H','BL1-H'], gates=[thelp_gate,treg_gate], cmap=cm.jet) #CD3, CD4
plt.title("CD3+CD4+ cells")
plt.xlabel("CD25")
plt.ylabel("CD127")
print("Percent gated as T helper: ", samplecd3cd4.gate(thelp_gate).counts/samplecd3cd4.counts*100., "%")
print("Percent gated as Treg: ", samplecd3cd4.gate(treg_gate).counts/samplecd3cd4.counts*100., "%")
```

The above shape is similar to S3(C), but the percentage of Tregs is lower than 6%.


```python
samplethelp = samplecd3cd4.gate(thelp_gate)
sampletreg = samplecd3cd4.gate(treg_gate)
```


```python
memthelpgate = ThresholdGate(6750.0, ['BL3-H'], region='below') & ThresholdGate(3800.0, ['BL3-H'], region='above')
naivethelpgate = ThresholdGate(7250.0, ['BL3-H'], region='above') & ThresholdGate(8700.0, ['BL3-H'], region='below')
_ = samplethelp.plot(['BL3-H'], gates=[memthelpgate,naivethelpgate])
plt.xlabel("CD45RA")
plt.title("T helper")
print("Percent gated as Memory T helper: ", samplethelp.gate(memthelpgate).counts/samplethelp.counts*100., "%")
print("Percent gated as Naive T helper: ", samplethelp.gate(naivethelpgate).counts/samplethelp.counts*100., "%")
```

The above shape and percentages are consistent with S3(C).


```python
memtreggate = ThresholdGate(6050.0, ['BL3-H'], region='below') & ThresholdGate(3600.0, ['BL3-H'], region='above')
naivetreggate = ThresholdGate(6700.0, ['BL3-H'], region='above') & ThresholdGate(8600.0, ['BL3-H'], region='below')
_ = sampletreg.plot(['BL3-H'], gates=[memtreggate,naivetreggate])
plt.xlabel("CD45RA")
plt.title("T reg")
print("Percent gated as Memory Treg: ", sampletreg.gate(memtreggate).counts/sampletreg.counts*100., "%")
print("Percent gated as Naive Treg: ", sampletreg.gate(naivetreggate).counts/sampletreg.counts*100., "%")
```

Shape of histogram differs from S3(C), so these gates based on percentages.

## D) NK Cells


```python
nk_gate = PolyGate([(3800, 5700), (6700, 8200), (7200, 8100), (5200, 5300)], ('VL4-H', 'BL3-H'), region='in', name='nk')
nkt_gate = PolyGate([(5500, 5500), (6000, 5000), (8000, 6500), (6900, 7400)], ('VL4-H', 'BL3-H'), region='in', name='nkt')
_ = panel2.plot(['VL4-H', 'BL3-H'], gates=[nk_gate,nkt_gate], gate_colors=['red','red'], cmap=cm.jet) #CD56 vs. CD3
samplenk = panel2.gate(nk_gate)
samplenkt = panel2.gate(nkt_gate)
print("Percent gated as NK: ", samplenk.counts/panel2.counts*100., "%")
print("Percent gated as NKT: ", samplenkt.counts/panel2.counts*100., "%")
plt.title("Singlet Lymphocytes")
plt.xlabel("CD3")
plt.ylabel("CD56")
```

Gated according to the projected shape of S3(D) without compensation. The percentage NK is slightly higher than 7%.

## E) CD8+ T cells


```python
cd8_gate = PolyGate([(5700, 5600), (5600, 7800), (7900, 7800), (7800, 6600)], ('VL4-H', 'VL6-H'), region='in', name='cd8')
_ = panel3.plot(['VL4-H', 'VL6-H'], gates=cd8_gate, cmap=cm.jet) #CD8 vs. CD3
plt.title("Singlet Lymphocytes")
plt.xlabel("CD3")
plt.ylabel("CD8")
print("Percent gated as CD3+CD8+: ", panel3.gate(cd8_gate).counts/panel3.counts*100., "%")
```

The shape and percentage are consistent with S3(E).


```python
cd3cd8sample = panel3.gate(cd8_gate)
memcd3cd8gate = ThresholdGate(7100.0, ['BL3-H'], region='below') & ThresholdGate(4200.0, ['BL3-H'], region='above')
naivecd3cd8gate = ThresholdGate(7500.0, ['BL3-H'], region='above') & ThresholdGate(8700.0, ['BL3-H'], region='below')
_ = cd3cd8sample.plot(['BL3-H'], gates=[memcd3cd8gate,naivecd3cd8gate])
plt.xlabel("CD45RA")
plt.title("CD3+CD8+ Cytotoxic T cells")
print("Percent gated as Cytotoxic Memory: ", cd3cd8sample.gate(memcd3cd8gate).counts/cd3cd8sample.counts*100., "%")
print("Percent gated as Cytotoxic Naive: ", cd3cd8sample.gate(naivecd3cd8gate).counts/cd3cd8sample.counts*100., "%")
```

The shape is similar to S3(E), but the percentage of Cytotoxic Memory is slightly greater than the reference (39%).


```python

```
