# Analysis of ESR Spectrum of Tooth Enamel

## Procedure

``normalize.r``

1. Normalization and calibration of magnetic flux density
2. Normalization of spectral intensity
 
``decompose.r``

3. Setting initial values
4. Optimization of parameters

### Installing R

- You need an environment with R. 


---
## Normalization and Calibration (JES X-320, JEOL)

Ensure that the measurement results (``measure.txt``) are stored in the following format:
```
waves=1 length=4096 data=CH1/2 
title= ~
sample= ~
comment= ~
date= 2023/09/03~15:43:03
cf336.800 st30 rmfirst fq100.00 md2.0x0.1
am1.00x1000 tc0.03 uF9441.960000 MHz uP5.00000    mW 
accumu= 60
===== CH1 data Wave No.1 ===== 
mT            Intensity 
    331.29150      -5.74963
    331.29394      -4.99963
    331.29637      -3.74963
    331.29880      -2.87463
    331.30123      -0.62463
    331.30366       1.68787
    331.30609       3.87537
    331.30852       4.81287
  ...
```
Additionally, confirm that the measurement conditions include:
- Both the third and fourth peaks of the Mn marker.


Download the spectrum formatting script ``normalize_JESX320.r`` and place it in the same folder as ``measure.txt``.
Then, execute the following command in the console (Mac terminal or Linux terminal):
```sh
Rscript normalize_JESX320.r measure.txt
```
Once the process is completed, two files will be generated:
- ``measure.txt_calnr.csv``
- ``measure.txt_calnr_ap.csv``

``measure.txt_calnr.csv`` converts the horizontal axis into a dimensionless quantity corresponding to the g-factor at resonance conditions, calculated as $h\nu/\beta H$, and calibrates it using the g-factor of the Mn marker. The vertical axis is normalized by the peak height of the Mn marker.

``measure.txt_calnr_ap.csv`` is essentially the same as ``measure.txt_calnr.csv``, but with a uniformly spaced horizontal axis using linear interpolation rather than the original measurement increments.

To convert multiple files at once, list them as arguments:
```sh
Rscript normalize_JESX320.r measure1.txt measure2.txt measure3.txt
```

---
## Required Spectral File for decompose.r

The formatted spectral file (``measure.txt_calnr.csv``) contains data in the following format:
```
g,a
2.03821022346621,-0.0222308264444024
2.03819480443265,-0.0217374900333399
2.03817944856056,-0.0217688431978222
2.03816409291465,-0.0217657913565938
2.0381487374949,-0.0217455370349838
2.03813338230131,-0.0217424852836534
2.03811796440415,-0.0217652248583653
2.0381026096638,-0.0217449706717819
2.0380872551496,-0.0217849253688014
2.03807190086154,-0.0218592851614268
...
```
The "g" column represents $h\nu/\beta H$ and the "a" column represents intensity.

---
## Running the Analysis

Prepare the following input files:
- Input file: ``input_decompose.in``
- Spectrum file for analysis: ``measure.txt_calnr.csv`` (arbitrary file name)

Download the analysis script ``decompose.r`` and run the following command:
```sh
Rscript decompose.r measure.txt_calnr.csv input_decompose.in output.out
```
Here, ``output.out`` (the name is arbitrary) will be generated, storing the optimized parameter set in the same format as ``input.in`` at the end of the optimization process.

- - - 
## Initial Parameter File
The file ``input_decompose.in`` should contain initial values and optimization control parameters.
It must follow a three-column format, separated by ``:``.

The first column specifies the type of parameter:

The first column specifies the type of parameter:
- ``CNST``：Control parameters for optimization or constants that remain fixed.
- ``FREE``： Initial values of parameters that will be optimized.
- ``LOCK``：Fixed values during optimization. Changing them to ``FREE`` makes them subject to optimization.
- ``NOTE``：Comment lines that are ignored in optimization. However, the three-column structure must be maintained.

Example of ``input_decompose.in``
```
KIND:CONTENT:VALUE
CNST:check initial guess:0
CNST:perturb initial guess:0
CNST:ep perturb:0.1
CNST:smoothing(D=0):2
CNST:number of hot batch:3
CNST:number of cooling batch:5
CNST:number of cold batch:5
CNST:Temperature hot:0.008
CNST:Temperature cold:1e-08
CNST:number of iteration:3000
CNST:g1 carbonate:2.0036
CNST:g2 carbonate:2.0023
CNST:g3 carbonate:1.9975
LOCK:width carbonate:3e-04
FREE:amplitude carbonate:0.749795010390121
CNST:g1 O2-:2.06
CNST:g2 O2-:2.001
CNST:g3 O2-:2.001
LOCK:width O2-:0.00015
FREE:amplitude O2-:0.038445050300495
CNST:g1 Hole center:2.0151
CNST:g2 Hole center:2.007
CNST:g3 Hole center:2.0032
LOCK:width Hole center:0.00018
LOCK:amplitude Hole center:0.26909462
CNST:g1 add1:2.0683
CNST:g2 add1:2.0683
CNST:g3 add1:2.0018
LOCK:width add1:2e-04
FREE:amplitude add1:0
FREE:cut(linear BG):0.00621488788068367
CNST:cut auto 1=auto 0=tune:1
CNST:LinearBG region min1:0.488
CNST:LinearBG region max1:0.489
CNST:LinearBG region min2:0.507
CNST:LinearBG region max2:0.508
LOCK:amplitude Mn:1
LOCK:width scale Mn:1
LOCK:amplitude Native:0
LOCK:broadening Native:1
LOCK:subamp1 Native:0.0
FREE:subamp2 Native:0.116473092366464
FREE:subamp3 Native:0.246247741979634
FREE:subamp4 Native:0.0659017187423335
FREE:subamp5 Native:0.0589587158632329
LOCK:subamp6 Native:0
FREE:subwid1 Native:0.011581297022604
FREE:subwid2 Native:0.00153921848694364
FREE:subwid3 Native:0.00666566029353311
FREE:subwid4 Native:0.00174533950823797
FREE:subwid5 Native:0.00108374968253501
LOCK:subwid6 Native:0.0010812664402403
LOCK:subgiv1 Native:0.495504260426887
LOCK:subgiv2 Native:0.499164153217957
LOCK:subgiv3 Native:0.49715487326045
LOCK:subgiv4 Native:0.498391373852862
LOCK:subgiv5 Native:0.498656382041184
LOCK:subgiv6 Native:0.497145764797701
LOCK:amplitude BG1:0
LOCK:g BG1:2.01711297742095
LOCK:width BG1:0.000783291226211285
LOCK:amplitude BG2:0
LOCK:g BG2:2.0051
LOCK:width BG2:0.001
CNST:refine ite:1000
CNST:initial ep:0.01
CNST:ymax(gross):0.7
CNST:ymax(close-up):0.12
CNST:ymax(2nd deriv):0.12
NOTE:Fitting Weight region min:0.495
NOTE:Fitting weight region max:0.503
CNST:Fitting Weight region min:0.494
CNST:Fitting weight region max:0.502
CNST:Weight inside fitting region:1
CNST:Weight outside fitting region:0
CNST:2nd deriv weight region min:0.498
CNST:2nd deriv weight region max:0.502
CNST:2nd deriv weight:1e-04
CNST:SINK:0
```

