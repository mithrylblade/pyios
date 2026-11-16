# pyios

Reference equations for **impulse oscillometry (IOS)**. Given a subject's age, height, weight and
sex, `pyios` returns the predicted value, % predicted, z-score and lower/upper limits of normal
(LLN/ULN) for an oscillometry parameter, using one of five published reference sets.

| `equation` | Reference | Population | Ages | Method |
|---|---|---|---|---|
| `moitra2020` | Moitra et al. 2020 | India | 18–88 | log-linear (age, weight, height) |
| `oostveen2013` | Oostveen et al. 2013 | Europe/Australia (forced oscillation) | 18–84 | log-linear (age, weight, height) |
| `schulz2013` | Schulz et al. 2013 | Germany | 45–89 | quantile regression (P5/P50/P95) |
| `berger2021` | Berger et al. 2021 | USA | 21–85 | log-linear (height, BMI) |
| `gochicoa2023` | Gochicoa-Rangel et al. 2023 | Mexico | 2.7–90 | segmented regression (sex, age, 1/height, 1/BMI) |

Full study metadata (population, devices, DOIs, validity ranges) is available from
`reference.metadata`.

## Install

```bash
pip install git+https://github.com/mithrylblade/pyios.git
```

This needs [git](https://git-scm.com/) installed, and it installs `numpy` and `pandas` automatically. To install from a local clone instead:

```bash
git clone https://github.com/mithrylblade/pyios.git
cd pyios
pip install .
```

## Check the installation

```bash
python -c "from pyios import OscillometryReference; print(OscillometryReference('gochicoa2023'))"
```

## Single subject

```python
from pyios import OscillometryReference

reference = OscillometryReference(equation="gochicoa2023")

result = reference.calculate(age=50, height=175, weight=80, sex="male", param="R5", value=0.30)
result.predicted, result.percent_predicted, result.z_score, result.lln, result.uln, result.flag
result.to_dict()
```

`value` (the measured result) is optional. Without it you get the predicted value and limits;
`percent_predicted` and `z_score` are `None` and this is not flagged.

## Whole DataFrame

```python
results = reference.compute(
    df,
    param="R5",
    sex_col="sex", age_col="age", height_col="height", weight_col="weight",  # defaults
    param_col="R5",       # column of measured values (defaults to a column named like `param`)
)
df = df.join(results)
```

Returns a DataFrame with the same index and columns `predicted`, `percent_predicted`, `z_score`,
`lln`, `uln`, `flag` (`prefix="R5_"` prefixes them). Non-numeric inputs become missing and are
flagged; bad rows never raise.

## Flexible inputs

* **Equation**: case, spaces, underscores and the year are optional – `"gochicoa"`, `"Gochicoa_2023"`,
  `"gochicoa 2023"`, `"GOCHICOA2023"` (`"shulz"` is accepted for Schulz).
* **Sex**: `0 "0" "f" "F" "female" "Female" "FEMALE"` → female; `1 "1" "m" "M" "male" "Male" "MALE"` → male.
  Anything else is flagged.
* **Parameter**: case-insensitive; `R5-R20`, `r5_r20`, `R5–R20` and `R5-20` are equivalent.
  `AX` means AX integrated from 5 Hz; Oostveen also offers `AX4`.

## Units

```python
OscillometryReference("oostveen2013", age_unit="years", height_unit="cm",
                      weight_unit="kg", result_unit="kPa")
```

| argument | default | options |
|---|---|---|
| `age_unit` | `years` | `years`, `months`, `weeks`, `days` |
| `height_unit` | `cm` | `cm`, `m`, `mm`, `in` |
| `weight_unit` | `kg` | `kg`, `g`, `lb` |
| `result_unit` | `kPa` | `kPa`, `hPa`, `cmH2O` |

Inputs are converted to whatever each equation requires (e.g. Oostveen/Berger use height in
metres, Berger works in cmH2O, Oostveen in hPa). Measured values passed in, and all returned
predicted values and limits, use `result_unit`. Non-pressure parameters (`Fres` in Hz, ratios in
%) are never scaled. Flag messages always quote age in years, height in cm and weight in kg.

## Flag

`flag` is `None` when everything is valid. Otherwise it is a `"; "`-joined message, for example:

* `parameter 'R10' is not available in reference set 'schulz2013' (available: ...)` – all outputs `None`
* `age 30 years below reference minimum 45 years` – values are still computed
* `height 220 cm above reference maximum 203 cm`, `weight ...`, `BMI ... above reference maximum 30 kg/m2`
* `sex missing or not recognised`, `height missing or not positive`, `age missing or not positive` – outputs `None`
* `z-score not available: schulz2013 equations are quantile-based ...`
* `measured value outside the domain of the log-transformed equation: z-score not calculated`

Age limits are the validated range of the specific equation (sex-specific for Oostveen). Height,
weight and BMI limits come from `reference_metadata.csv` (see caveats).

## How the statistics are computed

* **Log-linear / linear references** (Moitra, Oostveen, Berger): the regression is evaluated in its
  transformed space, `z = (T(observed) − T(predicted)) / SD` using the published RMSE/RSD, and
  LLN/ULN are the values at z = ∓1.6449, back-transformed. Offsets (Berger `+1`/`+4`,
  Oostveen `4 − X`) are handled, so reactance limits are correct (LLN of X is the more negative
  value). z-scores are signed so that z > 0 always means "observed above predicted".
* **Gochicoa**: the age segment is chosen from the published age bands (ages falling in the 0.01-year
  gaps between bands go to the next segment). The sex term uses male = 1, female = 0.
* **Schulz** (quantile regression): predicted = P50 equation, LLN = P5, ULN = P95. The paper provides
  no residual SD, so `z_score` is `None` (and flagged). % predicted is available.
* `% predicted = 100 × observed / predicted`.

## Parameters per reference

| | Parameters |
|---|---|
| `moitra2020` | R5, R20, R5-R20, Z5, X5, AX, Fres |
| `oostveen2013` | R4…R26 (14 frequencies), Rmean, X4, X5, X6, X8, X10, X12, X14, Fres, AX4, AX5 |
| `schulz2013` | R5, R20, R5-R20, (R5-R20)/R20, Z5, X5, AX, Fres |
| `berger2021` | R5, R10, R15, R20, R5-R15, R5-R20, X5, X10, AX, Fres |
| `gochicoa2023` | R5, R10, R15, R20, R5-R20, (R5-R20)/R5, X5, X10, X15, X20, AX, Fres |

`reference.parameters` lists them programmatically; `OscillometryReference.available_equations()` lists the sets.

## Data

The equations and metadata ship as CSV files in `src/pyios/data/`
(`reference_equations.csv`, `reference_metadata.csv`). Each equation row records coefficients,
transforms, units, validity ages, RMSE/RSD and the source table, so equations can be audited or
extended without touching the code.

## Caveats

* Research software; not validated for clinical decision-making.
* Reference sets differ in population, device and units; results from different sets are not
  interchangeable. Oostveen is forced-oscillation data from several devices, not strictly IOS.
* Height/weight/BMI limits are overall (both sexes) and, for Schulz, the authors' recommended
  limits; Moitra's BMI limit is an exclusion criterion; Gochicoa's ranges span all ages, so a value
  can be inside the limits yet implausible for a given age.
* Moitra's published coefficients are rounded to 1–3 significant figures, which limits precision.

## R package

An R port lives in [`r/`](r/README.md) (package `rios`, base R only). It has the same API
(`OscillometryReference(...)$calculate()` / `$compute()`), units handling, flag messages and
numerical results, and is tested against golden vectors generated from this Python package.
The CSV data files are duplicated in `r/inst/extdata/`; run `r/data-raw/sync_data.sh` after
changing the equations (the Python tests fail if the copies differ).

## Tests

```bash
python -m unittest discover -s tests      # or: pytest
```

## How to cite

If you use `pyios` (or the R port `rios`) in your work, please cite the software. GitHub's
**"Cite this repository"** button (right-hand sidebar) is generated from
[`CITATION.cff`](CITATION.cff).

> Moolla S. (2026). *pyios: Reference equations for impulse oscillometry* (Version 0.1.0)
> [Computer software]. https://github.com/mithrylblade/pyios
> *(DOI to be added: Zenodo DOI `10.5281/zenodo.XXXXXXX` &ndash; placeholder)*

BibTeX:

```bibtex
@software{moolla2026pyios,
  author  = {Moolla, Saadiq},
  title   = {pyios: Reference equations for impulse oscillometry},
  year    = {2026},
  version = {0.1.0},
  url     = {https://github.com/mithrylblade/pyios},
  license = {MIT}
}
```

**Please also cite the original paper(s) for the reference equation set(s) you used**
([References](#references) below). The equations are the work of those authors; this software
only implements them.

## References

The reference equations implemented in this package are taken from the following publications
(the equation id used by the package is shown after each reference). BibTeX entries are in
[`references.bib`](references.bib); their keys equal the equation ids.

1. Moitra S, Moitra S, Ghosh AK, Sengupta S, Das PK, Das A, et al. Reference values of impulse oscillometry (IOS) for healthy Indian adults. Int J Tuberc Lung Dis. 2020;24(5):536–539. [doi:10.5588/ijtld.19.0796](https://doi.org/10.5588/ijtld.19.0796) &mdash; `moitra2020`
2. Oostveen E, Boda K, van der Grinten CPM, James AL, Young S, Nieland H, et al. Respiratory impedance in healthy subjects: baseline values and bronchodilator response. Eur Respir J. 2013;42:1513–1523. [doi:10.1183/09031936.00126212](https://doi.org/10.1183/09031936.00126212) &mdash; `oostveen2013`
3. Schulz H, Flexeder C, Behr J, Heier M, Holle R, Huber RM, et al. Reference values of impulse oscillometric lung function indices in adults of advanced age. PLoS ONE. 2013;8(5):e63366. [doi:10.1371/journal.pone.0063366](https://doi.org/10.1371/journal.pone.0063366) &mdash; `schulz2013`
4. Berger KI, Wohlleber M, Goldring RM, Reibman J, Farfel MR, Friedman SM, et al. Respiratory impedance measured using impulse oscillometry in a healthy urban population. ERJ Open Res. 2021;7:00560-2020. [doi:10.1183/23120541.00560-2020](https://doi.org/10.1183/23120541.00560-2020) &mdash; `berger2021`
5. Gochicoa-Rangel L, Martínez-Briseño D, Guerrero-Zúñiga S, Contreras-Morales J, Arias-Jiménez D, Del-Río-Hidalgo R, et al. Reference equations using segmented regressions for impulse oscillometry in healthy subjects aged 2.7–90 years. ERJ Open Res. 2023;9:00503-2023. [doi:10.1183/23120541.00503-2023](https://doi.org/10.1183/23120541.00503-2023) &mdash; `gochicoa2023`

## Funding

This work was supported by the National Institute for Health and Care Research (NIHR)
(award 1556644) using UK international development funding from the UK Government to support
global health research, and by the Research Networks for Health Innovations in Sub-Saharan
Africa (RHISSA) funding initiative of the German Federal Ministry of Education and Research
(BMBF; funding code 01KA2223B).

The views expressed are those of the author(s) and not necessarily those of the NIHR, the
UK Government, or the BMBF.

## License

Released under the [MIT License](LICENSE).
