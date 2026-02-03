# BactoML gating: TCC, ICC, HNA (bactoml-master)

Below is how **TCC**, **ICC**, and **HNA** are computed in `bactoml-master`, based on the tutorial pipelines and tests.

## Where this is defined
- `bactoml-master/tutorial/MCCD_pipeline_example.py` (main reference pipeline)
- `bactoml-master/tutorial/pipeline_example.py` (simpler pipeline)
- `bactoml-master/tutorial/BactoML_Spiking_Experiment_Analysis.ipynb` (same gates + formulas, incl. ICC/DCC mode)
- `bactoml-master/bactoml/tests/test_df_pipeline.py` (unit-test definitions of gates)

## Common preprocessing
- **tlog transform** is applied before gating:
  - `transform('tlog', channels=['FL1','FL2','SSC'], th=1, r=1, d=1, auto_range=False)`

## TCC (Total Cell Count)
**Channels used:** `FL1` and `FL2`.

**Gate:** `PolyGate` in the FL1–FL2 plane. Two example definitions used in repo:
- `pipeline_example.py` and tests:
  - `PolyGate([[3.7,0],[3.7,3.7],[6.5,6],[6.5,0]], ['FL1','FL2'])`
- `MCCD_pipeline_example.py` and Spiking notebook (parameterized with `p_FL1=4`):
  - `PolyGate([[p_FL1,0.05],[p_FL1,3.2],[6.5,5.7],[6.5,0.05]], ['FL1','FL2'])`

**Count:**
- `TCC_abs` is computed by counting events after the TCC gate:
  - `event_counter_step = lambda x: x.shape[0]`

**Normalization (cells/L):**
- `VOL` is read from metadata: `float(x.get_meta()['$VOL']) * 1E-6`
- `TCC = TCC_abs / VOL`
  - In the LDC variant (see below): `TCC = (ICC_abs + DCC_abs) / VOL`

## ICC (Intact Cell Count)
**Channels used:** `FL1` and `FL2`.

**Gate:** `PolyGate` in the FL1–FL2 plane, defined as `icc_gate_step`.
- In the tutorial notebook, **ICC and TCC use the same polygon** (parameterized by `p_FL1`).

**Count and normalization:**
- The notebook defines the formulas for the LDC mode:
  - `ICC = ICC_abs / VOL`
  - `TCC = (ICC_abs + DCC_abs) / VOL`
- `ICC_abs` is the event count after applying the ICC gate (this is implied by the `icc_gate_step` + `event_counter_step` pattern used elsewhere).

## HNA (High Nucleic Acid)
**Channel used:** `FL1` only.

**Gate:** `ThresholdGate` on `FL1`, **above** a threshold.
- Examples used in repo:
  - `ThresholdGate(5.1, 'FL1', 'above')` (tests + `pipeline_example.py`)
  - `ThresholdGate(4.8, 'FL1', 'above')` (MCCD + Spiking notebook)

**Count:**
- `HNAC` is computed by counting events after the HNA gate:
  - `HNAC = event_count_step(hna_gate_step(...))`

**Normalization:**
- Standard mode: `HNAP = HNAC / TCC_abs * 100`
- LDC mode: `HNAP = HNAC / ICC_abs * 100`

**Order of gating:**
- In the Spiking notebook, **HNA gating is applied after TCC gating**:
  - `Pipeline([tcc_gate_step, hna_gate_step, event_counter_step])`
- In the simpler pipeline and tests, HNA is applied to the preprocessed measurement (after `tlog` and `tcc_gate` in the preprocessing pipeline).

## How gating is performed (mechanism)
- Gating uses **FlowCytometryTools** gates via `FCMeasurement.gate(...)`.
- In BactoML, each gate is wrapped in `DFLambdaFunction` and placed inside a `sklearn` `Pipeline`:
  - Example pattern: `DFLambdaFunction(lambda x: x.gate(PolyGate(...)))`
- Event counts are then extracted with `DFLambdaFunction(lambda x: x.shape[0])`.

## Notes on “gain” vs “gating”
- I did **not** find any “gain” calculation in `bactoml-master`.
- If you meant **gating**, it is performed exactly as described above with `PolyGate`/`ThresholdGate` on `FL1`/`FL2` (and `tlog` preprocessing).
