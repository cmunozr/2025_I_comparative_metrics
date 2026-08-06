# Modeling & Analysis Tracker

## August 4th. Task Checklist

### 1. Spatial Sensitivity Analysis (The Northern Region "Stress Test")
- [ ] **Configure Model 16 (M16) as Baseline**
  - Proceed with M16, which incorporates survey-effort variables to control for observational biases.
  - When running predictions, hold sampling effort variables (such as survey duration) constant at their mean so the final biodiversity metrics are independent of sampling biases .
- [ ] **Run North-South / East-West Performance Check**
  - Analyze how the baseline model's performance metrics vary geographically across Finland.
- [ ] **Execute the Southern Training / Northern Prediction Test**
  - Train the model exclusively on data from the South, then predict species occurrences/abundances in the North to verify if predictive relationships are indeed significantly worse up north.
- [ ] **Create the "Excluding the North" Model Filter**
  - Set up a filter to restrict the control-treatment dataset, removing the northern region (Lapland/specific Eli regions).
  - Verify that the treatment-control difference (the Metso signal) remains intact or is strengthened once the slower-growing, structurally distinct northern sites are removed .

### 2. Biodiversity Metrics & Artificial Scenarios
- [ ] **Run prediction over the no-north region**
  - Calculate complete XData for M016 
  - Calculate species richness and aggregate abundance across your 67 forest bird species . Do not run individual gradient analyses for each separate species to keep the final paper focused .
- [ ] **Run metrics**
  - Set up gradient simulations (e.g., varying canopy cover or stand age) .
  - Compare "Type 1" co-varying predictions (where other variables covary) against models where other covariates are locked at their mean to observe true structural effect directions .
 - [ ] **Build Artificial "Extreme Contrast" Scenarios**
  - Manually manipulate key structural variables (such as increasing or decreasing certain forest parameters) to create extreme simulated contrasts between Metso and control sites .

### 3. Species Classification Analysis
- [ ] **Integrate Species Classification List**
  - Categorize the 67 species.
- [ ] **Analyze Subgroup Responses**
  - Evaluate and compare the differences in treatment effects (Metso vs. control) for:
    - Forest specialists .
    - Forest generalists .
    - Declining bird species.

### 4. Methodological Draft & Asynchronous Sharing
- [ ] **Write Up Modeling Methodology**
  - Draft the methodological description of the HMSC modeling process .
  - Document the rationale and results for how the northern spatial constraints and geographic differences were handled .
- [ ] **Draft the Project Deliverable**
  - Instead of booking a live meeting, compile your newly calculated biodiversity metrics, gradient simulations, and drafted methodology into a **PowerPoint presentation or a paper draft** .
- [ ] **Share Deliverables Asynchronously**
  - Distribute the draft/PowerPoint to the team for review so everyone can look it over before the next follow-up call .