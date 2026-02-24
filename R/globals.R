# Suppress R CMD check notes for dplyr/tidyr column references
utils::globalVariables(c(
  # findrelation.R
  "pseudoR2.mcfadden",
  "pseudoR2.cox_snell",
  "pseudoR2.nagelkerke",
  "model_type",
  "model_eq",
  "R2",
  "value",
  "Resp",
  "Var",
  "Var_sq",
  "term",
  "statistic",
  "p.value",
  "isConv",
  "finTol",
  "n_breakpoints",
  "mc_fadden",
  "cox_and_snell_ml",
  "nagelkerke_cragg_and_uhler",
  # netSEMp1.R
  "adj.r.squared",
  "best_model",
  "r.squared",
  "x",
  "mfExt",
  "resp",
  # netSEMp2.R
  "IAdj-R2",
  "GAdj-R2",
  # plot files
  "adjRSqr",
  "endogenous",
  "exogenous",
  "Variable",
  "conp1.c",
  "regex",
  # paths.R
  "p_value",
  "R2adj",
  "var1",
  "var2",
  "model",
  "path",
  "min_R2",
  # var_selection.R
  "row_number"
))